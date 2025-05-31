#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import subprocess
import time
import logging
import sys
from datetime import datetime
import shutil

def check_path(path):
    """检查并标准化路径"""
    if path.startswith('~'):
        path = os.path.expanduser(path)
    return os.path.abspath(path)

def run_script(script_path, args=None, script_description=None):
    """运行Python脚本并捕获输出"""
    cmd = [sys.executable, script_path]
    if args:
        cmd.extend(args)
    
    if script_description:
        print(f"\n{'='*80}")
        print(f"开始执行: {script_description}")
        print(f"{'='*80}")
    
    start_time = time.time()
    
    try:
        subprocess.run(cmd, check=True)
        success = True
    except subprocess.CalledProcessError as e:
        print(f"错误: 脚本 {script_path} 执行失败，代码: {e.returncode}")
        success = False
    except Exception as e:
        print(f"错误: 执行 {script_path} 时遇到意外错误: {str(e)}")
        success = False
    
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    if script_description:
        status = "成功" if success else "失败"
        print(f"\n{script_description} {status}完成，耗时: {int(minutes)}分{int(seconds)}秒")
    
    return success

def setup_logging():
    """设置日志记录"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"custom_pipeline_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ]
    )
    return log_filename

def main():
    parser = argparse.ArgumentParser(description="运行完整的数据处理流程：生成口袋和自定义EHIGN格式数据")
    parser.add_argument("--train_pdb", default="train_pdb", help="包含蛋白质-配体数据的根目录")
    parser.add_argument("--no_pocket", action="store_true", help="跳过口袋生成步骤")
    parser.add_argument("--no_ehign", action="store_true", help="跳过EHIGN格式生成步骤")
    parser.add_argument("--no_rdkit", action="store_true", help="跳过RDKIT文件生成")
    parser.add_argument("--no_dgl", action="store_true", help="跳过DGL图文件生成")
    parser.add_argument("--no_parallel", action="store_true", help="禁用并行处理")
    parser.add_argument("--cutoff", type=float, default=5.0, help="口袋距离阈值(埃)")
    parser.add_argument("--dis_threshold", type=float, default=5.0, help="配体-口袋原子间距阈值(埃)")
    parser.add_argument("--retry_failed", action="store_true", help="尝试重新处理之前失败的样本")
    parser.add_argument("--error_report", default="error_report.txt", help="记录错误样本的报告文件")
    parser.add_argument("--max_attempts", type=int, default=1, help="尝试处理失败样本的最大次数")
    
    args = parser.parse_args()
    
    # 设置日志
    log_filename = setup_logging()
    logging.info(f"开始完整处理流程，日志文件: {log_filename}")
    
    # 记录执行参数
    for arg, value in vars(args).items():
        logging.info(f"参数 {arg}: {value}")
    
    # 标准化路径
    train_pdb_dir = check_path(args.train_pdb)
    logging.info(f"处理目录: {train_pdb_dir}")
    
    # 错误报告文件
    error_report_file = args.error_report
    
    # 如果启用了重试失败样本，但错误报告文件不存在，发出警告
    if args.retry_failed and not os.path.exists(error_report_file):
        logging.warning(f"重试模式启用，但错误报告文件不存在: {error_report_file}")
        args.retry_failed = False

    # 如果启用重试，从错误报告中加载失败样本列表
    failed_samples = set()
    if args.retry_failed and os.path.exists(error_report_file):
        logging.info(f"加载错误报告文件以重试失败样本: {error_report_file}")
        with open(error_report_file, 'r') as f:
            for line in f:
                if line.strip() and ':' in line:
                    sample_path = line.split(':', 1)[0].strip()
                    failed_samples.add(sample_path)
        logging.info(f"加载了 {len(failed_samples)} 个失败样本进行重试")

    # 使用临时文件记录当前运行的错误样本
    current_errors_file = f"current_errors_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    pipeline_start_time = time.time()
    
    # 步骤1: 生成口袋文件
    if not args.no_pocket:
        pocket_args = [
            "--base_dir", train_pdb_dir,
            "--cutoff", str(args.cutoff)
        ]
            
        if args.no_parallel:
            pocket_args.append("--no_parallel")
            
        if args.retry_failed:
            # 添加错误报告文件，用于记录失败样本
            pocket_args.extend(["--error_report", current_errors_file])
            # 仅处理之前失败的样本
            pocket_args.extend(["--only_process_failed", error_report_file])
        
        success = run_script(
            "./generate_protein_pockets.py", 
            pocket_args,
            "口袋文件生成"
        )
        
        if not success:
            logging.error("口袋文件生成失败，流程终止")
            return 1
    else:
        logging.info("跳过口袋文件生成步骤")
    
    # 步骤2: 生成自定义EHIGN格式文件
    if not args.no_ehign:
        ehign_args = [
            "--base_dir", train_pdb_dir,
            "--cutoff", str(args.cutoff),
            "--dis_threshold", str(args.dis_threshold)
        ]
        
        if args.no_rdkit:
            ehign_args.append("--no_rdkit")
        
        if args.no_dgl:
            ehign_args.append("--no_dgl")
            
        if args.no_parallel:
            ehign_args.append("--no_parallel")
            
        if args.retry_failed:
            # 添加错误报告文件，用于记录失败样本
            ehign_args.extend(["--error_report", current_errors_file])
            # 仅处理之前失败的样本
            ehign_args.extend(["--only_process_failed", error_report_file])
        else:
            # 记录错误样本到错误报告文件
            ehign_args.extend(["--error_report", current_errors_file])
        
        success = run_script(
            "./generate_ehign_custom.py", 
            ehign_args,
            "自定义EHIGN格式生成"
        )
        
        if not success:
            logging.error("自定义EHIGN格式生成失败，流程终止")
            return 1
    else:
        logging.info("跳过EHIGN格式生成步骤")
    
    # 合并当前错误报告到主错误报告文件
    if os.path.exists(current_errors_file):
        if os.path.exists(error_report_file):
            # 读取现有错误报告
            existing_errors = set()
            with open(error_report_file, 'r') as f:
                for line in f:
                    existing_errors.add(line.strip())
            
            # 读取当前错误
            current_errors = set()
            with open(current_errors_file, 'r') as f:
                for line in f:
                    current_errors.add(line.strip())
            
            # 合并错误报告
            combined_errors = existing_errors.union(current_errors)
            
            # 写回错误报告文件
            with open(error_report_file, 'w') as f:
                for error in sorted(combined_errors):
                    f.write(f"{error}\n")
                    
            logging.info(f"更新了错误报告文件 {error_report_file} 中的 {len(combined_errors)} 条记录")
        else:
            # 如果主错误报告不存在，直接使用当前错误报告
            shutil.copyfile(current_errors_file, error_report_file)
            logging.info(f"创建了新的错误报告文件: {error_report_file}")
        
        # 删除当前错误报告
        try:
            os.remove(current_errors_file)
        except:
            pass
    
    # 计算总处理时间
    pipeline_time = time.time() - pipeline_start_time
    hours, remainder = divmod(pipeline_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    logging.info(f"完整流程处理完成，总耗时: {int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}")
    print(f"\n{'='*80}")
    print(f"数据处理流程已完成！")
    print(f"现在您可以在DrugBAN中使用生成的数据文件")
    
    # 如果存在错误报告，输出统计信息
    if os.path.exists(error_report_file):
        failed_count = 0
        with open(error_report_file, 'r') as f:
            failed_count = sum(1 for _ in f)
        print(f"警告：有 {failed_count} 个样本处理失败，详见 {error_report_file}")
        print(f"可以使用 --retry_failed 参数重新处理失败的样本")
        
    print(f"总处理时间: {int(hours):02d}小时{int(minutes):02d}分钟{int(seconds):02d}秒")
    print(f"{'='*80}")
    
    return 0

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code) 