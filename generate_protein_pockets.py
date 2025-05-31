#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import glob
import subprocess
import multiprocessing
import logging
import time
from tqdm import tqdm
from rdkit import Chem
import numpy as np

# 设置日志
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    filename='pocket_generation.log',
                    filemode='w')
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

def create_pymol_script(protein_file, ligand_file, output_file, cutoff=5.0):
    """创建PyMOL脚本来提取蛋白质口袋"""
    pymol_script = f"""
import pymol
from pymol import cmd

# 加载蛋白质和配体
cmd.load("{protein_file}", "protein")
cmd.load("{ligand_file}", "ligand")

# 选择距离配体{cutoff}埃内的残基作为口袋
cmd.select("pocket", f"protein within {cutoff} of ligand")

# 保存口袋
cmd.save("{output_file}", "pocket")

# 退出PyMOL
cmd.quit()
"""
    return pymol_script

def generate_pocket(protein_dir, ligand_dir, cutoff=5.0, parallel=True):
    """为配体生成口袋PDB文件"""
    # 获取蛋白质ID和配体ID
    protein_id = os.path.basename(protein_dir)
    ligand_id = os.path.basename(ligand_dir)
    
    # 定位文件
    protein_file = os.path.join(protein_dir, f"{protein_id}.pdb")
    ligand_file = os.path.join(ligand_dir, f"{ligand_id}.sdf")
    pocket_file = os.path.join(ligand_dir, f"{ligand_id}_pocket_{cutoff}A.pdb")
    
    # 检查输入文件
    if not os.path.exists(protein_file):
        logging.error(f"蛋白质文件不存在: {protein_file}")
        return False
        
    if not os.path.exists(ligand_file):
        logging.error(f"配体文件不存在: {ligand_file}")
        return False
    
    # 检查输出文件是否已存在
    if os.path.exists(pocket_file):
        logging.info(f"口袋文件已存在: {pocket_file}")
        return True
    
    # 创建临时PyMOL脚本
    script_file = os.path.join(ligand_dir, "temp_pymol.pml")
    with open(script_file, "w") as f:
        f.write(create_pymol_script(protein_file, ligand_file, pocket_file, cutoff))
    
    try:
        # 使用PyMOL运行脚本
        subprocess.run(["pymol", "-cq", script_file], check=True, 
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # 检查生成的文件
        if os.path.exists(pocket_file) and os.path.getsize(pocket_file) > 0:
            logging.info(f"成功生成口袋文件: {pocket_file}")
            return True
        else:
            logging.error(f"口袋文件生成失败或为空: {pocket_file}")
            return False
    except Exception as e:
        logging.error(f"生成口袋时发生错误: {str(e)}")
        return False
    finally:
        # 清理临时文件
        if os.path.exists(script_file):
            os.remove(script_file)

def process_protein(protein_dir, cutoff=5.0, parallel=True):
    """处理单个蛋白质目录中的所有配体"""
    protein_id = os.path.basename(protein_dir)
    logging.info(f"处理蛋白质: {protein_id}")
    
    # 获取蛋白质目录下的所有配体目录
    ligand_dirs = [d for d in glob.glob(os.path.join(protein_dir, "*")) 
                  if os.path.isdir(d) and os.path.basename(d).startswith(protein_id)]
    
    if not ligand_dirs:
        logging.warning(f"蛋白质目录中没有找到配体子目录: {protein_dir}")
        return 0
    
    success_count = 0
    if parallel:
        # 并行处理每个配体
        with multiprocessing.Pool(processes=min(os.cpu_count(), len(ligand_dirs))) as pool:
            results = list(tqdm(
                pool.starmap(
                    generate_pocket, 
                    [(protein_dir, ligand_dir, cutoff, False) for ligand_dir in ligand_dirs]
                ),
                total=len(ligand_dirs),
                desc=f"生成 {protein_id} 的口袋"
            ))
            success_count = sum(results)
    else:
        # 顺序处理每个配体
        for ligand_dir in tqdm(ligand_dirs, desc=f"生成 {protein_id} 的口袋"):
            if generate_pocket(protein_dir, ligand_dir, cutoff, False):
                success_count += 1
    
    logging.info(f"蛋白质 {protein_id} 处理完成: {success_count}/{len(ligand_dirs)} 个口袋成功生成")
    return success_count

def process_all_proteins(base_dir, cutoff=5.0, parallel=True, error_report=None, only_failed=None):
    """处理基础目录中的所有蛋白质"""
    logging.info(f"开始在 {base_dir} 中生成口袋")
    
    # 获取所有蛋白质目录
    protein_dirs = [d for d in glob.glob(os.path.join(base_dir, "*")) 
                   if os.path.isdir(d)]
    
    if not protein_dirs:
        logging.error(f"在 {base_dir} 中没有找到蛋白质目录")
        return 0
    
    total_success = 0
    total_proteins = len(protein_dirs)
    
    logging.info(f"找到 {total_proteins} 个蛋白质目录")
    
    for i, protein_dir in enumerate(protein_dirs):
        logging.info(f"处理蛋白质 {i+1}/{total_proteins}: {os.path.basename(protein_dir)}")
        if only_failed and os.path.basename(protein_dir) in only_failed:
            logging.info(f"跳过失败样本: {os.path.basename(protein_dir)}")
            continue
        success_count = process_protein(protein_dir, cutoff, parallel)
        total_success += success_count
    
    logging.info(f"所有蛋白质处理完成。成功生成口袋: {total_success}")
    return total_success

def main():
    parser = argparse.ArgumentParser(description="根据蛋白质和配体生成口袋PDB文件")
    parser.add_argument("--base_dir", type=str, default="train_pdb", 
                        help="包含蛋白质和配体数据的基础目录")
    parser.add_argument("--cutoff", type=float, default=5.0, 
                        help="定义口袋的距离阈值(埃)")
    parser.add_argument("--no_parallel", action="store_true", 
                        help="禁用并行处理")
    parser.add_argument("--error_report", type=str, default="", 
                        help="记录错误样本的报告文件")
    parser.add_argument("--only_process_failed", type=str, default="", 
                        help="仅处理指定错误报告文件中的失败样本")
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    # 读取错误报告文件中的失败样本列表
    failed_samples = set()
    if args.only_process_failed and os.path.exists(args.only_process_failed):
        logging.info(f"读取失败样本列表: {args.only_process_failed}")
        with open(args.only_process_failed, 'r') as f:
            for line in f:
                if line.strip() and ':' in line:
                    sample_path = line.split(':', 1)[0].strip()
                    failed_samples.add(sample_path)
        logging.info(f"加载了 {len(failed_samples)} 个失败样本待处理")
    
    # 打开错误报告文件
    error_report = None
    if args.error_report:
        try:
            error_report = open(args.error_report, 'w')
            logging.info(f"将记录错误样本到文件: {args.error_report}")
        except Exception as e:
            logging.error(f"无法创建错误报告文件: {args.error_report}, 错误: {str(e)}")
    
    # 处理所有蛋白质
    success_count = process_all_proteins(
        args.base_dir, 
        cutoff=args.cutoff, 
        parallel=not args.no_parallel,
        error_report=error_report,
        only_failed=failed_samples if args.only_process_failed else None
    )
    
    # 关闭错误报告文件
    if error_report:
        error_report.close()
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    logging.info(f"任务完成! 耗时: {int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}")
    print(f"\n{'='*80}")
    print(f"口袋生成完成!")
    print(f"成功生成: {success_count} 个口袋文件")
    print(f"总处理时间: {int(hours):02d}小时{int(minutes):02d}分钟{int(seconds):02d}秒")
    print(f"{'='*80}")

if __name__ == "__main__":
    main() 