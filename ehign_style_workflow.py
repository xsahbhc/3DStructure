#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import time
import subprocess
import sys
from datetime import datetime

def run_command(command, description=None):
    """运行命令并捕获输出"""
    if description:
        print(f"\n{'='*80}")
        print(f"运行: {description}")
        print(f"{'='*80}")
        
    start_time = time.time()
    try:
        subprocess.run(command, check=True)
        success = True
    except subprocess.CalledProcessError as e:
        print(f"错误: 命令执行失败，代码: {e.returncode}")
        success = False
        
    elapsed_time = time.time() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    
    if description:
        status = "成功" if success else "失败"
        print(f"{description} {status}完成，耗时: {int(minutes)}分{int(seconds)}秒")
        
    return success

def main():
    parser = argparse.ArgumentParser(description="使用EHIGN风格的数据处理流程")
    parser.add_argument("--train_pdb", default="train_pdb", help="包含蛋白质-配体数据的根目录")
    parser.add_argument("--cutoff", type=float, default=5.0, help="口袋距离阈值(埃)")
    parser.add_argument("--dis_threshold", type=float, default=5.0, help="配体-口袋原子间距阈值(埃)")
    parser.add_argument("--processes", type=int, default=4, help="并行进程数")
    parser.add_argument("--skip_pocket", action="store_true", help="跳过口袋生成步骤")
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    # 步骤1: 生成口袋 (使用更简单的PyMOL方法)
    if not args.skip_pocket:
        pocket_cmd = [
            sys.executable, "generate_pocket_simple.py",
            "--base_dir", args.train_pdb,
            "--cutoff", str(args.cutoff)
        ]
        
        if not run_command(pocket_cmd, "口袋生成"):
            print("口袋生成失败，流程终止")
            return 1
    else:
        print("跳过口袋生成步骤")
        
    # 步骤2: 生成EHIGN格式数据 (使用EHIGN风格的代码)
    ehign_cmd = [
        sys.executable, "generate_ehign_simple.py",
        "--base_dir", args.train_pdb,
        "--cutoff", str(args.cutoff),
        "--dis_threshold", str(args.dis_threshold),
        "--processes", str(args.processes)
    ]
    
    if not run_command(ehign_cmd, "EHIGN格式生成"):
        print("EHIGN格式生成失败，流程终止")
        return 1
        
    # 计算总处理时间
    elapsed_time = time.time() - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    print(f"\n{'='*80}")
    print(f"EHIGN风格处理流程完成!")
    print(f"总处理时间: {int(hours):02d}小时{int(minutes):02d}分钟{int(seconds):02d}秒")
    print(f"{'='*80}")
    
    return 0
    
if __name__ == "__main__":
    sys.exit(main())
