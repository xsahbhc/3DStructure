#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import glob
import subprocess
import multiprocessing
from tqdm import tqdm
import pymol
from pymol import cmd

def generate_pocket_pymol(protein_dir, ligand_dir, cutoff=5.0):
    """使用PyMOL直接生成口袋，类似EHIGN方法"""
    # 获取蛋白质ID和配体ID
    protein_id = os.path.basename(protein_dir)
    ligand_id = os.path.basename(ligand_dir)
    
    # 定位文件
    protein_file = os.path.join(protein_dir, f"{protein_id}.pdb")
    ligand_file = os.path.join(ligand_dir, f"{ligand_id}.sdf")
    pocket_file = os.path.join(ligand_dir, f"{ligand_id}_pocket_{cutoff}A.pdb")
    
    # 检查输入文件
    if not os.path.exists(protein_file) or not os.path.exists(ligand_file):
        print(f"蛋白质或配体文件不存在: {protein_file}, {ligand_file}")
        return False
    
    # 检查输出文件是否已存在
    if os.path.exists(pocket_file):
        print(f"口袋文件已存在: {pocket_file}")
        return True
    
    # 使用PyMOL提取口袋
    try:
        cmd.reinitialize()
        cmd.load(protein_file, "protein")
        cmd.load(ligand_file, "ligand")
        cmd.remove('resn HOH')  # 去除水分子
        cmd.select('Pocket', f'byres protein within {cutoff} of ligand')
        cmd.save(pocket_file, 'Pocket')
        cmd.delete('all')
        
        if os.path.exists(pocket_file) and os.path.getsize(pocket_file) > 0:
            print(f"成功生成口袋文件: {pocket_file}")
            return True
        else:
            print(f"口袋文件生成失败或为空: {pocket_file}")
            return False
    except Exception as e:
        print(f"生成口袋时发生错误: {str(e)}")
        return False

def process_all_proteins_simple(base_dir, cutoff=5.0):
    """处理基础目录中的所有蛋白质"""
    # 获取所有蛋白质目录
    protein_dirs = [d for d in glob.glob(os.path.join(base_dir, "*")) 
                   if os.path.isdir(d)]
    
    total_success = 0
    total_ligands = 0
    for protein_dir in tqdm(protein_dirs, desc="处理蛋白质"):
        protein_id = os.path.basename(protein_dir)
        # 获取蛋白质目录下的所有配体目录
        ligand_dirs = [d for d in glob.glob(os.path.join(protein_dir, "*")) 
                      if os.path.isdir(d) and os.path.basename(d).startswith(protein_id)]
        
        for ligand_dir in ligand_dirs:
            total_ligands += 1
            if generate_pocket_pymol(protein_dir, ligand_dir, cutoff):
                total_success += 1
    
    print(f"所有蛋白质处理完成。成功生成口袋: {total_success}/{total_ligands}")
    return total_success

def main():
    parser = argparse.ArgumentParser(description="直接使用PyMOL生成口袋PDB文件")
    parser.add_argument("--base_dir", type=str, default="train_pdb", 
                        help="包含蛋白质和配体数据的基础目录")
    parser.add_argument("--cutoff", type=float, default=5.0, 
                        help="定义口袋的距离阈值(埃)")
    
    args = parser.parse_args()
    pymol.finish_launching(['pymol', '-qc'])  # 启动PyMOL
    
    success_count = process_all_proteins_simple(args.base_dir, args.cutoff)
    
    print(f"\n{'='*80}")
    print(f"口袋生成完成!")
    print(f"成功生成: {success_count} 个口袋文件")
    print(f"{'='*80}")
    
if __name__ == "__main__":
    main()
