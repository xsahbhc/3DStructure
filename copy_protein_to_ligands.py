#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import shutil
from tqdm import tqdm

# 蛋白质数据根目录
base_dir = "3D_structure/train_pdb"

# 找到所有蛋白质目录
protein_dirs = [d for d in glob.glob(os.path.join(base_dir, "*")) if os.path.isdir(d)]
print(f"找到 {len(protein_dirs)} 个蛋白质目录")

total_copied = 0
total_failed = 0

# 遍历每个蛋白质目录
for protein_dir in tqdm(protein_dirs, desc="处理蛋白质"):
    protein_id = os.path.basename(protein_dir)
    protein_file = os.path.join(protein_dir, f"{protein_id}.pdb")
    
    # 检查蛋白质PDB文件是否存在
    if not os.path.exists(protein_file):
        print(f"警告: 未找到蛋白质文件 {protein_file}")
        continue
        
    # 找到该蛋白质下的所有配体目录
    ligand_dirs = [d for d in glob.glob(os.path.join(protein_dir, "*")) 
                  if os.path.isdir(d) and os.path.basename(d).startswith(protein_id)]
    
    if not ligand_dirs:
        print(f"警告: 蛋白质 {protein_id} 没有关联的配体目录")
        continue
    
    # 复制蛋白质文件到每个配体目录
    for ligand_dir in ligand_dirs:
        target_file = os.path.join(ligand_dir, f"{protein_id}.pdb")
        try:
            if not os.path.exists(target_file):
                shutil.copy2(protein_file, target_file)
                total_copied += 1
        except Exception as e:
            print(f"复制到 {ligand_dir} 失败: {str(e)}")
            total_failed += 1

print(f"完成! 成功复制 {total_copied} 个蛋白质文件, {total_failed} 个复制失败") 