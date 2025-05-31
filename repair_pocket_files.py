#!/usr/bin/env python
import os
import glob
import argparse
from tqdm import tqdm

def repair_pdb_file(pdb_file):
    """修复PDB文件格式问题"""
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
            
        # 检查是否为空
        if not lines:
            print(f"文件为空: {pdb_file}")
            return False
            
        # 尝试修复常见的PDB格式问题
        fixed_lines = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # 确保行长度至少为54个字符（基本的坐标信息）
                if len(line) < 54:
                    continue
                    
                # 尝试修复原子名称和残基名称之间的间距
                if len(line) >= 27:
                    atom_name = line[12:16].strip()
                    rest_of_line = line[16:]
                    
                    # 重建ATOM行，确保正确的格式
                    new_line = f"ATOM  {line[6:11]} {atom_name:4s}{rest_of_line}"
                    fixed_lines.append(new_line)
                else:
                    fixed_lines.append(line)
            else:
                fixed_lines.append(line)
                
        # 确保文件以END结束
        if not any(line.startswith('END') for line in fixed_lines):
            fixed_lines.append('END\n')
            
        # 写回修复的文件
        with open(pdb_file, 'w') as f:
            f.writelines(fixed_lines)
            
        return True
        
    except Exception as e:
        print(f"修复文件 {pdb_file} 时出错: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="修复无法加载的口袋PDB文件")
    parser.add_argument("--base_dir", default="train_pdb", help="包含蛋白质-配体数据的根目录")
    parser.add_argument("--error_report", default="", help="错误报告文件，只修复报告中的文件")
    
    args = parser.parse_args()
    
    # 找到所有口袋文件
    pocket_files = []
    
    if args.error_report and os.path.exists(args.error_report):
        # 从错误报告中提取需要修复的蛋白质目录
        protein_dirs = set()
        with open(args.error_report, 'r') as f:
            for line in f:
                if ":" in line:
                    protein_dir = line.split(":", 1)[0].strip()
                    protein_dirs.add(protein_dir)
                    
        # 找出这些目录中的所有口袋文件
        for protein_dir in protein_dirs:
            if os.path.exists(protein_dir):
                for root, dirs, files in os.walk(protein_dir):
                    for file in files:
                        if "_pocket_" in file and file.endswith(".pdb"):
                            pocket_files.append(os.path.join(root, file))
    else:
        # 搜索所有口袋文件
        for root, dirs, files in os.walk(args.base_dir):
            for file in files:
                if "_pocket_" in file and file.endswith(".pdb"):
                    pocket_files.append(os.path.join(root, file))
    
    print(f"找到 {len(pocket_files)} 个口袋PDB文件")
    
    # 修复所有文件
    success_count = 0
    for pdb_file in tqdm(pocket_files, desc="修复口袋文件"):
        if repair_pdb_file(pdb_file):
            success_count += 1
    
    print(f"成功修复 {success_count}/{len(pocket_files)} 个口袋文件")

if __name__ == "__main__":
    main()
