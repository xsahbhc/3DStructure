#!/usr/bin/env python3
"""
清理BioSNAP数据集中不完整的配体文件夹
只保留完整可用的配体对
"""

import os
import shutil
from pathlib import Path
import csv

def check_pocket_file(pocket_file):
    """检查口袋文件是否有效"""
    if not os.path.exists(pocket_file):
        return False, "口袋文件不存在"
    
    size = os.path.getsize(pocket_file)
    if size < 100:
        return False, f"口袋文件太小 ({size} bytes)"
    
    try:
        with open(pocket_file, 'r') as f:
            content = f.read()
            if 'ATOM' not in content:
                return False, "口袋文件无ATOM记录"
    except:
        return False, "口袋文件读取失败"
    
    return True, "有效"

def is_complete_ligand_folder(ligand_dir):
    """检查配体文件夹是否完整"""
    ligand_id = ligand_dir.name
    
    # 检查必需文件
    required_files = {
        'sdf': ligand_dir / f"{ligand_id}.sdf",
        'pocket': ligand_dir / f"{ligand_id}_pocket_5.0A.pdb",
        'ligand_mol2': ligand_dir / f"{ligand_id}_ligand.mol2",
        'ligand_pdb': ligand_dir / f"{ligand_id}_ligand.pdb",
        'rdkit': ligand_dir / f"{ligand_id}.rdkit",
        'dgl': ligand_dir / f"{ligand_id}.dgl"
    }
    
    # 检查每个文件
    for file_type, file_path in required_files.items():
        if file_type == 'pocket':
            valid, _ = check_pocket_file(file_path)
            if not valid:
                return False
        else:
            if not file_path.exists() or file_path.stat().st_size < 100:
                return False
    
    return True

def cleanup_biosnap_data(dry_run=True):
    """清理BioSNAP数据，删除不完整的配体文件夹"""
    base_dir = Path("biosnap/train_pdb")
    
    if not base_dir.exists():
        print(f"❌ 错误: 目录 {base_dir} 不存在")
        return
    
    print("🧹 开始清理BioSNAP数据集...")
    print("=" * 80)
    
    if dry_run:
        print("🔍 DRY RUN模式 - 只显示将要删除的文件夹，不实际删除")
    else:
        print("⚠️  实际删除模式 - 将永久删除不完整的文件夹")
    
    print("=" * 80)
    
    # 统计变量
    total_ligand_folders = 0
    complete_folders = 0
    incomplete_folders = 0
    deleted_folders = []
    kept_folders = []
    
    # 遍历所有蛋白质文件夹
    for protein_dir in base_dir.iterdir():
        if not protein_dir.is_dir():
            continue
            
        protein_id = protein_dir.name
        protein_pdb = protein_dir / f"{protein_id}.pdb"
        
        # 检查蛋白质PDB文件
        if not protein_pdb.exists():
            print(f"⚠️  跳过 {protein_id}: 缺少蛋白质PDB文件")
            continue
        
        # 收集该蛋白质下的所有配体文件夹
        ligand_folders_to_delete = []
        ligand_folders_to_keep = []
        
        for ligand_dir in protein_dir.iterdir():
            if not ligand_dir.is_dir():
                continue
                
            ligand_id = ligand_dir.name
            
            # 跳过不符合命名规范的文件夹
            if not ligand_id.startswith(protein_id):
                continue
            
            total_ligand_folders += 1
            
            # 检查是否完整
            if is_complete_ligand_folder(ligand_dir):
                complete_folders += 1
                ligand_folders_to_keep.append(ligand_dir)
                kept_folders.append(str(ligand_dir))
            else:
                incomplete_folders += 1
                ligand_folders_to_delete.append(ligand_dir)
                deleted_folders.append(str(ligand_dir))
        
        # 删除不完整的配体文件夹
        for ligand_dir in ligand_folders_to_delete:
            if dry_run:
                print(f"🗑️  将删除: {ligand_dir}")
            else:
                try:
                    shutil.rmtree(ligand_dir)
                    print(f"✅ 已删除: {ligand_dir}")
                except Exception as e:
                    print(f"❌ 删除失败 {ligand_dir}: {e}")
        
        # 如果蛋白质文件夹下没有完整的配体了，也删除蛋白质文件夹
        if not ligand_folders_to_keep:
            if dry_run:
                print(f"🗑️  将删除空蛋白质文件夹: {protein_dir}")
            else:
                try:
                    shutil.rmtree(protein_dir)
                    print(f"✅ 已删除空蛋白质文件夹: {protein_dir}")
                except Exception as e:
                    print(f"❌ 删除蛋白质文件夹失败 {protein_dir}: {e}")
    
    # 输出统计结果
    print(f"\n📊 清理统计结果:")
    print(f"=" * 80)
    print(f"📊 总配体文件夹数量: {total_ligand_folders}")
    print(f"✅ 保留的完整配体文件夹: {complete_folders}")
    print(f"🗑️  {'将删除' if dry_run else '已删除'}的不完整文件夹: {incomplete_folders}")
    print(f"📈 保留率: {(complete_folders/total_ligand_folders*100):.2f}%")
    
    # 保存清理记录
    if not dry_run:
        # 保存删除记录
        with open('biosnap/deleted_folders.txt', 'w') as f:
            f.write("删除的不完整配体文件夹:\n")
            f.write("=" * 50 + "\n")
            for folder in deleted_folders:
                f.write(f"{folder}\n")
        
        # 保存保留记录
        with open('biosnap/kept_folders.txt', 'w') as f:
            f.write("保留的完整配体文件夹:\n")
            f.write("=" * 50 + "\n")
            for folder in kept_folders:
                f.write(f"{folder}\n")
        
        print(f"\n📝 清理记录已保存:")
        print(f"   - 删除记录: biosnap/deleted_folders.txt")
        print(f"   - 保留记录: biosnap/kept_folders.txt")
    
    print(f"\n🎉 清理{'预览' if dry_run else ''}完成!")
    print(f"✅ 最终可用于训练的配体对: {complete_folders}")
    
    return {
        'total': total_ligand_folders,
        'complete': complete_folders,
        'deleted': incomplete_folders,
        'dry_run': dry_run
    }

def main():
    print("🧹 BioSNAP数据清理工具")
    print("=" * 80)
    
    # 首先运行dry run
    print("第一步: 预览将要删除的文件夹")
    result = cleanup_biosnap_data(dry_run=True)
    
    if result['deleted'] == 0:
        print("✅ 所有配体文件夹都是完整的，无需清理！")
        return
    
    print(f"\n⚠️  将删除 {result['deleted']} 个不完整的配体文件夹")
    print(f"✅ 将保留 {result['complete']} 个完整的配体文件夹")
    
    # 确认是否执行实际删除
    print("\n" + "=" * 80)
    response = input("确认执行实际删除操作吗？(输入 'YES' 确认): ")
    
    if response == 'YES':
        print("\n🚀 开始实际清理...")
        cleanup_biosnap_data(dry_run=False)
    else:
        print("❌ 取消清理操作")

if __name__ == "__main__":
    main()
