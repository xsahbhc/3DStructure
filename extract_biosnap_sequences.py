#!/usr/bin/env python3
"""
为BioSNAP数据集提取成功转换为3D数据的样本对应的原始CSV序列数据
生成与BindingDB相同格式的1D/2D数据
"""

import pandas as pd
import os
import re
from pathlib import Path
from sklearn.model_selection import train_test_split

def extract_seqids_from_3d_data():
    """从BioSNAP 3D数据目录中提取所有成功转换的seqid"""
    train_pdb_dir = Path("biosnap/train_pdb")
    seqids = set()
    
    print("正在扫描BioSNAP 3D数据目录...")
    
    if not train_pdb_dir.exists():
        print(f"错误: 3D数据目录不存在: {train_pdb_dir}")
        return []
    
    # 遍历所有蛋白质目录
    for protein_dir in train_pdb_dir.iterdir():
        if protein_dir.is_dir():
            # 遍历每个蛋白质目录下的复合物文件
            for item in protein_dir.iterdir():
                if item.is_dir():
                    # 提取seqid (格式如: 1ACB_E_10127)
                    match = re.search(r'_(\d+)$', item.name)
                    if match:
                        seqid = int(match.group(1))
                        seqids.add(seqid)
    
    print(f"找到 {len(seqids)} 个成功转换的seqid")
    return sorted(list(seqids))

def create_biosnap_drugban_dataset():
    """创建BioSNAP的DrugBAN格式数据集"""
    
    # 1. 提取成功转换的seqid
    successful_seqids = extract_seqids_from_3d_data()
    
    if not successful_seqids:
        print("错误: 没有找到成功转换的3D数据")
        return None
    
    # 2. 读取原始train.csv
    print("读取原始BioSNAP train.csv...")
    train_csv_path = "biosnap/train_csv/train.csv"
    
    if not os.path.exists(train_csv_path):
        print(f"错误: 找不到原始train.csv文件: {train_csv_path}")
        return None
        
    train_df = pd.read_csv(train_csv_path)
    print(f"原始train.csv有 {len(train_df)} 条记录")
    
    # 3. 筛选成功转换的数据
    filtered_df = train_df[train_df['seqid'].isin(successful_seqids)].copy()
    print(f"筛选后有 {len(filtered_df)} 条记录")
    
    if len(filtered_df) == 0:
        print("错误: 筛选后没有数据")
        return None
    
    # 4. 转换为DrugBAN格式
    drugban_df = pd.DataFrame({
        'SMILES': filtered_df['SMILES'],
        'Protein': filtered_df['Protein'],
        'Y': filtered_df['Y']
    })
    
    # 5. 创建输出目录
    output_base = Path("3D_structure/biosnap_3d_sequences")
    output_dir = output_base / "random"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"创建输出目录: {output_dir}")
    
    # 6. 进行分层划分 (与BindingDB相同的比例)
    print("进行分层划分...")
    
    # 检查标签分布
    label_counts = drugban_df['Y'].value_counts()
    print(f"标签分布: {label_counts.to_dict()}")
    
    # 确保两个类别都有足够的样本
    min_class_count = label_counts.min()
    if min_class_count < 3:
        print(f"警告: 最小类别样本数过少 ({min_class_count})，可能影响分层划分")
    
    # 第一次划分: 训练集 vs (验证集+测试集)
    train_df, temp_df = train_test_split(
        drugban_df, 
        test_size=0.2,  # 20%用于验证和测试
        stratify=drugban_df['Y'],
        random_state=42
    )
    
    # 第二次划分: 验证集 vs 测试集
    val_df, test_df = train_test_split(
        temp_df,
        test_size=0.5,  # 一半用于验证，一半用于测试
        stratify=temp_df['Y'],
        random_state=42
    )
    
    print(f"数据划分结果:")
    print(f"  训练集: {len(train_df)} 条记录")
    print(f"  验证集: {len(val_df)} 条记录")
    print(f"  测试集: {len(test_df)} 条记录")
    
    # 7. 保存划分后的数据
    train_df.to_csv(output_dir / "train_stratified.csv", index=False)
    val_df.to_csv(output_dir / "val_stratified.csv", index=False)
    test_df.to_csv(output_dir / "test_stratified.csv", index=False)
    
    print(f"数据文件已保存到 {output_dir}/")
    
    # 8. 保存seqid映射信息
    seqid_mapping = filtered_df[['seqid', 'SMILES', 'Protein', 'Y']].copy()
    seqid_mapping.to_csv(output_dir / "seqid_mapping.csv", index=False)
    print(f"seqid映射信息已保存到 {output_dir}/seqid_mapping.csv")
    
    # 9. 验证数据完整性
    print("\n验证数据完整性...")
    total_generated = len(train_df) + len(val_df) + len(test_df)
    print(f"生成的总样本数: {total_generated}")
    print(f"原始筛选样本数: {len(filtered_df)}")
    
    if total_generated == len(filtered_df):
        print("✅ 数据完整性验证通过")
    else:
        print("❌ 数据完整性验证失败")
    
    # 10. 显示各集合的标签分布
    print("\n各数据集标签分布:")
    for name, df in [("训练集", train_df), ("验证集", val_df), ("测试集", test_df)]:
        counts = df['Y'].value_counts().sort_index()
        print(f"  {name}: {counts.to_dict()}")
    
    return output_dir

def move_to_3d_structure():
    """将生成的数据移动到3D_structure目录下"""
    source_dir = Path("3D_structure/biosnap_3d_sequences")
    target_dir = Path("3D_structure/biosnap/biosnap_3d_sequences")
    
    if source_dir.exists():
        print(f"移动数据从 {source_dir} 到 {target_dir}")
        
        # 创建目标目录
        target_dir.parent.mkdir(parents=True, exist_ok=True)
        
        # 移动目录
        import shutil
        if target_dir.exists():
            shutil.rmtree(target_dir)
        shutil.move(str(source_dir), str(target_dir))
        
        print(f"✅ 数据已移动到: {target_dir}")
        return target_dir
    else:
        print(f"错误: 源目录不存在: {source_dir}")
        return None

if __name__ == "__main__":
    try:
        print("🚀 开始为BioSNAP生成1D/2D序列数据...")
        
        # 生成数据
        output_dir = create_biosnap_drugban_dataset()
        
        if output_dir:
            print(f"\n✅ BioSNAP 1D/2D数据集创建成功!")
            print(f"📁 数据保存在: {output_dir}")
            
            # 移动到3D_structure目录
            final_dir = move_to_3d_structure()
            if final_dir:
                print(f"📁 最终位置: {final_dir}")
        else:
            print("❌ 数据集创建失败")
            
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()
