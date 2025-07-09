#!/usr/bin/env python3
"""
检查BindingDB和BioSNAP数据集的1D/2D和3D数据一致性
"""

import pandas as pd
import os
from pathlib import Path

def check_bindingdb_consistency():
    """检查BindingDB的1D/2D和3D数据一致性"""
    print("=" * 60)
    print("检查BindingDB数据一致性")
    print("=" * 60)
    
    # 1D/2D数据路径
    data_1d2d_base = Path("../drugban/datasets/bindingdb_3d_sequences/random")
    # 3D数据路径
    data_3d_base = Path("bindingdb/train_csv")
    
    splits = ['train_stratified', 'val_stratified', 'test_stratified']
    
    for split in splits:
        print(f"\n📊 检查 {split} 数据:")
        
        # 读取1D/2D数据
        file_1d2d = data_1d2d_base / f"{split}.csv"
        file_3d = data_3d_base / f"{split}.csv"
        
        if not file_1d2d.exists():
            print(f"  ❌ 1D/2D文件不存在: {file_1d2d}")
            continue
            
        if not file_3d.exists():
            print(f"  ❌ 3D文件不存在: {file_3d}")
            continue
            
        df_1d2d = pd.read_csv(file_1d2d)
        df_3d = pd.read_csv(file_3d)
        
        print(f"  📈 1D/2D数据: {len(df_1d2d)} 条记录")
        print(f"  📈 3D数据: {len(df_3d)} 条记录")
        
        # 检查标签分布
        if 'Y' in df_1d2d.columns:
            label_dist_1d2d = df_1d2d['Y'].value_counts().sort_index()
            print(f"  🏷️  1D/2D标签分布: {label_dist_1d2d.to_dict()}")
        
        if 'label' in df_3d.columns:
            label_dist_3d = df_3d['label'].value_counts().sort_index()
            print(f"  🏷️  3D标签分布: {label_dist_3d.to_dict()}")
        
        # 检查seqid一致性（如果可能）
        if 'seqid' in df_3d.columns:
            # 读取seqid映射
            seqid_mapping_file = data_1d2d_base / "seqid_mapping.csv"
            if seqid_mapping_file.exists():
                seqid_mapping = pd.read_csv(seqid_mapping_file)
                
                # 从3D数据中提取seqid
                seqids_3d = set(df_3d['seqid'].astype(str))
                seqids_1d2d = set(seqid_mapping['seqid'].astype(str))
                
                # 检查当前split的seqid是否在映射中
                if split == 'train_stratified':
                    # 对于训练集，检查1D/2D数据的seqid是否都在3D数据中
                    missing_in_3d = seqids_1d2d - seqids_3d
                    missing_in_1d2d = seqids_3d - seqids_1d2d
                    
                    print(f"  🔗 seqid映射检查:")
                    print(f"    - 1D/2D中有但3D中没有的seqid: {len(missing_in_3d)}")
                    print(f"    - 3D中有但1D/2D中没有的seqid: {len(missing_in_1d2d)}")
                    
                    if missing_in_3d:
                        print(f"    - 示例缺失seqid (1D/2D->3D): {list(missing_in_3d)[:5]}")
                    if missing_in_1d2d:
                        print(f"    - 示例缺失seqid (3D->1D/2D): {list(missing_in_1d2d)[:5]}")

def check_biosnap_consistency():
    """检查BioSNAP的1D/2D和3D数据一致性"""
    print("\n" + "=" * 60)
    print("检查BioSNAP数据一致性")
    print("=" * 60)
    
    # 1D/2D数据路径
    data_1d2d_base = Path("biosnap/biosnap_3d_sequences/random")
    # 3D数据路径
    data_3d_base = Path("biosnap/train_csv")
    
    splits = ['train_stratified', 'val_stratified', 'test_stratified']
    
    for split in splits:
        print(f"\n📊 检查 {split} 数据:")
        
        # 读取1D/2D数据
        file_1d2d = data_1d2d_base / f"{split}.csv"
        file_3d = data_3d_base / f"{split}.csv"
        
        if not file_1d2d.exists():
            print(f"  ❌ 1D/2D文件不存在: {file_1d2d}")
            continue
            
        if not file_3d.exists():
            print(f"  ❌ 3D文件不存在: {file_3d}")
            continue
            
        df_1d2d = pd.read_csv(file_1d2d)
        df_3d = pd.read_csv(file_3d)
        
        print(f"  📈 1D/2D数据: {len(df_1d2d)} 条记录")
        print(f"  📈 3D数据: {len(df_3d)} 条记录")
        
        # 检查标签分布
        if 'Y' in df_1d2d.columns:
            label_dist_1d2d = df_1d2d['Y'].value_counts().sort_index()
            print(f"  🏷️  1D/2D标签分布: {label_dist_1d2d.to_dict()}")
        
        if 'label' in df_3d.columns:
            label_dist_3d = df_3d['label'].value_counts().sort_index()
            print(f"  🏷️  3D标签分布: {label_dist_3d.to_dict()}")

def check_seqid_mapping_consistency():
    """检查seqid映射的一致性"""
    print("\n" + "=" * 60)
    print("检查seqid映射一致性")
    print("=" * 60)
    
    # BindingDB
    print("\n🔍 BindingDB seqid映射检查:")
    bindingdb_mapping = Path("../drugban/datasets/bindingdb_3d_sequences/random/seqid_mapping.csv")
    bindingdb_labels = Path("bindingdb/train_csv/labels.csv")
    
    if bindingdb_mapping.exists() and bindingdb_labels.exists():
        df_mapping = pd.read_csv(bindingdb_mapping)
        df_labels = pd.read_csv(bindingdb_labels)
        
        print(f"  📊 seqid映射文件: {len(df_mapping)} 条记录")
        print(f"  📊 3D标签文件: {len(df_labels)} 条记录")
        
        # 检查complex_id和seqid的对应关系
        if 'complex_id' in df_labels.columns:
            # 从complex_id中提取seqid
            seqids_from_complex = []
            for complex_id in df_labels['complex_id']:
                parts = complex_id.split('_')
                if len(parts) >= 3:
                    seqids_from_complex.append(int(parts[-1]))
            
            seqids_from_complex = set(seqids_from_complex)
            seqids_from_mapping = set(df_mapping['seqid'])
            
            print(f"  🔗 从complex_id提取的seqid数量: {len(seqids_from_complex)}")
            print(f"  🔗 映射文件中的seqid数量: {len(seqids_from_mapping)}")
            
            missing_in_mapping = seqids_from_complex - seqids_from_mapping
            missing_in_complex = seqids_from_mapping - seqids_from_complex
            
            print(f"  ❓ complex_id中有但映射中没有: {len(missing_in_mapping)}")
            print(f"  ❓ 映射中有但complex_id中没有: {len(missing_in_complex)}")
    
    # BioSNAP
    print("\n🔍 BioSNAP seqid映射检查:")
    biosnap_mapping = Path("biosnap/biosnap_3d_sequences/random/seqid_mapping.csv")
    biosnap_labels = Path("biosnap/train_csv/labels.csv")
    
    if biosnap_mapping.exists() and biosnap_labels.exists():
        df_mapping = pd.read_csv(biosnap_mapping)
        df_labels = pd.read_csv(biosnap_labels)
        
        print(f"  📊 seqid映射文件: {len(df_mapping)} 条记录")
        print(f"  📊 3D标签文件: {len(df_labels)} 条记录")
        
        # 检查complex_id和seqid的对应关系
        if 'complex_id' in df_labels.columns:
            # 从complex_id中提取seqid
            seqids_from_complex = []
            for complex_id in df_labels['complex_id']:
                parts = complex_id.split('_')
                if len(parts) >= 3:
                    seqids_from_complex.append(int(parts[-1]))
            
            seqids_from_complex = set(seqids_from_complex)
            seqids_from_mapping = set(df_mapping['seqid'])
            
            print(f"  🔗 从complex_id提取的seqid数量: {len(seqids_from_complex)}")
            print(f"  🔗 映射文件中的seqid数量: {len(seqids_from_mapping)}")
            
            missing_in_mapping = seqids_from_complex - seqids_from_mapping
            missing_in_complex = seqids_from_mapping - seqids_from_complex
            
            print(f"  ❓ complex_id中有但映射中没有: {len(missing_in_mapping)}")
            print(f"  ❓ 映射中有但complex_id中没有: {len(missing_in_complex)}")

if __name__ == "__main__":
    try:
        check_bindingdb_consistency()
        check_biosnap_consistency()
        check_seqid_mapping_consistency()
        
        print("\n" + "=" * 60)
        print("✅ 数据一致性检查完成")
        print("=" * 60)
        
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()
