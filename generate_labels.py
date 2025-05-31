import pandas as pd
import os
import sys

# 添加调试输出
print("开始执行标签生成脚本...")

# 路径定义
train_csv_path = '/home/work/workspace/shi_shaoqun/snap/3D_structure/train.csv'
root_dir = '/home/work/workspace/shi_shaoqun/snap/3D_structure/train_pdb'
output_path = '/home/work/workspace/shi_shaoqun/snap/3D_structure/labels.csv'

# 检查文件和目录是否存在
print(f"检查CSV文件: {train_csv_path}")
if not os.path.exists(train_csv_path):
    print(f"错误: 找不到CSV文件 {train_csv_path}")
    sys.exit(1)

print(f"检查数据目录: {root_dir}")
if not os.path.exists(root_dir):
    print(f"错误: 找不到数据目录 {root_dir}")
    sys.exit(1)

# 读取现有CSV文件
try:
    print("读取train.csv...")
    train_df = pd.read_csv(train_csv_path)
    print(f"CSV文件读取成功，共有{len(train_df)}行")
    print(f"CSV文件列名: {train_df.columns.tolist()}")
    
    # 检查必需的列
    required_columns = ['seqid', 'Y']
    for col in required_columns:
        if col not in train_df.columns:
            print(f"错误: CSV文件中缺少必需的列 '{col}'")
            sys.exit(1)
    
except Exception as e:
    print(f"错误: 读取CSV文件时出错: {e}")
    sys.exit(1)

# 创建seqid到标签的映射
print("创建seqid到标签的映射...")
seqid_to_label = dict(zip(train_df['seqid'].astype(str), train_df['Y']))
print(f"映射创建成功，包含{len(seqid_to_label)}个seqid")

# 创建labels.csv的数据
labels_data = []
found_count = 0
not_found_count = 0

print(f"开始遍历数据目录: {root_dir}")
try:
    # 检查数据目录是否有内容
    protein_folders = os.listdir(root_dir)
    print(f"在根目录中找到{len(protein_folders)}个文件/文件夹")
    
    # 遍历所有目录，找出所有复合物ID
    for protein_folder in protein_folders:
        protein_path = os.path.join(root_dir, protein_folder)
        if os.path.isdir(protein_path):
            print(f"处理蛋白质文件夹: {protein_folder}")
            complex_folders = os.listdir(protein_path)
            print(f"  在{protein_folder}中找到{len(complex_folders)}个文件/文件夹")
            
            for complex_folder in complex_folders:
                complex_path = os.path.join(protein_path, complex_folder)
                if os.path.isdir(complex_path):
                    # 找到了一个复合物目录
                    print(f"  处理复合物: {complex_folder}")
                    
                    # 从complex_id中提取seqid (最后一部分)
                    parts = complex_folder.split('_')
                    if len(parts) >= 3:  # 确保格式为 PDB_Chain_SeqID
                        seqid = parts[-1]  # 取最后一部分作为seqid
                        print(f"    提取的seqid: {seqid}")
                        
                        if seqid in seqid_to_label:
                            label = seqid_to_label[seqid]
                            labels_data.append({
                                'complex_id': complex_folder,
                                'label': label
                            })
                            print(f"    ✓ 成功匹配到标签: {label}")
                            found_count += 1
                        else:
                            print(f"    ✗ 警告: 在train.csv中找不到seqid={seqid}的对应标签")
                            not_found_count += 1
                    else:
                        print(f"    ✗ 警告: 复合物名称格式不正确: {complex_folder}")
except Exception as e:
    print(f"错误: 遍历目录时出错: {e}")
    sys.exit(1)

# 创建DataFrame并保存
try:
    labels_df = pd.DataFrame(labels_data)
    pos_count = labels_df['label'].sum() if not labels_df.empty else 0
    neg_count = len(labels_df) - pos_count if not labels_df.empty else 0
    
    print(f"\n统计信息:")
    print(f"- 找到{len(labels_df)}个有效的复合物ID")
    print(f"- 其中有{pos_count}个正样本(标签=1)和{neg_count}个负样本(标签=0)")
    print(f"- 未能匹配标签的复合物: {not_found_count}个")
    
    if not labels_df.empty:
        labels_df.to_csv(output_path, index=False)
        print(f"\n已成功生成labels.csv文件: {output_path}")
    else:
        print("\n错误: 没有找到有效的数据，未生成labels.csv文件")
except Exception as e:
    print(f"错误: 生成CSV文件时出错: {e}")
    sys.exit(1)

print("\n脚本执行完毕") 