import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split

print("开始执行分层数据集划分...")

# 配置路径
labels_file = '3D_structure/labels.csv'
output_dir = '3D_structure'

# 创建输出目录（如果不存在）
os.makedirs(output_dir, exist_ok=True)

# 读取标签文件
print(f"读取标签文件: {labels_file}")
try:
    df = pd.read_csv(labels_file)
    print(f"成功读取标签文件，共有 {len(df)} 条记录")
    print(f"正样本数量: {df['label'].sum()}, 正样本比例: {df['label'].mean():.2%}")
except Exception as e:
    print(f"读取文件出错: {e}")
    exit(1)

# 确保数据没有重复
df = df.drop_duplicates(subset=['complex_id'])

# 使用分层采样划分数据集 (8:1:1)
print("使用分层采样划分数据集...")
try:
    # 首先划分训练集和临时集(验证+测试)
    train_df, temp_df = train_test_split(
        df, 
        test_size=0.2,  # 20%用于验证+测试
        stratify=df['label'],  # 按标签分层
        random_state=42  # 固定随机种子，确保可重复性
    )
    
    # 再将临时集划分为验证集和测试集
    val_df, test_df = train_test_split(
        temp_df,
        test_size=0.5,  # 临时集的一半作为测试集
        stratify=temp_df['label'],  # 继续按标签分层
        random_state=42
    )
    
    # 打印统计信息
    print(f"训练集: {len(train_df)} 条记录, 正样本比例: {train_df['label'].mean():.2%}")
    print(f"验证集: {len(val_df)} 条记录, 正样本比例: {val_df['label'].mean():.2%}")
    print(f"测试集: {len(test_df)} 条记录, 正样本比例: {test_df['label'].mean():.2%}")
    
    # 保存为CSV文件
    train_file = os.path.join(output_dir, 'train_stratified.csv')
    val_file = os.path.join(output_dir, 'val_stratified.csv')
    test_file = os.path.join(output_dir, 'test_stratified.csv')
    
    train_df.to_csv(train_file, index=False)
    val_df.to_csv(val_file, index=False)
    test_df.to_csv(test_file, index=False)
    
    print(f"已保存训练集文件: {train_file}")
    print(f"已保存验证集文件: {val_file}")
    print(f"已保存测试集文件: {test_file}")
    
except Exception as e:
    print(f"划分数据集出错: {e}")
    exit(1)

# 生成DrugBAN项目配置建议
print("\n要在DrugBAN项目中使用这些分层划分的数据集，请修改以下配置:")
print("在main.py中，修改数据加载部分为:")
print("""
# 数据路径配置
train_path = os.path.join(dataFolder, "train_stratified.csv")
val_path = os.path.join(dataFolder, "val_stratified.csv")
test_path = os.path.join(dataFolder, "test_stratified.csv")

# 加载数据
df_train = pd.read_csv(train_path)
df_val = pd.read_csv(val_path)
df_test = pd.read_csv(test_path)
""")

print("或者，在DrugBAN3D.yaml中，添加数据集路径配置:")
print("""
DATA:
  TRAIN_FILE: "3D_structure/train_stratified.csv"
  VAL_FILE: "3D_structure/val_stratified.csv"
  TEST_FILE: "3D_structure/test_stratified.csv"
""")

print("\n另外，您可能需要修改dataloader_3d.py中的get_loader函数，使其读取这些预划分的数据集")
print("脚本执行完成!") 