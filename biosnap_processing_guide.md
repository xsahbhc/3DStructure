# BioSNAP数据集完整处理流程指南

## 🎯 概述

本指南详细描述了从原始BioSNAP数据到完整多模态训练数据的全过程，包括1D/2D序列数据和3D结构数据的生成。

## 📋 完整流程概览

### **阶段一：环境准备**
### **阶段二：1D/2D数据生成**
### **阶段三：3D数据生成（7步流程）**
### **阶段四：数据验证**

---

## 🔧 阶段一：环境准备

### 1. 环境配置脚本
- `3D_structure/setup_data_processing_env.sh` - 创建专用conda环境
- `3D_structure/activate_data_processing.sh` - 环境激活脚本
- `3D_structure/check_dependencies.py` - 依赖检查工具

```bash
# 设置数据处理环境
cd /home/work/workspace/shi_shaoqun/snap/3D_structure
./setup_data_processing_env.sh

## 📊 阶段二：1D/2D数据生成

### 2. 序列数据提取
**脚本**: `3D_structure/extract_biosnap_sequences.py` - 生成DrugBAN格式的1D/2D数据

```bash
# 生成1D/2D序列数据
cd /home/work/workspace/shi_shaoqun/snap
python drugban/extract_real_sequences.py
python 3D_structure/extract_biosnap_sequences.py
```

**输出**:
- `3D_structure/biosnap/biosnap_3d_sequences/random/train_stratified.csv` (6028条)
- `3D_structure/biosnap/biosnap_3d_sequences/random/val_stratified.csv` (753条)
- `3D_structure/biosnap/biosnap_3d_sequences/random/test_stratified.csv` (754条)
- `3D_structure/biosnap/biosnap_3d_sequences/random/seqid_mapping.csv` (7535条)

---

## 🧬 阶段三：3D数据生成（7步流程）

### 3. PDB ID获取
**脚本**: `3D_structure/get_pdbid.py`
**功能**: 使用DIAMOND工具从蛋白质序列获取PDB ID

**输入文件**: `biosnap/train_csv/train.csv`
**输出文件**:
- `biosnap/train_csv/with_pdbid.csv` (成功获取PDB ID的数据)
- `biosnap/train_csv/without_pdbid.csv` (未能获取PDB ID的数据)
- `biosnap/temp.fasta` (临时文件，会自动清理)

### 4. PDB文件下载
**脚本**: `3D_structure/pdb_download.py`
**功能**: 异步下载PDB结构文件

**输入文件**: `biosnap/train_csv/with_pdbid.csv`
**输出目录**: `biosnap/train_pdb/` (下载的PDB文件)

### 5. 配体SDF生成
**脚本**: `3D_structure/generate_ligand_sdf.py`
**功能**: 从SMILES生成3D配体结构

**输入文件**: `biosnap/train_csv/with_pdbid.csv`
**输入目录**: `biosnap/train_pdb/`
**输出**: 在各蛋白质目录下生成配体SDF文件
**日志文件**: `biosnap/sdf_generation_failures.log`

### 6. 蛋白质口袋生成
**脚本**: `3D_structure/generate_pocket_simple.py`
**功能**: 生成蛋白质口袋和多种格式文件

**默认输入目录**: `biosnap/train_pdb/`
**输出**: 在各蛋白质目录下生成口袋PDB文件和其他格式文件

### 7. 标签文件生成
**脚本**: `3D_structure/generate_labels.py`
**功能**: 生成完整的标签映射文件

**输入文件**: `biosnap/train_csv/train.csv`
**输入目录**: `biosnap/train_pdb/`
**输出文件**: `biosnap/train_csv/labels.csv`

### 8. 数据分层划分
**脚本**: `3D_structure/stratify_split.py`
**功能**: 生成训练/验证/测试集划分

**输出文件**:
- `biosnap/train_csv/train_stratified.csv` (6028条)
- `biosnap/train_csv/val_stratified.csv` (753条)
- `biosnap/train_csv/test_stratified.csv` (754条)

### 9. 图数据预处理
**脚本**: `DrugBAN-main/pre_cache.py`
**功能**: 生成DGL图和RDKit对象

**输出**: 各复合物目录下的`.dgl`和`.rdkit`文件

## 🚀 完整执行顺序

```bash
cd /home/work/workspace/shi_shaoqun/snap/3D_structure

# 激活环境
conda activate data_processing

# 1. 获取PDB ID
python get_pdbid.py

# 2. 下载PDB文件
python pdb_download.py

# 3. 生成配体SDF文件
python generate_ligand_sdf.py

# 4. 生成蛋白质口袋和图文件
python generate_protein_pockets.py

# 5. 生成标签文件
python generate_labels.py

# 6. 生成分层划分文件
python stratify_split.py biosnap

# 7. 预处理图数据 (在DrugBAN-main目录中)
cd ../DrugBAN-main
python pre_cache.py --data_dir ../3D_structure/biosnap/train_pdb

# 8. 验证结构一致性
cd ../3D_structure
python verify_biosnap_structure.py
```

---

## ✅ 阶段四：数据验证

### 10. 结构验证
**脚本**: `3D_structure/verify_biosnap_structure.py`
**功能**: 验证目录结构完整性

### 11. 数据一致性检查
**脚本**: `3D_structure/check_data_consistency.py`
**功能**: 验证1D/2D与3D数据一致性

```bash
# 验证数据一致性
cd /home/work/workspace/shi_shaoqun/snap/3D_structure
conda activate drugban
python check_data_consistency.py
```

---

## 📁 最终生成的完整数据结构

```
3D_structure/biosnap/
├── biosnap_3d_sequences/random/          # 1D/2D数据
│   ├── train_stratified.csv              # 训练集 (6028条)
│   ├── val_stratified.csv                # 验证集 (753条)
│   ├── test_stratified.csv               # 测试集 (754条)
│   └── seqid_mapping.csv                 # seqid映射文件 (7535条)
├── train_csv/                            # 3D数据标签
│   ├── train.csv                         # 原始BioSNAP数据
│   ├── with_pdbid.csv                    # 成功获取PDB ID的数据
│   ├── without_pdbid.csv                 # 未能获取PDB ID的数据
│   ├── labels.csv                        # 完整标签文件 (7535条)
│   ├── train_stratified.csv              # 3D训练集 (6028条)
│   ├── val_stratified.csv                # 3D验证集 (753条)
│   └── test_stratified.csv               # 3D测试集 (754条)
├── sdf_generation_failures.log           # SDF生成失败日志
└── train_pdb/                            # 3D结构数据
    └── [PDBID_CHAIN]/                    # 蛋白质目录 (如: 1ACB_E)
        ├── [PDBID_CHAIN].pdb             # 蛋白质结构文件
        └── [PDBID_CHAIN_SEQID]/          # 复合物子目录 (如: 1ACB_E_12171)
            ├── [COMPLEX_ID].sdf          # 配体SDF文件
            ├── [COMPLEX_ID]_ligand.pdb   # 配体PDB文件
            ├── [COMPLEX_ID]_ligand.mol2  # 配体MOL2文件
            ├── [COMPLEX_ID]_pocket_5.0A.pdb  # 蛋白质口袋文件
            ├── [COMPLEX_ID].dgl          # DGL图对象
            └── [COMPLEX_ID].rdkit        # RDKit分子对象
```

### 🏷️ 关键命名规则
1. **蛋白质目录**: `PDBID_CHAIN` (如: `1ACB_E`, `7RMG_R`)
2. **复合物子目录**: `PDBID_CHAIN_SEQID` (如: `1ACB_E_12171`)
3. **文件命名**:
   - 蛋白质: `PDBID_CHAIN.pdb`
   - 配体SDF: `PDBID_CHAIN_SEQID.sdf`
   - 配体PDB: `PDBID_CHAIN_SEQID_ligand.pdb`
   - 配体MOL2: `PDBID_CHAIN_SEQID_ligand.mol2`
   - 口袋: `PDBID_CHAIN_SEQID_pocket_5.0A.pdb`
   - DGL图: `PDBID_CHAIN_SEQID.dgl`
   - RDKit: `PDBID_CHAIN_SEQID.rdkit`

---

## 📊 数据统计与验证结果

### **BioSNAP数据完全一致性**：
- ✅ **训练集**: 1D/2D和3D都是6028条，标签分布完全匹配
- ✅ **验证集**: 1D/2D和3D都是753条，标签分布完全匹配
- ✅ **测试集**: 1D/2D和3D都是754条，标签分布完全匹配
- ✅ **seqid映射**: 7535条记录，完全一致（0个缺失）

---

## 🔧 使用的所有脚本文件

### **核心处理脚本**：
1. `3D_structure/get_pdbid.py` - PDB ID获取
2. `3D_structure/pdb_download.py` - PDB文件下载
3. `3D_structure/generate_ligand_sdf.py` - 配体SDF生成
4. `3D_structure/generate_protein_pockets.py` - 蛋白质口袋生成
5. `3D_structure/generate_labels.py` - 标签文件生成
6. `3D_structure/stratify_split.py` - 数据分层划分
7. `DrugBAN-main/pre_cache.py` - 图数据预处理
8. `3D_structure/extract_biosnap_sequences.py` - 1D/2D数据生成

### **验证和工具脚本**：
9. `3D_structure/verify_biosnap_structure.py` - 结构验证
10. `3D_structure/check_data_consistency.py` - 数据一致性检查
11. `3D_structure/setup_data_processing_env.sh` - 环境配置
12. `3D_structure/check_dependencies.py` - 依赖检查

### **配置和文档**：
13. `3D_structure/biosnap_processing_guide.md` - 本指南文件
14. `3D_structure/INSTALL_GUIDE.md` - 安装指南
15. `3D_structure/README_CUSTOM.md` - 自定义工具说明
16. `DrugBAN-main/configs.py` - BioSNAP路径配置

---

## ⚠️ 注意事项

1. **环境依赖**: 确保安装了rdkit, pandas, tqdm等依赖包
2. **DIAMOND工具**: get_pdbid.py需要DIAMOND工具和PDB数据库
3. **存储空间**: 3D数据处理需要大量存储空间
4. **处理时间**: 整个流程可能需要数小时到数天，取决于数据量
5. **错误处理**: 每个步骤都有日志记录，注意检查错误信息
6. **内存要求**: 图数据预处理需要足够的内存
7. **网络连接**: PDB文件下载需要稳定的网络连接

---

## 🎯 成果确认

✅ **数据完整性**: 所有必需的文件都已生成
✅ **格式一致性**: 与BindingDB数据格式完全一致
✅ **数据对应性**: 1D/2D与3D数据完美匹配
✅ **可用性**: 可直接用于DrugBAN多模态训练

## 🚀 后续使用

处理完成后，BioSNAP数据可以直接用于：
- **多模态DrugBAN训练**: 使用`run_multimodal_optimized.sh biosnap`
- **与BindingDB数据集的性能对比**: 相同格式便于对比分析
- **跨数据集泛化能力测试**: 在不同数据集间测试模型性能
- **消融实验**: 比较1D/2D、3D和多模态方法的效果

---

## 📞 故障排除

如果遇到问题：

1. **环境问题**: 运行`python check_dependencies.py`检查依赖
2. **数据一致性**: 运行`python check_data_consistency.py`验证数据
3. **结构完整性**: 运行`python verify_biosnap_structure.py`检查文件
4. **重新处理**: 可以从任何步骤重新开始，脚本会跳过已存在的文件
