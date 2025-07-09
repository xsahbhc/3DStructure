#!/bin/bash

# 3D数据处理环境安装脚本
# 用于创建专门的虚拟环境来处理BioSNAP和BindingDB数据

set -e  # 遇到错误时退出

echo "🚀 开始创建3D数据处理虚拟环境..."

# 环境名称
ENV_NAME="data"
PYTHON_VERSION="3.8"

# 检查conda是否可用
if ! command -v conda &> /dev/null; then
    echo "❌ 错误: conda未找到，请先安装Anaconda或Miniconda"
    exit 1
fi

# 检查环境是否已存在
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "⚠️  环境 ${ENV_NAME} 已存在"
    read -p "是否删除现有环境并重新创建? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "🗑️  删除现有环境..."
        conda env remove -n ${ENV_NAME} -y
    else
        echo "❌ 取消安装"
        exit 1
    fi
fi

# 创建新的conda环境
echo "📦 创建conda环境: ${ENV_NAME} (Python ${PYTHON_VERSION})"
conda create -n ${ENV_NAME} python=${PYTHON_VERSION} -y

# 激活环境
echo "🔄 激活环境..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ${ENV_NAME}


# 验证环境激活
if [[ "$CONDA_DEFAULT_ENV" != "$ENV_NAME" ]]; then
    echo "❌ 错误: 环境激活失败"
    exit 1
fi

echo "✅ 环境激活成功: $CONDA_DEFAULT_ENV"

# 更新pip
echo "📦 更新pip..."
pip install --upgrade pip

# 安装基础科学计算包
echo "📦 安装基础科学计算包..."
pip install numpy pandas matplotlib seaborn tqdm

# 安装化学信息学包
echo "🧪 安装化学信息学包..."
# 配置conda使用libmamba求解器加速
conda config --set solver libmamba
# 使用conda安装rdkit (使用libmamba求解器加速)
conda install -c conda-forge rdkit -y

# 安装网络请求包
echo "🌐 安装网络请求包..."
pip install aiohttp requests urllib3

# 安装机器学习包
echo "🤖 安装机器学习包..."
pip install scikit-learn

# 安装图处理包
echo "📊 安装图处理包..."
# 安装DGL (CPU版本)
pip install dgl -f https://data.dgl.ai/wheels/repo.html

# 安装PyTorch (CPU版本，用于基础图处理)
echo "🔥 安装PyTorch (CPU版本)..."
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# 安装生物信息学工具
echo "🧬 安装生物信息学工具..."
pip install biopython

# 检查系统工具
echo "🔧 检查系统工具..."

# 检查DIAMOND
if command -v ~/bin/diamond &> /dev/null; then
    echo "✅ DIAMOND已安装: $(~/bin/diamond --version | head -1)"
else
    echo "⚠️  DIAMOND未找到，请确保已安装并位于 ~/bin/diamond"
    echo "   下载地址: https://github.com/bbuchfink/diamond/releases"
fi

# 检查PyMOL
if command -v pymol &> /dev/null; then
    echo "✅ PyMOL已安装"
else
    echo "⚠️  PyMOL未找到，尝试安装..."
    # 尝试通过conda安装PyMOL (使用libmamba求解器)
    conda install -c conda-forge pymol-open-source -y || echo "❌ PyMOL安装失败，请手动安装"
fi

# 检查OpenBabel
if command -v obabel &> /dev/null; then
    echo "✅ OpenBabel已安装: $(obabel -V)"
else
    echo "⚠️  OpenBabel未找到，尝试安装..."
    conda install -c conda-forge openbabel -y || echo "❌ OpenBabel安装失败，请手动安装"
fi

# 创建环境激活脚本
echo "📝 创建环境激活脚本..."
cat > activate_data_processing.sh << 'EOF'
#!/bin/bash
# 激活数据处理环境的便捷脚本

echo "🔄 激活数据处理环境..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate data

if [[ "$CONDA_DEFAULT_ENV" == "data" ]]; then
    echo "✅ 环境激活成功: $CONDA_DEFAULT_ENV"
    echo "📍 当前Python路径: $(which python)"
    echo "📦 已安装的主要包:"
    echo "   - pandas: $(python -c 'import pandas; print(pandas.__version__)')"
    echo "   - rdkit: $(python -c 'from rdkit import rdBase; print(rdBase.rdkitVersion)')"
    echo "   - numpy: $(python -c 'import numpy; print(numpy.__version__)')"
    echo "   - scikit-learn: $(python -c 'import sklearn; print(sklearn.__version__)')"
    echo ""
    echo "🚀 现在可以运行数据处理脚本了！"
    echo "   cd /home/work/workspace/shi_shaoqun/snap/3D_structure"
    echo "   python get_pdbid.py"
else
    echo "❌ 环境激活失败"
    exit 1
fi
EOF

chmod +x activate_data_processing.sh

# 创建依赖检查脚本
echo "📝 创建依赖检查脚本..."
cat > check_dependencies.py << 'EOF'
#!/usr/bin/env python3
"""
检查数据处理环境的依赖包是否正确安装
"""

import sys
import importlib

def check_package(package_name, import_name=None):
    """检查包是否可以导入"""
    if import_name is None:
        import_name = package_name

    try:
        module = importlib.import_module(import_name)
        version = getattr(module, '__version__', 'unknown')
        print(f"✅ {package_name}: {version}")
        return True
    except ImportError as e:
        print(f"❌ {package_name}: 导入失败 - {e}")
        return False

def main():
    print("🔍 检查数据处理环境依赖...")
    print(f"🐍 Python版本: {sys.version}")
    print()

    # 必需的包
    required_packages = [
        ('pandas', 'pandas'),
        ('numpy', 'numpy'),
        ('rdkit', 'rdkit'),
        ('tqdm', 'tqdm'),
        ('aiohttp', 'aiohttp'),
        ('scikit-learn', 'sklearn'),
        ('dgl', 'dgl'),
        ('torch', 'torch'),
        ('biopython', 'Bio'),
    ]

    success_count = 0
    total_count = len(required_packages)

    for package_name, import_name in required_packages:
        if check_package(package_name, import_name):
            success_count += 1

    print()
    print(f"📊 依赖检查结果: {success_count}/{total_count} 包成功安装")

    if success_count == total_count:
        print("🎉 所有依赖包都已正确安装！")
        return True
    else:
        print("⚠️  部分依赖包安装失败，请检查上述错误信息")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
EOF

chmod +x check_dependencies.py

# 运行依赖检查
echo "🔍 运行依赖检查..."
python check_dependencies.py

# 创建使用说明
echo "📝 创建使用说明..."
cat > README_DATA_PROCESSING.md << 'EOF'
# 3D数据处理环境使用指南

## 环境激活

每次使用前，请先激活环境：

```bash
# 方法1: 使用便捷脚本
./activate_data_processing.sh

# 方法2: 手动激活
conda activate data_processing
```

## 依赖检查

检查所有依赖是否正确安装：

```bash
python check_dependencies.py
```

## 数据处理流程

激活环境后，按以下顺序运行脚本：

```bash
cd /home/work/workspace/shi_shaoqun/snap/3D_structure

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

# 7. 验证结构一致性
python verify_biosnap_structure.py
```

## 已安装的主要包

- **pandas**: 数据处理
- **rdkit**: 化学信息学
- **numpy**: 数值计算
- **tqdm**: 进度条
- **aiohttp**: 异步HTTP请求
- **scikit-learn**: 机器学习
- **dgl**: 图神经网络
- **torch**: 深度学习框架
- **biopython**: 生物信息学

## 外部工具

请确保以下工具已安装：

- **DIAMOND**: 序列比对工具 (~/bin/diamond)
- **PyMOL**: 分子可视化工具
- **OpenBabel**: 化学格式转换工具

## 故障排除

如果遇到问题：

1. 检查环境是否正确激活：`echo $CONDA_DEFAULT_ENV`
2. 运行依赖检查：`python check_dependencies.py`
3. 重新安装环境：`./setup_data_processing_env.sh`
EOF

echo ""
echo "🎉 数据处理环境安装完成！"
echo ""
echo "📋 下一步操作："
echo "1. 激活环境: ./activate_data_processing.sh"
echo "2. 检查依赖: python check_dependencies.py"
echo "3. 开始数据处理: cd /home/work/workspace/shi_shaoqun/snap/3D_structure"
echo ""
echo "📖 详细使用说明请查看: README_DATA_PROCESSING.md"
