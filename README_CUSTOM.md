# DrugBAN-EHIGN 自定义数据处理工具

这个工具集用于处理蛋白质-配体数据，生成口袋文件并转换为EHIGN需要的数据格式，同时**保留原始命名格式**。

## 背景

原始EHIGN使用特定的命名格式，但在我们的数据集中，一个蛋白质可能对应多个配体，因此我们需要保留原始文件名以便于区分和管理。此工具将生成EHIGN所需的所有文件，但使用我们自己的命名规则。

## 工具清单

- `generate_protein_pockets.py` - 生成蛋白质口袋文件
- `generate_ehign_custom.py` - 使用自定义命名格式生成EHIGN所需的数据文件
- `process_custom_pipeline.py` - 完整的数据处理流程

## 前提条件

使用前请确保已安装：

```bash
# 基本依赖
pip install rdkit dgl torch numpy pandas tqdm

# 外部工具
sudo apt-get install pymol openbabel
```

## 数据目录结构

工具假设数据结构如下：

```
train_pdb/
├── [蛋白质ID]/
│   ├── [蛋白质ID].pdb  # 蛋白质结构文件
│   ├── [蛋白质ID]_[配体ID]/
│   │   ├── [蛋白质ID]_[配体ID].sdf  # 配体结构文件
```

处理后，将在每个配体目录中生成以下文件：

```
[蛋白质ID]_[配体ID]/
├── [蛋白质ID]_[配体ID].sdf  # 原始配体SDF文件
├── [蛋白质ID]_[配体ID].mol2 # 转换的配体MOL2文件
├── [蛋白质ID]_[配体ID].pdb  # 转换的配体PDB文件
├── [蛋白质ID]_[配体ID]_pocket_5.0A.pdb  # 生成的口袋文件
├── [蛋白质ID]_[配体ID].rdkit  # RDKit对象序列化文件
├── [蛋白质ID]_[配体ID].dgl  # DGL异构图文件
```

## 使用方法

### 1. 生成口袋文件

```bash
cd 3D_structure
python generate_protein_pockets.py --base_dir train_pdb --cutoff 5.0
```

参数说明：
- `--base_dir`: 数据根目录
- `--cutoff`: 口袋区域的距离阈值(埃)
- `--no_parallel`: 禁用并行处理

### 2. 生成自定义EHIGN文件

```bash
cd 3D_structure
python generate_ehign_custom.py --base_dir train_pdb --cutoff 5.0 --dis_threshold 5.0
```

参数说明：
- `--base_dir`: 数据根目录
- `--cutoff`: 口袋区域距离阈值(埃)
- `--dis_threshold`: 配体-口袋原子间距阈值(埃)
- `--no_rdkit`: 跳过生成RDKit对象
- `--no_dgl`: 跳过生成DGL图
- `--no_parallel`: 禁用并行处理

### 3. 运行完整流程

```bash
cd 3D_structure
python process_custom_pipeline.py --train_pdb train_pdb --cutoff 5.0 --dis_threshold 5.0
```

参数说明：
- `--train_pdb`: 数据根目录
- `--cutoff`: 口袋区域距离阈值(埃)
- `--dis_threshold`: 配体-口袋原子间距阈值(埃)
- `--no_pocket`: 跳过口袋生成步骤
- `--no_ehign`: 跳过EHIGN格式生成步骤
- `--no_rdkit`: 跳过RDKIT文件生成
- `--no_dgl`: 跳过DGL图文件生成
- `--no_parallel`: 禁用并行处理

## 日志文件

- `pocket_generation.log` - 口袋生成日志
- `ehign_generation.log` - EHIGN数据生成日志
- `custom_pipeline_[时间戳].log` - 完整流程日志

## 在DrugBAN中使用

要在DrugBAN中使用这些生成的数据，您需要修改数据加载部分以适应您的命名规则。以下是一个示例代码片段：

```python
# 在DrugBAN中加载自定义格式的DGL图
def load_dgl_graph(complex_id):
    """加载特定蛋白质-配体复合物的DGL图
    
    Args:
        complex_id (str): 形如'8XFL_A_5538'的复合物ID
        
    Returns:
        dgl.DGLGraph: 加载的异构图对象
    """
    graph_path = os.path.join(data_dir, complex_id, f"{complex_id}.dgl")
    if not os.path.exists(graph_path):
        return None
        
    try:
        g, _ = torch.load(graph_path)
        return g
    except:
        return None

# 在DrugBAN中集成EHIGN
class DrugBAN3D(nn.Module):
    def __init__(self):
        super(DrugBAN3D, self).__init__()
        # 原始DrugBAN组件
        self.drugban_components = ...
        
        # EHIGN组件
        self.hgnn = CustomHeteroGNN(...)
        
    def forward(self, batch):
        # 处理原始2D输入
        drugban_output = ...
        
        # 处理3D异构图输入
        if hasattr(batch, 'dgl_graph'):
            graph_embedding = self.hgnn(batch.dgl_graph)
            # 融合2D和3D特征
            combined_features = torch.cat([drugban_output, graph_embedding], dim=1)
            output = self.final_predictor(combined_features)
            return output
            
        # 只有2D输入的情况
        return self.final_predictor(drugban_output)
```

## 注意事项

- 处理大量数据时，建议使用`--no_parallel`选项减少内存使用
- DGL图生成是计算密集型任务，确保有足够内存
- EHIGN使用异构图表示蛋白质-配体相互作用，与DrugBAN集成时需要进行适当修改 