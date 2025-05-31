#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import pickle
import multiprocessing
from tqdm import tqdm
import torch
import dgl
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from scipy.spatial import distance_matrix
import subprocess

def one_of_k_encoding_unk(x, allowable_set):
    """对原子特征进行独热编码"""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))

def atom_features(atom, atom_symbols=['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I']):
    """提取原子特征"""
    results = one_of_k_encoding_unk(atom.GetSymbol(), atom_symbols + ['Unknown']) + \
            one_of_k_encoding_unk(atom.GetDegree(),[0, 1, 2, 3, 4, 5, 6]) + \
            one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5, 6]) + \
            one_of_k_encoding_unk(atom.GetHybridization(), [
                Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
                Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D,
                Chem.rdchem.HybridizationType.SP3D2
                ]) + [atom.GetIsAromatic()]
    # 处理氢原子
    results = results + one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4])
    
    return np.array(results).astype(np.float32)

def cal_dist(pos1, pos2, ord=2):
    """计算两点之间的距离"""
    return np.linalg.norm(pos1 - pos2, ord=ord)

def area_triangle(vector1, vector2):
    """计算两个向量形成的三角形面积"""
    return 0.5 * np.linalg.norm(np.cross(vector1, vector2))

def angle(vector1, vector2):
    """计算两个向量之间的角度"""
    inner = np.inner(vector1, vector2)
    norms = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    cos = inner / (norms + 1e-6)
    cos = np.clip(cos, -1.0, 1.0)
    return np.arccos(cos)

def run_obabel(input_file, output_file, input_format=None, output_format=None):
    """运行OpenBabel转换格式"""
    cmd = ["obabel"]
    if input_format:
        cmd.extend(["-i", input_format])
    cmd.append(input_file)
    if output_format:
        cmd.extend(["-o", output_format])
    cmd.append("-O")
    cmd.append(output_file)
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"OpenBabel转换失败: {e}")
        print(f"错误信息: {e.stderr}")
        return False

def prepare_complex(ligand_dir, cutoff=5.0):
    """准备配体和口袋的RDKit对象，并转换配体格式"""
    ligand_id = os.path.basename(ligand_dir)
    protein_dir = os.path.dirname(ligand_dir)
    
    # 文件路径
    ligand_file = os.path.join(ligand_dir, f"{ligand_id}.sdf")
    pocket_file = os.path.join(ligand_dir, f"{ligand_id}_pocket_{cutoff}A.pdb")
    rdkit_file = os.path.join(ligand_dir, f"{ligand_id}.rdkit")
    
    # 目标文件路径
    ligand_mol2_file = os.path.join(ligand_dir, f"{ligand_id}_ligand.mol2")
    ligand_pdb_file = os.path.join(ligand_dir, f"{ligand_id}_ligand.pdb")
    
    # 检查文件是否存在
    if not os.path.exists(ligand_file) or not os.path.exists(pocket_file):
        print(f"配体或口袋文件不存在: {ligand_file}, {pocket_file}")
        return False
    
    # 转换配体文件到其他格式
    try:
        # 使用OpenBabel转换为MOL2
        if not os.path.exists(ligand_mol2_file):
            if not run_obabel(ligand_file, ligand_mol2_file, output_format="mol2"):
                print(f"转换配体为MOL2格式失败: {ligand_mol2_file}")
                return False
            print(f"已生成配体MOL2文件: {ligand_mol2_file}")
            
        # 使用OpenBabel转换为PDB
        if not os.path.exists(ligand_pdb_file):
            if not run_obabel(ligand_file, ligand_pdb_file, output_format="pdb"):
                print(f"转换配体为PDB格式失败: {ligand_pdb_file}")
                return False
            print(f"已生成配体PDB文件: {ligand_pdb_file}")
    except Exception as e:
        print(f"处理配体文件失败: {str(e)}")
        return False
    
    # 如果rdkit文件已存在，跳过
    if os.path.exists(rdkit_file):
        print(f"RDKit文件已存在: {rdkit_file}")
        return True
    
    try:
        # 加载配体和口袋 - 使用EHIGN风格
        ligand = Chem.SDMolSupplier(ligand_file)[0]
        if not ligand:
            print(f"无法加载配体文件: {ligand_file}")
            return False
            
        # EHIGN中的口袋处理方式 - 移除氢原子，更简单更宽容
        pocket = Chem.MolFromPDBFile(pocket_file, removeHs=True) 
        if not pocket:
            print(f"无法加载口袋文件: {pocket_file}")
            return False
            
        # 保存RDKit对象
        with open(rdkit_file, 'wb') as f:
            pickle.dump((ligand, pocket), f)
            
        print(f"成功生成RDKit文件: {rdkit_file}")
        return True
        
    except Exception as e:
        print(f"生成RDKit对象时出错: {str(e)}")
        return False

def edge_feats(mol, graph):
    """提取分子中的边特征"""
    geom = mol.GetConformers()[0].GetPositions()
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        angles_ijk = []
        areas_ijk = []
        dists_ik = []
        for neighbor in mol.GetAtomWithIdx(j).GetNeighbors():
            k = neighbor.GetIdx() 
            if mol.GetBondBetweenAtoms(j, k) is not None and i != k:
                vector1 = geom[j] - geom[i]
                vector2 = geom[k] - geom[i]

                angles_ijk.append(angle(vector1, vector2))
                areas_ijk.append(area_triangle(vector1, vector2))
                dists_ik.append(cal_dist(geom[i], geom[k]))

        angles_ijk = np.array(angles_ijk) if angles_ijk != [] else np.array([0.])
        areas_ijk = np.array(areas_ijk) if areas_ijk != [] else np.array([0.])
        dists_ik = np.array(dists_ik) if dists_ik != [] else np.array([0.])
        dist_ij1 = cal_dist(geom[i], geom[j], ord=1)
        dist_ij2 = cal_dist(geom[i], geom[j], ord=2)

        geom_feats = [
            angles_ijk.max()*0.1,
            angles_ijk.sum()*0.01,
            angles_ijk.mean()*0.1,
            areas_ijk.max()*0.1,
            areas_ijk.sum()*0.01,
            areas_ijk.mean()*0.1,
            dists_ik.max()*0.1,
            dists_ik.sum()*0.01,
            dists_ik.mean()*0.1,
            dist_ij1*0.1,
            dist_ij2*0.1,
        ]

        bond_type = bond.GetBondType()
        basic_feats = [
        bond_type == Chem.rdchem.BondType.SINGLE,
        bond_type == Chem.rdchem.BondType.DOUBLE,
        bond_type == Chem.rdchem.BondType.TRIPLE,
        bond_type == Chem.rdchem.BondType.AROMATIC,
        bond.GetIsConjugated(),
        bond.IsInRing()]

        graph.add_edge(i, j, feats=torch.tensor(basic_feats+geom_feats).float())

def mol2graph(mol):
    """将分子转换为图"""
    graph = nx.Graph()
    # 添加节点特征
    for i, atom in enumerate(mol.GetAtoms()):
        graph.add_node(i, feats=torch.from_numpy(atom_features(atom)))
    
    # 添加边特征
    edge_feats(mol, graph)

    graph = graph.to_directed()
    x = torch.stack([feats['feats'] for n, feats in graph.nodes(data=True)])
    edge_index = torch.stack([torch.LongTensor((u, v)) for u, v, feats in graph.edges(data=True)]).T
    edge_attr = torch.stack([feats['feats'] for u, v, feats in graph.edges(data=True)])

    return x, edge_index, edge_attr

def geom_feats_func(pos_i, pos_j, angles_ijk, areas_ijk, dists_ik):
    """计算几何特征"""
    angles_ijk = np.array(angles_ijk) if angles_ijk != [] else np.array([0.])
    areas_ijk = np.array(areas_ijk) if areas_ijk != [] else np.array([0.])
    dists_ik = np.array(dists_ik) if dists_ik != [] else np.array([0.])
    dist_ij1 = cal_dist(pos_i, pos_j, ord=1)
    dist_ij2 = cal_dist(pos_i, pos_j, ord=2)

    geom = [
        angles_ijk.max()*0.1,
        angles_ijk.sum()*0.01,
        angles_ijk.mean()*0.1,
        areas_ijk.max()*0.1,
        areas_ijk.sum()*0.01,
        areas_ijk.mean()*0.1,
        dists_ik.max()*0.1,
        dists_ik.sum()*0.01,
        dists_ik.mean()*0.1,
        dist_ij1*0.1,
        dist_ij2*0.1,
    ]

    return geom

def geom_feat_func(pos_i, pos_j, pos_k, angles_ijk, areas_ijk, dists_ik):
    """计算几何特征辅助函数"""
    vector1 = pos_j - pos_i
    vector2 = pos_k - pos_i
    angles_ijk.append(angle(vector1, vector2))
    areas_ijk.append(area_triangle(vector1, vector2))
    dists_ik.append(cal_dist(pos_i, pos_k))

def inter_graph(ligand, pocket, dis_threshold=5.0):
    """构建配体和口袋之间的相互作用图"""
    graph_l2p = nx.DiGraph()
    graph_p2l = nx.DiGraph()
    pos_l = ligand.GetConformers()[0].GetPositions()
    pos_p = pocket.GetConformers()[0].GetPositions()
    dis_matrix = distance_matrix(pos_l, pos_p)
    node_idx = np.where(dis_matrix < dis_threshold)
    
    for i, j in zip(node_idx[0], node_idx[1]):
        ks = node_idx[0][node_idx[1] == j]
        angles_ijk = []
        areas_ijk = []
        dists_ik = []
        for k in ks:
            if k != i:
                geom_feat_func(pos_l[i], pos_p[j], pos_l[k], angles_ijk, areas_ijk, dists_ik)
        geom = geom_feats_func(pos_l[i], pos_p[j], angles_ijk, areas_ijk, dists_ik)
        bond_feats = torch.FloatTensor(geom)
        graph_l2p.add_edge(i, j, feats=bond_feats)
        
        ks = node_idx[1][node_idx[0] == i]
        angles_ijk = []
        areas_ijk = []
        dists_ik = []
        for k in ks:
            if k != j:
                geom_feat_func(pos_p[j], pos_l[i], pos_p[k], angles_ijk, areas_ijk, dists_ik)     
        geom = geom_feats_func(pos_p[j], pos_l[i], angles_ijk, areas_ijk, dists_ik)
        bond_feats = torch.FloatTensor(geom)
        graph_p2l.add_edge(j, i, feats=bond_feats)
    
    edge_index_l2p = torch.stack([torch.LongTensor((u, v)) for u, v, feats in graph_l2p.edges(data=True)]).T if graph_l2p.edges() else torch.zeros((2, 0), dtype=torch.long)
    edge_attr_l2p = torch.stack([feats['feats'] for u, v, feats in graph_l2p.edges(data=True)]) if graph_l2p.edges() else torch.zeros((0, 11), dtype=torch.float)

    edge_index_p2l = torch.stack([torch.LongTensor((u, v)) for u, v, feats in graph_p2l.edges(data=True)]).T if graph_p2l.edges() else torch.zeros((2, 0), dtype=torch.long)
    edge_attr_p2l = torch.stack([feats['feats'] for u, v, feats in graph_p2l.edges(data=True)]) if graph_p2l.edges() else torch.zeros((0, 11), dtype=torch.float)

    return (edge_index_l2p, edge_attr_l2p), (edge_index_p2l, edge_attr_p2l)

def generate_dgl_graph(ligand_dir, cutoff=5.0, dis_threshold=5.0):
    """生成DGL格式的异构图"""
    ligand_id = os.path.basename(ligand_dir)
    rdkit_file = os.path.join(ligand_dir, f"{ligand_id}.rdkit")
    dgl_file = os.path.join(ligand_dir, f"{ligand_id}.dgl")
    
    # 检查文件是否存在
    if not os.path.exists(rdkit_file):
        print(f"RDKit文件不存在: {rdkit_file}")
        return False
    
    # 如果DGL文件已存在，跳过
    if os.path.exists(dgl_file):
        print(f"DGL图文件已存在: {dgl_file}")
        return True
        
    try:
        # 按照EHIGN风格处理
        with open(rdkit_file, 'rb') as f:
            ligand, pocket = pickle.load(f)
            
        atom_num_l = ligand.GetNumAtoms()
        atom_num_p = pocket.GetNumAtoms()
        
        if atom_num_l == 0 or atom_num_p == 0:
            print(f"配体或口袋原子数为零: {rdkit_file}")
            return False
            
        x_l, edge_index_l, edge_attr_l = mol2graph(ligand)
        x_p, edge_index_p, edge_attr_p = mol2graph(pocket)
        (edge_index_l2p, edge_attr_l2p), (edge_index_p2l, edge_attr_p2l) = inter_graph(ligand, pocket, dis_threshold=dis_threshold)
        
        # 创建异构图数据
        graph_data = {
            ('ligand', 'intra_l', 'ligand'): (edge_index_l[0], edge_index_l[1]),
            ('pocket', 'intra_p', 'pocket'): (edge_index_p[0], edge_index_p[1])
        }
        
        # 只有当边存在时才添加
        if edge_index_l2p.shape[1] > 0:
            graph_data[('ligand', 'inter_l2p', 'pocket')] = (edge_index_l2p[0], edge_index_l2p[1])
        if edge_index_p2l.shape[1] > 0:
            graph_data[('pocket', 'inter_p2l', 'ligand')] = (edge_index_p2l[0], edge_index_p2l[1])
            
        # 创建异构图
        g = dgl.heterograph(graph_data, num_nodes_dict={"ligand":atom_num_l, "pocket":atom_num_p})
        
        # 添加特征
        g.nodes['ligand'].data['h'] = x_l
        g.nodes['pocket'].data['h'] = x_p
        g.edges['intra_l'].data['e'] = edge_attr_l
        g.edges['intra_p'].data['e'] = edge_attr_p
        
        # 只有当边存在时才添加特征
        if 'inter_l2p' in g.etypes and edge_attr_l2p.shape[0] > 0:
            g.edges['inter_l2p'].data['e'] = edge_attr_l2p
        if 'inter_p2l' in g.etypes and edge_attr_p2l.shape[0] > 0:
            g.edges['inter_p2l'].data['e'] = edge_attr_p2l
            
        # 检查特征是否包含NaN
        if torch.isnan(x_l).any() or torch.isnan(x_p).any() or \
           torch.isnan(edge_attr_l).any() or torch.isnan(edge_attr_p).any():
            print(f"特征包含NaN值: {dgl_file}")
            return False
            
        # 保存为DGL图文件
        torch.save((g, torch.FloatTensor([0.0])), dgl_file)  # 使用0.0作为标签占位符
        print(f"成功生成DGL图文件: {dgl_file}")
        return True
        
    except Exception as e:
        print(f"生成DGL图时出错: {str(e)}")
        return False

def process_ligand(ligand_dir, cutoff=5.0, dis_threshold=5.0):
    """处理单个配体目录"""
    try:
        # 步骤1: 准备RDKit复合物
        if not prepare_complex(ligand_dir, cutoff):
            return False
            
        # 步骤2: 生成DGL图
        if not generate_dgl_graph(ligand_dir, cutoff, dis_threshold):
            return False
            
        return True
        
    except Exception as e:
        print(f"处理配体时出错: {ligand_dir}, {str(e)}")
        return False

def process_all_ligands(base_dir, cutoff=5.0, dis_threshold=5.0, num_processes=4):
    """处理所有配体"""
    # 收集所有配体目录
    all_ligand_dirs = []
    for protein_dir in glob.glob(os.path.join(base_dir, "*")):
        if not os.path.isdir(protein_dir):
            continue
            
        protein_id = os.path.basename(protein_dir)
        for ligand_dir in glob.glob(os.path.join(protein_dir, "*")):
            if not os.path.isdir(ligand_dir):
                continue
                
            ligand_id = os.path.basename(ligand_dir)
            if ligand_id.startswith(protein_id):
                all_ligand_dirs.append(ligand_dir)
                
    print(f"找到 {len(all_ligand_dirs)} 个配体目录")
    
    # 使用多进程处理所有配体
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = list(tqdm(
            pool.starmap(
                process_ligand, 
                [(ligand_dir, cutoff, dis_threshold) for ligand_dir in all_ligand_dirs]
            ),
            total=len(all_ligand_dirs),
            desc="处理配体"
        ))
        
    success_count = sum(results)
    print(f"成功处理 {success_count}/{len(all_ligand_dirs)} 个配体")
    return success_count

def main():
    parser = argparse.ArgumentParser(description="生成EHIGN格式数据文件（按照原始EHIGN逻辑）")
    parser.add_argument("--base_dir", type=str, default="train_pdb", 
                        help="包含蛋白质和配体数据的基础目录")
    parser.add_argument("--cutoff", type=float, default=5.0, 
                        help="口袋半径(埃)")
    parser.add_argument("--dis_threshold", type=float, default=5.0, 
                        help="配体-口袋相互作用距离(埃)")
    parser.add_argument("--processes", type=int, default=4, 
                        help="并行进程数")
                        
    args = parser.parse_args()
    
    process_all_ligands(args.base_dir, args.cutoff, args.dis_threshold, args.processes)
    
    print(f"\n{'='*80}")
    print(f"EHIGN格式数据处理完成!")
    print(f"{'='*80}")
    
if __name__ == "__main__":
    main()
