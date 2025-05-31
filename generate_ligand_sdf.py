import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import logging

# --- 配置 ---
# 基本目录，其中包含以PDB ID命名的蛋白质结构文件夹
TRAIN_PDB_BASE_DIR = 'train_pdb'
# 包含PDB_ID, SMILES, seqid信息的CSV文件路径
WITH_PDBID_CSV_PATH = 'with_pdbid.csv'

# CSV文件中的列名
PDB_ID_COL = 'PDB_ID'
SMILES_COL = 'SMILES'
SEQID_COL = 'seqid' # 用于唯一标识配体并用于文件夹/文件名
FAILURE_LOG_FILE = 'sdf_generation_failures.log'

# --- 设置日志记录 ---
# 配置日志记录器，用于记录失败信息
logging.basicConfig(filename=FAILURE_LOG_FILE, level=logging.INFO, 
                    format='%(asctime)s - PDB: %(pdb_id)s - SEQID: %(seqid)s - SMILES: %(smiles)s - REASON: %(reason)s',
                    filemode='w') # 'w' 模式会在每次运行时覆盖旧日志

logger = logging.getLogger(__name__)

# --- 工具函数 ---
def sanitize_name_component(name_component):
    """对用于文件名或目录名的部分进行基本清理。"""
    # 将组件转换为字符串，并替换可能引起问题的字符
    # 仅保留字母数字字符、下划线和连字符
    return "".join(char if char.isalnum() or char in ['_', '-'] else '_' for char in str(name_component))

def smiles_to_3d_sdf(smiles_string, sdf_file_path, pdb_id_for_log, seqid_for_log):
    """
    将SMILES字符串转换为带有3D坐标的SDF文件。
    成功返回True，失败返回False。
    """
    log_extra = {'pdb_id': pdb_id_for_log, 'seqid': seqid_for_log, 'smiles': smiles_string}
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            reason = f"无法从SMILES解析分子: '{smiles_string}'"
            logger.info(reason, extra={**log_extra, 'reason': "SMILES解析失败"})
            return False

        mol = Chem.AddHs(mol)  # 添加氢原子，这对于3D结构生成很重要

        # 生成3D坐标
        # ETKDG是一种较新的、通常效果更好的构象生成算法
        status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) 
        if status == -1:
            # 如果ETKDG失败 (例如，对于非常大或柔性的分子)，尝试使用随机坐标作为后备
            status_random = AllChem.EmbedMolecule(mol, useRandomCoords=True)
            if status_random == -1:
                reason = f"使用ETKDGv3和随机坐标为SMILES '{smiles_string}'生成3D坐标均失败"
                logger.info(reason, extra={**log_extra, 'reason': "3D坐标生成失败"})
                return False
        
        # 可选：优化分子几何结构 (对于大量分子可能较慢)
        # try:
        #    AllChem.UFFOptimizeMolecule(mol) # 或者 AllChem.MMFFOptimizeMolecule(mol)
        # except Exception as e:
        #    print(f"警告: 优化分子几何结构失败 for {sdf_file_path}: {e}")

        # 写入SDF文件
        writer = Chem.SDWriter(sdf_file_path)
        writer.write(mol)
        writer.close()
        return True
    except Exception as e:
        reason = f"处理SMILES '{smiles_string}' 到SDF ({sdf_file_path}) 时发生异常: {e}"
        logger.info(reason, extra={**log_extra, 'reason': f"异常: {str(e)}"})
        return False

# --- 主逻辑 ---
def main():
    # 确保日志文件在开始时是空的（如果使用 'w' filemode 则不需要这行）
    # if os.path.exists(FAILURE_LOG_FILE):
    #     os.remove(FAILURE_LOG_FILE)
    
    print(f"失败详情将记录在: {FAILURE_LOG_FILE}")

    # 检查基本PDB目录是否存在
    if not os.path.isdir(TRAIN_PDB_BASE_DIR):
        print(f"错误: 基础PDB目录 '{TRAIN_PDB_BASE_DIR}' 未找到。请确保脚本与此目录在同一级别，或者修改路径。")
        return

    # 检查CSV文件是否存在
    if not os.path.isfile(WITH_PDBID_CSV_PATH):
        print(f"错误: CSV文件 '{WITH_PDBID_CSV_PATH}' 未找到。请确保它与脚本在同一目录，或者修改路径。")
        return

    try:
        print(f"正在加载CSV文件: {WITH_PDBID_CSV_PATH}...")
        # 指定列的数据类型以避免潜在问题，特别是ID列
        df_with_pdbid = pd.read_csv(
            WITH_PDBID_CSV_PATH,
            dtype={PDB_ID_COL: str, SEQID_COL: str, SMILES_COL: str},
            on_bad_lines='warn' # 改为 'warn' 以便在控制台看到格式错误的行
        )
        print(f"从CSV加载了 {len(df_with_pdbid)} 条记录。")

        # 标准化CSV中的PDB_ID列以便匹配 (例如，转换为大写并去除首尾空格)
        df_with_pdbid[PDB_ID_COL] = df_with_pdbid[PDB_ID_COL].str.upper().str.strip()
        # 清理SMILES和seqid列中的潜在NaN值或非字符串值
        df_with_pdbid[SMILES_COL] = df_with_pdbid[SMILES_COL].astype(str).str.strip()
        df_with_pdbid[SEQID_COL] = df_with_pdbid[SEQID_COL].astype(str).str.strip()


    except Exception as e:
        print(f"错误: 加载或预处理CSV文件 '{WITH_PDBID_CSV_PATH}' 失败: {e}")
        return

    # 获取train_pdb目录下的所有PDB ID文件夹名称
    try:
        pdb_id_folders = [
            folder_name for folder_name in os.listdir(TRAIN_PDB_BASE_DIR)
            if os.path.isdir(os.path.join(TRAIN_PDB_BASE_DIR, folder_name))
        ]
    except FileNotFoundError:
        print(f"错误: 无法访问目录 '{TRAIN_PDB_BASE_DIR}'.")
        return
        
    if not pdb_id_folders:
        print(f"在 '{TRAIN_PDB_BASE_DIR}' 中没有找到PDB ID子文件夹。")
        return

    print(f"在 '{TRAIN_PDB_BASE_DIR}' 中找到 {len(pdb_id_folders)} 个PDB ID文件夹。开始处理配体...")

    successful_sdf_generations = 0
    failed_sdf_generations = 0
    skipped_existing_sdf = 0

    # 使用tqdm创建进度条
    for protein_folder_name in tqdm(pdb_id_folders, desc="处理PDB ID"):
        # protein_folder_name 是原始的文件夹名，例如 '1abc'
        # PDB ID用于在CSV中查找，应标准化 (例如，大写)
        pdb_id_for_lookup = protein_folder_name.upper().strip()

        # 在DataFrame中筛选当前PDB ID的配体记录
        ligand_entries_for_pdb = df_with_pdbid[df_with_pdbid[PDB_ID_COL] == pdb_id_for_lookup]

        if ligand_entries_for_pdb.empty:
            # print(f"信息: PDB ID {protein_folder_name} 在CSV中没有找到配体记录。")
            continue
        
        for _, ligand_row in ligand_entries_for_pdb.iterrows():
            smiles_str = ligand_row[SMILES_COL]
            seq_id_str = ligand_row[SEQID_COL]

            current_pdb_id_for_log = protein_folder_name # 原始文件夹名，更具信息性
            current_seqid_for_log = seq_id_str
            log_extra_main = {'pdb_id': current_pdb_id_for_log, 'seqid': current_seqid_for_log, 'smiles': str(smiles_str)}

            if not smiles_str or pd.isna(smiles_str) or smiles_str.lower() == 'nan' or smiles_str == '':
                reason = "SMILES为空或无效 (nan/'')"
                logger.info(reason, extra={**log_extra_main, 'reason': reason})
                failed_sdf_generations += 1
                continue
            
            if not seq_id_str or pd.isna(seq_id_str) or seq_id_str.lower() == 'nan' or seq_id_str == '':
                reason = f"seqid为空或无效 (nan/'')，无法创建唯一文件夹名"
                logger.info(reason, extra={**log_extra_main, 'reason': reason})
                failed_sdf_generations += 1
                continue

            # 清理seqid，使其适合用作文件名/目录名的一部分
            sanitized_seq_id = sanitize_name_component(seq_id_str)
            
            # 构建配体子文件夹名和路径
            # 格式: PDBID_SEQID (使用原始PDB文件夹名保持一致性)
            ligand_subfolder_name = f"{protein_folder_name}_{sanitized_seq_id}"
            ligand_subfolder_path = os.path.join(TRAIN_PDB_BASE_DIR, protein_folder_name, ligand_subfolder_name)

            # 创建配体子文件夹 (如果尚不存在)
            os.makedirs(ligand_subfolder_path, exist_ok=True)

            # 构建SDF文件名和完整路径
            # 格式: PDBID_SEQID.sdf
            sdf_file_name = f"{ligand_subfolder_name}.sdf"
            sdf_file_path = os.path.join(ligand_subfolder_path, sdf_file_name)

            # 检查SDF文件是否已存在
            if os.path.exists(sdf_file_path):
                # print(f"信息: SDF文件 {sdf_file_path} 已存在。跳过生成。")
                skipped_existing_sdf += 1
                successful_sdf_generations +=1 # 如果已存在，也算作成功处理
                continue
            
            # 生成SDF文件
            if smiles_to_3d_sdf(smiles_str, sdf_file_path, current_pdb_id_for_log, current_seqid_for_log):
                successful_sdf_generations += 1
            else:
                # 错误信息已在smiles_to_3d_sdf函数中打印
                failed_sdf_generations += 1

    print("\n--- 处理完成 ---")
    print(f"成功生成/已存在的SDF文件: {successful_sdf_generations}")
    print(f"跳过的已存在SDF文件: {skipped_existing_sdf}")
    print(f"生成失败的SDF文件 (详情见 {FAILURE_LOG_FILE}): {failed_sdf_generations}")

if __name__ == '__main__':
    # 确保脚本从`3D_structure`目录运行，或者调整相对路径
    # 例如，如果脚本在 3D_structure/scripts/ 中，则 TRAIN_PDB_BASE_DIR 和 WITH_PDBID_CSV_PATH 需要调整为 ../train_pdb 和 ../with_pdbid.csv
    main() 