import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import logging
import multiprocessing as mp
from functools import partial
import time
import threading

# --- 配置 ---
TRAIN_PDB_BASE_DIR = 'biosnap/train_pdb'
WITH_PDBID_CSV_PATH = 'biosnap/train_csv/with_pdbid.csv'

# CSV文件中的列名
PDB_ID_COL = 'PDB_ID'
SMILES_COL = 'SMILES'
SEQID_COL = 'seqid'
FAILURE_LOG_FILE = 'biosnap/sdf_generation_failures_fast.log'

# 性能优化配置
MAX_WORKERS = 1  # 使用单进程避免卡住
BATCH_SIZE = 100  # 批处理大小
SKIP_OPTIMIZATION = True  # 跳过几何优化，与之前数据集保持一致
TIMEOUT_SECONDS = 30  # 单个分子处理超时时间（秒）

# 设置日志记录
logging.basicConfig(filename=FAILURE_LOG_FILE, level=logging.INFO, 
                    format='%(asctime)s - PDB: %(pdb_id)s - SEQID: %(seqid)s - SMILES: %(smiles)s - REASON: %(reason)s',
                    filemode='w')
logger = logging.getLogger(__name__)

def run_with_timeout(func, args, timeout_seconds):
    """在指定时间内运行函数，超时则返回None"""
    result = [None]
    exception = [None]

    def target():
        try:
            result[0] = func(*args)
        except Exception as e:
            exception[0] = e

    thread = threading.Thread(target=target)
    thread.daemon = True
    thread.start()
    thread.join(timeout_seconds)

    if thread.is_alive():
        # 超时了
        return None, "timeout"
    elif exception[0]:
        return None, str(exception[0])
    else:
        return result[0], None

def sanitize_name_component(name_component):
    """对用于文件名或目录名的部分进行基本清理。"""
    return "".join(char if char.isalnum() or char in ['_', '-'] else '_' for char in str(name_component))

def smiles_to_3d_sdf_fast(smiles_string, sdf_file_path, pdb_id_for_log, seqid_for_log):
    """
    将SMILES字符串转换为带有3D坐标的SDF文件。
    保持与之前数据集一致的生成方法，确保质量。
    """
    log_extra = {'pdb_id': pdb_id_for_log, 'seqid': seqid_for_log, 'smiles': smiles_string}
    try:
        # 禁用RDKit的警告信息以减少输出
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')

        print(f"        🔄 开始处理SMILES: {smiles_string[:50]}...")

        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            reason = f"无法从SMILES解析分子: '{smiles_string}'"
            logger.info(reason, extra={**log_extra, 'reason': "SMILES解析失败"})
            print(f"        ❌ SMILES解析失败")
            return False

        print(f"        ➕ 添加氢原子...")
        mol = Chem.AddHs(mol)

        print(f"        🎯 生成3D构象...")
        # 使用超时机制的3D构象生成
        def generate_3d_conformation():
            params = AllChem.ETKDGv3()
            params.maxAttempts = 5  # 保持合理的尝试次数
            params.numThreads = 1   # 单线程避免冲突
            params.randomSeed = 42  # 固定随机种子保证可重现性

            status = AllChem.EmbedMolecule(mol, params)
            if status == -1:
                # 后备方案：使用随机坐标
                status_random = AllChem.EmbedMolecule(mol, useRandomCoords=True)
                return status_random
            return status

        # 使用30秒超时
        result, error = run_with_timeout(generate_3d_conformation, (), TIMEOUT_SECONDS)

        if error == "timeout":
            reason = f"3D构象生成超时 ({TIMEOUT_SECONDS}秒)"
            logger.info(reason, extra={**log_extra, 'reason': reason})
            print(f"        ⏰ 3D构象生成超时，跳过")
            return False
        elif error:
            reason = f"3D构象生成异常: {error}"
            logger.info(reason, extra={**log_extra, 'reason': reason})
            print(f"        ❌ 3D构象生成异常: {error}")
            return False
        elif result == -1:
            reason = f"3D坐标生成失败"
            logger.info(reason, extra={**log_extra, 'reason': "3D坐标生成失败"})
            print(f"        ❌ 3D坐标生成失败")
            return False

        # 进行几何优化（保持与之前数据集一致）
        if not SKIP_OPTIMIZATION:
            print(f"        ⚙️ 进行几何优化...")

            def optimize_geometry():
                return AllChem.UFFOptimizeMolecule(mol, maxIters=200)

            # 使用超时机制进行几何优化
            opt_result, opt_error = run_with_timeout(optimize_geometry, (), TIMEOUT_SECONDS)

            if opt_error == "timeout":
                print(f"        ⏰ 几何优化超时，跳过优化步骤")
                # 超时不影响整体流程，继续使用未优化的结构
            elif opt_error:
                print(f"        ⚠️ 几何优化异常: {opt_error}")
                # 优化失败不影响整体流程
            else:
                print(f"        ✅ 几何优化完成")
                # 优化成功

        print(f"        💾 写入SDF文件...")
        # 写入SDF文件
        writer = Chem.SDWriter(sdf_file_path)
        writer.write(mol)
        writer.close()
        print(f"        ✅ SDF文件生成成功")
        return True

    except Exception as e:
        reason = f"处理异常: {str(e)}"
        logger.info(reason, extra={**log_extra, 'reason': reason})
        print(f"        ❌ 处理异常: {str(e)}")
        return False

def process_ligand_batch(ligand_batch, batch_num):
    """处理一批配体数据"""
    results = {
        'success': 0,
        'failed': 0,
        'skipped': 0
    }

    for i, ligand_data in enumerate(ligand_batch):
        pdb_id, seqid, smiles, sdf_path = ligand_data

        # 显示当前处理的文件
        print(f"    [{batch_num}:{i+1}] 处理 {pdb_id}_{seqid}")

        # 检查文件是否已存在
        if os.path.exists(sdf_path):
            print(f"        ⏩ 跳过 (文件已存在): {sdf_path}")
            results['skipped'] += 1
            continue

        # 确保目录存在
        os.makedirs(os.path.dirname(sdf_path), exist_ok=True)

        # 生成SDF文件
        if smiles_to_3d_sdf_fast(smiles, sdf_path, pdb_id, seqid):
            print(f"        ✅ 成功生成: {sdf_path}")
            results['success'] += 1
        else:
            print(f"        ❌ 生成失败: {pdb_id}_{seqid} - SMILES: {smiles[:50]}...")
            results['failed'] += 1

    return results

def main():
    print(f"🚀 启动快速SDF生成脚本...")
    print(f"📊 配置: 最大进程数={MAX_WORKERS}, 批处理大小={BATCH_SIZE}")
    print(f"⚡ 性能优化: 跳过几何优化={SKIP_OPTIMIZATION}")
    print(f"📝 失败日志: {FAILURE_LOG_FILE}")
    
    # 检查基本目录
    if not os.path.isdir(TRAIN_PDB_BASE_DIR):
        print(f"❌ 错误: 基础PDB目录 '{TRAIN_PDB_BASE_DIR}' 未找到")
        return
    
    if not os.path.isfile(WITH_PDBID_CSV_PATH):
        print(f"❌ 错误: CSV文件 '{WITH_PDBID_CSV_PATH}' 未找到")
        return

    try:
        print(f"📖 正在加载CSV文件: {WITH_PDBID_CSV_PATH}...")
        df_with_pdbid = pd.read_csv(
            WITH_PDBID_CSV_PATH,
            dtype={PDB_ID_COL: str, SEQID_COL: str, SMILES_COL: str}
        )
        print(f"✅ 从CSV加载了 {len(df_with_pdbid)} 条记录")
        
        # 数据清理
        df_with_pdbid[PDB_ID_COL] = df_with_pdbid[PDB_ID_COL].str.upper().str.strip()
        df_with_pdbid[SMILES_COL] = df_with_pdbid[SMILES_COL].astype(str).str.strip()
        df_with_pdbid[SEQID_COL] = df_with_pdbid[SEQID_COL].astype(str).str.strip()
        
    except Exception as e:
        print(f"❌ 错误: 加载CSV文件失败: {e}")
        return

    # 获取PDB文件夹
    try:
        pdb_id_folders = [
            folder for folder in os.listdir(TRAIN_PDB_BASE_DIR)
            if os.path.isdir(os.path.join(TRAIN_PDB_BASE_DIR, folder))
        ]
        print(f"📁 找到 {len(pdb_id_folders)} 个PDB文件夹")
    except Exception as e:
        print(f"❌ 错误: 读取PDB目录失败: {e}")
        return

    # 准备所有需要处理的配体数据
    all_ligand_data = []
    
    print("🔍 准备配体数据...")
    for protein_folder_name in tqdm(pdb_id_folders, desc="扫描PDB文件夹"):
        pdb_id_for_lookup = protein_folder_name.upper().strip()
        ligand_entries = df_with_pdbid[df_with_pdbid[PDB_ID_COL] == pdb_id_for_lookup]
        
        if ligand_entries.empty:
            continue
        
        for _, ligand_row in ligand_entries.iterrows():
            smiles_str = ligand_row[SMILES_COL]
            seq_id_str = ligand_row[SEQID_COL]
            
            if pd.isna(smiles_str) or smiles_str.strip() == '' or smiles_str == 'nan':
                continue
            
            # 构建文件路径 - 使用与BindingDB一致的命名格式 PDBID_SEQID
            sanitized_seqid = sanitize_name_component(seq_id_str)
            ligand_folder_name = f"{protein_folder_name}_{sanitized_seqid}"  # PDBID_SEQID格式
            ligand_folder_path = os.path.join(TRAIN_PDB_BASE_DIR, protein_folder_name, ligand_folder_name)
            sdf_file_path = os.path.join(ligand_folder_path, f"{ligand_folder_name}.sdf")
            
            all_ligand_data.append((pdb_id_for_lookup, seq_id_str, smiles_str, sdf_file_path))
    
    print(f"📊 总共需要处理 {len(all_ligand_data)} 个配体")
    
    if not all_ligand_data:
        print("⚠️  没有找到需要处理的配体数据")
        return
    
    # 分批处理
    batches = [all_ligand_data[i:i + BATCH_SIZE] for i in range(0, len(all_ligand_data), BATCH_SIZE)]
    print(f"📦 分为 {len(batches)} 个批次处理")
    
    # 统计结果
    total_success = 0
    total_failed = 0
    total_skipped = 0
    
    start_time = time.time()
    
    # 使用单进程处理 (避免tmux中的多进程问题)
    print("🔄 使用单进程处理 (跳过几何优化，与之前数据集保持一致)...")
    results = []

    # 使用tqdm进度条
    with tqdm(total=len(batches), desc="处理批次", unit="batch") as pbar:
        for i, batch in enumerate(batches):
            # 更新进度条描述
            pbar.set_description(f"处理批次 {i+1}/{len(batches)}")

            result = process_ligand_batch(batch, i+1)
            results.append(result)

            # 更新进度条后缀信息
            temp_success = sum(r['success'] for r in results)
            temp_failed = sum(r['failed'] for r in results)
            temp_skipped = sum(r['skipped'] for r in results)
            elapsed = time.time() - start_time
            rate = (temp_success + temp_failed) / elapsed if elapsed > 0 else 0

            pbar.set_postfix({
                '成功': temp_success,
                '失败': temp_failed,
                '跳过': temp_skipped,
                '速度': f"{rate:.1f}/秒"
            })

            # 更新进度条
            pbar.update(1)

            # 每处理20个批次显示详细统计
            if (i + 1) % 20 == 0:
                print(f"\n📊 阶段统计 (批次{i+1}): 成功={temp_success}, 失败={temp_failed}, 跳过={temp_skipped}, 平均速度={rate:.2f}/秒")
    
    # 汇总结果
    for result in results:
        total_success += result['success']
        total_failed += result['failed']
        total_skipped += result['skipped']
    
    elapsed_time = time.time() - start_time
    
    print(f"\n🎉 处理完成!")
    print(f"⏱️  总耗时: {elapsed_time:.1f} 秒")
    print(f"✅ 成功生成: {total_success} 个SDF文件")
    print(f"⏩ 跳过已存在: {total_skipped} 个文件")
    print(f"❌ 生成失败: {total_failed} 个文件")
    print(f"📝 详细失败信息请查看: {FAILURE_LOG_FILE}")
    
    if total_success > 0:
        avg_time = elapsed_time / (total_success + total_failed)
        print(f"📊 平均处理速度: {avg_time:.2f} 秒/个")

if __name__ == '__main__':
    main()
