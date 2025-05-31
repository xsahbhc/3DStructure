import pandas as pd
import os

def get_failed_pdb_ids(failed_file_path):
    """从failed_downloads.txt中提取失败的PDB ID。"""
    failed_pdb_ids = set()
    try:
        with open(failed_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("❌"):
                    try:
                        # 格式: ❌ PDBID: 错误信息
                        pdb_id_with_chain = line.split(':')[0].split(' ')[1].strip()
                        failed_pdb_ids.add(pdb_id_with_chain)
                    except IndexError:
                        print(f"警告: 无法解析失败记录行: {line.strip()}")
    except FileNotFoundError:
        print(f"错误: 文件 {failed_file_path} 未找到。")
    return failed_pdb_ids

def get_seqids_from_failed_pdb_ids(with_pdbid_file_path, failed_pdb_ids):
    """从with_pdbid.csv中根据失败的PDB ID找到对应的seqid。"""
    seqids_from_failed = set()
    try:
        # 尝试自动检测分隔符
        try:
            df_with_pdbid = pd.read_csv(with_pdbid_file_path, usecols=['seqid', 'PDB_ID'])
        except Exception:
            print(f"尝试使用分号作为分隔符读取 {with_pdbid_file_path}")
            df_with_pdbid = pd.read_csv(with_pdbid_file_path, usecols=['seqid', 'PDB_ID'], sep=';')
        
        # 确保列名存在
        if 'PDB_ID' not in df_with_pdbid.columns or 'seqid' not in df_with_pdbid.columns:
            print(f"错误: {with_pdbid_file_path} 文件中缺少 'PDB_ID' 或 'seqid' 列。")
            return seqids_from_failed

        # PDB ID 可能需要处理，例如去除空格或统一大小写
        df_with_pdbid['PDB_ID'] = df_with_pdbid['PDB_ID'].astype(str).str.strip()
        
        # 筛选出失败PDB ID对应的行
        filtered_df = df_with_pdbid[df_with_pdbid['PDB_ID'].isin(failed_pdb_ids)]
        seqids_from_failed.update(filtered_df['seqid'].astype(str).str.strip().tolist())
        
        print(f"从 {with_pdbid_file_path} 中为 {len(failed_pdb_ids)} 个失败的PDB ID找到了 {len(seqids_from_failed)} 个对应的seqid。")
        if len(failed_pdb_ids) > 0 and len(seqids_from_failed) == 0:
            print(f"警告: 未能从 {with_pdbid_file_path} 中为任何失败的PDB ID找到对应的seqid。请检查PDB_ID列名和内容是否匹配。")
            print(f"文件中的PDB_ID示例 (前5个非空值): {df_with_pdbid['PDB_ID'].dropna().unique()[:5]}")
            print(f"失败的PDB_ID示例 (前5个): {list(failed_pdb_ids)[:5]}")


    except FileNotFoundError:
        print(f"错误: 文件 {with_pdbid_file_path} 未找到。")
    except Exception as e:
        print(f"处理 {with_pdbid_file_path} 时发生错误: {e}")
    return seqids_from_failed

def get_seqids_from_without_pdbid(without_pdbid_file_path):
    """从without_pdbid.csv中提取所有seqid。"""
    seqids_without = set()
    try:
        # 尝试自动检测分隔符
        try:
            df_without_pdbid = pd.read_csv(without_pdbid_file_path, usecols=['seqid'])
        except Exception:
            print(f"尝试使用分号作为分隔符读取 {without_pdbid_file_path}")
            df_without_pdbid = pd.read_csv(without_pdbid_file_path, usecols=['seqid'], sep=';')

        if 'seqid' not in df_without_pdbid.columns:
            print(f"错误: {without_pdbid_file_path} 文件中缺少 'seqid' 列。")
            return seqids_without
            
        seqids_without.update(df_without_pdbid['seqid'].astype(str).str.strip().tolist())
        print(f"从 {without_pdbid_file_path} 中找到了 {len(seqids_without)} 个seqid。")
    except FileNotFoundError:
        print(f"错误: 文件 {without_pdbid_file_path} 未找到。")
    except Exception as e:
        print(f"处理 {without_pdbid_file_path} 时发生错误: {e}")
    return seqids_without

def count_seqids_in_train_csv(train_csv_path, combined_seqids):
    """分块读取train.csv并统计combined_seqids中的seqid出现的次数。"""
    if not combined_seqids:
        print("没有需要统计的seqid。")
        return 0, 0

    total_matched_rows = 0
    found_seqids_in_train = set()
    chunk_size = 100000  # 根据内存调整块大小
    
    try:
        print(f"开始分块处理 {train_csv_path}...")
        iterator = pd.read_csv(train_csv_path, usecols=['seqid'], chunksize=chunk_size, low_memory=False)
        
        for i, chunk in enumerate(iterator):
            print(f"正在处理块 {i+1}...")
            if 'seqid' not in chunk.columns:
                print(f"错误: {train_csv_path} 的块中缺少 'seqid' 列。")
                return total_matched_rows, len(found_seqids_in_train)

            chunk['seqid'] = chunk['seqid'].astype(str).str.strip()
            
            # 筛选出当前块中与combined_seqids匹配的行
            matched_chunk = chunk[chunk['seqid'].isin(combined_seqids)]
            total_matched_rows += len(matched_chunk)
            found_seqids_in_train.update(matched_chunk['seqid'].tolist())
        
        print("处理完成。")
        return total_matched_rows, len(found_seqids_in_train)
        
    except FileNotFoundError:
        print(f"错误: 文件 {train_csv_path} 未找到。")
    except Exception as e:
        print(f"处理 {train_csv_path} 时发生错误: {e}")
    return 0, 0


if __name__ == "__main__":
    # 定义文件路径 (假设脚本与这些文件在同一目录下)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    failed_file = os.path.join(base_dir, "failed_downloads.txt")
    with_pdbid_file = os.path.join(base_dir, "with_pdbid.csv")
    without_pdbid_file = os.path.join(base_dir, "without_pdbid.csv")
    train_csv_file = os.path.join(base_dir, "train.csv")

    print("开始分析...")
    
    # 1. 从failed_downloads.txt获取失败的PDB ID
    failed_pdb_ids = get_failed_pdb_ids(failed_file)
    print(f"从 {failed_file} 中提取了 {len(failed_pdb_ids)} 个唯一的下载失败PDB ID。")
    if not failed_pdb_ids and os.path.exists(failed_file):
        print(f"注意: {failed_file} 中没有找到格式为 '❌ PDBID: ...' 的失败记录。")

    # 2. 从with_pdbid.csv找到失败PDB ID对应的seqid
    seqids_from_failed = get_seqids_from_failed_pdb_ids(with_pdbid_file, failed_pdb_ids)

    # 3. 从without_pdbid.csv获取seqid
    seqids_without = get_seqids_from_without_pdbid(without_pdbid_file)
    
    # 4. 合并所有目标seqid
    combined_seqids_to_check = seqids_from_failed.union(seqids_without)
    print(f"总共有 {len(combined_seqids_to_check)} 个唯一的seqid需要检查 (来自下载失败和无PDB ID的列表)。")

    # 5. 在train.csv中统计这些seqid
    if combined_seqids_to_check:
        total_rows, unique_seqids_found = count_seqids_in_train_csv(train_csv_file, combined_seqids_to_check)
        print("\n--- 统计结果 ---")
        print(f"在 {train_csv_file} 中:" )
        print(f"  - 找到了 {unique_seqids_found} 个来自合并列表的唯一seqid。")
        print(f"  - 这些seqid总共对应了 {total_rows} 条数据。")
    else:
        print("\n没有从任何来源获得需要检查的seqid，无法进行统计。")

    print("\n分析完成。") 