import pandas as pd
import subprocess
import os
import platform

print("开始运行脚本...")

# 读取 CSV 文件
try:
    print("正在读取 biosnap train.csv 文件...")
    df = pd.read_csv('biosnap/train_csv/train.csv')
    print(f"成功读取 {len(df)} 条序列数据")
except FileNotFoundError:
    print("错误：biosnap/train_csv/train.csv 文件不存在，请检查路径。")
    exit(1)

print("正在创建输出文件...")
# 确保biosnap目录存在
os.makedirs('biosnap', exist_ok=True)
# 创建输出文件
with open('biosnap/with_pdbid.csv', 'w') as f_with, open('biosnap/without_pdbid.csv', 'w') as f_without:
    # 写入表头（新增 seqid）
    f_with.write("seqid,Protein,PDB_ID,SMILES\n")
    f_without.write("seqid,Protein,SMILES\n")

    # 将所有序列写入临时 FASTA 文件
    try:
        print("正在创建临时 FASTA 文件...")
        with open('biosnap/temp.fasta', 'w') as temp_fasta:
            for idx, row in df.iterrows():
                protein_seq = row['Protein']
                temp_fasta.write(f">{idx}\n{protein_seq}\n")
        print("临时 FASTA 文件创建完成")
    except Exception as e:
        print("写入 FASTA 文件失败:", e)
        exit(1)

    # 使用 DIAMOND 进行批量查询
    print("开始运行 DIAMOND 比对...")
    diamond_cmd = [
        os.path.expanduser('~/bin/diamond'), 'blastp',
        '--query', 'biosnap/temp.fasta',
        '--db', 'pdb_db.dmnd',
        '--outfmt', '6',  # 使用制表符分隔的格式
        '--evalue', '1e-5',
        '--max-target-seqs', '1'
    ]
    try:
        result = subprocess.run(diamond_cmd, capture_output=True, text=True, check=True)
        print("DIAMOND 比对完成")
    except subprocess.CalledProcessError as e:
        print("DIAMOND 执行失败:", e.stderr)
        exit(1)

    # 解析 DIAMOND 输出
    print("正在处理比对结果...")
    if result.stdout:
        lines = result.stdout.strip().split('\n')
        results_dict = {}
        for line in lines:
            parts = line.split('\t')  # 确保输出格式正确
            if len(parts) < 2:
                continue
            query_id = parts[0]
            pdb_id = parts[1]
            results_dict[query_id] = pdb_id

        # 分类输出结果
        print("正在写入结果文件...")
        match_count = 0
        no_match_count = 0
        for idx, row in df.iterrows():
            seqid = row['seqid']         # 新增：获取 seqid 字段
            protein_seq = row['Protein']
            smiles = row['SMILES']
            if str(idx) in results_dict:
                f_with.write(f"{seqid},{protein_seq},{results_dict[str(idx)]},{smiles}\n")
                match_count += 1
            else:
                f_without.write(f"{seqid},{protein_seq},{smiles}\n")
                no_match_count += 1
        print(f"处理完成！找到PDB匹配的序列：{match_count}个，未找到匹配的序列：{no_match_count}个")
    else:
        # 如果无输出，全部写入无 PDB ID 文件
        print("未找到任何匹配结果，所有序列将写入without_pdbid.csv")
        for _, row in df.iterrows():
            seqid = row['seqid']
            protein_seq = row['Protein']
            smiles = row['SMILES']
            f_without.write(f"{seqid},{protein_seq},{smiles}\n")


# 清理临时文件
print("清理临时文件...")
if os.path.exists('biosnap/temp.fasta'):
    if platform.system() != 'Windows':
        os.system('rm biosnap/temp.fasta')
    else:
        os.system('del biosnap\\temp.fasta')

print("脚本执行完成！")