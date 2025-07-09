import os
import pandas as pd
import asyncio
import aiohttp
from tqdm import tqdm
import time

# 全局配置
MAX_CONCURRENT = 20  # 降低并发数
MAX_RETRIES = 3     # 最大重试次数
TIMEOUT = 60        # 增加超时时间
RETRY_DELAY = 1     # 重试间隔时间（秒）

async def download_pdb(session, pdb_id, full_pdb_id, pdb_dir, semaphore):
    # 检查文件是否已存在
    pdb_path = os.path.join(pdb_dir, f"{full_pdb_id}.pdb")
    if os.path.exists(pdb_path):
        return f"⏩ {full_pdb_id}: 文件已存在"

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    for retry in range(MAX_RETRIES):
        try:
            async with semaphore:
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=TIMEOUT)) as response:
                    if response.status == 200:
                        content = await response.read()
                        os.makedirs(pdb_dir, exist_ok=True)
                        with open(pdb_path, 'wb') as f:
                            f.write(content)
                        return f"✅ {full_pdb_id}: 下载成功"
                    elif response.status == 404:
                        return f"❌ {full_pdb_id}: PDB文件不存在 (404)"
                    else:
                        if retry < MAX_RETRIES - 1:
                            await asyncio.sleep(RETRY_DELAY)
                            continue # 进入下一次重试
                        return f"❌ {full_pdb_id}: HTTP {response.status} (已重试 {retry + 1} 次)"
        except asyncio.TimeoutError:
            if retry < MAX_RETRIES - 1:
                await asyncio.sleep(RETRY_DELAY)
                continue # 进入下一次重试
            return f"❌ {full_pdb_id}: 下载超时 (已重试 {retry + 1} 次)"
        except Exception as e:
            if retry < MAX_RETRIES - 1:
                await asyncio.sleep(RETRY_DELAY)
                continue # 进入下一次重试
            return f"❌ {full_pdb_id}: {str(e)} (已重试 {retry + 1} 次)"
    
    return f"❌ {full_pdb_id}: 所有 {MAX_RETRIES} 次重试均失败"

async def main():
    start_time = time.time()
    
    # 设置主目录
    main_dir = 'biosnap/train_pdb'
    os.makedirs(main_dir, exist_ok=True)

    # 读取CSV文件
    try:
        df = pd.read_csv('biosnap/train_csv/with_pdbid.csv')
        if 'PDB_ID' not in df.columns:
            print("❌ 错误：CSV文件中找不到 'PDB_ID' 列")
            return
    except FileNotFoundError:
        print("❌ 错误：找不到 biosnap/train_csv/with_pdbid.csv 文件")
        return
    except Exception as e:
        print(f"❌ 错误：读取CSV文件时出错 - {str(e)}")
        return

    # 检查CSV中的唯一PDB ID数量
    unique_pdb_ids = df['PDB_ID'].nunique()
    total_rows = len(df)
    print(f"ℹ️ CSV文件总行数: {total_rows}")
    print(f"ℹ️ CSV文件中唯一的PDB ID数量: {unique_pdb_ids}")
    if total_rows > unique_pdb_ids:
        print(f"⚠️ 警告: CSV文件中包含 {total_rows - unique_pdb_ids} 个重复的PDB ID条目。")

    print(f"📋 总共需要处理 {total_rows} 个PDB条目")
    
    semaphore = asyncio.Semaphore(MAX_CONCURRENT)
    connector = aiohttp.TCPConnector(limit=MAX_CONCURRENT, force_close=False, enable_cleanup_closed=True)
    client_timeout = aiohttp.ClientTimeout(total=TIMEOUT)
    
    async with aiohttp.ClientSession(connector=connector, timeout=client_timeout) as session:
        tasks = []
        for index, row in df.iterrows():
            # 去除PDB_ID前后可能存在的空格，并确保是字符串类型
            full_pdb_id = str(row['PDB_ID']).strip()
            
            if not full_pdb_id: # 跳过空的PDB_ID
                print(f"⚠️ 警告: 第 {index + 1} 行的 PDB_ID 为空，跳过")
                continue

            if len(full_pdb_id) < 4 or not full_pdb_id[:4].isalnum():
                print(f"⚠️ 警告: 第 {index + 1} 行的 PDB_ID '{full_pdb_id}' 格式不正确，跳过")
                continue
            
            pdb_id = full_pdb_id[:4].lower()
            pdb_dir = os.path.join(main_dir, full_pdb_id) 
            
            task = download_pdb(session, pdb_id, full_pdb_id, pdb_dir, semaphore)
            tasks.append(task)
        
        print("🚀 开始下载...")
        results = []
        failed_downloads = []
        
        with tqdm(total=len(tasks), desc="📥 下载进度") as pbar:
            for f_completed in asyncio.as_completed(tasks):
                result = await f_completed
                results.append(result)
                if "❌" in result:
                    failed_downloads.append(result)
                pbar.update(1)
                if "❌" in result:
                     # 避免在tqdm进度条刷新时破坏输出格式
                    pbar.write(f"{result}") 
    
    success = len([r for r in results if "✅" in r])
    skipped = len([r for r in results if "⏩" in r])
    failed = len([r for r in results if "❌" in r])
    
    end_time = time.time()
    duration = end_time - start_time
    
    print("\n📊 下载统计:")
    print(f"✅ 成功: {success}")
    print(f"⏩ 跳过 (文件已存在): {skipped}")
    print(f"❌ 失败: {failed}")
    print(f"⏱️ 总耗时: {duration:.2f} 秒")
    if len(results) > 0 and duration > 0:
        print(f"🚀 平均速度: {len(results)/duration:.2f} 个条目/秒")
    
    if failed_downloads:
        with open('failed_downloads.txt', 'w') as f_err:
            f_err.write('\n'.join(failed_downloads))
        print(f"\n❗ 失败的下载已保存到 failed_downloads.txt")

if __name__ == "__main__":
    # 确保在Windows上也能正常运行asyncio的tqdm
    if os.name == 'nt':
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    asyncio.run(main())