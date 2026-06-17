#!/usr/bin/env python3
"""
SRA智能下载工具 v4.7 - 不死鸟版 (双协议智能降级与防封控机制)
===============================================
特色升级：
  1. 解析层：打通 INSDC(国际同步) 与 CRA(国内自有) 网页，完美处理 HTML 边界。
  2. 传输层：[新增] 智能协议降级！FTP 被拒瞬间切 HTTP，永不中断。
  3. 防护层：[新增] 降低 aria2c 并发至 4 线程，防止触发 NGDC 服务器 IP 封控。
  4. 存储层：动态捕获真实文件后缀，断点续传完美支持。
"""

import os
import sys
import re
import csv

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import time
import argparse
import subprocess
import shutil
import threading
import random
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import wraps

import requests
from tqdm import tqdm

# ==================== 全局配置 ====================
DEFAULT_NGDC_METHOD = "aria2c"
DEFAULT_NGDC_CONCURRENCY = 1
DEFAULT_PREFETCH_CONCURRENCY = 4
DEFAULT_OUTPUT = "./sra_data"
TIMEOUT = 30
CHUNK_SIZE = 8192
RETRY_TIMES = 3                
RETRY_DELAY = 1               
TASK_DELAY = 0.5             

# 线程安全缓存与会话
_url_cache = {}
_cache_lock = threading.Lock()
_session = requests.Session()
_session.headers.update({
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'Connection': 'keep-alive',
})

class ProgressTracker:
    def __init__(self, total):
        self.total = total
        self.current = 0
        self.lock = threading.Lock()
        
    def increment(self):
        with self.lock:
            self.current += 1
            return self.current

def sanitize_dirname(name: str) -> str:
    if not name: return "Uncategorized"
    safe = re.sub(r'[\\/*?:"<>| ]', '_', name).strip()
    return safe if safe else "Uncategorized"

def retry(max_retries=RETRY_TIMES, delay=RETRY_DELAY):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            _delay = delay
            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    if attempt == max_retries: raise
                    sleep_time = _delay * (2 ** attempt) + random.uniform(0, 0.5)
                    print(f"⚠️ 请求异常 ({e}), 等待 {sleep_time:.1f}s 后第{attempt+1}次重试...")
                    time.sleep(sleep_time)
            return None
        return wrapper
    return decorator

def cached_get(key, func, *args, **kwargs):
    with _cache_lock:
        if key not in _url_cache:
            _url_cache[key] = func(*args, **kwargs)
        return _url_cache[key]

@retry(max_retries=RETRY_TIMES, delay=RETRY_DELAY)
def fetch_ngdc_url_from_gsa(accession):
    accession = accession.strip().upper()
    
    search_url = f"https://ngdc.cncb.ac.cn/gsa/search?searchTerm={accession}"
    resp = _session.get(search_url, timeout=TIMEOUT)
    resp.raise_for_status()

    def extract_download_urls(html_text):
        urls = {}
        ftp_pattern = rf'(ftp://download[0-9]*\.(?:cncb|big)\.ac\.cn/[^"\'<>\s]*{accession}[^"\'<>\s]*)'
        http_pattern = rf'(https?://download[0-9]*\.(?:cncb|big)\.ac\.cn/[^"\'<>\s]*{accession}[^"\'<>\s]*)'
        
        ftp_match = re.search(ftp_pattern, html_text, re.IGNORECASE)
        if ftp_match:
            urls['ftp'] = ftp_match.group(1).strip().replace('//SRR', '/SRR').replace('//ERR', '/ERR').replace('//CRR', '/CRR')
            
        http_match = re.search(http_pattern, html_text, re.IGNORECASE)
        if http_match:
            urls['http'] = http_match.group(1).strip().replace('//SRR', '/SRR').replace('//ERR', '/ERR').replace('//CRR', '/CRR')
            
        return urls if urls else None

    dl_urls = extract_download_urls(resp.text)
    if dl_urls: return dl_urls

    match = re.search(rf'href="((?:/gsa/)?browse/[^"]*{accession}[^"]*)"', resp.text)
    if not match:
        old_search_url = f"https://ngdc.cncb.ac.cn/gsa/search?db=GSA&term={accession}"
        resp_old = _session.get(old_search_url, timeout=TIMEOUT)
        match = re.search(rf'href="((?:/gsa/)?browse/[^"]*{accession}[^"]*)"', resp_old.text)
        if not match: return None

    path = match.group(1)
    if path.startswith('/gsa/'): detail_url = f"https://ngdc.cncb.ac.cn{path}"
    elif path.startswith('browse/'): detail_url = f"https://ngdc.cncb.ac.cn/gsa/{path}"
    else: detail_url = f"https://ngdc.cncb.ac.cn/gsa/{path.lstrip('/')}"

    detail_resp = _session.get(detail_url, timeout=TIMEOUT)
    detail_resp.raise_for_status()

    return extract_download_urls(detail_resp.text)

def download_with_requests(url, out_path):
    headers, mode, initial_pos = {}, 'wb', 0
    if os.path.exists(out_path):
        existing_size = os.path.getsize(out_path)
        headers['Range'] = f'bytes={existing_size}-'
        mode, initial_pos = 'ab', existing_size
        print(f"📌 检测到本地存在 {existing_size} 字节，启动断点续传...")

    try:
        with _session.get(url, headers=headers, stream=True, timeout=TIMEOUT) as r:
            if r.status_code == 416:
                print("✨ 本地文件已完整，完美跳过！")
                return True
            r.raise_for_status()
            
            total_size = int(r.headers.get('content-length', 0)) + initial_pos
            if r.status_code != 206:
                mode, initial_pos, total_size = 'wb', 0, int(r.headers.get('content-length', 0))

            with open(out_path, mode) as f, tqdm(
                total=total_size, unit='B', unit_scale=True, 
                desc=os.path.basename(out_path), initial=initial_pos, ascii=True
            ) as pbar:
                for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
        return True
    except Exception as e:
        print(f"❌ requests下载失败: {e}")
        return False

def download_with_aria2c(url, out_dir):
    # 【修复】将并发限制在 4 线程，避免触发 NGDC 服务器的反滥用(Anti-Abuse)机制
    cmd = ["aria2c", "-x", "3", "-s", "3", "-k", "1M", "-c", "--console-log-level=warn", "--summary-interval=10", "-d", out_dir, url]
    try:
        subprocess.run(cmd, check=True, timeout=3600)
        return True
    except subprocess.CalledProcessError: return False

def download_with_wget(url, out_dir):
    cmd = ["wget", "-c", "-q", "--show-progress", "-P", out_dir, url]
    try:
        subprocess.run(cmd, check=True, timeout=3600)
        return True
    except subprocess.CalledProcessError: return False

def download_single_ngdc(accession, method, base_out_dir, sample_name, tracker=None):
    try:
        out_dir = os.path.join(base_out_dir, sample_name)
        os.makedirs(out_dir, exist_ok=True)
        
        urls = cached_get(accession, fetch_ngdc_url_from_gsa, accession)
        if not urls: 
            if tracker: curr = tracker.increment()
            print(f"❌ [NGDC] 未找到相关下载链接: {accession}" + (f" (总体进度: {curr}/{tracker.total})" if tracker else ""))
            return False

        # === 【重磅更新】双协议兜底队列 ===
        protocols_to_try = []
        if method == "requests":
            if urls.get('http'): protocols_to_try.append(('HTTP', urls['http']))
        else:
            # 优先 FTP，如果 FTP 被拒绝，则备用 HTTP
            if urls.get('ftp'): protocols_to_try.append(('FTP', urls['ftp']))
            if urls.get('http'): protocols_to_try.append(('HTTP', urls['http']))

        if not protocols_to_try:
            if tracker: curr = tracker.increment()
            print(f"❌ [NGDC] 无可用协议链接: {accession}" + (f" (总体进度: {curr}/{tracker.total})" if tracker else ""))
            return False

        filename = protocols_to_try[0][1].split('/')[-1]
        out_path = os.path.join(out_dir, filename)
        success = False

        for proto_name, url in protocols_to_try:
            print(f"\n📦 [{sample_name}] 开始处理: {accession} [{proto_name}]")
            print(f"🔗 地址: {url}")
            print(f"📄 文件: {filename}")
            
            if method == "requests": success = download_with_requests(url, out_path)
            elif method == "aria2c" and shutil.which("aria2c"): success = download_with_aria2c(url, out_dir)
            elif method == "wget" and shutil.which("wget"): success = download_with_wget(url, out_dir)
            else: success = download_with_requests(url, out_path)

            if success and os.path.exists(out_path):
                break # 成功则直接跳出重试队列
            else:
                print(f"⚠️ [{proto_name}] 协议连接失败或中断，切换备用方案...")

        if success and os.path.exists(out_path):
            if tracker: curr = tracker.increment()
            print(f"✅ [NGDC] 成功完成: {accession}" + (f" (总体进度: {curr}/{tracker.total})" if tracker else ""))
            return True
        else:
            if tracker: curr = tracker.increment()
            print(f"❌ [NGDC] 所有协议下载均失败: {accession}" + (f" (总体进度: {curr}/{tracker.total})" if tracker else ""))
            return False
            
    except Exception as e:
        if tracker: curr = tracker.increment()
        print(f"❌ [NGDC] 处理 {accession} 时异常: {e}" + (f" (总体进度: {curr}/{tracker.total})" if tracker else ""))
        return False

def download_batch_ngdc(accession_list, acc_to_sample, method, base_out_dir, max_workers):
    results = {}
    tracker = ProgressTracker(len(accession_list))
    if max_workers <= 1:
        for idx, acc in enumerate(accession_list):
            results[acc] = download_single_ngdc(acc, method, base_out_dir, acc_to_sample.get(acc, "Uncategorized"), tracker)
            if idx < len(accession_list) - 1: time.sleep(TASK_DELAY)
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_acc = {
                executor.submit(download_single_ngdc, acc, method, base_out_dir, acc_to_sample.get(acc, "Uncategorized"), tracker): acc
                for acc in accession_list
            }
            for future in as_completed(future_to_acc):
                results[future_to_acc[future]] = future.result()
    return results

def download_with_prefetch(accession, base_out_dir, sample_name, tracker=None):
    out_dir = os.path.join(base_out_dir, sample_name)
    os.makedirs(out_dir, exist_ok=True)
    
    cmd = ["prefetch", "--max-size", "u", "-O", out_dir, accession]
    
    max_retries = 5 
    for attempt in range(max_retries):
        try:
            subprocess.run(cmd, check=True, timeout=86400, stdout=subprocess.DEVNULL)
            
            for root, _, files in os.walk(out_dir):
                if f"{accession}.sra" in files:
                    if tracker: curr = tracker.increment()
                    print(f"✅ [prefetch] 成功完成: {accession}" + (f" (总体进度: {curr}/{tracker.total})" if tracker else ""))
                    return True
            break
            
        except subprocess.TimeoutExpired:
            print(f"⚠️ [prefetch] {accession} 第 {attempt+1} 次下载超时 (超24小时)，准备断点续传...")
        except subprocess.CalledProcessError as e:
            print(f"⚠️ [prefetch] {accession} 网络连接断开 (尝试 {attempt+1}/{max_retries})，等待 10 秒后自动重试...")
            time.sleep(10)
        except Exception as e:
            print(f"⚠️ [prefetch] {accession} 出现未知错误: {e}")
            time.sleep(5)
            
    if tracker: curr = tracker.increment()
    print(f"❌ [prefetch] {accession} 失败: 历经多次重试仍无法完成下载" + (f" (总体进度: {curr}/{tracker.total})" if tracker else ""))
    return False

def download_batch_prefetch(accession_list, acc_to_sample, base_out_dir, max_workers):
    if not accession_list: return {}
    results = {}
    tracker = ProgressTracker(len(accession_list))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_acc = {
            executor.submit(download_with_prefetch, acc, base_out_dir, acc_to_sample.get(acc, "Uncategorized"), tracker): acc
            for acc in accession_list
        }
        for future in as_completed(future_to_acc):
            results[future_to_acc[future]] = future.result()
    return results

def main():
    parser = argparse.ArgumentParser(description="🌐 SRA智能下载工具 v4.7 - 不死鸟版")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--srr", help="单个SRA/CRR编号")
    group.add_argument("--list", help="批量编号列表文件 (TXT)")
    group.add_argument("--tsv", help="TSV矩阵文件 (第一行为样本名, 后续为编号)")

    parser.add_argument("--skip-list", help="已下载/需跳过的列表文件 (TXT)")
    
    parser.add_argument("--ngdc-method", choices=["requests", "aria2c", "wget"], default=DEFAULT_NGDC_METHOD)
    parser.add_argument("--ngdc-concurrency", type=int, default=DEFAULT_NGDC_CONCURRENCY)
    parser.add_argument("--prefetch-concurrency", type=int, default=DEFAULT_PREFETCH_CONCURRENCY)
    parser.add_argument("-o", "--output", default=DEFAULT_OUTPUT, help="基础输出目录")
    parser.add_argument("--no-fallback", action="store_true", help="禁用prefetch回退")
    args = parser.parse_args()

    skip_set = set()
    if args.skip_list and os.path.exists(args.skip_list):
        with open(args.skip_list, 'r', encoding='utf-8') as f:
            for line in f:
                acc = line.strip().upper()
                if acc and not acc.startswith('#'):
                    skip_set.add(acc)
        print(f"🛡️  已加载跳过清单: 发现 {len(skip_set)} 个预设的跳过样本。")
    elif args.skip_list:
        print(f"⚠️ 警告: 提供的跳过清单文件不存在 ({args.skip_list})，将正常下载所有项。")

    original_accessions = []
    acc_to_sample = {}

    if args.srr:
        acc = args.srr.strip().upper()
        original_accessions.append(acc)
        acc_to_sample[acc] = "Uncategorized"
    elif args.list:
        with open(args.list, 'r', encoding='utf-8') as f:
            for line in f:
                acc = line.strip().upper()
                if acc and not acc.startswith('#'):
                    original_accessions.append(acc)
                    acc_to_sample[acc] = "Uncategorized"
    elif args.tsv:
        print(f"📂 正在解析 TSV 实验设计表格: {args.tsv}")
        with open(args.tsv, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)
        if rows:
            headers = rows[0]
            num_cols = len(headers)
            for row in rows[1:]:
                row += [''] * (num_cols - len(row))
                for col_idx, cell in enumerate(row):
                    acc = cell.strip().upper()
                    if acc:
                        sample_name = sanitize_dirname(headers[col_idx])
                        original_accessions.append(acc)
                        acc_to_sample[acc] = sample_name
        else:
            sys.exit("❌ TSV文件为空！")

    original_accessions = list(dict.fromkeys(original_accessions))
    if not original_accessions: sys.exit("❌ 未找到有效的任务编号")
    print(f"\n📋 共识别到 {len(original_accessions)} 个有效任务项。")

    to_download = []
    skipped_accessions = []
    for acc in original_accessions:
        if acc in skip_set:
            skipped_accessions.append(acc)
        else:
            to_download.append(acc)

    if skipped_accessions:
        print(f"\n⏭️  触发过滤机制：自动跳过 {len(skipped_accessions)} 个已存在的样本。")
        for acc in skipped_accessions[:3]:
            print(f"   - 跳过: {acc} (归属: {acc_to_sample.get(acc, 'Unknown')})")
        if len(skipped_accessions) > 3:
            print(f"   ... 等共 {len(skipped_accessions)} 个")

    if not to_download:
        print("\n🎉 所有样本均在跳过清单中，无需下载，程序结束。")
        sys.exit(0)

    print(f"\n🚀 实际需要下载的样本数: {len(to_download)}")

    os.makedirs(args.output, exist_ok=True)
    
    print("\n" + "="*60)
    print(f"🔰 阶段一：NGDC并发下载 ({args.ngdc_method})")
    print("="*60)
    ngdc_results = download_batch_ngdc(to_download, acc_to_sample, args.ngdc_method, args.output, args.ngdc_concurrency)
    ngdc_success = [acc for acc, ok in ngdc_results.items() if ok]
    ngdc_failed = [acc for acc, ok in ngdc_results.items() if not ok]

    prefetch_success = []
    prefetch_failed = []
    
    if not args.no_fallback and ngdc_failed:
        print("\n" + "="*60)
        print(f"🔰 阶段二：自动触发 Prefetch 回退补漏")
        print("="*60)
        prefetch_results = download_batch_prefetch(ngdc_failed, acc_to_sample, args.output, args.prefetch_concurrency)
        prefetch_success = [acc for acc, ok in prefetch_results.items() if ok]
        prefetch_failed = [acc for acc, ok in prefetch_results.items() if not ok]
    elif args.no_fallback:
        prefetch_failed = ngdc_failed

    final_failed = prefetch_failed if not args.no_fallback else ngdc_failed

    date_str = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_csv = os.path.join(args.output, f"download_report_{date_str}.csv")
    failed_txt = os.path.join(args.output, f"failed_sra_{date_str}.txt")

    records = []
    for acc in original_accessions:
        sample = acc_to_sample.get(acc, "Unknown")
        if acc in skipped_accessions:
            status, method = "Skipped", "Pre-downloaded"
        elif acc in ngdc_success:
            status, method = "Success", f"NGDC ({args.ngdc_method})"
        elif acc in prefetch_success:
            status, method = "Success", "NCBI (prefetch)"
        else:
            status, method = "Failed", "None"
        records.append({"Accession": acc, "Sample": sample, "Method": method, "Status": status})

    with open(report_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["Accession", "Sample", "Method", "Status"])
        writer.writeheader()
        writer.writerows(records)

    if final_failed:
        with open(failed_txt, 'w', encoding='utf-8') as f:
            for acc in final_failed: f.write(f"{acc}\n")

    print("\n" + "="*60)
    print("🎯 下载流水线执行完毕")
    print(f"   ✓ 本次下载成功: {len(ngdc_success) + len(prefetch_success)} (NGDC:{len(ngdc_success)}, Prefetch:{len(prefetch_success)})")
    print(f"   ⏭️ 跳过已下载数: {len(skipped_accessions)}")
    print(f"   ✗ 最终失败任务: {len(final_failed)}")
    print(f"   📁 数据归档目录: {os.path.abspath(args.output)}")
    print(f"   📊 运行明细报表: {report_csv}")
    
    if final_failed:
        print(f"   ⚠️ 失败清单已导出: {failed_txt} (可直接作为下次的 --list 使用)")
    print("="*60)

if __name__ == "__main__":
    main()
