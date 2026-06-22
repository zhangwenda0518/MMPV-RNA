from transformers import AutoTokenizer, AutoModelForSequenceClassification
from Bio import SeqIO
import numpy as np
import argparse
import torch
import csv
import os
import shutil
import time
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import multiprocessing as mp

parser = argparse.ArgumentParser(description='ViraLM v1.0 - Multi-process Parallel')
parser.add_argument('--input', '-i', type=str, required=True, help='input FASTA file')
parser.add_argument('--output', '-o', type=str, default='result', help='output directory')
parser.add_argument('--database', '-d', type=str, required=True, help='model directory')
parser.add_argument('--processes', type=int, default=4, help='number of processes')
parser.add_argument('--batch_size', type=int, default=64, help='batch size for prediction')
parser.add_argument('--len', type=int, default=500, help='minimum sequence length')
parser.add_argument('--threshold', type=float, default=0.5, help='prediction threshold')
parser.add_argument('--force', '-f', action='store_true', help='force overwrite output directory')
parser.add_argument('--filename', '-n', type=str, default=None, help='custom output name')
parser.add_argument('--max_len', type=int, default=2000, help='maximum sequence chunk size')
parser.add_argument('--chunk_size', type=int, default=1000, help='sequences per process')
args = parser.parse_args()

# 禁用tokenizer并行以避免冲突
os.environ['TOKENIZERS_PARALLELISM'] = 'false'

input_pth = args.input
output_pth = args.output
model_pth = args.database
num_processes = min(args.processes, mp.cpu_count())
batch_size = args.batch_size
len_threshold = args.len
score_threshold = args.threshold
max_len = args.max_len
chunk_size = args.chunk_size
cache_dir = f'{output_pth}/cache'
filename = args.filename if args.filename else os.path.splitext(os.path.basename(input_pth))[0]

# 验证参数
if score_threshold < 0.5:
    print('Error: Threshold must be >= 0.5')
    exit(1)

if not os.path.exists(model_pth):
    print(f'Error: Model directory {model_pth} missing')
    exit(1)

# 创建输出目录
if os.path.isdir(output_pth):
    if args.force:
        shutil.rmtree(output_pth)
        os.makedirs(output_pth)
    else:
        print('Error: Output directory exists. Use -f to overwrite.')
        exit(1)
else:
    os.makedirs(output_pth)

os.makedirs(cache_dir, exist_ok=True)
print(f"Starting ViraLM prediction (Multi-process Mode)")
print(f"Input: {input_pth}")
print(f"Output: {output_pth}")
print(f"Number of processes: {num_processes}")
print(f"Chunk size: {chunk_size} sequences per process")

def is_valid_dna(seq):
    """验证DNA序列"""
    return set(seq.upper()).issubset(set('ACGT'))

def prepare_data_chunks(input_path, min_len=500, max_chunk=2000, chunk_size=1000):
    """准备数据块用于多进程处理，只处理长度>=min_len的序列"""
    print(f"Preparing data chunks (minimum length: {min_len}bp)...")
    all_chunks = []
    current_chunk = []
    current_accessions = []
    
    total_sequences = 0
    skipped_short = 0
    skipped_invalid = 0
    
    for record in tqdm(list(SeqIO.parse(input_path, "fasta")), desc="Reading sequences"):
        seq = str(record.seq).upper()
        seq_len = len(seq)
        total_sequences += 1
        
        # 跳过长度小于阈值的序列
        if seq_len < min_len:
            skipped_short += 1
            continue
        
        # 验证序列是否为有效DNA
        if not is_valid_dna(seq):
            skipped_invalid += 1
            continue
        
        # 长序列分割
        if seq_len > max_chunk:
            for i in range(0, seq_len, max_chunk):
                chunk = seq[i:i+max_chunk]
                # 确保分段长度至少为500bp
                if len(chunk) >= 500 and is_valid_dna(chunk):
                    current_chunk.append(chunk)
                    current_accessions.append(f"{record.id}_{i}_{i+len(chunk)}")
                    
                    if len(current_chunk) >= chunk_size:
                        all_chunks.append((current_chunk.copy(), current_accessions.copy()))
                        current_chunk = []
                        current_accessions = []
        else:
            current_chunk.append(seq)
            current_accessions.append(f"{record.id}_0_{seq_len}")
            
            if len(current_chunk) >= chunk_size:
                all_chunks.append((current_chunk.copy(), current_accessions.copy()))
                current_chunk = []
                current_accessions = []
    
    # 添加最后一个块
    if current_chunk:
        all_chunks.append((current_chunk, current_accessions))
    
    print(f"Preprocessing summary:")
    print(f"  Total sequences read: {total_sequences}")
    print(f"  Skipped (length < {min_len}bp): {skipped_short}")
    print(f"  Skipped (invalid DNA): {skipped_invalid}")
    print(f"  Processed sequences: {total_sequences - skipped_short - skipped_invalid}")
    print(f"  Created chunks for processing: {len(all_chunks)}")
    
    return all_chunks

# 准备数据
all_chunks = prepare_data_chunks(input_pth, len_threshold, max_len, chunk_size)

# 检查是否有数据
if len(all_chunks) == 0:
    print("Error: No sequences meet the length threshold or are valid DNA sequences.")
    exit(1)

def process_chunk(chunk_data, model_path, batch_size=64, score_threshold=0.5):
    """处理一个数据块"""
    sequences, accessions = chunk_data
    
    # 在每个进程中加载模型
    model = AutoModelForSequenceClassification.from_pretrained(
        model_path,
        num_labels=2,
        trust_remote_code=True,
        local_files_only=True,  # 避免重复下载
    )
    
    tokenizer = AutoTokenizer.from_pretrained(
        model_path,
        model_max_length=512,
        padding_side="right",
        use_fast=True,
        trust_remote_code=True,
    )
    
    model.eval()
    results = {}
    
    # 按批次处理
    for i in range(0, len(sequences), batch_size):
        batch_seqs = sequences[i:i+batch_size]
        batch_accs = accessions[i:i+batch_size]
        
        # 编码
        inputs = tokenizer(
            batch_seqs,
            truncation=True,
            padding=True,
            max_length=512,
            return_tensors='pt'
        )
        
        # 预测
        with torch.no_grad():
            outputs = model(**inputs)
            probs = torch.nn.functional.softmax(outputs.logits, dim=-1)
            virus_probs = probs[:, 1].cpu().numpy()
        
        # 收集结果
        for acc, prob in zip(batch_accs, virus_probs):
            seq_name = acc.rsplit('_', 2)[0]
            if seq_name not in results:
                results[seq_name] = []
            results[seq_name].append(prob)
    
    return results

# 多进程处理
print(f"Starting multi-process prediction with {num_processes} processes...")
start_time = time.time()

# 创建进程池
with ProcessPoolExecutor(max_workers=num_processes) as executor:
    # 提交任务
    future_to_chunk = {
        executor.submit(
            process_chunk, 
            chunk, 
            model_pth, 
            batch_size, 
            score_threshold
        ): i 
        for i, chunk in enumerate(all_chunks)
    }
    
    # 收集结果
    all_results = []
    for future in tqdm(as_completed(future_to_chunk), total=len(future_to_chunk), desc="Processing chunks"):
        try:
            result = future.result()
            all_results.append(result)
        except Exception as e:
            print(f"Error processing chunk: {e}")

# 合并结果
print("Merging results...")
final_results = {}
virus_set = set()

for chunk_result in all_results:
    for seq_name, probs in chunk_result.items():
        if seq_name not in final_results:
            final_results[seq_name] = []
        final_results[seq_name].extend(probs)

# 计算平均概率
processed_results = {}
for seq_name, probs in final_results.items():
    avg_prob = np.mean(probs)
    processed_results[seq_name] = avg_prob
    if avg_prob > score_threshold:
        virus_set.add(seq_name)

# 排序
processed_results = dict(sorted(processed_results.items(), key=lambda x: x[1], reverse=True))

# 保存结果
result_file = f'{output_pth}/result_{filename}.csv'
with open(result_file, 'w') as f:
    f.write('seq_name,prediction,virus_score\n')
    for seq_name, score in processed_results.items():
        pred = 'virus' if score > score_threshold else 'non-virus'
        f.write(f'{seq_name},{pred},{score:.6f}\n')

# 保存病毒序列
if virus_set:
    virus_records = []
    for record in SeqIO.parse(input_pth, "fasta"):
        if record.id in virus_set:
            virus_records.append(record)
    
    if virus_records:
        SeqIO.write(virus_records, f'{output_pth}/virus_{filename}.fasta', 'fasta')
        print(f"Found {len(virus_records)} virus sequences")

print(f"\nTotal processing time: {time.time() - start_time:.2f}s")
print(f"Results saved to: {result_file}")
