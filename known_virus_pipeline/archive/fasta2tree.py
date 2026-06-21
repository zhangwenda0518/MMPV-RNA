#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
import os
import logging

def setup_logger(log_file=None):
    """配置双轨日志输出系统 (屏幕看进度，文件存细节)"""
    logger = logging.getLogger("fasta2tree")
    logger.setLevel(logging.DEBUG)

    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger

logger = setup_logger()

def run_command(cmd, step_name):
    """拦截外部命令的输出并分发"""
    logger.info(f"========== 开始执行: {step_name} ==========")
    logger.info(f"运行命令: {cmd}")
    try:
        result = subprocess.run(
            cmd, shell=True, check=True, 
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        if result.stdout:
            logger.debug(f"\n[--- {step_name} 详细输出 ---]\n{result.stdout.strip()}\n")
            
        logger.info(f"========== {step_name} 执行成功 ==========\n")
    except subprocess.CalledProcessError as e:
        logger.error(f"{step_name} 执行失败! 退出码: {e.returncode}")
        if e.stdout:
            logger.error(f"\n[--- {step_name} 错误详情 ---]\n{e.stdout.strip()}\n")
        sys.exit(1)

def detect_seqtype(fasta_file):
    """自动检测序列类型 (nt / aa)"""
    logger.info("========== 正在自动检测序列类型 (nt/aa) ==========")
    valid_nt_chars = set("ACGTUNacgtun")
    total_chars = 0
    nt_chars = 0
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq = line.strip().replace("-", "").replace(" ", "").replace("?", "")
            if not seq:
                continue
                
            total_chars += len(seq)
            nt_chars += sum(1 for c in seq if c in valid_nt_chars)
            
            if total_chars > 10000:
                break
                
    if total_chars == 0:
        logger.error("FASTA 文件似乎没有包含任何有效的序列内容！")
        sys.exit(1)
        
    ratio = nt_chars / total_chars
    if ratio > 0.85:
        logger.info(f"🎯 检测结果: 核酸 (nt)。 [A/C/G/T/U/N 占比: {ratio:.1%}]\n")
        return 'nt'
    else:
        logger.info(f"🎯 检测结果: 蛋白质 (aa)。 [A/C/G/T/U/N 占比: {ratio:.1%}]\n")
        return 'aa'

def clean_fasta_headers(input_fasta, output_fasta):
    """清洗序列名并查重"""
    logger.info("========== 开始清洗 FASTA 序列名并查重 ==========")
    seen_ids = set()
    duplicate_found = False
    
    with open(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                clean_id = line.strip().split()[0]
                if clean_id in seen_ids:
                    logger.warning(f"发现重复的序列 ID: {clean_id}")
                    duplicate_found = True
                seen_ids.add(clean_id)
                fout.write(clean_id + '\n')
            else:
                fout.write(line)
                
    if not duplicate_found:
        logger.info("未发现重复的序列 ID。")
    logger.info(f"清洗后的序列已保存至: {output_fasta}\n")

def parse_best_model(modeltest_out_file, criterion="BIC"):
    """解析 modeltest-ng 的 .out 文件，提取最佳模型"""
    logger.info(f"========== 正在解析 ModelTest-NG 结果 ({criterion}) ==========")
    if not os.path.isfile(modeltest_out_file):
        logger.error(f"找不到 ModelTest-NG 的输出文件: {modeltest_out_file}")
        return None

    best_model = None
    in_summary = False
    
    with open(modeltest_out_file, 'r') as f:
        for line in f:
            clean_line = line.strip()
            if clean_line.startswith("Summary:"):
                in_summary = True
                continue
            
            if in_summary and clean_line.startswith(criterion):
                parts = clean_line.split()
                if len(parts) >= 2:
                    best_model = parts[1]
                    break

    if best_model:
        logger.info(f"成功提取最佳模型: {best_model}\n")
        return best_model
    else:
        logger.warning(f"解析失败！未能从 {modeltest_out_file} 中找到 {criterion} 模型。\n")
        return None

def main():
    parser = argparse.ArgumentParser(
        description="fasta2tree: 终极容错版 (MAFFT -> ClipKIT -> ModelTest -> RAxML-NG 或 IQ-TREE 智能切换)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input", required=True, help="输入的原始 FASTA 序列文件")
    parser.add_argument("-o", "--prefix", required=True, help="输出文件的前缀 (如: mydata)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="使用的 CPU 线程数")
    parser.add_argument("--seqtype", choices=['nt', 'aa', 'auto'], default='auto', 
                        help="序列类型: 'nt', 'aa', 或 'auto' (自动推断)")
    parser.add_argument("--trimmer", choices=['clipkit', 'trimal', 'gblocks', 'none'], default='clipkit', 
                        help="使用的比对修剪工具")
    parser.add_argument("--run_modeltest", action="store_true", 
                        help="开启此选项将运行 modeltest-ng 并自动解析最佳模型")
    parser.add_argument("--criterion", choices=['BIC', 'AIC', 'AICc'], default='BIC',
                        help="配合 --run_modeltest 使用，选择模型评估准则")
    parser.add_argument("--model", type=str, default="DEFAULT", 
                        help="手动指定建树模型。若为 DEFAULT，则 nt 默认 GTR+G4，aa 默认 LG+G4")

    args = parser.parse_args()

    # 重载全局 logger 并绑定到文件
    global logger
    log_filename = f"{args.prefix}.run_pipeline.log"
    logger = setup_logger(log_file=log_filename)
    logger.info(f"🗂️ 运行详情已被静默拦截，完整日志将实时写入: {log_filename}")

    if not os.path.isfile(args.input):
        logger.error(f"输入文件 {args.input} 不存在！")
        sys.exit(1)

    # 1. 检测序列类型
    seq_type = args.seqtype
    if seq_type == 'auto':
        seq_type = detect_seqtype(args.input)

    cleaned_fasta = f"{args.prefix}.cleaned.fasta"
    mafft_out = f"{args.prefix}.mafft.fasta"
    trimmed_out = f"{args.prefix}.mafft.trimmed.fasta"
    
    # 2. 清洗
    clean_fasta_headers(args.input, cleaned_fasta)

    # 3. 比对
    mafft_cmd = f"mafft --thread {args.threads} --adjustdirection --auto {cleaned_fasta} > {mafft_out}"
    run_command(mafft_cmd, "MAFFT Alignment")

    # 4. 修剪
    if args.trimmer == 'clipkit':
        trim_cmd = f"clipkit {mafft_out} -o {trimmed_out}"
        run_command(trim_cmd, "ClipKIT Trimming")
    elif args.trimmer == 'trimal':
        trim_cmd = f"trimal -in {mafft_out} -out {trimmed_out} -automated1"
        run_command(trim_cmd, "trimAl Trimming")
    elif args.trimmer == 'gblocks':
        t_type = 'd' if seq_type == 'nt' else 'p'
        trim_cmd = f"Gblocks {mafft_out} -t={t_type} -e=.gb && mv {mafft_out}.gb {trimmed_out}"
        run_command(trim_cmd, "Gblocks Trimming")
    elif args.trimmer == 'none':
        logger.info("跳过比对修剪步骤，直接使用 MAFFT 输出结果。")
        trimmed_out = mafft_out

    # 5. 模型测试与建树工具选择策略
    fallback_model = "GTR+G4" if seq_type == 'nt' else "LG+G4"
    final_model = args.model if args.model != "DEFAULT" else fallback_model
    use_iqtree_fallback = False # 控制是否切换到 IQ-TREE 的开关
    
    if args.run_modeltest:
        # 尝试运行 ModelTest-NG
        # 注意: 即使 modeltest 报错崩溃 (退出码非0), run_command 会终止程序。
        # 这里假设它正常跑完但没出有用的结果(比如全是警告)。
        modeltest_cmd = f"modeltest-ng -i {trimmed_out} -d {seq_type} -p {args.threads}"
        try:
            run_command(modeltest_cmd, "ModelTest-NG")
            modeltest_out_file = f"{trimmed_out}.out"
            parsed_model = parse_best_model(modeltest_out_file, criterion=args.criterion)
            
            if parsed_model:
                final_model = parsed_model
                logger.info(f"==> 将使用 ModelTest-NG 推荐模型: [{final_model}] (RAxML-NG 建树)")
            else:
                logger.warning("==> ⚠️ ModelTest-NG 解析不到合适模型！触发容错机制：自动移交 IQ-TREE 🚀")
                use_iqtree_fallback = True
        except SystemExit:
            # 如果 modeltest-ng 自身彻底崩溃导致 run_command 退出
            logger.warning("==> ⚠️ ModelTest-NG 运行崩溃！触发容错机制：自动移交 IQ-TREE 🚀")
            use_iqtree_fallback = True
    else:
        logger.info(f"==> 跳过模型测试，直接使用模型: [{final_model}] (RAxML-NG 建树)")

    # 6. 建树
    if use_iqtree_fallback:
        # 使用 IQ-TREE 的自适应模型测试 (MFP) 及你要求的参数
        iqtree_cmd = f"iqtree -s {trimmed_out} -m MFP -B 1000 -alrt 1000 -T {args.threads}"
        # 针对极短/极度保守序列可能出错的情况，让 IQ-TREE 尽力而为
        run_command(iqtree_cmd, "IQ-TREE 自动模型推断与建树")
        
        logger.info("🎉 所有的分析流程已顺利完成！")
        logger.info(f"🌲 最终的树文件 (Best Tree) : {trimmed_out}.treefile")
        logger.info(f"📊 完整的分析报告 (Report)  : {trimmed_out}.iqtree")
    else:
        # 传统 RAxML-NG 路线
        raxml_cmd = f"raxml-ng --msa {trimmed_out} --all --force --threads {args.threads} --model {final_model}"
        run_command(raxml_cmd, "RAxML-NG Tree Building")
        
        logger.info("🎉 所有的分析流程已顺利完成！")
        logger.info(f"🌲 最终的树文件 (Best Tree) : {trimmed_out}.raxml.bestTree")
        logger.info(f"📊 支持度树文件 (Support)  : {trimmed_out}.raxml.support")

    logger.info(f"💡 需要查看软件详细参数和报错细节，请检查日志: {log_filename}")

if __name__ == "__main__":
    main()
