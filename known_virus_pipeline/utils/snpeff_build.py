#!/usr/bin/env python3
import os
import sys
import shutil
import argparse
import subprocess
import urllib.request
import logging

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

class SnpEffBuilder:
    def __init__(self, jar, config, memory="4g"):
        self.jar = os.path.abspath(os.path.expanduser(jar))
        self.config = os.path.abspath(os.path.expanduser(config))
        self.memory = memory
        self.data_dir = os.path.join(os.path.dirname(self.config), "data")

    def _register_db(self, db_name):
        with open(self.config, "r") as f:
            if f"{db_name}.genome" in f.read():
                return
        with open(self.config, "a") as f:
            f.write(f"\n{db_name}.genome : {db_name}\n")

    def build_from_ncbi(self, db_name, accessions):
        db_path = os.path.join(self.data_dir, db_name)
        os.makedirs(db_path, exist_ok=True)
        gbk_file = os.path.join(db_path, "genes.gbk")
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(accessions)}&rettype=gbwithparts&retmode=text"
        
        logger.info(f"正在下载 GenBank 数据...")
        urllib.request.urlretrieve(url, gbk_file)
        self._register_db(db_name)
        
        cmd = ["java", f"-Xmx{self.memory}", "-jar", self.jar, "build", "-genbank", "-v", "-c", self.config, db_name]
        subprocess.run(cmd, check=True)

    def build_from_local(self, db_name, fa, gff):
        if not shutil.which("gffread"):
            logger.error("请先安装 gffread"); sys.exit(1)
        
        db_path = os.path.join(self.data_dir, db_name)
        os.makedirs(db_path, exist_ok=True)
        
        # 准备 SnpEff 标准命名文件
        try:
            subprocess.run(["gffread", gff, "-T", "-o", os.path.join(db_path, "genes.gtf")], check=True)
            shutil.copy(fa, os.path.join(db_path, "sequences.fa"))
            self._register_db(db_name)
            
            cmd = ["java", f"-Xmx{self.memory}", "-jar", self.jar, "build", "-gtf22", "-v", "-c", self.config, db_name]
            subprocess.run(cmd, check=True)
        except Exception as e:
            logger.error(f"构建失败: {e}")

def main():
    parser = argparse.ArgumentParser(description="SnpEff 数据库构建脚本")
    parser.add_argument("--db", required=True, help="数据库名称")
    parser.add_argument("--jar", default="~/snpEff/snpEff.jar")
    parser.add_argument("--config", default="~/snpEff/snpEff.config")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--acc", nargs="+", help="NCBI Accessions")
    group.add_argument("--local", nargs=2, metavar=("FASTA", "GFF"))
    
    args = parser.parse_args()
    builder = SnpEffBuilder(args.jar, args.config)
    
    if args.acc:
        builder.build_from_ncbi(args.db, args.acc)
    else:
        builder.build_from_local(args.db, args.local[0], args.local[1])

if __name__ == "__main__":
    main()
