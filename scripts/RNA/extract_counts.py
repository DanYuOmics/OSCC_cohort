import os
import glob
import re
import pandas as pd

quant_dir = "/data/OSCC2/fastq/RNA/ballgown_out"
gtf_files = sorted(glob.glob(f"{quant_dir}/*/*.gtf"))

if not gtf_files:
    # 兼容直接存放在 assemble 目录的情况
    gtf_files = sorted(glob.glob("/data/OSCC2/fastq/RNA/stringtie_assemble/*.gtf"))

print(f"Found {len(gtf_files)} GTF files.")

gene_counts = {}
transcript_counts = {}
samples = []

# 正则匹配属性字段
re_gene_id = re.compile(r'gene_id "([^"]+)"')
re_tx_id = re.compile(r'transcript_id "([^"]+)"')
re_cov = re.compile(r'cov "([^"]+)"')

# 平均 Read Length (PE150 设为 150)
READ_LEN = 150.0

for gtf in gtf_files:
    sample = os.path.basename(gtf).replace(".gtf", "")
    samples.append(sample)
    print(f"Parsing sample: {sample}")
    
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            
            info = fields[8]
            g_match = re_gene_id.search(info)
            t_match = re_tx_id.search(info)
            cov_match = re_cov.search(info)
            
            if not g_match or not t_match:
                continue
            
            gene_id = g_match.group(1)
            tx_id = t_match.group(1)
            
            # 转录本长度
            tx_len = int(fields[4]) - int(fields[3]) + 1
            cov = float(cov_match.group(1)) if cov_match else 0.0
            
            # 官方 prepDE read count 计算公式: (coverage * transcript_length) / read_length
            cnt = int(round((cov * tx_len) / READ_LEN))
            
            # 存入 Transcript 矩阵
            if tx_id not in transcript_counts:
                transcript_counts[tx_id] = {}
            transcript_counts[tx_id][sample] = cnt
            
            # 汇总到 Gene 矩阵
            if gene_id not in gene_counts:
                gene_counts[gene_id] = {}
            gene_counts[gene_id][sample] = gene_counts[gene_id].get(sample, 0) + cnt

# 构建 DataFrame 并将缺失值补 0
df_gene = pd.DataFrame.from_dict(gene_counts, orient='index', columns=samples).fillna(0).astype(int)
df_tx = pd.DataFrame.from_dict(transcript_counts, orient='index', columns=samples).fillna(0).astype(int)

df_gene.index.name = "gene_id"
df_tx.index.name = "transcript_id"

out_gene = "/data/OSCC2/fastq/RNA/gene_count_matrix.csv"
out_tx = "/data/OSCC2/fastq/RNA/transcript_count_matrix.csv"

df_gene.to_csv(out_gene)
df_tx.to_csv(out_tx)

print("\nMatrix extraction completed successfully:")
print(f"Gene matrix: {out_gene} | Dimensions: {df_gene.shape}")
print(f"Transcript matrix: {out_tx} | Dimensions: {df_tx.shape}")
