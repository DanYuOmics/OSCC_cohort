cd /data/OSCC2/fastq/RNA

# 1. 确认 48 个 GTF 均已就绪并生成清单
ls -1 /data/OSCC2/fastq/RNA/stringtie_assemble/*.gtf > mergelist.txt

# 2. 执行合并
stringtie --merge \
    -p 16 \
    -G /data/yudan_ref/RNA/genome.gtf \
    -o /data/OSCC2/fastq/RNA/stringtie_merged.gtf \
    mergelist.txt
