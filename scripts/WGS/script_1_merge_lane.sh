cd /data/OSCC2/fastq/WGS/tumour
mkdir -p raw_split_lanes
for s in S19TB S21TB S24TB S25TB; do
    cat ${s}_*_1.fq.gz > ${s}_combined_1.fq.gz
    cat ${s}_*_2.fq.gz > ${s}_combined_2.fq.gz
    mv ${s}_*L0*.fq.gz raw_split_lanes/
done
