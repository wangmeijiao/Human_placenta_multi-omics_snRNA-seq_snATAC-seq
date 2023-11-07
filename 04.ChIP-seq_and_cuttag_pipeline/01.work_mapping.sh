
cat sample.txt|grep -v "^#" |while read a b; do echo $a $b;

#mkdir unmapped_${a}
bowtie2 -p 50  --no-unal   -x /database/all_index/index_bowtie2_hg38/Homo_sapiens.GRCh38 -1 ../00.data/${a}_R1.fq.gz -2 ../00.data/${a}_R2.fq.gz -S ${a}.cuttag.sam  > log.bowtie2.$a 2>&1
/home/mjwang/progs/misc-tools/sam2bai ${a}.cuttag.sam

done > log.run 2>&1 &
