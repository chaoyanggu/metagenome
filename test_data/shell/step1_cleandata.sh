cd /datapool/pipeline/metagenome/upload/test_data/01.cleandata/sample1
/datapool/software/anaconda3/envs/RNAseq/bin/fastp  -i /datapool/bioinfo/ngs/sample1.R1.fastq.gz -o sample1_R1.clean.fq.gz  -I /datapool/bioinfo/ngs/sample1.R2.fastq.gz -O sample1_R2.clean.fq.gz -j sample1.json -h sample1.html
/datapool/software/anaconda3/bin/bowtie2  --very-sensitive  -p 30 -x /datapool/db/hg38/hg38  -1 sample1_R1.clean.fq.gz  -2 sample1_R2.clean.fq.gz --al-conc-gz sample1.map.fq.gz  --un-conc-gz  sample1.unmap.fq.gz  -S sample1.sam 2> bowtie2.log
/datapool/software/anaconda3/bin/samtools view -F 4 -Sb sample1.sam  > sample1.bam
/datapool/software/anaconda3/bin/samtools sort sample1.bam  -o sample1.sorted.bam
mv sample1.unmap.fq.1.gz sample1.unmap.1.fq.gz
mv sample1.unmap.fq.2.gz sample1.unmap.2.fq.gz
/datapool/software/anaconda3/envs/qiime2/bin/fastqc --extract sample1.unmap.1.fq.gz -o /datapool/pipeline/metagenome/upload/test_data/01.cleandata/nohostqc_result
/datapool/software/anaconda3/envs/qiime2/bin/fastqc --extract sample1.unmap.2.fq.gz -o /datapool/pipeline/metagenome/upload/test_data/01.cleandata/nohostqc_result
/datapool/software/anaconda3/bin/kraken2 -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  sample1.unmap.1.fq.gz  sample1.unmap.2.fq.gz  --threads 20 --output sample1.kraken2 --use-names --report sample1.label --use-mpa-style
###rm sample1.sam

cd /datapool/pipeline/metagenome/upload/test_data/01.cleandata/sample2
/datapool/software/anaconda3/envs/RNAseq/bin/fastp  -i /datapool/bioinfo/ngs/sample2.R1.fastq.gz -o sample2_R1.clean.fq.gz  -I /datapool/bioinfo/ngs/sample2.R2.fastq.gz -O sample2_R2.clean.fq.gz -j sample2.json -h sample2.html
/datapool/software/anaconda3/bin/bowtie2  --very-sensitive  -p 30 -x /datapool/db/hg38/hg38  -1 sample2_R1.clean.fq.gz  -2 sample2_R2.clean.fq.gz --al-conc-gz sample2.map.fq.gz  --un-conc-gz  sample2.unmap.fq.gz  -S sample2.sam 2> bowtie2.log
/datapool/software/anaconda3/bin/samtools view -F 4 -Sb sample2.sam  > sample2.bam
/datapool/software/anaconda3/bin/samtools sort sample2.bam  -o sample2.sorted.bam
mv sample2.unmap.fq.1.gz sample2.unmap.1.fq.gz
mv sample2.unmap.fq.2.gz sample2.unmap.2.fq.gz
/datapool/software/anaconda3/envs/qiime2/bin/fastqc --extract sample2.unmap.1.fq.gz -o /datapool/pipeline/metagenome/upload/test_data/01.cleandata/nohostqc_result
/datapool/software/anaconda3/envs/qiime2/bin/fastqc --extract sample2.unmap.2.fq.gz -o /datapool/pipeline/metagenome/upload/test_data/01.cleandata/nohostqc_result
/datapool/software/anaconda3/bin/kraken2 -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  sample2.unmap.1.fq.gz  sample2.unmap.2.fq.gz  --threads 20 --output sample2.kraken2 --use-names --report sample2.label --use-mpa-style
###rm sample2.sam

