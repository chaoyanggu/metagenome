/datapool/software/anaconda3/bin/metaspades.py -t 60 -m 500 -1 /datapool/test_data/01.cleandata/sample1/sample1.unmap.1.fq.gz -2 /datapool/test_data/01.cleandata/sample1/sample1.unmap.2.fq.gz  --only-assembler  -o /datapool/test_data/02.assembly/assembly/sample1
perl /datapool/pipeline/metagenome/bin/count_seq_len_GC.pl /datapool/test_data/02.assembly/assembly/sample1/contigs.fasta >/datapool/test_data/02.assembly/assembly/sample1/sample1.length.GC
/datapool/software/anaconda3/envs/seqkit/bin/seqkit  seq -m 1000 /datapool/test_data/02.assembly/assembly/sample1/contigs.fasta > /datapool/test_data/02.assembly/pick_1k_sequence/sample1.contig.fa
source activate metawrap
export PATH="/datapool/software/anaconda3/envs/metawrap/bin/:$PATH"
ln -s /datapool/test_data/01.cleandata/sample1/sample1.unmap.1.fq.gz /datapool/test_data/02.assembly/metawrap/sample1/fq/sample1_1.fastq
ln -s /datapool/test_data/01.cleandata/sample1/sample1.unmap.2.fq.gz /datapool/test_data/02.assembly/metawrap/sample1/fq/sample1_2.fastq
/datapool/software/anaconda3/envs/metawrap/bin/metawrap binning -o /datapool/test_data/02.assembly/metawrap/sample1/ -t 20 -a /datapool/test_data/02.assembly/pick_1k_sequence/sample1.contig.fa --metabat2 --maxbin2 --concoct /datapool/test_data/02.assembly/metawrap/sample1/fq/sample1_1.fastq /datapool/test_data/02.assembly/metawrap/sample1/fq/sample1_2.fastq
/datapool/software/anaconda3/envs/metawrap/bin/metawrap bin_refinement -o /datapool/test_data/02.assembly/metawrap/sample1/bin_refinement -t 12 -A /datapool/test_data/02.assembly/metawrap/sample1/metabat2_bins/ -B /datapool/test_data/02.assembly/metawrap/sample1/maxbin2_bins/ -C /datapool/test_data/02.assembly/metawrap/sample1/concoct_bins/ -c 50 -x 10
ls /datapool/test_data/02.assembly/metawrap/sample1/bin_refinement/metawrap_50_10_bins/*.fa > /datapool/test_data/02.assembly/metawrap/sample1/sample1.MG.list
less /datapool/test_data/02.assembly/metawrap/sample1/bin_refinement/metawrap_50_10_bins.stats  |sed '1d' |awk '$2>90&&$3<5' |awk -F "\t" '{print $1}'|perl -ne 'chomp;print qq(/datapool/test_data/02.assembly/metawrap/sample1/bin_refinement/metawrap_50_10_bins/$_.fa\n)' >/datapool/test_data/02.assembly/metawrap/sample1/sample1.HG.list
ls /datapool/test_data/02.assembly/metawrap/sample1/sample1.MG.list>>/datapool/test_data/03.annotation/binning.list
###rm

/datapool/software/anaconda3/bin/metaspades.py -t 60 -m 500 -1 /datapool/test_data/01.cleandata/sample2/sample2.unmap.1.fq.gz -2 /datapool/test_data/01.cleandata/sample2/sample2.unmap.2.fq.gz  --only-assembler  -o /datapool/test_data/02.assembly/assembly/sample2
perl /datapool/pipeline/metagenome/bin/count_seq_len_GC.pl /datapool/test_data/02.assembly/assembly/sample2/contigs.fasta >/datapool/test_data/02.assembly/assembly/sample2/sample2.length.GC
/datapool/software/anaconda3/envs/seqkit/bin/seqkit  seq -m 1000 /datapool/test_data/02.assembly/assembly/sample2/contigs.fasta > /datapool/test_data/02.assembly/pick_1k_sequence/sample2.contig.fa
source activate metawrap
export PATH="/datapool/software/anaconda3/envs/metawrap/bin/:$PATH"
ln -s /datapool/test_data/01.cleandata/sample2/sample2.unmap.1.fq.gz /datapool/test_data/02.assembly/metawrap/sample2/fq/sample2_1.fastq
ln -s /datapool/test_data/01.cleandata/sample2/sample2.unmap.2.fq.gz /datapool/test_data/02.assembly/metawrap/sample2/fq/sample2_2.fastq
/datapool/software/anaconda3/envs/metawrap/bin/metawrap binning -o /datapool/test_data/02.assembly/metawrap/sample2/ -t 20 -a /datapool/test_data/02.assembly/pick_1k_sequence/sample2.contig.fa --metabat2 --maxbin2 --concoct /datapool/test_data/02.assembly/metawrap/sample2/fq/sample2_1.fastq /datapool/test_data/02.assembly/metawrap/sample2/fq/sample2_2.fastq
/datapool/software/anaconda3/envs/metawrap/bin/metawrap bin_refinement -o /datapool/test_data/02.assembly/metawrap/sample2/bin_refinement -t 12 -A /datapool/test_data/02.assembly/metawrap/sample2/metabat2_bins/ -B /datapool/test_data/02.assembly/metawrap/sample2/maxbin2_bins/ -C /datapool/test_data/02.assembly/metawrap/sample2/concoct_bins/ -c 50 -x 10
ls /datapool/test_data/02.assembly/metawrap/sample2/bin_refinement/metawrap_50_10_bins/*.fa > /datapool/test_data/02.assembly/metawrap/sample2/sample2.MG.list
less /datapool/test_data/02.assembly/metawrap/sample2/bin_refinement/metawrap_50_10_bins.stats  |sed '1d' |awk '$2>90&&$3<5' |awk -F "\t" '{print $1}'|perl -ne 'chomp;print qq(/datapool/test_data/02.assembly/metawrap/sample2/bin_refinement/metawrap_50_10_bins/$_.fa\n)' >/datapool/test_data/02.assembly/metawrap/sample2/sample2.HG.list
ls /datapool/test_data/02.assembly/metawrap/sample2/sample2.MG.list>>/datapool/test_data/03.annotation/binning.list
###rm

