cd /datapool/pipeline/metagenome/upload/test_data/00.rawdata/
/datapool/software/anaconda3/envs/MetaPhage/bin/multiqc qc_result/ -o multiqc_output/ --filename multiqc_report.html
perl /datapool/pipeline/metagenome/bin/merge_fastqc.pl multiqc_output/multiqc_report_data/multiqc_general_stats.txt >multiqc.txt
perl /datapool/pipeline/metagenome/bin/sample_stats_extractor.pl /datapool/pipeline/metagenome/upload/test_data/01.cleandata/*/*.json >/datapool/pipeline/metagenome/upload/test_data/01.cleandata/sample.cleandata.stats
ls /datapool/pipeline/metagenome/upload/test_data/01.cleandata/*/bowtie2.log|while read f;do grep 'overall alignment rate' $f|awk -v f=$f '{print f"\t"$1}' ;done |perl -ne 'chomp;@a=split/\t/,$_;$id=(split/\//,$a[0])[-2];print"$id\t$a[-1]\n"' >/datapool/pipeline/metagenome/upload/test_data/01.cleandata/bowtie2.txt

cd /datapool/pipeline/metagenome/upload/test_data/01.cleandata/kraken_count
perl /datapool/pipeline/metagenome/bin/merge.kraken2.pl /datapool/pipeline/metagenome/upload/test_data/01.cleandata/*/*.label
ls *relative.txt|sed 's/.txt//'|perl -ne 'chomp;print"perl /datapool/pipeline/metagenome/bin/sortmat1-new.pl $_.txt > $_.sorted.txt\nperl /datapool/pipeline/metagenome/bin/toptable.pl $_.sorted.txt 10 > $_.top10.txt\nperl /datapool/pipeline/metagenome/bin/trantab.pl $_.top10.txt > $_.top10.trans.txt\nperl /datapool/pipeline/metagenome/bin/plot.top10-2.pl $_.top10.trans.txt $_.top10.svg\n/datapool/software/anaconda3/envs/DNBC4tools/bin/convert -density 300 $_.top10.svg $_.top10.png\n"' |sh

cd /datapool/pipeline/metagenome/upload/test_data/02.assembly/
find pick_1k_sequence -name "*.fa" | xargs /datapool/software/anaconda3/envs/seqkit/bin/seqkit stats -a | awk 'NR==1 || FNR>1 {print}' >seqkit.txt

