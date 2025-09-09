#cd /datapool/pipeline/metagenome/upload/test_data/03.annotation/gtdbtk/fa
#less /datapool/pipeline/metagenome/upload/test_data/03.annotation/binning.list |perl -ne 'chomp;print"ln -s $_ ./\n"'|sh
cd /datapool/pipeline/metagenome/upload/test_data/03.annotation/gtdbtk/
export GTDBTK_DATA_PATH=/datapool/software/anaconda3/envs/qiime2/share/gtdbtk-1.7.0/db/release202/
/datapool/software/anaconda3/envs/qiime2/bin/gtdbtk  classify_wf --genome_dir /datapool/pipeline/metagenome/upload/test_data/03.annotation/gtdbtk/fa   --out_dir ./new_gtdb --extension fa --cpus 60
