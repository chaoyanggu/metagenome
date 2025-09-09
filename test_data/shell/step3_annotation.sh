/datapool/software/anaconda3/bin/abricate /datapool/pipeline/metagenome/upload/test_data/02.assembly/pick_1k_sequence/sample1.contig.fa --db resfinder_new --minid=75 >/datapool/pipeline/metagenome/upload/test_data/03.annotation/resfinder/result/sample1_resfinder_new.tab
/datapool/software/anaconda3/bin/abricate /datapool/pipeline/metagenome/upload/test_data/02.assembly/pick_1k_sequence/sample1.contig.fa --db vfdb_new --minid=75 >/datapool/pipeline/metagenome/upload/test_data/03.annotation/vfdb/result/sample1_vfdb_new.tab
export PATH="/datapool/software/anaconda3/envs/prokka/bin/:$PATH"
/datapool/software/anaconda3/envs/prokka/bin/prokka --force --cpus 10 --outdir /datapool/pipeline/metagenome/upload/test_data/03.annotation/prokka/sample1 --prefix sample1 --locustag sample1 --metagenome --kingdom Bacteria /datapool/pipeline/metagenome/upload/test_data/02.assembly/pick_1k_sequence/sample1.contig.fa
###rm

/datapool/software/anaconda3/bin/abricate /datapool/pipeline/metagenome/upload/test_data/02.assembly/pick_1k_sequence/sample2.contig.fa --db resfinder_new --minid=75 >/datapool/pipeline/metagenome/upload/test_data/03.annotation/resfinder/result/sample2_resfinder_new.tab
/datapool/software/anaconda3/bin/abricate /datapool/pipeline/metagenome/upload/test_data/02.assembly/pick_1k_sequence/sample2.contig.fa --db vfdb_new --minid=75 >/datapool/pipeline/metagenome/upload/test_data/03.annotation/vfdb/result/sample2_vfdb_new.tab
export PATH="/datapool/software/anaconda3/envs/prokka/bin/:$PATH"
/datapool/software/anaconda3/envs/prokka/bin/prokka --force --cpus 10 --outdir /datapool/pipeline/metagenome/upload/test_data/03.annotation/prokka/sample2 --prefix sample2 --locustag sample2 --metagenome --kingdom Bacteria /datapool/pipeline/metagenome/upload/test_data/02.assembly/pick_1k_sequence/sample2.contig.fa
###rm

