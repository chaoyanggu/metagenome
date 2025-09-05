### step0_rawdata
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 0) step0_rawdata" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
perl /datapool/pipeline/metagenome/bin/dsub_gcy2.pl -thd 10 -mem 4 /datapool/pipeline/metagenome/upload/test_data/shell/step0_rawdata.sh 
date +"%D %T -> Finish 0) step0_rawdata" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step1_cleandata
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 1) step1_cleandata" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
perl /datapool/pipeline/metagenome/bin/dsub_gcy2.pl -thd 10 -mem 4 /datapool/pipeline/metagenome/upload/test_data/shell/step1_cleandata.sh 
date +"%D %T -> Finish 1) step1_cleandata" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step2_assembly
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 4) step2_assembly" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
perl /datapool/pipeline/metagenome/bin/dsub_gcy2.pl -thd 10 -mem 4 /datapool/pipeline/metagenome/upload/test_data/shell/step2_assembly.sh 
date +"%D %T -> Finish 4) step2_assembly" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step3_annotation
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 4) step3_annotation" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
perl /datapool/pipeline/metagenome/bin/dsub_gcy2.pl -thd 10 -mem 4 /datapool/pipeline/metagenome/upload/test_data/shell/step3_annotation.sh 
date +"%D %T -> Finish 4) step3_annotation" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step4_gtdbtk
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 4) step4_gtdbtk" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
nohup sh /datapool/pipeline/metagenome/upload/test_data/shell/step4_gtdbtk.sh >/datapool/pipeline/metagenome/upload/test_data/shell/step4_gtdbtk.log 
date +"%D %T -> Finish 4) step4_gtdbtk" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step0_rawdata
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 0) step0_rawdata" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
perl /datapool/pipeline/metagenome/bin/dsub_gcy2.pl -thd 10 -mem 4 /datapool/pipeline/metagenome/upload/test_data/shell/step0_rawdata.sh 
date +"%D %T -> Finish 0) step0_rawdata" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step1_cleandata
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 1) step1_cleandata" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
perl /datapool/pipeline/metagenome/bin/dsub_gcy2.pl -thd 10 -mem 4 /datapool/pipeline/metagenome/upload/test_data/shell/step1_cleandata.sh 
date +"%D %T -> Finish 1) step1_cleandata" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step2_assembly
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 4) step2_assembly" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
perl /datapool/pipeline/metagenome/bin/dsub_gcy2.pl -thd 10 -mem 4 /datapool/pipeline/metagenome/upload/test_data/shell/step2_assembly.sh 
date +"%D %T -> Finish 4) step2_assembly" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step3_annotation
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 4) step3_annotation" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
perl /datapool/pipeline/metagenome/bin/dsub_gcy2.pl -thd 10 -mem 4 /datapool/pipeline/metagenome/upload/test_data/shell/step3_annotation.sh 
date +"%D %T -> Finish 4) step3_annotation" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
### step5_stats_count
cd /datapool/pipeline/metagenome/upload/test_data/shell
date +"%D %T -> Start 5) step5_stats_count" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
nohup sh /datapool/pipeline/metagenome/upload/test_data/shell/step5_stats_count.sh >/datapool/pipeline/metagenome/upload/test_data/shell/step5_stats_count.log 
date +"%D %T -> Finish 4) step5_stats_count" >>/datapool/pipeline/metagenome/upload/test_data/shell/log
