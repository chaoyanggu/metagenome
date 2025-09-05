#!/usr/bin/perl 
=head1 Description
	metagenome
=head1 Usage
	perl metagenome.pl [options]
	general input and output:
	-samplelist  <a file>       list of fastq from rawdata(file suffix can be ".fq.gz"), one strain per line, PE read seperated by ",", and different by "\n"
	-host               host[human, mice or none]
	-o   <a directory>  output dir, default current directory [./]
	-thd <num>          thread for dsub
	-h                  show help
	-notrun		    only write the shell, but not run
=head1 Example
	perl metagenome.pl  -samplelist sample.list -host human  -o ./ 

=head1 Version
        Author: guchaoyang0826@163.com
        Date: 2025-08-22
=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use Cwd qw/abs_path/;

my ($samplelist,$host,$outdir,$annotation,$thd,$help,$notrun);

GetOptions (
	"samplelist:s"=>\$samplelist,
	"host:s"=>\$host,
	"o:s"=>\$outdir,
	"h"=>\$help,
	"notrun"=>\$notrun
);
die `pod2text $0` if ($help || !$samplelist);

unless ($outdir) {
	$outdir="./";
}

unless (-e $outdir) {
	mkdir $outdir;
}

unless ($thd){
	$thd=8;
}

$outdir=abs_path($outdir);

my $shelldir="$outdir/shell";
my $fastqc_00="$outdir/00.rawdata";
my $qc_01="$outdir/01.cleandata";
my $assembly_02="$outdir/02.assembly";
my $annotation_03="$outdir/03.annotation";
my $result_04="$outdir/04.result";
for($shelldir,$fastqc_00,$qc_01,$assembly_02,$annotation_03,$result_04){
    (-s $_) || mkdir $_;
    $_ = abs_path($_);
}
#============================================================================================================================================================
if (!$samplelist){
	print "No input read files\n";
	exit;
}

if (!$host){
        print "No input host type: 'human' or 'mice' or 'none'\n";
        exit;
}

my $fastqc="/datapool/software/anaconda3/envs/qiime2/bin/fastqc";
my $fastp="/datapool/software/anaconda3/envs/RNAseq/bin/fastp";
my $bowtie2="/datapool/software/anaconda3/bin/bowtie2";
my $kraken2="/datapool/software/anaconda3/bin/kraken2";
my $metaspades="/datapool/software/anaconda3/bin/metaspades.py";
my $metawrap="/datapool/software/anaconda3/envs/metawrap/bin/metawrap";
my $checkm="/datapool/software/anaconda3/envs/qiime2/bin/checkm";
my $samtools="/datapool/software/anaconda3/bin/samtools";
my $seqkit="/datapool/software/anaconda3/envs/seqkit/bin/seqkit";
my $abricate="/datapool/software/anaconda3/bin/abricate";
my $prokka="/datapool/software/anaconda3/envs/prokka/bin/prokka";
my $gtdbtk="/datapool/software/anaconda3/envs/qiime2/bin/gtdbtk";

open (OUT,">$shelldir/metagenome.sh")||die;

my $sample;
my ($fqlist,$genomelist);
if ($samplelist && $host eq "human"){


open (OUT0,">$shelldir/step0_rawdata.sh")||die;
open (OUT1,">$shelldir/step1_cleandata.sh")||die;
open (OUT2,">$shelldir/step2_assembly.sh")||die;
open (OUT3,">$shelldir/step3_annotation.sh")||die;
open (OUT4,">$shelldir/step4_gtdbtk.sh")||die;
open (OUT5,">$shelldir/step5_stats_count.sh")||die;

	(-d "$fastqc_00/qc_result") || `mkdir "$fastqc_00/qc_result"`;
	(-d "$qc_01/nohostqc_result") || `mkdir "$qc_01/nohostqc_result"`;

	(-d "$assembly_02/assembly") || `mkdir "$assembly_02/assembly"`;
	(-d "$assembly_02/pick_1k_sequence") || `mkdir "$assembly_02/pick_1k_sequence"`;
	(-d "$assembly_02/metawrap") || `mkdir "$assembly_02/metawrap"`;
	
	(-d "$annotation_03/vfdb/result") ||`mkdir -p "$annotation_03/vfdb/result"`;
	(-d "$annotation_03/resfinder/result") ||`mkdir -p "$annotation_03/resfinder/result"`;
	(-d "$annotation_03/prokka") ||`mkdir "$annotation_03/prokka"`;
	(-d "$annotation_03/gtdbtk/fa") ||`mkdir -p "$annotation_03/gtdbtk/fa"`;
	
	(-d "$qc_01/kraken_count") ||`mkdir "$qc_01/kraken_count"`;

	open ($fqlist,">$qc_01/cleanfq.list") ||die;
	open ($genomelist,">$assembly_02/genomes.list") ||die;
	open (IN,$samplelist)||die;
	while(my $line=<IN>){
		chomp($line);
		$sample=(split/\t/,$line)[0];
		my $fq=(split/\t/,$line)[1];
		my $fq1=(split/\,/,$fq)[0];
		my $fq2=(split/\,/,$fq)[1];
                
		my $cmd_0="$fastqc --extract $fq1 -o $fastqc_00/qc_result\n$fastqc --extract $fq2 -o $fastqc_00/qc_result\n###rm\n";
                print OUT0 $cmd_0."\n";
		
		my $dir1="$qc_01/$sample";
		my $dir2="$assembly_02/assembly/$sample/";
		my $dir3="$assembly_02/metawrap/$sample/";
		(-d $dir1) || `mkdir $dir1`;
		(-d $dir2) || `mkdir $dir2`;
		(-d "$dir3/fq") || `mkdir -p $dir3/fq`;

	my $cmd_1="cd $qc_01/$sample\n$fastp  -i $fq1 -o $sample\_R1.clean.fq.gz  -I $fq2 -O $sample\_R2.clean.fq.gz -j $sample.json -h $sample.html\n";
	$cmd_1.="$bowtie2  --very-sensitive  -p 30 -x /datapool/db/hg38/hg38  -1 $sample\_R1.clean.fq.gz  -2 $sample\_R2.clean.fq.gz --al-conc-gz $sample.map.fq.gz  --un-conc-gz  $sample.unmap.fq.gz  -S $sample.sam 2> bowtie2.log\n$samtools view -F 4 -Sb $sample.sam  > $sample.bam\n$samtools sort $sample.bam  -o $sample.sorted.bam\nmv $sample.unmap.fq.1.gz $sample.unmap.1.fq.gz\nmv $sample.unmap.fq.2.gz $sample.unmap.2.fq.gz\n";
	$cmd_1.="$fastqc --extract $sample.unmap.1.fq.gz -o $qc_01/nohostqc_result\n$fastqc --extract $sample.unmap.2.fq.gz -o $qc_01/nohostqc_result\n";
	$cmd_1.="$kraken2 -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  $sample.unmap.1.fq.gz  $sample.unmap.2.fq.gz  --threads 20 --output $sample.kraken2 --use-names --report $sample.label --use-mpa-style\n###rm $sample.sam\n";
	print OUT1 $cmd_1."\n";
	
	my $cmd_2="$metaspades -t 60 -m 500 -1 $qc_01/$sample/$sample.unmap.1.fq.gz -2 $qc_01/$sample/$sample.unmap.2.fq.gz  --only-assembler  -o $assembly_02/assembly/$sample\n";
	$cmd_2.="perl $Bin/count_seq_len_GC.pl $assembly_02/assembly/$sample/contigs.fasta >$assembly_02/assembly/$sample/$sample.length.GC\n$seqkit  seq -m 1000 $assembly_02/assembly/$sample/contigs.fasta > $assembly_02/pick_1k_sequence/$sample.contig.fa\n";
	$cmd_2.="source activate metawrap\nexport PATH=\"/datapool/software/anaconda3/envs/metawrap/bin/\:\$PATH\"\nln -s $qc_01/$sample/$sample.unmap.1.fq.gz $assembly_02/metawrap/$sample/fq/$sample\_1.fastq\nln -s $qc_01/$sample/$sample.unmap.2.fq.gz $assembly_02/metawrap/$sample/fq/$sample\_2.fastq\n$metawrap binning -o $assembly_02/metawrap/$sample/ -t 20 -a $assembly_02/pick_1k_sequence/$sample.contig.fa --metabat2 --maxbin2 --concoct $assembly_02/metawrap/$sample/fq/$sample\_1.fastq $assembly_02/metawrap/$sample/fq/$sample\_2.fastq\n$metawrap bin_refinement -o $assembly_02/metawrap/$sample/bin_refinement -t 12 -A $assembly_02/metawrap/$sample/metabat2_bins/ -B $assembly_02/metawrap/$sample/maxbin2_bins/ -C $assembly_02/metawrap/$sample/concoct_bins/ -c 50 -x 10\nls $assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins/*.fa > $assembly_02/metawrap/$sample/$sample.MG.list\nless $assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins.stats  |sed '1d' |awk '\$2>90&&\$3<5' |awk -F \"\\t\" '{print \$1}'|perl -ne 'chomp;print qq($assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins/\$_.fa\\n)' >$assembly_02/metawrap/$sample/$sample.HG.list\nls $assembly_02/metawrap/$sample/$sample.MG.list>>$annotation_03/binning.list\n###rm\n";
	print OUT2 $cmd_2."\n";

	print $fqlist "$sample\t$qc_01/$sample/$sample.unmap.1.fq.gz\t$qc_01/$sample/$sample.unmap.2.fq.gz\n";
	print $genomelist "$sample\t$assembly_02/assembly/$sample/contigs.fasta\t$assembly_02/pick_1k_sequence/$sample.contig.fa\n";
	
	my $cmd_3="$abricate $assembly_02/pick_1k_sequence/$sample.contig.fa --db resfinder_new --minid=75 >$annotation_03/resfinder/result/$sample\_resfinder_new.tab\n$abricate $assembly_02/pick_1k_sequence/$sample.contig.fa --db vfdb_new --minid=75 >$annotation_03/vfdb/result/$sample\_vfdb_new.tab\nexport PATH=\"/datapool/software/anaconda3/envs/prokka/bin/\:\$PATH\"\n$prokka --force --cpus 10 --outdir $annotation_03/prokka/$sample --prefix $sample --locustag $sample --metagenome --kingdom Bacteria $assembly_02/pick_1k_sequence/$sample.contig.fa\n###rm\n";
        print OUT3 $cmd_3."\n";

	}
	my $cmd_4="cd $annotation_03/gtdbtk/fa\nless $annotation_03/binning.list |perl -ne 'chomp;print\"ln -s \$_ ./\\n\"'|sh\nexport GTDBTK_DATA_PATH=/datapool/software/anaconda3/envs/qiime2/share/gtdbtk-1.7.0/db/release202/\n$gtdbtk  classify_wf --genome_dir $annotation_03/gtdbtk/fa   --out_dir ./new_gtdb --extension fa --cpus 60";
	print OUT4 $cmd_4."\n";
	my $cmd_5="cd $fastqc_00/\n/datapool/software/anaconda3/envs/MetaPhage/bin/multiqc qc_result/ -o multiqc_output/ --filename multiqc_report.html\nperl $Bin/merge_fastqc.pl multiqc_output/multiqc_report_data/multiqc_general_stats.txt >multiqc.txt\nperl $Bin/sample_stats_extractor.pl $qc_01/*/*.json >$qc_01/sample.cleandata.stats\nls $qc_01/*/bowtie2.log|while read f;do grep 'overall alignment rate' \$f|awk -v f=\$f '{print f\"\\t\"\$1}' ;done |perl -ne 'chomp;\@a=split/\\t/,\$_;\$id=(split/\\//,\$a[0])[-2];print\"\$id\\t\$a[-1]\\n\"' >$qc_01/bowtie2.txt\n\ncd $qc_01/kraken_count\nperl $Bin/merge.kraken2.pl $qc_01/*/*.label\nls *relative.txt|sed 's/.txt//'|perl -ne 'chomp;print\"perl $Bin/sortmat1-new.pl \$_.txt > \$_.sorted.txt\\nperl $Bin/toptable.pl \$_.sorted.txt 10 > \$_.top10.txt\\nperl $Bin/trantab.pl \$_.top10.txt > \$_.top10.trans.txt\\nperl $Bin/plot.top10-2.pl \$_.top10.trans.txt \$_.top10.svg\\n/datapool/software/anaconda3/envs/DNBC4tools/bin/convert -density 300 \$_.top10.svg \$_.top10.png\\n\"' |sh\n\ncd $assembly_02/\nfind pick_1k_sequence -name \"*.fa\" | xargs $seqkit stats -a | awk 'NR==1 || FNR>1 {print}' >seqkit.txt\ncp $Bin/plot_density.R $assembly_02/pick_1k_sequence/\ncd $assembly_02/pick_1k_sequence/\n/datapool/software/anaconda3/envs/R4.1/bin/Rscript plot_density.R\n";
	print OUT5 $cmd_5."\n";

        close IN;
        close OUT0;
        close OUT1;
        close OUT2;
        close OUT3;
        close OUT4;
	close OUT5;

	 my $cmd="### step0_rawdata\ncd $shelldir\ndate +\"\%D \%T -> Start 0) step0_rawdata\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step0_rawdata.sh \ndate +\"\%D \%T -> Finish 0) step0_rawdata\" >>$shelldir/log\n";
        $cmd.="### step1_cleandata\ncd $shelldir\ndate +\"\%D \%T -> Start 1) step1_cleandata\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step1_cleandata.sh \ndate +\"\%D \%T -> Finish 1) step1_cleandata\" >>$shelldir/log\n";
        $cmd.="### step2_assembly\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step2_assembly\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step2_assembly.sh \ndate +\"\%D \%T -> Finish 4) step2_assembly\" >>$shelldir/log\n";
        $cmd.="### step3_annotation\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step3_annotation\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step3_annotation.sh \ndate +\"\%D \%T -> Finish 4) step3_annotation\" >>$shelldir/log\n";
 

       print OUT $cmd. "### step4_gtdbtk\ncd $shelldir\n",
	"date +\"\%D \%T -> Start 4) step4_gtdbtk\" >>$shelldir/log\n", 
	"nohup sh $shelldir/step4_gtdbtk.sh >$shelldir/step4_gtdbtk.log \n",
	"date +\"\%D \%T -> Finish 4) step4_gtdbtk\" >>$shelldir/log\n";

	print OUT $cmd. "### step5_stats_count\ncd $shelldir\n",
        "date +\"\%D \%T -> Start 5) step5_stats_count\" >>$shelldir/log\n",
        "nohup sh $shelldir/step5_stats_count.sh >$shelldir/step5_stats_count.log \n",
          "date +\"\%D \%T -> Finish 4) step5_stats_count\" >>$shelldir/log\n";
}

if ($samplelist && $host eq "mice"){


open (OUT0,">$shelldir/step0_rawdata.sh")||die;
open (OUT1,">$shelldir/step1_cleandata.sh")||die;
open (OUT2,">$shelldir/step2_assembly.sh")||die;
open (OUT3,">$shelldir/step3_annotation.sh")||die;
open (OUT4,">$shelldir/step4_gtdbtk.sh")||die;
open (OUT5,">$shelldir/step5_stats_count.sh")||die;

        (-d "$fastqc_00/qc_result") || `mkdir "$fastqc_00/qc_result"`;
        (-d "$qc_01/nohostqc_result") || `mkdir "$qc_01/nohostqc_result"`;

        (-d "$assembly_02/assembly") || `mkdir "$assembly_02/assembly"`;
        (-d "$assembly_02/pick_1k_sequence") || `mkdir "$assembly_02/pick_1k_sequence"`;
        (-d "$assembly_02/metawrap") || `mkdir "$assembly_02/metawrap"`;

        (-d "$annotation_03/vfdb/result") ||`mkdir -p "$annotation_03/vfdb/result"`;
        (-d "$annotation_03/resfinder/result") ||`mkdir -p "$annotation_03/resfinder/result"`;
        (-d "$annotation_03/prokka") ||`mkdir "$annotation_03/prokka"`;
        (-d "$annotation_03/gtdbtk/fa") ||`mkdir -p "$annotation_03/gtdbtk/fa"`;

        (-d "$qc_01/kraken_count") ||`mkdir "$qc_01/kraken_count"`;

        open ($fqlist,">$qc_01/cleanfq.list") ||die;
        open ($genomelist,">$assembly_02/genomes.list") ||die;
        open (IN,$samplelist)||die;
        while(my $line=<IN>){
                chomp($line);
                $sample=(split/\t/,$line)[0];
                my $fq=(split/\t/,$line)[1];
                my $fq1=(split/\,/,$fq)[0];
                my $fq2=(split/\,/,$fq)[1];

                my $cmd_0="$fastqc --extract $fq1 -o $fastqc_00/qc_result\n$fastqc --extract $fq2 -o $fastqc_00/qc_result\n###rm\n";
                print OUT0 $cmd_0."\n";

                my $dir1="$qc_01/$sample";
                my $dir2="$assembly_02/assembly/$sample/";
                my $dir3="$assembly_02/metawrap/$sample/";
                (-d $dir1) || `mkdir $dir1`;
                (-d $dir2) || `mkdir $dir2`;
                (-d "$dir3/fq") || `mkdir -p $dir3/fq`;

	my $cmd_1="cd $qc_01/$sample\n$fastp  -i $fq1 -o $sample\_R1.clean.fq.gz  -I $fq2 -O $sample\_R2.clean.fq.gz -j $sample.json -h $sample.html\n";
        $cmd_1.="$bowtie2  --very-sensitive  -p 30 -x /datapool/db/mm39  -1 $sample\_R1.clean.fq.gz  -2 $sample\_R2.clean.fq.gz --al-conc-gz $sample.map.fq.gz  --un-conc-gz  $sample.unmap.fq.gz  -S $sample.sam 2> bowtie2.log\n$samtools view -F 4 -Sb $sample.sam  > $sample.bam\n$samtools sort $sample.bam  -o $sample.sorted.bam\nmv $sample.unmap.fq.1.gz $sample.unmap.1.fq.gz\nmv $sample.unmap.fq.2.gz $sample.unmap.2.fq.gz\n";
        $cmd_1.="$fastqc --extract $sample.unmap.1.fq.gz -o $qc_01/nohostqc_result\n$fastqc --extract $sample.unmap.2.fq.gz -o $qc_01/nohostqc_result\n";
        $cmd_1.="$kraken2 -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  $sample.unmap.1.fq.gz  $sample.unmap.2.fq.gz  --threads 20 --output $sample.kraken2 --use-names --report $sample.label --use-mpa-style\n###rm $sample.sam\n";
        print OUT1 $cmd_1."\n";

        my $cmd_2="$metaspades -t 60 -m 500 -1 $qc_01/$sample/$sample.unmap.1.fq.gz -2 $qc_01/$sample/$sample.unmap.2.fq.gz  --only-assembler  -o $assembly_02/assembly/$sample\n";
        $cmd_2.="perl $Bin/count_seq_len_GC.pl $assembly_02/assembly/$sample/contigs.fasta >$assembly_02/assembly/$sample/$sample.length.GC\n$seqkit  seq -m 1000 $assembly_02/assembly/$sample/contigs.fasta > $assembly_02/pick_1k_sequence/$sample.contig.fa\n";
        $cmd_2.="source activate metawrap\nexport PATH=\"/datapool/software/anaconda3/envs/metawrap/bin/\:\$PATH\"\nln -s $qc_01/$sample/$sample.unmap.1.fq.gz $assembly_02/metawrap/$sample/fq/$sample\_1.fastq\nln -s $qc_01/$sample/$sample.unmap.2.fq.gz $assembly_02/metawrap/$sample/fq/$sample\_2.fastq\n$metawrap binning -o $assembly_02/metawrap/$sample/ -t 20 -a $assembly_02/pick_1k_sequence/$sample.contig.fa --metabat2 --maxbin2 --concoct $assembly_02/metawrap/$sample/fq/$sample\_1.fastq $assembly_02/metawrap/$sample/fq/$sample\_2.fastq\n$metawrap bin_refinement -o $assembly_02/metawrap/$sample/bin_refinement -t 12 -A $assembly_02/metawrap/$sample/metabat2_bins/ -B $assembly_02/metawrap/$sample/maxbin2_bins/ -C $assembly_02/metawrap/$sample/concoct_bins/ -c 50 -x 10\nls $assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins/*.fa > $assembly_02/metawrap/$sample/$sample.MG.list\nless $assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins.stats  |sed '1d' |awk '\$2>90&&\$3<5' |awk -F \"\\t\" '{print \$1}'|perl -ne 'chomp;print qq($assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins/\$_.fa\\n)' >$assembly_02/metawrap/$sample/$sample.HG.list\nls $assembly_02/metawrap/$sample/$sample.MG.list>>$annotation_03/binning.list\n###rm\n";
        print OUT2 $cmd_2."\n";

        print $fqlist "$sample\t$qc_01/$sample/$sample.unmap.1.fq.gz\t$qc_01/$sample/$sample.unmap.2.fq.gz\n";
        print $genomelist "$sample\t$assembly_02/assembly/$sample/contigs.fasta\t$assembly_02/pick_1k_sequence/$sample.contig.fa\n";

        my $cmd_3="$abricate $assembly_02/pick_1k_sequence/$sample.contig.fa --db resfinder_new --minid=75 >$annotation_03/resfinder/result/$sample\_resfinder_new.tab\n$abricate $assembly_02/pick_1k_sequence/$sample.contig.fa --db vfdb_new --minid=75 >$annotation_03/vfdb/result/$sample\_vfdb_new.tab\nexport PATH=\"/datapool/software/anaconda3/envs/prokka/bin/\:\$PATH\"\n$prokka --force --cpus 10 --outdir $annotation_03/prokka/$sample --prefix $sample --locustag $sample --metagenome --kingdom Bacteria $assembly_02/pick_1k_sequence/$sample.contig.fa\n###rm\n";
        print OUT3 $cmd_3."\n";

        }
        my $cmd_4="cd $annotation_03/gtdbtk/fa\nless $annotation_03/binning.list |perl -ne 'chomp;print\"ln -s \$_ ./\\n\"'|sh\nexport GTDBTK_DATA_PATH=/datapool/software/anaconda3/envs/qiime2/share/gtdbtk-1.7.0/db/release202/\n$gtdbtk  classify_wf --genome_dir $annotation_03/gtdbtk/fa   --out_dir ./new_gtdb --extension fa --cpus 60";
        print OUT4 $cmd_4."\n";
        my $cmd_5="cd $fastqc_00/\n/datapool/software/anaconda3/envs/MetaPhage/bin/multiqc qc_result/ -o multiqc_output/ --filename multiqc_report.html\nperl $Bin/merge_fastqc.pl multiqc_output/multiqc_report_data/multiqc_general_stats.txt >multiqc.txt\nperl $Bin/sample_stats_extractor.pl $qc_01/*/*.json >$qc_01/sample.cleandata.stats\nls $qc_01/*/bowtie2.log|while read f;do grep 'overall alignment rate' \$f|awk -v f=\$f '{print f\"\\t\"\$1}' ;done |perl -ne 'chomp;\@a=split/\\t/,\$_;\$id=(split/\\//,\$a[0])[-2];print\"\$id\\t\$a[-1]\\n\"' >$qc_01/bowtie2.txt\n\ncd $qc_01/kraken_count\nperl $Bin/merge.kraken2.pl $qc_01/*/*.label\nls *relative.txt|sed 's/.txt//'|perl -ne 'chomp;print\"perl $Bin/sortmat1-new.pl \$_.txt > \$_.sorted.txt\\nperl $Bin/toptable.pl \$_.sorted.txt 10 > \$_.top10.txt\\nperl $Bin/trantab.pl \$_.top10.txt > \$_.top10.trans.txt\\nperl $Bin/plot.top10-2.pl \$_.top10.trans.txt \$_.top10.svg\\n/datapool/software/anaconda3/envs/DNBC4tools/bin/convert -density 300 \$_.top10.svg \$_.top10.png\\n\"' |sh\n\ncd $assembly_02/\nfind pick_1k_sequence -name \"*.fa\" | xargs $seqkit stats -a | awk 'NR==1 || FNR>1 {print}' >seqkit.txt\ncp $Bin/plot_density.R $assembly_02/pick_1k_sequence/\ncd $assembly_02/pick_1k_sequence/\n/datapool/software/anaconda3/envs/R4.1/bin/Rscript plot_density.R\n";
        print OUT5 $cmd_5."\n";

	close IN;
        close OUT0;
        close OUT1;
        close OUT2;
        close OUT3;
        close OUT4;
        close OUT5;

         my $cmd="### step0_rawdata\ncd $shelldir\ndate +\"\%D \%T -> Start 0) step0_rawdata\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step0_rawdata.sh \ndate +\"\%D \%T -> Finish 0) step0_rawdata\" >>$shelldir/log\n";
        $cmd.="### step1_cleandata\ncd $shelldir\ndate +\"\%D \%T -> Start 1) step1_cleandata\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step1_cleandata.sh \ndate +\"\%D \%T -> Finish 1) step1_cleandata\" >>$shelldir/log\n";
        $cmd.="### step2_assembly\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step2_assembly\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step2_assembly.sh \ndate +\"\%D \%T -> Finish 4) step2_assembly\" >>$shelldir/log\n";
        $cmd.="### step3_annotation\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step3_annotation\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step3_annotation.sh \ndate +\"\%D \%T -> Finish 4) step3_annotation\" >>$shelldir/log\n";


       print OUT $cmd. "### step4_gtdbtk\ncd $shelldir\n",
        "date +\"\%D \%T -> Start 4) step4_gtdbtk\" >>$shelldir/log\n",
        "nohup sh $shelldir/step4_gtdbtk.sh >$shelldir/step4_gtdbtk.log \n",
          "date +\"\%D \%T -> Finish 4) step4_gtdbtk\" >>$shelldir/log\n";

        print OUT $cmd. "### step5_stats_count\ncd $shelldir\n",
        "date +\"\%D \%T -> Start 5) step5_stats_count\" >>$shelldir/log\n",
        "nohup sh $shelldir/step5_stats_count.sh >$shelldir/step5_stats_count.log \n",
          "date +\"\%D \%T -> Finish 4) step5_stats_count\" >>$shelldir/log\n";
}


if ($samplelist && $host eq "none"){


open (OUT0,">$shelldir/step0_rawdata.sh")||die;
open (OUT1,">$shelldir/step1_cleandata.sh")||die;
open (OUT2,">$shelldir/step2_assembly.sh")||die;
open (OUT3,">$shelldir/step3_annotation.sh")||die;
open (OUT4,">$shelldir/step4_gtdbtk.sh")||die;
open (OUT5,">$shelldir/step5_stats_count.sh")||die;

        (-d "$fastqc_00/qc_result") || `mkdir "$fastqc_00/qc_result"`;
        (-d "$qc_01/nohostqc_result") || `mkdir "$qc_01/nohostqc_result"`;

        (-d "$assembly_02/assembly") || `mkdir "$assembly_02/assembly"`;
        (-d "$assembly_02/pick_1k_sequence") || `mkdir "$assembly_02/pick_1k_sequence"`;
        (-d "$assembly_02/metawrap") || `mkdir "$assembly_02/metawrap"`;

        (-d "$annotation_03/vfdb/result") ||`mkdir -p "$annotation_03/vfdb/result"`;
        (-d "$annotation_03/resfinder/result") ||`mkdir -p "$annotation_03/resfinder/result"`;
        (-d "$annotation_03/prokka") ||`mkdir "$annotation_03/prokka"`;
        (-d "$annotation_03/gtdbtk/fa") ||`mkdir -p "$annotation_03/gtdbtk/fa"`;

        (-d "$qc_01/kraken_count") ||`mkdir "$qc_01/kraken_count"`;

        open ($fqlist,">$qc_01/cleanfq.list") ||die;
        open ($genomelist,">$assembly_02/genomes.list") ||die;
        open (IN,$samplelist)||die;
        while(my $line=<IN>){
                chomp($line);
                $sample=(split/\t/,$line)[0];
                my $fq=(split/\t/,$line)[1];
                my $fq1=(split/\,/,$fq)[0];
                my $fq2=(split/\,/,$fq)[1];

                my $cmd_0="$fastqc --extract $fq1 -o $fastqc_00/qc_result\n$fastqc --extract $fq2 -o $fastqc_00/qc_result\n###rm\n";
                print OUT0 $cmd_0."\n";

                my $dir1="$qc_01/$sample";
                my $dir2="$assembly_02/assembly/$sample/";
                my $dir3="$assembly_02/metawrap/$sample/";
                (-d $dir1) || `mkdir $dir1`;
                (-d $dir2) || `mkdir $dir2`;
                (-d "$dir3/fq") || `mkdir -p $dir3/fq`;


	my $cmd_1="cd $qc_01/$sample\n$fastp  -i $fq1 -o $sample\_R1.clean.fq.gz  -I $fq2 -O $sample\_R2.clean.fq.gz -j $sample.json -h $sample.html\n";
        $cmd_1.="$fastqc --extract $sample\_R1.clean.fq.gz  -o $qc_01/nohostqc_result\n$fastqc --extract $sample\_R2.clean.fq.gz  -o $qc_01/nohostqc_result\n";
        $cmd_1.="$kraken2 -db /datapool/db/Kraken2_db/minikraken2_v2_8GB_201904_UPDATE/ --paired --gzip-compressed  $sample\_R1.clean.fq.gz  $sample\_R2.clean.fq.gz  --threads 20 --output $sample.kraken2 --use-names --report $sample.label --use-mpa-style\n###rm $sample.sam\n";
        print OUT1 $cmd_1."\n";

        my $cmd_2="$metaspades -t 60 -m 500 -1 $qc_01/$sample/$sample\_R1.clean.fq.gz -2 $qc_01/$sample/$sample\_R2.clean.fq.gz  --only-assembler  -o $assembly_02/assembly/$sample\n";
        $cmd_2.="perl $Bin/count_seq_len_GC.pl $assembly_02/assembly/$sample/contigs.fasta >$assembly_02/assembly/$sample/$sample.length.GC\n$seqkit  seq -m 1000 $assembly_02/assembly/$sample/contigs.fasta > $assembly_02/pick_1k_sequence/$sample.contig.fa\n";
        $cmd_2.="source activate metawrap\nexport PATH=\"/datapool/software/anaconda3/envs/metawrap/bin/\:\$PATH\"\nln -s $qc_01/$sample/$sample\_R1.clean.fq.gz $assembly_02/metawrap/$sample/fq/$sample\_1.fastq\nln -s $qc_01/$sample/$sample\_R2.clean.fq.gz  $assembly_02/metawrap/$sample/fq/$sample\_2.fastq\n$metawrap binning -o $assembly_02/metawrap/$sample/ -t 20 -a $assembly_02/pick_1k_sequence/$sample.contig.fa --metabat2 --maxbin2 --concoct $assembly_02/metawrap/$sample/fq/$sample\_1.fastq $assembly_02/metawrap/$sample/fq/$sample\_2.fastq\n$metawrap bin_refinement -o $assembly_02/metawrap/$sample/bin_refinement -t 12 -A $assembly_02/metawrap/$sample/metabat2_bins/ -B $assembly_02/metawrap/$sample/maxbin2_bins/ -C $assembly_02/metawrap/$sample/concoct_bins/ -c 50 -x 10\nls $assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins/*.fa > $assembly_02/metawrap/$sample/$sample.MG.list\nless $assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins.stats  |sed '1d' |awk '\$2>90&&\$3<5' |awk -F \"\\t\" '{print \$1}'|perl -ne 'chomp;print qq($assembly_02/metawrap/$sample/bin_refinement/metawrap_50_10_bins/\$_.fa\\n)' >$assembly_02/metawrap/$sample/$sample.HG.list\nls $assembly_02/metawrap/$sample/$sample.MG.list>>$annotation_03/binning.list\n###rm\n";
        print OUT2 $cmd_2."\n";

        print $fqlist "$sample\t$qc_01/$sample/$sample.unmap.1.fq.gz\t$qc_01/$sample/$sample.unmap.2.fq.gz\n";
        print $genomelist "$sample\t$assembly_02/assembly/$sample/contigs.fasta\t$assembly_02/pick_1k_sequence/$sample.contig.fa\n";

        my $cmd_3="$abricate $assembly_02/pick_1k_sequence/$sample.contig.fa --db resfinder_new --minid=75 >$annotation_03/resfinder/result/$sample\_resfinder_new.tab\n$abricate $assembly_02/pick_1k_sequence/$sample.contig.fa --db vfdb_new --minid=75 >$annotation_03/vfdb/result/$sample\_vfdb_new.tab\nexport PATH=\"/datapool/software/anaconda3/envs/prokka/bin/\:\$PATH\"\n$prokka --force --cpus 10 --outdir $annotation_03/prokka/$sample --prefix $sample --locustag $sample --metagenome --kingdom Bacteria $assembly_02/pick_1k_sequence/$sample.contig.fa\n###rm\n";
        print OUT3 $cmd_3."\n";

        }
        my $cmd_4="cd $annotation_03/gtdbtk/fa\nless $annotation_03/binning.list |perl -ne 'chomp;print\"ln -s \$_ ./\\n\"'|sh\nexport GTDBTK_DATA_PATH=/datapool/software/anaconda3/envs/qiime2/share/gtdbtk-1.7.0/db/release202/\n$gtdbtk  classify_wf --genome_dir $annotation_03/gtdbtk/fa   --out_dir ./new_gtdb --extension fa --cpus 60";
        print OUT4 $cmd_4."\n";
        my $cmd_5="cd $fastqc_00/\n/datapool/software/anaconda3/envs/MetaPhage/bin/multiqc qc_result/ -o multiqc_output/ --filename multiqc_report.html\nperl $Bin/merge_fastqc.pl multiqc_output/multiqc_report_data/multiqc_general_stats.txt >multiqc.txt\nperl $Bin/sample_stats_extractor.pl $qc_01/*/*.json >$qc_01/sample.cleandata.stats\nls $qc_01/*/bowtie2.log|while read f;do grep 'overall alignment rate' \$f|awk -v f=\$f '{print f\"\\t\"\$1}' ;done |perl -ne 'chomp;\@a=split/\\t/,\$_;\$id=(split/\\//,\$a[0])[-2];print\"\$id\\t\$a[-1]\\n\"' >$qc_01/bowtie2.txt\n\ncd $qc_01/kraken_count\nperl $Bin/merge.kraken2.pl $qc_01/*/*.label\nls *relative.txt|sed 's/.txt//'|perl -ne 'chomp;print\"perl $Bin/sortmat1-new.pl \$_.txt > \$_.sorted.txt\\nperl $Bin/toptable.pl \$_.sorted.txt 10 > \$_.top10.txt\\nperl $Bin/trantab.pl \$_.top10.txt > \$_.top10.trans.txt\\nperl $Bin/plot.top10-2.pl \$_.top10.trans.txt \$_.top10.svg\\n/datapool/software/anaconda3/envs/DNBC4tools/bin/convert -density 300 \$_.top10.svg \$_.top10.png\\n\"' |sh\n\ncd $assembly_02/\nfind pick_1k_sequence -name \"*.fa\" | xargs $seqkit stats -a | awk 'NR==1 || FNR>1 {print}' >seqkit.txt\ncp $Bin/plot_density.R $assembly_02/pick_1k_sequence/\ncd $assembly_02/pick_1k_sequence/\n/datapool/software/anaconda3/envs/R4.1/bin/Rscript plot_density.R\n";
        print OUT5 $cmd_5."\n";

	close IN;
        close OUT0;
        close OUT1;
        close OUT2;
        close OUT3;
        close OUT4;
        close OUT5;

        my $cmd="### step0_rawdata\ncd $shelldir\ndate +\"\%D \%T -> Start 0) step0_rawdata\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step0_rawdata.sh \ndate +\"\%D \%T -> Finish 0) step0_rawdata\" >>$shelldir/log\n";
        $cmd.="### step1_cleandata\ncd $shelldir\ndate +\"\%D \%T -> Start 1) step1_cleandata\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step1_cleandata.sh \ndate +\"\%D \%T -> Finish 1) step1_cleandata\" >>$shelldir/log\n";
        $cmd.="### step2_assembly\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step2_assembly\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step2_assembly.sh \ndate +\"\%D \%T -> Finish 4) step2_assembly\" >>$shelldir/log\n";
        $cmd.="### step3_annotation\ncd $shelldir\ndate +\"\%D \%T -> Start 4) step3_annotation\" >>$shelldir/log\nperl $Bin/dsub_gcy2.pl -thd 10 -mem 4 $shelldir/step3_annotation.sh \ndate +\"\%D \%T -> Finish 4) step3_annotation\" >>$shelldir/log\n";

	print OUT $cmd. "### step4_gtdbtk\ncd $shelldir\n",
        "date +\"\%D \%T -> Start 4) step4_gtdbtk\" >>$shelldir/log\n",
        "nohup sh $shelldir/step4_gtdbtk.sh >$shelldir/step4_gtdbtk.log \n",
        "date +\"\%D \%T -> Finish 4) step4_gtdbtk\" >>$shelldir/log\n";
	
	print OUT $cmd. "### step5_stats_count\ncd $shelldir\n",
	"date +\"\%D \%T -> Start 5) step5_stats_count\" >>$shelldir/log\n",
	"nohup sh $shelldir/step5_stats_count.sh >$shelldir/step5_stats_count.log \n",
	"date +\"\%D \%T -> Finish 4) step5_stats_count\" >>$shelldir/log\n";
}


close OUT;
	
$notrun && exit;
