#!/usr/bin/perl


=head1 Description
	parallelly run by dsub
=head1 Usage
	perl dsub_batch.pl [options] <shell file>
		--mem maximum memory to use [Gb], default 2
		--thd maximum thread to run, default 2 in each node;
		--queue Queue name, default: batch
		--walltime Maximum walltime [HH:MM:SS], default: 10000:00:00
		--env the environment of software, current default base
		--h show help
=head1 Example

=head1 Version
	Author: Du Pengcheng dupengcheng@icdc.cn
	Date: 2015-12-15

=cut

use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ($mem, $thd, $queue, $walltime, $env, $help);
GetOptions (
	"mem:s"=>\$mem,
	"thd:s"=>\$thd,
	"queue:s"=>\$queue,
	"walltime:s"=>\$walltime,
	"env:s"=>\$env,
	"h"=>\$help,
);
die `pod2text $0` if ($help || @ARGV==0);

my $sh = $ARGV[0];
$sh = abs_path($sh);
my $cur_dir = abs_path("./");

$mem = 10 unless ($mem);
$thd = 1 unless ($thd);
$queue = 'batch' unless ($queue);
$walltime = 10000 unless ($walltime);
$env = "base" unless ($env);


my $master_script = $ARGV[0];
die "Error: Input file '$master_script' does not exist!\n" unless -e $master_script;

my $timestamp = sub_format_datetime(localtime(time()));
my $shell_dir = $cur_dir."/".basename($sh)."\.".$timestamp;

mkdir $shell_dir or die "Cannot create output directory: $!\n";

#读取并分割成子程序
open my $in, '<', $master_script or die "Cannot open $master_script: $!";

my @lines = <$in>;
close $in;

# 查找第一个###rm标记的行号
my $split_line = -1;
for (my $i = 0; $i < @lines; $i++) {
    if ($lines[$i] =~ /^###rm/) {
        $split_line = $i;
        last;  # 只取第一个匹配到的行号
    }
}

# 如果没有找到###rm，直接退出
die "Error: No ###rm marker found in the script.\n" unless $split_line >= 0;

my $lines_per_job = $split_line + 1 + 1;

# # 分割任务
my @jobs;
for (my $i = 0; $i < @lines; $i += $lines_per_job) {
    my $end = ($i + $lines_per_job - 1 < @lines) ? $i + $lines_per_job - 1 : $#lines;
    push @jobs, join("", @lines[$i..$end]);
}

# 提交每个任务
my $job_count = 0;
my @job_ids;

foreach my $job (@jobs) {
    $job_count++;
    my $job_file = sprintf("$shell_dir/work_%03d.sh", $job_count);

    open my $out, '>', $job_file or die "Cannot create $job_file: $!";
    print $out "#!/bin/bash\n";
    print $out "#PBS -q $queue\n";
    print $out "#PBS -l mem=${mem}gb,nodes=1:ppn=$thd,walltime=$walltime\n\n";
    print $out "source activate $env\n\n";
    print $out $job;
    close $out;

    # 提交作业
    my $qsub_cmd = "cd $shell_dir && qsub $job_file";
    my $qsub_output = `$qsub_cmd`;
    if ($qsub_output =~ /(\d+)\.\S+/) {
        push @job_ids, $1;
        print "Submitted work $job_count with ID $1\n";
    } else {
        warn "Failed to submit work $job_count\n";
    }
}

# 监控作业运行
if (@job_ids) {
    print "\nWaiting for works to complete...\n";
    my %active_jobs = map { $_ => 1 } @job_ids;

    while (keys %active_jobs) {
        sleep 30;
        my $qstat = `qstat`;

        foreach my $id (keys %active_jobs) {
            delete $active_jobs{$id} unless $qstat =~ /\b$id\b/;
        }

        printf "Works remaining: %d\r", scalar keys %active_jobs;
    }

    print "\nAll works completed successfully!\n";
}

exit 0;

sub sub_format_datetime
{
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
        $wday = $yday = $isdst = 0;
        sprintf("%4d-%02d-%02d.%02d.%02d.%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

