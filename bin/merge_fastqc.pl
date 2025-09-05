#!/usr/bin/perl
use strict;
use warnings;

unless (@ARGV) {
    die "Usage: $0 <input.fasta>\n";
}
my $input_file = $ARGV[0];

open(my $fh, '<', $input_file) or die "Cannot open file $input_file: $!";

my %sample_data;

while (<$fh>) {
    chomp;
    next if /^Sample/;  # 跳过表头

    # 解析每一行
    my ($sample, $percent_dup, $percent_gc, $avg_len, $percent_fails, $total_seqs) = split /\t/;
    
    # 提取样本名称（去掉.R1/.R2）
    my $base_sample = $sample;
    $base_sample =~ s/(?:\.R|_)[12]$//;
    
    # 存储数据
    push @{$sample_data{$base_sample}{'duplicates'}}, $percent_dup;
    push @{$sample_data{$base_sample}{'gc'}}, $percent_gc;
    push @{$sample_data{$base_sample}{'length'}}, $avg_len;
    push @{$sample_data{$base_sample}{'sequences'}}, $total_seqs;
}

 # 打印表头
print "Sample\tTotal Sequences\tAverage sequence length\tTotal Bases\tGC content\tDuplicates\n";

 # 处理每个样本
foreach my $sample (sort keys %sample_data) {
    my $data = $sample_data{$sample};

    # 计算统计值
    my $total_seqs = sum(@{$data->{'sequences'}});
    my $avg_length = average(@{$data->{'length'}});
    my $total_bases = $total_seqs * $avg_length;
    my $avg_gc = average(@{$data->{'gc'}});
    my $avg_dup = average(@{$data->{'duplicates'}});

    # 输出结果
    print join("\t", 
        $sample,
        $total_seqs,
        sprintf("%.1f", $avg_length),
        $total_bases,
        sprintf("%.1f", $avg_gc),
        sprintf("%.2f", $avg_dup)
    ), "\n";
}

# 辅助函数：计算平均值
sub average {
    my @values = @_;
    my $sum = 0;
    $sum += $_ for @values;
    return $sum / scalar(@values);
}

# 辅助函数：求和
sub sum {
    my @values = @_;
    my $sum = 0;
    $sum += $_ for @values;
    return $sum;
}
