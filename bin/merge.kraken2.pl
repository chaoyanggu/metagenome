#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
#
## 用法：perl merge_kraken_by_level_strict.pl *.label
## 输入：Kraken结果文件（格式：分类层级\t计数）
## 输出：各层级的独立表格（门/纲/目/科/属/种）
#
my %data;          # 原始数据 {taxonomy}{sample} = count
my %level_data;    # 分级数据 {level}{taxonomy}{sample} = count
my %level_totals; # 记录每个样本所有层级的累加值（用于归一化）
my @samples;       # 样本顺序

my %level_pattern = (
    domain  => qr/^d__[^|]+$/,  # domain层级匹配所有d__开头的分类
    phylum  => qr/^d__Bacteria.*?\|p__[^|]*($|\|)/,  # 必须包含d__Bacteria和p__
    class   => qr/^d__Bacteria.*?\|c__[^|]*($|\|)/,   # 必须包含d__Bacteria和c__
    order   => qr/^d__Bacteria.*?\|o__[^|]*($|\|)/,   # 必须包含d__Bacteria和o__
    family  => qr/^d__Bacteria.*?\|f__[^|]*($|\|)/,   # 必须包含d__Bacteria和f__
    genus   => qr/^d__Bacteria.*?\|g__[^|]*($|\|)/,   # 必须包含d__Bacteria和g__
    species => qr/^d__Bacteria.*?\|s__[^|]*$/    # 必须包含d__Bacteria和s__
);

# 排除模式：用于确保不包含更高级别的标记
my %exclude_pattern = (
    phylum  => qr/(c__|o__|f__|g__|s__)/,
    class   => qr/(o__|f__|g__|s__)/,
    order   => qr/(f__|g__|s__)/,
    family  => qr/(g__|s__)/,
    genus   => qr/s__/,
    species => qr//,  # species是最后一级，无需排除
);

# 读取样本文件
foreach my $file (@ARGV) {
    my $sample = basename($file);
    $sample =~ s/\.label$//;  # 从文件名提取样本名
    push @samples, $sample;
    
    $level_totals{$sample} = 0;  # 初始化样本总计数

    open(my $fh, '<', $file) or die "Cannot open $file: $!";
    while (<$fh>) {
        chomp;
        my ($taxonomy, $count) = split /\t/;
        $data{$taxonomy}{$sample} = $count;
        $level_totals{$sample} += $count;  # 累加样本总计数
        
        # 优先处理species（最后一级）
        if ($taxonomy =~ $level_pattern{species}) {
            $level_data{species}{$taxonomy}{$sample} = $count;
            $level_totals{species}{$sample} += $count;
            next;  # 跳过后续层级检查
        }        

        # 其他层级只处理细菌数据
        if ($taxonomy =~ /^d__Bacteria/) {
            foreach my $level (qw(phylum class order family genus)) {
               if (($taxonomy =~ $level_pattern{$level}) && ($taxonomy !~ $exclude_pattern{$level})) {
                    $level_data{$level}{$taxonomy}{$sample} = $count;
                    $level_totals{$level}{$sample} += $count;
                    last;  # 匹配到第一个符合条件的层级后退出
                }
            }
        }
        # Domain层级处理所有数据
        if ($taxonomy =~ $level_pattern{domain}) {
               $level_data{domain}{$taxonomy}{$sample} = $count;
               $level_totals{domain}{$sample} += $count;
        }
   
   }

    close($fh);
}

# 为每个层级生成独立表格
foreach my $level (qw(domain phylum class order family genus species)) {
    next unless exists $level_data{$level};  # 跳过未出现的层级
 

    # 绝对丰度表格
    my $abs_file = "${level}_level_absolute.txt";
    open(my $abs_out, '>', $abs_file) or die "Cannot write $abs_file: $!";
    print $abs_out join("\t", "Taxonomy", @samples), "\n";
    
    # 相对丰度表格
    my $rel_file = "${level}_level_relative.txt";
    open(my $rel_out, '>', $rel_file) or die "Cannot write $rel_file: $!";
    print $rel_out join("\t", "Taxonomy", @samples), "\n";
    
    foreach my $taxonomy (sort keys %{$level_data{$level}}) {
         print $abs_out $taxonomy;
         print $rel_out $taxonomy;
    
         foreach my $sample (@samples) {
              my $count = $level_data{$level}{$taxonomy}{$sample} || 0;
              print $abs_out "\t$count";
    
    #计算相对丰度
    my $total_count = $level_totals{$level}{$sample} || 0; 
    my $rel_abundance = $total_count >0 ? $count / $total_count :0;
    print $rel_out "\t",$rel_abundance;
         }
        print $abs_out "\n";
        print $rel_out "\n";
    }
    close ($abs_out);
    close ($rel_out);
    
    print "Generated: $abs_file\n";
    print "Generated: $rel_file\n";
}

