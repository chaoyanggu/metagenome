#!/usr/bin/perl
use strict;
use warnings;
use JSON;

# 检查是否提供了输入文件
if (@ARGV == 0) {
    die "Usage: $0 <json_file1> [<json_file2> ...]\n";
}


# 打印表头
print join("\t",
    "Sample",
    "Reads(Before)", "Bases(Before)",
    "Reads(After)", "Bases(After)",
    "Q20(After)", "Q30(After)",
    "GC Content", "Dup Rate"
), "\n";


# 处理每个JSON文件
foreach my $file (@ARGV) {
    open my $fh, '<', $file or die "Cannot open $file: $!";
    local $/;
    my $json_text = <$fh>;
    close $fh;
    
    my $json = JSON->new;
    my $data = eval { $json->decode($json_text) } or next;
 
    # 从文件名提取样本名（去掉路径和扩展名）
    my $sample = $file;
    $sample =~ s/^.*\///;  # 去掉路径
    $sample =~ s/\.json$//; # 去掉.json扩展名
    
    my $before = $data->{summary}{before_filtering};
    my $after = $data->{summary}{after_filtering};
    my $dup_rate = $data->{duplication}{rate};
    
    # 添加数据到表格

    print join("\t",
        $sample,
        format_number($before->{total_reads}),
        format_number($before->{total_bases}),
        format_number($after->{total_reads}),
        format_number($after->{total_bases}),
        sprintf("%.2f%%", $after->{q20_rate} * 100),
        sprintf("%.2f%%", $after->{q30_rate} * 100),
        sprintf("%.2f%%", $after->{gc_content} * 100),
        sprintf("%.2f%%", $dup_rate * 100)
    ), "\n";
}


# 辅助函数：格式化大数字
sub format_number {
    my $num = shift;
    $num =~ s/(\d)(?=(\d{3})+(\D|$))/$1,/g;
    return $num;
}

