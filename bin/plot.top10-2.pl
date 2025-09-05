#!/usr/bin/perl
use strict;
use warnings;

BEGIN {
        unshift(@INC,"/datapool/stu/yuejl/bin/dupc/SVG_script/");
}

use SVG;;


## 检查参数
die "Usage: $0 <input.txt> <output.svg>\n" unless @ARGV == 2;
my ($input_file, $output_file) = @ARGV;

## 读取输入数据
open my $fh, '<', $input_file or die "Cannot open $input_file: $!";

## 读取表头（细菌分类）
my $header = <$fh>;
chomp $header;
my @taxa = split /\t/, $header;
shift @taxa;  # 去掉第一列的"Taxonomy"

## 读取样本数据
my (@samples, @data);
while (<$fh>) {
    chomp;
    my @fields = split /\t/;
    push @samples, shift @fields;  # 样本名
    push @data, [@fields];         # 相对丰度值
}
close $fh;

## 创建SVG图形
my $num_samples = scalar @samples;
my $bar_width = 40;
my $space_between_bars = 30;

# 动态计算图形宽度
my $graph_width = $num_samples * ($bar_width + $space_between_bars) + $space_between_bars;

my $graph_height = 450;
my $margin = 50;

# 计算图例所需宽度（基于最长的分类名称）
my $longest_taxon = 0;
foreach my $taxon (@taxa) {
    $longest_taxon = length($taxon) if length($taxon) > $longest_taxon;
}
my $legend_width = 30 + $longest_taxon * 5;  # 动态计算图例宽度
$legend_width = 150 if $legend_width < 150;  # 最小宽度
$legend_width = 400 if $legend_width > 400;  # 最大宽度


my $total_width = $graph_width + 2 * $margin + $legend_width; # 额外空间给图例
my $total_height = $graph_height + 2 * $margin + 100; # 增加高度以容纳样本名

my $svg = SVG->new(
    width => $total_width,
    height => $total_height,
);

# 定义颜色方案
my @colors = ('#DC143C', '#0000FF', '#20B2AA', '#FFA500', '#9370DB', '#98FB98', '#F08080', '#1E90FF', '#7CFC00', '#808000', '#B2B2B2');

# # 计算坐标轴位置
my $x0 = $margin;
my $y0 = $margin;
my $x1 = $x0 + $graph_width;
my $y1 = $y0 + $graph_height;
#

# 添加标题
$svg->text(
    x => $x0 + $graph_width/2,
    y => $margin - 20,
    'text-anchor' => 'middle',
    style => {
        'font-size' => '16px',
        'font-weight' => 'bold',
    },
)->cdata('Annotation');



# 绘制X轴
$svg->line(
    x1 => $x0, y1 => $y1,
    x2 => $x1, y2 => $y1,
    style => {
        'stroke' => 'black',
        'stroke-width' => 2,
    }
);

# 添加X轴刻度线（缩短的刻度线）
for my $i (0..$#samples) {
    my $x_pos = $x0 + $space_between_bars + $i * ($bar_width + $space_between_bars) + $bar_width/2;
    $svg->line(
        x1 => $x_pos, y1 => $y1,
        x2 => $x_pos, y2 => $y1 + 5,  # 缩短刻度线长度（5像素）
        style => {
            'stroke' => 'black',
            'stroke-width' => 1,
        }
    );
}


# 绘制Y轴
$svg->line(
    x1 => $x0, y1 => $y0,
    x2 => $x0, y2 => $y1,
    style => {
        'stroke' => 'black',
        'stroke-width' => 2,
    }
);

# 添加Y轴标签
$svg->text(
    x => $x0-40,
    y => $y0 + $graph_height/2,
    'writing-mode' => 'tb',
    'text-anchor' => 'middle',
    transform => "rotate(90,$x0-40,$y0+$graph_height/2)",
    style => {
        'font-size' => '12px',
    },
)->cdata('Relative Abundance');

# 计算每个样本的位置
my %sample_pos;
my $current_x = $x0 + $space_between_bars;
for my $i (0..$#samples) {
    $sample_pos{$i} = $current_x;
    $current_x += $bar_width + $space_between_bars;
}

# 绘制堆叠柱状图
for my $i (0..$#samples) {
    my $current_y = $y1;
    for my $j (0..$#taxa) {
        my $value = $data[$i][$j] || 0;
        my $height = $value * $graph_height;
        $svg->rect(
            x => $sample_pos{$i},
            y => $current_y - $height,
            width => $bar_width,
            height => $height,
            style => {
                'fill' => $colors[$j % @colors],
                'stroke' => 'black',
                'stroke-width' => 0.5,
            }
        );
        $current_y -= $height;
    }
}



# 添加X轴标签（样本名）
for my $i (0..$#samples) {
   my $x_pos = $sample_pos{$i} + $bar_width/2; 
   $svg->text(
        x => $x_pos,
        y => $y1,
        'text-anchor' => 'end',
        transform => "rotate(-45,$x_pos,$y1)",
        style => {
            'font-size' => '10px',
            'dominant-baseline' => 'hanging',  # 上对齐
        },
    )->cdata($samples[$i]);
}




# 添加图例
my $legend_x = $x1 + 30;
my $legend_y = $y0;
for my $i (0..$#taxa) {
    $svg->rect(
        x => $legend_x,
        y => $legend_y + $i * 20,
        width => 15,
        height => 15,
        style => {
            'fill' => $colors[$i % @colors],
            'stroke' => 'black',
            'stroke-width' => 0.5,
        }
    );
    $svg->text(
        x => $legend_x + 20,
        y => $legend_y + $i * 20 + 12,
        style => {
            'font-size' => '10px',
        },
    )->cdata($taxa[$i]);
}

# 添加Y轴刻度
my $max_value = 1;
my $step = 0.2;

for (my $val = 0; $val <= $max_value; $val += $step) {
    my $y = $y1 - ($val * $graph_height);
    
# 绘制刻度线
    $svg->line(
        x1 => $x0 - 5,
        y1 => $y,
        x2 => $x0,
        y2 => $y,
        style => {
            'stroke' => 'black',
            'stroke-width' => 1,
        }
    );
    
# 添加刻度标签
     $svg->text(
        x => $x0 - 10,
        y => $y+4,
        'text-anchor' => 'end',
        style => {
            'font-size' => '10px',
        },
    )->cdata(sprintf("%.1f", $val));
}


# 保存SVG文件
open my $out, '>', $output_file or die "Cannot write to $output_file: $!";
print $out $svg->xmlify;
close $out;

