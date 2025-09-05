#!/usr/bin/perl
use strict;
use warnings;
my $usage=<<USAGE;
perl $0 <TAB> > <TranTab>
USAGE
die $usage unless @ARGV == 1;
open IN,"<",$ARGV[0] or die $!;
my $fl = <IN>;
chomp $fl;
my @arr = split /\t/,$fl;
my $tit = shift @arr;
my @data;
my @head;
while (<IN>) {
    chomp;
    my @a = split /\t/;
    my $lab = shift @a;
    push @head,$lab;
    push @data,\@a;
}
close IN;
print join("\t",$tit,@head),"\n";
for my $i (0..$#arr) {
    print $arr[$i];
    foreach my $j (0..$#data) {
	print "\t${$data[$j]}[$i]";
    }
    print "\n";
}
exit 0;
__END__
