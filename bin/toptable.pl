#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 <table> <lines> > out\n" unless @ARGV == 2;
open IN,"<",$ARGV[0] or die $!;
my $idx = 0;
my @others;
while (<IN>) {
    chomp;
    next unless $_;
    if ($idx > $ARGV[1]) {
	my @arr = split /\t/,$_;
	shift @arr;
	for (0..$#arr) {
	    $others[$_]+= $arr[$_];
	}
	next;
    }
    $idx++;
    print "$_\n";
}
if (@others) {
    print join("\t","Others",@others),"\n";
}
close IN;
__END__
