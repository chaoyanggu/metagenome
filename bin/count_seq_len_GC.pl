#! /usr/bin/perl
use strict;
use warnings;

if (@ARGV!=1) {
	print "Usage:perl $0 <seq_file>\n";
	exit;
}
my $seq_file=$ARGV[0];
my (%seq,$total_len,$total_GC);
$/="\>";
open (IN,$seq_file) || die;
<IN>;
while (<IN>) {
	my ($id,$seq)=split /\n/,$_,2;
	$id=(split /\s+/,$id)[0];
	$seq=~s/[\r\n\>]//g;
	my $GC=countNT("G",$seq)+countNT("C",$seq);
	$total_len+=length($seq);
	$total_GC+=$GC;
	print $id."\t";
	print length($seq);
	print "\t";
	printf "%.2f",$GC/length($seq)*100;
	print "\n";
}

if ($total_len>0) {
	print "Total"."\t".$total_len."\t";
	printf "%.2f",$total_GC/$total_len*100;
	print"\n";
}


sub countNT {
	my ($nt,$seq)=@_;
	return $seq=~s/$nt/$nt/gi;
}
