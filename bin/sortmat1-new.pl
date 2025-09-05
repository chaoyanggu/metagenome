#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 <input.mat> > <sort.mat>\n" unless @ARGV == 1;
open IN,"<",$ARGV[0] or die $!;
my %best;
my %info;
my $fl = <IN>;
print $fl;
while (<IN>) {
    chomp;
    my $l = $_;
    $l =~ s/\t+$//;
    my @a = split /\t/,$l;
    my $id = shift @a;
    my $max = 0;
    foreach (@a) {
	$max = $_ if ($_ > $max);
    }
    $best{$id} = $max;
    $info{$id} = $l;
}
close IN;
foreach (sort {$best{$b}<=>$best{$a}} keys %best) {
    next if $_ eq 'others';
    print "$info{$_}\n";
}
if ($info{others}) {
    print "$info{others}\n";
}
exit 0;
__END__
