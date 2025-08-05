#!/usr/bin/perl
##Author: Jia Lei
##过滤基因集中重复序列基因、短基因

die("Usage: perl $0 <repeat list> <invalid list(length)> <gff3> <output file>") if @ARGV != 4;

open FH_IN_R, "$ARGV[0]" or die;
open FH_IN_L, "$ARGV[1]" or die;
open GFF3, "$ARGV[2]" or die;
open FH_OUT, ">$ARGV[3]" or die;

#my(%repeat, %invalid, %high, %keep);
my(%repeat, %invalid, %highu);
while(<FH_IN_R>){
	chomp;
	$repeat{$_}++;
}

while(<FH_IN_L>){	#长度不足基因
	chomp;
	$invalid{$_} ++;
}

#while(<FH_IN_H>){
#	chomp;
#	my($id, $judge) = split/\t/;
#	$high{$id} ++;
#	$keep{$id} = $judge;
#}

while(<GFF3>){
    chomp;
    my @tmp = split/\t/;
    if($tmp[2] =~ /mRNA/){
        if($tmp[-1] =~ /ID\=(.+?)\;Parent/){
            $high{$1}++;
        }
    }
}

my($cinv, $crepeat);
foreach my $name(sort keys %high){
	if($invalid{$name} == 1){
		$cinv++;
		next;
	}
	if($repeat{$name} == 1){
		$crepeat++;
		next;
	}
	print FH_OUT "$name\n";
}

print "Filter stat:\nInvalid: $cinv\nRepeat: $crepeat\n";
