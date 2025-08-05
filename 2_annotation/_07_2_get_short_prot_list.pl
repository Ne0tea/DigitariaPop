#!/usr/bin/perl
#Author: Jia Lei
#基于蛋白长度进行过滤，去除基因ID保存到invalid_genelist.txt中，保留蛋白序列保存到.validprotseq为后缀的文件中。

die("Usage: perl $0 <protein file> <prot lenght>") if @ARGV != 2;

open FH_IN, $ARGV[0];
open FH_OUT1, ">invalid_genelist.txt";
open FH_OUT2,">$ARGV[0].validprotseq";

my $length = $ARGV[1];

my($title, %seq);
while(<FH_IN>){
	chomp;
	if(/>(\S+)/){
		$title = $1;
	}else{
		$seq{$title} .= $_;
	}
}

my($cstart, $clen, $cx);
foreach (sort keys %seq){
	$cstart ++ if $seq{$_} !~ /^M.*\*$/;
	$len = length $seq{$_};
	$clen ++ if $len < $length;
	$cx ++ if $seq{$_} =~ /X/;
	if($seq{$_} =~ /^M.*\*$/ and $len >= $length and $seq{$_} !~ /X/){
		print FH_OUT2 ">$_\n$seq{$_}\n";
	}else{
		print FH_OUT1 "$_\n";
	}
}

print "Filter stat:\nStart without M: $cstart\nShort length: $clen\nSeq with X: $cx\n";
