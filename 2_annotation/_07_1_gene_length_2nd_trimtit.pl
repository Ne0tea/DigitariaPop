#!/usr/bin/perl
#获取蛋白长度, 去除 gene ID 后续无效信息

die("Usage: perl $0 <prot file>\n") if @ARGV != 1;

open FH_IN, "$ARGV[0]";
open FH_OUT1, ">$ARGV[0].len";
open FH_OUT2, ">trimtit.$ARGV[0]";

my(@titles, %seq);
while(<FH_IN>){
	chomp;
	if(/>(.*?)\s+/){
		push @titles, $1;
	}else{
		$seq{$1} .= $_;	
	}	
}

foreach $title (@titles){
	$seqlength = length $seq{$title};
	$seq{$title} =~ s/\*//g;
	print FH_OUT1 "$title\t$seqlength\n";
	print FH_OUT2 ">$title\n$seq{$title}\n";
}
