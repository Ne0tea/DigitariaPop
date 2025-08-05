#!/usr/bin/perl
##Author: Jia Lei
##根据提供的基因列表提取GFF3文件中基因

die("Usage: perl $0 <reserved gene list > <GFF3>") if @ARGV != 2;

open FH_IN_L, "$ARGV[0]" or die;
open FH_IN_G, "$ARGV[1]" or die;
open FH_OUT, ">filtered.$ARGV[1]";

my %genes;	#genes存放保留geneID
while(<FH_IN_L>){
	chomp;
	$genes{$1}++ if /evm\.model\.(.+)/;
}

my $output;
while(<FH_IN_G>){
	chomp;
	next if /^$/;	#如果第一行为空行, 跳过
	my @tmp = split/\t/;
	if($tmp[2] eq 'gene'){	#判断一整个基因信息
		(my $ID = $tmp[-1]) =~ s/ID=evm\.TU\.|;.*//g;
		if($genes{$ID}){	#保留整个基因结构
			$output .= "$_\n";
			while(<FH_IN_G>){
				$output .= $_;
				last if /^$/;
			}
		}else{	#属于去除的基因, 去除整个基因结构
			while(<FH_IN_G>){
				last if /^$/;
			}
		}
	}else{
		print "Logicality Error, Check script or GFF3 file!\n";
	}
}

print FH_OUT "$output";
