#!/usr/bin/perl
##Author: Jia Lei

die("Usage: perl $0 <GFF3 file> <evidence file>") if @ARGV != 2;

open FH_IN_GFF, "$ARGV[0]";
open FH_IN_EV, "$ARGV[1]";
open FH_OUT_L, ">abinitio_low_confidence_gene.list";
open FH_OUT_H, ">abinitio_high_confidence_gene.list";

my %gene_infor;
while(<FH_IN_GFF>){
	chomp;
	my @tmp = split/\t/;
	if($tmp[2] =~ /mRNA/){
		if($tmp[-1] =~ /ID\=(.+?)\;Parent/){
			$gene = $1;
			print $tmp[-1] unless $1;
			$gene_infor{"$tmp[0]\t$tmp[3]\t$tmp[4]"} = $gene;
		}
	}
}

my %ab_type = (	#三种从头预测类型，无需改动
	'Augustus' => 1,
	'AUGUSTUS' => 1,
	'Fgenesh' => 1,
	'GeneMark.hmm' => 1,
);

while(<FH_IN_EV>){
	chomp;
	my @tmp = split/\t/;
	my $loc = "$tmp[0]\t$tmp[1]\t$tmp[2]";
	my($ab_score, $score);
	for(my $i = 3; $i < @tmp; $i++){
		my $type = $tmp[$i];
		if($ab_type{$type} == 1){
			$ab_score ++;
		}elsif($type ne '-'){
			$score ++;
		}
	}
	if($score >= 1 and $gene_infor{$loc}){
		print FH_OUT_H "$gene_infor{$loc}\tAll_keep\n";
	}elsif($score == 0 and $ab_score >= 2 and $gene_infor{$loc}){
		print FH_OUT_H "$gene_infor{$loc}\tNeed_remove\n";
	}else{
		print FH_OUT_L "$gene_infor{$loc}\n";
	}
}

