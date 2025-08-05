#!/usr/bin/perl

die("Usage: perl $0 <gff3 file>") if @ARGV != 1;

open FH_IN, "$ARGV[0]";
open FH_OUT, ">ordered.$ARGV[0]";

my(%count, $output, $num);
while(<FH_IN>){
	unless(/^$/){
		my @tmp = split/\t/;
		$count{$tmp[0]}++ if $tmp[2] eq 'gene';
		s/=(.*?)$tmp[0]\.\d+/=$tmp[0].$count{$tmp[0]}/g;
		if($tmp[2] eq 'mRNA'){
			s/ID=(.*?);/ID=$1.mRNA1;/;
			s/Name=(.*)/Name=$1.mRNA1/;
		}elsif($tmp[2] eq 'exon'){
			s/Parent=(.*)/Parent=$1.mRNA1/;
			m/exon(\d+)/;
			$num = $1;
		}elsif($tmp[2] eq 'CDS'){
			s/ID=(.*?);/ID=$1.CDS$num;/;
			s/Parent=(.*)/Parent=$1.mRNA1/;
		}
	}
	$output .= $_;
}

print FH_OUT "$output";
