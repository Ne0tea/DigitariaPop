#!/usr/bin/perl

my ($blast_result,$gene_len);

if(@ARGV ne 3){
	die "Usage: perl $0 <blast_result> <gene_len> <gene or protein>";	
}else{
	($blast_result,$gene_len,$seq_type)	= @ARGV;
}

if($seq_type =~ /protein/){
	$identity_standard = 30;
}elsif($seq_type =~ /gene/){
	$identity_standard = 80;
}

open OUT, ">$blast_result.repeat_overlap";
my %gene2len;
open GL, $gene_len;
while(<GL>){
	chomp;
	$gene2len{$1} = $2 if /(.+?)\t(.+)/;	
}
close GL;

open BR, $blast_result;
my %basecount;
my %loci;
my %count;
my %genenames;
while(<BR>){
	chomp;
	my @tmp = split/\t/;
	my ($genename,$similar,$start,$end,$evalue) = @tmp[0,2,-6,-5,-2];
	$genenames{$genename}++;
	if($similar > $identity_standard and $evalue < 1e-5){
		for my $i ($start..$end){
			$count{$genename."\t".$i}++;
			if($count{$genename."\t".$i} == 1){
				push @{$loci{$genename}},$i;    		
			}
		}
	}
}
close BR;


foreach (sort keys%genenames){
	if($loci{$_}){
		my $overlap_base = scalar @{$loci{$_}};
		my $overlap_percent = $overlap_base/$gene2len{$_};
		#print OUT "$_\t$overlap_base\t$gene2len{$_}\t$overlap_percent\n" if $overlap_percent>=0.25;
		#print OUT "$_\t$overlap_base\t$gene2len{$_}\t$overlap_percent\n" if $overlap_percent>=0.2;
		print OUT "$_\t$overlap_base\t$gene2len{$_}\t$overlap_percent\n" if $overlap_percent>=0.1;
	}
}

close OUT;
