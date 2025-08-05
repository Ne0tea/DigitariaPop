#!/usr/bin/perl

die("Usage: perl $0 <tophat ref genome> <RNAseq ID list> <RNAseq path> <pbs template>") if @ARGV != 4;

open FH_IN, "$ARGV[1]";
my $genome = $ARGV[0];
my $path = $ARGV[2];
my $pbs = $ARGV[3];

my @ID;
while(<FH_IN>){
	chomp;
	push @ID, $_;
}

foreach my $id (@ID){
	print "Solving $id\n";
	`qsub -N "RNAseq_$id" -o "$id.log" -v genome=$genome,path=$path,sp_id=$id $pbs`;
}

#perl 07_tophat_cufflinks.pl /public4/home/huangyj/Digitaria/assembly/04gap_filler/Dsan_chr.V2.fa rna_list /public4/home/huangyj/Digitaria/Dsan_ref_data/RNA-seq/ _07_RNAseq_template.pbs 
