#!/usr/bin/perl

die("Usage: perl $0 <part prefix> <genome> <separate part>") if @ARGV != 3;

my($prefix, $genome, $part) = @ARGV;
my $output="${prefix}_combined_evm_out";

`rm -f "${output}.gff3"`;

print "### Combine all separate gff files\n";
for(my $nu = 1; $nu <= $part; $nu++){
	`cat $prefix.$nu/evm_out.gff3 >>$output.gff3`;
}

print "### Generate protain file\n";
`/public/software/apps/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl $output.gff3 $genome prot >$output.prot.fa`;
#`/public/software/apps/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl $output.gff3 $genome CDS >$output.CDS.fa`
