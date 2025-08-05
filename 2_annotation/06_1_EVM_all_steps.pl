#!/usr/bin/perl 
#Author: Jia Lei
#No transcript

die("Usage: perl $0 <Fasta file> <Separate num> <Outdir prefix> <evm template PBS> <weight file(abs path)> <GFF3 gene prediction> <GFF3 protein> <GFF3 transcripts>") if @ARGV != 8;

my($fasta, $sep_num, $out_pre, $pbs, $weight) = @ARGV;

my(@anno_type, %anno_saved);
for(my $nu = 5; $nu < @ARGV; $nu++){
	(my $gff = $ARGV[$nu]) =~ s/\.gff3$//;
	push @anno_type, $gff;
	open FH_IN_T, $ARGV[$nu];
	while(<FH_IN_T>){
		next if /^#/;
		my($contig,) = split/\t/;
		$anno_saved{$gff."\t".$contig} .= $_;
	}
}
print "Finish GFF3 infor saving\n";

open FH_IN, "$fasta";
my($title, $len_sum, %saved, @orders);
while(<FH_IN>){
	chomp;
	if(/^>/){
		($title = $_) =~ s/^>//;
		push @orders, $title;
	}else{
		$saved{$title} .= $_;
		$len_sum += length $_;
	}
}
my $sep_len = int($len_sum/$sep_num);
print "Finish Fasta file saving\nData len: $len_sum\nSeparate len: $sep_len\n";

my($tmp_len, $tmp_seq, $tmp_id, @id_list);
my $part = 1;
foreach my $id (@orders){
	$tmp_len += length $saved{$id};
	$tmp_seq .= ">$id\n$saved{$id}\n";
	$tmp_id .= "$id\n";
	push @id_list, $id;
	if($tmp_len >= $sep_len){	#进行输出
		my $tmp_dir = "$out_pre.$part";
		mkdir $tmp_dir;
		open FH_OUT_S, ">./$tmp_dir/SeparateFile.$part";
		open FH_OUT_C, ">./$tmp_dir/Contigs.$part";
		print FH_OUT_S $tmp_seq;
		print FH_OUT_C $tmp_id;
		foreach my $type(@anno_type){
			open FH_OUT_A, ">./$tmp_dir/$type.$part";
			foreach(@id_list){
				print FH_OUT_A $anno_saved{$type."\t".$_};
			}
			close FH_OUT_A;
		}
		`qsub -N "EVM_$out_pre" -o "log.$out_pre.$part" -v part="$out_pre.$part",weight="$weight",sep="SeparateFile.$part",abini="$anno_type[0].$part",homolog="$anno_type[1].$part",transcripts="$anno_type[2].$part" $pbs`;
		$tmp_len = 0;
		$tmp_seq = '';
		$tmp_id = '';
		$part ++;
		@id_list = ();
	}else{
		next;
	}
}
if($tmp_len != 0){
	my $tmp_dir = "$out_pre.$part";
	mkdir $tmp_dir;
	open FH_OUT1, ">./$tmp_dir/SeparateFile.$part";
	open FH_OUT2, ">./$tmp_dir/Contigs.$part";
	print FH_OUT1 $tmp_seq;
	print FH_OUT2 $tmp_id;
	foreach my $type(@anno_type){
		open FH_OUT3, ">./$tmp_dir/$type.$part";
		foreach(@id_list){
			print FH_OUT3 $anno_saved{$type."\t".$_};
		}
		close FH_OUT3;
	}
	`qsub -N "EVM_$out_pre" -o "log.$out_pre.$part" -v part="$out_pre.$part",weight="$weight",sep="SeparateFile.$part",abini="$anno_type[0].$part",homolog="$anno_type[1].$part",transcripts="$anno_type[2].$part" $pbs`;
}
print "$sep_num parts defined\n$part parts actually divided\n";
