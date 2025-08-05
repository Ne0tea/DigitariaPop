#!/usr/bin/perl 

if(@ARGV ne 3){
        die "Usage: perl $0 <Fa_file> <Separate_num> <outdir_prefix>\n";
}else{
        ($file,$separate_num,$prefix) = @ARGV;
}

open IN,"<$file";
###读取基因组，按照title储存title与序列信息，并且把title与第几个contig信息匹配
while(<IN>){
        chomp;
        if(/>(.+)/){
                $contig_num++;
                $tit = $1;
                $contigname2num{$tit} = $contig_num;  #因为后面要根据位置引用title，所以储存了两次哈希
                $contignum2name{$contig_num} = $tit;
                push @tits,$tit;
        }else{
                $allseq .= $_;
                $seq{$tit} .= $_;
        }
}

$genomelen = length $allseq;
print "Genome length: $genomelen\n";

$bases_per_file = int($genomelen/$separate_num); #用于下面基因组分割时的碱基计数
print "Bases per file: $bases_per_file\n";

$bases = 0;
$content = '';

$firstcontig = $tits[0];
$lastcontig = $tits[$#tits];

foreach $tit (@tits){

        $allbases += length($seq{$tit});

        $bases += length($seq{$tit});
        $content .= ">$tit\n$seq{$tit}\n";  #储存了分割后输出的序列信息
        if($bases >= $bases_per_file){

                $count_file++;
                $contig_flag{$tit} = $count_file;
                $last1contig = $tit;  #分割位置之前的最后一个contig
                $first1contig = $contignum2name{$contigname2num{$last1contig} + 1}; #分隔位置之后第二个文件的第一个contig

                mkdir "$prefix.$count_file";
                open OUT, '>'."./$prefix.$count_file/SeparateFile.$count_file";
                open OUT1, '>'."./$prefix.$count_file/Contigs.$count_file";
                print OUT "$content";
                print OUT1 "$tit\n";
                $content = '';
                $bases = 0;
        }
}

#统计最后一个分的文件（会有比较多的contigs）
mkdir "$prefix.$separate_num";
open OUT, '>'."./$prefix.$separate_num/SeparateFile.$separate_num";
open OUT1, '>'."./$prefix.$separate_num/Contigs.$separate_num";
$final_start_num = $contigname2num{$last1contig};  #该文件的第一个contig（来自上一个文件的最后一个contig）
for $i ($final_start_num..$#tits){
        $contigname = $tits[$i];
        print OUT1 "$contigname\n";
        print OUT ">$contigname\n$seq{$contigname}\n";
}
