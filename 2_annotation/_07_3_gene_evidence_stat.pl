#!/usr/bin/perl
##Author: Jia Lei
##统计整理EVM整合后各个基因证据支持情况

die("Usage: perl $0 <evm out file> <scaffold num> <evm weight> <output file>") if @ARGV != 4;

open IN_EVM, "$ARGV[0]";	#evm注释输出
open IN_W, "$ARGV[2]";	#权重文件，获取所有注释类型
open OUT, ">>$ARGV[3]";	#注释支持证据输出

my $scafname = $ARGV[1];

my %EV_type;
while(<IN_W>){
	my @tmp = split/\h+/;
	$EV_type{$tmp[1]} = 1;
}

my(%EV_count, @position, %All_infor);
while(<IN_EVM>){
	chomp;
	next if /^#/;
	if(/^$/){
		next unless @position;	#跳过第一个空行
		my @pos_sort = sort {$a <=> $b} @position;
		my($start, $end) = ($pos_sort[0], $pos_sort[-1]);
		$All_infor{$start} .= "$scafname\t$start\t$end";
		foreach(sort keys %EV_type){
			if($EV_count{$_} > 0){
				$All_infor{$start} .= "\t$_";
			}else{
				$All_infor{$start} .= "\t-";
			}
		}
		$All_infor{$start} .= "\n";
		@position = ();
		%EV_count = ();
	}else{
		my @tmp = split/\t/;
		push @position, $tmp[0];
		push @position, $tmp[1];
		my @EV_record = split/\,/, $tmp[-1];
		foreach my $per (@EV_record){
			$EV_count{$1}++ if $per =~ /;(.+?)\}/;
		}
	}
}

foreach(sort {$a <=> $b} keys %All_infor){
	print OUT "$All_infor{$_}";
}
