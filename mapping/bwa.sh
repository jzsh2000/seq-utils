#! /bin/bash
reference=$1
fq1=$2
fq2=$3
bam=$4
/WPS2/zhangxn/software/bwa mem -t 16 $reference $fq1 $fq2 > $bam.sam
/WPS2/zhangxn/software/samtools-1.1/samtools view -b -u -S $bam.sam >$bam.bam
/WPS2/zhangxn/software/samtools-1.1/samtools sort -@ 16 $outdir/$bam.bam $bam.sort
rm $bam.sam
rm $bam.bam
mv $bam.sort.bam $bam.bam
/WPS2/zhangxn/software/samtools-1.1/samtools view -F 4 $bam.bam | perl -ne 'BEGIN{my %h}@a=split;$t=substr $a[2],0,15;if($h{$t}){$h{$t}++}else{$h{$t}=1}END{foreach my $key (keys %h)
{print "$key\t$h{$key}\n"}}' >$bam.counts.txt
    /WPS2/zhangxn/software/samtools-1.1/samtools flagstat $bam.bam >$bam.infos.txt
