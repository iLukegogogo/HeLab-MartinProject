#! /usr/bin/perl

use strict;
my $exp_file      = $ARGV[0];


my $BOWTIE2_INDEX        = "/srv/persistent/keliu/genomes/mm10/mm10.index";
my $RT_GTF_ANNOTATION    = "./annotation/mm10.repbase.gtf";


############# make the necessary directory #########

my @dir_list =qw /tmp alignment fastq RData/;
grep {
    chomp $_;
    `mkdir $_` if not -e $_;
} @dir_list;


################## parse exp file to extract SRA id ########
my %h;
open(INFILE,"< $exp_file");
while(my $line=<INFILE>){
    chomp $line;
    while($line=~/([SE]RR\d+)/g){$h{$1}=1;}
}

close INFILE;
my @sra_no = keys %h;


#################### download data from SRA ###################
my @fastq_dump_cmd_list;
my @wget_cmd_list;
grep{
    my $line=$_;
    $line=~s/\r//g;
    chomp $line;
    if($line ne "" ){
        $line=~/(([SE]RR)\d\d\d)/;
        my   $wget_cmd=" wget  ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/$2/$1/$line/$line.sra";
        push @wget_cmd_list,      $wget_cmd."\n"                                 if  (($line!~/^#/) && (! -e "alignment/${line}.all.rmdup.bam"));
        push @fastq_dump_cmd_list,"fastq-dump  --split-files --gzip $line.sra\n" if  (($line!~/^#/) && (! -e "alignment/${line}.all.rmdup.bam"));
    }
} @sra_no;

chdir "fastq";

open (OUTFILE,"> wget_cmd_list");
print OUTFILE join("",@wget_cmd_list);
close OUTFILE;


open(OUTFILE,"> fastq_dump_cmd_list");
print OUTFILE join("",@fastq_dump_cmd_list);
close OUTFILE;


`parallel -j 1 --no-notice  < wget_cmd_list`;
`parallel -j 10 --no-notice < fastq_dump_cmd_list`;
`rm wget_cmd_list`;
`rm fastq_dump_cmd_list`;

chdir "../";

#################### fastq read mapping ####################
my $bowtie2_option = "--very-sensitive -p 8";

my $map_cmd;
my $picard_cmd;
my $filter_cmd;
my $unique_cmd;
my @map_cmd_list;
my @picard_cmd_list;
my @filter_cmd_list;
my @unique_cmd_list;

foreach my $line (@sra_no){
    chomp $line;
    next if $line eq "";
    next if $line=~/^#/;
    $line=~s/\r//;
    my $sra_no = $line;
    
    if ( -e "fastq/${sra_no}_2.fastq.gz" ){
        $map_cmd    = "bowtie2 $bowtie2_option --no-discordant --no-mixed -X 1000 -x $BOWTIE2_INDEX -1 fastq/${sra_no}_1.fastq.gz -2 fastq/${sra_no}_2.fastq.gz | samtools view -f 2 -bS - | samtools sort -o tmp/${sra_no}.all.bam -T sorted.${sra_no}.all \n";
        $picard_cmd = "java -jar /users/keliu/bin/picard/picard.jar MarkDuplicates INPUT=tmp/${sra_no}.all.bam  OUTPUT=tmp/${sra_no}.all.rmdup.bam METRICS_FILE=tmp/$sra_no.all.metrics REMOVE_DUPLICATES=true\n";
        $filter_cmd = "samtools view -b -F 1804    tmp/${sra_no}.all.rmdup.bam        > alignment/${sra_no}.all.rmdup.bam\n";
        $unique_cmd = "samtools view -b -q 10      alignment/${sra_no}.all.rmdup.bam  > alignment/${sra_no}.unique.rmdup.bam\n";
        
    }else{
        $map_cmd    = "bowtie2 $bowtie2_option  -x $BOWTIE2_INDEX -U fastq/${sra_no}_1.fastq.gz                                                              | samtools view      -bS - | samtools sort -o tmp/${sra_no}.all.bam -T sorted.${sra_no}.all  \n";
        $picard_cmd = "java -jar /users/keliu/bin/picard/picard.jar MarkDuplicates INPUT=tmp/${sra_no}.all.bam  OUTPUT=tmp/${sra_no}.all.rmdup.bam METRICS_FILE=tmp/$sra_no.all.metrics REMOVE_DUPLICATES=true\n";
        $filter_cmd = "samtools view  -b -F 1804   tmp/${sra_no}.all.rmdup.bam        > alignment/${sra_no}.all.rmdup.bam\n";
        $unique_cmd = "samtools view  -b -q 10     alignment/${sra_no}.all.rmdup.bam  > alignment/${sra_no}.unique.rmdup.bam\n";
    }

    if (! -e "alignment/${sra_no}.all.rmdup.bam"){
        push @map_cmd_list,$map_cmd;
        push @picard_cmd_list,$picard_cmd;
        push @filter_cmd_list,$filter_cmd;
    }
    if (! -e "alignment/${sra_no}.unique.rmdup.bam"){
        push @unique_cmd_list,$unique_cmd;
    }

}


open(OUTFILE,"> tmp/map_cmd_list");
print OUTFILE join("",@map_cmd_list);
close OUTFILE;

open(OUTFILE,"> tmp/picard_cmd_list");
print OUTFILE join("",@picard_cmd_list);
close OUTFILE;

open(OUTFILE,"> tmp/filter_cmd_list");
print OUTFILE join("",@filter_cmd_list);
close OUTFILE;


open(OUTFILE,"> tmp/unique_cmd_list");
print OUTFILE join("",@unique_cmd_list);
close OUTFILE;

`parallel --no-notice -j 8 < tmp/map_cmd_list`;
`parallel --no-notice -j 8 < tmp/picard_cmd_list`;
`parallel --no-notice -j 8 < tmp/filter_cmd_list`;


############### use featureCounts to do read counting ############
chdir 'alignment';
my $out_put       = `ls |grep all.rmdup.bam`;
my $bam_file      = join("  ",split /\n+/,$out_put);
my $featureCounts_cmd;


$featureCounts_cmd = "featureCounts -a $RT_GTF_ANNOTATION -B -C -p --primary -T 10 -O -M -g repeat -o ../RData/rt.read.count.matrix  $bam_file";
`$featureCounts_cmd`;



exit();
########### macs2 peak calling ###########
my @macs2_call_peak_cmd;
my @macs2_bdgcmp_cmd;
my @sort_cmd;
my @bigwig_cmd;
my @mv_peak_cmd;


my @sort_peak_cmd;
my @bigbed_cmd;


open(INFILE,"< $exp_file") || die "Can not find $exp_file\n";
<INFILE>;
while(my $line=<INFILE>){
    $line=~s/\r//;
    chomp $line;
    next if $line=~/^#/;
    next if $line eq "";
    
    my @tmp = split /,/,$line;
    my $ip       = $tmp[1];
    my $input    = $tmp[2];
    my $exp_name = $tmp[0];
    my $layout   = $tmp[3];
    my $target   = $tmp[4];
    my @ip_list     = map {"alignment/$_.unique.rmdup.bam"}  split /:/,$ip;
    my @input_list  = map {"alignment/$_.unique.rmdup.bam"}  split /:/,$input;
    $ip    = join(" ",@ip_list);
    $input = join(" ",@input_list);
    
    my $macs2_cmd;
    if($layout eq "SINGLE"){
        $macs2_cmd = " macs2 callpeak -t $ip -c $input --outdir macs2.output -g $SPECIES -n $exp_name -B --SPMR -q 0.05 --keep-dup all ";
        $macs2_cmd = $macs2_cmd . " --nomodel --extsize 150 " if $target eq "HISTONE";
        $macs2_cmd = $macs2_cmd . " \n";
    }else{
        $macs2_cmd = " macs2 callpeak -t $ip -c $input --outdir macs2.output -g $SPECIES -n $exp_name -B --SPMR -q 0.05 --keep-dup all  -f BAMPE \n";
    }
    my $ppois_bdg   = " macs2.output/${exp_name}_ppois.bdg";
    my $bigwig_file = " signal/${exp_name}_ppois.bigwig";

    my $bdgcmp_cmd = " macs2  bdgcmp -t macs2.output/${exp_name}_treat_pileup.bdg  -c macs2.output/${exp_name}_control_lambda.bdg -o $ppois_bdg -m ppois \n";
    my $sort_cmd   = " sort -k1,1 -k2,2n -o $ppois_bdg $ppois_bdg \n";
    my $bigwig_cmd = " bedGraphToBigWig $ppois_bdg  $GLOBAL_CHR_SIZE_FILE $bigwig_file \n";
    my $mv_peak_cmd  = " cat macs2.output/".$exp_name ."_peaks.*Peak | grep -v '#' |cut -f 1,2,3 > peaks/".$exp_name ."_peaks.bed\n";

    push @macs2_call_peak_cmd, $macs2_cmd;
    push @macs2_bdgcmp_cmd,    $bdgcmp_cmd;
    push @sort_cmd,            $sort_cmd;
    push @bigwig_cmd,          $bigwig_cmd;
    push @mv_peak_cmd,         $mv_peak_cmd;

}
close INFILE;


open(OUTFILE, "> macs2.output/macs2_cmd_list");
print OUTFILE join("",@macs2_call_peak_cmd);
close OUTFILE;


open(OUTFILE,"> macs2.output/bdgcmp_cmd_list");
print OUTFILE join("",@macs2_bdgcmp_cmd);
close OUTFILE;

open(OUTFILE,"> macs2.output/sort_cmd_list");
print OUTFILE join("",@sort_cmd);
close OUTFILE;

open(OUTFILE,"> macs2.output/bigwig_cmd_list");
print OUTFILE join("",@bigwig_cmd);
close OUTFILE;

open(OUTFILE,"> macs2.output/mv_peak_cmd_list");
print OUTFILE join("",@mv_peak_cmd);
close OUTFILE;


`parallel --no-notice -j 5 < macs2.output/macs2_cmd_list`;
`parallel --no-notice -j 5 < macs2.output/bdgcmp_cmd_list`;
`parallel --no-notice -j 5 < macs2.output/sort_cmd_list`;
`parallel --no-notice -j 5 < macs2.output/bigwig_cmd_list`;
`parallel --no-notice -j 5 < macs2.output/mv_peak_cmd_list`;





