#! /usr/bin/perl

use strict;
my $exp_file      = $ARGV[0];



my $STAR_INDEX           = "/srv/persistent/keliu/genomes/mm10.STAR";
my $RT_GTF_ANNOTATION    = "/users/keliu/project/HeLin/annotation/mm10.repbase.gtf";
my $GENE_GTF_ANNOTATION  = "/srv/persistent/keliu/genomes/mm10/annotation/gencode.VM6.annotation.gtf";


############# make the necessary directory #########

my @dir_list =qw / tmp alignment fastq   RData/;
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
chdir "fastq";
my @fastq_dump_cmd_list;
grep{
    my $line=$_;
    $line=~s/\r//g;
    chomp $line;
    if($line ne "" ){
        push @fastq_dump_cmd_list,"fastq-dump --gzip --split-files $line \n" if  (($line!~/^#/) && (! -e "../alignment/${line}.bam"));
    }
} @sra_no;

open(OUTFILE,"> fastq_dump_cmd_list");
print OUTFILE join("",@fastq_dump_cmd_list);
close OUTFILE;

`parallel -j 5 --no-notice < fastq_dump_cmd_list`;
`rm fastq_dump_cmd_list`;
chdir "../";


#################### fastq read mapping ####################
my $star_option = "--runThreadN 8 --genomeLoad NoSharedMemory --limitBAMsortRAM 0   --outSAMtype BAM SortedByCoordinate  --outSAMstrandField intronMotif   --outSAMattributes NH   HI   NM   MD   AS   XS --outFilterMultimapNmax 20000   --outFilterMultimapScoreRange 1   --outFilterScoreMinOverLread 0.66   --outFilterMatchNminOverLread 0.33   --outFilterMismatchNmax 10   --alignIntronMax 500000   --alignMatesGapMax 1000000   --alignSJDBoverhangMin 1   --sjdbOverhang 100   --sjdbScore 2  --readFilesCommand zcat";

my $map_cmd;
my $filter_cmd;
my @map_cmd_list;
my @filter_cmd_list;

foreach my $line (@sra_no){
    chomp $line;
    next if $line eq "";
    next if $line=~/^#/;
    $line=~s/\r//;
    my $sra_no = $line;
    
    if ( -e "fastq/${sra_no}_2.fastq.gz" ){
        $map_cmd    = "STAR --genomeDir $STAR_INDEX $star_option --outFileNamePrefix  ./tmp/${sra_no}/  --readFilesIn fastq/${sra_no}_1.fastq.gz  fastq/${sra_no}_2.fastq.gz  \n";
        $filter_cmd = "samtools view -b -F 1804    ./tmp/${sra_no}/Aligned.sortedByCoord.out.bam        > alignment/${sra_no}.bam\n";
        
    }else{
        $map_cmd    = "STAR --genomeDir $STAR_INDEX $star_option --outFileNamePrefix  ./tmp/${sra_no}/  --readFilesIn fastq/${sra_no}_1.fastq.gz  \n";
        $filter_cmd = "samtools view -b -F 1804    ./tmp/${sra_no}/Aligned.sortedByCoord.out.bam        > alignment/${sra_no}.bam\n";
    }

    if (! -e "alignment/${sra_no}.bam"){
        push @map_cmd_list,$map_cmd;
        push @filter_cmd_list,$filter_cmd;
    }
}


open(OUTFILE,"> tmp/map_cmd_list");
print OUTFILE join("",@map_cmd_list);
close OUTFILE;



open(OUTFILE,"> tmp/filter_cmd_list");
print OUTFILE join("",@filter_cmd_list);
close OUTFILE;



`parallel --no-notice -j 8 < tmp/map_cmd_list`;
`parallel --no-notice -j 8 < tmp/filter_cmd_list`;


############### use featureCounts to do read counting ############
chdir 'alignment';
my $out_put       = `ls |grep bam`;
my $bam_file      = join("  ",split /\n+/,$out_put);
my $featureCounts_cmd;




$featureCounts_cmd = "featureCounts -a $RT_GTF_ANNOTATION   -B -C -p   --primary -T 10 -O -M -g repeat  -o ../RData/rt.read.count.matrix    $bam_file";
`$featureCounts_cmd`;

$featureCounts_cmd = "featureCounts -a $GENE_GTF_ANNOTATION -B -C -p   --primary -T 10 -O    -g gene_id -o ../RData/gene.read.count.matrix  $bam_file";
`$featureCounts_cmd`;
chdir '../';









