############################################
# Before running this script, please ensure:
# 1. Bowtie2, R, PICARD,featureCounts have been instatlleld
# 2. Edit the file code/chip.seq.preprocess.pl to  specify the correct path for mm10 bowtie2 index and other tools 
# 3. Edit the file code/rna.seq.preprocess.pl  to  specify the correct path for mm10 STAR    index
#
# After running this script, all .RData files needed to reproduce the client-side results could be downloaded under
# the two folders: chip.seq/RData  nad rna.seq/Rdata 
#########################

cd chip.seq

cd annotation
wget http://merlot.lbl.gov/keliu/HeLab-MartinProject.support.files/mm10.annotation/*.gtf
wait
cd ../

perl code/chip.seq.preprocess.pl  exp.file/mouse.chip.seq.csv 
wait

Rscript code/chip.seq.MA.R  exp.file/mouse.chip.seq.csv 
wait

Rscript code/chip.seq.MA.R  exp.file/otx2.csv 
wait

Rscript code/chip.seq.MA.R  exp.file/oct4.csv 
wait

Rscript code/chip.seq.MA.R  exp.file/Arid3a.csv 
wait



cd ../rna.seq

cd annotation
wget http://merlot.lbl.gov/keliu/HeLab-MartinProject.support.files/mm10.annotation/*.gtf
wait
cd ../

perl code/rna.seq.preprocess.pl  exp.file/Dux.dox.csv 
wait

perl code/rna.seq.preprocess.pl  exp.file/otx2.csv 
wait

Rscript code/rna.seq.MA.R exp.file/Dux.dox.csv 
wait

Rscript code/rna.seq.MA.R exp.file/otx2.csv 
wait

