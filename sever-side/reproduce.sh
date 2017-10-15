############################################
# Before running this script, please ensure:
# 1. Bowtie2, R, PICARD,featureCounts have been instatlleld
# 2. Edit the file code/chip.seq.preprocess.pl to  fill correct mm10 bowtie2 index
# 3. Edit the file code/rna.seq.preprocess.pl  to  fill correct mm10  STAR index
#
#########################

cd chip.seq

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

perl code/rna.seq.preprocess.pl  exp.file/Dux.dox.csv 
wait

perl code/rna.seq.preprocess.pl  exp.file/otx2.csv 
wait

Rscript code/rna.seq.MA.R exp.file/Dux.dox.csv 
wait

Rscript code/rna.seq.MA.R exp.file/otx2.csv 
wait

