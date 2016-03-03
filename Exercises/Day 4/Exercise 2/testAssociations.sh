#!/bin/bash
#--------------------------------------
#run spaunge
#--------------------------------------
rm *.sam *.bam *.bai
../bin/spaunge task=simulate ind=100 coverage=25 explVar=0.05 SNPFreq=0.5 verbose seqLen=1000


#--------------------------------------
#make bam and index files
#--------------------------------------
rm simulated.bams.list
for file in spaunge_Ind*.sam
do
	samtools view -bSh ${file} > ${file}.bam
	samtools index ${file}.bam
	echo ${file}.bam >> simulated.bams.list
done

samtools faidx spaunge_refSeq.fasta 

#-------------------------------------
#Run ANGSD
#-------------------------------------
angsd -GL 2 -doGlf 2 -doMaf 1 -doMajorMinor 1 -SNP_pval 0.001 -bam simulated.bams.list -out spaunge_ANGSD

#--------------------------------------
#Run Spaunge
#-------------------------------------
../bin/spaunge zbeagle=spaunge_ANGSD.beagle.gz pheno=spaunge_phenotypes.txt verbose
