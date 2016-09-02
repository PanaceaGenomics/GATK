#!/bin/bash
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR

module load samtools/0.1.19 # Required to contain bgzip & Tabix
module load GATK/3.6


################################################
# Include variables from file called "variables"
################################################

. variables


#########################
## Copy reference files
#########################

cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa.fai .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.dict .


#########################
# Fix issue with header
#########################

echo -e '@RG\tID:'${sampleID}'_lane1\tPL:ILLUMINA\tLB:Library\tSM:'${sampleID}'' > rg

rm ${sampleID}.GATK.recal*bai
samtools view -H ${sampleID}.GATK.recal.bam | grep -v '@RG' | cat - rg > header
samtools reheader header ${sampleID}.GATK.recal.bam > ${sampleID}_gready.bam && mv ${sampleID}_gready.bam ${sampleID}.GATK.recal.bam
samtools index ${sampleID}.GATK.recal.bam


#########################
# GATK variants calling
#########################

java -jar /local/software/GATK/3.6/source/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	--emitRefConfidence GVCF \
	-R hs38.fa \
	-I ${sampleID}.GATK.recal.bam \
	-o ${sampleID}.g.vcf


############################################################################
# Get alternative genotypes (include multiple input GVCFs here if needed)
############################################################################

java -jar /local/software/GATK/3.6/source/GenomeAnalysisTK.jar \
	-T GenotypeGVCFs \
	-R hs38.fa \
	-V ${sampleID}.g.vcf \
	-o ${sampleID}.vcf


#############################################################
## Qsub script 3 if final VCF file exists and is non-empty
#############################################################

if [ -s ${sampleID}.vcf ]; then
        qsub 3_ANNOTATE_V4.0.sh
fi


######################################
# BGzip and index GVCF for future use
######################################

bgzip ${sampleID}.g.vcf
tabix -p vcf ${sampleID}.g.vcf.gz


###########
# Clean up
###########

rm hs38.fa hs38.fa.fai hs38.dict rg header
