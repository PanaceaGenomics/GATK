#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=60:00:00
#PBS -l mem=20gb
cd $PBS_O_WORKDIR

module load samtools/0.1.19
module load picard/1.97
module load bwa/0.7.12
module load GATK/3.6


#######################
#include variables from file called "variables"
#######################

. variables


#######################
#copying reference idx files
#######################

cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa.fai .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.dict .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa.amb .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa.ann .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa.bwt .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa.pac .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/hs38.fa.sa .
cp /scratch/EXOME_DATA/REF_GENOMES/HG38/BROAD/hg38bundle/dbsnp_144.hg38.vcf .

##########################
#Concatentate Multilane FQ Files if required
##########################

if [ $no_lanes -eq "1" ]; then
	fq_align_1=$lane1_pair1
	fq_align_2=$lane1_pair2
fi

if [ $no_lanes -eq "2"  ]; then
	cat $lane1_pair1 $lane2_pair1 > "$sampleID"_cat1.fq
	cat $lane1_pair2 $lane2_pair2 > "$sampleID"_cat2.fq
	fq_align_1="$sampleID"_cat1.fq
	fq_align_2="$sampleID"_cat2.fq
fi

if [ $no_lanes -eq "3" ]; then
	cat $lane1_pair1 $lane2_pair1 $lane3_pair1 > "$sampleID"_cat1.fq
	cat $lane1_pair2 $lane2_pair2 $lane3_pair3 > "$sampleID"_cat2.fq
	fq_align_1="$sampleID"_cat1.fq
	fq_align_2="$sampleID"_cat2.fq
fi


###########################################
# aligning fastq files  with the reference genome using bwa-mem
# -t = Threads -M = flag shorter split hits as secondary
# -R = Readgroups -O = Gap open penalty -E =Gap extension penalty
###########################################

bwa mem \
	-t 16 \
	-M \
	-R '@RG\tID:'$sampleID'_lane1\tSM:'$sampleID'\tPL:ILLUMINA\tLB:Library' \
	-O 65 \
	-E 7 \
	hs38.fa \
	$fq_align_1 $fq_align_2 \
	> "$sampleID"_aligned.sam

## Remove contactenated FastQ files if they were generated. 
rm "$sampleID"_cat1.fq "$sampleID"_cat2.fq


######################
# convert SAM to BAM
######################

samtools view -bS -o "$sampleID"_aligned.bam "$sampleID"_aligned.sam


######################
# Picard tools 
######################

# sort bam file
picard SortSam \
	INPUT="$sampleID"_aligned.bam \
	OUTPUT="$sampleID"_aligned_sorted.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=tmp1 \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=2000000

# mark duplicates
picard MarkDuplicates \
	INPUT="$sampleID"_aligned_sorted.bam \
	METRICS_FILE="$sampleID"_dup_metrics \
	OUTPUT="$sampleID"_marked_dups_sorted.bam \
	TMP_DIR=tmp1 \
	VALIDATION_STRINGENCY=SILENT

# Sort BAM file
picard SortSam \
	INPUT="$sampleID"_marked_dups_sorted.bam \
	OUTPUT="$sampleID".DelDup.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=tmp2 \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=2000000

# Index BAM file
picard BuildBamIndex \
	INPUT="$sampleID".DelDup.bam \
	TMP_DIR=tmp3 \
	VALIDATION_STRINGENCY=SILENT

# Fix mate pair information by picard
picard FixMateInformation \
	INPUT="$sampleID".DelDup.bam \
	OUTPUT="$sampleID".GATK.fixedmateinfo.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=tmp2 \
	VALIDATION_STRINGENCY=SILENT \
	MAX_RECORDS_IN_RAM=500000 \
	CREATE_INDEX=true


#############################################
#Remove SAM, Intermediate BAMs(Retaining original *_aligned.bam)
#############################################

rm "$sampleID"_aligned.sam "$sampleID"_aligned_sorted.bam "$sampleID"_marked_dups_sorted.bam "$sampleID".DelDup.bam "$sampleID".DelDup.bai


##########################################################################################
# GATK BQRS; NB: Indel realignment not required in using HaplotypeCaller downstream
##########################################################################################

# Recalibratiing base quality (longer: more than 60 minutes)
java -jar /local/software/GATK/3.6/source/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-I "$sampleID".GATK.fixedmateinfo.bam \
	-R  hs38.fa \
	-knownSites dbsnp_144.hg38.vcf \
	-o "$sampleID".recal_data.grp \
	-nct 16

# Shorter: less than 15 minutes
java -jar /local/software/GATK/3.6/source/GenomeAnalysisTK.jar \
	-T PrintReads \
	-R hs38.fa \
	-I "$sampleID".GATK.fixedmateinfo.bam \
	-BQSR "$sampleID".recal_data.grp \
	-l INFO -log "$sampleID".newQual.log \
	-o "$sampleID".GATK.recal.bam \
	-nct 16

# Index final BAM file
picard BuildBamIndex \
	INPUT="$sampleID".GATK.recal.bam \
	TMP_DIR=tmp1 \
	VALIDATION_STRINGENCY=SILENT


#############################################################
## Qsub script 2 if final BAM file exists and is non-empty
#############################################################

if [ -s "$sampleID".GATK.recal.bam ]; then
	qsub 2_CALLING_V4.0.sh
fi


##############################
# Remove reference files, Penultimate BAM. Picard tmp directories
##############################

rm hs38.fa.sa hs38.fa.amb hs38.fa.ann hs38.fa.pac hs38.fa.bwt hs38.fa.fai hs38.fa hs38.dict dbsnp_144.hg38.vcf dbsnp_144.hg38.vcf.idx "$sampleID".GATK.fixedmateinfo.bam "$sampleID".GATK.fixedmateinfo.bai
rm -rf tmp1 tmp2 tmp3
