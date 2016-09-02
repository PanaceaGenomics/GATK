#/bin/bash
#PBS -l walltime=2:00:00
cd $PBS_O_WORKDIR

module load bcftools/1.2.1
module load samtools/1.2.1
module load annovar/2015Dec14


####################################################################################################################################
# load variables from file "variables" need: "$sample_ID" "$diseaseDB" "$nondiseaseDB" "$somaticDB"
# soton databases should be in /scratch/GENE_DATABASES/humandb/HG38/ and named as hg38_bloodcancer.txt for $diseaseDB=bloodcancer
####################################################################################################################################

. variables


#####################################
# CONVERT VCF file                  #
#####################################

# filter depth>=4, convert vcf to annovar format
vcfutils.pl varFilter -d 4 -1 0 -2 0 -3 0 -4 0 "$sampleID".vcf > "$sampleID"_filtered.vcf
convert2annovar.pl -format vcf4 "$sampleID"_filtered.vcf -outfile "$sampleID"_annovar.input -allsample -withfreq -include 2>annovar.log

#Check number of samples and print to STDOUT
sample_no=`awk '{print NF-17}' "$sampleID"_annovar.input | head -1`
echo "number of samples in vcf: $sample_no"


#######################
#  ANNOVAR ANNOTATION; have removed in house databases #
#######################

table_annovar.pl "$sampleID"_annovar.input /scratch/GENE_DATABASES/humandb/HG38/ -buildver hg38 -out "$sampleID" -remove -protocol refGene,knownGene,avsnp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_eas,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,exac03,cosmic70,dbnsfp30a,dbnsfp31a_interpro,dbscsnv11,clinvar_20160302,hrcr1,kaviar_20150923,nci60 -operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . 2>>annovar.log


##################
# ADD Extra info #
##################

#GET GENOMIC SEQUENCE +/- 50bp OF VARIANT AND FLAG  - output to seqcontext
cp /scratch/GENE_DATABASES/humandb/HG38/EXTRA/getseq.c .
cc -o getseq getseq.c
./getseq "$sampleID"_annovar.input

#add VCF INFO fields as separate columns  (see end of this file for descriptions)
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"AC=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_AC"}' "$sampleID"_annovar.input > AC
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"AF=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_AF"}' AC > AF
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"AN=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_AN"}' AF > AN
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"BaseQRankSum=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_BaseQRankSum"}' AN > BaseQRankSum
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"ClippingRankSum=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_ClippingRankSum"}' BaseQRankSum > ClippingRankSum
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"DP=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_DP"}' ClippingRankSum > DP
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"DS=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_DS"}' DP > DS
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"END=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_END"}' DS > END
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"ExcessHet=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_ExcessHet"}' END > ExcessHet
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"FS=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_FS"}' ExcessHet > FS
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"HaplotypeScore=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_HaplotypeScore"}' FS > HaplotypeScore
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"InbreedingCoeff=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_InbreedingCoeff"}' HaplotypeScore > InbreedingCoeff
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"MLEAC=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_MLEAC"}' InbreedingCoeff > MLEAC
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"MLEAF=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_MLEAF"}' MLEAC > MLEAF
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"MQ=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_MQ"}' MLEAF > MQ
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"MQRankSum=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_MQRankSum"}' MQ > MQRankSum
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"QD=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_QD"}' MQRankSum > QD
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"RAW_MQ=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_RAW_MQ"}' QD > RAW_MQ
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"ReadPosRankSum=.*;",a)) {split(a[0],b,"="); split(b[2],c,";"); print $0,c[1];} else print $0,"no_ReadPosRankSum"}' RAW_MQ > ReadPosRankSum
awk 'BEGIN {FS=OFS="\t"} {if(match($16,"SOR=.*",a)) {split(a[0],b,"="); print $0,b[2];} else print "no_SOR",$0}' ReadPosRankSum > SOR


############
# COLLATE  #
############

#add 6 column keys
awk 'BEGIN{FS=OFS="\t"}{print $1":"$2":"$3":"$4":"$5,$0}' "$sampleID"_annovar.input >input.key
awk 'BEGIN{FS=OFS="\t"}{print $1":"$2":"$3":"$4":"$5,$0}' "$sampleID".hg38_multianno.txt | grep -v "Chr:Start:End:Ref:Alt" >multianno.key
awk 'BEGIN{FS=OFS="\t"}{print $1":"$2":"$3":"$4":"$5,$0}' SOR >SOR_file
cut -f 1,20- SOR_file > SOR_short
awk 'BEGIN{FS=OFS="\t"}{print $1":"$2":"$3":"$4":"$5,$6,$7}' seqcontext >seqcontext_file

#append matched annotation
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0;next}{print $0,a[$1]?a[$1]:"not annotated"}' multianno.key input.key >t1
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0;next}{print $0,a[$1]?a[$1]:"no_INFO"}' SOR_short t1 >t2
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2;b[$1]=$3;next}{print $0,a[$1]?a[$1]"\t"b[$1]:".\t."}' seqcontext_file t2 > t3  

#Add Grantham based on refGene
awk -v smp_no="$sample_no" 'BEGIN{FS=OFS="\t"} {if($(28+smp_no)=="nonsynonymous SNV") {split($(29+smp_no),a,","); split(a[1],b,":"); print $0,substr(b[5],3,1)substr(b[5],length(b[5]),1);} else print $0,"XX"}' t3 > t4
awk -v smp_no="$sample_no" 'BEGIN{FS=OFS="\t"} NR==FNR{a[$5]=$4;b[$5]=$3;next}{print $0,a[$(132+smp_no)]?a[$(132+smp_no)]"_"b[$(132+smp_no)]:"."}' /scratch/GENE_DATABASES/humandb/HG38/EXTRA/GS2 t4 > t5

#format gene name for cross referencing
awk -v smp_no="$sample_no" 'BEGIN{FS=OFS="\t"}{if($(25+smp_no) != "intergenic") {split($(26+smp_no),a,";"); split(a[1],b,"("); split(b[1],c,","); print c[1],$0} else print $(26+smp_no),$0}' t5 > t6

# add mutability flag >= 0.001697 bp per novel potentially pathogenic variant mean + 2*sd
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$1;next}{print $0,a[$1]?"1":"0"}' /scratch/GENE_DATABASES/humandb/HG38/EXTRA/highly_mutable t6 > t7
# add false positive Genes list from Fuentes et al 2011
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$1;next}{print $0,a[$1]?"1":"0"}' /scratch/GENE_DATABASES/humandb/HG38/EXTRA/fuentes_falsePos t7 > t8
# add ACMG gene list flag
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$1;next}{print $0,a[$1]?"1":"0"}' /scratch/GENE_DATABASES/humandb/HG38/EXTRA/ACMG t8 > t9

#ADD headers
vcf_header=`grep "#CHROM" "$sampleID"_filtered.vcf`
multianno_header=`head -1 "$sampleID".hg38_multianno.txt`
annovar_header="annovar_chr\tannovar_bp1\tannovar_bp2\tannovar_ref\tannovar_alt\tGroupAlleleFreq\tGenotypeQuality?last_sample_or_average\tGenotypeDepth?last_sample_or_average"
extra_header="KEY3\tINFO_AC\tINFO_AF\tINFO_AN\tINFO_BaseQRankSum\tINFO_ClippingRankSum\tINFO_DP\tINFO_DS\tINFO_END\tINFO_ExcessHet\tINFO_FS\tINFO_HaplotypeScore\tINFO_InbreedingCoeff\tINFO_MLEAC\tINFO_MLEAF\tINFO_MQ\tINFO_MQRankSum\tINFO_QD\tINFO_RAW_MQ\tINFO_ReadPosRankSum\tINFO_SOR\tseq\tseq_flag\tAAchange\tGrantham\tMutability\tFuentesFalsePositive\tACMG"

printf "GENE\tKEY1\t$annovar_header\t$vcf_header\tKEY2\t$multianno_header\t$extra_header\n" > t10
cat t9 >> t10

#remove repeated columns
cut -f 11-20,22-110,112- t10 > ${sampleID}_annotated.var 

#extract exonic|splicing
awk -v smp_no="$sample_no" 'BEGIN{OFS=FS="\t"}{if(($(15+smp_no) ~ "exonic") || ($(15+smp_no) ~ "splicing") || ($1 == "#CHROM")) print $0}' ${sampleID}_annotated.var > ${sampleID}_annotated_coding.var


########################
# TEMP FILE CLEAN UP   #
########################

rm "$sampleID"_annovar.input t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 seqcontext seqcontext_file multianno.key input.key getseq ${sampleID}.hg38_multianno.txt AC AF AN BaseQRankSum ClippingRankSum DP DS END ExcessHet FS HaplotypeScore InbreedingCoeff MLEAC MLEAF MQ MQRankSum QD RAW_MQ ReadPosRankSum SOR SOR_file SOR_short getseq.c getseq


#############################################################
# Give appropriate file permissions to HGIG group
#############################################################

chgrp -c hgig . ; chgrp -Rc hgig * ; chmod -c g+rwX . ; chmod -Rc g+rwX *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##VCF INFO fields
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinb erg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in t he same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
