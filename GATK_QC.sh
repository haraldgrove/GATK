#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script performs the post-calling process quality reports upon one sample.
##BQSR Plot and Deduplication Metrics have been integrated into the main variant-calling flow.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Check if the individual GATK process for every sample has been accomplished.
##iii>	Run the command 'bash /path/to/GATK_QC.sh sample_name'
##-------------





##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
picard_dir=/home/harald/software/picard-tools-2.1.1
gatk_dir=/home/harald/software/GATK-3.5
samtools_dir=/home/harald/software/samtools # 1.3-20-gd49c73b
java8_dir=/usr/lib/jvm/jdk1.8.0_73/bin
home_dir=/home/harald/Projects/GATKtest/WES
ref_dir=/home/harald/Rawdata/hg38bundle
out_dir=${home_dir}/Output
bed_dir=${home_dir}/bed
target_bed=${bed_dir}/SureSelect5.target.bed # xx
bait_bed=${bed_dir}/SureSelect5.bait.bed # xx
exon_bed=${bed_dir}/exon_hg19.bed # xx
gene_list=${home_dir}/gene_list.txt # xx

##-------------
##Step0-2: References
##-------------
ref_genome=${ref_dir}/Homo_sapiens_assembly38.fasta

##-------------
##Step0-3: Other Parametres
##-------------
java_mem=4g

##-------------
##Step0-4: Sample Name
##-------------
sample_name=$1

##-------------
##Step0-5: Output Folders Creation
##-------------
mkdir -p ${out_dir}/${sample_name}/{Report,Script/QC}

##-------------
##Step0-6: Picard Bed Files
##-------------
(${samtools_dir}/samtools view -H ${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam | grep "@SQ" ; sed 's/\r//g' ${target_bed} | awk '{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g' ) >| ${bed_dir}/OnTarget.picard.bed
(${samtools_dir}/samtools view -H ${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam | grep "@SQ" ; sed 's/\r//g' ${bait_bed} | awk '{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g' ) >| ${bed_dir}/OnBait.picard.bed



##-------------
##Step1: Ti/Tv 
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/01_${sample_name}_TiTv.sh
#!/bin/bash
##-------------
##Step1-1: Prepare SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--excludeFiltered \
--variant ${out_dir}/${sample_name}/QC/FILTERED_ON_BAIT/${sample_name}_CLEAN_SNV+INDEL.vcf \
-o ${out_dir}/${sample_name}/TEMP/${sample_name}_QC_FILTERED_TiTv_SNV.vcf

##-------------
##Step1-2-1: Select Known Variants
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/QC/FILTERED_ON_BAIT/${sample_name}_CLEAN_SNV+INDEL.vcf \
--excludeFiltered \
--concordance ${ref_dir}/dbsnp_138.hg19.excluding_sites_after_129.vcf \
-o ${out_dir}/${sample_name}/TEMP/${sample_name}_QC_FILTERED_TiTv_SNV_KNOWN.vcf

##-------------
##Step1-2-2: Select Novel Variants
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/QC/FILTERED_ON_BAIT/${sample_name}_CLEAN_SNV+INDEL.vcf \
--excludeFiltered \
--discordance ${ref_dir}/dbsnp_138.hg19.excluding_sites_after_129.vcf \
-o ${out_dir}/${sample_name}/TEMP/${sample_name}_QC_FILTERED_TiTv_SNV_NOVEL.vcf

##-------------
##Step1-3-1: Calculate TiTv of Known Variants
##-------------
${samtools_dir}/bcftools/vcfutils.pl qstats ${out_dir}/${sample_name}/TEMP/${sample_name}_QC_FILTERED_TiTv_SNV_KNOWN.vcf \
> ${out_dir}/${sample_name}/Report/${sample_name}_TiTv_Known.txt

##-------------
##Step1-3-2: Calculate TiTv of Novel Variants
##-------------
${samtools_dir}/bcftools/vcfutils.pl qstats ${out_dir}/${sample_name}/TEMP/${sample_name}_QC_FILTERED_TiTv_SNV_NOVEL.vcf \
> ${out_dir}/${sample_name}/Report/${sample_name}_TiTv_Novel.txt

EOL



##-------------
##Step2: Depth of Coverage
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/02_${sample_name}_depth_of_coverage.sh
#!/bin/bash
##-------------
##Step2-1: Depth of Coverage on Exon
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R ${ref_genome} \
-I ${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
-geneList:REFSEQ ${gene_list} \
-L ${exon_bed} \
-mmq 20 \
-mbq 10 \
--outputFormat csv \
-omitBaseOutput \
-o ${out_dir}/${sample_name}/Report/${sample_name}_coverage.exon \
-ct 8 \
-ct 15 \
-ct 20

##-------------
##Step2-2: Depth of Coverage on Bait
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R ${ref_genome} \
-I ${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
-geneList:REFSEQ ${gene_list} \
-L ${bait_bed} \
-mmq 20 \
-mbq 10 \
--outputFormat csv \
-omitBaseOutput \
-o ${out_dir}/${sample_name}/Report/${sample_name}_coverage.bait \
-ct 8 \
-ct 15 \
-ct 20

##-------------
##Step2-3: Adding Suffixes
##-------------
rename 's/$/.csv/' ${out_dir}/${sample_name}/Report/${sample_name}_coverage*

EOL



##-------------
##Step3: Alignment Summary Metrics
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/03_${sample_name}_alignment_summary_metrics.sh
#!/bin/bash
##-------------
##Step3: Alignment Summary Metrics
##-------------
java -jar ${picard_dir}/CollectAlignmentSummaryMetrics.jar \
INPUT=${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
OUTPUT=${out_dir}/${sample_name}/Report/${sample_name}_alignment_summary_metrics.txt \
R=${ref_genome} \
VALIDATION_STRINGENCY=SILENT

EOL



##-------------
##Step4: Quality Score Distribution
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/04_${sample_name}_quality_score_distribution.sh
#!/bin/bash
##-------------
##Step4: Quality Score Distribution
##-------------
java -jar ${picard_dir}/QualityScoreDistribution.jar \
INPUT=${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
OUTPUT=${out_dir}/${sample_name}/Report/${sample_name}_quality_score_distribution.txt \
CHART=${out_dir}/${sample_name}/Report/${sample_name}_quality_score_distribution_chart.pdf \
R=${ref_genome} \
VALIDATION_STRINGENCY=SILENT \
INCLUDE_NO_CALLS=true

EOL



##-------------
##Step5: GC Bias Metrics
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/05_${sample_name}_gc_bias_metrics.sh
#!/bin/bash
##-------------
##Step5: GC Bias Metrics
##-------------
java -jar ${picard_dir}/CollectGcBiasMetrics.jar \
INPUT=${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
OUTPUT=${out_dir}/${sample_name}/Report/${sample_name}_gc_bias_metrics.txt \
CHART_OUTPUT=${out_dir}/${sample_name}/Report/${sample_name}_gc_bias_metrics.pdf \
SUMMARY_OUTPUT=${out_dir}/${sample_name}/Report/${sample_name}_gc_bias_summary.txt \
R=${ref_genome} \
VALIDATION_STRINGENCY=SILENT

EOL



##-------------
##Step6: Insert Size Metrics
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/06_${sample_name}_insert_size_metrics.sh
#!/bin/bash
##-------------
##Step6: Insert Size Metrics
##-------------
java -jar ${picard_dir}/CollectInsertSizeMetrics.jar \
INPUT=${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
OUTPUT=${out_dir}/${sample_name}/Report/${sample_name}_insert_size_metrics.txt \
H=${out_dir}/${sample_name}/Report/${sample_name}_insert_size_metrics_histogram.pdf \
R=${ref_genome} \
VALIDATION_STRINGENCY=SILENT

EOL



##-------------
##Step7: Mean Quality by Cycle
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/07_${sample_name}_mean_quality_by_cycle.sh
#!/bin/bash
##-------------
##Step7: Mean Quality by Cycle
##-------------
java -jar ${picard_dir}/MeanQualityByCycle.jar \
INPUT=${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
OUTPUT=${out_dir}/${sample_name}/Report/${sample_name}_mean_quality_by_cycle.txt \
CHART=${out_dir}/${sample_name}/Report/${sample_name}_mean_quality_by_cycle_chart.pdf \
R=${ref_genome} \
VALIDATION_STRINGENCY=SILENT

EOL



##-------------
##Step8: Hybridization Selection Metrics
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/08_${sample_name}_hybridization_selection_metrics.sh
#!/bin/bash
##-------------
##Step8: Hybridization Selection Metrics
##-------------
java -jar ${picard_dir}/CalculateHsMetrics.jar \
INPUT=${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
OUTPUT=${out_dir}/${sample_name}/Report/${sample_name}_hybridization_selection_metrics.txt \
PER_TARGET_COVERAGE=${out_dir}/${sample_name}/Report/${sample_name}_per_target_coverage.txt \
R=${ref_genome} \
BI=${bed_dir}/OnBait.picard.bed \
TI=${bed_dir}/OnTarget.picard.bed \
VALIDATION_STRINGENCY=SILENT

EOL





##-------------
##MASTER SCRIPT
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/QC/${sample_name}_QC.sh
#!/bin/bash
##-------------
##${sample_name}'s Post-Calling Quality Control
##-------------
bash ${out_dir}/${sample_name}/Script/QC/01_${sample_name}_TiTv.sh
bash ${out_dir}/${sample_name}/Script/QC/02_${sample_name}_depth_of_coverage.sh
bash ${out_dir}/${sample_name}/Script/QC/03_${sample_name}_alignment_summary_metrics.sh
bash ${out_dir}/${sample_name}/Script/QC/04_${sample_name}_quality_score_distribution.sh
bash ${out_dir}/${sample_name}/Script/QC/05_${sample_name}_gc_bias_metrics.sh
bash ${out_dir}/${sample_name}/Script/QC/06_${sample_name}_insert_size_metrics.sh
bash ${out_dir}/${sample_name}/Script/QC/07_${sample_name}_mean_quality_by_cycle.sh
bash ${out_dir}/${sample_name}/Script/QC/08_${sample_name}_hybridization_selection_metrics.sh

EOL





##-------------
##EXECUTION
##-------------
bash ${out_dir}/${sample_name}/Script/QC/${sample_name}_QC.sh
