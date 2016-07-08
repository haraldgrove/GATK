#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##Updated APR 2016, Harald Grove
##Replaced software with newest version
##Replaced GATK bundle-files with version hg38
##Added region information from exome capture kit
##-------------
##This script performs the entire variant-calling process upon one sample, following the Genome Analysis Toolkit (GATK)'s pipeline.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Run the command 'bash /path/to/GATK_individual.sh [options] sample_name'
##-------------


## Make sure the script stop at first error
set -e


##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
bwa_dir=/home/harald/software/bwa-0.7.13
picard_dir=/home/harald/software/picard-tools-2.1.1
gatk_dir=/home/harald/software/GATK-3.5
fastq_dir=/home/harald/Rawdata/ExomeSeq
ref_dir=/home/harald/Rawdata/hg38bundle
out_dir=/home/harald/Projects/GATKtest/wes
java8_dir=/usr/lib/jvm/jdk1.8.0_73/bin
java7_dir=/usr/lib/jvm/jdk1.7.0_79/bin

##-------------
##Step0-2: References
##-------------
ref_genome=${ref_dir}/Homo_sapiens_assembly38.fasta
indel_1=${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
indel_2=${ref_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz
DBSNP=${ref_dir}/dbsnp_144.hg38.vcf.gz
exon_bed=~/Rawdata/ExomeSeq/kit-info/remapped_agilentV5_regions_onlychr.bed

##-------------
##Step0-3: Other Parametres
##-------------
java_mem=4g
cores=4
gatk_num_threads=1
gatk_num_cpu_threads=4

##-------------
##Step0-4: Input Arguments
##-------------
while test $# -gt 0 ; do
        case "$1" in
                -h|--help)
						echo ""
                        echo "Usage: bash $0 [options] sample_name"
                        echo ""
                        echo "This script performs the entire variant-calling process upon one sample, following the Genome Analysis Toolkit (GATK)'s pipeline."
                        echo ""
                        echo "Options:"
                        echo "-h, --help				display this help and exit"
						echo "-v, --version				display version of this script and exit"
						echo "-XS, --no-summary			suppress the command summary before execution"
						echo "-XP, --no-prompt			suppress the user prompt before execution, only when the command summary is displayed"
						echo "-XX, --no-exec				suppress automatic execution, generating only script files"
                        echo "-e, --exome				call only exonic variants, drastically accelerating the Call Haplotype and the Genotype processes"
						echo ""
                        exit 0
                        ;;
				-v|--version)
						echo ""
						echo "GATK_individual.sh"
                        echo ""
						echo "Created MAR 2016"
						echo "Updated JUL 2016"
						echo "by"
						echo "PURIN WANGKIRATIKANT [purin.wan@mahidol.ac.th]"
						echo "Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)"
						echo "Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand"
						echo ""
						exit 0
						;;
				-XS|--no-summary)
						no_summary=1
						shift
						;;
				-XP|--no-prompt)
						no_prompt=1
						shift
						;;
				-XX|--no-exec)
						no_exec=1
						shift
						;;
				-XG|--no-genotyping)
						no_geno=1
						shift
						;;
                -e|--exome)
						seq_type='EXOME'
						bed_argument='-L '${exon_bed}
                        shift
                        ;;
				*)
						sample_name=$1
						shift
						;;
		esac
done

##-------------
##Step0-5: Default Value Setting
##-------------
if [[ ! -v seq_type ]] ; then
		seq_type='GENOME'
fi

##-------------
##Step0-6: Input Verification
##-------------
if [[ ! -e ${fastq_dir}/${sample_name}_1.fastq.gz ]] ; then
		echo
		echo 'Invalid SAMPLE NAME: '${sample_name}
		echo ${fastq_dir}/${sample_name}_1.fastq.gz not found.
		echo 'Terminated.'
		echo
		exit 1
fi
if [[ ! -e ${fastq_dir}/${sample_name}_2.fastq.gz ]] ; then
		echo
		echo 'Invalid SAMPLE NAME: '${sample_name}
		echo ${fastq_dir}/${sample_name}_2.fastq.gz not found.
		echo 'Terminated.'
		echo
		exit 1
fi

##-------------
##Step0-7: Summarisation & User's Confirmation Prompt
##-------------
if [[ ${no_summary} != 1 ]] ; then
		echo
		echo '---------------------------------------'
		echo 'INDIVIDUAL VARIANT CALLING PROCESS'
		echo 'SAMPLE NAME =			'${sample_name}
		echo 'SEQUENCED DATA =		'${seq_type}
		echo '---------------------------------------'
		echo

		if [[ ${no_prompt} != 1 ]] ; then
				while true ; do
						read -p "Are all the input arguments correct? (Y/N): " confirm
						case ${confirm} in
								Y|y)
										echo "Confirmed. Initiating..."
										echo
										break
										;;
								N|n)
										echo "Terminated."
										echo
										exit 1
										;;
								* )
										echo "Please enter Y or N."
										echo
										;;
						esac
				done
		fi
fi

##-------------
##Step0-8: Output Folders Creation
##-------------
mkdir -p ${out_dir}
mkdir -p ${out_dir}/${sample_name} ; mkdir -p ${out_dir}/${sample_name}/{Scripts,LOG,TEMP,SAM,BAM,BQSR,GVCF,VCF,QC,QC/FILTERED,Report}



##-------------
##Step1: Align
##-------------
cat << EOL > ${out_dir}/${sample_name}/Scripts/1_${sample_name}_align.sh
#!/bin/bash
set -e
##-------------
##Step1: Align
##-------------
${bwa_dir}/bwa mem -t ${cores} -R "@RG\tID:DM_${sample_name}\tSM:${sample_name}\tPL:Illumina\tLB:WES\tPU:unit1" ${ref_genome} \
${fastq_dir}/${sample_name}_1.fastq.gz ${fastq_dir}/${sample_name}_2.fastq.gz > ${out_dir}/${sample_name}/SAM/${sample_name}_aligned.sam

EOL



##-------------
##Step2: Sort
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/2_${sample_name}_sort.sh
#!/bin/bash
set -e
##-------------
##Step2: Sort
##-------------
${java8_dir}/java -Xmx${java_mem} -jar ${picard_dir}/picard.jar SortSam \
INPUT=${out_dir}/${sample_name}/SAM/${sample_name}_aligned.sam \
OUTPUT=${out_dir}/${sample_name}/BAM/${sample_name}_sorted.bam \
SORT_ORDER=coordinate \
TMP_DIR=${out_dir}/${sample_name}/TEMP

EOL



##-------------
##Step3: Deduplicate
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/3_${sample_name}_deduplicate.sh
#!/bin/bash
set -e
##-------------
##Step3: Deduplicate
##-------------
${java8_dir}/java -Xmx${java_mem} -jar ${picard_dir}/picard.jar MarkDuplicates \
INPUT=${out_dir}/${sample_name}/BAM/${sample_name}_sorted.bam \
OUTPUT=${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam \
METRICS_FILE=${out_dir}/${sample_name}/Report/${sample_name}_deduplication_metrics.txt \
CREATE_INDEX=TRUE TMP_DIR=${out_dir}/${sample_name}/TEMP

EOL



##-------------
##Step4: Build Index
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/4_${sample_name}_build_index.sh
#!/bin/bash
ser -e
##-------------
##Step4: Build Index
##-------------
${java8_dir}/java -Xmx${java_mem} -jar ${picard_dir}/picard.jar BuildBamIndex \
INPUT=${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam \
TMP_DIR=${out_dir}/${sample_name}/TEMP

EOL



##-------------
##Step5: Indel Realignment
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/5_${sample_name}_realign_indels.sh
#!/bin/bash
set -e
##-------------
##Step5-1: Create Aligner Target
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
--disable_auto_index_creation_and_locking_when_reading_rods \
-known ${indel_1} \
-known ${indel_2} \
-R ${ref_genome} \
-I ${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam \
-dt NONE \
-nt ${gatk_num_threads} \
-o ${out_dir}/${sample_name}/BAM/${sample_name}_indel_target_intervals.list \
-log ${out_dir}/${sample_name}/LOG/5-1_${sample_name}_indel_target_intervals.log

##-------------
##Step5-2: Realign Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T IndelRealigner \
--disable_auto_index_creation_and_locking_when_reading_rods \
-known ${indel_1} \
-known ${indel_2} \
-I ${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam \
-R ${ref_genome} \
-targetIntervals ${out_dir}/${sample_name}/BAM/${sample_name}_indel_target_intervals.list \
-dt NONE \
-o ${out_dir}/${sample_name}/BAM/${sample_name}_realigned.bam \
-log ${out_dir}/${sample_name}/LOG/5-2_${sample_name}_indel_realigned.log

EOL



##-------------
##Step6: Base Quality Score Recalibration
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/6_${sample_name}_recalibrate_base.sh
#!/bin/bash
set -e
##-------------
##Step6-1: Perform Base Recalibration
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R ${ref_genome} \
-knownSites ${indel_1} \
-knownSites ${indel_2} \
-knownSites ${DBSNP} \
-I ${out_dir}/${sample_name}/BAM/${sample_name}_realigned.bam \
-nct ${gatk_num_cpu_threads} \
-o ${out_dir}/${sample_name}/BQSR/${sample_name}_perform_bqsr.table \
-log ${out_dir}/${sample_name}/LOG/6-1_${sample_name}_perform_bqsr.log

##-------------
##Step6-2: Generate Post-BQSR Table
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R ${ref_genome} \
-knownSites ${indel_1} \
-knownSites ${indel_2} \
-knownSites ${DBSNP} \
-I ${out_dir}/${sample_name}/BAM/${sample_name}_realigned.bam \
-nct ${gatk_num_cpu_threads} \
-BQSR ${out_dir}/${sample_name}/BQSR/${sample_name}_perform_bqsr.table \
-o ${out_dir}/${sample_name}/BQSR/${sample_name}_after_bqsr.table \
-log ${out_dir}/${sample_name}/LOG/6-2_${sample_name}_after_bqsr.log

##-------------
##Step6-3: Plot Base Recalibration
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R ${ref_genome} \
-before ${out_dir}/${sample_name}/BQSR/${sample_name}_perform_bqsr.table \
-after ${out_dir}/${sample_name}/BQSR/${sample_name}_after_bqsr.table \
-plots ${out_dir}/${sample_name}/Report/${sample_name}_bqsr.pdf \
-log ${out_dir}/${sample_name}/LOG/6-3_${sample_name}_plot_bqsr.log ;

##-------------
##Step6-4: Print Reads
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T PrintReads \
-R ${ref_genome} \
--disable_auto_index_creation_and_locking_when_reading_rods \
-I ${out_dir}/${sample_name}/BAM/${sample_name}_realigned.bam \
-BQSR ${out_dir}/${sample_name}/BQSR/${sample_name}_perform_bqsr.table \
-dt NONE \
-EOQ \
-nct ${gatk_num_cpu_threads} \
-o ${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
-log ${out_dir}/${sample_name}/LOG/6-4_${sample_name}_final_bam.log

EOL



##-------------
##Step7: Call Haplotype
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/7_${sample_name}_call_haplotype.sh
#!/bin/bash
set -e
##-------------
##Step7: Call Haplotype
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ${ref_genome} \
--input_file ${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
--genotyping_mode DISCOVERY # Default\
-stand_emit_conf 30 # Default \
-stand_call_conf 30 # Default \
${bed_argument} \
-A DepthPerSampleHC \
-A ClippingRankSumTest \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A FisherStrand \
-A GCContent \
-A AlleleBalanceBySample \
-A AlleleBalance \
-A QualByDepth \
-pairHMM VECTOR_LOGLESS_CACHING \
<<<<<<< HEAD
-nct ${gatk_num_cpu_threads} \
-o ${out_dir}/${sample_name}/GVCF/${sample_name}_GATK.gvcf \
-log ${out_dir}/${sample_name}/LOG/07_${sample_name}_haplotype_caller.log
=======
-o ${out_dir}/${sample_name}/GVCF/${sample_name}_GATK.g.vcf \
-log ${out_dir}/${sample_name}/LOG/7_${sample_name}_haplotype_caller.log
>>>>>>> 4f948e132910f1a7000b00340852ba4ef531a29b

EOL



##-------------
##Step8: Genotype
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/8_${sample_name}_genotype_gvcf.sh
#!/bin/bash
set -e
##-------------
##Step8: Genotype
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/GVCF/${sample_name}_GATK.g.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-nt 1 \
${bed_argument} \
-o ${out_dir}/${sample_name}/VCF/${sample_name}_RAW.vcf \
-log ${out_dir}/${sample_name}/LOG/8_${sample_name}_genotype_gvcf.log

EOL



##-------------
##Step9: SNV Quality Control
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/9_${sample_name}_SNV_quality_control.sh
#!/bin/bash
set -e
##-------------
##Step9-1-1: Extract SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VCF/${sample_name}_RAW.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType SNP \
--excludeFiltered \
-nt 1 \
-o  ${out_dir}/${sample_name}/QC/${sample_name}_RAW_SNV.vcf \
-log ${out_dir}/${sample_name}/LOG/9-1-1_${sample_name}_QC_select_snv.log

##-------------
##Step9-1-2: Extract Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VCF/${sample_name}_RAW.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType INDEL \
-selectType MNP \
-selectType MIXED \
-selectType SYMBOLIC \
--excludeFiltered \
-nt 1 \
-o ${out_dir}/${sample_name}/QC/${sample_name}_RAW_INDEL.vcf \
-log ${out_dir}/${sample_name}/LOG/9-1-2_${sample_name}_QC_select_INDEL.log

##-------------
##Step9-2-1: Annotate SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/QC/${sample_name}_RAW_SNV.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--dbsnp ${DBSNP} \
-L ${out_dir}/${sample_name}/QC/${sample_name}_RAW_SNV.vcf \
-A GCContent \
-A VariantType \
-dt NONE \
-nt 1 \
-o ${out_dir}/${sample_name}/QC/${sample_name}_RAW_SNV_ANNOTATED.vcf \
-log ${out_dir}/${sample_name}/LOG/9-2-1_${sample_name}_QC_snv_annotation.log

##-------------
##Step9-3-1: Filter SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/QC/${sample_name}_RAW_SNV_ANNOTATED.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--filterExpression 'QD < 2.0' \
--filterName 'QD' \
--filterExpression 'MQ < 30.0' \
--filterName 'MQ' \
--filterExpression 'FS > 40.0' \
--filterName 'FS' \
--filterExpression 'MQRankSum < -12.5' \
--filterName 'MQRankSum' \
--filterExpression 'ReadPosRankSum < -8.0' \
--filterName 'ReadPosRankSum' \
--filterExpression 'DP < 8.0' \
--filterName 'DP' \
--logging_level ERROR \
-o ${out_dir}/${sample_name}/QC/${sample_name}_FILTERED_SNV.vcf \
-log ${out_dir}/${sample_name}/LOG/9-3-1_${sample_name}_QC_filter_snv.log

##-------------
##Step9-4-1: Clean SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/QC/${sample_name}_FILTERED_SNV.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--excludeFiltered \
-nt 1 \
-o  ${out_dir}/${sample_name}/QC/FILTERED/${sample_name}_CLEAN_SNV.vcf \
-log ${out_dir}/${sample_name}/LOG/9-4-1_${sample_name}_QC_clean_snv.log

##-------------
##Step9-5: Combine SNVs + Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T CombineVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/QC/FILTERED/${sample_name}_CLEAN_SNV.vcf \
--variant ${out_dir}/${sample_name}/QC/${sample_name}_RAW_INDEL.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--genotypemergeoption UNSORTED \
-o ${out_dir}/${sample_name}/QC/FILTERED/${sample_name}_CLEAN_SNV+INDEL.vcf \
-log ${out_dir}/${sample_name}/LOG/9-5_${sample_name}_QC_combine_variants.log

EOL





##-------------
##MASTER SCRIPT
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Scripts/${sample_name}_GATK.sh
#!/bin/bash
set -e
##-------------
##${sample_name}'s Vaiant Calling
##-------------
(bash ${out_dir}/${sample_name}/Scripts/1_${sample_name}_align.sh) 2>&1 | tee ${out_dir}/${sample_name}/LOG/1_${sample_name}_alignment.log
(bash ${out_dir}/${sample_name}/Scripts/2_${sample_name}_sort.sh) 2>&1 | tee ${out_dir}/${sample_name}/LOG/2_${sample_name}_sort.log
(bash ${out_dir}/${sample_name}/Scripts/3_${sample_name}_deduplicate.sh) 2>&1 | tee ${out_dir}/${sample_name}/LOG/3_${sample_name}_deduplication.log
(bash ${out_dir}/${sample_name}/Scripts/4_${sample_name}_build_index.sh) 2>&1 | tee ${out_dir}/${sample_name}/LOG/4_${sample_name}_building_index.log

if [[ -e ${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam ]] ; then
		rm ${out_dir}/${sample_name}/SAM/${sample_name}_aligned.sam
fi

bash ${out_dir}/${sample_name}/Scripts/5_${sample_name}_realign_indels.sh
bash ${out_dir}/${sample_name}/Scripts/6_${sample_name}_recalibrate_base.sh

if [[ -e ${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam ]]; then
  rm -f ${out_dir}/${sample_name}/BAM/${sample_name}_{deduplicated,sorted,realigned}.{bam,bai}
fi

bash ${out_dir}/${sample_name}/Scripts/7_${sample_name}_call_haplotype.sh

if [[ ${no_geno} != 1 ]] ; then
		bash ${out_dir}/${sample_name}/Scripts/8_${sample_name}_genotype.sh
		bash ${out_dir}/${sample_name}/Scripts/9_${sample_name}_SNV_quality_control.sh
fi

EOL





##-------------
##EXECUTION
##-------------
if [[ ${no_exec} != 1 ]] ; then
		bash ${out_dir}/${sample_name}/Scripts/${sample_name}_GATK.sh
fi
