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
version=0.1

##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
map_dir=/
fastq_dir=tiger/harald/SCD/fastq
ref_dir=tiger/harald/resources/hg38bundle
out_dir=tiger/harald/SCD

##-------------
##Step0-2: References
##-------------
ref_genome=${ref_dir}/Homo_sapiens_assembly38.fasta
indel_1=${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
indel_2=${ref_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz
DBSNP=${ref_dir}/dbsnp_144.hg38.vcf.gz
exon_bed=${ref_dir}/remapped_agilentV5_regions_onlychr.bed

##-------------
##Step0-3: Other Parametres
##-------------
cores=16
gatk_num_threads=4
gatk_num_cpu_threads=8

##-------------
##Step0-4: Input Arguments
##-------------

no_geno=0
while test $# -gt 0 ; do
    case "$1" in
        -h|--help)
            echo ""
            echo "Usage: bash $0 [options] sample_name"
            echo ""
            echo "This script performs the entire variant-calling process upon one sample,"
            echo "following the Genome Analysis Toolkit (GATK)'s pipeline."
            echo ""
            echo "Options:"
            echo "-h, --help          display this help and exit"
            echo "-v, --version       display version of this script and exit"
            echo "-XS, --no-summary   suppress the command summary before execution"
            echo "-XP, --no-prompt    suppress the user prompt before execution, only when the command summary is displayed"
            echo "-XX, --no-exec      suppress automatic execution, generating only script files"
            echo "-e, --exome         limit variant detection to certain regions, requires bed-file"
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
            bed_argument='-L '/data/${exon_bed}
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
if [[ ! -e /${fastq_dir}/${sample_name}_1.fastq.gz ]] ; then
    echo
    echo 'Invalid SAMPLE NAME: '${sample_name}
    echo /${fastq_dir}/${sample_name}_1.fastq.gz not found.
    echo 'Terminated.'
    echo
    exit 1
fi
if [[ ! -e /${fastq_dir}/${sample_name}_2.fastq.gz ]] ; then
    echo
    echo 'Invalid SAMPLE NAME: '${sample_name}
    echo /${fastq_dir}/${sample_name}_2.fastq.gz not found.
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
mkdir -p /${out_dir}
mkdir -p /${out_dir}/${sample_name} ; mkdir -p /${out_dir}/${sample_name}/{Scripts,LOG,TEMP,SAM,BAM,BQSR,GVCF,VCF,QC,QC/FILTERED,Report}
chmod -R 777 /${out_dir}/${sample_name}


##-------------
##Step1: Align
##-------------
cat << EOL > /${out_dir}/${sample_name}/Scripts/1_${sample_name}_align.sh
#!/bin/bash
set -e
##-------------
##Step1: Align
##-------------
docker run --rm -v /:/data biodckr/bwa bwa mem \
-t ${cores} \
-R "@RG\tID:DM_${sample_name}\tSM:${sample_name}\tPL:Illumina\tLB:WES\tPU:unit1" \
/data/${ref_genome}.gz \
/data/${fastq_dir}/${sample_name}_1.fastq.gz \
/data/${fastq_dir}/${sample_name}_2.fastq.gz \
> /${out_dir}/${sample_name}/SAM/${sample_name}_aligned.sam

EOL



##-------------
##Step2: Sort
##-------------
cat <<EOL > /${out_dir}/${sample_name}/Scripts/2_${sample_name}_sort.sh
#!/bin/bash
set -e
##-------------
##Step2: Sort
##-------------
docker run --rm -v /:/data biodckr/picard /opt/conda/jre/bin/java -Xmx30G -jar /opt/conda/share/picard-2.3.0-0/picard.jar \
SortSam \
INPUT=/data/${out_dir}/${sample_name}/SAM/${sample_name}_aligned.sam \
OUTPUT=/data/${out_dir}/${sample_name}/BAM/${sample_name}_sorted.bam \
SORT_ORDER=coordinate \
TMP_DIR=/tmp

EOL



##-------------
##Step3: Deduplicate
##-------------
cat <<EOL > /${out_dir}/${sample_name}/Scripts/3_${sample_name}_deduplicate.sh
#!/bin/bash
set -e
##-------------
##Step3: Deduplicate
##-------------
docker run --rm -v /:/data biodckr/picard /opt/conda/jre/bin/java -Xmx30G -jar /opt/conda/share/picard-2.3.0-0/picard.jar \
MarkDuplicates \
INPUT=/data/${out_dir}/${sample_name}/BAM/${sample_name}_sorted.bam \
OUTPUT=/data/${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam \
METRICS_FILE=/data/${out_dir}/${sample_name}/BAM/${sample_name}_deduplication_metrics.txt \
CREATE_INDEX=TRUE \
TMP_DIR=/tmp

EOL



##-------------
##Step4: Build Index
##-------------
cat <<EOL > /${out_dir}/${sample_name}/Scripts/4_${sample_name}_build_index.sh
#!/bin/bash
set -e
##-------------
##Step4: Build Index
docker run --rm -v /:/data biodckr/picard /opt/conda/jre/bin/java -Xmx30G -jar /opt/conda/share/picard-2.3.0-0/picard.jar \
BuildBamIndex \
INPUT=/data/${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam \
TMP_DIR=/tmp

EOL



##-------------
##Step6: Base Quality Score Recalibration
##-------------
cat <<EOL > /${out_dir}/${sample_name}/Scripts/6_${sample_name}_recalibrate_base.sh
#!/bin/bash
set -e
##-------------
##Step6-1: Perform Base Recalibration
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T BaseRecalibrator \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R /data/${ref_genome} \
-knownSites /data/${indel_1} \
-knownSites /data/${indel_2} \
-knownSites /data/${DBSNP} \
${bed_argument} \
--interval_padding 100 \
-I /data/${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam \
-nct ${gatk_num_cpu_threads} \
-o /data/${out_dir}/${sample_name}/BQSR/${sample_name}_perform_bqsr.table \
-log /data/${out_dir}/${sample_name}/LOG/6-1_${sample_name}_perform_bqsr.log
#--fix_misencoded_quality_scores

##-------------
##Step6-4: Print Reads
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T PrintReads \
-R /data/${ref_genome} \
--disable_auto_index_creation_and_locking_when_reading_rods \
-I /data/${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam \
-BQSR /data/${out_dir}/${sample_name}/BQSR/${sample_name}_perform_bqsr.table \
-dt NONE \
-EOQ \
-nct ${gatk_num_cpu_threads} \
-o /data/${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
-log /data/${out_dir}/${sample_name}/LOG/6-4_${sample_name}_final_bam.log

EOL


##-------------
##Step7: Call Haplotype
##-------------
cat <<EOL > /${out_dir}/${sample_name}/Scripts/7_${sample_name}_call_haplotype.sh
#!/bin/bash
set -e
##-------------
##Step7: Call Haplotype
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T HaplotypeCaller \
-R /data/${ref_genome} \
--input_file /data/${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
--emitRefConfidence GVCF \
--genotyping_mode DISCOVERY \
-stand_emit_conf 30 \
-stand_call_conf 30 \
${bed_argument} \
--interval_padding 100 \
-o /data/${out_dir}/${sample_name}/GVCF/${sample_name}_GATK.g.vcf \
-log /data/${out_dir}/${sample_name}/LOG/7_${sample_name}_haplotype_caller.log \
-A DepthPerSampleHC \
-pairHMM VECTOR_LOGLESS_CACHING \
-nct ${gatk_num_cpu_threads}

EOL



##-------------
##Step8: Genotype
##-------------
cat <<EOL > /${out_dir}/${sample_name}/Scripts/8_${sample_name}_genotype_gvcf.sh
#!/bin/bash
set -e
##-------------
##Step8: Genotype
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T GenotypeGVCFs \
-R /data/${ref_genome} \
--variant /data/${out_dir}/${sample_name}/GVCF/${sample_name}_GATK.g.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-nt 1 \
-o /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW.vcf \
-log /data/${out_dir}/${sample_name}/LOG/8_${sample_name}_genotype_gvcf.log

docker run --rm -v /:/data alexcoppe/gatk \
-T VariantAnnotator \
-R /data/${ref_genome} \
--variant /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-I /data/${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam \
--dbsnp /data/${DBSNP} \
-L /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW.vcf \
-dt NONE \
-nt 1 \
-o /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW_ANNOTATED.vcf \
-log /data/${out_dir}/${sample_name}/LOG/8-1_${sample_name}_QC_snv_annotation.log \
-A ClippingRankSumTest \
-A ReadPosRankSumTest \
-A MappingQualityRankSumTest \
-A GCContent \
-A AlleleBalanceBySample \
-A AlleleBalance \
-A VariantType

EOL

##-------------
##Step9: SNV Quality Control
##-------------
cat <<EOL > /${out_dir}/${sample_name}/Scripts/9_${sample_name}_SNV_quality_control.sh
#!/bin/bash
set -e
##-------------
##Step9-1-1: Extract SNPs
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T SelectVariants \
-R /data/${ref_genome} \
--variant /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW_ANNOTATED.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType SNP \
--excludeFiltered \
-nt 1 \
-o /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW_SNV.vcf \
-log /data/${out_dir}/${sample_name}/LOG/9-1-1_${sample_name}_QC_select_snv.log

##-------------
##Step9-1-2: Filter SNPs
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T VariantFiltration \
-R /data/${ref_genome} \
--variant /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW_SNV.vcf \
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
-o /data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.vcf \
-log /data/${out_dir}/${sample_name}/LOG/9-1-2_${sample_name}_QC_filter_snv.log

##-------------
##Step9-2-1: Extract INDELs
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T SelectVariants \
-R /data/${ref_genome} \
--variant /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW_ANNOTATED.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType INDEL \
-selectType MNP \
-selectType MIXED \
-selectType SYMBOLIC \
--excludeFiltered \
-nt 1 \
-o /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW_INDEL.vcf \
-log /data/${out_dir}/${sample_name}/LOG/9-2-1_${sample_name}_QC_select_INDEL.log

##-------------
##Step9-2-2: Filter INDELs
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T VariantFiltration \
-R /data/${ref_genome} \
--variant /data/${out_dir}/${sample_name}/VCF/${sample_name}_RAW_INDEL.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--filterExpression 'QD < 2.0' \
--filterName 'QD' \
--filterExpression 'FS > 40.0' \
--filterName 'FS' \
--filterExpression 'ReadPosRankSum < -8.0' \
--filterName 'ReadPosRankSum' \
--filterExpression 'DP < 8.0' \
--filterName 'DP' \
--logging_level ERROR \
-o /data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.vcf \
-log /data/${out_dir}/${sample_name}/LOG/9-2-2_${sample_name}_QC_filter_indel.log

EOL

##-------------
##MASTER SCRIPT
##-------------
cat <<EOL > /${out_dir}/${sample_name}/Scripts/${sample_name}_GATK.sh
#!/bin/bash
set -e
##-------------
##${sample_name}'s Variant Calling
## Script version: ${version}
##-------------
(bash /${out_dir}/${sample_name}/Scripts/1_${sample_name}_align.sh) 2>&1 | tee /${out_dir}/${sample_name}/LOG/1_${sample_name}_alignment.log
(bash /${out_dir}/${sample_name}/Scripts/2_${sample_name}_sort.sh) 2>&1 | tee /${out_dir}/${sample_name}/LOG/2_${sample_name}_sort.log
(bash /${out_dir}/${sample_name}/Scripts/3_${sample_name}_deduplicate.sh) 2>&1 | tee /${out_dir}/${sample_name}/LOG/3_${sample_name}_deduplication.log
(bash /${out_dir}/${sample_name}/Scripts/4_${sample_name}_build_index.sh) 2>&1 | tee /${out_dir}/${sample_name}/LOG/4_${sample_name}_building_index.log

if [[ -e /${out_dir}/${sample_name}/BAM/${sample_name}_deduplicated.bam ]] ; then
    rm /${out_dir}/${sample_name}/SAM/${sample_name}_aligned.sam
fi

bash /${out_dir}/${sample_name}/Scripts/6_${sample_name}_recalibrate_base.sh

if [[ -e /${out_dir}/${sample_name}/BAM/${sample_name}_GATK.bam ]]; then
    rm -f /${out_dir}/${sample_name}/BAM/${sample_name}_{deduplicated,sorted,realigned}.{bam,bai}
fi

bash /${out_dir}/${sample_name}/Scripts/7_${sample_name}_call_haplotype.sh

if [[ ${no_geno} != 1 ]] ; then
    bash /${out_dir}/${sample_name}/Scripts/8_${sample_name}_genotype_gvcf.sh
    bash /${out_dir}/${sample_name}/Scripts/9_${sample_name}_SNV_quality_control.sh
fi

EOL




##-------------
##EXECUTION
##-------------
if [[ ${no_exec} != 1 ]] ; then
    bash /${out_dir}/${sample_name}/Scripts/${sample_name}_GATK.sh
fi
