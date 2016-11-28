#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script performs the entire joint variant-calling process upon all samples, following the Genome Analysis Toolkit (GATK)'s pipeline.
##If no \"sample_name\" is specified in the command, all the folders in \"out_dir\" are automatically taken as the complete set of samples.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Check if the individual GATK process for every sample has been accomplished.
##iii>	Run the command 'bash /path/to/GATK_joint.sh [options] [sample_name]'
##-------------


version=0.1


##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
map_dir=/
ref_dir=tiger/harald/resources/hg38bundle
out_dir=tiger/harald/N-S-T

##-------------
##Step0-2: References
##-------------
ref_genome=${ref_dir}/Homo_sapiens_assembly38.fasta
indel_1=${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
indel_2=${ref_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz
DBSNP=${ref_dir}/dbsnp_144.hg38.vcf.gz
exon_bed=${ref_dir}/remapped_agilentV5_regions_onlychr.bed
file_list=${out_dir}/GVCF_list.list

##-------------
##Step0-3: Other Parametres
##-------------
java_mem=10G
cores=16
gatk_num_threads=4
gatk_num_cpu_threads=8

##-------------
##Step0-4: Input Arguments
##-------------
while test $# -gt 0 ; do
    case "$1" in
        -h|--help)
            echo ""
            echo "Usage: bash $0 [options] [sample_name] [sample_name] [...]"
            echo ""
            echo "This script performs the entire joint variant-calling process upon all samples,"
            echo "following the Genome Analysis Toolkit (GATK)'s pipeline."
            echo ""
            echo "Options:"
            echo "-h, --help           display this help and exit"
            echo "-v, --version        display version of this script and exit"
            echo "-XS, --no-summary    suppress the command summary before execution"
            echo "-XP, --no-prompt     suppress the user prompt before execution, only when the command summary is displayed"
            echo "-XX, --no-exec       suppress automatic execution, generating only script files"
            echo "-e, --exome          call only exonic variants, drastically accelerating the Joint Genotype process"
            echo ""
            exit 0
            ;;
        -v|--version)
            echo ""
            echo "GATK_joint.sh"
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
        -e|--exome)
            seq_type='EXOME'
            bed_argument='-L /data/'${exon_bed}
            shift
            ;;
        *)
            joint_name=$1
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
if [[ ! -v joint_name ]] ; then
    joint_name='COMBINED'
fi

##-------------
##Step0-7: Summarisation & User's Confirmation Prompt
##-------------
if [[ ${no_summary} != 1 ]] ; then
    echo
    echo '---------------------------------------'
    echo 'JOINT VARIANT CALLING PROCESS'
    echo 'SEQUENCED DATA =     '${seq_type}
    echo 'PREFIX =             '${joint_name}
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
mkdir -p /${out_dir}/${joint_name}
mkdir -p /${out_dir}/${joint_name}/{Script,LOG,VCF,VQSR}
chmod -R 777 /${out_dir}/${joint_name}


##-------------
##Step2: Joint Genotype
##-------------
cat <<EOL > /${out_dir}/${joint_name}/Script/2_${joint_name}_joint_genotype.sh
#!/bin/bash

set -e

##-------------
##Step2: Joint Genotype
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T GenotypeGVCFs \
-R /data/${ref_genome} \
--variant /data/${out_dir}/GVCF_list.list \
-nt ${gatk_num_threads} \
${bed_argument} \
-dt NONE \
-o /data/${out_dir}/${joint_name}/VCF/${joint_name}_RAW.vcf \
-log /data/${out_dir}/${joint_name}/LOG/2-1_${joint_name}_genotype_gvcf.log

docker run --rm -v /:/data alexcoppe/gatk \
-T VariantAnnotator \
-R /data/${ref_genome} \
--variant /data/${out_dir}/${joint_name}/VCF/${joint_name}_RAW.vcf \
--dbsnp /data/${DBSNP} \
-L /data/${out_dir}/${joint_name}/VCF/${joint_name}_RAW.vcf \
-nt ${gatk_num_threads} \
-o /data/${out_dir}/${joint_name}/VCF/${joint_name}_RAW_ANNOTATED.vcf \
-log /data/${out_dir}/${joint_name}/LOG/2-2_${joint_name}_QC_snv_annotation.log \
-A ClippingRankSumTest \
-A ReadPosRankSumTest \
-A MappingQualityRankSumTest \
-A GCContent \
-A AlleleBalanceBySample \
-A AlleleBalance \
-A VariantType

EOL


##-------------
##Step3: Variant Quality Score Recalibration
##-------------
cat <<EOL > /${out_dir}/${joint_name}/Script/3_${joint_name}_recalibrate_variant.sh
#!/bin/bash

set -e

##-------------
##Step3-1: Select SNVs
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T VariantRecalibrator \
-R /data/${ref_genome} \
-input /data/${out_dir}/${joint_name}/VCF/${joint_name}_RAW_ANNOTATED.vcf \
-recalFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.recal \
-tranchesFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.tranches \
-nt ${gatk_num_threads} \
-dt NONE \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/${ref_dir}/hapmap_3.3.hg38.vcf.gz \
-resource:omni,known=false,training=true,truth=true,prior=12.0 /data/${ref_dir}/1000G_omni2.5.hg38.vcf.gz \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 /data/${ref_dir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/${DBSNP} \
-an QD \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
-an InbreedingCoeff \
-mode SNP \
-rscriptFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.R \
-log /data/${out_dir}/${joint_name}/LOG/3-1-1_${joint_name}_VQSR_snv_recalibration.log

docker run --rm -v /:/data alexcoppe/gatk \
-T ApplyRecalibration \
-R /data/${ref_genome} \
-input /data/${out_dir}/${joint_name}/VCF/${joint_name}_RAW_ANNOTATED.vcf \
-tranchesFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.tranches \
-recalFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.recal \
-o /data/${out_dir}/${joint_name}/VCF/${joint_name}_SNV.recalibrated.filtered.vcf \
--ts_filter_level 99.5 \
-mode SNP \
-log /data/${out_dir}/${joint_name}/LOG/3-1-2_${joint_name}_VQSR_apply_snv_recalibration.log

##-------------
##Step3-2: Select Indels
##-------------
docker run --rm -v /:/data alexcoppe/gatk \
-T VariantRecalibrator \
-R /data/${ref_genome} \
-input /data/${out_dir}/${joint_name}/VCF/${joint_name}_RAW_ANNOTATED.vcf \
-recalFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.recal \
-tranchesFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.tranches \
-nt ${gatk_num_threads} \
-dt NONE \
--maxGaussians 4 \
-resource:mills,known=false,training=true,truth=true,prior=12.0 /data/${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/${DBSNP} \
-an QD \
-an FS \
-an SOR \
-an ReadPosRankSum \
-an MQRankSum \
-an InbreedingCoeff \
-mode INDEL \
-rscriptFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.R \
-log /data/${out_dir}/${joint_name}/LOG/3-2-1_${joint_name}_VQSR_indel_recalibration.log

docker run --rm -v /:/data alexcoppe/gatk \
-T ApplyRecalibration \
-R /data/${ref_genome} \
-input /data/${out_dir}/${joint_name}/VCF/${joint_name}_RAW_ANNOTATED.vcf \
-tranchesFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.tranches \
-recalFile /data/${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.recal \
-o /data/${out_dir}/${joint_name}/VCF/${joint_name}_INDEL.recalibrated.filtered.vcf \
--ts_filter_level 99.0 \
-mode INDEL \
-log /data/${out_dir}/${joint_name}/LOG/3-2-2_${joint_name}_VQSR_apply_indel_recalibration.log

EOL


##-------------
##MASTER SCRIPT
##-------------
cat <<EOL > /${out_dir}/${joint_name}/Script/${joint_name}_GATK.sh
#!/bin/bash
##-------------
##Joint Variant Calling
## Version: ${version}
##-------------
#bash /${out_dir}/${joint_name}/Script/1_${joint_name}_combine_vcfs.sh
bash /${out_dir}/${joint_name}/Script/2_${joint_name}_joint_genotype.sh
bash /${out_dir}/${joint_name}/Script/3_${joint_name}_recalibrate_variant.sh
#bash /${out_dir}/${joint_name}/Script/4_${joint_name}_SNV_quality_control.sh

EOL





##-------------
##EXECUTION
##-------------
if [[ ${no_exec} != 1 ]] ; then
    bash /${out_dir}/${joint_name}/Script/${joint_name}_GATK.sh
fi
