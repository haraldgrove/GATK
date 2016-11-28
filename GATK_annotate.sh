#!/bin/bash
##-------------
##Original development:
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-----------
##Addition of docker images and use of hg38 reference
##Harald Grove, JUN 2016
##NOTICE:
##  The use of docker images creates some points to be aware of:
##      By default the root folder ('/') is being mapped to /data in the container
##      This means that root should be excluded from the given paths (e.g. out_dir)
##-------------
##This script annotates a given VCF file with the following databases:
##i>	dbSNP
##ii>	dbNSFP
##iii>	gwasCat
##iv>	PhastCons
##v>	ClinVar
##vi>	SnpEff
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Run the command 'bash /path/to/GATK_annotate.sh [options] input_file'
##
##snpEff produces two output files that is placed in the folder where the command is run.
##	snpEff_summary.html
##	snpEff_genes.txt
##-------------



set -e

##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
snpeff_dir=/home/biodocker/snpEff
out_dir=tiger/harald/SCD

##-------------
##Step0-2: References
##-------------
DBSNP=tiger/harald/resources/hg38bundle/dbsnp_144.hg38.vcf.gz
DBNSFP=${snpeff_dir}/data/dbNSFP3.2a.txt.gz
GWASCATALOG=${snpeff_dir}/data/gwas_catalog_v1.0.1-associations_e84_r2016-07-10.tsv
PHASTCONS=${snpeff_dir}/data/phastCons100way
CLINVAR=${snpeff_dir}/data/clinvar.vcf 

##-------------
##Step0-3: Other Parametres
##-------------
java_mem=10G

##-------------
##Step0-4: Input Arguments
##-------------
while test $# -gt 0 ; do
    case "$1" in
        -h|--help)
            echo ""
            echo "Usage: bash $0 [options] input_file"
            echo ""
            echo "This script annotates a given VCF file with the following databases:"
            echo "i>	dbSNP"
            echo "ii>	dbNSFP"
            echo "iii>	gwasCat"
            echo "iv>	PhastCons"
            echo "v>	ClinVar"
            echo "vi>	SnpEff"
            echo ""
            echo "Options:"
            echo "-h, --help           display this help and exit"
            echo "-v, --version        display version of this script and exit"
            echo "-XS, --no-summary    suppress the command summary before execution"
            echo "-XP, --no-prompt     suppress the user prompt before execution, only when the command summary is displayed"
            # echo "-O, --out-dir  OUT_DIR   specify output directory, the same as input's by default"
            # echo "-o, --out-file  OUT_FILE  specify output file, by default the input file's name with the suffix \"annotated\" appended"
            echo "-r, --replace        overwrite input VCF file after finished annotation"
            echo ""
            exit 0
            ;;
        -v|--version)
            echo ""
            echo "GATK_annotate.sh"
            echo ""
            echo "Created MAR 2016"
            echo "Updated JUL 2016"
            echo "by"
            echo "PURIN WANGKIRATIKANT [purin.wan@mahidol.ac.th]"
            echo "Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)"
            echo ""
            echo "This version uses a Docker image to store snpEff, SnpSift and most of the databases (except DBSNP)"
            echo "Updated NOV 2016"
            echo "by"
            echo "HARALD GROVE [harald.gro@mahidol.ac.th]"
            echo "Bioinformatics and Data Management for Research"
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
        -o|--out-file)
            shift
            out_file=$1
            out_dir=$( echo $1 | sed 's/\/[^\/]*$//')/.
            shift
            ;;
        -r|--replace)
            replacement=YES
            shift
            ;;
        # -O|--out-dir)
            # shift
            # out_dir=$( echo $1 | sed 's/\/$//' )
            # shift
            # ;;
        *)
            sample_name=$1
            shift
            ;;
    esac
done

##-------------
##Step0-5: Default Value Setting
##-------------
#if [[ ! -v out_dir ]] ; then
#    out_dir=$( pwd )
#fi
#if [[ ! -v out_file ]] ; then
#    out_file=$( echo ${in_file} | sed 's/.vcf$/_annotated.vcf/' )
#fi
if [[ ! -v replacement ]] ; then
    replacement=NO
fi

##-------------
##Step0-6: Input Verification
##-------------
in_file=/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.vcf
if [[ ! -e ${in_file} ]] ; then
    echo
    echo 'Invalid INPUT FILE: '${in_file}
    echo ${in_file} not found.
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
    echo 'VARIANT ANNOTATION'
    echo 'INPUT FILE =			'${in_file}
    echo 'OUTPUT FILE =			'${out_file}
    echo 'REPLACEMENT =			'${replacement}
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
##Step0-8: Output Folder Creation
##-------------
#mkdir -p /${out_dir}
#cp ${in_file} ${out_file}


cat << EOL > /${out_dir}/${sample_name}/Scripts/${sample_name}_annotate.sh
#!/bin/bash
# set -e

##-------------
##Step1: dbSNP
##-------------
echo '1/6 dbSNP Annotation Started'
docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
annotate \
/data/${DBSNP} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp1.vcf

docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
annotate \
/data/${DBSNP} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp1.vcf

echo '1/6 dbSNP Annotation Completed'



##-------------
##Step2: dbNSFP
##-------------
echo '2/6 dbNSFP Annotation Started'
docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
dbnsfp \
-db ${DBNSFP} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp1.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp2.vcf

docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
dbnsfp \
-db ${DBNSFP} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp1.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp2.vcf

echo '2/6 dbNSFP Annotation Completed'



##-------------
##Step3: gwasCat
##-------------
echo '3/6 gwasCat Annotation Started'
docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
gwasCat \
-db ${GWASCATALOG} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp2.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp3.vcf

docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
gwasCat \
-db ${GWASCATALOG} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp2.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp3.vcf

echo '3/6 gwasCat Annotation Completed'



##-------------
##Step4: PhastCons
##-------------
echo '4/6 PhastCons Annotation Started'
docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
phastCons \
${PHASTCONS} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp3.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp4.vcf

docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
phastCons \
${PHASTCONS} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp3.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp4.vcf

echo '4/6 PhastCons Annotation Completed'



##-------------
##Step5: ClinVar
##-------------
echo '5/6 ClinVar Annotation Started'
docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
annotate \
${CLINVAR} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp4.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp5.vcf

docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
annotate \
${CLINVAR} \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp4.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp5.vcf

echo '5/6 ClinVar Annotation Completed'

##-------------
##Step6: ESP6500
##-------------
#echo '6/7 ESP6500 Annotation Started'
#docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
#annotate \
#/home/biodocker/snpEff/data/ESP6500/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf \
#/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp5.vcf \
#> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp6.vcf
#
#docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx${java_mem} -jar /home/biodocker/snpEff/SnpSift.jar \
#annotate \
#/home/biodocker/snpEff/data/ESP6500/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf \
#/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp5.vcf \
#> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp6.vcf
#echo '6/7 ESP6500 Annotation Completed'

##-------------
##Step6: SnpEff
##-------------
echo '6/6 SnpEff Annotation Started'
docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx10G -jar /home/biodocker/snpEff/snpEff.jar \
GRCh38.82 \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp5.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp6.vcf

docker run --rm -v /:/data haraldgrove/snpeff:v2 java -Xmx10G -jar /home/biodocker/snpEff/snpEff.jar \
GRCh38.82 \
/data/${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp5.vcf \
> /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp6.vcf

echo '6/6 SnpEff Annotation Completed'



##-------------
##Step8: Clean-up
##-------------

mv /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp6.vcf /${out_dir}/${sample_name}/VCF/${sample_name}_SNV.vcf
mv /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp6.vcf /${out_dir}/${sample_name}/VCF/${sample_name}_INDEL.vcf

if [[ -e /${out_dir}/${sample_name}/VCF/${sample_name}_INDEL.vcf ]] ; then
    rm /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.temp?.vcf
fi

if [[ -e /${out_dir}/${sample_name}/VCF/${sample_name}_SNV.vcf ]] ; then
    rm /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.temp?.vcf
fi

if [[ ${replacement} == 'YES' ]] ; then
    rm /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_INDEL.vcf
    rm /${out_dir}/${sample_name}/VCF/${sample_name}_FILTERED_SNV.vcf
fi

EOL

##-------------
##EXECUTION
##-------------
if [[ ${no_exec} != 1 ]] ; then
    bash /${out_dir}/${sample_name}/Scripts/${sample_name}_annotate.sh
fi
