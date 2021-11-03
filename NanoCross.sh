#!/bin/bash -e

usage() { echo "Usage: $0 -I <ONT.fastq> -R <reference_file> -O <Output> -C <Chromosome> [-wt]" 1>&2; 
	echo "  -I   Bulk-gamete sequencing use Oxford Nanopore" 1>&2;
	echo "  -R   Reference genome file" 1>&2;
	echo "  -O   Output file containing recombinant molecules " 1>&2;
	echo "  -C   Chromosome information" 1>&2;
	echo "------------optional--------";
	echo "  -w   Work Directory; default=Current Directory" 1>&2;
	echo "  -t   Number of threads; default=10" 1>&2;
	exit 1; }
	
while getopts ":I:R:O:C:w:t:" option; do
    case "${option}" in
	I)
            sample=${OPTARG}
            ;;
	R)
            reference=${OPTARG}
            ;;
	O)
			output=${OPTARG}
            ;;
	C)
            chromosome=${OPTARG}
            ;;
        w)
            work_dir=${OPTARG}
            ;;
		t)
            threads=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if  [[ -z "${sample}" || -z "${reference}" || -z "${output}" || -z "${chromosome}" ]] ; then
        usage
fi
if [ -z "${work_dir}" ] ; then
        work_dir=./
fi
if [ -z "${threads}" ] ; then
	threads=10
fi

#
# Step 1 Canu correct raw ONT Reads
#

Corr_Reads_dir=$work_dir/corr_reads
mkdir -p $Corr_Reads_dir
cd $Corr_Reads_dir 
${canu} -correct -p canu -d $sample  useGrid=false  -nanopore-raw $sample  genomeSize=2.5g \
corMinCoverage=0 corMhapSensitivity=high minReadLength=5000 minOverlapLength=500 corOutCoverage=120 \
 > canu.$sample.txt 2> error.canu.$sample.txt
 
# 
# Step 2 dehomopolymerate compression corrected ONT reads and Reference genome
#
Deho_reads_dir=$work_dir/Deho_dir
mkdir -p $Deho_dir
cd $Deho_dir
dehomopolymerate -f  $Corr_Reads_dir/$sample/canu.correctedReads.fasta.gz > $Deho_dir/$sample.deho.fasta \
> output_deho_$sample.txt 2> error_deho_$sample.txt

awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' $reference > $Deho_dir/reference.fasta
dehomopolymerate -f  $Deho_dir/reference.fasta > $Deho_dir/reference.deho.fasta \
> output_deho_$reference.txt 2> error_deho_$reference.txt

#
# Step 3 Align reads & Variants call
#

Bam_dir=$work_dir/Bam_dir
mkdir -p $Bam_dir
cd $Bam_dir
minimap2 -t 20  -ax map-ont -a  $Deho_dir/reference.deho.fasta $Deho_dir/$sample.deho.fasta | samtools view  -Shb  > $bam_dir/$sample.bam
samtools sort -o $Bam_dir/$sample.sort.bam  -@ $threads $Bam_dir/$sample.bam
samtools index $sample.sort.bam 
samtools view -bF 2308 -q 60 $sample.sort.bam > $sample.uni.bam 
samtools index $sample.uni.bam

Vcf_dir=$work_dir/Vcf_dir
mkdir -p $Vcf_dir
cd $Vcf_dir
model=/your_ont_model_path
python clair.py callVarBam --chkpnt_fn $model/model --bam_fn $Bam_dir/$sample.uni.bam \
--ref_fn $Deho_dir/reference.deho.fasta --call_fn $Vcf_dir/clair.$sample.vcf \
> calir.var.$sample.txt >2 clair.error.$sample.txt 

igvtools index $Vcf_dir/clair.$sample.vcf
gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantFiltration -O clair.$sample.filtered.vcf  \
-V  clair.$sample.vcf --filter-expression 'QUAL < 625' --filter-name lowQualFilter --cluster-window-size 10 \
 --cluster-size 3 --missing-values-evaluate-as-failing

#
# Step 4 Haplotype construction
#

Phase_dir=$work_dir/Phase_dir
mkdir -p $Phase_dir
cd $Phase_dir

whatshap phase -o $Phase_dir/whatshap.sample.phase.vcf --max-coverage 20 --reference=$Deho_dir/reference.deho.fasta \
$Vcf_dir/clair.$sample.filtered.vcf  $Bam_dir/$sample.uni.bam  \
> whatshap_phase.$sample.txt 2> error_phase.$sample.txt 
# Extraction of haplotype block information 
Rscript Phase_Extract.R  $phase_dir/whatshap.sample.phase.vcf $chromosome

#
# Step 5 Detection of recombinant molecules
#
Recom_dir=$work_dir/Recom_dir
mkdir -p $Recom_dir
cd $Recom_dir

Rscript Detect_Recombination.R $Bam_dir/$sample.uni.bam $Phase_dir/$chromosome.phase.txt 20 ./$chromosome.recom.txt 
