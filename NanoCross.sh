#!/bin/bash -e

usage() { echo "Usage: $0 -I <ONT.fastq> -R <reference_file> -O <Output> -C <Chromosome> [-wt]" 1>&2; 
	echo "  -I   Bulk-gamete sequencing use Oxford Nanopore" 1>&2;
	echo "  -N   sample name of input file" 1>&2;
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
            input=${OPTARG}
            ;;
    N)
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

if  [[ -z "input" || -z "${sample}" ||-z "${reference}" || -z "${output}" || -z "${chromosome}" ]] ; then
        usage
fi
if [ -z "${work_dir}" ] ; then
        work_dir=./
fi
if [ -z "${threads}" ] ; then
	threads=10
fi


# 
# Step 1 dehomopolymerate compression corrected ONT reads and Reference genome
#

NanoCross_dir=./
Deho_reads_dir=$work_dir/Deho_dir
mkdir -p $Deho_dir
cd $Deho_dir

dehomopolymerate -f  $input > $Deho_dir/$sample.deho.fasta \
> output_deho.txt 2> error_deho.txt

genome_dir=$work_dir/genome
mkdir -p $genome_dir
cd $genome_dir
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' $reference > $genome_dir/reference.upper.fasta

dehomopolymerate -f  $genome_dir/reference.upper.fasta > $genome_dir/reference.upper.deho.fasta \
> output_deho_$reference.txt 2> error_deho_$reference.txt
samtools faidx reference.upper.deho.fasta

#
# Step 2 Align reads & Variants call
#

Bam_dir=$work_dir/Bam_dir
mkdir -p $Bam_dir
cd $Bam_dir
minimap2 -t 20  -ax map-ont -a  $genome_dir/reference.upper.deho.fasta  $Deho_dir/$sample.deho.fasta | samtools view  -Shb  > $bam_dir/$sample.deho.bam
samtools sort -o $Bam_dir/$sample.deho.sort.bam  -@ $threads $Bam_dir/deho.bam
samtools index $sample.deho.sort.bam 
samtools view -bF 2308 -q 60 $sample.deho.sort.bam $chromosome > $sample.$chromosome.deho.filtered.bam 
samtools index $sample.deho.filtered.bam

Vcf_dir=$work_dir/Vcf_dir
mkdir -p $Vcf_dir
cd $Vcf_dir
model=/your_ont_model_path
python clair.py callVarBam --chkpnt_fn $model/model --bam_fn $bam_dir/$sample.$chromosome.deho.filtered.bam \
--ref_fn $genome_dir/reference.upper.deho.fasta --call_fn $Vcf_dir/$sample.$chromosome.vcf \
> calir.var.$sample.$chromosome.txt >2 clair.error.$sample.$chromosome.txt 

igvtools index $Vcf_dir/clair.$sample.vcf
gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantFiltration -O $sample.$chromosome.filtered.vcf  \
-V  $sample.$chromosome.vcf --filter-expression 'QUAL < 650' --filter-name lowQualFilter --cluster-window-size 10 \
 --cluster-size 3 --missing-values-evaluate-as-failing
#Extract the variation through filtering
awk '/^#/{print $0}'  $sample.$chromosome.filtered.vcf  > $sample.$chromosome.passed.vcf 
awk '$7=="PASS"'  $sample.$chromosome.filtered.vcf  >> $sample.$chromosome.passed.vcf 

#
# Step 3 Haplotype construction
#

Phase_dir=$work_dir/Phase_dir
mkdir -p $Phase_dir
cd $Phase_dir

whatshap phase -o $Phase_dir/$sample.$chromosome.phase.vcf --max-coverage 20 --reference=$genome_dir/reference.upper.deho.fasta \
$Vcf_dir/$sample.$chromosome.passed.vcf  $Bam_dir/$sample.$chromosome.deho.filtered.bam  \
> whatshap_phase.$sample.txt 2> error_phase.$sample.txt 

# Extraction of haplotype block information 
Rscript $NanoCross_dir/Phase_Extract.R  $Phase_dir/$sample.$chromosome.phase.vcf $chromosome

#
# Step 4 Detection of recombinant molecules
#

Recom_dir=$work_dir/Recom_dir
mkdir -p $Recom_dir
cd $Recom_dir

Rscript $NanoCross_dir/Detect_Recombination.R $Bam_dir/$sample.$chromosome.deho.filtered.bam $Phase_dir/$chromosome.phase.txt $threads ./$sample$chromosome.recom.txt 
