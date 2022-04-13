# nanoCross
>A pipeline that Detecting recombination crossover using ONT sequence libraries.

Author: Zuoquan Chen

Email: zuoquanchen@outlook.com

Draft date: Feb. 1, 2022
## Dependencies

- g++ >= 4.8.5
- [dehomopolymerate](https://github.com/tseemann/dehomopolymerate)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [clair](https://github.com/HKU-BAL/Clair)
- [igvtools](https://software.broadinstitute.org/software/igv/home)
- [GATK](https://github.com/broadinstitute/gatk)
- [Whatshap](https://github.com/whatshap/whatshap)
- [R >= 3.6.0](https://www.r-project.org/)

## Installation

```
git clone https://github.com/zuoquanchen/nanoCross.git
```
Then, To configure the other tools, add the following directories to your PATH:
```
export PATH=$path_to_dehomopolymerate_bin/dehomopolymerate:$PATH
export PATH=$path_to_minimap2_bin/minimap2:$PATH
export PATH=$path_to_samtools_bin/samtools:$PATH
export PATH=$path_to_clair_bin/clair:$PATH
export PATH=$path_to_igvtools_bin/igvtools:$PATH
export PATH=$path_to_whatshap_bin/whatshap:$PATH
export PATH=$path_to_Rscript_bin/Rscript:$PATH
```

## Usage

NanoCross requires bam files and genomic haplotype files as input files and output of molecules where recombination exists.
```
NanoCross.sh -I <ONT.fastq> -R <reference_file> -O <Output> -C <Chromosome> [-wt]
  -I   Bulk-gamete sequencing use Oxford Nanopore
  -N   sample name of input file
  -R   Reference genome file
  -O   Output file containing recombinant molecules
  -C   Chromosome information
------------optional--------
  -w   Work Directory; default=Current Directory
  -t   Number of threads; default=10
```

## Output

NanoCross will create these directories: 
- Deho_reads contains dehomopolymerate sequence 
- Bam_dir contains align file 
- Vcf_dir contains the  variant file 
- Phase_dir contains phased variant
- Recom_dir contains recombinant molecules and location information

## License

This project is licensed under the GNU General Public License - see the [LICENSE.md](LICENSE.md) file for details
