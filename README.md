# nanoCross
>A pipeline that Detecting recombination breakpoints using ONT sequence libraries.

Author: Zuoquan Chen

Email: zuoquanchen@outlook.com

Draft date: Nov. 1, 2021
## Dependencies

- g++ >= 4.8.5
- [dehomopolymerate](https://github.com/tseemann/dehomopolymerate)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [Canu](https://github.com/marbl/canu)
- [clair](https://github.com/HKU-BAL/Clair)
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
export PATH=$path_to_canu_bin/canu:$PATH
export PATH=$path_to_clair_bin/clair:$PATH
export PATH=$path_to_whatshap_bin/whatshap:$PATH
```

## Usage

NanoCross requires bam files and genomic haplotype files as input files and output of molecules where recombination exists.
```
NanoCross.sh -I <ONT.fastq> -R <reference_file> -O <Output> -C <Chromosome> [-wt]
  -I   Bulk-gamete sequencing use Oxford Nanopore
  -R   Reference genome file
  -O   Output file containing recombinant molecules
  -C   Chromosome information
------------optional--------
  -w   Work Directory; default=Current Directory
  -t   Number of threads; default=10
```

## Output

NanoCross will create these directories:
- Corr_Reads contains Corrected sequence 
- deho_reads contains dehomopolymerate sequence 
- Bam_dir contains align file 
- vcf_dir contains the  variant file 
- phase_dir contains phased variant
- recom_dir contains recombinant molecules information

## License

This project is licensed under the GNU General Public License - see the [LICENSE.md](LICENSE.md) file for details
