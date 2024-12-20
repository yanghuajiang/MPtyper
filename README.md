# MPtyper
A tool for identifying the genotypes of Mycoplasma pneumoniae contained in short-read data.
## MPtyper was developed and tested in Python 3.9. MPtyper depends on several Python libraries:
~~~~~~~~~~~~~~
click
glob
shutil
subprocess
collections
re
numpy
~~~~~~~~~~~~~~
## MPtyper also calls the following 3rd party programs for different pipelines:
~~~~~~~~~~~~~~
minimap2
samtools
~~~~~~~~~~~~~~
# Quick Start (with examples)
~~~~~~~~~~~~~~
python MPtyper.py -db examples/MP -c -b -o examples/test_out -r examples/r1.fq.gz -r examples/r2.fq.gz
~~~~~~~~~~~~~~
# USAGE
~~~~~~~~~~~~~~
Usage: MPtyper.py [OPTIONS]

Options:
  Internally available: <>

  -db, --database TEXT   dirname for the database. [default: MP]
  -o, --prefix TEXT  prefix for the outputs.  [required]
  -r, --reads TEXT   files for short reads, can be specified at most twice. 
                     [required]
  -c, --consensus    flag to generate consensus sequences. (for phylogenetic
                     analysis)
  --min_depth FLOAT  minimum read depth for high quality bases. [only for
                     consensus sequences, default:3]
  --min_cons FLOAT   minimum proportion of consensus reads for high quality
                     bases. [only for consensus sequences, default:0.8]
  -b, --bam          flag to keep intermediate BAM file.
  --help             Show this message and exit.
~~~~~~~~~~~~~~
# Notes
## Scope of use
MPtyper is capable of analyzing single-end and paired-end sequencing data of single-clone bacteria or next-generation metagenome. However, because of the predefined parameters of minimap2, it is currently applicable only to short-read sequencing data.
## Composition of the database
The database is a folder containing a reference genome in fasta format, SNP table(s) with ".def" as suffix in format as "genotype  contig  site  point-mutation" as follows: [keep sites in order]
~~~~~~~~~~~~~~
......
EC1	LR214945.1	435750	C->A
EC1	LR214945.1	488617	C->A
EC1	LR214945.1	754986	C->T
EC1	LR214945.1	785483	C->T
EC2	LR214945.1	152818	G->A
EC2	LR214945.1	243011	C->T
EC2	LR214945.1	423818	C->T
EC2	LR214945.1	555564	C->A
......
~~~~~~~~~~~~~~
## Ootputs
### There are 1-3 outputs stored genotype information, consensus sequences and matched reads in BAM format:
~~~~~~~~~~~~~~
examples/test.genotypes
examples/test.consensus.fa [if "-c" was given]
examples/test.bam [if "-b" was given]
~~~~~~~~~~~~~~
### The number and the order of columns in the ".genotypes" file is determined by the number and the name's order of SNVs tables in the database folder, and multiple genotypes can be stored within the same table.
### Output will give two kind of information for each genotype: Genotype name, Liklihood. Each SNPs table presents two genotypes. The genotype represented by the original base and represented by the substituted base are directly opposing. However, in the genotyping of Mycoplasma pneumoniae genomes, only P1-1 and P1-2 are directly opposing. EC1 is merely a prevalent clone within P1-1, and EC2 is merely a prevalent clone within P1-2. Therefore, EC1 and EC2 are not directly opposing. In other words, the direct opposite of EC1 is "non-EC1." This interpretation can be extended to all other SNP-based genotyping scenarios. In summary, the magnitude of the probability value represents the ratio of the number of reads supporting that genotype to the total number of matching reads for both genotypes, with a maximum value of 1. The positive and negative signs indicate the genotype, where a positive sign represents the original genotype and a negative sign represents its opposing genotype. For instance: 
~~~~~~~~~~~~~~
P1-1(-1.00) means "non P1-1" with liklihood of 1
EC2(0.98) means "EC2" with liklihood of 0.98
~~~~~~~~~~~~~~

