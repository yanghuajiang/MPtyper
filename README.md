# MPtyper
A tool for identifying the genotypes of Mycoplasma pneumoniae contained in sequence data.
## MPtyper was developed and tested in Python 3.9 and calls the following 3rd party programs for different pipelines:
~~~~~~~~~~~~~~
minimap2
samtools
~~~~~~~~~~~~~~
# Quick Start (with examples)
~~~~~~~~~~~~~~
python MPtyper.py -db examples/MP -c -b -o examples/test -r examples/r1.fq.gz -r examples/r2.fq.gz
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
  --cutoff FLOAT     cutoff of total depth in type-specific sites for validation of genotypes.
                     [default:1, to assign genotypes for low-quality data]
  -b, --bam          flag to keep intermediate BAM file.
  --help             Show this message and exit.
~~~~~~~~~~~~~~
# Notes
## Scope of use
MPtyper is capable of analyzing single-end and paired-end sequencing data of single-clone bacteria or next-generation metagenome. However, because of the predefined parameters of minimap2 (-ax sr), it is currently applicable better to short-read sequencing data.
## Composition of the database
The database is a folder containing a reference genome in fasta format, SNPs table(s) with ".def" as suffix in format as "genotype  contig  site  point-mutation" as follows: [please keep sites in order for each genotype to reduce the caculating time]
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
### The core information of genotypes with its possibility will print on the screen, including [prefix], [P1 and EC type], [barcode type], [MR mutation]:
~~~~~~~~~~~~~~
test	P1-1(1.00);EC1(1.0)	Barcode_1.3.1(1.0)	MR_A2063G(1.0)
~~~~~~~~~~~~~~
### There are 1-3 outputs stored detailed genotype information, consensus sequences and matched result in BAM format:
~~~~~~~~~~~~~~
examples/test.genotypes
examples/test.consensus.fa [if "-c" was given]
examples/test.bam (.gz format) [if "-b" was given]
~~~~~~~~~~~~~~



