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
[output print: test	P1-1(1.00);EC1(1.00)	Barcode_1.3.1(1.00)	MR_A2063G]
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
MPtyper is capable of analyzing single-end and paired-end sequencing data of single-clone bacteria or next-generation metagenome. However, because of the predefined parameters of minimap2 (-ax sr), it is currently applicable better to short-read sequencing data.
## Composition of the database
The database is a folder containing a reference genome in fasta format, SNVs table(s) with ".def" as suffix in format as "genotype  contig  site  point-mutation" as follows: [please keep sites in order for each genotype to reduce the caculating time]
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
test	P1-1(1.00);EC1(1.00)	Barcode_1.3.1(1.00)	MR_A2063G
~~~~~~~~~~~~~~
### There are 1-3 outputs stored detailed genotype information, consensus sequences and matched result in BAM format:
~~~~~~~~~~~~~~
examples/test.genotypes
examples/test.consensus.fa [if "-c" was given]
examples/test.bam [if "-b" was given]
~~~~~~~~~~~~~~
### The number and the order of columns in the ".genotypes" file is determined by the number and the name's order of SNVs tables in the database folder, and multiple genotypes can be stored within the same table.
### It gives two kind of information for each genotype: genotype name, matched reads and matched rate. Each type-specific SNVs present two genotypes: the genotype represented by the original base and represented by the substituted base, which are directly opposite. 
However, in the genotyping of Mycoplasma pneumoniae genomes, only P1-1 and P1-2 are directly opposite. EC1 is merely a prevalent clone within P1-1, and EC2 is merely a prevalent clone within P1-2. Therefore, EC1 and EC2 are not directly opposite. In other words, the direct opposite of EC1 is "non-EC1." This interpretation can be extended to all other SNP-based genotyping scenarios. 
In summary, the magnitude of the probability value represents the ratio of the number of reads supporting that genotype to the total number of matched reads for both genotypes, with a maximum value of 1. 
All detailed information for each genotype will saved in the ".genotypes" file: 
~~~~~~~~~~~~~~
P1-2(15,0.00) means the possibility of "P1-2" is 0, 15 reads supports the opposite genotype (P1-1).
P1-2(445,0.15) means the possibility of "P1-2" is 0.15, 445 reads supports this genotype, which also means the possibility of the opposite genotype (P1-1) is 0.85. So it is more likely to be "P1-1". The core information printed on the screen will simply show "P1-1(0.85)".
EC1(138,1.00) means the possibility of "EC1" is 1, 138 reads supports this genotype.
EC2(0,-1) means there is no reads supporting both "EC2" and "non-EC2", the possibility of "-1" means lack of information, which is set to distinct against the possibility of "0".
~~~~~~~~~~~~~~

