# Coalescent Simulator for Tumor Evolution (CSiTE)

CSiTE is a Python program for simulating short reads for tumor samples.
It takes a coalescent tree (in the newick format) and jointly simulates Single 
Nucleotide Variants (SNVs) and Copy Number Variants (CNVs) along the history of 
the genealogy. By integrating those somatic variants into the normal genome, 
perturbed genome of each tumor cell can be built, which are then piped to ART to 
simulate NGS short reads. Germline variants can also be integrated to the genome 
to make the data more realistic. With those charactors, CSiTE can be used to 
generate benchmark data for both bulk and single cell sample of tumor.

## 1. Installing

CSiTE is written in Python3 (>=3.5), and requires numpy, pyfaidx and yaml. 
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
is also required if you want to simulate short reads with CSiTE.

CSiTE can be downloaded from github:

    git clone https://github.com/hchyang/CSiTE.git

## 2. Running CSiTE

There are five modules in CSiTE. 

    Program: csite.py (a Coalescent Simulator for Tumor Evolution)
    Version: 0.9.0

    Usage:   csite.py <command> [options]

    Command: vcf2fa     build normal genome from input (germline) vcf file
             phylovar   simulate somatic variations on a phylogeny
             chain2fa   build tumor genomes from somatic mutations (chain file)
             fa2ngs     simulate short reads from normal and tumor fasta
             allinone   a wrapper for short reads simulation

### 2.1 Module vcf2fa 

The module `vcf2fa` can be used to build the genome sequences of normal cells. 
It takes a reference genome in fasta format (e.g. hg19.fa) and a list of SNPs 
in VCF format as input (only SNPs are acceptable). 
For the fasta file, only the chromosomes which will be simulated later 
should be included. For example, if you want to simulate short reads for human 
chromosome 1-5, all other chromosomes should not appear in the fasta file. 
For the VCF file, there should be only one sample's genotype information in it 
and the genotypes should be phased. 

After the running of `vcf2fa`, two fasta files will be generated. Each contains 
the sequences of all chromosomes of that haplotype, and the variants loci are 
replaced by according alleles.
By default, each chromosome will be treated as autosome and has one copy in 
each haplotype's fasta. User should explicitly specify sex chromosomes by 
`--sex_chr` if needed. For example, for simulating a male's genome, user should 
use `--sex_chr X,Y`. And for simulating a female's genome, user should use 
`--sex_chr X,X` or without any setting for `--sex_chr`. (Note: X,Y should be 
exactly the same string as in the fasta file.) 

### 2.2 Module phylovar

This module is the core module of CSiTE. It jointly simulates SNVs and CNVs of 
a sample of tumor cells along the history of a genealogy.
To run the module, a coalescent tree which captures the ancestral relationships
of the tumor cells in a sample is required. The input tree file should be in the
newick format, which can be generated by the [ms
program](http://home.uchicago.edu/rhudson1/source/mksamples.html) with the `-T`
option. ms program has the full-assemblage of the options needed to generate
complex demographic histories of the sample. 

With `--name`, users can specify the name of the sequence to simulate. 
And with `--length`, users can set the length of the sequence to simulate. 
The two most important parameters in the simulation are the mutation rates of 
SNVs (`--snv_rate`) and CNVs (`--cnv_rate`). These two parameters specify the 
rates of the whole sequence, no matter how long it is. The mutational events 
in CSiTE are simulated according to a Poisson process with parameters specified 
by user (see [Notes](https://github.com/hchyang/CSiTE#notes) for extra 
discussions).  

Other than the rate of SNVs, `--tstv` also affects the simulation of SNVs. This
parameter specifies the ratio of transition/transversion of SNVs. For
`phylovar`, the variants are in conceptual. So we use 0/1/2 to represent SNVs.
0 stands for transition and 1/2 stand for transversions. Those SNVs will be 
translated to concrete nucleotides by the module `chain2fa`. The map of the 
translation is:

reference | 0 | 1 | 2
----------|---|---|---
 N | N | N | N
 A | G | C | T
 G | A | C | T
 C | T | A | G
 T | C | A | G

Except for the rate of CNVs, there are five other parameters guiding CNVs
simulation. 

`--cnv_length_beta` can be used to specify the beta/mean parameter of the
exponential distribution (We assume the length of CNVs follows an exponential
distribution. This parameter specifies the mean length of the CNV events).

`--cnv_length_max` can be used to set the upper limit of CNVs' length (This will
effectively truncate the exponential distribution at this limit). 

`--del_prob` A CNV event can be a deletion or an amplification. Using this
parameter, we can specify the probability that a CNV event is a deletion.

`--copy_max` can be used to set the upper bound of an amplification event.  When
an amplification event happens, it will randomly pick a copy number according to
a geometric like distribution with Pr(n+1)=p\*Pr(n). The p parameter can be
specified by `-c/--copy_parameter`. The overall distribution will be normalized,
s.t. the total probability is 1. 

Chromosome level aneuploidy is widespread in tumor cells. In order to simulate 
this phenomenon of tumor genome, `phylovar` supplies `--parental` parameter. The 
default value of this parameter is '01', which tells `phylovar` to simulate two 
copies of the chromosome, one copy from each of the normal haplotypes. With the 
value of '001', `phylovar` can simulate a tumor chromosome with 3 haplotypes, 
and two of them are based on the haploype '0' of the normal chromosome, and the 
third one are based on the haplotype '1' of the normal chromosome.  

All of the previous parameters can be set in a config file (`--config`). This 
is quite convenient if you want to simulate variants in a genome, which 
contains multiple chromosomes. The format of the config file will be described 
in the following section.

The coalescent tree only describes the polymorphic part of the sample. You can
also simulate the truncal part of the tumor evolution (tumorigenesis) through
two different approaches: a) specify the trunk length e.g. using
`--trunk_length`. `--trunk_length` accepts a single float number that can be
used as the branch length leading to the root of the coalescent tree.  b)
specify the events on the trunk explicitly in a file through `--trunk_vars
filename`. The format of the trunk variants file is described in the later
section.

When the number of simulated cells is very large, the computational load can be
very heavy. Given the fact that most of the mutational events are extremely
rare, we implemented a pruning algorithm to trim the tree using `--prune` and
`--prune_proportion` options. For example, if you want to simulate the somatic
variants for a population containing 1,000,000 cells, after you setting `--prune
10000` or `--prune_proportion 0.01`, all subtrees with <=10,000 tips will be
trimmed into a tip node, which means there will be no polymorphic variants on
those subtrees with <=10,000 tips. So the tips belonging to the same subtree 
(with <=10,000 tips) will show the same genotypes.

The `phylovar` module of CSiTE can be used alone to simulte the conceptual somatic 
variants, or be used with other modules to simulate the concrete short reads of 
the tumor samples. To accomplish the first mode, `phylovar` accepts the 
`-D/--depth` and `-P/--purity` parameters. With these setting, `phylovar` simulates 
the read depth on each somatic variant sites and the frequency of the mutant 
allele. Those information will be saved in SNVs file, which will be described in 
the later section. 

###### Notes
By default, the branch length from the ms program is measured in 4N
generations. `phylovar` will simulate SNVs/CNVs events following a Poisson 
process. For example, if we set -r to 100 in `phylovar` (CSiTE), this is 
equivalent to set the population rescaled parameter -t to 100 in ms (see ms 
manual for details). 

#### 2.2.1 Input files

##### Tree file (-t/--tree)

Tree file should be in newick format, for example:

    ((2:0.083,4:0.083):0.345,(5:0.322,(1:0.030,3:0.030):0.292):0.105);

##### Trunk variants file (--trunk\_vars)

You can simulate truncal variants by specify `--trunk_length`, or import trunk
variants directly using `--trunk_vars` option.  The file format of truncal
variants is like:
    
    #chr hap start end var
    1 0 364645 364646 0
    1 1 464646 466646 +3
    2 0 465636 465637 1 
    2 1 468637 472633 -1

- **chr**:       the chromosome of the variant
- **hap**:       which haplotype the variant locates in
- **start**:     the start position of the variant (0-based, inclusive)
- **end**:       the end position of the variant (0-based, exclusive)
- **var**:       the type of the variant. 0/1/2: SNV, -1: deletion, 
+int: amplification

P.S. start and end are 0 based. And the region of each variant is similar to the
bed file (inclusive on the left, but exclusive on the right end,
i.e.\[start,end)).

##### Configure file (--config)

The command line paramaters are only suitable for simulating a single sequence.
When simulate multiple sequences, e.g. all chromosomes in a genome, user can 
supply a configure file. A typical contains two sections, genome section and 
chromosomes section. In the genome section, user can specify genome-wide 
variants settings. And in chromsomes section, user can customise the variants 
settings for each chromosome. It's very flexible. Different chromosomes can have 
different lengths, different SNV/CNV muations rates, and even different haplotype 
composition. 

The configure file should be specified in YAML format. Here is an example of the 
configure file:

    genome:
        snv_rate: 2.0
        cnv_rate: 0.6
        del_prob: 0.5
        cnv_length_beta: 200
        cnv_length_max: 400
        copy_parameter: 0.5
        copy_max: 5
        parental: '01'
        tstv: 2.0
        length: 109000
    chromosomes:
        - '1':
            length: 100000
            parental: '00'
        - '2':
            snv_rate: 2.0
            length: 9000

All of the parameters under genome section have to be specified in the configure 
file. The snv\_rate/cnv\_rate/length in the genome section, is the sum of those 
values of all chromosomes. The parental and the name of each chromosome should be 
quoted to make them in the type of string.

#### 2.2.2 Output files

##### SNVs file (-S/--snv)

This file contains the frequency information of simulated SNVs. There are six 
columns in this file. 
- **chr**:         the chromosome of SNVs
- **pos**:         the position (0-based) of SNVs
- **form**:        the form of SNVs (0/1/2)
- **true\_freq**:  the true frequency of alternative allele across the cell 
population
- **sim\_dp**:     the simulated total coverage of tumor and normal cells in 
the tumor sample at the position of SNV
- **sim\_freq**:   the observed frequency of alternative allele across the cell 
population

P.S. Do not mix the information of simulated depth and frequency here up with 
the coverage and frequency of simulated short reads. They are two independent 
processes.

##### CNVs file (-V/--cnv)

This file contains the information of simulated CNVs. There are five columns in
this file. 
- **chr**:         the chromosome of CNVs
- **start**:       the start position of CNVs (0-based, inclusive)
- **end**:         the end position of CNVs (0-based, exclusive)
- **copy**:        the copy changes of CNVs
- **carrier**:     the number of tumor cells in the sample carring the CNVs

##### CNV profile file (--cnv\_profile)

This file contains the CNV profile across each chromosome. There are four 
columns in this file. 
- **chr**:         the chromosome of regions
- **start**:       the start position of regions (0-based, inclusive)
- **end**:         the end position of regions (0-based, exclusive)
- **local\_cp**:   the copy number of the local region

##### Nodes variants file (-N/--nodes\_vars)

`phylovar` can output the variants (SNVs/CNVs) occured on the branch leading 
to each node. There are six columns in this file:
- **node**:        the id of nodes (we name each node in the format 'nodeXXX', 
in which XXX is an integer starting from 1)
- **chr**:         the chromosome of the variant
- **hap**:         the haplotype the variant locates in
- **start**:       the start position of the variant (0-based, inclusive)
- **end**:         the end position of the variant (0-based, exclusive)
- **var**:         the type of the variant. 0/1/2: SNV, -1: deletion, 
+int: amplification

##### Variants tree file (-T/--vars\_tree)

`phylovar` can output a [NHX
file](https://sites.google.com/site/cmzmasek/home/software/forester/nhx) with
each node's id and all variants attached. The variants are encoded in the form 
of 'chr#hap#start#end#var'. As the 'var' cloumn in nodes variants file, 0/1/2 stand 
for SNV, -1 stands for deletion and +int stands for amplification.

##### SNV genotype file (--snv\_genotype)

This file contains the snv\_genotype of each tumor cell at SNV loci. Each SNV
has one record. The first two column are the coordinate of the SNV. 
Subsequently, there is one column for each tumor cell, which encodes the 
genotype of each tumor cells. The snv genotype is in the form of ‘M:N’. M 
denotes the number of alternative allele and N denotes the number of reference 
allele.
- **chr**:         the chromosome of SNVs
- **pos**:         the position (0-based) of SNVs
- **form**:        the form of SNVs (0/1/2)
- **cell1**:       the genotype of cell 1
- **cell2**:       the genotype of cell 2
- **...**:         ...

##### Individual CNVs file (--ind\_cnvs)

This file contains the CNVs on each parental copy of each single cell in the
sample. There are six columns in this file:
    
    #cell parental  start     end       copy
    1       0       7912422   7930111   2
    1       1       43110140  43341629  1
    2       0       2255734   2299608   -1
    2       0       22660687  22788472  -1
    2       1       59756841  61142076  3

- **cell**:        the id of cells in the sample
- **chr**:         the chromosome of SNVs
- **parental**:    which parental copy the variant locates in (0 means one 
of the parental copy, 1 means another copy, 2 means the third copy if your 
sample is triploid...)
- **start**:       the start position of the CNV
- **end**:         the end position of the CNV
- **copy**:        an integer. -1: deletion; positive integer: amplification

P.S. start and end are 0 based. And the region of each variant is similar to the
bed file (inclusive on the left, but exclusive on the right end,
i.e.\[start,end)).

##### parental copy file (--parental\_copy)

This file contains the information of parental copy for each SNV. The first two
columns are the coordinate of the SNV, and followed by N columns if the ploidy 
is N.

##### chain file and nodes profile file (--chain)

There will be two kind of files generated under the folder specified by this
parameter. One is chain file and the other is nodes profile file.

There is a chain file for each tip node. This file contains the information of 
perturbed genome of simulated tumor cell. Each haplotype of each chromosome of 
the tumor genome will have a block in the chain file. Each block has a header 
section and a body section. For example:

    >1_Haplotype0 parental:0
    1       0       33      ref
    1       33      55      -1
    1       55      99      ref
    1       55      100000  ref

The first line, which starts with '>' is the header line. There are two fields 
in the header line. The first is the sequence name and the second is the parental 
copy of that sequence. And then there four columns in the body section: 
- **chr**:         the chromosome of SNVs
- **start**:       the start position of the variant (0-based, inclusive)
- **end**:         the end position of the variant (0-based, exclusive)
- **var**:         the type of the variant. 0/1/2: SNV, -1: deletion, 'ref': 
reference

P.S. There is no +int in the fourth column of body section to indicate a CNV event.
Instead, the amplified regions will have multiple record to represent that event.
e.g. the last two lines in the example, indicate there a amplification in region 
1:55-59 in that sequence.

There are three columns in a nodes profile file:
- **node**:        the id of tip nodes (after pruning if the tree is pruned)
- **l\_count**:    the count of original leaves under this node 
- **l\_name**:     the names of original leaves under this node 

##### Log file (-g/--log)

This file contains logging information, e.g. the command line parameters and the
random seed used. You can use these information to replicate the simulation.
After setting the `--loglevel` as DEBUG, `phylovar` will output detailed information
of the simulation process. 

##### all options of module phylovar

Will be filled later.

### 2.3 Module chain2fa

After generating genomes of normal cells with `vcf2fa` and chain files of tumor
cells with `phylovar`, we can build the genome sequences for each tumor cell. 
This work is done by the module `chain2fa`. And the genomes of tumor cells and 
normal cells can then be used by ART to simulate short reads. For the 
`--reference` parameter, all reference fasta files should list sequentially and
seperated by comma, like `--reference normal_hap0.fa,normal_hap1.fa`.  

### 2.4 Module fa2ngs

After the running of previous three modules, we will get the genome sequences of 
normal/tumor cells generated. By supplying the purity and depth settings we want 
to simulate, module `fa2ngs` will figure out the fold of covarge for each genome 
to simulate and employ ART to simulate the short reads for the genomes in the 
tumor sample. `fa2ngs` uses an very flexible way to cooperate with ART. You can
pass all parameters to ART by `--art` except the fcov, in, out, id, and rndSeed
parameters of ART, as those parameters are handled by `ngs2fa` itself. For 
example, you can use `--art '/path/to/ART/art_illumina --noALN --quiet --paired 
--len 100 --mflen 500 --sdev 20'` if ART are not installed system-wild. Or you 
can use `--art 'echo art_illumina --noALN --quiet --paired --len 100 --mflen 
500 --sdev 20'` if you just want the `ngs2fa` print out the art command it used.

This module will also generate a meta file for Wessim, which is a targeted 
re-sequencing simulator that generates synthetic exome sequencing reads from a 
given sample genome. With this file and probe information, Wessim can simulate 
the exome-sequencing data of the simulated tumor sample.

### 2.5 Module allinone

This is a convenient wrapper for simulating reads for tumor sample. It calls the
module `vcf2fa` to build the genome of normal cells, and calls the modules 
`phylovar` and `chain2fa` to build the genomes of tumor cells. And then according 
the purity and coverage settings customised by user, `fa2ngs` are called to 
simulate short reads for the tumor sample. Besides the short reads for tumor 
sample, `allinone` will also simulate the short reads for nomal sample.

Most of the parameters of this module will be passed to specific modules, except 
the parameters `-o/--output`, `-g/--log` and `--start`. `--output` takes a string
as the ouput folder name, and all output of the simulation will be put in this 
folder. `allinone` will create 4 subfolders under that folder. These subfolders 
are `normal_fa`, `tumor_chain`, `tumor_fa` and `art_reads`, and their content 
should be self-explantory. The `--log` parameter specify the name of log file, 
which will records the random seed and all of the sub-commands and their 
parameters called by `allinone`. Sometimes, `allinone` will fail in the middle 
if user specified an inappropriate parameter. In this case, user can use the 
`--start` parameter to skip some steps which are finished appropriately. Only 
1-4 are acceptable:

- **1**:     vcf2fa
- **2**:     phylovar
- **3**:     chain2fa
- **4**:     fa2ngs

## 3. Examples

* Simulate the coalescent tree of a sample of 1000 tumor cells, which are
sampled from a exponetially growing tumor. 
(consult the manual of ms for more information)

    `ms 1000 1 -T -G 1 |tail -n1 > ms_tree.txt`

* Simulate the somatic SNVs of this sample. We assume the sequencing depth is
60, and the purity of the sample is 0.8, which means there are 250 normal
cells other than these 1000 tumor cells. Other settings are: 
a) the mutation rate of SNVs and CNVs are 10 and 0.1 respectively; 
b) the sequence length is 135534747 (chr10 of hg19); 
c) the cells of the sample are diploid. We save the
frequncy of the simulated SNVs to file 'snvs\_freq.txt'.

    `csite.py -t ms_tree.txt -P 0.8 --length 135534747 -r 10 -R 0.1 -D 60 -S
    snvs_freq.txt`

* There are no truncal muations in the simulation above. If you want to simulate
truncal muations, use the option `--trunk_length`.

    `csite.py -t ms_tree.txt -P 0.8 --length 135534747 -r 10 -R 0.1 -D 60 -S
    snvs_freq.txt --trunk_length 2.0`

* If you want to ignore the variants with the frequency <=0.01, you can use
`--prune 20` or `--prune_proportion 0.02` (we use `--prune 20` instead of
`--prune 10` for the cells are diploid in the our simulation). These two
options can be used to accelerate the simulation when your tree is huge.

    `csite.py -t ms_tree.txt -P 0.8 --length 135534747 -r 10 -R 0.1 -D 60 -S
    snvs_freq.txt --trunk_length 2.0 --prune 20`

    or

    `csite.py -t ms_tree.txt -P 0.8 --length 135534747 -r 10 -R 0.1 -D 60 -S
    snvs_freq.txt --trunk_length 2.0 --prune_proportion 0.02`

* If you want to save the SNVs genotypes for each single cell for exactly the
same simulation as above, use the options `--snv_genotype` and
`--random_seed`.

    `csite.py -t ms_tree.txt -P 0.8 --length 135534747 -r 10 -R 0.1 -D 60 -S
    snvs_freq.txt --trunk_length 2.0 --prune_proportion 0.02 --snv_genotype
    snvs_genotype.txt --random_seed xxxx`

    P.S. The random seed xxxx can be found in the log file of the previous
    simulation.

## Authors

* [Hechuan Yang](https://github.com/hchyang)

## License

This project is licensed under the GNU GPLv3 License - see the
[LICENSE](LICENSE) file for details.

