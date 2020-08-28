This pipeline is designed to take Long-read and illumina sequencing data of a heterozygous diploid and provide a phased diploid assembly <br/>
This has been tested on *S.cerevisiae* with nanopore data and appears to work well with strains containing a at least 0.2% heterozygosity, average long-read length of 8kb and >60X coverage. <br/>
**_NOTE: Example commands below use S.cerevisiae and nanopore specifications if required by the tool_** <br/>

The essential pipeline has N steps
1. *De-novo* genome assembly using long-reads <br/>
2. Genome polishing using long (Racon and Medaka) and short (Pilon) reads <br/>
OPTIONAL STEP: Scaffolding and negative-gap closing <br/>
3. Alignment of long-reads to assembly to generate bam file <br/>
4. Variant calling using illumina reads to generate a vcf containing heterozygous variants <br/>
5. Input of bam and vcf into whatshap to phase variants <br/>
6. Input of bam and phased-vcf into whathap to phase reads <br/>
7. Assessing phasing results through custom R script and statistics such as v90 <br/>
8. Arbitrarily combining blocks and unphased reads to generate haplotypes read sets <br/>
9. *De-novo* assembly and long-read polishing of each haplotype seperately <br/>
10. Illumina read polishing of diploid genome <br/>

long-read data = long_reads.fq.gz (do not have test set currently MAYBE CAN LOOK AT USING THE HYBRID MADE IN THE LAB) <br/>
Illumina data = illumina.1.fq.gz and illumina.2.fq.gz

**1. *De-novo* genome assembly using long-reads** <br/>
The assembler of choice is at to your discretion. <br/>
For this example I will use canu (best results with latest, at the time, version v2) at is has proved to work very well in the case of *S.cerevisiae* <br/>

    ## Canu assembly
    # predicted assembly size is based on the haploid genome size of S.cerevisiae
    COMMANDS

**2. Genome polishing using long and short reads** <br/>
This is only a recommended polishing scheme <br/>
After benchmarking different polishing schemes with a variety of assemblers on *S. cerevisiae*, those with read correction steps benefitted most from 1xRacon-2xMedaka-3xPilon <br/>
Assemblers without read-correction benefit from two more additional rounds of Racon <br/>

    ## 1 round of racon poloshing
    COMMANDS
    ## 2 rounds of Medaka Polishing
    COMMANDS
    ## 3 rounds of Pilon
    COMMANDS

**OPTIONAL STEP: Scaffolding and gap closing** <br/>
Scaffolding assemblies using a reference or optical data can help with phasing <br/>
A simple scaffolding procedure uses ragout or RaGOO and a reference genome <br/>

    ## RaGOO example
    COMMANDS
    WATCH OUT FOR OVER ASSEMBLING AND CONTIG 0 GENERATION

Additionally it has been found that some assemblies can fail to connect overlapping contigs, called negative-gaps <br/>
Tools such as gap5 and TGS-gapcloser can be used it order to find end-to-end overlaps in contigs in order to join contigs and create a consensus at the overlap position <br/>
An additional round of pilon can help with correcting any errors introduced during the generation of a consensus

**3. Alignment of long-reads to assembly to generate bam file** <br/>

    ## Minimap alignment
    COMMANDS

**4. Variant calling using illumina reads to generate a vcf containing heterozygous variants** <br/>

    ## GATK variant calling
    COMMANDS
    
**5. Input of bam and vcf into whatshap to phase variants** <br/>
Whatshap will take long read mapping and variants in order to produce a vcf containing variant haplotype and block information <br/>
This can be seen in the last vcf column; the first value 1|0 meaning haplotype 1 (and vice-versa) and the last value XXXX denoting the block number.

    ## Whatshap variant phasing
    COMMANDS
     
**6. Input of bam and phased-vcf into whatshap to phase reads** <br/>
Whatshap only uses a subset of the most informative reads in order to phase <br/>
Therefore it is necessary to run another script in order to phase the entire long-read dataset

    ## Whatshap read phasing
    COMMANDS
    ## Bash commands to generate fastq.gz files for each block haplotype
    COMMANDS

**7. Assessing phasing results through custom R script and statistics such as v90** <br/>
To begin a few simple bash commands can provide a summary file of phasing statistics <br/>
However, these usual phasing statitics, such as the number of blocks, do not provide an intuitive image of phasing quality <br/>
Here I introduce a new statistic V90 which is the minimum number of blocks to cover 90% of all phased variants (similar to statistic such as N90) <br/>
Simply the lower V90 the more contiguous the phasing. However, in respect to the assembly, this does not tell you if two large blocks exist within a single large contig <br/>
The normV90 goes one step further and divides the number of blocks in the V90 by the number of contigs involved within these respective blocks <br/>
Therefore a normV90 value of 1 signifies that for all contigs containing 90% of phased variants there is only one significant contiguous block

    ## Bash script commands to get statistics
    COMMANDS
    
Additionally another route is to visually inspect how well your genome is phased <br/>
Whatshap contains an IGV route for visualisation here: https://whatshap.readthedocs.io/en/latest/guide.html#visualizing-phasing-results
Additionally, below is an R script which will visualise the position of heterozygous variants within the *de-novo* assembly and the phased blocks

    ## R script needs to become automatised
    COMMANDS
    
    
**8. Arbitrarily combining blocks and unphased reads to generate haplotypes read sets** <br/>
On the assumption that the phasing is contiguous and that normV90 is close to 1, arbitrarily combining reads from blocks denoted haplotype 1 will combine a minimal number of read-blocks **within** a contig. Additionally one read-block will contain the vast majority of all variants within this contig. <br/>
However due to LOH, only taking phased variants would eliminate reads within homozygous regions <br/>
Therefore, assuming that a significant number of reads relative to the total number have been phased, the unphased reads can be added to both haplotype read sets in order to generate more complete haplotype assemblies <br/>

    ## Bash commands to join haplotypes 
    COMMANDS
    
    
**9. *De-novo* assembly and long-read polishing of each haplotype seperately** <br/>
As in section 1 and 2 each haplotype read set must be assembled individually and then polished by their respective reads

    ## Canu Assembly
    COMMANDS
    ## 1 round of racon poloshing
    COMMANDS
    ## 2 rounds of Medaka Polishing
    COMMANDS

**10. Illumina read polishing of diploid genome** <br/>
Illumina reads have not been phased and therefore it appears better to feed Pilon the combined haplotype assemblies for polishing <br/>

    ## Join haplotype assemblies into diploid assembly (need to rename contigs)
    COMMANDS
    ## 3 rounds of Pilon polishing
    COMMANDS
    
