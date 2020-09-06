**PHASED DIPLOID GENOME ASEMBLY PIPELINE**
This pipeline is designed to take Long-read and illumina sequencing data of a heterozygous diploid and provide a phased diploid assembly <br/>
This has been tested on *S.cerevisiae* with nanopore data
It appears to work well with strains containing at least 0.2% heterozygosity, an average long-read length of 8kb and >60X coverage. <br/>
**_NOTE: Example commands below use S.cerevisiae and nanopore specifications if required by the tool_** <br/>
The test data (strain SO002) comes from a cross between stable haploids YLF161 and YLF132 (https://doi.org/10.1002/yea.2991) from West African (DVBPG6044) and North American (YPS128) backgrounds respectively <br/>
https://drive.google.com/drive/folders/1IkuQMK5FLbndHcrX8vhoCPzKWePJleAl?usp=sharing <br/>
Provided are nanopore reads (long_reads.fq.gz), basecalled by guppy (v3.4.5) and with adapter and barcodes removed using porechop <br/>
176,003 reads <br/>
106X coverage <br/>
4,402 median length <br/>
7,277 mean length <br/>
12,519 N50 <br/>
Both parents have reference quality assemblies for comparison (https://doi.org/10.1038/ng.3847)



The essential pipeline has 10 steps: <br/>
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
Illumina data = illumina_1.fq.gz and illumina_2.fq.gz

**1. *De-novo* genome assembly using long-reads** <br/>
The assembler of choice is at to your discretion. <br/>
For this example I will use canu (best results with latest, at the time, version v2) at is has proved to work very well in the case of *S.cerevisiae* <br/>
    
    ##downsample entire reads to 40x using filtlong, selecting for the longest reads
    filtlong long_reads.fq.gz --length_weight 10 --min_length 1000 -t 480000000 | gzip > ./long_reads_fl.fq.gz
    ## Canu assembly
    # predicted assembly size is based on the haploid genome size of S.cerevisiae
    # perhaps you need to add useGrid=false unless running on a configured cluster
    canu -p SO002 -d SO002_haploid_canu_dir \
    genomeSize=12m \
    -nanopore long_reads_fl.fq.gz
    # assembly is in SO002_haploid_canu_dir/SO002_haploid.contigs.fasta
    mv SO002_haploid_canu_dir/SO002_haploid.contigs.fasta ./SO002_haploid.canu.fa

**2. Genome polishing using long and short reads** <br/>
This is only a recommended polishing scheme <br/>
After benchmarking different polishing schemes with a variety of assemblers on *S. cerevisiae*, those with read correction steps benefitted most from 1xRacon-2xMedaka-3xPilon <br/>
Assemblers without read-correction benefit from two more additional rounds of Racon <br/>

    ## 1 round of racon poloshing
    # minimap mapping first
    minimap2 -ax map-ont SO002_haploid.canu.fa long_reads.fq.gz > SO002_haploid.sam
    # feed sam file into racon
    racon long_reads.fq.gz SO002_haploid.sam SO002_haploid.canu.fa > SO002_haploid.canu_r1.fa
    # clean up, remove sam file
    rm SO002_haploid.sam
    
    ## 2 rounds of Medaka Polishing
    # may need to change model (-m) depending on organism
    medaka_consensus -i long_reads.fq.gz -d SO002_haploid.canu_r1.fa -o SO002_haploid.canu_r1_m1 -v -m r941_min_high_g351
    # get polished genome
    mv SO002_haploid.canu_r1_m1/consensus.fasta ./SO002_haploid.canu_r1_m1.fa
    # clean up, remove other medaka output
    rm -r SO002_haploid.canu_r1_m1/
    # repeat for second round of polishing using new polished genome as input
    medaka_consensus -i long_reads.fq.gz -d SO002_haploid.canu_r1_m1.fa -o SO002_haploid.canu_r1_m2 -v -m r941_min_high_g351
    mv SO002_haploid.canu_r1_m2/consensus.fasta ./SO002_haploid.canu_r1_m2.fa
    rm -r SO002_haploid.canu_r1_m2/
    # clean up, remove index files
    rm *.mmi
    rm *.fai
    
    ## 3 rounds of Pilon
    # indexing of genome to polish   
    bwa index SO002_haploid.canu_r1_m2.fa
    picard CreateSequenceDictionary R= SO002_haploid.canu_r1_m2.fa O= SO002_haploid.canu_r1_m2.dict
	samtools faidx SO002_haploid.canu_r1_m2.fa
    # mapping of illumina to long-read polished genome to generate a bam file
    bwa mem SO002_haploid.canu_r1_m2.fa illumina_1.fq.gz illumina_2.fq.gz | samtools sort -o SO002_haploid.bwamem.sorted.bam -
    # add some read groups
    gatk AddOrReplaceReadGroups \
    -I SO002_haploid.bwamem.sorted.bam \
    -O SO002_haploid.bwamem.sorted_RG.bam \
    --RGID SO002 --RGPL ILLUMINA --RGSM SO002 --RGLB 1 --RGPU 1
    # remove duplicates using GATK4
    gatk MarkDuplicatesSpark \
 	-I SO002_haploid.bwamem.sorted_RG.bam \
 	-O SO002_haploid.bwamem.sorted_RG_markdup.bam \
 	--create-output-bam-index \
    --remove-all-duplicates TRUE
    # clean up,remove older bams
    rm SO002_haploid.bwamem.sorted.bam
    rm SO002_haploid.bwamem.sorted_RG.bam
    # provide bam and long-read polished assembly to pilon
    pilon \
    --genome SO002_haploid.canu_r1_m2.fa \
    --bam SO002_haploid.bwamem.sorted_RG_markdup.bam  \
    --fix snps,indels --vcf --changes --output SO002_haploid.pilon 
    # clean up, remove pilon output, genome indexes, previous bam files and change name of polished genome file
    rm SO002_haploid.pilon.changes
    rm SO002_haploid.pilon.vcf
    rm SO002_haploid.canu_r1_m2.fa.*
    rm SO002_haploid.canu_r1_m2.dict
    rm SO002_haploid.bwamem.sorted_RG_markdup.bam
    mv SO002_haploid.pilon.fasta SO002_haploid.canu_r1_m2_p1.fa
    # repeat polishing steps 2 more times
    bwa index SO002_haploid.canu_r1_m2_p1.fa
    picard CreateSequenceDictionary R= SO002_haploid.canu_r1_m2_p1.fa O= SO002_haploid.canu_r1_m2_p1.dict
	samtools faidx SO002_haploid.canu_r1_m2_p1.fa
    bwa mem SO002_haploid.canu_r1_m2_p1.fa illumina_1.fq.gz illumina_2.fq.gz | samtools sort -o SO002_haploid.bwamem.sorted.bam -
    gatk AddOrReplaceReadGroups \
    -I SO002_haploid.bwamem.sorted.bam \
    -O SO002_haploid.bwamem.sorted_RG.bam \
    --RGID SO002 --RGPL ILLUMINA --RGSM SO002 --RGLB 1 --RGPU 1
    gatk MarkDuplicatesSpark \
 	-I SO002_haploid.bwamem.sorted_RG.bam \
 	-O SO002_haploid.bwamem.sorted_RG_markdup.bam \
 	--create-output-bam-index \
    --remove-all-duplicates TRUE
    rm SO002_haploid.bwamem.sorted.bam
    rm SO002_haploid.bwamem.sorted_RG.bam
    #provide bam and long-read polished assembly to pilon
    pilon \
    --genome SO002_haploid.canu_r1_m2_p1.fa \
    --bam SO002_haploid.bwamem.sorted_RG_markdup.bam  \
    --fix snps,indels --vcf --changes --output SO002_haploid.pilon 
    rm SO002_haploid.pilon.changes
    rm SO002_haploid.pilon.vcf
    rm SO002_haploid.canu_r1_m2_p1.fa.*
    rm SO002_haploid.canu_r1_m2_p1.dict
    rm SO002_haploid.bwamem.sorted_RG_markdup.bam
    mv SO002_haploid.pilon.fasta SO002_haploid.canu_r1_m2_p2.fa
    # third and final round
    bwa index SO002_haploid.canu_r1_m2_p2.fa
    picard CreateSequenceDictionary R= SO002_haploid.canu_r1_m2_p2.fa O= SO002_haploid.canu_r1_m2_p2.dict
	samtools faidx SO002_haploid.canu_r1_m2_p2.fa
    bwa mem SO002_haploid.canu_r1_m2_p2.fa illumina_1.fq.gz illumina_2.fq.gz | samtools sort -o SO002_haploid.bwamem.sorted.bam -
    gatk AddOrReplaceReadGroups \
    -I SO002_haploid.bwamem.sorted.bam \
    -O SO002_haploid.bwamem.sorted_RG.bam \
    --RGID SO002 --RGPL ILLUMINA --RGSM SO002 --RGLB 1 --RGPU 1
    gatk MarkDuplicatesSpark \
 	-I SO002_haploid.bwamem.sorted_RG.bam \
 	-O SO002_haploid.bwamem.sorted_RG_markdup.bam \
 	--create-output-bam-index \
    --remove-all-duplicates TRUE
    rm SO002_haploid.bwamem.sorted.bam
    rm SO002_haploid.bwamem.sorted_RG.bam
    pilon \
    --genome SO002_haploid.canu_r1_m2_p2.fa \
    --bam SO002_haploid.bwamem.sorted_RG_markdup.bam  \
    --fix snps,indels --vcf --changes --output SO002_haploid.pilon 
    rm SO002_haploid.pilon.changes
    rm SO002_haploid.pilon.vcf
    rm SO002_haploid.canu_r1_m2_p2.fa.*
    rm SO002_haploid.canu_r1_m2_p2.dict
    rm SO002_haploid.bwamem.sorted_RG_markdup.bam
    # viola, the polished assembly
    mv SO002_haploid.pilon.fasta SO002_haploid.canu_r1_m2_p3.fa

**OPTIONAL STEP: Scaffolding and negative gap closing** <br/>
Scaffolding assemblies using a reference or optical data can help with phasing <br/>
A simple scaffolding procedure uses ragout, RaGOO or RagTag and a reference genome <br/>

    ## RagTag example, increasing minimum unique alignment length (-f) to reduce overlaps
    ragtag.py -f 5000 yeast.tidy.fa SO002_haploid.canu_r1_m2_p3.fa
    mv ragtag_output/ragtag.scaffolds.fasta SO002_haploid.canu_polished_scaffolded.fa
    #clean up, remove ragtag output
    rm -r ragtag_output/
    ##check scaffolding against reference using dotplot in order to see scaffolding job
    #make sure that all the reference is covered in the renamed and orientated contigs as we want to remove haplotigs and reduce assembly to complete haploid merged assembly. Any lost haplotig variation is recovered after re-assembling.
    

Additionally it has been found that some assemblies can fail to connect overlapping contigs, called negative-gaps <br/>
Tools such as gap5 and TGS-gapcloser can be used in order to find end-to-end overlaps in contigs in order to join contigs and create a consensus at the overlap position <br/>
An additional round of pilon can help with correcting any potential errors introduced during the generation of a consensus

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
    
