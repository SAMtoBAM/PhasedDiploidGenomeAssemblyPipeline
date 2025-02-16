##downsample entire reads to 40x using filtlong, selecting for the longest reads and removing any reads shorter than 1kb
filtlong long_reads.fq.gz --length_weight 10 --min_length 1000 -t 480000000 | gzip > ./long_reads_fl.fq.gz
## Canu assembly
# predicted assembly size is based on the haploid genome size of S.cerevisiae
canu -p SO002 \
	-d SO002_haploid_canu_dir \
	genomeSize=12m \
	useGrid=false \
	-nanopore long_reads_fl.fq.gz
# assembly is in SO002_haploid_canu_dir/SO002_haploid.contigs.fasta
mv SO002_haploid_canu_dir/SO002_haploid.contigs.fasta ./SO002_haploid.canu.fa


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
rm SO002_haploid.bwamem.sorted_RG_markdup.ba*
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
rm SO002_haploid.bwamem.sorted_RG_markdup.ba*
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
rm SO002_haploid.bwamem.sorted_RG_markdup.ba*
# et viola, the polished assembly
mv SO002_haploid.pilon.fasta SO002_haploid.canu_r1_m2_p3.fa

##ragout example
#create recipe file with names of reference prefix, new assembly prefix, reference assembly, new assembly, prefix of genome for scaffold names
echo "./reference = ref" > ragout.recipe
echo "target = SO002_haploid" >> ragout.recipe
echo "ref.fasta = S288c.masked.fa" >> ragout.recipe
echo "SO002_haploid.fasta = SO002_haploid.canu_r1_m2_p3.fa" >> ragout.recipe
echo ".naming_ref = ref" >> ragout.recipe
#run ragout
ragout -o SO002_haploid_ragout --solid-scaffolds ragout.recipe
cat  SO002_haploid_ragout/SO002_haploid_scaffolds.fasta SO002_haploid_ragout/SO002_haploid_unplaced.fasta > SO002_haploid.canu_polished_scaffolded.fa
# clean up, remove ragout output file
rm -r SO002_haploid_ragout/
## check scaffolding against reference using dotplot (an easy sway to do this is using a tool called D-GENIES) in order to see scaffolding results
# make sure that all the reference is covered in the renamed and orientated contigs as we want to remove haplotigs and reduce assembly to complete haploid merged assembly
## select only the scaffolds/contigs associated with the reference by name (this assembly is within the google drive folder)
cat SO002_haploid.canu_polished_scaffolded.fa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' | grep "\S" | sed 's/>chr_/>chr/g' | grep -A 1 'chr' > SO002_haploid.FINAL.fa


## Minimap alignment
minimap2 -ax map-ont SO002_haploid.FINAL.fa long_reads.fq.gz | samtools sort -o SO002_haploid.LR_FINAL.sorted.bam -
samtools index SO002_haploid.LR_FINAL.sorted.bam



## GATK variant calling
#As before for polishing..
bwa index SO002_haploid.FINAL.fa
picard CreateSequenceDictionary R= SO002_haploid.FINAL.fa O= SO002_haploid.FINAL.dict
samtools faidx SO002_haploid.FINAL.fa
#map illumina and can just put the read group information directly whilst mapping
bwa mem SO002_haploid.FINAL.fa \
-R '@RG\tID:SO002\tPL:ILLUMINA\tPI:0\tSM:SO002\tLB:1' \
illumina_1.fq.gz illumina_2.fq.gz | samtools sort -o SO002_haploid.SR_FINAL.sorted.bam -
gatk MarkDuplicatesSpark \
-I SO002_haploid.SR_FINAL.sorted.bam \
-O SO002_haploid.SR_FINAL.sorted.markdup.bam \
--create-output-bam-index
gatk HaplotypeCaller \
-R SO002_haploid.FINAL.fa \
-I SO002_haploid.SR_FINAL.sorted.markdup.bam \
-O SO002_haploid.SR_FINAL_raw_SNPs_indels.vcf
# clean up, no longer need bam files etc, can remove
rm SO002_haploid.SR_FINAL.sorted.bam
rm SO002_haploid.SR_FINAL.sorted.markdup.ba*
# can apply some crude filters to the vcf in order to improve the calls
# I don't apply anything sophisticated, just removing those with low coverage
# here I remove everything with less than 30x coverage and it appears to work well with most strains with high coverage
vcftools --vcf SO002_haploid.SR_FINAL_raw_SNPs_indels.vcf \
--min-meanDP 30 \
--recode --recode-INFO-all \
--out SO002_haploid.SR_FINAL_SNPs_indels
mv SO002_haploid.SR_FINAL_SNPs_indels.recode.vcf SO002_haploid.SR_FINAL.vcf
# we can check on the number of variants present
# the total number of variants
grep -e "1/1" -e "0/1" -e "1/2" SO002_haploid.SR_FINAL.vcf -c
# the total number of heterozygous variants
grep -e "0/1" -e "1/2" SO002_haploid.SR_FINAL.vcf -c
# there should be minimal difference between the two due to the polishing correcting other differences
# doing the same on the file before filtering will tell if the filtering has been potentially too conservative
# During my final test run, I received 69076 and 67787 respectively
# clean up
rm SO002_haploid.SR_FINAL_raw_SNPs_indels.vc*
rm SO002_haploid.FINAL.dict
rm SO002_haploid.SR_FINAL_SNPs_indels.log

## Whatshap variant phasing
# run the 'phase' script from whatshap to generate a phased vcf
# here I also have indels included, this appears to improve phasing substantially
# some tool output summary data is placed in the SO002_haploid.whatshap_output.txt file
whatshap phase --ignore-read-groups --indels \
--reference SO002_haploid.FINAL.fa \
-o SO002_haploid.phased.vcf  \
SO002_haploid.SR_FINAL.vcf SO002_haploid.LR_FINAL.sorted.bam
# check out phasing results
whatshap stats --gtf SO002_haploid.phased.gtf SO002_haploid.phased.vcf >> SO002_haploid.whatshap_output.txt
## in google drive is the output from my last run "SO002_haploid.whatshap_output.txt"
## in step 7 we will simplify this output to make it more directly interpretable

## Whatshap read phasing
# compress and index vcf
bgzip SO002_haploid.phased.vcf
tabix -p vcf SO002_haploid.phased.vcf.gz
# phase reads using haplotag
whatshap haplotag --ignore-read-groups \
-o SO002_haploid.phased.bam \
--reference SO002_haploid.FINAL.fa \
SO002_haploid.phased.vcf.gz SO002_haploid.LR_FINAL.sorted.bam
## Bash commands to generate fastq.gz files for each block haplotype, perhaps there are better ways to do all this
# what I essentially do is split all the reads into blocks and then haplotypes, incase it is important to have them seperated by block
# for example in this case we could imagine pairing block haplotypes by their original parents
# get all the phase block numbers
samtools view SO002_haploid.phased.bam | grep 'PS:i:' | awk -F ":" '{print $NF}' | sort | uniq > SO002_haploid.phased_PS.txt
# go block by block extracting reads into fastq files
mkdir SO002_PS_blocks/
cat SO002_haploid.phased_PS.txt | while read block 
do
	mkdir SO002_PS_blocks/$block
	# extract all alignments for each block
	samtools view SO002_haploid.phased.bam |\
	awk '{if($1 ~ /^@/) print $0; else if ($0 ~ /PS:i:'$line'/) print $0}' > SO002_PS_blocks/$block/${block}.sam
	# extract only HP1 reads
	cat SO002_PS_blocks/$block/${block}.sam |\
	awk '{if($1 ~ /^@/) print $0; else if ($0 ~ /HP:i:1/) print $0}' |\
	grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip > SO002_PS_blocks/$block/${block}.HP1.fq.gz
	# same for HP2
	cat SO002_PS_blocks/$block/${block}.sam |\
	awk '{if($1 ~ /^@/) print $0; else if ($0 ~ /HP:i:2/) print $0}' |\
	grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip > SO002_PS_blocks/$block/${block}.HP2.fq.gz
	
	# clean up
	rm SO002_PS_blocks/$block/${block}.sam
done
# now there should be two gzipped fastq files per block in 'SO002_PS_blocks/*'


## Bash script to get statistics
# in the google drive there is a bash script "SO002_calculating_all_stats.sh"
# it's geared up to run for the test but the input files are easily modified at the top
bash SO002_calculating_all_stats.sh
# now you will have a directory "SO002_haploid_stats" with a summary file "haploid_phasing_stats.csv"
# my summary file is in google drive

## R script needs to become automatised
# for this the script needs some data we can produce in bash
# first create a phased vcf containing only the phased variants which should be the majority
zcat SO002_haploid.phased.vcf.gz | awk '{if($1 ~ /^#/) print $0; else if($9 ~ /:PS/) print $0}' | gzip > SO002_haploid.phased_only.vcf.gz
# file containing information on the variants
zcat SO002_haploid.phased_only.vcf.gz |\
awk '{if($1 !~ /^#/) {print $1"\t"$2"\t"$4"\t"$5"\t"$10}}' |\
awk '{if(length($3)>1 || length($4)>1) {print $1"\t"$2"\t"$5"\tINDEL"} else {print $1"\t"$2"\t"$5"\tSNP"}}' |\
awk -F ":" '{print $1"\t"$6}' |\
awk 'BEGIN{print "contig\tpos\tblock\ttype"} {print $1"\t"$2"\t"$4"\t"$5}' > SO002_haploid.phased_only_GRAPHING.csv
# a bed file for the genome used to phase
zcat SO002_haploid.phased_only.vcf.gz |\
awk '{if($1 ~ /^##contig/) {print $0}}' |\
awk -F "=" 'BEGIN {print "contig\tend"} {print $3"\t"$4}' |\
sed 's/,length//g' |\
sed 's/>//g' > SO002_haploid.bed
# a bed file for the phased blocks (uses the whatshap stats output)
cat SO002_haploid.phased.gtf |\
awk 'BEGIN{print "contig\tblock\tstart\tend"}{print $1"\t"$4"\t"$4"\t"$5}' > SO002_haploid.phased_block.bed
# now you can run the rscript in the same directory
Rscript Phasing_plots_SO002.R



## Bash commands to join haplotypes
# get names of reads for both haplotypes
zcat SO002_PS_blocks/*/*.HP1.fq.gz | paste - - - - | sed 's/^@/>/g'|\
cut -f1-2 | tr '\t' '\n' | grep '>' |\
sed 's/>//g' > SO002_HP1.readnamesList.txt
zcat SO002_PS_blocks/*/*.HP2.fq.gz | paste - - - - | sed 's/^@/>/g'|\
cut -f1-2 | tr '\t' '\n' | grep '>' |\
sed 's/>//g' > SO002_HP2.readnamesList.txt
# for each haplotype, remove from all reads the reads associated with the opposite haplotype
# this creates our haplotype specific set of reads!
cat long_reads.fq.gz | paste - - - - | grep -v -F -f SO002_HP1.readnamesList.txt | tr "\t" "\n" | gzip > SO002_HP2more.fq.gz
cat long_reads.fq.gz | paste - - - - | grep -v -F -f SO002_HP2.readnamesList.txt | tr "\t" "\n" | gzip > SO002_HP1more.fq.gz


## Join haplotype assemblies into diploid assembly (need to rename contigs)
ls SO002_HP1.canu_r1_m2.fa SO002_HP2.canu_r1_m2.fa | while read genomes
do
hp=$( echo $genomes | awk -F "." '{print $1}' | awk -F "_" '{print $2}' )
cat $genomes | awk -v hp="$hp" '{if($0 ~ ">") print $1"_"hp ; else print $0}' >> SO002_diploid.canu_r1_m2.fa
done


##split the assembly into seperate haplotypes again
cat SO002_diploid.canu_r1_m2_p3.fa |\
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' |\
grep "\S" |\
sed 's/_pilon_pilon_pilon//g' |\
grep -A 1 'HP1' > SO002_HP1.canu_r1_m2_p3.fa
cat SO002_diploid.canu_r1_m2_p3.fa |\
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' |\
grep "\S" |\
sed 's/_pilon_pilon_pilon//g' |\
grep -A 1 'HP2' > SO002_HP2.canu_r1_m2_p3.fa

