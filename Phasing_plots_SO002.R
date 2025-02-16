
#libraries
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('scales')) install.packages('scales'); library('scales')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('data.table')) install.packages('data.table'); library('data.table')

#read in files produced previously for this script
strain_bed=read.csv(file="SO002_haploid.bed", header=T, sep='\t')
strain_variants=read.csv(file="SO002_haploid.phased_only_GRAPHING.csv", header=T, sep='\t')
strain_blocks=read.csv(file="SO002_haploid.phased_block.bed", sep='\t', header = T)

##reformulate the datatable for better alphanumeric ordering
strain_bed=data.table(strain_bed)
strain_variants=data.table(strain_variants)
strain_blocks=data.table(strain_blocks)
#reorder contigs
strain_bed=strain_bed[gtools::mixedorder(as.character(strain_bed$contig))]
strain_variants=strain_variants[gtools::mixedorder(as.character(strain_variants$contig))]
strain_blocks=strain_blocks[gtools::mixedorder(as.character(strain_blocks$contig))]
#create list of reordered contigs for plotting
chrom_order= as.character(strain_bed$contig)
chrom_key = setNames(object=as.character(strain_bed$contig), nm=chrom_order)

##create the chromosome order in the files
strain_bed[["contig"]] = factor(x=strain_bed[["contig"]], levels=chrom_order)
strain_variants[["contig"]] = factor(x=strain_variants[["contig"]], levels=chrom_order)
strain_blocks[["contig"]] = factor(x=strain_blocks[["contig"]], levels=chrom_order)

##plot
pdf("SO002_haploid.phasing.pdf")
ggplot(data=strain_bed)+
  geom_rect(aes(xmin = as.numeric(contig) - 0.2,
                xmax = as.numeric(contig) + 0.2,
                ymax = end, ymin = 0),
                colour="black", fill = "white") + 
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "black", size=12),
        axis.text.y = element_text(colour = "black", size=11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = expression("SO002 "~italic(de-novo)~"assembly"), limits = names(chrom_key))+
  # add bands to the blocks for each SNP or INDEL
  geom_rect(data = strain_variants, aes(xmin = as.numeric(contig) - 0.2, 
                                   xmax = as.numeric(contig) + 0.2, 
                                   ymax = pos+100, ymin = pos-100, fill = type))+
  scale_fill_manual(values = c("#D55E00", "#56B4E9")) +
  geom_rect(data = strain_blocks, aes(xmin = as.numeric(contig) + 0.4, 
                                   xmax = as.numeric(contig) + 0.25, 
                                   ymax = end, ymin = start),fill="grey50", alpha=0.75)+
  ylab("position (bp)")

dev.off() 
