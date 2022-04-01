#' ---
#' title: "Cullen_et_al_2022_animal_open_space"
#' author: "James Cullen"
#' date: "March 7th, 2022"
#' ---

#' This script relates to Cullen et al. (2022) published in Animal - Open Space.
#' The 'DNA quantity analysis' scripts were used for plotting DNA extraction quantites.
#' The 'Bacterial (16S) and fungal (ITS2) sequencing analysis' scripts were used for 
#' importing and analysing QIIME2 artifacts generated during using the script available
#' at https://github.com/JamesTCullen/Cullen_et_al_2022_animal_open_space/blob/main/QIIME2_analysis.sh


#' DNA quantity analysis

#load ggplot2 for plotting
library(ggplot2)

#Create data.frame of data from .txt file 
dna.df <- read.table(file = 'meth_dev_quantity_metadata.txt', header = TRUE)

#Define levels of BB variable for plotting
dna.df$Bead_beating <- factor(dna.df$Bead_beating, levels=c("BB0","BB3","BB10", "BB15", "BB20"))

DNA_plot <- ggplot(dna.df, aes(x=Bead_beating, y=DNA_concentration_ng_ul)) + 
  geom_boxplot(aes(fill=Bead_beating)) + stat_summary(fun = "mean", colour = "black", size = 1.5, geom = "point") +
  facet_wrap(.~Sample_type, scales = "free") + theme_bw(base_size = 16) +
  labs(x="Bead beating duration",
       y="DNA concentration (ng/µL)",
       fill="Bead beating") + 
  theme(axis.text=element_text(size=16), 
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 40,
                                   hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        legend.position = "top")

#' Bacterial (16S) and fungal (ITS2) sequencing analysis

#Load packages
library(qiime2R)#to import QIIME2 artifacts into phyloseq
library(phyloseq)
library(dplyr)

#Build phyloseq object with artifacts exported from QIIME2 for 16S data (ASV table, 
#taxonomy of associated ASVs, phylogenetic tree and sample metadata) 
p<-qza_to_phyloseq(features="meth-dev-16S-trimmed-table.qza",
                   taxonomy="meth-dev-16S-trimmed-taxonomy.qza", 
                   metadata="meth-dev-metadata-16S.tsv", 
                   tree = "meth-dev-16S-trimmed-rooted-tree.qza")

#Inspect 16S phyloseq object
p

#Build phyloseq object with artifacts exported from QIIME2 for ITS data 
#(ASV table, taxonomy of associated ASVs, phylogenetic tree and sample metadata)
i<-qza_to_phyloseq(features="meth-dev-its-trimmed-table.qza",
                   taxonomy="meth-dev-its-trimmed-taxonomy.qza", 
                   metadata="meth-dev-metadata-its.tsv", 
                   tree = "meth-dev-its-trimmed-rooted-tree.qza")

#Inspect ITS phyloseq object
i

#Inspect feature tables
#16S
tail(tax_table(p),10)

#ITS
tail(tax_table(i),10)

#Check assignment of ASVs to different taxonomic levels
#16S
apply(tax_table(p)[,2:7],2,function(x){1-mean(is.na(x))})

#ITS
apply(tax_table(i)[,2:7],2,function(x){1-mean(is.na(x))})

#Check metadata
#16S
colnames(sample_data(p))

#ITS
colnames(sample_data(i))

#Renaming ASVs  
#16S
taxa_names(p) <- paste0("ASV16S", seq(ntaxa(p)))

#ITS
taxa_names(i) <- paste0("ASVITS", seq(ntaxa(i)))

#Use decontam to identify contaminant ASVs
library(decontam)

#Check library sizes as a function of true samples and controls
#16S
df <- as.data.frame(sample_data(p)) 
df$LibrarySize <- sample_sums(p)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample.or.Control)) + geom_point()

#ITS
df <- as.data.frame(sample_data(i)) 
df$LibrarySize <- sample_sums(i)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample.or.Control)) + geom_point()

#Identify contaminants in both data sets (frequency method)
#16S
contamdf.freq.p <- isContaminant(p, method="frequency", conc="indexed.lib.conc")

#ITS
contamdf.freq.i <- isContaminant(i, method="frequency", conc="indexed.lib.conc")

#This produces a data.frame: $p is the probability for classifying contaminants
#$contaminant=TRUE if $p < 0.1 (statistical evidence ASV is a contaminant) 
#Check if contaminants were ID'd
#16S
table(contamdf.freq.p$contaminant) #27 contaminant ASVs

#ITS
table(contamdf.freq.i$contaminant) #3 contaminant ASVs

#Check which ASVs were ID'd
#16S
head(which(contamdf.freq.p$contaminant))

#ITS
head(which(contamdf.freq.i$contaminant))

#Remove contaminant ASVs
#16S
p.noncontam <- prune_taxa(!contamdf.freq.p$contaminant, p)

#ITS
i.noncontam <- prune_taxa(!contamdf.freq.i$contaminant, i)

#Check decontaminated phyloseq objects
#16S
p.noncontam
#ITS
i.noncontam

#Define levels of bead beating variable for correct order when plotting
#16S
sample_data(p.noncontam)$bead.beating<-factor(sample_data(p.noncontam)$bead.beating, levels=c("BB0","BB3","BB10", "BB15", "BB20"))

#ITS
sample_data(i.noncontam)$bead.beating<-factor(sample_data(i.noncontam)$bead.beating, levels=c("BB0","BB3","BB10", "BB15", "BB20"))

#Check levels of beat beating variable
#16S
levels(sample_data(p.noncontam)$bead.beating)
#ITS
levels(sample_data(i.noncontam)$bead.beating)

#Removing unwanted taxa from 16S
p.1 <- p.noncontam %>% subset_taxa(Genus!="Chloroplast" | is.na(Genus))

p.2 <- p.1 %>% subset_taxa(Genus!="Mitochondria" | is.na(Genus))

p.3 <- subset_taxa(p.2, Kingdom != "Unassigned")

p.4 <- subset_taxa(p.3, Kingdom != "d__Eukaryota")

p.5 <- subset_taxa(p.4, Kingdom != "d__Archaea")

#Post-QC library stats
#16S
mean(sample_sums(p.5))
#16S
range(sample_sums(p.5))

#ITS
mean(sample_sums(i.noncontam))
#ITS
range(sample_sums(i.noncontam))

#Make a data.frame of read depths
#16S
reads.p<-data.frame(reads=sample_sums(p.5))
#ITS
reads.i<-data.frame(reads=sample_sums(i.noncontam))

#Add on the sample IDs
reads.p$Sample<-rownames(reads.p) #16S
reads.i$Sample<-rownames(reads.i) #ITS

#Extract the metadata from the phyloseq object 
meta.p<-data.frame(sample_data(p.5)) #16S
meta.p$Sample<-rownames(meta.p) #16S

meta.i<-data.frame(sample_data(i.noncontam)) #ITS
meta.i$Sample<-rownames(meta.i) #ITS

#Join on the Metadata
reads.p<-left_join(reads.p,meta.p,"Sample") #16S
reads.i<-left_join(reads.i,meta.i,"Sample") #ITS

#Some boxplots to visualise reads by sample type and bead beating time
#Plotting reads by bead beating time
bb.reads.p <- ggplot(reads.p,aes(x=bead.beating,y=reads)) + geom_boxplot(aes(fill=bead.beating))#16S

bb.reads.i <- ggplot(reads.i,aes(x=bead.beating,y=reads)) + geom_boxplot(aes(fill=bead.beating))#ITS

#Plotting reads by sample type
bb.sample.p <- ggplot(reads.p,aes(x=sample.type,y=reads)) + geom_boxplot(aes(fill=sample.type))#16S

bb.sample.i <- ggplot(reads.i,aes(x=sample.type,y=reads)) + geom_boxplot(aes(fill=sample.type))#ITS

#What are the mean number of reads?
mean(sample_sums(p.5)) #16S

mean(sample_sums(i.noncontam)) #ITS

#Which samples have a low number of reads?
names(sample_sums(p.5))[which(sample_sums(p.5)<15000)] #16S

names(sample_sums(i.noncontam))[which(sample_sums(i.noncontam)<10000)] #ITS

#Remove sample TLF-BB3-3-16S and TLF-BB10-1-ITS due to low number of reads 
p.new = subset_samples(p.5, sample_names(p.5) != "TLF-BB3-3-16S") #16S

i.new = subset_samples(i.noncontam, sample_names(i.noncontam) != "TLF-BB10-1-ITS") #ITS

#Number of reads for final datasets 
mean(sample_sums(p.new)) #16S
mean(sample_sums(i.new)) #ITS

#Rarefaction curves for both datasets
library(vegan)
rarecurve.p <- rarecurve(t(otu_table(p.new)), step=50, cex=0.5) #16S

rarecurve.i <- rarecurve(t(otu_table(i.new)), step=50, cex=0.5) #ITS

#' Alpha diversity analysis

#Prune phyloseq objects to subset sample types for both datasets
#16S - Faeces
p.f<-prune_samples(sample_data(p.new)$sample.type =="Faeces",p.new) #15 samples

#16S - Liquid Feed
p.lf<-prune_samples(sample_data(p.new)$sample.type =="Liquid-feed",p.new) #14 samples

#ITS - Faeces
i.f<-prune_samples(sample_data(i.new)$sample.type =="Faeces",i.new) #15 samples

#ITS - Liquid Feed
i.lf<-prune_samples(sample_data(i.new)$sample.type =="Liquid Feed",i.new) #14 samples

#Create data.frames of alpha diversity measures with metadata for each object
#16S - Faeces
adiv.p.f <- data.frame(
  "Observed" = estimate_richness(p.f, measures = "Observed"),
  "Shannon" = estimate_richness(p.f, measures = "Shannon"),
  "Simpson" = estimate_richness(p.f, measures = "Simpson"),
  "SampleType" = phyloseq::sample_data(p.f)$sample.type, 
  "Bead_Beating" = phyloseq::sample_data(p.f)$bead.beating)

#16S - Liquid Feed
adiv.p.lf <- data.frame(
  "Observed" = estimate_richness(p.lf, measures = "Observed"),
  "Shannon" = estimate_richness(p.lf, measures = "Shannon"),
  "Simpson" = estimate_richness(p.lf, measures = "Simpson"),
  "SampleType" = phyloseq::sample_data(p.lf)$sample.type, 
  "Bead_Beating" = phyloseq::sample_data(p.lf)$bead.beating)

#ITS - Faeces
adiv.i.f <- data.frame(
  "Observed" = estimate_richness(i.f, measures = "Observed"),
  "Shannon" = estimate_richness(i.f, measures = "Shannon"),
  "Simpson" = estimate_richness(i.f, measures = "Simpson"),
  "SampleType" = phyloseq::sample_data(i.f)$sample.type, 
  "Bead_Beating" = phyloseq::sample_data(i.f)$bead.beating)

#ITS - Liquid Feed
adiv.i.lf <- data.frame(
  "Observed" = estimate_richness(i.lf, measures = "Observed"),
  "Shannon" = estimate_richness(i.lf, measures = "Shannon"),
  "Simpson" = estimate_richness(i.lf, measures = "Simpson"),
  "SampleType" = phyloseq::sample_data(i.lf)$sample.type, 
  "Bead_Beating" = phyloseq::sample_data(i.lf)$bead.beating)

#Use alpha_boxplot function to create individual alpha diversity plots
library(amplicon)

#16S - Faeces
p.f.alpha <- alpha_boxplot(adiv.p.f, adiv.p.f, index = "Shannon", groupID = "Bead_Beating") +
  theme_bw() +
  labs(x="Bead beating duration",
       y="Shannon",
       fill="Bead beating") + 
  theme(axis.text=element_text(size=16), 
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 16)) + 
  ggtitle(label = "Bacteria in faeces") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_blank(), legend.position = "none")

#16S - Liquid Feed
p.lf.alpha <- alpha_boxplot(adiv.p.lf, adiv.p.lf, index = "Shannon", groupID = "Bead_Beating") +
  theme_bw() +
  labs(x="Bead beating duration",
       y="Shannon",
       fill="Bead beating") + 
  theme(axis.text=element_text(size=16), 
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 16)) + 
  ggtitle(label = "Bacteria in liquid feed") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        legend.title = element_blank(), legend.position = "none")

#ITS - Faeces
i.f.alpha <- alpha_boxplot(adiv.i.f, adiv.i.f, index = "Shannon", groupID = "Bead_Beating") +
  theme_bw() +
  labs(x="Bead beating duration",
       y="Shannon",
       fill="Bead beating") + 
  theme(axis.text=element_text(size=16), 
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 16)) + 
  ggtitle(label = "Fungi in faeces") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        legend.title = element_blank(), legend.position = "none")

#ITS - Liquid Feed
i.lf.alpha <- alpha_boxplot(adiv.i.lf, adiv.i.lf, index = "Shannon", groupID = "Bead_Beating") +
  theme_bw() +
  labs(x="Bead beating duration",
       y="Shannon",
       fill="Bead beating") + 
  theme(axis.text=element_text(size=16), 
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 16)) + 
  ggtitle(label = "Fungi in liquid feed") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        legend.title = element_blank(), legend.position = "none")

library(patchwork)

#Combine plots
p.f.alpha + i.f.alpha + p.lf.alpha + i.lf.alpha

#' Beta diversity analysis 

#Agglomerate to genus level
p.f.genus<-tax_glom(p.f, taxrank = "Genus") #16S - Faeces

p.lf.genus<-tax_glom(p.lf, taxrank = "Genus") #16S - Liquid Feed

i.f.genus<-tax_glom(i.f, taxrank = "Genus") #ITS - Faeces

i.lf.genus<-tax_glom(i.lf, taxrank = "Genus") #ITS - Liquid Feed


###Ordination using built in functions in phyloseq (calls vegan)    
ord.bray.p.f <- ordinate(p.f.genus, method="PCoA",k=2, 
                         distance="bray",trymax=50) #16S - Faeces
ord.bray.p.lf <- ordinate(p.lf.genus, method="PCoA",k=2, 
                          distance="bray",trymax=50) #16S - Liquid Feed
ord.bray.i.f <- ordinate(i.f.genus, method="PCoA",k=2, 
                         distance="bray",trymax=50) #ITS - Faeces
ord.bray.i.lf <- ordinate(i.lf.genus, method="PCoA",k=2, 
                          distance="bray",trymax=50) #ITS - Liquid Feed


#Make ordination plots and store them
#16S - Faeces
ord.p.f<-plot_ordination(p.f.genus, ord.bray.p.f, color="bead.beating", 
                         shape = "sample.type", title="Bacteria in faeces")
#16S - Liquid Feed
ord.p.lf<-plot_ordination(p.lf.genus, ord.bray.p.lf, color="bead.beating", 
                          shape = "sample.type", title="Bacteria in liquid feed")
#ITS - Faeces
ord.i.f<-plot_ordination(i.f.genus, ord.bray.i.f, color="bead.beating", 
                         shape = "sample.type", title="Fungi in faeces")
#ITS - Liquid Feed
ord.i.lf<-plot_ordination(i.lf.genus, ord.bray.i.lf, color="bead.beating", 
                          shape = "sample.type", title="Fungi in liquid feed")

#Plot, format and combine plots using patchwork 
#16S - Faeces
plot1 <- ord.p.f + theme_bw(base_size = 16) + geom_point(size=4) + 
  labs(color="Bead beating") + theme(legend.position = "none", 
                                     legend.title = element_text(), 
                                     title = element_text(size = 14))
#16S - Liquid Feed
plot2 <- ord.p.lf + theme_bw(base_size = 16) + geom_point(size=4) + 
  labs(color="Bead beating") + theme(legend.position = "none", 
                                     legend.title = element_text(), 
                                     title = element_text(size = 14))
#ITS - Faeces
plot3 <- ord.i.f + theme_bw(base_size = 16) + geom_point(size=4) + 
  labs(color="Bead beating") + theme(legend.position = "none", 
                                     legend.title = element_text(),
                                     title = element_text(size = 14))
#ITS - Liquid Feed
plot4 <- ord.i.lf + theme_bw(base_size = 16) + geom_point(size=4) + 
  labs(color="Bead beating") + theme(legend.position = "none", 
                                     legend.title = element_text(), 
                                     title = element_text(size = 14))

#Combine all plots for figure
bdivplot <- plot1 + plot2 + plot3 + plot4

#Add legend
bdivplot + theme(legend.position = "bottom", legend.direction = "horizontal")


#' Differential abundance analysis 

#Need to identify which bacterial and fungal genera are >1% mean relative 
#abundance at each BB time for each sample type. 
#Then manually filter those from the phyloseq objects before running DESeq2
#(phyloseq objects already agglomerated to genus level from beta diversity 
#analysis above)

p.f.genus #16S - Faeces
p.lf.genus #16S - Liquid Feed
i.f.genus #ITS - Faeces
i.lf.genus #ITS - Liquid Feed

#Transform to relative abundance and create a data.frame of mean relative 
#abundance with respect to genus and BB procedure

#16S - Faeces
p.f.genus.rel <- transform_sample_counts(p.f.genus, function(x){x/sum(x)})
p.f.genus.rel.df = psmelt(p.f.genus.rel)
p.f.genus.rel.df.agr = aggregate(Abundance~bead.beating+Genus, data=p.f.genus.rel.df, FUN=mean)
write.csv(p.f.genus.rel.df.agr,'p.f.genus.rel.csv') #export the dataframe as .csv

#16S - Liquid Feed
p.lf.genus.rel <- transform_sample_counts(p.lf.genus, function(x){x/sum(x)})
p.lf.genus.rel.df = psmelt(p.lf.genus.rel)
p.lf.genus.rel.df.agr = aggregate(Abundance~bead.beating+Genus, data=p.lf.genus.rel.df, FUN=mean)
write.csv(p.lf.genus.rel.df.agr,'p.lf.genus.rel.csv') #export the dataframe as .csv

#ITS - Faeces
i.f.genus.rel <- transform_sample_counts(i.f.genus, function(x){x/sum(x)})
i.f.genus.rel.df = psmelt(i.f.genus.rel)
i.f.genus.rel.df.agr = aggregate(Abundance~bead.beating+Genus, data=i.f.genus.rel.df, FUN=mean)
write.csv(i.f.genus.rel.df.agr,'i.f.genus.rel.csv') #export the dataframe as .csv

#ITS - Liquid Feed
i.lf.genus.rel <- transform_sample_counts(i.lf.genus, function(x){x/sum(x)})
i.lf.genus.rel.df = psmelt(i.lf.genus.rel)
i.lf.genus.rel.df.agr = aggregate(Abundance~bead.beating+Genus, data=i.lf.genus.rel.df, FUN=mean)
write.csv(i.lf.genus.rel.df.agr,'i.lf.genus.rel.csv') #export the dataframe as .csv

#Creating objects containing only genera >1% mean relative abundance across 
#BB procedures for each sample type

#16S - Faeces
p.f.keep=c("Prevotella","Clostridium_sensu_stricto_1","Prevotellaceae_NK3B31_group","Muribaculaceae",
           "Lactobacillus","Treponema","Alloprevotella","Rikenellaceae_RC9_gut_group","Phascolarctobacterium",
           "UCG-005","UCG-010","Anaerovibrio","Terrisporobacter","Christensenellaceae_R-7_group",
           "Clostridia_vadinBB60_group","Streptococcus","Blautia","Succinivibrio","Clostridia_UCG-014",
           "WCHB1-41","Prevotellaceae_UCG-003","UCG-002","Faecalibacterium","NK4A214_group","Campylobacter",
           "Ruminococcus","Subdoligranulum","Prevotellaceae_UCG-003","Gastranaerophilales","Megasphaera",
           "[Eubacterium]_coprostanoligenes_group","Parabacteroides")
gen.p.f = subset_taxa(p.f, Genus %in% p.f.keep)

#Agglomerate at genus level
gen.p.f.g <- tax_glom(gen.p.f, taxrank = "Genus")

#16S - Liquid Feed
p.lf.keep=c("Clostridium_sensu_stricto_1","Lactobacillus","Weissella","Leuconostoc","Pediococcus")
gen.p.lf = subset_taxa(p.lf, Genus %in% p.lf.keep)

#Agglomerate at genus level
gen.p.lf.g <- tax_glom(gen.p.lf, taxrank = "Genus")

#ITS - Faeces
i.f.keep=c("Kazachstania","Mucor","Gamsia","Pichia","Saccharomyces","Wallemia","Cladosporium","Monascus",
           "Scopulariopsis","Malassezia","Acaulium","Penicillium","Pilaira","Aspergillus","Naganishia",
           "Alternaria","Debaryomyces","Sporobolomyces","Skeletocutis","Vishniacozyma","Papiliotrema",
           "Peniophora","Rickenella","Filobasidium","Trichosporon") 
gen.i.f = subset_taxa(i.f, Genus %in% i.f.keep)

#Agglomerate at genus level
gen.i.f.g <- tax_glom(gen.i.f, taxrank = "Genus")

#ITS - Liquid Feed
i.lf.keep=c("Kazachstania","Alternaria","Saccharomyces","Neoascochyta","Monographella","Cladosporium",
            "Vishniacozyma","Fusarium","Diaporthe","Gibberella","Pichia")
gen.i.lf = subset_taxa(i.lf, Genus %in% i.lf.keep)

#Agglomerate at genus level
gen.i.lf.g <- tax_glom(gen.i.lf, taxrank = "Genus")

#Define function to handle results objects and annotate with taxonomy
taxo<-function(resultsobject,physeqobject,alpha){
  sigtab<-resultsobject[which(resultsobject$padj<alpha),]
  sigtab<- cbind(as(sigtab, "data.frame"), as(tax_table(physeqobject)[rownames(sigtab), ], "matrix"))
  colnames(sigtab)[7:12]<-c("Kingdom","Phylum","Class","Order","Family","Genus")
  return(sigtab)
}          

#Define function to calculate geometric means
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}   


#Function to parse significant data from DESeq2 results
deseqplot_data<-function(sigtab){
  # Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  
#Copy across genus labels and fill in any unassigned
  sigtab$Genus.long<-as.character(sigtab$Genus)
  sigtab$Genus.long[grep("unclassified",sigtab$Genus)] <-paste0("[",as.character(sigtab$Family[grep("unclassified",sigtab$Genus)]),"]")
  
  #Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus.long, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus.long = factor(as.character(sigtab$Genus.long), levels=names(x))
  return(sigtab)
}

library(DESeq2)

#Fit The Model We Are Interested in Using The function in phyloseq - gets the data in a format DESeq can read
p.f.mod<-phyloseq_to_deseq2(gen.p.f.g, ~ bead.beating) #16S - Faeces
p.lf.mod<-phyloseq_to_deseq2(gen.p.lf.g, ~ bead.beating) #16S - Liquid Feed
i.f.mod<-phyloseq_to_deseq2(gen.i.f.g, ~ bead.beating) #ITS - Faeces
i.lf.mod<-phyloseq_to_deseq2(gen.i.lf.g, ~ bead.beating) #ITS - Liquid Feed


# Calculate geometric mean counts   
#16S - Faeces
p.f_ge<-apply(DESeq2::counts(p.f.mod),1,gm_mean) 
p.f_size<-estimateSizeFactors(p.f.mod,geoMeans=p.f_ge) 

#16S - Liquid Feed
p.lf_ge<-apply(DESeq2::counts(p.lf.mod),1,gm_mean) 
p.lf_size<-estimateSizeFactors(p.lf.mod,geoMeans=p.lf_ge) #ITS - Liquid Feed

#ITS - Faeces
i.f_ge<-apply(DESeq2::counts(i.f.mod),1,gm_mean) 
i.f_size<-estimateSizeFactors(i.f.mod,geoMeans=i.f_ge) 

#ITS - Liquid Feed
i.lf_ge<-apply(DESeq2::counts(i.lf.mod),1,gm_mean) 
i.lf_size<-estimateSizeFactors(i.lf.mod,geoMeans=i.lf_ge) 

#Fit the DESeq model   
#16S - Faeces
p.f.deseq = DESeq(p.f_size, test="Wald", fitType="mean")

#16S - Liquid Feed
p.lf.deseq = DESeq(p.lf_size, test="Wald", fitType="mean")

#ITS - Faeces
i.f.deseq = DESeq(i.f_size, test="Wald", fitType="mean")

#ITS - Liquid Feed
i.lf.deseq = DESeq(i.lf_size, test="Wald", fitType="mean")


#Extract the results for the contrasts we want
#We will compare BB3 vs BB0, BB10 vs BB3, BB15 vs BB10 & BB20 vs BB15 for each sample type for 16S and ITS

#16S - Faeces
p.f_results1<-results(p.f.deseq,contrast=c("bead.beating","BB3","BB0"))
p.f_results2<-results(p.f.deseq,contrast=c("bead.beating","BB10","BB3"))
p.f_results3<-results(p.f.deseq,contrast=c("bead.beating","BB15","BB10"))
p.f_results4<-results(p.f.deseq,contrast=c("bead.beating","BB20","BB15"))

#16S - Liquid Feed
p.lf_results1<-results(p.lf.deseq,contrast=c("bead.beating","BB3","BB0"))
p.lf_results2<-results(p.lf.deseq,contrast=c("bead.beating","BB10","BB3"))
p.lf_results3<-results(p.lf.deseq,contrast=c("bead.beating","BB15","BB10"))
p.lf_results4<-results(p.lf.deseq,contrast=c("bead.beating","BB20","BB15"))

#ITS - Faeces
i.f_results1<-results(i.f.deseq,contrast=c("bead.beating","BB3","BB0"))
i.f_results2<-results(i.f.deseq,contrast=c("bead.beating","BB10","BB3"))
i.f_results3<-results(i.f.deseq,contrast=c("bead.beating","BB15","BB10"))
i.f_results4<-results(i.f.deseq,contrast=c("bead.beating","BB20","BB15"))

#ITS - Liquid Feed
i.lf_results1<-results(i.lf.deseq,contrast=c("bead.beating","BB3","BB0"))
i.lf_results2<-results(i.lf.deseq,contrast=c("bead.beating","BB10","BB3"))
i.lf_results3<-results(i.lf.deseq,contrast=c("bead.beating","BB15","BB10"))
i.lf_results4<-results(i.lf.deseq,contrast=c("bead.beating","BB20","BB15"))

#Annotate the taxa with the taxonomy  
#16S - Faeces
p.f_results_taxa1<-taxo(p.f_results1,p.f.genus,0.05) ###Cutoff for P value 
p.f_results_taxa2<-taxo(p.f_results2,p.f.genus,0.05) ###Cutoff for P value 
p.f_results_taxa3<-taxo(p.f_results3,p.f.genus,0.05) ###Cutoff for P value 
p.f_results_taxa4<-taxo(p.f_results4,p.f.genus,0.05) ###Cutoff for P value 

#Write results to .csv files
#16S - Faeces
write.csv(p.f_results_taxa1, 'p.f_results_taxa1.csv')
write.csv(p.f_results_taxa2, 'p.f_results_taxa2.csv')
write.csv(p.f_results_taxa3, 'p.f_results_taxa3.csv')
write.csv(p.f_results_taxa4, 'p.f_results_taxa4.csv')

#16S - Liquid Feed
p.lf_results_taxa1<-taxo(p.lf_results1,p.lf.genus,0.05) #Cutoff for P value 
p.lf_results_taxa2<-taxo(p.lf_results2,p.lf.genus,0.05) 
#p.lf_results_taxa3<-taxo(p.lf_results3,p.lf.genus,0.05) #not significant
#p.lf_results_taxa4<-taxo(p.lf_results4,p.lf.genus,0.05) #not significant

#Write results to .csv files
write.csv(p.lf_results_taxa1, 'p.lf_results_taxa1.csv')
write.csv(p.lf_results_taxa2, 'p.lf_results_taxa2.csv')
#write.csv(p.lf_results_taxa3, 'p.lf_results_taxa3.csv') #not significant
#write.csv(p.lf_results_taxa4, 'p.lf_results_taxa4.csv') #not significant

#ITS - Faeces
i.f_results_taxa1<-taxo(i.f_results1,i.f.genus,0.05)  
i.f_results_taxa2<-taxo(i.f_results2,i.f.genus,0.05) 
#i.f_results_taxa3<-taxo(i.f_results3,i.f.genus,0.05) #not significant
i.f_results_taxa4<-taxo(i.f_results4,i.f.genus,0.05) 

#Write results to .csv files
write.csv(i.f_results_taxa1, 'i.f_results_taxa1.csv')
write.csv(i.f_results_taxa2, 'i.f_results_taxa2.csv')
#write.csv(i.f_results_taxa3, 'i.f_results_taxa3.csv') #not significant
write.csv(i.f_results_taxa4, 'i.f_results_taxa4.csv')

#ITS - Liquid Feed
i.lf_results_taxa1<-taxo(i.lf_results1,i.lf.genus,0.05)  
i.lf_results_taxa2<-taxo(i.lf_results2,i.lf.genus,0.05)  
i.lf_results_taxa3<-taxo(i.lf_results3,i.lf.genus,0.05)  
i.lf_results_taxa4<-taxo(i.lf_results4,i.lf.genus,0.05) 

#Write results to .csv files
write.csv(i.lf_results_taxa1, 'i.lf_results_taxa1.csv')
write.csv(i.lf_results_taxa2, 'i.lf_results_taxa2.csv')
write.csv(i.lf_results_taxa3, 'i.lf_results_taxa3.csv')
write.csv(i.lf_results_taxa4, 'i.lf_results_taxa4.csv')

################################### END ########################################












