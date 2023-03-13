## Meredith Persico: mxp1191@psu.edu
## Code for: "The age of absorptive roots impacts root-adjacent microbial composition in grapevines"
## Post dada2 Data Analysis
## Initial exploration of data (PCoA)
## ITS

## Set directory & check files
setwd("/Users/mjp23/Desktop/RootAgePostDada2")

list.files("/Users/mjp23/Desktop/RootAgePostDada2")

### Load required packages ###
library(ade4)
library(RColorBrewer)
library(phyloseq)
library(plyr)
library(gdata)
library(ggplot2)
library(vegan)
library(agricolae)
library(tidyverse)
library(gtools)

##install.packages("janitor")
library(janitor)

## Import data
tax <- read.table(file="Ph_RootAge_its_ESV_Taxonomy.txt", header=T, row.names=1)
tax ## looks good

asv_tab <- read.table(file="Ph_RootAge_its_ESV_Abund_Table.txt", header=T, row.names=1) 
asv_tab ## asv don't need to transpose, but need to remove blanks (95,96,112)
data.frame(rownames(asv_tab)) ## want to remove its87,88,112 (blanks in metadata file)
asv <- asv_tab[-c(86,87,103),] ## these rows correspond to the above samples
asv ## looks good; no blank samples

meta <- read.table("Ph_RootAge_its_metadata.txt", header=T, row.names=1)
meta ## need to remove blanks
data.frame(rownames(meta)) ## blanks correspond to same rows as above
metadat <- meta[-c(86,87,103),]
metadat ## looks good; blanks are removed
metadat$Root.Age <- as.numeric(as.character(metadat$Root.Age)) ## changes root age to numeric value
is.numeric(metadat$Root.Age) ## confirms yes root age is numeric

## Create new columns in metadata
metadat$Vine <- paste(metadat$Row, metadat$Block) ## created a new column, "Vine"
metadat$YoungvsOld <- ifelse(metadat$Root.Age<11,"Young", "Old") ## Creates a new, old vs. young based on Volder et al.
metadat$YoungvsOld

metadat ## looks good

## Remove problematic taxa
taxon<- tax
taxon[is.na(taxon)]<-"Unclassified"
taxon ## Now has "Unclassified" instead of NA

unclass <- subset(taxon, Phylum == "Unclassified")

rownames(unclass) ##279 unclassified rows

un.0 <- rownames(unclass)

red.taxon <- taxon[-which(rownames(taxon) %in% un.0),]
red.taxon 

## Match asv table with new taxonomy file
asv.filt<-asv[,which(colnames(asv) %in% rownames(red.taxon))]
asv.filt ## just looking at what comes up

## Calculate reads per sample
rowSums(asv.filt) ## several have no samples

## Determine minimum number of reads in your data set
min(rowSums(asv.filt))

## Cut samples with too few reads
n2 <- names(which(rowSums(asv.filt) > 1000)) 
asv.re<-asv.filt[which(rownames(asv.filt) %in% n2),]
min(rowSums(asv.re)) ## min now 1158 reads

metadat.red <-metadat[which(rownames(metadat) %in% rownames(asv.re)),] ## match metadata with removed samples

## Rarefy to obtain even numbers of reads by sample 
set.seed(336)
asv.r<-rrarefy(asv.re, 1158)
asv.r

## Convert asvs to percentages
asv.perc<-asv.r/rowSums(asv.r)*100
asv.perc

## Set factors as factors before analysis
# For whole dataset
metadat.red$Block <-as.factor(metadat.red$Block)
metadat.red$Row <- as.factor(metadat.red$Row)
metadat.red$Root.type <- as.factor(metadat.red$Root.type)
metadat.red$Age.testing.or.variation.testing. <- as.factor(metadat.red$Age.testing.or.variation.testing.)
metadat.red$Vine <- as.factor(metadat.red$Vine)
metadat.red$Mass..mg.<- as.factor(metadat.red$Mass..mg.)
metadat.red$Color <-as.factor(metadat.red$Color)
metadat.red$YoungvsOld<-as.factor(metadat.red$YoungvsOld)

## Counting root and vine distribution of full data set
metadat.red %>% 
  group_by(metadat.red$Vine) %>% 
  summarize(number_rows=n())
metadat.red %>%
  group_by(metadat.red$Root.Age) %>%
  summarize(number_rows=n())

tabyl(metadat.red, Root.Age, Vine)
tabyl(age.met, Root.Age, Vine)

## Group the asv, metadata, and taxonomy files
asv.perc<-as.data.frame(asv.perc)
asv.perc

## create phyloseq object
ITS.ASV <- otu_table(asv.perc, taxa_are_rows=F)
ITS.Meta <- sample_data(metadat.red,errorIfNULL=TRUE)
ITS.Taxo <- tax_table(as.matrix(red.taxon), errorIfNULL=TRUE)
ITS.Phylo <- phyloseq(ITS.ASV, ITS.Meta,ITS.Taxo)

### PCoA and NMDS USING PHYLOSEQ
RootageITS.ord <- ordinate(ITS.Phylo, "PCoA", "bray")

p1 = plot_ordination(ITS.Phylo, RootageITS.ord, type="taxa", color="Phylum", title="All Roots: Fungal Phylums")
print(p1)

p2 =plot_ordination(ITS.Phylo, RootageITS.ord, type="samples", color= "YoungvsOld") + geom_point(size=4) +
  labs(title="All Roots ITS: Young (< 11 days) vs. Old (11.5-40 days)") 
print(p2)

## Incorporating vine as color; root age as shape
p3= plot_ordination(ITS.Phylo, RootageITS.ord, type="samples", color="Vine", shape = "YoungvsOld") + geom_point(size=4)+ 
  labs(title="PCoA ITS: All sampled roots (n = 85)") +
  theme_bw()
print(p3) 


## Data shows separation of roots by color (vine); therefore, subsequent analysis comparing young vs. old roots will 
## factor in the variation caused by vine; furthermore, vine 21A only has one root, and it will be removed.

