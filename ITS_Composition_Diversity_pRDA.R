## Meredith Persico: mxp1191@psu.edu
## Code for: "The age of absorptive roots impacts root-adjacent microbial composition in grapevines"
### This code covers ITS beta composition, alpha diversity, and partial RDA


## Data import and cleaning
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
library(microbiome)

##install.packages("janitor")
library(janitor)

## Import data
tax <- read.table(file="Ph_RootAge_its_ESV_Taxonomy.txt", header=T, row.names=1)
tax ## looks good

asv_tab <- read.table(file="Ph_RootAge_its_ESV_Abund_Table.txt", header=T, row.names=1) 
asv_tab ## asv don't need to transpose, but need to remove blanks (95,96,112)
data.frame(rownames(asv_tab)) ## want to remove its87,88,112 (blanks in metadata file)
asv <- asv_tab[-c(12,33,40,47,52, 68, 86,87,91, 103),] ## these rows correspond to the above samples (and the second order roots from the metadat file
asv ## looks good; no blank samples

meta <- read.csv("Ph_RootAge_its_metadata.csv", header=T, row.names=1) ### added color to fungal data 3.17.22, make sure to use csv because .txt not necessarily changed
meta ## need to remove blanks & second order roots
data.frame(rownames(meta)) ## blanks correspond to same rows as above, 2nd order roots slightly different
metadat <- meta[-c(12,33,40,47,52, 68, 86,87,91, 103),] ## also removed 21A (1 root vine)
metadat ## looks good; blanks and second orders are removed
metadat$Root.Age <- as.numeric(as.character(metadat$Root.Age)) ## changes root age to numeric value
is.numeric(metadat$Root.Age) ## confirms yes root age is numeric

## Create new columns in metadata
## Create new columns in metadata
metadat$Vine <- paste(metadat$Row, metadat$Block) ## created a new column, "Vine"
metadat$AgeCat <- ifelse(metadat$Root.Age<13,"Young", "MiddleAge") ## Create new column for age category ("young" and "middle age")
metadat$AgeCat <- ifelse(metadat$Root.Age<27, metadat$AgeCat, "Old") ## Create new column for age category ("old")
metadat$AgeCat ## check categories are correct
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

## Skip transpose step

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
## ^ MP Note: Use this as metadata file

## Rarefy to obtain even numbers of reads by sample 
set.seed(336)
asv.r<-rrarefy(asv.re, 1158)
asv.r

## Convert asvs to percentages
asv.perc<-asv.r/rowSums(asv.r)*100
asv.perc

## Split up age and variation roots from metadata and asv.perc files
asv.fac <- cbind(metadat.red, asv.perc)
asv.age <- subset(asv.fac, Age.testing.or.variation.testing.%in% c("Age")) ## split into just age-related roots
asv.var <- subset(asv.fac, Age.testing.or.variation.testing.%in% c("Variation")) ## split to just variation-related roots

asv.age.perc <- asv.age[,12:1122] ## new ASV
age.met <- asv.age[,1:11] ## new metadata

## removing the residual zero asvs from when age and var datasets were combined
test <- t(asv.age.perc)
rowSums(test)
no_zeroes <- test[rowSums(test)>0,] 
no.zero.taxa <- red.taxon[which(rownames(red.taxon) %in% rownames(no_zeroes)),] ## remove zero taxa
asv.nozero <- t(no_zeroes) ## asv table
meta.nozeroes <- metadat.red[which(rownames(metadat.red) %in% rownames(asv.nozero)),] ## metadata no zeros

## "variation" roots dataset
asv.var.perc <- asv.var[,12:1122] ## new ASV
asv.var.met <- asv.var[,1:11] ## new metadata

## Just age roots data-set, 
## with residual info(e.g.,taxa not in just-age roots) removed from the full data-set
## w/variation roots 
test <- t(asv.age.perc)
rowSums(test)
no_var <- test[rowSums(test)>0,] ## anything with ALL zeroes removed (due to residual info from the "variation" dataset)

no.var.taxa <- red.taxon[which(rownames(red.taxon) %in% rownames(no_var)),] ## taxonomy
asv.novar <- t(no_var) ## asv
meta.novar <- metadat.red[which(rownames(metadat.red) %in% rownames(asv.novar)),] ## metadat

## Set factors as factors before analysis
# For no zero age roots
meta.novar$Vine <- as.factor(meta.novar$Vine)
meta.novar$YoungvsOld<-as.factor(meta.novar$YoungvsOld)

## create phyloseq object
novar.ITS.ASV.age <- otu_table(asv.novar, taxa_are_rows=F)
novar.ITS.Meta.age <- sample_data(meta.novar,errorIfNULL=TRUE)
novar.ITS.Taxo.age <- tax_table(as.matrix(no.var.taxa, errorIfNULL=TRUE))
novar.ITS.Phylo.age <- phyloseq(novar.ITS.ASV.age, novar.ITS.Meta.age, novar.ITS.Taxo.age) 
novar.ITS.Phylo.age ## 539 total taxa in just-age dataset

saveRDS(novar.ITS.Phylo.age, file = "NovarITS.rds")
readRDS("NovarITS.rds")

## Aggregate to Phylum and Class levels
phylum_data = aggregate_taxa(novar.ITS.Phylo.age, "Phylum")
class_data = aggregate_taxa(novar.ITS.Phylo.age, "Class")

phylum_data@tax_table
class_data@tax_table

df.t = data.frame(sample_data(novar.ITS.Phylo.age))
df.p = data.frame(sample_data(phylum_data))
df.c = data.frame(sample_data(class_data))


## 1. Beta Composition
### a. PERMANOVA
#### ASV
bc.tax = distance(novar.ITS.Phylo.age, "bray")
set.seed(10)
tax.perm = adonis(bc.tax ~ df.t$YoungvsOld, strata = df.t$Vine, data= df.t, permuations = 999)
tax.perm

#### Phylum
bc.phy = distance(phylum_data, "bray")
set.seed(10)
phy.perm = adonis(bc.phy~ YoungvsOld, data= df.p, permutations = 999)
phy.perm

#### Class
bc.class = distance(class_data, "bray")
set.seed(10)
class.perm = adonis(bc.class~ YoungvsOld, data= df.c, permutations = 999)
class.perm

## b. Graphs
### PCoA
#### ASV
ASV.ord.pc = ordinate(novar.ITS.Phylo.age, "PCoA", "bray")
pc.ASV <- plot_ordination(novar.ITS.Phylo.age, ASV.ord.pc, "samples", color="YoungvsOld")
pc.ASV <- pc.ASV + geom_point(size = 2) + 
  labs(title = "PCoA ITS: Selection of roots (n = 39)") +
  stat_ellipse(type="norm", level = 0.90, linetype=2) +
  theme_bw()
print(pc.ASV)

## with vine included
pc.ASV2 <- plot_ordination(novar.ITS.Phylo.age, ASV.ord.pc, "samples", color="Vine", shape = "YoungvsOld")
pc.ASV2 <- pc.ASV2 + geom_point(size = 5) + 
  labs(title = "PCoA ITS: Selection of roots (n = 39)") +
  theme_bw()
print(pc.ASV2)

#### Phylum
phy.ord.pc = ordinate(phylum_data, "PCoA", "bray")
pc.phy <- plot_ordination(phylum_data, phy.ord.pc, "samples", color="YoungvsOld")
pc.phy <- pc.phy + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(pc.phy)

#### class
class.ord.pc = ordinate(class_data, "PCoA", "bray")
pc.class <- plot_ordination(class_data, class.ord.pc, "samples", color="YoungvsOld")
pc.class <- pc.class + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(pc.class)

## CAP
#### ASV
tax.ord = ordinate(novar.ITS.Phylo.age, "CAP", "bray", ~YoungvsOld)
cap.tax <- plot_ordination(novar.ITS.Phylo.age, tax.ord, "samples", color="YoungvsOld")
cap.tax <- cap.tax + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(cap.tax)

#### Phylum
phy.ord = ordinate(phylum_data, "CAP", "bray", ~YoungvsOld)
cap.phy <- plot_ordination(phylum_data, phy.ord, "samples", color="YoungvsOld")
cap.phy <- cap.phy + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(cap.phy)

#### Class 
class.ord = ordinate(class_data, "CAP", "bray", ~YoungvsOld)
cap.class <- plot_ordination(class_data, class.ord , "samples", color="YoungvsOld")
cap.class <- cap.class + geom_point(size = 2) + 
  stat_ellipse(type="norm", level = 0.90, linetype=2) 
print(cap.class)

## 2. Beta Dispersion
### ASV
disp.yo.tax <- betadisper(bc.tax,(as.factor(df.t$YoungvsOld)))
disp.yo.tax
anova(disp.yo.tax) 

plot(disp.yo.tax, ellipse = TRUE, hull=FALSE, main= 'title' ) # 1 sd data ellipse

### Phylum
disp.yo.phy <- betadisper(bc.phy,(as.factor(df.p$YoungvsOld)))
disp.yo.phy
anova(disp.yo.phy) 

plot(disp.yo.phy, ellipse = TRUE, hull=FALSE, main= 'title' ) # 1 sd data ellipse

### Class 
disp.yo.class <- betadisper(bc.class,(as.factor(df.c$YoungvsOld)))
disp.yo.class
anova(disp.yo.class)

## 3. Alpha Diversity
## Alpha diversity measures on rareified but not changed to percentages data
## need to make a different phyloseq object; first split age dataset from variation dataset
comp.asv.fac <- cbind(metadat.red, asv.r)
comp.asv.age <- subset(comp.asv.fac, Age.testing.or.variation.testing.%in% c("Age")) ## split into just age-related roots
comp.asv.var <- subset(comp.asv.fac, Age.testing.or.variation.testing.%in% c("Variation")) ## split to just variation-related roots
comp.asv.age.noperc <- comp.asv.age[,12:1122] ## new ASV
comp.age.met <- comp.asv.age[,1:11] ## new metadata
comp.asv.var.noperc <- comp.asv.var[,12:1122] ## new ASV
comp.asv.var.met <- comp.asv.var[,1:11] ## new metadata

comp.ITS.ASV.age <- otu_table(comp.asv.age.noperc, taxa_are_rows=F)
comp.ITS.Meta.age <- sample_data(comp.age.met,errorIfNULL=TRUE)
comp.ITS.Taxo.age <- tax_table(as.matrix(red.taxon), errorIfNULL=TRUE)
comp.ITS.Phylo.age <- phyloseq(comp.ITS.ASV.age, comp.ITS.Meta.age,comp.ITS.Taxo.age) ## Phyloseq object for just age-selected roots and no percentages
comp.ITS.Phylo.age ## tax will be higher because still includes taxa names from variation data-set

## aggregate at phylum and class levels
phylum_data2 = aggregate_taxa(comp.ITS.Phylo.age, "Phylum")
class_data2 = aggregate_taxa(comp.ITS.Phylo.age, "Class")

phylum_data2@tax_table
class_data2@tax_table

df.t2 = data.frame(sample_data(comp.ITS.Phylo.age))
df.p2 = data.frame(sample_data(phylum_data2))
df.c2 = data.frame(sample_data(class_data2))

## Taxonomic level
rich = estimate_richness(comp.ITS.Phylo.age, measures=c("Observed", "Shannon"))
meta<-as.data.frame(comp.ITS.Phylo.age@sam_data)
richtable<-cbind(meta,rich)

## observed species qq plot 
qqnorm(richtable$Observed, pch = 1, frame = FALSE)
qqline(richtable$Observed, col = "steelblue", lwd = 2)

## shannon diversity qq plot (to check if normally distributed)
qqnorm(richtable$Shannon, pch = 1, frame = FALSE)
qqline(richtable$Shannon, col = "steelblue", lwd = 2)

anova.obs <-aov(rich$Observed~comp.ITS.Phylo.age@sam_data$YoungvsOld)
summary(anova.obs)

anova.shan <-aov(rich$Shannon~comp.ITS.Phylo.age@sam_data$YoungvsOld)
summary(anova.shan)

young <- subset(richtable, YoungvsOld %in% c("Young"))
old <- subset(richtable, YoungvsOld %in% c("Old"))
mean(young$Observed) ## 31.25
mean(young$Shannon) ## 2.28

mean(old$Observed) ##  22.052
mean(old$Shannon) ## 1.96

## Phylum level
rich.phy = estimate_richness(phylum_data2, measures=c("Observed", "Shannon"))
meta.phy <-as.data.frame(phylum_data2@sam_data)
richtable.phy <-cbind(meta.phy,rich.phy)

## observed species qq plot
qqnorm(richtable.phy$Observed, pch = 1, frame = FALSE)
qqline(richtable.phy$Observed, col = "steelblue", lwd = 2)

## shannon diversity qq plot
qqnorm(richtable.phy$Shannon, pch = 1, frame = FALSE)
qqline(richtable.phy$Shannon, col = "steelblue", lwd = 2)

anova.obs.phy <-aov(rich.phy$Observed~phylum_data2@sam_data$YoungvsOld) 
summary(anova.obs.phy)

anova.shan.phy <-aov(rich.phy$Shannon~phylum_data@sam_data$YoungvsOld)
summary(anova.shan.phy) 

## Class level
rich.class = estimate_richness(class_data2, measures=c("Observed", "Shannon"))
meta.class <-as.data.frame(class_data2@sam_data)
richtable.class <-cbind(meta.class,rich.class)

## observed species qq plot
qqnorm(richtable.class$Observed, pch = 1, frame = FALSE)
qqline(richtable.class$Observed, col = "steelblue", lwd = 2)

## shannon diversity qq plot 
qqnorm(richtable.class$Shannon, pch = 1, frame = FALSE)
qqline(richtable.class$Shannon, col = "steelblue", lwd = 2)

anova.obs.class <-aov(rich.class$Observed~class_data2@sam_data$YoungvsOld) 
summary(anova.obs.class)

anova.shan.class <-aov(rich.class$Shannon~class_data2@sam_data$YoungvsOld)
summary(anova.shan.class) 

### 4. Additional composition Information
## Young and old dominant taxa
young.tax <- subset_samples(comp.ITS.Phylo.age, YoungvsOld %in% c("Young"))
old.tax <- subset_samples(comp.ITS.Phylo.age, YoungvsOld %in% c("Old"))

## sort top young
top50.young <- names(sort(taxa_sums(young.tax), decreasing=TRUE))[1:50]
ps.top50.young <- transform_sample_counts(young.tax, function(OTU) OTU/sum(OTU))
ps.top50.young <- prune_taxa(top50.young, ps.top50.young)

## sort bottom young
bottom50.young <- names(sort(taxa_sums(young.tax), decreasing=TRUE))[1:50]
ps.bottom50.young <- transform_sample_counts(young.tax, function(OTU) OTU/sum(OTU))
ps.bottom50.young <- prune_taxa(bottom50.young, ps.bottom50.young)

## what phyla did the top young taxa belong to
topPhyla.young <- tax_table(ps.top50.young)[, "Phylum"]
topPhyla.young.view <- as(topPhyla.young, "vector")
topyoung <- as.data.frame(topPhyla.young.view)

## what phyla did the bottom young taxa belong to
bottomPhyla.young <- tax_table(ps.bottom50.young)[, "Phylum"]
bottomPhyla.young.view <- as(bottomPhyla.young, "vector")
bottomyoung <- as.data.frame(bottomPhyla.young.view)
print(bottomyoung)

## sort top old
top50.old <- names(sort(taxa_sums(old.tax), decreasing=TRUE))[1:50]
ps.top50.old <- transform_sample_counts(young.tax, function(OTU) OTU/sum(OTU))
ps.top50.old <- prune_taxa(top50.old, ps.top50.old)

## what phyla did the top old taxa belong to?
topPhyla.old <- tax_table(ps.top50.old)[, "Phylum"]
topPhyla.old.view <- as(topPhyla.old, "vector")
topold <- as.data.frame(topPhyla.old.view)


## 5. Partial RDA
counts.ja <- as.data.frame(asv.age.perc) ### Just age roots, percent asvs
meta.rda.ja <-as.data.frame(age.met) ## again just age roots
meta.rda.env.ja <- meta.rda.ja[ , -which(names(meta.rda.ja) %in% c("OLD.ID","Age.testing.or.variation.testing.","AgeCat"))]

as.factor(meta.rda.env.ja$Row)
is.factor(meta.rda.env.ja$Row) 
class(meta.rda.env.ja$Row)

counts.hel.ja <-decostand(counts.ja, "hellinger")
decorana(counts.hel.ja)

ja.rda.p <- rda(counts.hel.ja ~ YoungvsOld + Condition(Vine), meta.rda.env.ja) 
summary(ja.rda.p)
set.seed(1)
anova(ja.rda.p)
(ja.r2.root.p <- RsquareAdj(ja.rda.p)$r.squared) ##  
(ja.r2.root.p.adj <-RsquareAdj(ja.rda.p)$adj.r.squared) ## 