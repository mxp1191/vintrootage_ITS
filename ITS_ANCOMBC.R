### Meredith Persico: mxp1191@psu.edu
### Code for: "The age of absorptive roots impacts root-adjacent microbial composition in grapevines"
### This code covers differential abundance analysis using ancombc for ITS data

## Set directory & check files
setwd("/Users/mjp23/Desktop/RootAgePostDada2")

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
library(ANCOMBC)
library(microbiome)
library(scales)
library(pheatmap)

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
asv.filt 

## Calculate reads per sample
rowSums(asv.filt) ## several have no samples

## Determine minimum number of reads in your data set
min(rowSums(asv.filt))

## Cut samples with too few reads
n2 <- names(which(rowSums(asv.filt) > 1000)) 
asv.re<-asv.filt[which(rownames(asv.filt) %in% n2),]
min(rowSums(asv.re)) ## min now 1158 reads

metadat.red <-metadat[which(rownames(metadat) %in% rownames(asv.re)),] ## match metadata with removed samples

## Split up age and variation roots from metadata and asv.perc files
ancom.asv.fac <- cbind(metadat.red, asv.re)
ancom.asv.age <- subset(ancom.asv.fac, Age.testing.or.variation.testing.%in% c("Age")) ## split into just age-related roots
ancom.asv.var <- subset(ancom.asv.fac, Age.testing.or.variation.testing.%in% c("Variation")) ## split to just variation-related roots

ancom.asv.age.noperc <- ancom.asv.age[,12:1122] ## new ASV
ancom.age.met <- ancom.asv.age[,1:11] ## new metadata

## removing the residual zero asv's from when age and var datasets were combined
test <- t(ancom.asv.age.noperc)
rowSums(test)
no_zeroes <- test[rowSums(test)>0,] 
no.zero.taxa <- red.taxon[which(rownames(red.taxon) %in% rownames(no_zeroes)),] ## remove zero taxa
asv.nz <- t(no_zeroes) ## asv
meta.nzes <- metadat.red[which(rownames(metadat.red) %in% rownames(asv.nz)),]

ancom.asv.var.noperc <- ancom.asv.var[,12:1122] ## new ASV
ancom.asv.var.met <- ancom.asv.var[,1:11] ## new metadata

ancom.age.met$Vine <- as.factor(ancom.age.met$Vine)
ancom.age.met$YoungvsOld <-as.factor(ancom.age.met$YoungvsOld)


#### DID NOT NORMALIZE READS COUNTS OR REMOVE LOW ABUNDANT TAXA ####

## create phyloseq object of non-normalized data
ancom.ITS.ASV.age <- otu_table(ancom.asv.age.noperc, taxa_are_rows=F)
ancom.ITS.Meta.age <- sample_data(ancom.age.met,errorIfNULL=TRUE)
ancom.ITS.Taxo.age <- tax_table(as.matrix(red.taxon), errorIfNULL=TRUE)
ancom.ITS.Phylo.age <- phyloseq(ancom.ITS.ASV.age, ancom.ITS.Meta.age,ancom.ITS.Taxo.age) 
ancom.ITS.Phylo.age


#### phyloseq object of non-normalized data and no zeroes
nz.ITS.ASV.age <- otu_table(asv.nz, taxa_are_rows=F)
nz.ITS.Meta.age <- sample_data(meta.nzes,errorIfNULL=TRUE)
nz.ITS.Taxo.age <- tax_table(as.matrix(no.zero.taxa, errorIfNULL=TRUE))
nz.ITS.Phylo.age <- phyloseq(nz.ITS.ASV.age, nz.ITS.Meta.age, nz.ITS.Taxo.age) 
nz.ITS.Phylo.age

view(nz.ITS.Taxo.age)

########### ANCOMBC: ASV LEVEL

### Step 1: zero_cut threshold is ".99" which means it includes all samples
### and even taxa that are rare.

out = ancombc(phyloseq = nz.ITS.Phylo.age , formula = "YoungvsOld", 
              p_adj_method = "holm", zero_cut = 0.99, lib_cut = 1000, 
              group = "YoungvsOld", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

res = out$res
res_global = out$res_global
tab_p = res$p_val ## p value
tab_q = res$q ## adjusted p values


tab_diff = res$diff_abn
as.data.frame(tab_diff)
count(tab_diff, tab_diff$YoungvsOld == "TRUE")
true <- subset(tab_diff, tab_diff$YoungvsOldYoung == "TRUE")
as.data.frame(true)
rownames(tab_diff) 

acomb.true.tax<-red.taxon[which(rownames(red.taxon) %in% rownames(true)), 2:7]
as.data.frame(acomb.true.tax)
tab_w = res$W

out$feature_table
out$samp_frac
out$zero_ind
out$samp_frac
out$resid
out$delta_wls


### Step 2: zero_cut threshold is now "0.75" which means it includes only
### taxa appearing in > 25% of samples

out.25 = ancombc(phyloseq = nz.ITS.Phylo.age , formula = "YoungvsOld", 
                 p_adj_method = "holm", zero_cut = 0.75, lib_cut = 1000, 
                 group = "YoungvsOld", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                 max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

res.25 = out.25$res
tab_p.25 = res.25$p_val ## p value
tab_q.25 = res.25$q ## adjusted p values


tab_diff.25 = res.25$diff_abn
as.data.frame(tab_diff.25)
count(tab_diff.25, tab_diff.25$YoungvsOld == "TRUE")
true.25 <- subset(tab_diff.25, tab_diff.25$YoungvsOldYoung == "TRUE")
as.data.frame(true.25)
rownames(tab_diff.25) 

acomb.true.tax.25<-red.taxon[which(rownames(red.taxon) %in% rownames(true.25)), 2:7]
as.data.frame(acomb.true.tax.25)
tab_w = res$W

out.25$feature_table
out.25$samp_frac
out.25$zero_ind
out.25$samp_frac
out.25$resid
out.25$delta_wls

# Phylum              Class          Order                 Family             Genus        Species
# F_asv13 p__Basidiomycota c__Tremellomycetes o__Tremellales f__Trimorphomycetaceae      g__Saitozyma   s__podzolica
# F_asv16    p__Ascomycota c__Sordariomycetes o__Hypocreales        f__Hypocreaceae    g__Trichoderma s__aggressivum
# F_asv18    p__Ascomycota c__Sordariomycetes o__Hypocreales      f__Bionectriaceae   g__Clonostachys       s__rosea
# F_asv54    p__Ascomycota c__Sordariomycetes o__Sordariales       f__Chaetomiaceae g__Dichotomopilus s__subfunicola
# F_asv61    p__Ascomycota c__Sordariomycetes o__Sordariales       f__Chaetomiaceae     g__Chaetomium s__homopilatum
# F_asv72 p__Glomeromycota  c__Glomeromycetes   Unclassified           Unclassified      Unclassified   Unclassified

#### Calculate which samples have more/less of the different taxa
view(out.25$feature_table)
feat.flip <- t(out.25$feature_table) ## flip feature table
meta.filt <-nz.ITS.Phylo.age@sam_data[which(rownames(nz.ITS.Phylo.age@sam_data) %in% rownames(feat.flip)),]
meta.filt ## just looking at what comes up

taxa.tab <- cbind(feat.flip, meta.filt) ## merge metadata and feature table


young.tab.tax <- subset(taxa.tab, YoungvsOld %in% c("Young"))
old.tab.tax <- subset(taxa.tab, YoungvsOld %in% c("Old"))

mean(young.tab.tax$F_asv13) ## 617.7857
mean(old.tab.tax$F_asv13) ## 171.1875

mean(young.tab.tax$F_asv16) ## 926.9286
mean(old.tab.tax$F_asv16) ## 80.8125

mean(young.tab.tax$F_asv18) ## 682.9286
mean(old.tab.tax$F_asv18)## 103.5

mean(young.tab.tax$F_asv54) ## 213.4286
mean(old.tab.tax$F_asv54) ## 43.1875

mean(young.tab.tax$F_asv61) ## 91.85714
mean(old.tab.tax$F_asv61) ##  49.1875

mean(young.tab.tax$F_asv72) ## 125.1429
mean(old.tab.tax$F_asv72) ## 7.5

###### REORGANIZE DATA FOR HEATMAP
## need to average the samples by young vs. old 
taxa.tab2 <- taxa.tab[,-c(21:30)]
avgs <- aggregate(. ~ YoungvsOld, data=taxa.tab2, FUN= mean)
str(avgs)
row.names(avgs) <- avgs$YoungvsOld
avgs <- avgs[,-1]

## make dataframe numeric
avgs2 <- avgs %>% mutate_at(1:20, as.numeric)
str(avgs2) ## sums2 now numeric
avgs2 <- as.matrix(avgs2)

avgs2 <- t(avgs2)

## change ASVs to genus/species
heat.tax <-red.taxon[which(rownames(red.taxon) %in% rownames(avgs2)), 6:7]
heat.tax$Genus <- gsub('g__', '', as.character(heat.tax$Genus))
heat.tax$Species <- gsub('s__','', as.character(heat.tax$Species))
heat.tax$con <- paste(heat.tax$Genus, heat.tax$Species)

# ##rename the unclassified, asvs with same species label, and indicate * for sig diff
heat.tax[3, 3] = "Fusarium acutatum 1"
heat.tax[5, 3] = "Fusarium acutatum 2"
heat.tax[7, 3] = "Saitozyma podzolica *"
heat.tax[9, 3] = "Trichoderma aggressivum *"
heat.tax[11, 3] = "Clonostachys rosea 1 *"
heat.tax[14, 3] = "Clonostachys rosea 2"
heat.tax[16, 3] = "Unclassified (Sordariomycetes)"
heat.tax[17, 3] = "Dichotomopilus subfunicola *"
heat.tax[19,3] = "Chaetomium homopilatum *"
heat.tax[20,3] = "Unclassified (Glomeromycetes) *"


row.names(avgs2) <- heat.tax$con
view(avgs2)
rownames(avgs2)

col_order <- c("Young", "Old")
avgs2 <- avgs2[, col_order] ## re order young vs. old

## order rows alphabetical
rr <- avgs2[order(row.names(avgs2)), ]

## heatmap (> 25% abundance)
pheatmap(rr, 
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = hcl.colors(100, "Tofino"))


## At class level (99 cutoff)
class_data = aggregate_taxa(nz.ITS.Phylo.age, "Class")
# The taxonomy table
tax_mat = as(tax_table(class_data), "matrix")
# Run ancombc function
out4.1 = ancombc(phyloseq = class_data, formula = "YoungvsOld",
                 p_adj_method = "holm", zero_cut = 0.99, lib_cut = 0,
                 group = "YoungvsOld", struc_zero = TRUE, neg_lb = FALSE,
                 tol = 1e-5, max_iter = 100, conserve = TRUE,
                 alpha = 0.05, global = FALSE)
res4.1 = out4.1$res
tab_diff4.1 = res4.1$diff_abn
as.data.frame(tab_diff4.1)
count(tab_diff4.1, tab_diff4.1$YoungvsOld == "TRUE") ## 5 different

## at 25% prevalence
# Run ancombc function
out4.2 = ancombc(phyloseq = class_data, formula = "YoungvsOld",
                 p_adj_method = "holm", zero_cut = 0.75, lib_cut = 0,
                 group = "YoungvsOld", struc_zero = TRUE, neg_lb = FALSE,
                 tol = 1e-5, max_iter = 100, conserve = TRUE,
                 alpha = 0.05, global = FALSE)
res4.2 = out4.2$res
tab_diff4.2 = res4.2$diff_abn
as.data.frame(tab_diff4.2)
count(tab_diff4.2, tab_diff4.2$YoungvsOld == "TRUE") 

## at 50% prevalence
# Run ancombc function
out4.3 = ancombc(phyloseq = class_data, formula = "YoungvsOld",
                 p_adj_method = "holm", zero_cut = 0.50, lib_cut = 0,
                 group = "YoungvsOld", struc_zero = TRUE, neg_lb = FALSE,
                 tol = 1e-5, max_iter = 100, conserve = TRUE,
                 alpha = 0.05, global = FALSE)
res4.3 = out4.3$res
tab_diff4.3 = res4.3$diff_abn
as.data.frame(tab_diff4.3)
count(tab_diff4.3, tab_diff4.3$YoungvsOld == "TRUE") 