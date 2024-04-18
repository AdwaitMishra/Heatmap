#Beta Diversity heatmap
BiocManager::install("phyloseq")
install.packages(c("RColorBrewer", "patchwork"))
library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
setwd("C:/Users/anant/OneDrive/Desktop") 
merged_metagenomes <- import_biom("sycamore_cuatroc.biom")

class(merged_metagenomes)
merged_metagenomes@tax_table@.Data <-
substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

View(merged_metagenomes@otu_table@.Data)
merged_metagenomes_bac <- subset_taxa(merged_metagenomes, Kingdom %in% c("Bacteria", "Fungi", "Viruses"))
ref_data <- merged_metagenomes_bac@otu_table@.Data

colnames(merged_metagenomes_bac@otu_table@.Data) <- gsub("_.*", "", colnames(merged_metagenomes_bac@otu_table@.Data))
# c ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48")
merged_metagenomes_bac@sam_data@row.names <- colnames(merged_metagenomes_bac@otu_table@.Data)
# c ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48")
sample_sums(merged_metagenomes_bac)
summary(merged_metagenomes_bac@otu_table@.Data)
summary(merged_metagenomes_bac@tax_table@.Data== "")
head(merged_metagenomes_bac@otu_table@.Data)

percentages <- transform_sample_counts(merged_metagenomes_bac, function(x) x*100 / sum(x) )
percent_1 <- as.data.frame(percentages@otu_table@.Data)

library(dplyr)
z <- rownames(filter_all(percent_1, any_vars(. > 1)))
percentages@otu_table@.Data <-percentages@otu_table@.Data[z, ]
percent_1 <- as.data.frame(percentages@otu_table@.Data)
percent_1$Genus<- merged_metagenomes_bac@tax_table@.Data[rownames(percent_1), 6]
percent_1$Species<- merged_metagenomes_bac@tax_table@.Data[rownames(percent_1), 7]
for (i in 1:nrow(percent_1)) {
  percent_1$name[i] <-paste(percent_1$Genus[i],percent_1$Species[i], sep = " ")}
percent_1 <- subset(percent_1, name!=" ")


rownames(percent_1) <- percent_1$name
percent_1 <- percent_1[, -c(41,42,40)]
percent_1 <- t(percent_1)
percent_syc <- as.data.frame(percent_1)
pheatmap(log(percent_1+0.01), cluster_cols = F, cluster_rows = F, display_numbers =   T, fontsize_number = 3)
