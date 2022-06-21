# 20191128
# BJL
# merged FPKM

# ------------------------------------------------------------------------------
# this script is to calculate FPKM from the featureCounts output of RNA-seq data
# having first pooled the data from biological replicates
# this script filters for fpkm > 1 and excludes pseudogenes 
# before it writes the final object
# ------------------------------------------------------------------------------



# if run interactive
# srun
# ml Anaconda2
# source activate bRy
# R + copy and paste

# if run script
# Rscript this script using 
# BJL submitR script


# libraries

library(Rsubread)
library(DESeq2)
library(tidyverse)
library(lubridate)
library(refGenome)
library(rtracklayer)


# directories

dataDir <- "/camp/home/leekeb//working/Bryony/opossum_adult/allele-specific/data/rna-seq"
plotDir <- "/camp/home/leekeb/working/Bryony/opossum_adult/allele-specific/plots"
annotDir <- "/camp/home/leekeb/working/shared_projects/OOPs/genome"

# load data

setwd(dataDir)
load("myFeatureCounts.Rdata")  # contains featureCounts output 
sample_data <- read.csv("LEE63-sample-sheet.csv")  # sample metadata


# gtf of opossum genes - modified by JZ to match changes made to reference genome - includsion of pseudoY and filling gaps in Rsx locus and to contain both copies of DNMT1
setwd(annotDir)
gtf <- import(file.path(annotDir, "Monodelphis_domestica.monDom5.97_mod_shifted_XY_DNMT1.gtf"))


# featureCounts output
cts <- as.data.frame(myFeatureCounts$counts)  # extract count matrix
cts <- cts %>% select(matches("genome1|genome2|allele_flagged")) # exclude columns for unassigned, conflicting
cts <- rownames_to_column(cts, var = "gene_id") # make rownames a column for gene id



# merge samples by condition
longCts <- gather(cts, key = "bam_name", value = "count", -gene_id) %>%
  separate(bam_name, into = c("limsid", "genome"), sep = "[.]", remove = F) %>%
  separate(limsid, into = "limsid", sep = "_") # makes cts table long, adds and tidies column for limsid and genome while keeping bam_name. FYI throws warning 'Too many values..' but this is OK, it's due to using separate() to split and extract only some fields from 'bam_name'. 

longCts$sex <- ifelse (longCts$bam_name %in% sample_data$bam_name[sample_data$sex == "male"], "male", "female") # adds column for sample sex

# adds column for sample tissue
longCts <- longCts %>% mutate(tissue = case_when(
  (bam_name %in% sample_data$bam_name[sample_data$tissue == "brain"]) ~ "brain",
  (bam_name %in% sample_data$bam_name[sample_data$tissue == "liver"]) ~ "liver",
  (bam_name %in% sample_data$bam_name[sample_data$tissue == "spleen"]) ~ "spleen"
))

longCts <- longCts %>% mutate(condition = paste0(sex, ".", tissue, ".", genome)) # adds column for sample condition

mergedCts <- longCts %>% group_by(gene_id, condition) %>% summarise(merged_count = sum(count)) # creates new table containing the sum of reads per condition

mergedCts <- mergedCts %>% spread(key = condition, value = merged_count) %>% remove_rownames() %>% column_to_rownames('gene_id') %>% as.matrix() # formats to DESseqData compatible


# construct coldata
coldata <- data.frame(condition = colnames(mergedCts)) %>% 
  separate(condition, into = c("sex", "tissue", "genome"), remove = FALSE, sep = "[.]") %>%
  mutate(rownames = condition)
coldata <- column_to_rownames(coldata, var = "rownames") # rowids from condition


# construct DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = mergedCts,
                              colData = coldata,
                              design = ~ condition) 

# add feature lengths as rowData
mcols(dds) <- myFeatureCounts$annotation$Length
names(mcols(dds)) <- "basepairs"

# estimate size factors
dds <- estimateSizeFactors(dds)

# calculate FPKM
fpkms <- fpkm(dds, robust = T)

# exclude pseudogenes
pseudogenes <- unique(gtf$gene_id[str_detect(gtf$gene_biotype, "pseudogene")])

fpkms <- as.data.frame(fpkms) %>% 
  rownames_to_column(var = "gene_id") 

fpkms <- fpkms[!(fpkms$gene_id %in% pseudogenes),]

# exclude columns from allele-specific files
fpkms_allele_flagged <- fpkms %>% select(gene_id, as.vector(coldata$condition[coldata$genome == "allele_flagged"]))


# filter for genes with FPKM >1  per tissue 
fpkms_brain <- select(fpkms_allele_flagged, gene_id, as.vector(coldata$condition[coldata$tissue == "brain" & coldata$genome == "allele_flagged"])) %>%
  filter_all(all_vars(. > 1))

fpkms_liver <- select(fpkms_allele_flagged, gene_id, as.vector(coldata$condition[coldata$tissue == "liver" & coldata$genome == "allele_flagged"])) %>%
  filter_all(all_vars(. > 1))


fpkms_spleen <- select(fpkms_allele_flagged, gene_id, as.vector(coldata$condition[coldata$tissue == "spleen" & coldata$genome == "allele_flagged"])) %>%
  filter_all(all_vars(. > 1))

# save objects
setwd(dataDir)
save(coldata, file = "coldata_merged.R")
save(mergedCts, file = "mergedCts.R")
save(fpkms, file = "merged_fpkms.R")
save(fpkms_allele_flagged, file = "merged_fpkms_allele_flagged.R")
save(fpkms_liver, file = "merged_liver_fpkms.R")    
save(fpkms_brain, file = "merged_brain_fpkms.R")   
save(fpkms_spleen, file = "merged_spleen_fpkms.R")