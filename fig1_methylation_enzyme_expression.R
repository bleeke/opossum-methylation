# 20220120
# BJL
# expression of methylation enzymes in opossum embryos
# single cell RNA-seq data Mahadevaiah 2020

# functions 
# from https://chrischizinski.github.io/SNR_R_Group/2016-10-05-Themes_Facets
theme_mine <- function(base_size = 18, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 24),
      strip.text.y = element_text(size = 18),
      axis.text.x = element_text(size=18),
      axis.text.y = element_text(size=24,hjust=1),
      axis.ticks =  element_line(colour = "black"), 
      axis.title.x= element_text(size=16),
      axis.title.y= element_text(size=16,angle=90),
      panel.background = element_blank(), 
      panel.border =element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.margin = unit(1.0, "lines"), 
      plot.background = element_blank(), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1)
    )
}

# libraries
library(tidyverse)
library(DESeq2)
library(scater)
library(SingleCellExperiment)
library(plotrix) # for function std.error

# directories

datadir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/mahadevaiah_2020/"
plotdir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots"
annotdir <- "/camp/lab/turnerj/working/Bryony/annotations"


# data
setwd(datadir)
load("myFeatureCounts.Rdata")
counts <- as.data.frame(myFeatureCounts$counts)  # extract count matrix


# sample metadata
setwd(datadir)
sample_data <- read.delim("sample_identifiers.txt") %>% # sample metadata 
  distinct() %>% # remove dups
  dplyr::filter(inferred.cell.type != "not available") %>% # whole-embryos removed 
  separate(single.cell.identifier, into = c("sample_id"), sep = "_", remove = F) %>%
  #make new column for numerical age
  separate(age, into = c("age_num"), sep = " ") %>%
  # make new column with identifier for age + cell type
  mutate(age_cell_type = paste0(age_num, "_", inferred.cell.type))


# sample conditions to keep
keep <- c("7.5_hypoblast", "7.5_epiblast", "7.5_trophectoderm")

sample_data_keep <- sample_data %>% 
  filter(age_cell_type %in% keep | age_num != 7.5) %>%
  # construct a new column for categorisation when plotting
  mutate(age_cell_type_plot = case_when(age_num != 7.5 ~ age_num,
                                        age_num == 7.5 & inferred.cell.type == "hypoblast" ~ "7.5_hypo",
                                        age_num == 7.5 & inferred.cell.type == "epiblast" ~ "7.5_epi",
                                        age_num == 7.5 & inferred.cell.type == "trophectoderm" ~ "7.5_te"
  ))


# build objects compatible with scater
#---------------------------------------

#build annot data object
cell_ans <- data.frame(sample_id = sample_data$sample_id,
                       sex = sample_data$sex,
                       inferred_cell_type = sample_data$inferred.cell.type,
                       age = sample_data$age_num,
                       individual = sample_data$individual)
rownames(cell_ans) <- sample_data$sample_id

# tidy counts object
names <- colnames(counts)
names <- as.data.frame(names)
names <- names %>% separate(names,into=c("names"), sep = "[.]", remove = T)

colnames(counts) <- names$names

# order cell_ans to match colnames(counts)
cell_ans <- cell_ans[match(colnames(counts), rownames(cell_ans)),]

#build singlecellexperiment object
reads <-SingleCellExperiment(list(counts = as.matrix(counts)),
                             colData = cell_ans,
                             rowData = list(gene_id = rownames(counts)))


# some QC
#---------------------------------------
#
#setwd(annotdir)
#mito_genes <- read.delim("20201023_mondom5_mito_genes.txt", header = F, col.names = "gene_id") # mondom5 mitochondrial genes
#
#per.cell <- perCellQCMetrics(reads, subsets =list(Mito=mito_genes$gene_id))
#colData(reads) <- cbind(colData(reads), per.cell)
#
#setwd(plotdir)
#plotColData(reads, x = "sum", y="detected", colour_by="age") 
#ggsave(filename = "detected_sum_age.png")
#
#
#plotColData(reads, x = "sum", y="detected", colour_by="individual") 
#ggsave(filename = "detected_sum_individual.png")
#
#
#plotColData(reads, x = "sum", y = "subsets_Mito_percent",
#            other_fields = "age", colour_by = "individual") + facet_wrap(~age)
#ggsave("sum_Mitoperc_age.png")
#
#
#per.feat <- perFeatureQCMetrics(reads)
#
#
#plotHighestExprs(reads, exprs_values = "counts")
#ggsave("highest_expressed.png")

# examine expression of genes of interest
#---------------------------------------

#calculate expression
reads <- logNormCounts(reads) # log norm counts
fpkm(reads) <- calculateFPKM(reads, lengths = myFeatureCounts$annotation$Length) # fpkm
assay(reads, "log2fpkm") <- log2((fpkm(reads)+1))

# get fpkm values as dataframe
fpkm_table <- fpkm(reads) %>% as.data.frame() %>% rownames_to_column(var = "gene_id")

# genes of interest
setwd(datadir)
goi <- read_delim("opossum_methylation_enzymes.txt", delim = "\t")

# filter fpkm table for genes of interest
fpkm_table_goi <- fpkm_table %>% 
  filter(gene_id %in% goi$gene_id)


#  plot 
setwd(plotdir)

# colours
cols <- c("gray38", # oocyte
          "gray38", # E1.5
          "gray38", # E2.5
          "gray38", # E3.5
          "gray38", # E4.5
          "gray38", # E5.5
          "gray38", # E6.5
          "#CD3333", # E7.5
          "#f66e6e",
          "#8B2323")

# for splitting lines
lineage <- c("7.5_epi", "7.5_hypo", "7.5_te")

types <- unique(goi$gene_name)

for (i in 1:length(types)){
  type <- fpkm_table_goi %>%
    # extract genes of interest
    filter(gene_id %in% goi$gene_id[goi$gene_name == types[i]]) %>%
    # make long
    gather(-gene_id, key = sample_id, value = fpkm) %>%
    # exclude timepoints/cell types not of interest
    filter(sample_id %in% sample_data_keep$sample_id) %>%
    # join with metadata
    left_join(select(sample_data_keep, age_num, age_cell_type_plot, sample_id), by = "sample_id") %>%
    # calculate log2(fpkm+1)
    mutate(log2_fpkm = log2(fpkm+1)) %>%
    # group by cell type/timepoint and get mean and se
    group_by(gene_id, age_cell_type_plot, age_num) %>%
    summarise(mean = mean(log2_fpkm),
              se = std.error(log2_fpkm))
  
  first_line <- type %>%
    filter(!age_cell_type_plot %in% lineage)
  
  second_line <- type %>%
    filter(age_cell_type_plot %in% c("6.5", "7.5_epi"))
  
  third_line <- type %>%
    filter(age_cell_type_plot %in% c("6.5", "7.5_hypo"))
  
  fourth_line <- type %>%
    filter(age_cell_type_plot %in% c("6.5", "7.5_te"))
  
  ggplot(type,
         aes(x = age_num,
             y = mean,
             colour = age_cell_type_plot))+
    geom_pointrange(aes(ymin = mean - (1.96*se), ymax = mean + (1.96*se)), 
                    #position = position_jitter(width = 0.2, seed = 123)
                    ) +
    geom_line(data=first_line, aes(x = age_num, y = mean, group = gene_id), size = 1)+
    geom_line(data=second_line, aes(x = age_num, y = mean, group = gene_id), colour = "#CD3333", size = 1)+
    geom_line(data=third_line, aes(x = age_num, y = mean, group = gene_id),  colour = "#f66e6e", size = 1)+
    geom_line(data=fourth_line, aes(x = age_num, y = mean, group = gene_id), colour = "#8B2323", size = 1)+
    scale_colour_manual(values = cols)+
    ggtitle(types[i], )+
    ylab("")+
    xlab("")+
    ylim(-1, 12)+
    theme_mine()+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank())
  
  ggsave(paste0(types[i], ".pdf"), width = 6, height = 6, useDingbats=FALSE)
}

