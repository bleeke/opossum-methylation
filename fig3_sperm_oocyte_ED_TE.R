# BJL
# differential methylation
# annotate diffMeth regions with genomic features
# calculate for all
# make plots for:
# oocyte vs sperm
# ED vs TE

# libraries
library(tidyverse)
library(methylKit)
library(rtracklayer)
library(genomation)

# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots"

# load data
setwd(dataDir)
# list of methylDiff objects for pairwise timecourse comparisons
diffs <- readRDS("diffMeths_20.RDS")

# load annotations 

setwd(annotDir)
# load from gtf to GRanges
genes <- import("mondom5.97_Xshifted_Y_DNMT1.gtf")
# filter out entries for transcripts etc and include pseudoY chromosome CDS 'genes'
genes <- genes[genes$type == "gene" | seqnames(genes) == "chrY"]
# exclude pseudogenes
pseudogenes <- unique(genes$gene_id[str_detect(genes$gene_biotype, "pseudogene")])
genes <- genes[!genes$gene_id %in% pseudogenes]
# exclude duplicated loci
genes_granges <- genes[!paste0(seqnames(genes), ranges(genes)) %>% duplicated,]
# tidy up
rm(genes, pseudogenes)

#intergenic
# set seqlengths of genes_granges
chrom_lengths <- chrom_lengths <- read_delim("genome/chrom_lengths.txt", delim = "\t")
seqlengths(genes_granges) <- chrom_lengths$length
intergenic_granges <- setdiff(as(seqinfo(genes_granges), "GRanges"), genes_granges, ignore.strand = TRUE)
#tidy up
rm(chrom_lengths)

#promoters
promoters_granges <- unique(promoters(genes_granges))

# repeats
repeats <- read_delim("monDom5_UCSC_RepeatMasker_Xshifted_ClassFamily.txt",
                      delim = "\t",
                      col_names = TRUE)
repeats_granges <- GRanges(seqnames = repeats$chrom ,
                           IRanges(start = repeats$start,
                                   end = repeats$end),
                           strand = repeats$strand)

# cgis
cgis_granges <- readRDS("mondom5_CGI_Xshifted.RDS")

# list of all features
features <- list(repeats_granges, genes_granges, intergenic_granges, promoters_granges, cgis_granges)
features_names <- c( "repeatRegions", "geneRegions", "intergenicRegions", "promoterRegions", "cgiRegions")



# annotate hyper and hypo diffMeth regions with genomic features

# to store loop output
features_hyper <- matrix(NA, ncol = length(features), nrow = length(diffs), 
                         dimnames = list(names(diffs), features_names))

features_hypo <- matrix(NA, ncol = length(features), nrow = length(diffs), 
                        dimnames = list(names(diffs), features_names))
features_same <- matrix(NA, ncol = length(features), nrow = length(diffs), 
                        dimnames = list(names(diffs), features_names))
features_all <- matrix(NA, ncol = length(features), nrow = length(diffs), 
                        dimnames = list(names(diffs), features_names))

n_sites_compared <- list()
n_sites_hyper <- list()
n_sites_hypo <- list()

for(i in 1:length(diffs)){
  pair_diffMeth <- diffs[[i]]
  n_sites_compared[i] <- nrow(pair_diffMeth)
  pair_hyper <- getMethylDiff(pair_diffMeth, type = "hyper", difference = 20)
  n_sites_hyper[i] <- nrow(pair_hyper)
  pair_hypo <- getMethylDiff(pair_diffMeth, type = "hypo", difference = 20)
  n_sites_hypo[i] <- nrow(pair_hypo)
  pair_same <- pair_diffMeth[pair_diffMeth$meth.diff >=-20 & pair_diffMeth$meth.diff <=20 & pair_diffMeth$qvalue < 0.01 | pair_diffMeth$qvalue >= 0.01]
  for(k in 1:length(features)){
    hyper <- annotateWithFeature(features[[k]],
                                 as(pair_hyper, "GRanges"))
                                    
    hypo <- annotateWithFeature(features[[k]],
                                as(pair_hypo, "GRanges"))
    
    same <- annotateWithFeature(features[[k]], as(pair_same, "GRanges")
                                )
    all <-  annotateWithFeature(features[[k]], as(pair_diffMeth, "GRanges"))
    features_hyper[i,k] <- hyper@perc.of.OlapFeat
    features_hypo[i,k] <- hypo@perc.of.OlapFeat
    features_same[i,k] <- same@perc.of.OlapFeat
    features_all[i,k] <- all@perc.of.OlapFeat
  }
}

names(n_sites_compared) <- names(diffs)
names(n_sites_hyper) <- names(diffs)
names(n_sites_hypo) <- names(diffs)

# save data
setwd(dataDir)
saveRDS(features_hyper, "diffMethCpGs_hyper_annot_genomic_feature.RDS")
saveRDS(features_hypo, "diffMethCpGs_hypo_annot_genomic_feature.RDS")
saveRDS(features_same, "diffMethCpGs_same_annot_genomic_feature.RDS")s()
saveRDS(features_all, "diffMethCpGs_all_annot_genomic_feature.RDS")
saveRDS(n_sites_compared, "diffMethCpGs_n_sites_compared.RDS")
saveRDS(n_sites_hyper, "diffMethCpGs_n_sites_hyper.RDS")
saveRDS(n_sites_hypo, "diffMethCpGs_n_sites_hypo.RDS")

# make plots

setwd(plotDir)


features_data <- list(as.data.frame(features_hyper),
                      as.data.frame(features_hypo),
                       as.data.frame(features_all))
names(features_data) <- c("hyper", "hypo", "all")


plot_data <- features_data %>%
  map(rownames_to_column, var = "comparison") %>%
  bind_rows(.id = "type") %>%
  pivot_longer(!c(type, comparison), names_to = "region", values_to = "perc") %>%
  mutate(comparison_fact = factor(comparison,
                                  levels = c("1->2",
                                             "1->3",
                                             "2->3",
                                             "3->4",
                                             "4->5",
                                             "5->6",
                                             "6->7",
                                             "7->8",
                                             "8->9",
                                             "7->10",
                                             "7->11",
                                             "10->11",
                                             "10->12",
                                             "10->13",
                                             "10->14")))

labels_plot = c("1->2" = "sp -> oo",
                "1->3" = "sp -> 1.5",
                "2->3" = "oo -> 1.5",
                "3->4" = "1.5 -> 2.5",
                "4->5" = "2.5 -> 3.5",
                "5->6" = "3.5 -> 4.5",
                "6->7" = "4.5 -> 5.5",
                "7->8" = "5.5 -> 6.5",
                "8->9" = "6.5 -> 7.5",
                "7->10" = "5.5 -> 7.5_ED",
                "7->11" = "5.5 -> 7.5_TE",
                "10->11" =  "7.5_ED -> 7.5_TE",
                "10->12" = "7.5_ED -> brain",
                "10->13" = "7.5_ED -> liver",
                "10->14" = "7.5_ED -> spleen")


all_plot <- plot_data %>%
  ggplot(aes(x = comparison_fact,
             y = perc,
             group = type,
             fill = type))+
  geom_col(position = "dodge")+
  facet_grid(region~.)+
  scale_x_discrete(labels = labels_plot)
ggsave("test.pdf", all_plot, width = 28)
  

# try a normalisation to all detected

features_hyper_norm <- features_hyper/features_all
features_hypo_norm <- features_hypo/features_all

features_data_norm <- list(as.data.frame(features_hyper_norm),
                      as.data.frame(features_hypo_norm))
names(features_data_norm) <- c("hyper", "hypo")



plot_data <- features_data_norm %>%
  map(rownames_to_column, var = "comparison") %>%
  bind_rows(.id = "type") %>%
  pivot_longer(!c(type, comparison), names_to = "region", values_to = "perc") %>%
  mutate(comparison_fact = factor(comparison,
                                  levels = c("1->2",
                                             "1->3",
                                             "2->3",
                                             "3->4",
                                             "4->5",
                                             "5->6",
                                             "6->7",
                                             "7->8",
                                             "8->9",
                                             "7->10",
                                             "7->11",
                                             "10->11",
                                             "10->12",
                                             "10->13",
                                             "10->14")))

labels_plot = c("1->2" = "sp -> oo",
                "1->3" = "sp -> 1.5",
                "2->3" = "oo -> 1.5",
                "3->4" = "1.5 -> 2.5",
                "4->5" = "2.5 -> 3.5",
                "5->6" = "3.5 -> 4.5",
                "6->7" = "4.5 -> 5.5",
                "7->8" = "5.5 -> 6.5",
                "8->9" = "6.5 -> 7.5",
                "7->10" = "5.5 -> 7.5_ED",
                "7->11" = "5.5 -> 7.5_TE",
                "10->11" =  "7.5_ED -> 7.5_TE",
                "10->12" = "7.5_ED -> brain",
                "10->13" = "7.5_ED -> liver",
                "10->14" = "7.5_ED -> spleen")


all_plot <- plot_data %>%
  ggplot(aes(x = comparison_fact,
             y = perc,
             group = type,
             fill = type))+
  geom_col(position = "dodge")+
  facet_grid(region~.)+
  scale_x_discrete(labels = labels_plot)
ggsave("test_norm.pdf", all_plot, width = 28)
