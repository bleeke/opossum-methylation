# BJL
# differential methylation

# libraries
library(tidyverse)
library(methylKit)

# directories

dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/"
adultDir  <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out/"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"


# variables

cutoff <- 20

# load data
# make inputs for methRead function
# sample names
samples_emb <- as.list(c("sperm", "oocyte", "1.5_embryo",
                         "2.5_embryo", "3.5_embryo", "4.5_embryo",
                  "5.5_embryo", "6.5_embryo", "7.5_embryo",
                  "7.5_ED", "7.5_TE"))

samples_ad <- as.list(c("brain", "liver", "spleen")) 

samples <- c(samples_emb, samples_ad)

# files to import from
files_emb <- as.list(paste0(dataDir,
                            samples_emb,
                       "_destranded_pooled_df_filt5.csv"))

files_ad <- as.list(paste0(adultDir,
                           samples_ad,
                           "_destranded_pooled_df_filt5.csv"))

files <- (c(files_emb, files_ad))

# 'treatment' vector
condition <- c(1:length(files))

# read data
meth_dat <- methRead(location = files, 
                     sample.id = samples, 
                     assembly = "MonDom5", 
                     pipeline = list(fraction = TRUE, 
                                     chr.col = 1, 
                                     start.col = 2, 
                                     end.col = 3, 
                                     coverage.col = 5, 
                                     strand.col = 4, 
                                     freqC.col = 8
                                     ), 
                     context = "CpG", 
                     resolution = "base", 
                     treatment = condition, 
                     mincov = 5)

# differential methylation
# unite to keep only rows present in both

# previous comparisons
#pair1 <- c(1, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 10, 9, 9, 9)
#pair2 <- c(2, 3, 3, 4, 5, 6 ,7, 8, 9, 10, 11, 11, 12, 13, 14)


pair1 <- c(1, 1, 2, 3, 4, 5, 6, 7, 8, 7, 7, 10, 10, 10, 10)
pair2 <- c(2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 12, 13, 14)



counts <- matrix(NA, ncol = length(pair1), nrow = 4, 
                 dimnames = list(c("hypo", "same", "hyper", "total"), 
                                 paste(pair1, pair2, sep = "->")
                 )
)
diffs <- list()

for(i in 1:length(pair1)) {
  a <- pair1[i]
  b <- pair2[i]
  pair_list <- new("methylRawList", 
                   meth_dat@.Data[c(a, b)], 
                   treatment = c(a, b)
  )
  pair_united <- methylKit::unite(pair_list)
  pair_diffMeth <- calculateDiffMeth(pair_united)
  pair_prop <- diffMethPerChr(pair_diffMeth, plot = F, meth.cutoff = cutoff)
  pair_counts <- data.frame(
    hypo = pair_prop$diffMeth.all$percentage.of.hypomethylated, 
    hyper = pair_prop$diffMeth.all$percentage.of.hypermethylated
  )
  pair_counts$same <- 100 - pair_counts$hypo - pair_counts$hyper
  
  counts[, i] <- c(unlist(pair_counts)[c(1, 3, 2)], nrow(pair_united))
  diffs[[i]] <- pair_diffMeth
}

names(diffs) <- paste(pair1, pair2, sep = "->")

setwd(dataDir)
#saveRDS(counts, file = paste0("diffMethCounts_", cutoff, ".RDS"))
#saveRDS(diffs, file = paste0("diffMeths_", cutoff, ".RDS"))


changing <- apply(counts, 2, function(x) x[c(1,3)] * x[4] / 100)

# plot proportion changing sites 

# labels


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


# filter adult data out

split_adult <- c("10->12",
                 "10->13",
                 "10->14")

# main figure
setwd(plotDir)
counts %>%
  as.data.frame() %>%
  rownames_to_column("diffMeth") %>%
  gather(-diffMeth, key = comparison, value = perc) %>%
  filter(diffMeth != "total") %>%
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
                                             "10->14"))) %>%
  filter(!comparison_fact %in% split_adult) %>%
  mutate(diffMeth = factor(diffMeth, levels = c("hyper", "same", "hypo")))%>%
  filter(diffMeth != "same")%>%
  ggplot(aes(x = comparison_fact,
             y = perc,
             fill = diffMeth,
             colour = diffMeth))+
  scale_fill_manual(values = c("#C02324", "#3A5FCD"),
                    labels = c("higher","lower"),
                    name = "methylation\ndifference")+
  scale_colour_manual(values = c("#861819", "#22397b"),
                      guide = FALSE)+
  scale_x_discrete(labels = labels_plot)+
  ylab("CpGs (%)")+
  xlab("timepoint comparison")+
  geom_col(size = 1.2,
           alpha = 0.6)+
  theme_bw()+
  ylim(0,10)+
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
ggsave(paste0("diffMeth_", cutoff, ".pdf"), width = 12, height = 4)


# adult figure
setwd(plotDir)
counts %>%
  as.data.frame() %>%
  rownames_to_column("diffMeth") %>%
  gather(-diffMeth, key = comparison, value = perc) %>%
  filter(!diffMeth %in% c("total", "same")) %>%
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
                                             "10->14"))) %>%
  filter(comparison_fact %in% split_adult) %>%
  mutate(diffMeth = factor(diffMeth, levels = c("hyper", "same", "hypo")))%>%
  
  ggplot(aes(x = comparison_fact,
             y = perc,
             fill = diffMeth,
             colour = diffMeth))+
  scale_fill_manual(values = c("#C02324", "#3A5FCD"),
                    labels = c("higher","lower"),
                    name = "methylation\ndifference")+
  scale_colour_manual(values = c("#861819", "#22397b"),
                      guide = FALSE)+
  scale_x_discrete(labels = labels_plot)+
  ylab("CpGs (%)")+
  xlab("timepoint comparison")+
  geom_col(size = 1.2,
           alpha = 0.6)+
  theme_bw()+
  ylim(0,10)+
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
ggsave(paste0("diffMeth_BrainLiverSpleen", cutoff, ".pdf"), width = 6, height = 4)


# get n_sites_compared per timepoint comparison for adding to the plot manually

n_sites_compared <- vector()

for(i in 1:length(diffs)){
  n_sites_compared[i] <- nrow(diffs[[i]])/1000
}

