# BJL
# differential methylation
# gametes vs E4.5

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
samples <- as.list(c("sperm", "oocyte", "4.5_embryo"))



# files to import from
files <- as.list(paste0(dataDir,
                            samples,
                            "_destranded_pooled_df_filt5.csv"))


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

pair1 <- c(1, 2)
pair2 <- c(3, 3)



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

# not run
#changing <- apply(counts, 2, function(x) x[c(1,3)] * x[4] / 100)

# plot proportion changing sites 
# labels
labels_plot = c("1->3" = "sp -> 4.5",
                "2->3" = "oo -> 4.5")

#  figure
setwd(plotDir)
counts %>%
  as.data.frame() %>%
  rownames_to_column("diffMeth") %>%
  gather(-diffMeth, key = comparison, value = perc) %>%
  filter(diffMeth != "total") %>%
  mutate(comparison_fact = factor(comparison,
                                  levels = c("1->3",
                                             "2->3"))) %>%
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
ggsave(paste0("diffMeth_", cutoff, "gametes_4.5.pdf"), width = 5, height = 4)



# get n_sites_compared per timepoint comparison for adding to the plot manually

n_sites_compared <- vector()

for(i in 1:length(diffs)){
  n_sites_compared[i] <- nrow(diffs[[i]])/1000
}

