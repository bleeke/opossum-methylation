# BJL
# library CpG coverage plot
# Figure 1S

# libraries
library(tidyverse)
library(GenomicRanges)

# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/"
adultDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out/"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots"
outDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/bs_stats"


# read in embryo data
setwd(dataDir)
files <- list.files(pattern = "destranded_pooled_df.csv")
files <- set_names(files, str_remove(files, pattern = "_destranded_pooled_df.csv"))

embryo_data <- files %>%
  map_dfr(read_delim,
          delim = "\t",
          .id = "file")

# read in adult data
setwd(adultDir)
files <- list.files(pattern = "destranded_pooled_df.csv")
files <- set_names(files, str_remove(files, pattern = "_destranded_pooled_df.csv"))

adult_data <- files %>%
  map_dfr(read_delim,
          delim = "\t",
          .id = "file")

data <- rbind(embryo_data, adult_data)

# count number of CpG sites captured at different coverage thresholds
coverages <- data.frame()

for(i in 1:5){
  tmp <- data %>%
    filter(coverage1 > i-1) %>%
    group_by(file) %>%
    summarise(n = n()) %>%
    select(file, n) %>%
    mutate(cov = i)
  coverages <- rbind(coverages, tmp)
}


# get file with all mondom5 CpGs
setwd(annotDir)
all_cpgs <- readRDS("all_mondom_cpgs.RDS")
count_all_cpgs <- length(all_cpgs) # number of 'destranded' CpGs in the monDom5 genome

# tidy up 
rm(all_cpgs, embryo_data, adult_data) 


# calculcate percentage of all mondom5 CpGs captured at each coverage threshold
coverages <- coverages %>%
  mutate(perc_cover = (n/count_all_cpgs)*100)

# write coverage values to file
setwd(outDir)
write_delim(coverages, "coverages_pooled_libraries.txt", delim = "\t")

# plot
setwd(plotDir)


# labels for plots
timepoint_labels <- c("sperm" = "sperm",
                      "oocyte" = "oocyte",
                      "1.5_embryo" = "E1.5",
                      "2.5_embryo" = "E2.5",
                      "3.5_embryo" = "E3.5",
                      "4.5_embryo" = "E4.5",
                      "5.5_embryo" = "E5.5",
                      "6.5_embryo" = "E6.5",
                      "7.5_embryo" = "E7.5",
                      "7.5_TE" = "E7.5 TE",
                      "7.5_ED" = "E7.5 ED",
                      "brain" = "brain",
                      "liver" = "liver",
                      "spleen" = "spleen")
# labels
coverage_labels <- c(">0", ">1", ">2", ">3", ">4")

plot_cols <- c("#9467BD",
          "#BCBD22",
          "#AEC7E8",
          "#98DF8A",
          "#DBDB8D",
          "#FFBB78",
          "#FF9896",
          "#E377C2",
          "#CD3333",
          "#f66e6e", # ED
          "#8B2323", # TE
          "#27408B",
          "#3A5FCD", 
          "#436EEE")

plot_data <- coverages %>%
  mutate(sample = factor(file, 
                         levels = c("sperm",
                         "oocyte",
                         "1.5_embryo",
                         "2.5_embryo",
                         "3.5_embryo",
                         "4.5_embryo",
                         "5.5_embryo",
                         "6.5_embryo",
                         "7.5_embryo",
                         "7.5_ED",
                         "7.5_TE",
                         "brain",
                         "liver",
                         "spleen")))

cov_plot <- plot_data %>%
  ggplot(aes(x = as.factor(coverage_level),
             y = perc_cover,
             colour = sample))+
  geom_line(aes(group = sample),
            size = 2)+
  scale_colour_manual(values = plot_cols)+
  labs(x = "depth of coverage",
       y = "CpGs covered (%)")+
  theme_bw()+
  ylim(0,100)+
  facet_grid(.~sample,
             labeller = labeller(sample = timepoint_labels))+
  theme(legend.position = "none",
        axis.text = element_text(size = 26),
        axis.title = element_text(size = 30),
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 22))
cov_plot
ggsave(file = paste0("opossum_embryo_BSseq_coverages.pdf"), width = 26, height = 4.1)

  





