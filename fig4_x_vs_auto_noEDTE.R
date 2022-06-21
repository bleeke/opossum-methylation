# BJL
# distribution of methylation at X chromosome vs autosomes
# embryo data


# libraries
library(tidyverse)
library(scales)
library(ggridges)
library(methylKit)



# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/sex_specific"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"


# import data
setwd(dataDir)

# data preparation
# did once

## define samples to work on
#samples_female <- c(
#  paste0(
#    c("oocyte",
#      "1.5_embryo",
#      "2.5_embryo",
#      "3.5_embryo",
#      "4.5_embryo",
#      "5.5_embryo",
#      "6.5_embryo",
#      "7.5_embryo"),
#    "_",
#    "female")
#)
#
#samples_male <- c(
#  paste0(
#    c("sperm",
#      "1.5_embryo",
#      "2.5_embryo",
#      "3.5_embryo",
#      "4.5_embryo",
#      "5.5_embryo",
#      "6.5_embryo",
#      "7.5_embryo"),
#    "_",
#    "male")
#)

## make one samples df
#samples_df <- data.frame(female = samples_female, male = samples_male)
## tidy up
#rm(samples_female, samples_male)
#
## to put loop output in
#df <- data.frame()
#
## read in male and female files for each timepoint
## find sites in common
## calculate percentage methylation
## reorganise and save to dataframe
#for(i in 1:nrow(samples_df)){
#  
#  samples <- as.list(c(as.character(samples_df$female[i]),
#                       as.character(samples_df$male[i])))
#  
#  condition <- c(1:length(samples))
#  
#  files <- as.list(paste0(samples,
#                          "_destranded_pooled_df_filt5.csv")) 
#  
#  # tidy the names to remove 'embryo'
#  samples <- samples %>%
#    str_remove(pattern = "embryo_") %>%
#    as.list()
#  
#  meth_data <- methRead(location = files,
#                        sample.id = samples,
#                        assembly = "monDom5",
#                        pipeline = list(fraction = TRUE,
#                                        chr.col = 1,
#                                        start.col = 2,
#                                        end.col = 3,
#                                        coverage.col = 5,
#                                        strand.col = 4,
#                                        freqC.col = 8),
#                        context = "CpG",
#                        resolution = "base",
#                        treatment = condition,
#                        mincov = 5)
#  
#  # unite to generate methylBase object containing only those sites present in all samples
#  data_united <- methylKit::unite(meth_data)
#  # extract percentage methylation
#  data_pc <- percMethylation(data_united,
#                             rowids = TRUE)
#  data_pc <- data_pc %>%
#    as.data.frame() %>%
#    rownames_to_column("position") %>%
#    pivot_longer(cols = !position,
#                 names_to = c("timepoint", "sex"),
#                 names_sep = "_") %>%
#    separate(position, into = c("chr", "start", "end"),
#             sep = "[.]")
#  
#  df <- rbind(df, data_pc)
#}

setwd(dataDir)
#saveRDS(df, "figure3_embryo_data_noEDTE.RDS")
df <- readRDS("figure3_embryo_data_noEDTE.RDS")
df <- df %>% 
  # a variable for creating facets on the plot
  mutate(timepoint_plot = factor(str_replace(timepoint, "oocyte|sperm", "gamete"),
                                 levels = c("gamete",
                                            "1.5",
                                            "2.5",
                                            "3.5",
                                            "4.5",
                                            "5.5",
                                            "6.5",
                                            "7.5"))) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("sperm",
                                       "oocyte",
                                       "1.5",
                                       "2.5",
                                       "3.5",
                                       "4.5",
                                       "5.5",
                                       "6.5",
                                       "7.5"))) %>%
  # adjust levels for order of plotting
  mutate(sex = factor(sex,
                      levels = c("male", "female"))) %>%
  # remove data from chromosomes MT and pseudoY
  dplyr::filter(chr != "chrMT" & chr != "chrY") %>%
  # categorise data as chrX or autosome - note autosomes includes chrUn
  mutate(chromosome = case_when(chr == "chrX" ~ "chrX",
                                chr != "chrX" ~ "autosome")) %>%
  mutate(chromosome = factor(chromosome,
                             levels = c("autosome",
                                        "chrX")))



# plot
setwd(plotDir)

# fills
fills_plot <- c("#9467BD",
                "#BCBD22",
                "#AEC7E8",
                "#98DF8A",
                "#DBDB8D",
                "#FFBB78",
                "#FF9896",
                "#E377C2",
                "#CD3333")

labs_plot <- c("gamete" = "gametes",
               "1.5" = "E1.5",
               "2.5" = "E2.5",
               "3.5" = "E3.5",
               "4.5" = "E4.5",
               "5.5" = "E5.5",
               "6.5" = "E6.5",
               "7.5" = "E7.5")
               

ridge_plot <- df  %>%
  ggplot(aes(x = value,
             y = 1,
             colour = timepoint,
             fill = timepoint,
             group = sex, 
             linetype = sex)) +
  geom_density_ridges(alpha = 0.2,
                      size = 1.2,
                      show.legend = F)+
  theme_bw()+
  xlab("CpG methylation")+
  ylab("density")+
  facet_grid(chromosome ~ timepoint_plot,
             labeller = labeller(timepoint_plot = labs_plot))+
  scale_fill_manual(values = fills_plot)+
  scale_colour_manual(values = fills_plot)+
  scale_linetype_manual(values = c(3, 1))+
  scale_x_continuous(breaks = c(0, 50, 100),
                     labels = c("0", "0.5", "1"))+
  theme(axis.text = element_text(size = 8),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  # blank layer with toy data to get a linetype legend
  geom_line(data = df[1:2,],
            aes(linetype=sex),
            size = 1.7,
            show.legend = TRUE)+
  guides(linetype = guide_legend(override.aes=list(linetype = c(3,1))),
         colour = FALSE,
         fill = FALSE)

setwd(plotDir)
ggsave(paste0("chr_ridgeplot_embryo.pdf"), plot = ridge_plot, width = 28, height = 5)

# count number of CpGs for legend
count_CpG <- df %>% 
  group_by(sex, timepoint_plot, chromosome) %>%
  summarise(count = n())
setwd(dataDir)
write.table(count_CpG, "Count_CpG_embryo_ridges.txt")