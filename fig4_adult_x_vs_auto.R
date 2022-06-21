# BJL
# ridgeplots comparing autosome and X chromosome methylation level in opossum adult tissues

# libraries
library(tidyverse)
library(scales)
library(ggridges)
library(methylKit)

# directories
dataDir  <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out/sex_specific"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"

# import data
setwd(dataDir)


# data preparation
# did once

#tissues <- c("brain",
#             "liver",
#             "spleen")
#
## to put loop output in
#df <- data.frame()
#
## read in male and female files for each tissue
## find sites in common
## calculate percentage methylation
## reorganise and save to dataframe
#for(i in 1:length(tissues)){
#  
#samples <- as.list(paste0(c("female","male"),
#                          "_",
#                          tissues[i]))
#                     
#condition <- c(1:length(samples))
#
#files <- as.list(paste0(samples,
#                        "_destranded_pooled_df_filt5.csv")) 
#
#adult_data <- methRead(location = files,
#                              sample.id = samples,
#                              assembly = "monDom5",
#                              pipeline = list(fraction = TRUE,
#                                              chr.col = 1,
#                                              start.col = 2,
#                                              end.col = 3,
#                                              coverage.col = 5,
#                                              strand.col = 4,
#                                              freqC.col = 8),
#                              context = "CpG",
#                              resolution = "base",
#                              treatment = condition,
#                              mincov = 5)
#
## unite to generate methylBase object containing only those sites present in all samples
#adult_united <- methylKit::unite(adult_data)
## extract percentage methylation
#adult_pc <- percMethylation(adult_united,
#                            rowids = TRUE)
#adult_pc <- adult_pc %>%
#  as.data.frame() %>%
#  rownames_to_column("position") %>%
#  pivot_longer(cols = !position,
#               names_to = c("sex", "tissue"),
#               names_sep = "_") %>%
#  separate(position, into = c("chr", "start", "end"),
#           sep = "[.]")
#
#df <- rbind(df, adult_pc)
#}
#
##setwd(dataDir)
##saveRDS(df, "figure3_adult_data.RDS")

df <- readRDS("figure3_adult_data.RDS")
df <- df %>%   
  # remove data from pseudoY and mitochondrial chromosomes
  dplyr::filter(chr != "chrMT" & chr != "chrY") %>%
  # categorise data as chrX or autosome - note autosomes includes chrUn
  mutate(chromosome = case_when(chr == "chrX" ~ "chrX",
                                chr != "chrX" ~ "autosome")) %>%
  mutate(chromosome = factor(chromosome,
                             levels = c("autosome",
                                        "chrX")))


# plot
setwd(plotDir)

# fill colour for plot
plot_fills <- c("royalblue4",
                "royalblue3", 
                "royalblue2")

# make  ridgeplot brain - main figure
ridge_plot <- df %>%
  filter(tissue == "brain") %>%
  ggplot(aes(x = value,
             y = sex,
             fill = tissue,
             colour = tissue,
             group = sex, 
             linetype = sex))+
  # make ridgeplots quite see-through
  geom_density_ridges(size = 1.2,
                      alpha=0.2,
                      show.legend = F)+
  theme_bw()+
  xlab("CpG methylation")+
  ylab("density")+
  facet_grid(chromosome ~ tissue)+
  scale_fill_manual(values = plot_fills)+
  scale_linetype_manual(values = c(1, 3))+
  scale_colour_manual(values = plot_fills)+
  # make labels from 0 to 1 not percentage
  scale_x_continuous(breaks = c(0, 50, 100),
                     labels = c("0", "0.5", "1"))+
  # no labels on y axis
  theme(axis.text = element_text(size = 8),
        strip.text = element_text(size = 18),
        legend.position = "top",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18))+
  # blank layer to get a linetype legend
  geom_line(data = df[1:2,],
            aes(linetype=sex, colour = NA),
            size = 1.2,
            show.legend = TRUE)+
  guides(linetype = guide_legend(override.aes=list(linetype = c(1,3))),
         colour = FALSE,
         fill = FALSE)
setwd(plotDir)

ggsave(paste0("chr_ridgeplot_adultbrain.pdf"), plot = ridge_plot, width = 6, height = 3)


# make  ridgeplot liver and spleen - supp figure
df_subset <- df %>%
  filter(tissue %in% c("liver", "spleen"))
  ridge_plot <- df_subset %>%
    ggplot(aes(x = value,
             y = sex,
             fill = tissue,
             colour = tissue,
             group = sex, 
             linetype = sex))+
  # make ridgeplots quite see-through
  geom_density_ridges(size = 1.2,
                      alpha=0.2,
                      show.legend = F)+
  theme_bw()+
  xlab("CpG methylation")+
  ylab("density")+
  facet_grid(chromosome ~ tissue,
             drop = TRUE)+
  scale_fill_manual(values = plot_fills[2:3])+
  scale_linetype_manual(values = c(1, 3))+
  scale_colour_manual(values = plot_fills[2:3])+
  # make labels from 0 to 1 not percentage
  scale_x_continuous(breaks = c(0, 50, 100),
                     labels = c("0", "0.5", "1"))+
  # no labels on y axis
  theme(axis.text = element_text(size = 8),
        strip.text = element_text(size = 18),
        legend.position = "top",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
  axis.title = element_text(size = 18))+
  # blank layer to get a linetype legend
  geom_line(data = df_subset[1:2,],
            aes(linetype=sex, colour = NA),
            size = 1.2,
            show.legend = TRUE)+
  guides(linetype = guide_legend(override.aes=list(linetype = c(1,3))),
         colour = FALSE,
         fill = FALSE)
setwd(plotDir)

ggsave(paste0("chr_ridgeplot_adultliverspleen.pdf"), plot = ridge_plot, width = 7, height = 3)


# count number of CpGs for legend
count_CpG <- df %>% 
  group_by(sex, tissue, chromosome) %>%
  summarise(count = n())
setwd(dataDir)
write.table(count_CpG, "Count_CpG_adult_ridges.txt")
