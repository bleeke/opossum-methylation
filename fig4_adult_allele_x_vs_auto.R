# BJL
# distribution of methylation at X chromosome vs autosomes 
# adult data allele - specific


# libraries
library(tidyverse)
library(scales)
library(ggridges)
library(methylKit)


# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/allele_specific_methylKit_out"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"


# import data
setwd(dataDir)

# data preparation
# did once

## for autosomes
#sex_alleles <- c("_male_genome1",
#                 "_male_genome2",
#                 "_female_genome1",
#                 '_female_genome2')
#
#tissues <- c("brain",
#             "liver",
#             "spleen")
#
## to put loop output in
#df_auto <- data.frame()
#
## read in male and female files for each tissue
## find sites in common
## calculate percentage methylation
## reorganise and save to dataframe
#for(i in 1:length(tissues)){
#  
#  samples <- as.list(paste0(tissues[i],
#                            sex_alleles))
#  
#  condition <- c(1:length(samples))
#  
#  files <- as.list(paste0(samples,
#                          "_destranded_pooled_df_filt5.csv")) 
#  
#  adult_data <- methRead(location = files,
#                         sample.id = samples,
#                         assembly = "monDom5",
#                         pipeline = list(fraction = TRUE,
#                                         chr.col = 1,
#                                         start.col = 2,
#                                         end.col = 3,
#                                         coverage.col = 5,
#                                         strand.col = 4,
#                                         freqC.col = 8),
#                         context = "CpG",
#                         resolution = "base",
#                         treatment = condition,
#                         mincov = 5)
#  
#  # unite to generate methylBase object containing only those sites present in all samples
#  adult_united <- methylKit::unite(adult_data)
#  # extract percentage methylation
#  adult_pc <- percMethylation(adult_united,
#                              rowids = TRUE)
#  adult_pc <- adult_pc %>%
#    as.data.frame() %>%
#    rownames_to_column("position") %>%
#    pivot_longer(cols = !position,
#                 names_to = c("tissue", "sex", "genome"),
#                 names_sep = "_") %>%
#    separate(position, into = c("chr", "start", "end"),
#             sep = "[.]") %>%
#  # remove data from pseudoY, chrX and mtDNA
#  dplyr::filter(chr != "chrMT" & chr != "chrY" & chr != "chrX")
#  
#  df_auto <- rbind(df_auto, adult_pc)
#}
#
#
## for chrX
#sex_alleles <- c("_male_genome1",
#                 # no paternal allele for male chrX
#                 "_female_genome1",
#                 '_female_genome2')
#
#tissues <- c("brain",
#             "liver",
#             "spleen")
#
#
## to put loop output in
#df_X <- data.frame()
#
## read in male and female files for each tissue
## find sites in common
## calculate percentage methylation
## reorganise and save to dataframe
#for(i in 1:length(tissues)){
#  
#  samples <- as.list(paste0(tissues[i],
#                            sex_alleles))
#  
#  condition <- c(1:length(samples))
#  
#  files <- as.list(paste0(samples,
#                          "_destranded_pooled_df_filt5.csv")) 
#  
#  adult_data <- methRead(location = files,
#                         sample.id = samples,
#                         assembly = "monDom5",
#                         pipeline = list(fraction = TRUE,
#                                         chr.col = 1,
#                                         start.col = 2,
#                                         end.col = 3,
#                                         coverage.col = 5,
#                                         strand.col = 4,
#                                         freqC.col = 8),
#                         context = "CpG",
#                         resolution = "base",
#                         treatment = condition,
#                         mincov = 5)
#  
#  # unite to generate methylBase object containing only those sites present in all samples
#  adult_united <- methylKit::unite(adult_data)
#  # extract percentage methylation
#  adult_pc <- percMethylation(adult_united,
#                              rowids = TRUE)
#  adult_pc <- adult_pc %>%
#    as.data.frame() %>%
#    rownames_to_column("position") %>%
#    pivot_longer(cols = !position,
#                 names_to = c("tissue", "sex", "genome"),
#                 names_sep = "_") %>%
#    separate(position, into = c("chr", "start", "end"),
#             sep = "[.]") %>%
#  # get data from only chrX
#  dplyr::filter(chr == "chrX")
#  
#  df_X <- rbind(df_X, adult_pc)
#}
#
## combine
#df <- rbind(df_auto, df_X) %>%   
#  # categorise data as chrX or autosome - note autosomes includes chrUn
#  mutate(chromosome = case_when(chr == "chrX" ~ "chrX",
#                                chr != "chrX" ~ "autosome")) %>%
#  mutate(chromosome = factor(chromosome,
#                             levels = c("autosome",
#                                        "chrX"))) %>%
#  tidyr::unite(col = "plot_type",
#               sex, genome,
#               sep = "_",
#               remove = FALSE) %>%
#  mutate(plot_type = case_when(plot_type == "male_genome1" ~ "male_maternal_allele",
#                               plot_type == "male_genome2" ~ "male_paternal_allele",
#                               plot_type == "female_genome1" ~ "female_maternal_allele",
#                               plot_type == "female_genome2" ~ "female_paternal_allele")) %>%
#  mutate(plot_type = factor(plot_type,
#         levels = c("female_paternal_allele",
#                    "female_maternal_allele",
#                    "male_maternal_allele",
#                    "male_paternal_allele")))
#
#setwd(dataDir)
#saveRDS(df, "figure3_adult_allele_data.RDS")

df <- readRDS("figure3_adult_allele_data.RDS")
# plot
setwd(plotDir)

# fill colour for plot
plot_fills <- c("royalblue4",
                "royalblue3", 
                "royalblue2")



# make plot brain - main fig3
ridge_plot <- df %>%
  filter(tissue == "brain") %>%
  ggplot(aes(x = value,
             y = plot_type,
             fill = tissue,
             colour = tissue,
             group = plot_type, 
             linetype = plot_type))+
  # make ridgeplots quite see-through
  geom_density_ridges(size = 1.2,
                      alpha=0.2,
                      show.legend = F)+
  theme_bw()+
  xlab("CpG methylation")+
  ylab("density")+
  facet_grid(chromosome ~ tissue,
             drop = TRUE)+
  scale_fill_manual(values = plot_fills,
                    guide = "none")+
  scale_colour_manual(values = plot_fills,
                      guide = "none")+
  # make labels from 0 to 1 not percentage
  scale_x_continuous(breaks = c(0, 50, 100),
                     labels = c("0", "0.5", "1"))+
  scale_linetype_manual(values = c(1, 3, 2, 4))+
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "top")+
  # blank layer to get a linetype legend
  geom_line(aes(linetype=plot_type, colour = NA),
            size = 1,
            show.legend = TRUE)+
  guides(linetype = guide_legend(
  
                                 override.aes=list(linetype = c(1, 3, 2, 4))),
         colour = FALSE,
         fill = FALSE)
setwd(plotDir)

ggsave(paste0("chr_ridgeplot_allele_adult_brain.pdf"), plot = ridge_plot, width = 6, height = 3.3)




# make plot liver and spleen - supp fig
ridge_plot <- df %>%
  filter(tissue %in% c("liver", "spleen")) %>%
  ggplot(aes(x = value,
             y = plot_type,
             fill = tissue,
             colour = tissue,
             group = plot_type, 
             linetype = plot_type))+
  # make ridgeplots overlap quite a lot 
  # and be quite see-through
  geom_density_ridges(size = 1.2,
                      alpha=0.2,
                      show.legend = F)+
  theme_bw()+
  xlab("CpG methylation")+
  ylab("density")+
  facet_grid(chromosome ~ tissue,
             drop = TRUE)+
  scale_fill_manual(values = plot_fills[2:3],
                    guide = "none")+
  scale_colour_manual(values = plot_fills[2:3],
                      guide = "none")+
  # make labels from 0 to 1 not percentage
  scale_x_continuous(breaks = c(0, 50, 100),
                     labels = c("0", "0.5", "1"))+
  scale_linetype_manual(values = c(1, 3, 2, 4))+
  theme(strip.text = element_text(size = 12),
        legend.position = "top")+
  # blank layer to get a linetype legend
  geom_line(aes(linetype=plot_type, colour = NA),
            size = 1.2,
            show.legend = TRUE)+
  guides(linetype = guide_legend(nrow = 2,
                                 override.aes=list(linetype = c(1, 3, 2, 4))),
         colour = FALSE,
         fill = FALSE)

ggsave(paste0("chr_ridgeplot_allele_adult_liverspleen.pdf"), plot = ridge_plot, width = 7, height = 3.3)


# count number of CpGs for legend
count_CpG <- df %>% 
  group_by(sex, tissue, plot_type, chromosome) %>%
  summarise(count = n())
setwd(dataDir)
write.table(count_CpG, "Count_CpG_adult_allele_ridges.txt")
