# 20211031
# BJL
# repeat expression summary plot


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
library(rtracklayer) 
library(ComplexHeatmap)
library(lubridate)
library(ggrepel)
library(plotrix)



# directories
datadir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/repeat_expression/"
annotdir <- "/camp/lab/turnerj/working/Bryony/annotations/"
plotdir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots"


# data

# repeat annot
setwd(annotdir)
repeat_annot <- read_delim("20210126_UCSC_monDom5_repeatmasker_sorted_shifted_classfamily.txt", delim = "\t") # import the repeat track you built

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
  

# data
setwd(datadir)

# did once
# tidied Telescope output of rna-seq reads overlapping repeats
# cpm
#counts <- readRDS("telescope_report_tidied_cpm_merged.RDS") %>% 
#  # remove other summary counts e.g. "__no_feature" 
#  filter(transcript %in% repeat_annot$locus) %>%
#  # make sample_names column to match metadata
#  separate(source, sep = "[.]", into = c("sample_id")) %>%
#  # for consistency with repeat annotation object
#  dplyr::rename(locus= transcript) %>%
#  #make column per sample
#  spread(key = "sample_id", value = "cpm")
#counts[is.na(counts)] <- 0 # restore zeroes for NAs
#
#saveRDS(counts, "20210201_telescope_report_tidied_merged_cpm_reorganised.RDS")

# get expression data
# cpm object previously generated
counts <- readRDS("20210201_telescope_report_tidied_merged_cpm_reorganised.RDS")

# join with repeat annotation info
counts <- counts %>% left_join(select(repeat_annot, locus, repFamily, repClass, repName), by ="locus")

# plot
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


repeat_classes <- c("LTR", "LINE", "SINE")
for(i in 1:length(repeat_classes)){
  
  tmp <- counts %>%
    # get repeat class to work on
    filter(repClass == repeat_classes[i]) %>% 
    # reorganise longer
    pivot_longer(!c(locus, repName, repClass, repFamily),
                 names_to = "sample_id",
                 values_to = "count") %>%
    # exclude timepoints/cell types not of interest
    filter(sample_id %in% sample_data_keep$sample_id) %>%
    # join with metadata
    left_join(select(sample_data_keep, age_num, age_cell_type_plot, sample_id), by = "sample_id") %>%
    # group by cell type/timepoint and get total cpm for repeat
    group_by(repFamily, age_cell_type_plot, sample_id, age_num) %>%
    summarise(total = sum(count)) %>%
    # group by cell type/timepoint and get mean and se
    group_by(repFamily, age_cell_type_plot, age_num) %>%
    summarise(mean = mean(total),
              se = std.error(total))
  #write_delim(tmp, path = paste0(datadir, "_", repeat_classes[i], "summarised_counts.txt"), delim = "\t") 

  first_line <- tmp %>%
    filter(!age_cell_type_plot %in% lineage)
  
  second_line <- tmp %>%
    filter(age_cell_type_plot %in% c("6.5", "7.5_epi"))
  
  third_line <- tmp %>%
    filter(age_cell_type_plot %in% c("6.5", "7.5_hypo"))
  
  fourth_line <- tmp %>%
    filter(age_cell_type_plot %in% c("6.5", "7.5_te"))
  
  
   # plot timepoint on x axis against total cpm on y axis          
  ggplot(data = tmp,
         aes(x = age_num,
             y = mean,
             group = repFamily,
             colour = age_cell_type_plot)) +
    geom_pointrange(aes(ymin = mean - (1.96*se), ymax = mean + (1.96*se)), 
                    #position = position_jitter(width = 0.2, seed = 123)
    ) +
    geom_line(data=first_line, aes(x = age_num, y = mean, group = repFamily), size = 1)+
    geom_line(data=second_line, aes(x = age_num, y = mean, group = repFamily), colour = "#CD3333", size = 1)+
    geom_line(data=third_line, aes(x = age_num, y = mean, group = repFamily),  colour = "#f66e6e", size = 1)+
    geom_line(data=fourth_line, aes(x = age_num, y = mean, group = repFamily), colour = "#8B2323", size = 1)+
    scale_colour_manual(values = cols)+
    ggtitle(repeat_classes[i], )+
    ylab("")+
    xlab("")+
    #ylim()+
    theme_mine()+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank())
  
  ggsave(paste0(repeat_classes[i], ".pdf"), width = 6, height = 6, useDingbats=FALSE)
}