# BJL
# 20211018
# xy average methylation plot 


# libraries
library(tidyverse)
library(scales)



# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/"
adultDir  <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"




# get data


# embryo data
setwd(dataDir)

embryo_data <- list.files(pattern = "_destranded_pooled_df_filt") %>% 
  set_names(str_replace(.,"_d.*", "")) %>%
  map_dfr(read_delim,
          delim = "\t",
          .id = "file") %>%
  mutate(percent = ratio * 100)


# adult data
setwd(adultDir)

adult_data <- list.files(pattern = "_destranded_pooled_df_filt") %>% 
  set_names(str_replace(.,"_d.*", "")) %>%
  map_dfr(read_delim,
          delim = "\t",
          .id = "file") %>%
  mutate(percent = ratio * 100) 

# summarise data
summarised_data <- bind_rows(embryo_data, adult_data) %>%
  group_by(file) %>%
  summarise(mean_percent = mean(percent),
            sd = sd(percent)) %>%
  mutate(file = factor(file,
                       levels = c("oocyte",
                                  "sperm",
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
                                  "spleen"))) %>%
  mutate(timepoint_plot_placement = factor(c("E1.5",
                                      "E2.5",
                                      "E3.5",
                                      "E4.5",
                                      "E5.5",
                                      "E6.5",
                                      "E7.5",
                                      "E7.5",
                                      "E7.5",
                                      "adult",
                                      "adult",
                                      "gamete",
                                      "gamete",
                                      "adult"),
                                      levels = c("gamete",
                                                 "E1.5",
                                                 "E2.5",
                                                 "E3.5",
                                                 "E4.5",
                                                 "E5.5",
                                                 "E6.5",
                                                 "E7.5",
                                                 "adult"))) %>%
  arrange(file)

write.table(summarised_data, "embryo_methylation_percents.txt")
# make plot

# colours
cols <- c("#BCBD22",
          "#9467BD",
          "#AEC7E8",
          "#98DF8A",
          "#DBDB8D",
          "#FFBB78",
          "#FF9896",
          "#E377C2",
          "#CD3333",
          "#f66e6e",
          "#8B2323",
          "#27408B",
          "#3A5FCD", 
          "#436EEE")

point_plot <- ggplot(data = summarised_data,
                    aes(x = timepoint_plot_placement,
                        y = mean_percent/100,
                        colour = file))+
  geom_point(size = 3)+
  geom_text(aes(label = file),
            hjust=-0.25,
            vjust=-0.1,
            colour = "black",
            size = 2)+
  scale_colour_manual(values = cols,
                      )+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("average CpG methylation")+
  xlab("")
point_plot
ggsave(file = paste0(plotDir, "point_plot.pdf"), width = 7, height = 3, useDingbats=FALSE)



