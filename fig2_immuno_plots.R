# 20210322
# BJL
# immuno quantification graphs



# packages
library(tidyverse)

# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/immuno_quant"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots"

# import data
setwd(dataDir)
data <- read.delim("5meC_H3K9me3_29.5_hpm.txt")


# plot1
immuno_plot <- data %>%
  ggplot(aes(x = antibody,y = ratio_pat_mat, 
             fill = antibody)) +
  scale_fill_manual(values = c("chartreuse3",
                               "mediumorchid2"))+
  geom_boxplot(alpha = 0.6)+
  geom_point()+
  ylab("signal intensity (paternal/maternal)")+
  theme_bw()+
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylim(0,1.5)
  
ggsave("pat_mat_ratio_immuno_5meC_H3K9me3.pdf",
       path = plotDir,
       plot = immuno_plot,
       width = 5.5, height = 5)


# import data
setwd(dataDir)
data <- read.delim("H3_H3K9me3_29.5_hpm.txt")


# plot2
immuno_plot <- data %>%
  ggplot(aes(x = antibody,y = ratio_pat_mat, 
             fill = antibody)) +
  scale_fill_manual(values = c("chartreuse3",
                               "mediumorchid2"))+
  geom_boxplot(alpha = 0.6)+
  geom_point()+
  ylab("signal intensity (paternal/maternal)")+
  theme_bw()+
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylim(0,1.5)

ggsave("pat_mat_ratio_immuno_H3_H3K9me3.pdf",
       path = plotDir,
       plot = immuno_plot,
       width = 5.5, height = 5)


# import data
setwd(dataDir)
data <- read.delim("5hmC_H3K9me3_29.5_hpm.txt")


# plot2
immuno_plot <- data %>%
  ggplot(aes(x = antibody,y = ratio_pat_mat, 
             fill = antibody)) +
  scale_fill_manual(values = c("chartreuse3",
                               "mediumorchid2"))+
  geom_boxplot(alpha = 0.6)+
  geom_point()+
  ylab("signal intensity (paternal/maternal)")+
  theme_bw()+
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylim(0,1.5)

ggsave("pat_mat_ratio_immuno_5hmC_H3K9me3.pdf",
       path = plotDir,
       plot = immuno_plot,
       width = 5.5, height = 5)
