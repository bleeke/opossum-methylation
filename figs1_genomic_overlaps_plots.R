# BJL
# plot genomic distributions of CpGs captured in BSseq libraries

# libraries
library(tidyverse)
library(RColorBrewer)

# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out"
adultDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"

# get data
setwd(dataDir)
files <- list.files(pattern = "genomic_overlaps")
data <- files %>%
  set_names() %>%
  map_dfr(readRDS,
          .id = "file") %>%
  separate(file, sep = "_g", into = c("file"))

setwd(adultDir)
files <- list.files(pattern = "genomic_overlaps")
adult_data <- files %>%
  set_names() %>%
  map_dfr(readRDS,
          .id = "file") %>%
  separate(file, sep = "_g", into = c("file"))

setwd(annotDir)
theoretical_data <- readRDS("all_mondom_cpgs_genomic_overlaps.RDS") %>% 
  mutate(file = "all_CpGs")

df <- rbind(data, adult_data, theoretical_data) %>%
  mutate(sample = factor(file,
                         levels = c("all_CpGs",
                                    "sperm", "oocyte",
                                    "1.5_embryo","2.5_embryo", "3.5_embryo",
                                    "4.5_embryo", "5.5_embryo", "6.5_embryo",
                                    "7.5_embryo", "7.5_ED", "7.5_TE",
                                    "brain", "liver", "spleen")))

# plot proportions of CpGs in each genomic category

# make plot colours: require palette of 14 so add extras manually
plot_colours <- c(brewer.pal(12, "Set3"), "#CD3333", "#483D8B")


setwd(plotDir)

df %>% 
  filter(context == "chr_df") %>%
  mutate(type = case_when(type == "chrX" ~ "chrX",
                          type == "chrUn" ~ "chrUn",
                          type != "chrX" & type != "chrUn" ~ "autosome")) %>%
  mutate(type = factor(type, levels = c("chrX", "chrUn", "autosome"))) %>%
  group_by(type, sample) %>%
  summarise(perc= sum(perc)) %>%
  ggplot(aes(x = sample,
             y = perc))+
  geom_col(aes(fill = type),
           colour="black")+
  scale_fill_manual(values = plot_colours,
                    name = "percentage of CpG sites by chromosome:")+
  theme_bw()+
  theme(axis.text = element_text(size = 26),
        axis.title= element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position="top")+
  scale_x_discrete(labels = c("all\nCpGs",
                              "sperm", "oocyte",
                              "E1.5","E2.5", "E3.5",
                              "E4.5", "E5.5", "E6.5",
                              "E7.5", "E7.5\nED", "E7.5\nTE",
                              "brain", "liver", "spleen"))
ggsave("chrom_CpGs.tiff", width = 26, height = 4)


df %>% 
  filter(context == "CGI_df") %>%
  mutate(type = factor(type, levels = c("intergenic_CGI", "intragenic_CGI", "prom_CGI", "no_CGI")))%>%
  ggplot(aes(x = sample,
             y = perc))+
  geom_col(aes(fill = type),
           colour = "black")+
  scale_fill_manual(values = plot_colours,
                    name = "percentage of CpG sites by CGI type:",
                    labels = c("intergenic", "intragenic", "promoter", "no CGI"))+
  theme_bw()+
  theme(axis.text = element_text(size = 26),
        axis.title= element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position="top")+
  scale_x_discrete(labels = c("all\nCpGs",
                              "sperm", "oocyte",
                              "E1.5","E2.5", "E3.5",
                              "E4.5", "E5.5", "E6.5",
                              "E7.5", "E7.5\nED", "E7.5\nTE",
                              "brain", "liver", "spleen"))
ggsave("CGI_CpGs.tiff", width = 26, height = 4)

df %>% 
  filter(context == "genes_df") %>%
  mutate(type = factor(type, levels = c("promoter", "intragenic","intergenic")))%>%
  ggplot(aes(x = sample,
             y = perc))+
  geom_col(aes(fill = type),
           colour = "black")+
  scale_fill_manual(values = plot_colours,
                    name = "percentage of CpG sites by gene region:",
                    labels = c("promoter", "gene", "intergenic"))+
  theme_bw()+
  theme(axis.text = element_text(size = 26),
        axis.title= element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position="top")+
  scale_x_discrete(labels = c("all\nCpGs",
                              "sperm", "oocyte",
                              "E1.5","E2.5", "E3.5",
                              "E4.5", "E5.5", "E6.5",
                              "E7.5", "E7.5\nED", "E7.5\nTE",
                              "brain", "liver", "spleen"))
ggsave("gene_CpGs.tiff", width = 26, height = 4)


df %>% 
  filter(context == "repeat_df") %>%
  mutate(type = factor(type,
                       levels = c( "Satellite", "RC", "Unknown", "snRNA",
                                   "tRNA", "scRNA", "Low_complexity",
                                   "Simple_repeat", "rRNA", "DNA",
                                   "LTR", "SINE", "LINE", "no_repeat"))) %>%
  ggplot(aes(x = sample,
             y = perc))+
  geom_col(aes(fill = type),
           colour = "black")+
  scale_fill_manual(values = rev(plot_colours),
                    name = "percentage of CpG sites by repeat class:",
                    labels = c( "Satellite", "RC", "Unknown", "snRNA",
                                "tRNA", "scRNA", "Low complexity",
                                "Simple repeat", "rRNA", "DNA",
                                "LTR", "SINE", "LINE", "no repeat"))+
  theme_bw()+
  theme(axis.text = element_text(size = 26),
        axis.title= element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position="top")+
  scale_x_discrete(labels = c("all\nCpGs",
                              "sperm", "oocyte",
                              "E1.5","E2.5", "E3.5",
                              "E4.5", "E5.5", "E6.5",
                              "E7.5", "E7.5\nED", "E7.5\nTE",
                              "brain", "liver", "spleen"))
ggsave("repeat_CpGs.tiff", width = 26, height = 4.7)