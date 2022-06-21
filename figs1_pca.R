# 20210203
# BJL
# PCA plot for opossum BS-seq 

# libraries
library(tidyverse)
library(methylKit)
library(lubridate)
library(ggfortify)

# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methyl_extract/bismark_CpG_report/"
adultDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methyl_extract/bismark_CpG_report/"
sampleDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots"

# load data 

setwd(sampleDir)

# embryo data
# import sample metadata
emb_sample_data <- read_delim( "PM18153_sample_data.txt", delim = "\t") %>%
  # filter for samples that passed mapping rate QC
  filter(include_in_analysis == "include")

# define sample names
emb_samples <- as.list(emb_sample_data$sample_id)
  
emb_files <- as.list(paste0(dataDir,
                              "/",
                              emb_samples,
                              "_merged_fastq.gz_trimmed_bismark_bt2.deduplicated.CpG_report.txt"))
  
# define sample treatment
emb_design <- emb_sample_data$sample_design

# import sample metadata

adult_sample_data <- read_delim("adult_sample_data.txt", delim = "\t")

# adult data
# define sample names
adult_samples <- as.list(adult_sample_data$sample_id)

# define files
adult_files <- as.list(paste0(adultDir,
                              "/",
                              adult_samples,
                       "_merged_fastq.gz_trimmed_bismark_bt2.allele_flagged.deduplicated_sorted.CpG_report.txt"))

# define sample treatment
adult_design <- adult_sample_data$sample_design


# combine vectors for embryo and adult
files <- c(emb_files, adult_files)
samples <- c(emb_samples, adult_samples)
design <- c(emb_design, adult_design)

# make unified sample metadata for embryo and adult

sample_data <- bind_rows(emb_sample_data,
                         adult_sample_data) %>%
  dplyr::select(sample_id, sample_design) %>% 
  mutate(sample_design_fact = as.factor(sample_design))

#import
dat <- methRead(location = files,
                sample.id = samples,
                assembly = "monDom5",
                pipeline="bismarkCytosineReport",
                treatment = design,
                mincov = 1)

# unite with min.per.group = 0L to create methylBase object while keeping all data

united_dat <- unite(dat, min.per.group = 0L)

# PCA 
setwd(plotDir)
PCA <- PCASamples(united_dat,
                  obj.return = T)

# plot PCA
autoplot(PCA, data = sample_data, colour = 'sample_design_fact', size = 4) +
  theme_bw() +
  theme(axis.text = element_text(size =18),
        axis.title = element_text(size =20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  scale_color_manual(values = c("#BCBD22",
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
                                "#436EEE"),
                     labels = c("oocyte",
                                "sperm",
                                "E1.5",
                                "E2.5",
                                "E3.5",
                                "E4.5",
                                "E5.5",
                                "E6.5",
                                "E7.5",
                                "E7.5 ED",
                                "E7.5 TE",
                                "brain",
                                "liver",
                                "spleen"),
                     name = "samples")+
  labs(colour = "library")
ggsave("oposs_pca_0L.pdf", useDingbats = FALSE) 

# save screeplot
pdf(file = "oposs_screeplot.pdf")
PCASamples(united_dat,
           screeplot=T)
dev.off()
