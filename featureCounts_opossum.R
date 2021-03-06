# 20191007
# BJL

#  this script is launched from a bash script using Rscript with input arguments

# takes input sam files and uses featureCounts from RsubRead to do read summarisation
# here is applied to rna-seq data from opossum adult JvDB corss
# produces a gene expression matrix with dimensions: nrows = ngenes, ncols = nsamples

# libraries
library(Rsubread)

# arguments as defined in the bash submission script
args <- commandArgs(trailingOnly = T)
indir <- args[1]
annot <- args[2]
outmatrix <- args[3]
outfile <- args[4]

# move to input directory
setwd(indir)
print(paste("We have moved to", getwd()))

# define a character vector of inout files
infiles <- list.files(pattern = ".bam")

# run featureCounts command to quantify transcripts

myFeatureCounts <- featureCounts(files = infiles,
                                 annot.ext = annot,
                                 isGTFAnnotationFile = TRUE,
                                 isPairedEnd=TRUE, 
                                 useMetaFeatures = TRUE,
                                 countMultiMappingReads = FALSE) 

# save the whole output list object to file
save(myFeatureCounts, file = outfile)

# Extract counts matrix
# save as table
countsMatrix <- myFeatureCounts$counts

write.table(countsMatrix, file = outmatrix)




sessionInfo()                              