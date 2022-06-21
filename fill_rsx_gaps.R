    library(seqinr)
    setwd("/camp/lab/turnerj/working/jasmin/data/MonDom5")
    
    # MonDom5 X chromosome sequence from public download
    
    chrX <- read.fasta("mdm_ref_MonDom5_chrX.fa", forceDNAtolower = F, as.string = T)
    
    # BAC file equals Supp File S3 in Sprague et al., RNA, 2019
    # RSX file equals Supp File S2 in Sprague et al., RNA, 2019
    
    bac <- read.fasta("/camp/lab/turnerj/working/resources/opossum/rsx_bac_nanopore.fasta", 
        forceDNAtolower = F, as.string = T)
    
    rsx <- read.fasta("/camp/lab/turnerj/working/resources/opossum/rsx_with_introns.fasta", 
        forceDNAtolower = F, as.string = T)
    
    # positions for gap and fill seq in MonDom5 and RSX_BAC_nanopore
    # from Supp Tab S3 in Sprague et al., RNA, 2019
    
    part1 <- substr(chrX, 1, 35622568)
    gap1 <- substr(bac, 105578, 105938)
    part2 <- substr(chrX, 35623970, 35636736)
    gap2 <- substr(bac, 118744, 118832)
    part3 <- substr(chrX, 35637246, 35638910)
    gap3 <- substr(bac, 120502, 134312)
    part4 <- substr(chrX, 35650055, 35690787)
    gap4 <- substr(bac, 174820, 175832)
    part5 <- substr(chrX, 35691813, 35710775)
    gap5 <- substr(bac, 194845, 196043)
    part6 <- substr(chrX, 35711970, nchar(chrX))
    
    new_seq <- paste0(part1, gap1, part2, gap2, part3, gap3, part4, gap4, part5, gap5, part6)
    
    write.fasta(new_seq, name = "chrX gaps filled with bac in MonDom5|JZ|220819", 
        file.out = "mdm_ref_MonDom5_chrX_gaps_filled.fa", as.string = T)
    
    chrX.filled <- read.fasta("mdm_ref_MonDom5_chrX_gaps_filled.fa", 
        forceDNAtolower = F, as.string = T)
    
    
    > sessionInfo()
    R version 3.3.1 (2016-06-21)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: CentOS Linux 7 (Core)
    
    locale:
     [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
     [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
     [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
     [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
     [9] LC_ADDRESS=C               LC_TELEPHONE=C
    [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C
    
    attached base packages:
    [1] stats4    parallel  stats     graphics  grDevices utils     datasets
    [8] methods   base
    
    other attached packages:
    [1] Biostrings_2.40.2   XVector_0.12.0      IRanges_2.6.1
    [4] S4Vectors_0.10.2    BiocGenerics_0.18.0 seqinr_3.2-0
    
    loaded via a namespace (and not attached):
    [1] zlibbioc_1.18.0 tools_3.3.1     ade4_1.7-4

