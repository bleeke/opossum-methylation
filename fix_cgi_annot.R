    setwd("/camp/lab/turnerj/working/Bryony/annotations")
    
    ### read in the data
    cgi <- read.table("20180614_monDom5_CGI.bed", skip = 1)
    colnames(cgi) <- c("chr", "start", "end", "id")
    cgi <- cgi[order(cgi$chr, cgi$start),]
    
    
    ### check for CGIs in RSX locus
    
    head(cgi[grepl("chrX", cgi$chr) & (cgi$start >= 35603406) & (cgi$start <= 35652800), ])
    # revealed one cgi within RSX locus:
    # chrX 35626390 35626640 CpG:_27
    # substrat 1040 bp from either end (diff of gap 1)
    # new coordinates: 35625350-35625600
    
    ### fix annotation of CGI within RSX locus
    
    cgi[grepl("chrX", cgi$chr) & grepl(35626390, cgi$start), ]$start <- 35625350
    cgi[grepl("chrX", cgi$chr) & grepl(35626640, cgi$end), ]$end <- 35625600
    
    
    ### shift coordinates of CGIs downstream of RSX by 1200 bp
    
    cgi[grepl("chrX", cgi$chr) & (cgi$start > 35625350), ]$start <- 
        cgi[grepl("chrX", cgi$chr) & (cgi$start > 35625350), ]$start + 1200
    cgi[grepl("chrX", cgi$chr) & (cgi$end > 35625600), ]$end <- 
        cgi[grepl("chrX", cgi$chr) & (cgi$end > 35625600), ]$end + 1200
    
    
    write.table(cgi, file = "20180614_monDom5_CGI_shifted.bed", 
                sep = "\t", quote = F, col.names = F, row.names = F)
