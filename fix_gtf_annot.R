    ### delete extended Hprt1 exon and shift positions post-RSX ###
    
    gtf <- read.table("Monodelphis_domestica.monDom5.97_mod.gtf", header = F, sep = "\t")
    colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", 
            "frame", "attribute")
            
    c1 <- grepl("chrX", gtf$seqname)
    c2 <- grepl(35567216, gtf$start)
    c3 <- grepl("gene", gtf$feature)
    gtf[c1 & c2 & c3, ]$end <- 35598004
    
    c1 <- grepl("chrX", gtf$seqname)
    c2 <- grepl(35589163, gtf$start)
    c3 <- grepl("transcript", gtf$feature)
    gtf[c1 & c2 & c3, ]$end <- 35590230
    
    gtf <- gtf[!(grepl(35817805, gtf$start) & grepl("chrX", gtf$seqname)), ] 
    
    genes_post_rsx <- gtf[(gtf$start >= 35734441) & grepl("chrX", gtf$seqname), ]
    genes_post_rsx$end <- genes_post_rsx$end + 1200
    genes_post_rsx$start <- genes_post_rsx$start + 1200
    
    write.table(gtf, file = "Monodelphis_domestica.monDom5.97_mod_shifted.gtf", 
                sep = "\t", quote = F, col.names = F, row.names = F)
    rm(list = ls())
    
    
    ### change coordinates of RSX/XSR based on location of gaps ###
    # calculations of new coordinates were done in IGV and Excel
    
    gtf <- read.table("RSX_XSR.tmp.gtf", header = F, sep = "\t")
    colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", 
            "frame", "attribute")
    rsx <- gtf[1:4, ]
    xsr <- gtf[5:10, ]
    
    rsx$start <- c(35603406, 35603406, 35635154, 35652336) 
    rsx$end <- c(35654007, 35603519, 35635231, 35654007) 
    
    xsr$start <- c(35605416, 35652322, 35634737, 35606592, 35605416, 35605421) 
    xsr$end <- c(35652816, 35652816, 35634872, 35633831, 35605715, 35605423) 
    
    gtf_new <- rbind(rsx, xsr)
    gtf_new$source <- rep("Jasmin", 10)
    gtf_new$attribute <- gsub("Mahesh", "Turner lab", gtf_new$attribute)
    gtf_new$seqname <- rep("chrX", 10)
    gtf_new$freq <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 1)
    gtf_out <- gtf_new[rep(row.names(gtf_new), gtf_new$freq), 1:9]
    gtf_out$feature <- c("gene", "transcript", "exon", "CDS", "exon", "CDS", "exon", 
                        "CDS", "gene", "transcript", "exon", "CDS", "exon",
                        "CDS", "exon", "CDS", "exon", "CDS", "stop_codon") 
    
    write.table(gtf_out, file = "RSX_XSR_mod.gtf", sep = "\t", quote = F, 
                col.names = F, row.names = F)
