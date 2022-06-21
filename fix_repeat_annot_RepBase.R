    setwd("/camp/lab/turnerj/working/Bryony/annotations")
    
    ### read in the data and sort
    annot <- get(load("20191101_opossum_RepeatMasker.R"))
    annot <- annot[order(annot$query_sequence, annot$position_in_query_begin), ]
    
    
    ### check for annots in RSX locus
    
    get_chrX <- grepl("chrX", annot$query_sequence)
    rsx_start <- annot$position_in_query_begin >= 35603406
    rsx_end <- annot$position_in_query_begin <= 35652800
    
    annot[get_chrX & rsx_start & rsx_end, 5:7]
    # revealed 16 annot within RSX locus
    # 1-10 are upstream of gap 1 and can be ignored
    # 11+12 are between gaps 1 and 2 (+1040)
    # 13+14 are between gaps 2 and 3 (+1460)
    # 15+16 are between gaps 3 and 4 (-1207)
    
    
    ### fix annotation of annot within RSX locus
    
    rep11 <- grepl(35633089, annot$position_in_query_begin)
    rep12 <- grepl(35636074, annot$position_in_query_begin)
    annot[get_chrX & (rep11 | rep12), 6:7] <- annot[get_chrX & (rep11 | rep12), 6:7] - 1040
    
    rep13 <- grepl(35637964, annot$position_in_query_begin)
    rep14 <- grepl(35638686, annot$position_in_query_begin)
    annot[get_chrX & (rep13 | rep14), 6:7] <- annot[get_chrX & (rep13 | rep14), 6:7] - 1460
    
    rep15 <- grepl(35650512, annot$position_in_query_begin)
    rep16 <- grepl(35651385, annot$position_in_query_begin)
    annot[get_chrX & (rep15 | rep16), 6:7] <- annot[get_chrX & (rep15 | rep16), 6:7] + 1207
    
    
    
    ### shift coordinates of annots downstream of RSX by 1200 bp
    
    annot[get_chrX & (annot$position_in_query_begin > 35625350), 6:7] <-
         annot[get_chrX & (annot$position_in_query_begin > 35625350), 6:7] + 1200
    
    
    save(annot, file = "20191101_opossum_RepeatMasker_shifted.RDS")
