### adding "chr" to chromosome names
    head -n 5 Monodelphis_domestica.monDom5.97.gtf > head.gtf
    tail -n +6 Monodelphis_domestica.monDom5.97.gtf | awk '$1="chr"$1' | \
    cat head.gtf - > Monodelphis_domestica.monDom5.97_mod.gtf
    

### add new DNMT1 annotation from release 98 to modified release 97
wget ftp://ftp.ensembl.org/pub/release-98/gtf/monodelphis_domestica/ \
        Monodelphis_domestica.ASM229v1.98.chr.gtf.gz
    gzip -d Monodelphis_domestica.ASM229v1.98.chr.gtf.gz
    grep "ENSMODG00000048820" Monodelphis_domestica.ASM229v1.98.chr.gtf > dnmt1a
    grep "ENSMODG00000005463" Monodelphis_domestica.ASM229v1.98.chr.gtf > dnmt1b
    cat dnmt1a dnmt1b > dnmt1
    cut -f1-8 dnmt1 > cols
    cut -f9 dnmt1 > att
    awk '{ $1="chr" $1; } 1' OFS='\t' cols > cols_mod
    paste cols_mod att > dnmt1_mod
    mv dnmt1_mod /camp/lab/turnerj/working/shared_projects/OOPs/genome
    grep -v "ENSMODG00000005463" Monodelphis_domestica.monDom5.97_mod_shifted_XY.gtf \
        > gtf_tmp
    cat gtf_tmp dnmt1_mod > Monodelphis_domestica.monDom5.97_mod_shifted_XY_DNMT1.gtf
    bedtools sort -i Monodelphis_domestica.monDom5.97_mod_shifted_XY_DNMT1.gtf
    rm gtf_tmp dnmt1_mod
