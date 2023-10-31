

##post figR  process script (compare dorc <Ma et al 2021> with p2g <Granja et al 2019> method, overlap dar etc. )###


library('magrittr')



#readin DAR region 

regionSets <- readRDS('/home/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/cisTarget/RcisTarget/regionSets.rds')

sapply(regionSets, FUN = length)
#c6 7000 c3 7000 c9 5244 c5 6188 c7 6140 c1 5276 c8 2722 c2 4929 c4 1158


###read in figR cisCorr.filt.rds and p2g.filter.rds to try to link cistrome to gene################

cisCorr.filt <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/liger/figR_Kartha/FigR/PLA-8w-combined/cisCorr.filt.rds') #19674

length(unique(cisCorr.filt$Gene)) #10359

cisCorr <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/liger/figR_Kartha/FigR/PLA-8w-combined/cisCorr.rds') #99609

length(unique(cisCorr$Gene)) #15678


cisCorr.filt0.1 <- cisCorr %>% dplyr::filter(pvalZ <= 0.1)
#31481


# cistrome.test <- motif_tf.highConf.list[['cisbp__M4634']][['region']]

# table(cistrome.test %in% cisCorr.filt$PeakRanges)

# idx <- match(cistrome.test,cisCorr.filt$PeakRanges)
# idx <- idx[!is.na(idx)]

# cisCorr.filt[idx,'Gene']


p2g.filter <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/liger/peak_gene_link/result_good/gene_peak_links/out.p2g.filter.rds') #93053

p2g.filter.uniq <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/liger/peak_gene_link/result_good/gene_peak_links/out.p2g.filter.uniq.rds') #46558


length(unique(p2g.filter$gene))
#9922


length(unique(p2g.filter.uniq$gene))
#7678

p2g.filter.split <- split(p2g.filter, p2g.filter$gene)
len.dorc.p2g <- sapply(p2g.filter.split, nrow)
plot(1:length(len.dorc.p2g),sort(len.dorc.p2g), pch = 19, cex = 0.1 )

p2gGenes <- head(names(sort(len.dorc.p2g,decreasing = TRUE)),800)
p2gGenes.count <- head(sort(len.dorc.p2g,decreasing = TRUE),800)




grep('PAPPA',p2gGenes,value=TRUE)
'PAPPA2''PAPPA'

grep('STAT',p2gGenes,value=TRUE)
'STAT5B''STAT5A' 'STAT4'

grep('PPAR',p2gGenes,value=TRUE) #no
'PPARG'

grep('ESR',p2gGenes,value=TRUE)

grep('ESR',names(sort(len.dorc.p2g,decreasing = TRUE)),value=TRUE)
'ESRP2''ESR1''ESRRB''ESRRA'


dorcGenes <- readRDS('dorcGenes.rds')


length(intersect(p2gGenes,dorcGenes)) #144 if top 800 #100 if top 500

'INHBA''FLT1''PAPPA''SDC1''KANK1''SLC6A6''IL2RB''LINC00474''POMP''TRIB1''KCTD17''ITGA6''LAMB1''ADHFE1''TCF7L2''ANKH''RBPMS''SELL''ADRB1''BMP1''CDH1''PVT1''TMEM229B''GPD1L''MSI2''NLRC5''HSPB1''ILDR2''LAMC2''LHB''CBX4''SH3BP5''SMG6''SYDE1''TIMP2''F5''HERPUD1''KANK4''KCNA7''C1QTNF6''CTNNAL1''EPS8L1''HOPX''MFSD2A''NRP2''PARD3''THRB''DUSP14''GALNT18''MDFI''OAF''RPS6KA1''SMYD2''TMEM150C''CMTM7''COL4A2''DUSP1''KCNK3''RNF144B''SLC2A5''VWA2''ADGRG1''BMP7''CELA2B''CPXM2''LHPP''LRP5''PARP1''PCDH1''RNF19B''ADORA1''CTDSP1''DNMT1''LAMB4''TMEM245''TNRC18''ABHD17C''ARNTL''DGKZ''EHBP1''FOXI3''ITPK1''MAPRE3''PROSER2-AS1''ZNF431''AES''AGPAT4''ELFN2''FBN2''FGFR2''GALNT2''GNG7''LINC01119''LINC01509''MFSD12''MPP7''MTSS1''POU2F3''PPP1R14C''SLC22A11''TNFRSF1B''WLS''YPEL2''CACNA2D3''CARD16''CASD1''CCDC113''CHSY1''DCBLD1''EFNA5''ICE1''OSBPL9''SLC15A1''SLC43A2''SRGAP1''TFDP2''TNFAIP2''UBE4B''ASAP1''ATP2B4''EFHD2''MKLN1''PDLIM1''PLCG2''PSD4''SLC29A1''SSBP3''STAT4''UBA5''ALPL''CTNNBIP1''FHL2''HDAC11''ITGA5''KRT7''KRT80''LINC01411''NFATC2''PHACTR2''RASSF3''SIPA1L3''SLC12A9''TMPRSS3''TNIK'

#'INHBA''FLT1''PAPPA''SDC1''KANK1''SLC6A6''IL2RB''LINC00474''POMP''TRIB1''KCTD17''ITGA6''LAMB1''ADHFE1''TCF7L2''ANKH''RBPMS''SELL''ADRB1''BMP1''CDH1''PVT1''TMEM229B''GPD1L''MSI2''NLRC5''HSPB1''ILDR2''LAMC2''LHB''CBX4''SH3BP5''SMG6''SYDE1''TIMP2''F5''HERPUD1''KANK4''KCNA7''C1QTNF6''CTNNAL1''EPS8L1''HOPX''MFSD2A''NRP2''PARD3''THRB''DUSP14''GALNT18''MDFI''OAF''RPS6KA1''SMYD2''TMEM150C''CMTM7''COL4A2''DUSP1''KCNK3''RNF144B''SLC2A5''VWA2''ADGRG1''BMP7''CELA2B''CPXM2''LHPP''LRP5''PARP1''PCDH1''RNF19B''ADORA1''CTDSP1''DNMT1''LAMB4''TMEM245''TNRC18''ABHD17C''ARNTL''DGKZ''EHBP1''FOXI3''ITPK1''MAPRE3''PROSER2-AS1''ZNF431''AES''AGPAT4''ELFN2''FBN2''FGFR2''GALNT2''GNG7''LINC01119''LINC01509''MFSD12''MPP7''MTSS1''POU2F3''PPP1R14C''SLC22A11'



grep('STAT',dorcGenes,value=TRUE)
'STAT4'

grep('PPAR',dorcGenes,value=TRUE)


grep('ESR',dorcGenes,value=TRUE)



#############



##how many dars in dorc or p2g#####
table(mcols(regionSets[['c2']])$name %in% cisCorr.filt$PeakRanges )
FALSE  TRUE 
 4631   298

table(mcols(regionSets[['c2']])$name %in% cisCorr$PeakRanges )
FALSE  TRUE 
 3538  1391 

table(mcols(regionSets[['c2']])$name %in% p2g.filter$peak )
FALSE  TRUE 
 4201   728

table(mcols(regionSets[['c2']])$name %in% p2g.filter.uniq$peak )
FALSE  TRUE 
 4201   728


sapply( regionSets, FUN = function(x){ nrow(mcols(x))  }  ) 
c6 7000 c3 7000 c9 5244 c5 6188 c7 6140 c1 5276 c8 2722 c2 4929 c4 1158

sapply( regionSets, FUN = function(x){ table(mcols(x)$name %in% p2g.filter$peak)['TRUE']   }  ) 
#c6.TRUE 3805 c3.TRUE 3399 c9.TRUE 1495 c5.TRUE 1591 c7.TRUE 1205 c1.TRUE 938 c8.TRUE 510 c2.TRUE 728 c4.TRUE 141

sapply( regionSets, FUN = function(x){ table(mcols(x)$name %in% p2g.filter.uniq$peak)['TRUE']   }  ) 
#c6.TRUE 3805 c3.TRUE 3399 c9.TRUE 1495 c5.TRUE 1591 c7.TRUE 1205 c1.TRUE 938 c8.TRUE 510 c2.TRUE 728 c4.TRUE 141

sapply( regionSets, FUN = function(x){ table(mcols(x)$name %in% cisCorr.filt$PeakRanges)['TRUE']   }  )  #0.01
#c6.TRUE 912 c3.TRUE 871 c9.TRUE 479 c5.TRUE 562 c7.TRUE 541 c1.TRUE 485 c8.TRUE 236 c2.TRUE 298 c4.TRUE 91

sapply( regionSets, FUN = function(x){ table(mcols(x)$name %in% cisCorr.filt0.1$PeakRanges)['TRUE']   }  )
#c6.TRUE 1302 c3.TRUE 1223 c9.TRUE 724 c5.TRUE 894 c7.TRUE 778 c1.TRUE 715 c8.TRUE 343 c2.TRUE 482 c4.TRUE 145

#sapply( regionSets, FUN = function(x){ table(mcols(x)$hg19name %in% cisCorr.filt0.1$PeakRanges)['TRUE']   }  )
#all NA

###################get dar overlap with dorc , then do cistarget##########

regionSets_filter_dorc0.01 <- lapply( regionSets, FUN = function(x){  x[mcols(x)$name %in% cisCorr.filt$PeakRanges] } ) 

sapply(regionSets_filter_dorc0.01, length)
#c6 912 c3 871 c9 479 c5 562 c7 541 c1 485 c8 236 c2 298 c4 91

saveRDS(regionSets_filter_dorc0.01,'~/pwdex/placenta_10X_combine/02.snapATAC_harmony/cisTarget/i-cisTarget/DARs_overlap_dorc/regionSets_filter_dorc0.01.rds')


regionSets_filter_dorc0.01.df <- lapply(regionSets_filter_dorc0.01, function(x){
       x.df <- as.data.frame(x)
       x.df$id <- paste(x.df$name, '|' ,x.df$hg19name,sep='')
    
      x.df[,c('seqnames','start','end','id')]
    
   }
)


for(i in names(regionSets_filter_dorc0.01.df)){
    cat(i,'\n')
    write.table(regionSets_filter_dorc0.01.df[[i]],  file =  paste('~/pwdex/placenta_10X_combine/02.snapATAC_harmony/cisTarget/i-cisTarget/DARs_overlap_dorc/dar_',i,'.bed',sep=''),  col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
    
}




#############dar overlap with p2g.filter(p2g.filter.uniq the same) (h38 peak id)##################

regionSets_filter_p2g <- lapply(regionSets, FUN = function(x){  x[mcols(x)$name %in% p2g.filter.uniq$peak] } ) 

sapply(regionSets_filter_p2g, length)
#c6 3805 c3 3399 c9 1495 c5 1591 c7 1205 c1 938 c8 510 c2 728 c4 141


saveRDS(regionSets_filter_p2g,'~/pwdex/placenta_10X_combine/02.snapATAC_harmony/cisTarget/i-cisTarget/DARs_overlap_p2g/regionSets_filter_p2g.rds')


regionSets_filter_p2g.df <- lapply(regionSets_filter_p2g, function(x){
       x.df <- as.data.frame(x)
       x.df$id <- paste(x.df$name, '|' ,x.df$hg19name,sep='')
    
      x.df[,c('seqnames','start','end','id')]
    
   }
)


for(i in names(regionSets_filter_p2g.df)){
    cat(i,'\n')
    write.table(regionSets_filter_p2g.df[[i]],  file =  paste('~/pwdex/placenta_10X_combine/02.snapATAC_harmony/cisTarget/i-cisTarget/DARs_overlap_p2g/dar_',i,'.overlap_p2g.bed',sep=''),  col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
    
}





####################overlap peaks and links of dorc method and p2g method#########################

length(intersect(cisCorr.filt$PeakRanges,p2g.filter$peak ))
#7213

length(intersect(cisCorr$PeakRanges,p2g.filter$peak ))
#25476


share_peak <- intersect(cisCorr.filt$PeakRanges,p2g.filter$peak )


length(intersect(cisCorr.filt$Gene,p2g.filter$gene ))
#6346

sample(intersect(cisCorr.filt$Gene,p2g.filter$gene ),100)


idx1 <- match(share_peak,cisCorr.filt$PeakRanges)
idx2 <- match(share_peak,p2g.filter$peak)

cisCorr.filt.share <- cisCorr.filt[idx1,]
p2g.filter.share <- as.data.frame(p2g.filter[idx2,])

length(unique(cisCorr.filt[idx1,'Gene'])) #4268
length(unique(p2g.filter[idx2,'gene'])) #3438

share_gene <- intersect(unique(cisCorr.filt[idx1,'Gene']),unique(p2g.filter[idx2,'gene']))
#1877



grep('PAPPA',share_gene,value=TRUE)
'PAPPA2''PAPPA'
grep('FLT1',share_gene,value=TRUE)
'FLT1'
grep('CSH',share_gene,value=TRUE)

grep('LEP',share_gene,value=TRUE)
'LEP'

grep('STAT4',share_gene,value=TRUE)
'STAT4'

grep('STAT5A',share_gene,value=TRUE) #no

grep('STAT5B',share_gene,value=TRUE) #no

grep('AR',share_gene,value=TRUE)
'AR'



grep('NPAS2',share_gene,value=TRUE) #no

grep('FOSL2',share_gene,value=TRUE)
'FOSL2'

grep('JUN',share_gene,value=TRUE)
'JUNB''JUND'

grep('BACH',share_gene,value=TRUE) #no



grep('MYCN',share_gene,value=TRUE)
'MYCNUT'

grep('POU2F3',share_gene,value=TRUE) #no


grep('INHBA',share_gene,value=TRUE)
'INHBA'

grep('PSG8',share_gene,value=TRUE)
'PSG8'


grep('GCM1',share_gene,value=TRUE) #no
grep('SH3TC2',share_gene,value=TRUE)
'SH3TC2'



grep('ESRRG',share_gene,value=TRUE) #no

grep('PPARD',share_gene,value=TRUE)
'PPARD'

grep('TBX3',share_gene,value=TRUE)

grep('REL',share_gene,value=TRUE)#no


grep('ERVFRD-1',share_gene,value=TRUE)
'ERVFRD-1'

grep('PPARG',share_gene,value=TRUE)
'PPARG'


grep('CDH1',share_gene,value=TRUE)
'CDH1'

grep('DNMT1',share_gene,value=TRUE)
'DNMT1'

grep('TEAD4',share_gene,value=TRUE) #no


##peaks linked to marker gene

getLinks_gene <- function(gene.sel = NULL){
    share_n <- length(intersect(subset(cisCorr.filt.share,Gene == gene.sel)$PeakRanges, subset(p2g.filter.share,gene == gene.sel)$peak))

    figR_n <- nrow(subset(cisCorr.filt.share,Gene == gene.sel) )
    p2g_n <- nrow(subset(p2g.filter.share,gene == gene.sel))

    cat("shared:",share_n,', figR_method:',figR_n,' p2g_method:',p2g_n,'\n')

}

getLinks_gene(gene.sel = 'SH3TC2') 
#shared: 2 , figR_method: 2  p2g_method: 2




#compare all genes






#####
