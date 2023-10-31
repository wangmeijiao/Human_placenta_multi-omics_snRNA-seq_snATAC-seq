
##################


##cellpairing with iNMF mat -> doDORC -> doGRN -> plot driver TFs and build network

####result table:
##peak to gene linking
#cisCorr.filt : peak to gene link table
#dorcTab.split: split cisCorr.filt by gene and only kept dorcGenes

##tf mining 
#figR.d : tf scaning result table




#devtools::install_github("buenrostrolab/FigR")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#FigR_stim_tutorial.html
#workingDir <- "/home/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/liger/figR_Kartha/FigR/example_independent"
# stimZip <- "https://s3.us-east-1.amazonaws.com/vkartha/FigR/FigR_stim.zi
# #tfname <- 'ESRRA'p"
# download.file(url = stimZip,
#               destfile = paste0(workingDir,basename(stimZip)))
# unzip(paste0(workingDir,basename(stimZip)),exdir = workingDir,overwrite = FALSE)

# setwd(paste0(workingDir))



library(doParallel)
library(BuenColors)
library(FigR) #0.1.0
#library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

library(tidyr)
library(ComplexHeatmap)
library(pheatmap)


library(ggplot2)

library(igraph)
library(ggraph)
library(graphlayouts)


color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" )

color_good <- c("#E7D654", "#6F1482" ,"#DC7035", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,"#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,
 "#CA362E" ,"#2B3918","#1E1E1E" )


#############ArchR colors############
#hmcols <- colorRamps::blue2green2red(length(bks) ) #colors
color_peak <- ArchR::paletteContinuous(set = 'solarExtra',n=256,reverse=FALSE)  
color_tfdev = ArchR::paletteContinuous(set = 'blueYellow',n=257,reverse=FALSE)                       
#color_ga <- paletteContinuous(set='solarExtra',n=257,reverse=FALSE) 
#color_ga <- paletteContinuous(set='horizon',n=257,reverse=FALSE)                      
#color_ga <- paletteContinuous(set='horizonExtra',n=257,reverse=FALSE)  #good        
#color_ga <- paletteContinuous(set='greenBlue',n=257,reverse=FALSE)
#color_ga <- paletteContinuous(set='blueYellow',n=257,reverse=FALSE)
color_ga <- ArchR::paletteContinuous(set='greyMagma',n=257,reverse=FALSE)


color_rna <- colorRampPalette(c('grey','red'))(10) 



color_gradient_my <- c(
    rgb(5,48,97,maxColorValue = 255),
    rgb(42,113,178,maxColorValue = 255),
    rgb(147,198,222,maxColorValue = 255),
    rgb(239,243,245,maxColorValue = 255),
    rgb(253,219,199,maxColorValue = 255),
    rgb(214,96,77,maxColorValue = 255),
    rgb(121,5,34,maxColorValue = 255)

)



#########set colorset
color_gradient_biasedBR  <- readRDS('/home/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/chromVAR_TF_specific/full12526/color_gradient_biasedBR.rds')
#plot
barplot(1:length(color_gradient_biasedBR),col=rev(color_gradient_biasedBR) )


color_set_yellowbrick.flat <- readRDS('/home/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/liger/color_set_yellowbrick.flat.rds') 
#saveRDS(color_set_yellowbrick.flat,'complexHeatmap_interactive/color_set_yellowbrick.flat.rds')

colorset1 <- color_set_yellowbrick.flat[['RdYlGn.11']]
colorset2 <- rev(color_set_yellowbrick.flat[['YlGnBu.9']]) #or this
colorset3 <- rev(color_set_yellowbrick.flat[['RdPu.9']]) #use this !
colorset4 <- color_set_yellowbrick.flat[['Spectral.8']]
colorset5 <- color_set_yellowbrick.flat[['RdGy.11']]
colorset5 <- color_set_yellowbrick.flat[['Oranges.5']]
#colorset <- colorset1

colorset <- color_set_yellowbrick.flat[['Reds.9']]

colorset <- colorset5
colorset <- colorset1
barplot(1:length(colorset),col=colorset )

colorset_go <- rev(color_set_yellowbrick.flat[['YlGnBu.9']])
colorset_pathway <- rev(color_set_yellowbrick.flat[['RdPu.9']])

barplot(1:length(colorset_go),col=colorset_go )
barplot(1:length(colorset_pathway),col=colorset_pathway )

##color set done





#ATAC.se <- readRDS("/home/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/chromVAR_TF_specific/se.raw.rds")

#ATAC.se <- readRDS("/home/mjwang/pwdex/placenta_10X_term_combine/02.snapATAC_harmony/chromVAR_TF_specific/TE8334/se.raw.rds")

ATAC.se <-readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_late_combine/02.snapATAC_harmony/chromVAR_TF_specific/se.raw.rds')
175965 x 24692
#299537 x 11093, raw se obj with SnapATAC cluster, umap, smat.dmat, meta etc.

#rename id to placenta.atac_ to match matINMF or GLUE or cca?
#placenta_donor1#AAACGAAAGGAAGGTA-1 to placenta_atac_AAACGAAAGGAAGGTA-9
cellid <- colnames(ATAC.se) 


cellid_mod <- sapply(stringr::str_split(cellid, "donor|#|-", n = 4), function(x){ paste0( 'placenta_atac_',x[3],'-',x[2]  )   } )



# idy <- which(grepl(pattern = "^placenta_donor2",x = cellid    ))
# cellid[idy] <- gsub( pattern = '-1$',replacement = '-2' ,x =  cellid [idy]  )
# cellid  <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = cellid )

# cellid <- paste0('placenta.atac_',cellid)

colnames(ATAC.se) <- cellid_mod 


renameID_atac <- function(cellid=NULL){
    
    
    cellid_bk <- cellid
    cellid_mod <- sapply(stringr::str_split(cellid, "donor|#|-", n = 4), function(x){ paste0( 'placenta_atac_',x[3],'-',x[2]  )   } )


    
    
#     cellid_bk <- cellid
#     #cellid <- colnames(ATAC.se) 
#     idy <- which(grepl(pattern = "^placenta_donor2",x = cellid    ))
#     cellid[idy] <- gsub( pattern = '-1$',replacement = '-2' ,x =  cellid [idy]  )
#     cellid  <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = cellid )

#     cellid <- paste0('placenta.atac_',cellid)

    #colnames(ATAC.se) <- cellid 
    return (data.frame(cellid = cellid_bk, cellid_new = cellid_mod, stringsAsFactors = FALSE))
    
    
}

# #renamed <- renameID_atac(colnames(readRDS("/home/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/chromVAR_TF_specific/se.raw.rds")))
# renamed <- renameID_atac(colnames(readRDS("/home/mjwang/pwdex/placenta_10X_term_combine/02.snapATAC_harmony/chromVAR_TF_specific/TE8334/se.raw.rds")))


# all.equal(renamed$cellid_new, colnames(ATAC.se),check.attributes =FALSE) #TRUE
# #



#RNAmat <- readRDS("/home/mjwang/pwdex/placenta_10X_combine/02.seurat_harmony/exprMat.count.rds")
#24307 x 10198 #Seurat raw count mat

#RNAmat <- readRDS('/home/mjwang/pwdex/placenta_10X_term_combine/02.seurat_harmony/exprMat_raw_count.rds')
#22355 x 9331

RNAmat <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_late_combine/02.seurat_harmony/exprMat.count.rds')
28686  x 23981


##rename PLA-rna-late1_AAACCCAGTGTAAATG-1 to -1 -2 .. -9: placenta_rna_TTTGTTGTCAGACCGC-9

cellid <- colnames(RNAmat)
cellid_mod <- sapply(stringr::str_split(cellid, "late|_|-", n = 6), function(x){ paste0( 'placenta_rna_',x[5],'-',x[4]  )   } )


colnames(RNAmat) <- cellid_mod

#ATAC.se <- readRDS("./control1h_PBMC_atac_SE.rds") #219136 x 5352
#RNAmat <- readRDS("./control1h_PBMC_RNAnorm.rds") #15584 x 3508

#all.equal(RNAmat, readRDS("/home/mjwang/pwdex/placenta_10X_combine/02.seurat_harmony/placenta.final.final.rds")@assays$RNA@counts  )#TRUE, exprMat.count.rds is the counts slot

# Load CCA components (used for pairing)
# This was derived by running CCA (using all cells) between the original ATAC/RNA data using the top variable scATAC-seq gene scores and scRNA-seq gene expression
#CCA_PCs <- readRDS('./control1h_PBMC_atac_rna_CCA_l2.rds')

#CCA_PCs <- readRDS('/home/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/liger/peak_gene_link/data/matINMF.rds')

#GLUE method, harmony corrected
#CCA_PCs <- read.table('/home/mjwang/pwdex/placenta_10X_term_combine/03.snRNA_snATAC/GLUE/x_glue_harmony.txt')
#17665 x 50

CCA_PCs <- read.table('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_late_combine/03.snRNA_snATAC/GLUE/x_glue_harmony.txt')

#rownames(CCA_PCs) <- gsub('^placenta_','placenta.',rownames(CCA_PCs))


dim(CCA_PCs) # ATAC + RNA (rows), n components (columns)
#48673 x 50

#17665 x 50
#21268 x 20
#8860 x 50

table (c(colnames(ATAC.se),colnames(RNAmat))  %in% rownames(CCA_PCs) )
TRUE 
48673 

# TRUE 
# 17665

# FALSE  TRUE 
#    23 21268


isATAC <- grepl("placenta_atac",rownames(CCA_PCs))

table(isATAC)

FALSE  TRUE 
23981 24692

# FALSE  TRUE 
#  9331  8334

# FALSE  TRUE 
# 10194 11074

# FALSE  TRUE 
#  3508  5352 

ATACcells <- rownames(CCA_PCs)[isATAC]
RNAcells <- rownames(CCA_PCs)[!isATAC]


#visualize UMAP of the ATAC-RNA cell co-embedding
#nPCs <- 20 # Num CCA PCs to use when running UMAP / pairing

nPCs <- 50

# Run UMAP
set.seed(123)
umap.out <- uwot::umap(CCA_PCs[,1:nPCs],
                       metric="cosine",
                       n_neighbors=30)

umap.d <- as.data.frame(umap.out)
colnames(umap.d) <- c("UMAP1","UMAP2")
rownames(umap.d) <- rownames(CCA_PCs)
  
umap.d$Assay <- ifelse(isATAC,"ATAC","RNA")

BuenColors::shuf(umap.d) %>% 
  ggplot(aes(UMAP1,UMAP2,color=Assay)) + 
  geom_point(size=0.1) + 
  theme_classic() + 
  scale_color_manual(values = c("cadetblue","darkorange"))+
  guides(colour = guide_legend(override.aes = list(size=3)))


##pairing cells
ATAC_PCs <- CCA_PCs[isATAC,] #24692 × 50
RNA_PCs <- CCA_PCs[!isATAC,] #23981 x 50


#Maximum matching problem may have only 1e+07 - (nrows + ncols + 2) finite entries; 16388489 too many. Set 'options("optmatch_max_problem_size" = Inf)' to disable this check.
options("optmatch_max_problem_size" = Inf)
pairing <- pairCells(ATAC = as.matrix(ATAC_PCs),
                     RNA = as.matrix(RNA_PCs),
                     keepUnique = TRUE)

source('pairCells_new.r')

environment(pairCells_new) <- environment(pairCells) #clone the enviroment
pairing <- pairCells_new(ATAC = as.matrix(ATAC_PCs), #24692
                         RNA = as.matrix(RNA_PCs), #23981
                         keepUnique = TRUE,
                         chunk.size = 5000
                        )
24505 × 3

Running geodesic pairing in after chunking data ..
Number of cells in bigger dataset:  24692 
Number of cells in smaller dataset:  23981 
Difference in cell count between 2 datasets:  711 
Chunking larger dataset to match smaller datset ..
Chunk size n =  5000  cells
Total number of chunks:  5 

Chunk #  1 
No. cells in chunk:  5000 

...

Barcodes in the larger dataset will now be unique ..



#chunkList <- chunk_CCA(CCA_1 = as.matrix(ATAC_PCs), CCA_2 = as.matrix(RNA_PCs), seed = 123,chunkSize = 5000)
#5 chunk list 


#9227 x 3
#10722 × 3
#
sum(duplicated( pairing$ATAC) ) #0 #997 dups
sum(duplicated( pairing$RNA) ) #8239  #0 dups


plotPairs(ATAC = pairing$ATAC,
          RNA=pairing$RNA,
          max.show = 200,
          umap.df = umap.d)


ATAC.se.paired <- ATAC.se[,pairing$ATAC] #24505 #9227 sum(duplicated(colnames(ATAC.se.paired))) 997 duplicated
RNAmat.paired <- RNAmat[,pairing$RNA] #24505 (8239 duplicated)  #9227




#########get atac_cluster vs rna_cluster confusion map (from granja method step 0)############

##add cluster to pairing table

sum(duplicated(pairing$ATAC)) #0 #997
sum(duplicated(pairing$RNA))  #8239 #0


#cluster.df.add #atac se obj, 8334 x 33, see below

cluster.df.add.atac <- readRDS('../03.snRNA_snATAC/GLUE/cluster.df.add.atac.rds') #24962 x 37 #8334 x 33

#all.equal(cluster.df.add.atac,cluster.df.add,check.attributes = FALSE) #TRUE, but rowname diff


cluster.df.add.rna <- readRDS('../03.snRNA_snATAC/GLUE/cluster.df.add.rna.rds') #23981 x 18 #9331 x 22

#all.equal(cluster.df.add.rna[,2:3],umap.d.rna,check.attributes = FALSE) #TRUE, but rowname diff

#rownames(cluster.df.add.rna) <- gsub('^placenta_rna','placenta.rna',rownames(cluster.df.add.rna))
#rownames(cluster.df.add.atac) <- gsub('^placenta_atac','placenta.atac',rownames(cluster.df.add.atac))

table(pairing$ATAC %in%  rownames(cluster.df.add.atac))
 TRUE 
24505

# TRUE 
# 9227

table(pairing$RNA %in%  rownames(cluster.df.add.rna))
TRUE 
24505

# TRUE 
# 9227


cluster.atac <- as.character(cluster.df.add.atac[pairing$ATAC,'cluster'])
cluster.rna <- as.character(cluster.df.add.rna[pairing$RNA,'cluster'])

pairing$cluster.atac <- cluster.atac
pairing$cluster.rna <- cluster.rna

confuse.df <- as.matrix(table(pairing$cluster.atac,pairing$cluster.rna))

cM <- as.matrix(ArchR::confusionMatrix(pairing$cluster.atac, pairing$cluster.rna))

preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

rna  atac
5	2
3	1
1	4
4	5
2	3
8	8
4	6
7	7
6	9

#rna atac
# 8	7
# 5	5
# 10	2
# 1	1
# 4	4
# 6	6
# 3	3
# 9	8

pheatmap::pheatmap(as.data.frame(cM))


options(repr.plot.width=5.5,repr.plot.height=5)
#pheatmap::pheatmap(as.data.frame(log(cM+1))
res.p <- pheatmap::pheatmap(
    as.data.frame(log2(cM+1))[
                      #c('8','2','1','4','5','3','6','7'),
                      #c('9','7','10','2','4','1','3','5','6','8')
                      #c('8','7','6','5','4','3','2','1'),
                      #c('1','2','3','4','5','6','7','8','9','10')
                      c('9','5','6','4','1','3','7','2','8'),
                      c('10','9','11','6','4','1','2','3','7','5','8')
                     ],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_set_yellowbrick.flat[['Oranges.9']],#color_tfdev,
    border_color = 'white',fontsize = 24,
    display_numbers = TRUE,
    fontsize_number = 8
)

pdf('pdfs/rna-atac-match.pdf',width=5.5,height=5)
print(res.p)
dev.off()



##Peak-gene association testing

cisCorr <- runGenePeakcorr(ATAC.se = ATAC.se.paired,
                           RNAmat = RNAmat.paired,
                           genome = "hg38", # One of hg19, mm10 or hg38 
                           nCores = 10, #slow, need many cores, no use if too many
                           p.cut = NULL, # Set this to NULL and we can filter later
                           n_bg = 100)


saveRDS(cisCorr,'cisCorr.rds')
64309 x 5

#46737
#99609
#72438


timetag <- Sys.time() #12:34, ~10h

#Determining DORCs
head(cisCorr)
boxplot(cisCorr$pvalZ,main = 'cisCorr pvalZ')



#cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05) 
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.1) 
#19360 # pvalZ <= 0.1
#12125 # pvalZ <= 0.05

#9038
#19647



#filtered
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
#numDorcs <- cisCorr %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs #max 16 #max 12 #max 9
10624 x 2

#####plot numDorc distribution for cutoff decision

hist(numDorcs$n,breaks = 20);
abline(v=5)

##better include genes of interest
idx <- grep('FLT1',numDorcs$Gene ) #7 #has 7 links, is DORC
idx <- grep('FLT4',numDorcs$Gene ) #8 #has 5
idx <- grep('LVRN',numDorcs$Gene ) #1 #has 6

idx <- grep('PAPPA$',numDorcs$Gene ) #4 #has 15 
idx <- grep('PAPPA2',numDorcs$Gene ) #3 #has 4
idx <- grep('CSH',numDorcs$Gene ) #CSHL1 has 3, CSH1 has 2, CSH2 has 1 #CSHL1 has 5, CSH1 CSH2 has 1 each 
idx <- grep('GH2',numDorcs$Gene ) #3 #has 1

idx <- grep('ERV',numDorcs$Gene ) #ERVFRD-1 has 1, ERVV-1 has 2, ERVV-2	has 3  #ERVFRD-1 has 7, ERVV-2 has 4
idx <- grep('ANX',numDorcs$Gene ) #ANXA1 has 3, ANXA2 has 2, ANXA2 has 2 #ANXA1 has 4, ANXA3 has 3

idx <- grep('SH3TC',numDorcs$Gene ) #SH3TC2 has 9, SH3TC2 has 3 #SH3TC1 has 3, SH3TC2 has 2


##regulator of regulator?

idx <- grep('STAT4',numDorcs$Gene ) #has 7 #has 6
idx <- grep('STAT',numDorcs$Gene ) #STAT5A	3  STAT5B not found #STAT5B has 2, STAT5A not found?

idx <- grep('ESRR',numDorcs$Gene )#ESRRA 1 ESRRB 1  #has 1 links, not DORC
idx <- grep('PPARG',numDorcs$Gene ) #has 2 #has 4
idx <- grep('PPARD',numDorcs$Gene ) #has 3 #has 3
idx <- grep('REL',numDorcs$Gene ) #REL,RELA,RELB each has 2 #RELA/B has 1
idx <- grep('OVOL',numDorcs$Gene ) #OVOL1-AS1 1, OVOL2 1, OVOL3 1 #OVOL1 has 2, OVOL2 has 1, OVOL3 has 1
idx <- grep('GCM',numDorcs$Gene ) #has 4 #not found

idx <- grep('MYCN',numDorcs$Gene ) #has 3 #not found

idx <- grep('SLC19A1',numDorcs$Gene) #has 3

numDorcs[idx,]


####

library(ggplot2)
options(repr.plot.width=7.5,repr.plot.height=7.5)
dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                         #cutoff = 10, # No. sig peaks needed to be called a DORC
                         cutoff = 5,
                         labelTop = 20,
                         returnGeneList = TRUE, # Set this to FALSE for just the plot
                         force=5 #geom_repl_text, repulsion degree
                      )

length(dorcGenes) 
#2082

#377 cutoff 5
#54
#481

numDorcs[numDorcs$n >= 3,'Gene'] #2082 (>=3) #377 (>=5) #137 #481
all.equal(dorcGenes,  unlist(numDorcs[numDorcs$n >= 5,'Gene']), check.attributes = FALSE  )
#TRUE

###let more dorcs been defined
dorcGenes.cutoff3 <- dorcJPlot(dorcTab = cisCorr.filt,
                         #cutoff = 10, # No. sig peaks needed to be called a DORC
                         cutoff = 3,
                         labelTop = 20,
                         returnGeneList = TRUE, # Set this to FALSE for just the plot
                         force=5)

length(dorcGenes.cutoff3)
#2082 cutoff 3
#505 
#2217

#head 100 of cutoff = 3
head(dorcGenes.cutoff3,100)
'SLCO2A1''ADORA1''TMEM229B''CITED4''CLMN''COL27A1''EPB41L1''TNFAIP2''ADCY5''CDK14''CERS4''LAMA3''LCP1''PSTPIP2''PTGFRN''SCCPDH''SLC6A6''ZC3H12C''ACOXL''CPXM2''CRH''DLGAP4''LAMB4''RGCC''ATP11A''FLNB''ITPK1''LMCD1''MME''OSMR''PLEKHH1''POU2F3''SH3TC1''SLC37A1''TMEM74B''TMPRSS6''AGAP1''ANGPT2''AOX1''ATP2B4''ATP8A2''CAB39''CDC42SE2''CEBPB''CROT''FAM167A''FLT4''FRZB''HIVEP3''HK2''HPCAL1''MAP3K13''MTMR4''MYO16-AS1''MYO1B''NDRG2''NPAS2''PLA2G2F''PTPRJ''SLC40A1''SLC4A4''TP73''VPS16''ACPP''ADAMTSL1''ASAP3''ATG7''BCL3''CGNL1''COL4A2''CTDSP1''EPHA1-AS1''ERP44''FLT1''FOLR1''GFPT2''GPR155''GTF2A1L''H2AFY''IER5''IL17REL''INSIG2''ITPR3''KCNK3''KRT37''LRP11''LYN''MRPL1''MUC20''OSMR-AS1''PDCD1LG2''PSG6''PTGES''SLC20A1''SLC25A37''SLC7A1''SNX10''SPEF2''STAT4''SULT2B1'


#'DLGAP4''CBLC''CLMN''OSMR''SLC15A1''ACOXL''ADCY5''CFLAR-AS1''COL4A2''EPB41L1''HM13-AS1''ITPK1''LAMA3''MRPL28''PCBP1''PDCD4-AS1''TCHH''TNFAIP2''ABHD1''ACTN1''ADORA1''ANGPT2''AQP3''ATP12A''CD44''CD59''CITED4''CPXM2''CROT''CYP17A1-AS1''EPHA1-AS1''F5''FXYD5''GSTM2''HIST1H3A''HPCAL1''ISM1''KIF6''LINC01629''MDM2''MYO1B''PAQR9''PLD2''PLEKHH1''RAB11FIP3''RAB30''RGCC''SCCPDH''SLC20A1''SLCO2A1''SREBF2''TGFB1''TMEM212''ZC3HC1''ABCC3''ABCF3''ADAMTSL1''AGMAT''ALX3''ANKDD1A''ANKLE1''AOX1''ATP13A1''ATP13A3''BCO2''BTBD16''C12orf57''C9orf139''C9orf16''CBR3''CCDC129''CCDC30''CD53''CD93''CDC42SE2''CDKL3''CENPA''CERS4''CIRBP''CLASP1''CLIP4''COL16A1''CRIM1''DNMT1''DOCK5''EIF6''EMP1''EVA1B''FAM167A''FLNB''FLT1''GBP3''GCNT2''GDF15''GNA15''GRHL1''HIPK4''HOXB7''HRASLS5''HS6ST1'



numDorcs[numDorcs$n >= 3,'Gene'] #505 #2217
all.equal(dorcGenes.cutoff3,  unlist(numDorcs[numDorcs$n >= 3,'Gene']), check.attributes = FALSE  )
#TRUE



dorcGenes.cutoff4 <- dorcJPlot(dorcTab = cisCorr.filt,
                         #cutoff = 10, # No. sig peaks needed to be called a DORC
                         cutoff = 4,
                         labelTop = 20,
                         returnGeneList = TRUE, # Set this to FALSE for just the plot
                         force=5)

length(dorcGenes.cutoff4) #863 #156

#'DLGAP4''CBLC''CLMN''OSMR''SLC15A1''ACOXL''ADCY5''CFLAR-AS1''COL4A2''EPB41L1''HM13-AS1''ITPK1''LAMA3''MRPL28''PCBP1''PDCD4-AS1''TCHH''TNFAIP2''ABHD1''ACTN1''ADORA1''ANGPT2''AQP3''ATP12A''CD44''CD59''CITED4''CPXM2''CROT''CYP17A1-AS1''EPHA1-AS1''F5''FXYD5''GSTM2''HIST1H3A''HPCAL1''ISM1''KIF6''LINC01629''MDM2''MYO1B''PAQR9''PLD2''PLEKHH1''RAB11FIP3''RAB30''RGCC''SCCPDH''SLC20A1''SLCO2A1''SREBF2''TGFB1''TMEM212''ZC3HC1''ABCC3''ABCF3''ADAMTSL1''AGMAT''ALX3''ANKDD1A''ANKLE1''AOX1''ATP13A1''ATP13A3''BCO2''BTBD16''C12orf57''C9orf139''C9orf16''CBR3''CCDC129''CCDC30''CD53''CD93''CDC42SE2''CDKL3''CENPA''CERS4''CIRBP''CLASP1''CLIP4''COL16A1''CRIM1''DNMT1''DOCK5''EIF6''EMP1''EVA1B''FAM167A''FLNB''FLT1''GBP3''GCNT2''GDF15''GNA15''GRHL1''HIPK4''HOXB7''HRASLS5''HS6ST1''IL1RL2''IQSEC1''KCNK17''KCNN4''KRT8''KRTAP27-1''LHB''LINC00518''LINC00857''LINC00862''LINC01136''LINC01569''LPL''LRP11''MAP3K13''MAP3K4''MAP4K3''MARC2''MFSD2A''MIF-AS1''MPEG1''MRPL52''NCOA4''PCP2''PDLIM4''PFKFB2''PHYHD1''POR''PRKCE''PSTPIP2''PTPRU''QPCTL''RDH12''RMND5B''SAMD14''SCAMP4''SDC1''SH3TC1''SLC11A1''SLC12A9''SLC24A4''SLC25A25''SLC6A6''SNCG''SNHG15''SNX10''STC2''TALDO1''TEPP''TRIB1''TSPO''XDH''ZC2HC1A''ZC3H12C''ZNF441''ZNF701'




grep('PAPPA$',dorcGenes.cutoff3)

subset(cisCorr.filt,Gene == 'PAPPA') #4 (filtered by pvalZ 0.1 or 0.05)

subset(cisCorr,Gene == 'PAPPA') #18

	     Peak	PeakRanges	                Gene	rObs	        pvalZ
30625	167769	chr9:116112455-116113340	PAPPA	0.0263458447	0.21674306
30626	167770	chr9:116115840-116116401	PAPPA	0.0059112476	0.58805273
30627	167771	chr9:116123027-116123855	PAPPA	0.0148325745	0.33360163
30628	167772	chr9:116124287-116125252	PAPPA	0.0245705040	0.16836212
30629	167773	chr9:116136899-116137706	PAPPA	0.0046722513	0.52167696
30630	167774	chr9:116142184-116143260	PAPPA	0.0069349113	0.61428042
30631	167775	chr9:116150043-116150354	PAPPA	0.0024051364	0.72194653
30632	167776	chr9:116152050-116152802	PAPPA	0.0009478672	0.74298422
30633	167777	chr9:116153248-116154200	PAPPA	0.0248800348	0.14931778
30634	167779	chr9:116163398-116163951	PAPPA	0.0098932673	0.55744741
30635	167781	chr9:116166900-116167138	PAPPA	0.0177198513	0.27198216
30636	167782	chr9:116174090-116174823	PAPPA	0.0067172722	0.66314104
30637	167784	chr9:116187327-116188114	PAPPA	0.0244264309	0.04079142
30638	167785	chr9:116189729-116190066	PAPPA	0.0262705941	0.04111042
30639	167786	chr9:116192567-116192918	PAPPA	0.0315067073	0.01043951
30640	167787	chr9:116196222-116196691	PAPPA	0.0222493583	0.10528511
30641	167788	chr9:116199738-116200068	PAPPA	0.0348615416	0.02537932
30642	167789	chr9:116203188-116203491	PAPPA	0.0060215550	0.62848717


# 21751	122261	chr9:116112477-116113448	PAPPA	0.024482389	0.1579252
# 21752	122262	chr9:116115917-116116391	PAPPA	0.004755223	0.3705540
# 21753	122263	chr9:116123226-116123800	PAPPA	0.000607471	0.6246751
# 21754	122264	chr9:116124298-116125218	PAPPA	0.024378202	0.2236505
# 21755	122265	chr9:116136909-116137698	PAPPA	0.008212099	0.3285438
# 21756	122266	chr9:116142296-116143254	PAPPA	0.021306259	0.1637092
# 21757	122268	chr9:116152114-116152666	PAPPA	0.014077567	0.5955555
# 21758	122269	chr9:116153249-116154203	PAPPA	0.020373479	0.2538970
# 21759	122270	chr9:116161073-116161313	PAPPA	0.022888228	0.1125954
# 21760	122271	chr9:116162155-116162872	PAPPA	0.006864127	0.6462784
# 21761	122273	chr9:116174031-116174820	PAPPA	0.009739544	0.5594841
# 21762	122275	chr9:116187321-116188084	PAPPA	0.019388759	0.1360272
# 21763	122277	chr9:116192480-116192920	PAPPA	0.003335690	0.5377771
# 21764	122278	chr9:116196201-116196641	PAPPA	0.018696808	0.1692644


subset(cisCorr.filt,Gene == 'LAMA3') #11
subset(cisCorr,Gene == 'LAMA3') #17


# 7336	52032	chr18:23675873-23676548	LAMA3	0.02683339	0.038059173
# 7337	52035	chr18:23689210-23690507	LAMA3	0.05030469	0.002328512
# 7338	52037	chr18:23695579-23696316	LAMA3	0.02496487	0.037199321
# 7339	52038	chr18:23699149-23699456	LAMA3	0.03605418	0.005045246
# 7340	52039	chr18:23699680-23700057	LAMA3	0.03490288	0.005726948
# 7341	52040	chr18:23701849-23702136	LAMA3	0.04363785	0.001046833


############

saveRDS(dorcGenes,'dorcGenes.rds')
saveRDS(dorcGenes.cutoff3,'dorcGenes.cutoff3.rds')
saveRDS(dorcGenes.cutoff4,'dorcGenes.cutoff4.rds') #use this?



dorcGenes <- readRDS('dorcGenes.rds')
dorcGenes.cutoff3 <- readRDS('dorcGenes.cutoff3.rds')


dorcGenes <- dorcGenes.cutoff3# dorcGenes.cutoff4 

#######get dorc mat by aggregate accessiblity of peak regions of dorc

source('bin/utils.R')
source('getDORCScores.r')

#environment(getDORCScores) <- environment('FigR::getDORCScores')

dorcMat <- getDORCScores(ATAC.se = ATAC.se.paired, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 4)

dim(dorcMat) #16 min , aggregate cell-peak mat to gene_dorc x cell mat by gene-dorc peak table 
#will smooth by default
2082 x 24505

#54 x 9227
#481 x 10722
#100 x 4914

##row names of dorcMat
rownames(dorcMat)[1:100] #total 2081
'ABALON''ABCA1''ABCC10''ABCC3''ABCD3''ABHD12''ABHD3''ABI2''ABI3''ABTB2''ACAA1''ACAD9''ACBD6''ACHE''ACIN1''ACKR1''ACKR2''ACKR3''ACO1''ACOT11''ACOT13''ACOT7''ACOXL''ACPP''ACSBG2''ACSS2''ACTL10''ACTN4''ACVRL1''ADAM12''ADAM15''ADAM9''ADAMTS3''ADAMTS4''ADAMTSL1''ADAMTSL5''ADAT2''ADCY10''ADCY5''ADCY7''ADGRF4''ADGRG5''ADGRL1''ADGRL3''ADH4''ADIRF''ADM''ADORA1''ADPRHL1''AES''AFF1''AGA''AGAP1''AGPAT3''AGRN''AHCY''AHDC1''AHNAK''AHRR''AICDA''AIMP2''AK3''AKAP12''AKAP7''AKR1B15''AKR1D1''AKR7A2''ALAD''ALAS1''ALDOA''ALG1L2''ALG9''ALOX5AP''ALOXE3''ALPP''ALX3''AMBRA1''AMN''AMN1''AMZ2''ANGPT2''ANGPTL4''ANKDD1A''ANKH''ANKLE1''ANKRD1''ANKRD11''ANKRD13A''ANKRD13D''ANKRD23''ANKRD24''ANKRD27''ANKRD28''ANKRD33''ANKRD9''ANO6''ANP32A''ANTXRL''ANXA1''AOC1'

#'ABCC3''ABCF3''ABHD1''ACOXL''ACTN1''ADAMTSL1''ADCY5''ADORA1''AGMAT''ALX3''ANGPT2''ANKDD1A''ANKLE1''AOX1''AQP3''ATP12A''ATP13A1''ATP13A3''BCO2''BTBD16''C12orf57''C9orf139''C9orf16''CBLC''CBR3''CCDC129''CCDC30''CD44''CD53''CD59''CD93''CDC42SE2''CDKL3''CENPA''CERS4''CFLAR-AS1''CIRBP''CITED4''CLASP1''CLIP4''CLMN''COL16A1''COL4A2''CPXM2''CRIM1''CROT''CYP17A1-AS1''DLGAP4''DNMT1''DOCK5''EIF6''EMP1''EPB41L1''EPHA1-AS1''EVA1B''F5''FAM167A''FLNB''FLT1''FXYD5''GBP3''GCNT2''GDF15''GNA15''GRHL1''GSTM2''HIPK4''HIST1H3A''HM13-AS1''HOXB7''HPCAL1''HRASLS5''HS6ST1''IL1RL2''IQSEC1''ISM1''ITPK1''KCNK17''KCNN4''KIF6''KRT8''KRTAP27-1''LAMA3''LHB''LINC00518''LINC00857''LINC00862''LINC01136''LINC01569''LINC01629''LPL''LRP11''MAP3K13''MAP3K4''MAP4K3''MARC2''MDM2''MFSD2A''MIF-AS1''MPEG1''MRPL28''MRPL52''MYO1B''NCOA4''OSMR''PAQR9''PCBP1''PCP2''PDCD4-AS1''PDLIM4''PFKFB2''PHYHD1''PLD2''PLEKHH1''POR''PRKCE''PSTPIP2''PTPRU''QPCTL''RAB11FIP3''RAB30''RDH12''RGCC''RMND5B''SAMD14''SCAMP4''SCCPDH''SDC1''SH3TC1''SLC11A1''SLC12A9''SLC15A1''SLC20A1''SLC24A4''SLC25A25''SLC6A6''SLCO2A1''SNCG''SNHG15''SNX10''SREBF2''STC2''TALDO1''TCHH''TEPP''TGFB1''TMEM212''TNFAIP2''TRIB1''TSPO''XDH''ZC2HC1A''ZC3H12C''ZC3HC1''ZNF441''ZNF701'


# 'ABHD1''ACOXL''ACTN1''ADCY5''ADORA1''ANGPT2''AQP3''ATP12A''CBLC''CD44''CD59''CFLAR-AS1''CITED4''CLMN''COL4A2''CPXM2''CROT''CYP17A1-AS1''DLGAP4''EPB41L1''EPHA1-AS1''F5''FXYD5''GSTM2''HIST1H3A''HM13-AS1''HPCAL1''ISM1''ITPK1''KIF6''LAMA3''LINC01629''MDM2''MRPL28''MYO1B''OSMR''PAQR9''PCBP1''PDCD4-AS1''PLD2''PLEKHH1''RAB11FIP3''RAB30''RGCC''SCCPDH''SLC15A1''SLC20A1''SLCO2A1''SREBF2''TCHH''TGFB1''TMEM212''TNFAIP2''ZC3HC1'


#ABCD3,ABCG1,ABHD17C,ACOT11,ACSS1,ACTG1,ADCY4,ADCY7,ADGRG1,ADHFE1,ADORA1,ADRB1,AES,AGL,AGPAT2,AGPAT4,ALOX5,ALPL,ANKH,ANKRD33,ANKRD33B,ANKRD37,ANKRD44,ANKRD55,ANPEP,ANTXRL,AP1S3,APOC1,ARHGAP18,ARHGAP22,ARHGAP23,ARHGAP24,ARHGAP26-AS1,ARHGAP28,ARHGEF10L,ARHGEF16,ARHGEF28,ARNTL,ARSI,ASAP1,ASB2,ASPH,ASPHD2,ATG7,ATP10D,ATP13A4,ATP2B1,ATP2B4,ATP8A2,ATP8B1,AXL,B3GNT6,BAZ2A,BCAT1,BIN3,BMP1,BMP7,BOC,BRD4,BTBD3,C10orf91,C11orf45,C15orf39,C17orf80,C1QL1,C1QTNF6,CACNA2D3,CALML3-AS1,CAMSAP1,CAPZB,CARD16,CASD1,CASKIN2,CBX4,CCDC113,CD36,CD9,CDH1,CDKN1C,CEL,CELA2B,CERS4,CFL2,CFLAR-AS1,CHMP4B,CHST15,CHSY1,CLEC1A,CLIC3,CLIC5,CMTM7,COL17A1,COL27A1,COL4A2,COMTD1,CORO6,CPD,CPXM2,CSE1L-AS1,CSF2RA,CSHL1,CSPG4,CTDSP1,CTNNAL1,CTNNBIP1,CYB5B,DCAF12,DCBLD1,DENND4C,DGKD,DGKZ,DHRS3,DKK3,DLC1,DLX4,DNMT1,DOCK9,DOT1L,DPEP1,DRAM2,DUSP1,DUSP14,EDAR,EFHD2,EFNA5,EGLN3,EHBP1,ELAVL1,ELF3,ELFN2,EPS8L1,ERVFRD-1,ETV3,F5,FABP3,FADD,FAM167A,FAM184A,FAM20A,FAM225B,FAM71F1,FAM83D,FAM83F,FBLN1,FBN2,FCGR3B,FER1L6-AS2,FGFR2,FGFRL1,FHL2,FLNB,FLT1,FLT4,FMO5,FOXI3,FOXK1,FRAT1,FSTL3,FXYD3,GALNT18,GALNT2,GAS6-AS1,GCNT2,GCOM1,GFPT2,GJC2,GNB1L,GNG7,GPC3,GPD1L,GPR146,GPR17,GPR32,GPR87,GPRC5C,GPSM2,GPX3,GPX4,GRAMD1A,GRAMD4,GRK3,GSN-AS1,GYPC,HACD4,HBM,HDAC11,HDAC5,HERPUD1,HES4,HEY1,HMGCS2,HOPX,HOTTIP,HSD11B2,HSPB1,HTRA1,ICE1,ID2,IL12RB1,IL1A,IL27,IL2RB,ILDR2,INHBA,IQCJ-SCHIP1-AS1,IRF5,ISG20,ISM1,ISM2,ITGA5,ITGA6,ITGB1,ITGB6,ITGB8,ITPK1,JAG1,JAM3,KANK1,KANK4,KCNA7,KCNK3,KCTD17,KIF6,KLHL5,KRCC1,KRT7,KRT80,LAMA3,LAMA4,LAMA5,LAMB1,LAMB4,LAMC2,LCLAT1,LCP1,LDLRAP1,LEP,LGALS2,LGALSL,LHB,LHPP,LIMD1,LINC00474,LINC00482,LINC00882,LINC01095,LINC01119,LINC01136,LINC01267,LINC01270,LINC01307,LINC01411,LINC01509,LPCAT1,LPL,LRP5,LURAP1,LVRN,LY6E,LYL1,LZTS1,LZTS1-AS1,MAB21L3,MAFA-AS1,MAFK,MAPRE3,MAST1,MCF2L2,MDFI,MET,MFSD12,MFSD2A,MGAT3,MGST1,MICAL2,MIDN,MINK1,MISP,MKLN1,MKNK1,MME,MOCS1,MPEG1,MPP7,MRPL45,MS4A10,MSANTD1,MSI2,MSRB1,MST1R,MTCL1,MTSS1,MTSS1L,MUC4,MYH10,MYO18B,N4BP3,NCMAP,NECTIN3,NEDD4L,NEDD9,NEURL1B,NFATC2,NID1,NIPAL1,NLRC5,NLRP2,NLRP7,NPB,NRP2,NT5E,NUAK2,OAF,OGFRP1,OSBPL9,OTUB2,PADI1,PAGE4,PAPPA,PARD3,PARP1,PARP12,PCDH1,PDGFA,PDLIM1,PEAR1,PEBP4,PEG10,PFKP,PHACTR2,PHF20,PIM3,PITX2,PKD2L1,PKP1,PLA2G2F,PLA2G4D,PLAC4,PLCG2,PLEC,PLEKHG3,PLEKHG6,PLIN5,PLXND1,POMP,POU2F3,PPL,PPP1R14C,PPP2R2B,PPP5C,PROSER2-AS1,PRR5,PRSS12,PSD4,PSTPIP2,PVT1,RAB17,RALB,RALBP1,RASGEF1B,RASSF3,RBBP6,RBPMS,RELT,RHBDL3,RHPN2,RNASE9,RNF144B,RNF165,RNF19B,RNF43,RPL32,RPS6KA1,SCIN,SCRN2,SDC1,SELL,SEMA3F,SEMA5B,SFTPB,SGPP2,SGSM1,SH2D5,SH3BP2,SH3BP5,SIGLEC6,SIPA1L3,SLC12A9,SLC13A4,SLC15A1,SLC19A1,SLC19A3,SLC1A6,SLC22A11,SLC22A5,SLC24A3,SLC26A7,SLC29A1,SLC2A1,SLC2A3,SLC2A5,SLC35B4,SLC40A1,SLC43A2,SLC48A1,SLC6A6,SLC7A1,SLC9A5,SLCO2A1,SLCO4A1,SMAD7,SMAGP,SMG6,SMYD2,SPR,SPSB1,SRGAP1,SRSF4,SSBP3,STARD8,STAT4,STON1-GTF2A1L,SYDE1,SYNPO2L,TAF4B,TC2N,TCF7,TCF7L2,TCHH,TEX33,TFDP2,TGFB1,THRB,TIMP2,TM4SF1,TMEM150C,TMEM189-UBE2V1,TMEM229B,TMEM236,TMEM245,TMEM45A,TMEM74B,TMPRSS13,TMPRSS3,TNFAIP2,TNFRSF1B,TNFSF10,TNIK,TNRC18,TOMM40L,TP53I11,TPST2,TRIB1,TRIM40,TRIM71,TSPAN18,TTC19,UBA5,UBAC1,UBE4B,UNC13C,USP31,VAPA,VAV1,VDR,VWA2,VWA5B1,VWCE,WDR25,WLS,WNT6,WNT7A,XDH,YPEL2,ZBTB10,ZFPM1,ZFX,ZNF280B,ZNF395,ZNF431,ZNF488


saveRDS(dorcMat,'dorcMat.rds')


###aggregate docrMat by cluster then hclust (or NMF?)


##plot umap of the paired umap
meta <- metadata(ATAC.se.paired)
names(meta)
#'dr.umap''cluster''smat.dmat''cluster.df.add''meta'
#'dr.umap''cluster''metadata''smat.dmat''cluster.df.add''meta'


dr.umap <- meta$dr.umap
cluster <- meta$cluster
cluster.df.add <- meta$cluster.df.add
smat.dmat <- meta$smat.dmat

rownames(dr.umap) <- renameID_atac(rownames(dr.umap) )$cellid_new
names(cluster) <- renameID_atac( names(cluster) )$cellid_new
rownames(cluster.df.add ) <- renameID_atac(rownames(cluster.df.add ) )$cellid_new
rownames(smat.dmat ) <- renameID_atac(rownames(smat.dmat ) )$cellid_new




#colData(ATAC.se.paired) %>% 
cluster.df.add %>%
  as.data.frame() %>% 
  ggplot(aes(UMAP_1,UMAP_2)) + 
  geom_point(color="gray",size=0.8)+ theme_classic()




#quick plot umap function
quickDimPlot <- function(data = cluster.df.add_TM_ok, feature = 'cluster', title= '8w all'){
  options(repr.plot.height=7.5,repr.plot.width=8.5)
  centers <- data %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))   

  p <- ggplot(data=data, aes_string(x = "UMAP_1", y = "UMAP_2", col = `feature`)) +
  #ggplot(data=cluster.df.add_TM_ok, aes(x = UMAP_1, y = UMAP_2, col = library)) +
  geom_point(size = .1) +
  geom_text(data = centers, 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            size = 10) +
  scale_colour_manual(values = color_good)  +
  guides(col = guide_legend(override.aes = list(size = 6))) + 
  ggtitle(paste(title,', total nuclei: ',nrow(data),sep='') ) +
  #theme_bw() +
  theme(legend.position = 'right',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )
  
  print(p)
  return('done')
}

options(repr.plot.width = 7.5, repr.plot.height = 7.5)
#plot(cluster.df.add[,c('UMAP_1','UMAP_2')])
quickDimPlot(data = cluster.df.add, feature = 'cluster', title= 'late pregnancy atac')



aggregateClusters <- function(idents,ids,data){
   data.res <- data.frame(row.names = row.names(data))
   if(ids == "all"){ids = levels(idents)}
   for(id in ids){
     barcodes.sel <- names(idents[idents == id])
     ##data.sel <- Matrix::rowSums(data[,barcodes.sel])
     #table(barcodes.sel %in% colnames(data))
     barcodes.sel <- barcodes.sel[barcodes.sel %in% colnames(data)]
     data.sel <- Matrix::rowMeans(data[,barcodes.sel]) #use this for tfdev
     data.res <- cbind(data.res,data.sel)
     
   }
   colnames(data.res) <- ids
   return(data.res)
  
}



table(colnames(dorcMat) %in% rownames(cluster.df.add) )
 TRUE 
24505

#9227 TRUE

table(rownames(cluster.df.add)  %in% colnames(dorcMat) )
FALSE  TRUE 
  187 24505


FALSE  TRUE 
  104  8230

sum(duplicated(colnames(dorcMat))) #0 #997


#TRUE 
#10722 of 11093


cluster.atac <- as.character(cluster.df.add$cluster)
names(cluster.atac) <- rownames(cluster.df.add)

table(cluster.atac)
 1    2    3    4    5    6    7    8    9 
3457 3583 3186 5910 2983 2651 1449 1050  423 

cluster.atac <- factor(cluster.atac, levels =  c('9','5','6','4','1','3','7','2','8'))#c('8','2','1','4','5','3','6','7') )

table(cluster.atac)
 9    5    6    4    1    3    7    2    8 
 423 2983 2651 5910 3457 3186 1449 3583 1050


#cluster.atac <- cluster.atac [names(cluster.atac) %in% colnames(dorcMat.s) ]

dorcMat.aggre <- aggregateClusters(idents = cluster.atac,ids='all',data=dorcMat)

options(repr.plot.width = 7.5, repr.plot.height = 25)
pheatmap::pheatmap(dorcMat.aggre,scale = 'row',border_color = 'white',fontsize_row = 12,fontsize_col = 20,show_rownames = TRUE)



##smooth the DORC

lsi <- smat.dmat

#lsi <- readRDS("./control1h_PBMC_atac_lsi.rds")
dim(lsi)
24692 x 25

#8334 x 16
#11093 x 16

stopifnot(all(colnames(dorcMat) %in% rownames(lsi)))
#TRUE

#TRUE 
#9227

# Subset to paired ATAC
lsi <- lsi[colnames(dorcMat),]

dim(lsi)
24505 x 25

#9227 x 16
#10722 x 16

sum(duplicated(rownames(lsi))) #0 #997

# Get cell KNNs
cellkNN <- FNN::get.knn(lsi,k=30)$nn.index #24505 × 30 
rownames(cellkNN) <- colnames(dorcMat)

all.equal(rownames(lsi),rownames(cellkNN) ) #TRUE

# Smooth dorc scores using cell KNNs (k=30)
library(doParallel)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN,mat = dorcMat,nCores = 4) #9 min, running window smooth guided by  knn?
2082 x 24505

#156 x 9227
#481 x 10722


##rename dorcMat.s to dorcMat.s (have duplication rowid)

all.equal(rownames(dorcMat),rownames(dorcMat.s)) #TRUE

#colid <- colnames(dorcMat.s)

#colid_rename <- gsub('\\.{3}[0-9]+$','',colid)

#all.equal(colid_rename,colnames(dorcMat)) #TRUE


#colnames(dorcMat.s) <- colid_rename

sum(duplicated(colnames(dorcMat.s))) #0 #997

all.equal(colnames(dorcMat.s),colnames(ATAC.se.paired) ) #TRUE



# smoothScoresNN <- function (NNmat, mat, geneList = NULL, barcodesList = NULL, nCores = 1) 
# {
#     stopifnot(all.equal(nrow(NNmat), ncol(mat)))
#     if (is.null(rownames(NNmat))) 
#         stop("NN matrix has to have matching cell IDs as rownames\n")
#     if (!all.equal(rownames(NNmat), colnames(mat))) 
#         stop("Nearest-neighbor matrix and cell data matrix don't have matching cells barcodes ..\n")
#     cat("Number of cells in supplied matrix: ", ncol(mat), "\n")
#     cat("Number of genes in supplied matrix: ", nrow(mat), "\n")
#     cat("Number of nearest neighbors being used per cell for smoothing: ", 
#         ncol(NNmat), "\n")
#     if (!is.null(geneList)) {
#         if (!(all(geneList %in% rownames(mat)))) {
#             cat("One or more of the gene names supplied is not present in the matrix provided: \n")
#             cat(geneList[!geneList %in% rownames(mat)], sep = ", ")
#             cat("\n")
#             stop()
#         }
#         cat("Running smoothing only on genes:", geneList, sep = "\n")
#         cat("........\n")
#         mat <- mat[rownames(mat) %in% geneList, ]
#     }
#     else {
#         if (nrow(mat) > 10000) {
#             cat("Running smoothing for all genes in matrix! (n = ", 
#                 nrow(mat), ") This is bound to take more time than querying specific markers ..\n", 
#                 sep = "")
#         }
#     }
#     if (Sys.info()["sysname"] %in% "Windows") {
#         message("Windows OS detected .. Cannot support parallilzation using mclapply for mc.cores > 1")
#         message("Using 1 core instead ..\n")
#         nCores <- 1
#     }
#     opts <- list()
#     pb <- txtProgressBar(min = 0, max = ncol(mat), style = 3)
#     progress <- function(n) setTxtProgressBar(pb, n)
#     opts <- list(progress = progress)
#     time_elapsed <- Sys.time()
#     cl <- parallel::makeCluster(nCores)
#     doSNOW::registerDoSNOW(cl)
#     if (!is.null(barcodesList)) {
#         cat("Subsetting to ", length(barcodesList), " barcodes in dataset..\n")
#         NNmat <- NNmat[barcodesList, ]
#     }
#     cat("Running in parallel using ", nCores, "cores ..\n")
#     matL <- foreach::foreach(x = 1:nrow(NNmat), .options.snow = opts, 
#         .packages = c("Matrix", "data.table", "dplyr")) %dopar% 
#         {
#             smoothedScore <- data.table(Matrix::rowMeans(mat[, 
#                 NNmat[x, ]]))
#             rownames(smoothedScore) <- rownames(mat)
#             colnames(smoothedScore) <- rownames(NNmat)[x]
#             smoothedScore
#         }
#     parallel::stopCluster(cl)
#     close(pb)
#     cat("Merging results ..\n")
#     smoothedMat <- dplyr::bind_cols(matL) %>% data.matrix() %>% 
#         Matrix(sparse = TRUE)
#     rownames(smoothedMat) <- rownames(mat)
#     time_elapsed <- Sys.time() - time_elapsed
#     cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed), 
#         "\n"))
#     return(smoothedMat)
# }



###aggregate docrMat.s by cluster then hclust (or NMF?) again


table(colnames(dorcMat.s) %in% rownames(cluster.df.add) )

 TRUE 
24505

# 9227 with duplications
# TRUE 
# #10722 of 11093


cluster.atac <- as.character(cluster.df.add$cluster)
names(cluster.atac) <- rownames(cluster.df.add)

cluster.atac <- factor(cluster.atac, levels =  c('9','5','6','4','1','3','7','2','8') )#c('8','2','1','4','5','3','6','7') )

table(cluster.atac)
 9    5    6    4    1    3    7    2    8 
 423 2983 2651 5910 3457 3186 1449 3583 1050


cluster.atac <- cluster.atac [names(cluster.atac) %in% colnames(dorcMat.s) ]

all.equal(sort(names(cluster.atac)),sort(colnames(dorcMat.s))) #TRUE

dorcMat.s.aggre <- aggregateClusters(idents = cluster.atac,ids='all',data=dorcMat.s)

options(repr.plot.width = 7.5, repr.plot.height = 25)
#options(repr.plot.width = 7.5, repr.plot.height = 65)
res.dorcmat.s.aggre <- pheatmap::pheatmap(dorcMat.s.aggre,scale = 'row',border_color = 'white',fontsize_row = 12,fontsize_col = 20, show_rownames = TRUE)


##plot dorc for TF gene only ?






# Smooth RNA using cell KNNs
# This takes longer since it's all genes


colnames(RNAmat.paired) <- colnames(ATAC.se.paired) # Just so that the smoothing function doesn't throw an error (matching cell barcodes in the KNN and the matrix)
RNAmat.s <- smoothScoresNN(NNmat = cellkNN,mat = RNAmat.paired,nCores = 4) #17 min
28686 x 24505

##rename RNAmat.s to RNAmat (have duplication rowid)

all.equal(rownames(RNAmat.paired),rownames(RNAmat.s)) #TRUE

# colid <- colnames(RNAmat.s)
# colid_rename <- gsub('\\.{3}[0-9]+$','',colid)
# all.equal(colid_rename,colnames(RNAmat.paired)) #TRUE


# colnames(RNAmat.s) <- colid_rename

sum(duplicated(colnames(RNAmat.s))) #0 #997

all.equal(colnames(RNAmat.s),colnames(ATAC.se.paired) ) #TRUE



# Visualize on pre-computed UMAP
# This is the ATAC UMAP shown in the paper (based on ATAC LSI)

#umap.d <- as.data.frame(colData(ATAC.se.paired)[,c("UMAP1","UMAP2")])

umap.d <- cluster.df.add[,c('UMAP_1','UMAP_2')]
colnames(umap.d) <- c('UMAP_1','UMAP_2')




table(rownames(umap.d) %in% colnames(ATAC.se.paired))
FALSE  TRUE 
  187 24505

# ALSE  TRUE 
#   104  8230

# FALSE  TRUE 
#   371 10722

table(colnames(ATAC.se.paired) %in% rownames(umap.d) )
 TRUE 
24505

# TRUE 
# 9227

# TRUE 
# 10722 

umap.d <- umap.d[colnames(ATAC.se.paired),]
sum(duplicated(rownames(umap.d))) #0

# DORC scores for top DORC(s)

myDORCs <- c('PAPPA','FLT1','CSHL1','STAT4')
#myDORCs <- c('ERVFRD-1','INHBA','DNMT1','CSHL1')
myDORCs <- c('BMP1','LAMA3','STAT4')
myDORCs <- head(dorcGenes,10)
#myDORCs <- c("SEMA7A","IL7R","CD93")

myDORCs <- myDORCs[myDORCs %in% rownames(dorcMat.s) & myDORCs %in% rownames(RNAmat.s)  ]


options(repr.plot.width = 7.5, repr.plot.height = 7.5)
dorcGGlist <- lapply(myDORCs,function(x) { 
  plotMarker2D(umap.d,
               dorcMat.s,
               markers = x,
               maxCutoff = "q0.99",
               colorPalette = "brewer_heat"
  ) + ggtitle(paste0(x," DORC"))
})


# Paired RNA expression for top DORC(s)
# Plot on the same reference ATAC UMAP
rnaGGlist <- lapply(myDORCs,function(x) { 
  plotMarker2D(umap.d,
               RNAmat.s,
               markers = x,
               maxCutoff = "q0.99",
               colorPalette = "brewer_purple"
  ) + ggtitle(paste0(x," RNA"))
})

options(repr.plot.width = 7.5, repr.plot.height = 7.5)
library(patchwork)
(dorcGGlist[[1]] + dorcGGlist[[2]] + dorcGGlist[[3]]) /  (rnaGGlist[[1]] + rnaGGlist[[2]] + rnaGGlist[[3]])



##TF-gene associations

#load the modified runFigRGRN (below code) first!

figR.d <- runFigRGRN(ATAC.se = ATAC.se.paired, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     #genome = "hg19",
                     genome = "hg38",
                     dorcMat = dorcMat.s,
                     dorcK = 5, 
                     rnaMat = RNAmat.s, 
                     nCores = 10)
                     #nCores = 4)


timetag <- Sys.time()
saveRDS(figR.d,'figR.d.rds') 


#######

Getting peak x motif matches ..
Getting peak x motif match positions ..
Determining background peaks ..
Using  50  iterations ..

Testing  948  TFs
Testing  2082  DORCs

# Getting peak x motif matches ..
# Determining background peaks ..
# Using  50  iterations ..

# Testing  860  TFs
# Testing  156  DORCs


# figR.d.new <- runFigRGRN(ATAC.se = ATAC.se.paired, # Must be the same input as used in runGenePeakcorr()
#                      dorcTab = cisCorr.filt, # Filtered peak-gene associations
#                      #genome = "hg19",
#                      genome = "hg38",
#                      dorcMat = dorcMat.s,
#                      dorcK = 5, 
#                      rnaMat = RNAmat.s, 
#                      nCores = 10)

# figR.d.newnew <- runFigRGRN(ATAC.se = ATAC.se.paired, # Must be the same input as used in runGenePeakcorr()
#                      dorcTab = cisCorr.filt, # Filtered peak-gene associations
#                      #genome = "hg19",
#                      genome = "hg38",
#                      dorcMat = dorcMat.s,
#                      dorcK = 5, 
#                      rnaMat = RNAmat.s, 
#                      nCores = 10)


saveRDS(ATAC.se.paired,'ATAC.se.paired.rds')
saveRDS(RNAmat.paired,'RNAmat.paired.rds')
saveRDS(dorcMat.s,'dorcMat.s.rds')
saveRDS(RNAmat.s,'RNAmat.s.rds')
saveRDS(cisCorr.filt,'cisCorr.filt.rds')


# saveRDS(figR.d,'figR.d.rds') 
# saveRDS(figR.d.new,'figR.d.new.rds') #return more intermediate results
#saveRDS(figR.d.newnew,'figR.d.newnew.rds') #add motif position

###
ATAC.se.paired <- readRDS('ATAC.se.paired.rds')
RNAmat.paired <- readRDS('RNAmat.paired.rds')
dorcMat.s <- readRDS('dorcMat.s.rds')
RNAmat.s <- readRDS('RNAmat.s.rds')
cisCorr.filt <- readRDS('cisCorr.filt.rds')


#figR.d <- readRDS('figR.d.rds') 
#figR.d.new <- readRDS('figR.d.new.rds') #return more intermediate results
#figR.d.newnew <- readRDS('figR.d.newnew.rds')

#all.equal(figR.d, figR.d.new$TFenrich.d   ) #TRUE

#all.equal(figR.d, figR.d.newnew$TFenrich.d) #TRUE

#all.equal(figR.d.new$motif_ix, figR.d.newnew$motif_ix) #TRUE

#all.equal(figR.d.new$bg, figR.d.newnew$bg) #TRUE

#all.equal(figR.d.new$DORC.knn, figR.d.newnew$DORC.knn) #TRUE


#all.equal(motif_ix,figR.d.new$motif_ix ) #motif column diffs
#colnames(figR.d.new$motif_ix)
#colnames(motif_ix)
#all.equal(motif_ix[,colnames(figR.d.new$motif_ix)],figR.d.new$motif_ix ) #TRUE

#all.equal(bg, figR.d.new$bg) #TRUE
#all.equal(DORC.knn,figR.d.new$DORC.knn)#TRUE



####################figR functions to create GRN in 'step by step' way#########################

##load this slightly  modified version

runFigRGRN <- function(ATAC.se, # SE of scATAC peak counts. Needed for chromVAR bg peaks etc.
                       dorcK=30, # How many dorc kNNs are we using to pool peaks
                       dorcTab, # peak x DORC connections (should contain indices relative to peaks in ATAC.se)
                       n_bg=50, # No. of background peaks to use for motif enrichment Z test
                       genome, # One of mm10, hg19, hg38, with no default
                       dorcMat, # Expect smoothed
                       rnaMat, # Expect smoothed
                       dorcGenes=NULL, # If only running on a subset of genes
                       nCores=1
){
  # Must be matched data
  stopifnot(all.equal(ncol(dorcMat),ncol(rnaMat)))

  stopifnot(all.equal(colnames(dorcMat),colnames(rnaMat))) #rna cellid substitude by atac cellid after pairing; and used to do correlation of tf expression with dorc overall accessibility (tf-to-peak links)
    
  # Expects "Gene" / "Peak" in dorcTab
  if(!all(c("Peak","Gene") %in% colnames(dorcTab)))
    stop("Expecting fields Peak and Gene in dorcTab data.frame .. see runGenePeakcorr function in BuenRTools")

  if(all(grepl("chr",dorcTab$Peak,ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")

    if(!(all(grepl("chr",rownames(ATAC.se),ignore.case = TRUE))))
      stop("Peak regions provided in dorcTab data.frame but not found as rownames in input SE")

    if(!all(dorcTab$Peak %in% rownames(ATAC.se)))
      stop("Found DORC peak region not present in input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  } else{
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
  # If using index, make sure no indices are outside range of SE
    if(max(dorcTab$Peak) > nrow(ATAC.se))
      stop("Found DORC peak index outside range of input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }


  if(is.null(dorcGenes)) {
    dorcGenes <- rownames(dorcMat)
  } else {
    cat("Using specified list of dorc genes ..\n")
    if (!(all(dorcGenes %in% rownames(dorcMat)))) {
      cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
      cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], sep = ", ")
      cat("\n")
      stop()
    }
  }

  DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))),k = dorcK)$nn.index # Scaled
  rownames(DORC.knn) <- rownames(dorcMat)

  if (is.null(SummarizedExperiment::rowData(ATAC.se)$bias)) {
    if (genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if (genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
  }

  # Set data subfolder path
  packagePath <- find.package("FigR", lib.loc=NULL, quiet = TRUE)

  if(grepl("hg",genome)){
    pwm <- readRDS(paste0(packagePath,"/data/cisBP_human_pfms_2021.rds"))
  } else {
    pwm <- readRDS(paste0(packagePath,"/data/cisBP_mouse_pfms_2021.rds"))
  }

  # Old motif naming convention
  if(all(grepl("_",names(pwm),fixed = TRUE)))
     names(pwm) <- FigR::extractTFNames(names(pwm))

  message("Removing genes with 0 expression across cells ..\n")
  rnaMat <- rnaMat[Matrix::rowSums(rnaMat)!=0,]
  myGeneNames <- gsub(x = rownames(rnaMat),pattern = "-",replacement = "") # NKX2-1 to NKX21 (e.g.)
  rownames(rnaMat) <- myGeneNames

  # Only non-zero expression TFs (also found in rnaMat)
  motifsToKeep <- intersect(names(pwm),myGeneNames)

  # This has to be done on the full SE (same peakset used as input to dorc calling)
  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se,pwms = pwm[motifsToKeep],genome=genome)
  #299935 x 926
  # Keep TFs with some peak x motif match
  motif_ix <- motif_ix[,Matrix::colSums(assay(motif_ix))!=0]
  #TRUE 
  #926
    
  # This has to be done on the full SE (same peakset used as input to dorc calling)
  cat("Getting peak x motif match positions ..\n")
  motif_pos <- motifmatchr::matchMotifs(subject = ATAC.se,pwms = pwm[motifsToKeep],genome=genome, out = 'positions')

    
  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat,1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=ATAC.mat),rowRanges = granges(ATAC.se))
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
    #299537 x 50
  } else {
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }

   
   stopifnot(all.equal(rownames(assay(motif_ix))[dorcTab$Peak], dorcTab$PeakRanges)) #add by mjwang
    
    
  # For each DORC, do motif enrichment among dorc sig Peaks, and correlation of DORC accessibility (smoothed) to TF RNA levels

  cat("Testing ",length(motifsToKeep)," TFs\n")
  cat("Testing ",nrow(dorcMat)," DORCs\n")
  library(doParallel)
  if(nCores > 1)
    message("Running FigR using ",nCores," cores ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(nCores)
  #doSNOW::registerDoSNOW(cl)
  doParallel::registerDoParallel(cl)
  mZtest.list <- foreach(g=dorcGenes,
                         #.options.snow = opts, #modified
                         .packages = c("FigR", "dplyr","Matrix","Rmpfr")) %dopar%   {
                           # Take peaks associated with gene and its k neighbors
                           # Pool and use union for motif enrichment
                           DORCNNpeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% c(g,dorcGenes[DORC.knn[g,]])]) 

#c(g,rownames(dorcMat)[DORC.knn[g,]])]) #in fact, dorcGenes== rownames(dorcMat), see above if statement 
                           DORC.knn_genes <- c(g,dorcGenes[DORC.knn[g,]])
                           DORCNNpeaks_names <- unique(dorcTab$PeakRanges[dorcTab$Gene %in% c(g,dorcGenes[DORC.knn[g,]])])
                           if(usePeakNames)
                             DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks) # Convert to index relative to input

                           mZ <- FigR::motifPeakZtest(peakSet = DORCNNpeaks,
                                                bgPeaks = bg,
                                                tfMat = assay(motif_ix))

                           mZ <- mZ[,c("gene","z_test")]
                           colnames(mZ)[1] <- "Motif"
                           colnames(mZ)[2] <- "Enrichment.Z"
                           mZ$Enrichment.P <- 2*pnorm(abs(mZ$Enrichment.Z),lower.tail = FALSE) # One-tailed
                           mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
                           mZ <- cbind("DORC"=g,mZ)
                           # Correlate smoothed dorc with smoothed expression, with spearman
                           corr.r <- cor(dorcMat[g,],t(as.matrix(rnaMat[mZ$Motif,])),method = "spearman")
                           stopifnot(all.equal(colnames(corr.r),as.character(mZ$Motif)))
                           #stopifnot(all.equal(colnames(corr.r),mZ$Motif)) #modified
      
                           mZ$Corr <- corr.r[1,] # Correlation coefficient
                           mZ$Corr.Z <- scale(mZ$Corr,center = TRUE,scale = TRUE)[,1] # Z-score among all TF correlations
                           mZ$Corr.P <- 2*pnorm(abs(mZ$Corr.Z),lower.tail = FALSE) # One-tailed
                           mZ$Corr.log10P <- sign(mZ$Corr.Z)*-log10(mZ$Corr.P)
                           return(mZ)
                         }
  parallel::stopCluster(cl)

  cat("Finished!\n")
  cat("Merging results ..\n")
  # Merge and save table for downstream filtering and plotting (network)
  TFenrich.d <- do.call('rbind',mZtest.list)
  dim(TFenrich.d)
  rownames(TFenrich.d) <- NULL

  # Make combined score based on multiplication
  # Here, we only sign by corr
  # Since sometimes we lose digit precision (1 - v small number is 1, instead of 0.9999999..)
  # Use Rmpfr, increase precision limits above default (100 here)
  TFenrich.d <- TFenrich.d %>% dplyr::mutate("Score"=sign(Corr)*as.numeric(-log10(1-(1-Rmpfr::mpfr(Enrichment.P,100))*(1-Rmpfr::mpfr(Corr.P,100)))))
  TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
  TFenrich.d
    
  return(list(TFenrich.d= TFenrich.d,  DORC.knn = DORC.knn ,motif_ix=motif_ix, bg=bg, motif_pos = motif_pos )    )
    
}


# motifPeakZtest <- function(peakSet,
#                            bgPeaks,
#                            tfMat
# ) {

#   if(nrow(tfMat)!=nrow(bgPeaks))
#     stop("Reference peak set used for TF and background peaks matrix must match..\n")

#   if(!all(peakSet %in% 1:nrow(bgPeaks)))
#     stop("One or more of the provided peak indices are out of the background peak set range ..\n")


#   # Filter out TFs whose motifs that overlap no peaks (if they exist)
#   tfMat <- tfMat[,Matrix::colSums(tfMat)!=0]

#   # get selected peak motif frequencies
#   cat("Getting selected peak motif frequencies ..\n")

#   # get frequency of motifs in test set (observed)
#   p.tab <- Matrix::colMeans(tfMat[peakSet, ])

#   # get the background frequency in peak sets of the same size
#   cat("Getting background peak motif frequencies ..\n")
#   # extract relevant rows (i.e. peakset being tested) from background peak matrix
#   bg.f <- as.matrix(bgPeaks[peakSet, ])

#   # calculate (background) motif frequencies in each iteration of background peaks corresponding to peakset
#   bg.tab <- apply(bg.f[, c(1:ncol(bgPeaks))], 2, function(bg_iter) {

#     b.i <- Matrix::colMeans(tfMat[bg_iter, ])
#     return(b.i)

#   })

#   cat("Calculating empirical p values and z score p values ..\n")

#   # loop over each motif and generate enrichment statistics compared to background
#   m.p <- dplyr::bind_rows(lapply(names(p.tab), function(motif) {

#     # calculate sd and mean frequencies for bg and selected peaks
#     s <- sd(bg.tab[motif, ])
#     bg_freq <- mean(bg.tab[motif, ])

#     z_score <- (p.tab[motif] - bg_freq) / s

#     if(is.nan(z_score))
#       z_score <- 0

#     # generate data.frame object of relevant statistics
#     d <- data.frame(
#       motifID = motif,
#       gene = extractTFNames(motif),
#       motif_obs_freq = p.tab[motif],
#       motif_bg_freq = mean(bg.tab[motif, ]),
#       motif_counts = p.tab[motif] * length(peakSet),
#       emp_pval = 1 - (sum(bg.tab[motif, ] < p.tab[motif]) / ncol(bg.tab)),
#       z_test = z_score,
#       pval.z = 2 * pnorm(-abs(z_score)),
#       signed.log10p = -log10(2 * pnorm(-abs(z_score))) * sign(z_score)
#     )
#     return(d)
#   }))

#   # sort by enrichment pval, motif observed frequency
#   # Note: this returned an error saying pval.z not found in the arrange, look into it
#   #m.p <- dplyr::arrange(m.p,pval.z, motif_obs_freq)

#   # return df of enrichment scores
#   return(m.p)
# }



###########Visualizing FigR results####################


TFenrich.d= figR.d$TFenrich.d
DORC.knn = figR.d$DORC.knn
motif_ix=figR.d$motif_ix
bg=figR.d$bg
motif_pos = figR.d$motif_pos

figR.d <- TFenrich.d
1973736 × 10
#441558 x 10


#Global regulation profile scatter plot
options(repr.plot.width = 7.5, repr.plot.height = 7.5)
figR.d %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=0.01,shape=16) + 
  
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

ggsave("pdfs/tf_driver_all_scatter_plot.pdf",width = 7.5, height = 7.5)


all.equal(as.character(unique(figR.d$DORC)), rownames(dorcMat) ) #TRUE



##Scatter view
options(repr.plot.width = 7.5, repr.plot.height = 7.5)


plotDrivers(figR.d,score.cut = 1.2,marker = "PAPPA", label = TRUE) #STAT5A, MITF, ATF3?
#plotDrivers(figR.d,score.cut = 1.2,marker = "PAPPA2", label = TRUE)

ggsave("pdfs/tf_driver_PAPPA.pdf",width = 7.5, height = 7.5)


plotDrivers(figR.d,score.cut = 1.2,marker = "LAMA3", label = TRUE)

ggsave("pdfs/tf_driver_LAMA3.pdf",width = 7.5, height = 7.5)


plotDrivers(figR.d,score.cut = 1.2,marker = "CSHL1", label = TRUE)


plotDrivers(figR.d,score.cut = 1.1,marker = "FLT1", label = TRUE) #JUND, BACH1, ZBTB43, FOSL2 (score.cut = 1)
ggsave("pdfs/tf_driver_FLT1.pdf",width = 7.5, height = 7.5)

plotDrivers(figR.d,score.cut = 1.2,marker = "FLT4", label = TRUE) #GRHL1, TWIST1
plotDrivers(figR.d,score.cut = 1.1,marker = "FSTL3", label = TRUE) #MYCN, MNT, ZNF586



plotDrivers(figR.d,score.cut = 1.1,marker = "INHBA", label = TRUE) #HMG20B
plotDrivers(figR.d,score.cut = 1.1,marker = "LVRN", label = TRUE) #POU2F3, RUNX1



# plotDrivers <- function (figR.d, marker, score.cut = 1, label = TRUE) 
# {
#     if (!marker %in% figR.d$DORC) 
#         stop("Marker specified is not a valid DORC symbol found in the data.frame")
#     d <- figR.d %>% filter(DORC %in% marker) %>% mutate(isSig = ifelse(abs(Score) >= 
#         score.cut, "Yes", "No"))
#     if (label) {
#         d$Label <- d$Motif
#         d$Label[d$isSig %in% "No"] <- ""
#     }
#     else {
#         d$Label <- ""
#     }
#     gScatter <- d %>% ggplot(aes(x = Corr.log10P, y = Enrichment.log10P, 
#         color = isSig, label = Label)) + geom_hline(yintercept = 0, 
#         color = "gray60", linetype = "dashed") + geom_vline(xintercept = 0, 
#         color = "gray60", linetype = "dashed") + geom_point(size = 0.8) + 
#         theme_classic() + scale_color_manual(values = c("gray66", 
#         "firebrick3")) + scale_x_continuous(breaks = scales::pretty_breaks()) + 
#         scale_y_continuous(breaks = scales::pretty_breaks()) + 
#         labs(y = "Enrichment log10 P", x = "Correlation log10 P", 
#             title = marker) + ylim(-ceiling(max(abs(d$Enrichment.log10P))), 
#         ceiling(max(abs(d$Enrichment.log10P)))) + xlim(-ceiling(max(abs(d$Corr.log10P))), 
#         ceiling(max(abs(d$Corr.log10P)))) + theme(legend.position = "none", 
#         axis.text = element_text(color = "black"), plot.title = element_text(hjust = 0.5, 
#             face = "italic"), panel.background = element_rect(fill = NA)) + 
#         geom_text(hjust = 1.1, fontface = "italic", color = "black", 
#             size = 3)
#     gScatter
# }




#Ranking TF drivers
options(repr.plot.width = 12.5, repr.plot.height = 7.5)
rankDrivers(figR.d,rankBy = "meanScore")

ggsave("pdfs/tf_driver_rank_by_meanScore.pdf",width = 12.5, height = 7.5)


rankDrivers(figR.d,score.cut = 1.5,rankBy = "nTargets",interactive = TRUE)
rankDrivers(figR.d,score.cut = 1.2,rankBy = "nTargets",interactive = TRUE)




# rankDrivers <- function (figR.d, rankBy = c("meanScore", "nTargets"), myLabels = NULL, 
#     score.cut = NULL, interactive = FALSE) 
# {
#     if (!rankBy %in% c("meanScore", "nTargets")) 
#         stop("rankBy parameter has to be one of meanScore or nTargets to rank drivers using ..\n")
#     if (rankBy %in% "meanScore") {
#         message("Ranking TFs by mean regulation score across all DORCs ..\n")
#         figR.summ <- figR.d %>% group_by(Motif) %>% dplyr::summarise(Score = mean(Score)) %>% 
#             arrange(desc(Score)) %>% mutate(Motif = factor(Motif, 
#             levels = as.character(Motif)))
#         figR.summ$TF <- as.character(figR.summ$Motif)
#         if (is.null(myLabels)) {
#             figR.summ$TF[figR.summ$Score >= quantile(figR.summ$Score, 
#                 0.05) & figR.summ$Score <= quantile(figR.summ$Score, 
#                 0.95)] <- ""
#         }
#         else {
#             figR.summ$TF[!figR.summ$TF %in% myLabels] <- ""
#         }
#         library(ggrepel)
#         gAll <- ggplot(figR.summ, aes(x = Motif, y = Score, label = TF)) + 
#             geom_bar(size = 0.1, stat = "identity", fill = "darkorange", 
#                 color = NA) + theme_classic() + theme(axis.text.x = element_blank(), 
#             axis.text = element_text(color = "black")) + ggrepel::geom_text_repel(size = 3, 
#             min.segment.length = 0.1, segment.size = 0.2, max.overlaps = 20) + 
#             geom_hline(yintercept = 0) + labs(x = "TF Motifs", 
#             y = "Regulation Score")
#     }
#     else {
#         message("Ranking TFs by total number of associated DORCs ..\n")
#         if (is.null(score.cut)) {
#             message("Regulation score cut-off not specified ..\n")
#             score.cut <- 1
#         }
#         message("Using absolute score cut-off of: ", score.cut, 
#             " ..\n")
#         figR.summ <- figR.d %>% filter(abs(Score) >= score.cut) %>% 
#             group_by(Motif) %>% dplyr::select(-DORC) %>% dplyr::summarize(numActivated = sum(Score > 
#             0), numRepressed = sum(Score < 0)) %>% dplyr::mutate(diff = numActivated - 
#             numRepressed) %>% mutate(numActivatedY = ifelse(diff > 
#             0, numActivated, -numActivated), numRepressedY = ifelse(diff > 
#             0, numRepressed, -numRepressed)) %>% dplyr::arrange(desc(diff)) %>% 
#             mutate(Motif = factor(Motif, levels = as.character(Motif))) %>% 
#             dplyr::select(-diff) %>% reshape2::melt(id.vars = c("Motif", 
#             "numActivated", "numRepressed"))
#         gAll <- figR.summ %>% ggplot(aes(x = Motif, y = value, 
#             fill = variable, numActivated = numActivated, numRepressed = numRepressed)) + 
#             geom_bar(stat = "identity", color = "lightgray", 
#                 size = 0.1) + theme_classic() + geom_hline(yintercept = 0) + 
#             scale_fill_manual(values = c("firebrick3", "steelblue4"), 
#                 labels = c("Activated", "Repressed")) + theme(axis.text.x = element_text(angle = 90, 
#             vjust = 0.5, hjust = 1, size = 6), axis.text = element_text(color = "black")) + 
#             labs(x = "Ranked TF Motifs", y = paste0("# Associated genes \nabs(Score) >= ", 
#                 score.cut), fill = "Class") + scale_y_continuous(labels = abs)
#     }
#     if (!interactive) {
#         gAll
#     }
#     else {
#         if (rankBy %in% "meanScore") {
#             plotly::ggplotly(gAll)
#         }
#         else {
#             plotly::ggplotly(gAll + theme(legend.position = "none", 
#                 axis.text.x = element_blank()), tooltip = c("Motif", 
#                 "numActivated", "numRepressed"))
#         }
#     }
# }

#Heatmap view

source('plotfigRHeatmap_new.r')

#options(repr.plot.width = 15, repr.plot.height = 55)
set.seed(123)
options(repr.plot.width = 15, repr.plot.height = 20)
#res.h <- plotfigRHeatmap_new(figR.d = figR.d, #positive only
res.h <- plotfigRHeatmap(figR.d = figR.d,
                #score.cut = 1.2,
                score.cut = 1.5,
                #score.cut = 1.8,
                #score.cut = 2.0,
                #DORCs = genes.to.label[genes.to.label %in% figR.d$DORC],
                #TFs = tf.to.label[tf.to.label %in% figR.d$Motif],
                #column_names_gp = gpar(fontsize=8,fontface = 'bold'), # from ComplexHeatmap
                #row_names_gp = gpar(fontsize=10),
                show_row_dend = FALSE # from ComplexHeatmap
                )

##will arrange the 'DORC' ~ 'TF' column to build a matrix, with value of the last column (reshape2::dcast use reshape2::guess_value to decide the data of the matrix, if not given, use the last column of figR.d)

# DORCsToKeep <- figR.d %>% filter(Score >= score.cut) %>% 
#         pull(DORC) %>% unique()
# TFsToKeep <- figR.d %>% filter(Score >= score.cut) %>% 
#         pull(Motif) %>% unique()

# net.d <- figR.d %>% filter(DORC %in% DORCsToKeep & Motif %in% 
#         TFsToKeep) %>% reshape2::dcast(DORC ~ Motif) %>% tibble::column_to_rownames("DORC") %>% 
#         as.matrix()


#options(repr.plot.width = 15, repr.plot.height = 105)

#options(repr.plot.width = 25, repr.plot.height = 15)
options(repr.plot.width = 7.5, repr.plot.height = 15)


options(repr.plot.width = 9.5, repr.plot.height = 35)
#options(repr.plot.width = 25, repr.plot.height = 15) #for select tf
set.seed(123);draw(res.h,heatmap_legend_side = "right",gap = unit(0, "cm"))



res.h.mat <- res.h@matrix
set.seed(123);row_idx <- row_order(res.h)
set.seed(123);col_idx <- column_order(res.h)

rowid_dorc_list <- rownames(res.h.mat)[row_idx]
colid_tf_list <- colnames(res.h.mat)[col_idx]

regMat <- res.h.mat[row_idx,col_idx]



saveRDS(regMat,'regMat.score_cutoff_2.0.rds')
saveRDS(regMat,'regMat.score_cutoff_1.8.rds')
saveRDS(regMat,'regMat.score_cutoff_1.5.rds')#use this?
saveRDS(regMat,'regMat.score_cutoff_1.2.rds')

saveRDS(regMat,'regMat.score_cutoff_1.8.positive_reg_only.rds')
saveRDS(regMat,'regMat.score_cutoff_1.5.positive_reg_only.rds')
saveRDS(regMat,'regMat.score_cutoff_1.2.positive_reg_only.rds')



all.equal(regMat,readRDS('regMat.score_cutoff_1.5.rds')) #TRUE



regMat <- readRDS('regMat.score_cutoff_2.0.rds')
regMat <- readRDS('regMat.score_cutoff_1.5.rds') #use this?
regMat <- readRDS('regMat.score_cutoff_1.2.rds')

regMat <- readRDS('regMat.score_cutoff_1.5.positive_reg_only.rds')
regMat <- readRDS('regMat.score_cutoff_1.2.positive_reg_only.rds')



################extract matrix and replot with heatmap annotation##################


# #create mat and row labels


# row_sel = c('CPD','ARL6IP5')
# col_sel = c('NFE2L1','ZNF548')

# rowid <- rownames(regMat)
# colid <- colnames(regMat)

# rowid_start <- grep(paste0("^",row_sel[1],"$"), rowid)
# rowid_end <- grep(paste0("^",row_sel[2],"$"), rowid)
# stopifnot(rowid_start < rowid_end)

# row_slice <- rowid_start:rowid_end

# colid_start <- grep(paste0("^",col_sel[1],"$"), colid ) #exact match
# colid_end <- grep(paste0("^",col_sel[2],"$"), colid)
# stopifnot(colid_start < colid_end)

# col_slice <- colid_start :colid_end


# regMat.sel <- regMat[row_slice,col_slice]



# ##column cluster annotation

# colnames(cluster) <- 'Cell type'
# ha_top <- columnAnnotation(
#      df = cluster, #the dataframe, each column will annotated
#      col = list('Cell type'=c(unlist(map_cellcolor_rna),'na'='white') ), #a color list with name
#      #simple_anno_size = unit(0.1, "mm")
#      show_annotation_name = FALSE
# )


# ###left deg group belonging annotation

# DEG_type <- as.data.frame(marker.genes.de[,c('cluster','gene')])
# colnames(DEG_type) <- c('DEGs','gene')

# ##for ID_map add dar type
# #DAR_type.bk <- DAR_type
# #sum(duplicated(DAR_type.bk$peak )) #0
# #rownames(DAR_type.bk) <- DAR_type.bk$peak

# ##DAR_type <- DAR_type[,2,drop=FALSE]

# ha_left <- rowAnnotation(
#      df = DEG_type[,1,drop=FALSE],
#      col = list('DEGs'=c(unlist(map_cellcolor_rna),'na'='white') )
# )



###right gene label annotation
#genes.to.label <- c('MKI67','TEAD4','DNMT1','ERVFRD-1','SH3TC2','LAMA3','PAPPA','CSHL1','FLT1','ENG','MYCN')
# genes.to.label <- c('DNMT1','ITGA6','CDH1',
#                     'ANXA1','ERVFRD-1',
#                     'SH3TC2',
#                     'CGA','PSG8','PSG5','PSG1', 'SEMA3B','EBF2','CSF3R','PDE4D','BACE2','BMP1','CSHL1','PAPPA','LAMA3','CDKN1A','INPP5D','GDF15','FLT1','LVRN','ENG',,,'OSMR','ACDY5','MDM2','CCDC30','PTCHD4','CDKN1A','CROT','CAMTA1','OSMR','GH2','JAK1','BDNF','CSH1','CSH2','CSHL1','BACE2','GDF15','SASH1','MAP3K13','POSTN','CYP19A1','ADAM12','INHA','TGFB1','MAP4K4','AOX1','SLC25A25','XDH')


genes.to.label <- list(#marker gene and deg gene,and XDH-like gene, FLT1-like gene manual selected
       'CTB' = c('DNMT1','CDH1','ITGA6'),#CROT
       'CTB_fusion' = c('ERVFRD-1','ANXA1','ERVV-1','ERVV-2'),
       'STB_nascent' = c('SH3TC2','PDE4D','BACE2','SEMA3B','OSMR','TGFB1','CSF3R'),
       'STB_general' = c('PSG8','PSG5','PSG1','CGA','SDC1'),
       'STB_premature1' = c('BMP1','CSH1','CSH2','CSHL1'),
       'STB_Mature1_LAMA3' = c('LAMA3','HK2','ADCY5','EBF2','TMEM74B'),
       'STB_Mature1' = c('PAPPA','TMEM108','ELL2','ADAMTSL1','XDH','AOX1','SLC25A25','CLIC5','ADAM10','JAK2'), 
       'STB_Mature1_apoptosis' = c('INPP5D','CDKN1A','CCDC30','GDF15','PTCHD4','MDM2','MAP4K4','BDNF'),
       'STB_Mature2_residual' = c('FLT1','CAMTA1','MAP3K13','LVRN','DTNB','ADAM12','PSOTN'),
       'STB_Mature2_apoptosis' = c('DDX60','DDX58','ENG','TP73')
)



genes.to.label.sel <- list(#marker gene and deg gene,and XDH-like gene, FLT1-like gene manual selected
       'CTB' = c('CDH1','ITGA6'),#CROT
       'CTB_fusion' = c('ERVFRD-1'),
       'STB_nascent' = c('SH3TC2','PDE4D','BACE2','SEMA3B','OSMR','TGFB1','CSF3R'),
       'STB_general' = c('PSG8','PSG5','PSG1','CGA','SDC1'),
       'STB_premature1' = c('BMP1','CSH1','CSH2','CSHL1'),
       'STB_Mature1_LAMA3' = c('LAMA3','HK2','EBF2','TMEM74B'),
       'STB_Mature1' = c('PAPPA','TMEM108','ELL2','ADAMTSL1','XDH','CLIC5','ADAM10','JAK2'), 
       'STB_Mature1_apoptosis' = c('INPP5D','CDKN1A','CCDC30','PTCHD4','MDM2','MAP4K4'),
       'STB_Mature2_residual' = c('FLT1','CAMTA1','MAP3K13','LVRN','DTNB','PSOTN'),
       'STB_Mature2_residual_apoptosis' = c('DDX60','DDX58','ENG','TP73')
)


#genes.to.label <- unique(unlist(genes.to.label))
genes.to.label <- unique(unlist(genes.to.label.sel))
#genes.to.label <- unique(unlist(genes.to.label.filter))

# marker.genes.late <- list(
#     #'Quality control' = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.ribo', 'percent.hb', 'percent.chrY', 'percent.xist'),
#     #'Sex gene' = c('XIST','RPS4Y1'),
#     'Trophoblast' = c('KRT7', 'GATA3', 'TFAP2A'),
#     'CTB proliferation' = c('TEAD4','ITGA6','MKI67','TOP2A'),
#     'CTB' = c('DNMT1', 'CDH1', 'PPARG',   'TP53','TP63', 'TP73', 'BCAM'),
#     'CTB fusion' = c('ERVFRD-1','OVOL1','PPARD', 'TEAD3'), 
#     'STB nascent' = c('SH3TC2', 'BACE2','PDE4D','SEMA3B','GCM1', 'ESRRG'), 
#     'STB general' = c('PSG8',  'PSG2', 'PSG5','CGA','SDC1', 'LEP'), 
#     'STB Mature 1' = c('PAPPA', 'ADAMTSL1','ADAMTS6','ANGPTL4', 'GH2', 'JAK2', 'STAT5A', 'STAT5B', 'LAMA3','STAT4','AR', 'VDR'), 
#      'STB Premature 1' = c('BMP1','CSHL1', 'CSH1', 'CSH2','MAFK','EBF2','CSF3R','JUND'),
#     'STB Mature 2?' = c('FLT1', 'ENG', 'ANGPTL4', 'FSTL3', 'INHBA', 'INHA','MYCN', 'POU2F3', 'LVRN', 'TGFB1', 'FOSL2', 'FOS', 'FOSB','JUND', 'JUNB', 'JUN'), 
#      'STB FLT1 residual' = c('FLT1','LVRN','CAMTA1','BAZ2B','GLDC','MAP4K4'),
#     'STB apoptosis' = c('DDX60', 'DDX58', 'SASH1','SERPINE1','JUP','MAP4K4', 'SPATA5', 'GDF15', 'CROT', 'CDKN1A', 'ADCY5'),
#     'STB syncytial-knot-like' = c('CDKN1A','INPP5D','PTCHD4','ADCY5','MDM2','GDF15','CCDC30','CROT','BDNF') #CROT, BDNF two neuropeptide
# #     'EVT' = c('HLA-G', 'LAIR2', 'PLAC8', 'MMP2'), 
# #     'STR general' = c('VIM', 'DLK1'), 
# #     'Vascular Endothelial Cell' = c('PECAM1'), 
# #     'STR' = c('HIVEP3', 'HLA-A', 'HLA-DPA1', 'HLA-DPB1'), 
# #     'Mesenchymal STR'= c('THY1'), 
# #     'Hofbauer Cell'= c('CD68', 'CD14'), 
# #     'Red blood' = c('HBA1', 'HBZ')
# )

#genes.to.label <- unique(c(genes.to.label, unique(unlist(marker.genes.late))  ))


#genes.to.label <- c('GCM1','PPARD','DNMT1','CDH1','BCAM','SH3TC2','CGA','PAPPA','ADAMTSL1','ANGPTL4','GH2','LAMA3','VDR','EBF2','FLT1','INHA','MYCN','POU2F3','TGFB1','MAP4K4','SERPINE1','CROT','CDKN1A','PTCHD4','MDM2','BDNF')


# top.use <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_late_combine/02.seurat_harmony/trajectory_infer/top.use.fullcell.rds') #monocle2 deg ftsb only 
# sum(duplicated(top.use)) #0

# genes.to.label <- top.use


# deg.top100 <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_late_combine/02.seurat_harmony//DEGs/result_do_GO_Pathway_quickonestep_top100_veryloose/geneList.veryloose_use_all.top100.rds')


# genes.to.label <- names(deg.top100[['up']][['c3']])



idx <- vector()
for(i in genes.to.label){
    idx <- c(idx,grep(pattern = paste0('^',i,'$') , x = rownames(regMat) ) )
    
}

idx.id <- rownames(regMat)[idx]

ha_right <- rowAnnotation(foo = anno_mark(at = idx,  #HeatmapAnnotation(... which = 'row')
                                    labels = idx.id,
                                    #labels = as.character(ID_map.use$Gene),
                                    labels_gp = gpar(fontsize = 12),
                                    link_width = unit(2, "mm"),
                                    extend = unit(10, "mm"),
                                    padding = unit(0.8, "mm") #important
                                    #annotation_width= unit(4, "cm")
                                   )

) 


####tf label

tf.to.label <- c(
    'trophoblast' = c('GATA3','TFAP2C','TFAP2A','TFAP2C'),
    'CTB' = c('TP53','TP63','TEAD4','RXRA','PPARG','TCF7L2'),
    'CTB fusion' = c('REL','RELB','RELA','FOSL1','ZNF217','GRHL3'),
    'STB Nascent' = c('GCM1','ESRRG','OVOL1','MAFK','MAFG','MAFF'),
    'STB Mature 1 PAPPA' = c('STAT5A','STAT4','MITF','AR','VDR',
                            'FOSL2','CEBPB','CEBPD','CEBPG','JUND','JUNB','SP3'),
    'STB Mature 1 PAPPA LAMA3' = c('STAT6','EBF2'),
    'STB Mature 2 residual' = c('ARNT','SMAD3','GCM1'),
    'STB Mature 2 residua apoptosis' = c('TP73')
   
)


tf.to.label.sel <- c( #for score cutoff 1.5, both active and repressive regulator
    'trophoblast' = c('GATA3','TFAP2C','TFAP2A','TFAP2C'),
    'CTB' = c('TP53','TP63','TEAD4','RXRA','TCF7L2'),
    'CTB fusion' = c('REL','RELB','RELA','FOSL1','ZNF217','GRHL3'),
    'STB Nascent' = c('GCM1','ESRRG','OVOL1','MAFK','MAFG'),
    'STB Mature 1 PAPPA' = c('STAT5A','STAT4','MITF','AR',
                            'FOSL2','CEBPB','CEBPG','JUND','JUNB','SP3'),
    'STB Mature 1 PAPPA LAMA3' = c('STAT6'),
    'STB Mature 2 residual' = c('ARNT','GCM1')
    #'STB Mature 2 residua apoptosis' = c('TP73')
   
)



#colid.sel <- c('SMAD3','NPAS2','RELB','GLIS2','MAFK',,,)
#tf.to.label <- unique(c(tf.to.label,colid.sel))

#tf.to.label <- unique(unlist(tf.to.label))
tf.to.label <- unique(unlist(tf.to.label.sel))

idx <- vector()
for(i in tf.to.label){
    idx <- c(idx,grep(pattern = paste0('^',i,'$') , x = colnames(regMat) ) )
    
}

idx.id <- colnames(regMat)[idx]

ha_top <- columnAnnotation(foo = anno_mark(at = idx,  #HeatmapAnnotation(... which = 'row')
                                    labels = idx.id,
                                    #labels = as.character(ID_map.use$Gene),
                                    labels_gp = gpar(fontsize = 12),
                                    #link_width = unit(2, "mm"),
                                    #extend = unit(5, "mm"),
                                    #padding = unit(1, "mm") #important
                                    #annotation_width= unit(4, "cm")
                                   )

) 


#color_use <- rev(colorset_pathway)
#color_use <- color_gradient_my
color_use <- color_peak #use this

#color_use <- circlize::colorRamp2(seq(-2, 2, length.out = 9),  #inside plotfigRHeatmap
#        colors = BuenColors::jdb_palette("solar_flare"))

#str(color_use)

res.hp = Heatmap(#regMat.sel, name = "TF-target", 
                 regMat, name = "Score", 
         cluster_rows = FALSE, 
         cluster_columns = FALSE, 
         show_row_names = TRUE,#FALSE,
         show_column_names = TRUE,#FALSE,
         show_row_dend = FALSE,
         show_column_dend = FALSE,
         use_raster = FALSE,#will use raster if >2000 row or cols, however rstudio do not support raster
         ##col = circlize::colorRamp2(seq(-1.5,1.5,by=3/10), viridis(n = 11,option = "C")),
         col = circlize::colorRamp2(seq(-1.5,1.5,by=3/(length(color_use)-1)), color_use),
         #col = color_use,
         #col = circlize::colorRamp2(seq(-1,1,by=2/(length(color_gradient_my)-1)), color_gradient_my),
         #col = circlize::colorRamp2(seq(-1.5,1.5,by=3/(length(color_tfdev1)-1)), color_tfdev1),
         #col = circlize::colorRamp2(seq(-1.5,1.5,by=3/255), color_peak),
         na_col = 'white',
         #column_km = 3,
         #row_km = 3,
         #heatmap_legend_param = list(color_bar = "continuous"),
         #right_annotation = ha#,heatmap_width=unit(8, "cm"),
         ##clustering_distance_rows  = 'pearson', 
         ##clustering_distance_columns  = 'pearson',
         #column_split = cluster[,1],
         column_gap = unit(.1,'cm'),
         row_gap = unit(.1,'cm'),
         #column_labels = levels(cluster[,1]),
         column_names_side = 'bottom',
         top_annotation = ha_top,#trajectory bin cluster belonging
         #bottom_annotation = , #trajectory arrow and text, no use decorate_heatmap_body
         ##left_annotation = ha_left, #peak dar cluster belonging
         right_annotation = ha_right #peak annotation gene label
)


pdf(file = 'pdfs/TF-mining-regulation.heatmap.withlable.withArialMT.pdf',width = 7.5,height = 8.5,useDingbats = FALSE,fonts = NULL) #family to get linux fonts, fonts to get extral fonts file 
#pdf(file = 'DEGs/de.marker.top50.heatmap.nolegend.pdf',width = 5.5,height = 7.5,useDingbats = FALSE)

#"AvantGarde", "Bookman", "Courier", "Helvetica","Helvetica-Narrow", "NewCenturySchoolbook", "Palatino", "Times", Arial is a microsoft font, not availabel in linux


#options(repr.plot.height = 5.5, repr.plot.width = 19)
#options(repr.plot.height = 15, repr.plot.width = 7.5)
options(repr.plot.height = 8.5, repr.plot.width = 7.5)

options(repr.plot.height = 155, repr.plot.width = 35)
draw(res.hp, heatmap_legend_side = "right",gap = unit(0.1, "cm")) 


dev.off()




# res.h.mat <- res.h@matrix

# set.seed(123);row_idx <- row_order(res.h)
# set.seed(123);col_idx <- column_order(res.h)

# rowid_dorc_list <- rownames(res.h.mat)[row_idx]
# colid_tf_list <- colnames(res.h.mat)[col_idx]

# regMat_new <- res.h.mat[row_idx,col_idx]

# all.equal(regMat_new,regMat) #plotfigRHeatmap mat == customize heatmap mat



# #Heatmap for given TF
# plotfigRHeatmap(figR.d = figR.d,
#                 score.cut = 1.5,
#                 TFs = c('CEBPG','VDR','TBX3','ZNF266','MAFF','MYCN'),
#                 column_names_gp = gpar(fontsize=12), # from ComplexHeatmap
#                 #row_names_gp = gpar(fontsize=12),
#                 show_row_dend = FALSE # from ComplexHeatmap
#                 )





# score.cut <- 1.5
# column_names_gp = gpar(fontsize=10)
# show_row_dend = FALSE
# DORCs = NULL
# TFs = NULL

# plotfigRHeatmap <- function (figR.d, score.cut = 1, DORCs = NULL, TFs = NULL, ...) 
# {
#     message("Using absolute score cut-off of: ", score.cut, " ..\n")
#     DORCsToKeep <- figR.d %>% filter(abs(Score) >= score.cut) %>% 
#         pull(DORC) %>% unique()
#     TFsToKeep <- figR.d %>% filter(abs(Score) >= score.cut) %>% 
#         pull(Motif) %>% unique()
#     if (!is.null(DORCs)) {
#         if (!all(DORCs %in% figR.d$DORC)) 
#             stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
#         DORCsToKeep <- intersect(DORCsToKeep, DORCs)
#         TFsToKeep <- figR.d %>% filter(abs(Score) >= score.cut & 
#             DORC %in% DORCsToKeep) %>% pull(Motif) %>% unique()
#     }
#     if (!is.null(TFs)) {
#         if (!all(TFs %in% figR.d$Motif)) 
#             stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
#         TFsToKeep <- intersect(TFsToKeep, TFs)
#         DORCsToKeep <- figR.d %>% filter(abs(Score) >= score.cut & 
#             Motif %in% TFsToKeep) %>% pull(DORC) %>% unique()
#     }
#     net.d <- figR.d %>% filter(DORC %in% DORCsToKeep & Motif %in% 
#         TFsToKeep) %>% reshape2::dcast(DORC ~ Motif) %>% tibble::column_to_rownames("DORC") %>% 
#         as.matrix()
#     message("Plotting ", nrow(net.d), " DORCs x ", ncol(net.d), 
#         "TFs\n")
#     myCols <- circlize::colorRamp2(seq(-2, 2, length.out = 9), 
#         colors = BuenColors::jdb_palette("solar_flare"))
#     myHeat <- ComplexHeatmap::Heatmap(net.d, col = myCols, clustering_distance_rows = "pearson", 
#         clustering_distance_columns = "pearson", name = "Score", 
#         border = TRUE, row_names_gp = gpar(fontsize = 10, fontface = "bold"), 
#         ...)
#     myHeat
# }


# grep('ESRRG',unique(figR.d$Motif),value=TRUE)
# unique(subset(figR.d,Motif == 'ESRRG' & Score > 0 & Enrichment.log10P > 0.5)$DORC )


# peaks_withmotif <- rownames(motif_ix)[as.matrix(assay(motif_ix[,'MITF']))[,1]]
# length(peaks_withmotif)
# 9249

# table(subset(cisCorr.filt,Gene == 'ADCY7')$PeakRanges %in% peaks_withmotif)
# FALSE  TRUE 
#     8     1

# peak_withmotif_idx <- which(subset(cisCorr.filt,Gene == 'ADCY7')$PeakRanges %in% peaks_withmotif)

# subset(cisCorr.filt,Gene == 'ADCY7')[peak_withmotif_idx,]



####plot expression of heatmap row (DORC) and column (TF) gene

exprMat <- readRDS("../02.seurat_harmony/exprMat.data.rds")
exprMat.z <- readRDS("../02.seurat_harmony/exprMat.scale.rds")
exprMat.z.aggre <- readRDS('../02.seurat_harmony/exprMat.ave.z.rds')

# exprMat <- readRDS('/home/mjwang/pwdex/placenta_10X_term_combine/02.seurat_harmony/exprMat_data.rds')
# #22355 x 9331

# exprMat.z <- readRDS('/home/mjwang/pwdex/placenta_10X_term_combine/02.seurat_harmony/exprMat_scale_data.rds')

# exprMat.z.aggre <- readRDS('/home/mjwang/pwdex/placenta_10X_term_combine/02.seurat_harmony/exprMat.z.aggre.rds')
# #22355 x 10

#rename id
colnames(exprMat) <- paste0('placenta.rna_',colnames(exprMat))
colnames(exprMat.z) <- paste0('placenta.rna_',colnames(exprMat.z))


##

gene_list <- rowid_dorc_list #1763 #135
title <- 'DORC'

gene_list <- colid_tf_list #402 #140
title <- 'TF '

gene_list <- gene_list[gene_list %in% rownames(exprMat.z.aggre) ] #135

marker.mat.z <- exprMat.z.aggre[gene_list,]



grep("PAPPA",rowid_dorc_list,value=TRUE) #1139
rowid_dorc_list[1130:1145]


grep("LAMA3",rowid_dorc_list,value=FALSE) #1208

grep("FLT1",rowid_dorc_list,value=FALSE) #949

grep("DNMT1",rowid_dorc_list,value=FALSE) #276

# make_bold_names <- function(mat, rc_fun, rc_names) { #different font boldness and color for rowname/colnames
#   bold_names <- rc_fun(mat)
#   ids <- rc_names %>% match(rc_fun(mat))
#   ids %>%
#     walk(
#       function(i)
#         bold_names[i] <<-
#         bquote(bold(.(rc_fun(mat)[i]))) %>%
#         as.expression()
#     )
#   bold_names
# }


# #gene_list_deg
# #marker.mat.z.sel <- marker.mat.z[gene_list_deg,]


# rowid_hi <- gene_list_deg

# qs.value <- quantile(unlist(marker.mat.z), probs = seq(0,1,0.01) , na.rm = TRUE)
# cutoff.low <- qs.value['5%'] 
# cutoff.high <- qs.value['99%'] 


# ##plot the all target gene in eGRN
# ##options(repr.plot.height = 25, repr.plot.width = 5.5)
# options(repr.plot.height = 20, repr.plot.width = 5)
# #options(repr.plot.width = width, repr.plot.height = height)
# res.p <- pheatmap(
#                marker.mat.z, 
#                cluster_cols = TRUE,
#                cluster_rows = TRUE, 
#                treeheight_col = 2.5,
#                treeheight_row = 2.5,
#                #cellwidth = 20, 
#                #cellheight = 2.8,
#                na_col = 'white', 
#                color = rev(colorset_pathway),
#     #c(colorset_go,'white','white',rev(colorset_pathway)),#colorset,#rev(colorset_go),#rev(colorset_pathway),
#                border =TRUE,
#                border_color = 'black',
#                #labels_col = colid, 
# #                labels_row = make_bold_names(marker.mat.z,
# #                                             rownames, 
# #                                             rowid_hi
# #                                            ),

#                angle_col = 315,
#                #display_numbers = data.df.text,
#                #number_color = 'white',
#                #fontsize_number = 10,
#                fontsize_col = 20,
#                fontsize_row = 12,
#                main = paste(title,' gene expression zscore',sep=''),
#                silent = FALSE,
#                scale = 'row'
#                #legend_breaks = c(2,4,6,8,10),
#                #legend_labels = c(2,4,6,8,10),
#                #legend = TRUE
# #                        height = height,
# #                        width = width,
# #                        filename = paste('result_do_GO_Pathway_quickonestep/ck.GO_BP.pvalue_',
# #                                         pvalue,'.qvalue_', qvalue,'.heatmap.pdf' ,sep = '' 
# #                                        )

#               )



###manually select  regulation score matrix modules#####




sliceModule <- function(mat = regMat, row_sel = c('LPL','MPEG1'), col_sel = c('STAT5A','SHOX'), data = NULL, plot = FALSE ){
    
    rowid <- rownames(mat)
    colid <- colnames(mat)
    
    rowid_start <- grep(paste0("^",row_sel[1],"$"), rowid)
    rowid_end <- grep(paste0("^",row_sel[2],"$"), rowid)
    stopifnot(rowid_start < rowid_end)
    
    row_slice <- rowid_start:rowid_end
    
    colid_start <- grep(paste0("^",col_sel[1],"$"), colid ) #exact match
    colid_end <- grep(paste0("^",col_sel[2],"$"), colid)
    stopifnot(colid_start < colid_end)
    
    col_slice <- colid_start :colid_end
    
    
    ##plot if TRUE
    if(plot){
        gene_list1 <- rowid[row_slice]
        title1 <- 'DORC'

        gene_list2 <- colid[col_slice]
        title2 <- 'TF '

        gene_list1 <- gene_list1[gene_list1 %in% rownames(data) ] #135
        marker.mat1 <- data[gene_list1,]

        gene_list2 <- gene_list2[gene_list2 %in% rownames(data) ] #135
        marker.mat2 <- data[gene_list2,]
        
        stopifnot(length(gene_list1) > 2)
        stopifnot(length(gene_list2) > 2)
        
        #qs.value <- quantile(unlist(marker.mat), probs = seq(0,1,0.01) , na.rm = TRUE)
        #cutoff.low <- qs.value['5%'] 
        #cutoff.high <- qs.value['99%'] 


        ##plot the all target gene in eGRN
        ##options(repr.plot.height = 25, repr.plot.width = 5.5)
        options(repr.plot.height = 20, repr.plot.width = 5.5)
        #options(repr.plot.width = width, repr.plot.height = height)
        res.p1 <- pheatmap(
                       marker.mat1, 
                       cluster_cols = TRUE,
                       cluster_rows = TRUE, 
                       treeheight_col = 2.5,
                       treeheight_row = 2.5,
                       #cellwidth = 20, 
                       #cellheight = 2.8,
                       na_col = 'white', 
                       color = rev(colorset_pathway),
            #c(colorset_go,'white','white',rev(colorset_pathway)),#colorset,#rev(colorset_go),#rev(colorset_pathway),
                       border =TRUE,
                       border_color = 'black',
                       #labels_col = colid, 
        #                labels_row = make_bold_names(marker.mat.z,
        #                                             rownames, 
        #                                             rowid_hi
        #                                            ),

                       angle_col = 315,
                       #display_numbers = data.df.text,
                       #number_color = 'white',
                       #fontsize_number = 10,
                       #fontsize = 10,
                       fontsize_col = 13,
                       fontsize_row = 12,
                       main = paste(title1,' gene expression zscore',sep=''),
                       silent = TRUE,
                       scale = 'row'
                       #legend_breaks = c(2,4,6,8,10),
                       #legend_labels = c(2,4,6,8,10),
                       #legend = TRUE
        #                        height = height,
        #                        width = width,
        #                        filename = paste('result_do_GO_Pathway_quickonestep/ck.GO_BP.pvalue_',
        #                                         pvalue,'.qvalue_', qvalue,'.heatmap.pdf' ,sep = '' 
        #                                        )

                      )
        res.p2 <- pheatmap(
               marker.mat2, 
               cluster_cols = TRUE,
               cluster_rows = TRUE, 
               treeheight_col = 2.5,
               treeheight_row = 2.5,
               #cellwidth = 20, 
               #cellheight = 2.8,
               na_col = 'white', 
               color = rev(colorset_pathway),
    #c(colorset_go,'white','white',rev(colorset_pathway)),#colorset,#rev(colorset_go),#rev(colorset_pathway),
               border =TRUE,
               border_color = 'black',
               #labels_col = colid, 
#                labels_row = make_bold_names(marker.mat.z,
#                                             rownames, 
#                                             rowid_hi
#                                            ),

               angle_col = 315,
               #display_numbers = data.df.text,
               #number_color = 'white',
               #fontsize_number = 10,
               fontsize_col = 13,
               fontsize_row = 12,
               main = paste(title2,' gene expression zscore',sep=''),
               silent = TRUE,
               scale = 'row'
               #legend_breaks = c(2,4,6,8,10),
               #legend_labels = c(2,4,6,8,10),
               #legend = TRUE
#                        height = height,
#                        width = width,
#                        filename = paste('result_do_GO_Pathway_quickonestep/ck.GO_BP.pvalue_',
#                                         pvalue,'.qvalue_', qvalue,'.heatmap.pdf' ,sep = '' 
#                                        )

              )
        
    }
    
    return(list(matSlice = mat[row_slice,col_slice], 
                row_slice=rowid[row_slice],
                col_slice = colid[col_slice],
                res.p1 = res.p1,
                res.p2 = res.p2
               ) 
          )
    
    
}


##module1
module1 <- sliceModule(mat = regMat, row_sel = c('FHL2','ANKRD24'), col_sel = c('NPAS2','CXXC1'), data = exprMat.z.aggre, plot = TRUE )

pdf('plot_network/module1.heatmap.dorc.gene.expr.pdf',width = 4, height = 7.5)
options(repr.plot.width = 4, repr.plot.height = 7.5)
module1$res.p1
dev.off()

pdf('plot_network/module1.heatmap.tf.gene.expr.pdf',width = 4, height = 7.5)
options(repr.plot.width = 4, repr.plot.height = 7.5)
module1$res.p2
dev.off()






##module2
module2 <- sliceModule(mat = regMat, row_sel = c('PTGFRN','LRIG3'), col_sel = c('STAT6','ZNF467'), data = exprMat.z.aggre, plot = TRUE )

pdf('plot_network/module2.heatmap.dorc.gene.expr.pdf',width = 4, height = 7.5)
options(repr.plot.width = 4, repr.plot.height = 7.5)
module2$res.p1
dev.off()

pdf('plot_network/module2.heatmap.tf.gene.expr.pdf',width = 3.5, height = 7.5)
options(repr.plot.width = 3.5, repr.plot.height = 7.5)
module2$res.p2
dev.off()



##module 3
module3 <- sliceModule(mat = regMat, row_sel = c('LHB','LITAF'), col_sel = c('NFIA','PBX1'), data = exprMat.z.aggre, plot = TRUE )

pdf('plot_network/module3.heatmap.dorc.gene.expr.pdf',width = 3.5, height = 12.5)
options(repr.plot.width = 3.5, repr.plot.height =12.5)
module3$res.p1
dev.off()

pdf('plot_network/module3.heatmap.tf.gene.expr.pdf',width = 3.5, height = 3.5)
options(repr.plot.width = 3.5, repr.plot.height = 3.5)
module3$res.p2
dev.off()


##module 4
module4 <- sliceModule(mat = regMat, row_sel = c('TMEM218','CDH1'), col_sel = c('NR1H4','ZNF300'), data = exprMat.z.aggre, plot = TRUE )

pdf('plot_network/module4.heatmap.dorc.gene.expr.pdf',width = 3.8, height = 3.5)
options(repr.plot.width = 3.8, repr.plot.height = 5.0)
module4$res.p1
dev.off()

pdf('plot_network/module4.heatmap.tf.gene.expr.pdf',width = 3.5, height = 3.5)
options(repr.plot.width = 3.5, repr.plot.height = 3.5)
module4$res.p2
dev.off()

##module 5
module5 <- sliceModule(mat = regMat, row_sel = c('TRMT61B','BHLHE40'), col_sel = c('PPARG','ETV5'), data = exprMat.z.aggre, plot = TRUE )

pdf('plot_network/module5.heatmap.dorc.gene.expr.pdf',width = 3.9, height = 3.5)
options(repr.plot.width = 3.9, repr.plot.height = 15.5)
module5$res.p1
dev.off()

pdf('plot_network/module5.heatmap.tf.gene.expr.pdf',width = 3.5, height = 3.5)
options(repr.plot.width = 3.5, repr.plot.height = 4.5)
module5$res.p2
dev.off()

##module 6
module6 <- sliceModule(mat = regMat, row_sel = c('FREM2','CMC2'), col_sel = c('ZSCAN1','REST'), data = exprMat.z.aggre, plot = TRUE )

pdf('plot_network/module6.heatmap.dorc.gene.expr.pdf',width = 3.9, height = 4.5)
options(repr.plot.width = 3.9, repr.plot.height = 4.5)
module6$res.p1
dev.off()

pdf('plot_network/module6.heatmap.tf.gene.expr.pdf',width = 3.6, height = 3.5)
options(repr.plot.width = 3.6, repr.plot.height = 3.5)
module6$res.p2
dev.off()

##module 7
module7 <- sliceModule(mat = regMat, row_sel = c('CFAP45','M1AP'), col_sel = c('JUNB','JUND'), data = exprMat.z.aggre, plot = TRUE )

pdf('plot_network/module7.heatmap.dorc.gene.expr.pdf',width = 3.9, height = 4.5)
options(repr.plot.width = 3.9, repr.plot.height = 4.5)
module7$res.p1
dev.off()

pdf('plot_network/module7.heatmap.tf.gene.expr.pdf',width = 3.6, height = 3.5)
options(repr.plot.width = 3.6, repr.plot.height = 3.5)
module7$res.p2
dev.off()



saveRDS(module1,'module1.rds')
saveRDS(module2,'module2.rds')
saveRDS(module3,'module3.rds')
saveRDS(module4,'module4.rds')
saveRDS(module5,'module5.rds')
saveRDS(module6,'module6.rds')
saveRDS(module6,'module7.rds')

##plot tf regulator and target gene expression correlation: gene DORC, gene real expression, gene imputated expression

#realRNAmat.s <- readRDS('../02.seurat_harmony/knn_smooth_exprMat/RNAmat.s.rds')
realRNAmat.s <- readRDS('../02.seurat_harmony/placenta.smooth.rds')@assays$RNA@data
28686 x 23981

#realRNAmat.s <- realRNAmat.s@assays$RNA@data

#umap.d.rna <- readRDS('../02.seurat_harmony/knn_smooth_exprMat/umap.d.rna.rds')
cluster.df.add <- readRDS('../02.seurat_harmony/cluster.df.add.cstb.rds')
23981 × 18

all.equal(colnames(realRNAmat.s),rownames(cluster.df.add)) #TRUE

umap.d.rna <- cluster.df.add[,c('UMAP_1','UMAP_2')]




# myDORCs <- c('STAT4','STAT5A','EBF2','LAMA3') #3 TF + 1 target , use this?
# myDORCs <- c('LAMA3','SLC25A25','AOX1') #3 DORCs (in module1  network) include LAMA3, use this ?
# #myDORCs <- c('AHRR','MLX','TP73','JUND')
# myDORCs <- c('GCM1','CEBPB','ADCY5','JUND')
# myDORCs <- c('OSMR','SLC6A6','FAM167A','TGFB1')
# myDORCs <- c('LAMA3','ADCY5','OSMR','STAT5A')



myDORCs <- c('STAT5A','STAT4','PAPPA','XDH') #module1
myDORCs <- c('STAT6','AR','LAMA3','ADAMTSL1') #module2



table(myDORCs %in% rownames(dorcMat.s))
table(myDORCs %in% rownames(RNAmat.s))
table(myDORCs %in% rownames(realRNAmat.s))

myDORCs_atac <- myDORCs[myDORCs %in% rownames(dorcMat.s)]
myDORCs_rna <- myDORCs[myDORCs %in% rownames(RNAmat.s)]
myDORCs_realrna <- myDORCs[myDORCs %in% rownames(realRNAmat.s)]



#########plot dorc on atac umap
#options(repr.plot.width = 7.5, repr.plot.height = 7.5)
dorcGGlist <- lapply(myDORCs_atac,function(x) { 
  plotMarker2D(umap.d,
               dorcMat.s,
               markers = x,
               minCutoff = 'q0.3',
               maxCutoff = "q0.99",
               colorPalette = "brewer_heat",
               legend.position = 'right'
  ) + ggtitle(paste0(x," DORC")) + theme(plot.title=element_text(size = 18))
})

#options(repr.plot.width = 15, repr.plot.height = 7.5)
#dorcGGlist[[1]] + dorcGGlist[[2]]
#dorcGGlist[[3]] + dorcGGlist[[4]]


###########Paired RNA expression for top DORC(s) (impute to ATAC space?)
# Plot on the same reference ATAC UMAP
rnaGGlist <- lapply(myDORCs_rna,function(x) { 
  plotMarker2D(umap.d, #24505
               RNAmat.s,
               markers = x,
               #minCutoff = 'q0.2',
               maxCutoff = "q0.99",
               colorPalette = "brewer_purple",
               legend.position = 'right'
  ) + ggtitle(paste0(x," RNA imputation")) + theme(plot.title=element_text(size = 18))
})

#options(repr.plot.width = 15, repr.plot.height = 7.5)
#rnaGGlist[[1]] + rnaGGlist[[2]]
#rnaGGlist[[3]] + rnaGGlist[[4]]

#######real RNA expression for given DORC(s) (load from work_pretty_plot_knn_smooth.r saved rds)
# Plot on real RNA UMAP

#all.equal(rownames(umap.d.rna),colnames(realRNAmat.s)) #TRUE

realGGlist <- lapply(myDORCs_rna,function(x) { 
  plotMarker2D(umap.d.rna,
               realRNAmat.s,
               markers = x,
               minCutoff = 'q0.3',
               maxCutoff = "q0.99",
               colorPalette = "brewer_green",
               legend.position = 'right'
  ) + ggtitle(paste0(x," RNA expression")) + theme(plot.title=element_text(size = 18))
})

#options(repr.plot.width = 15, repr.plot.height = 15)
#(realGGlist[[1]] + realGGlist[[2]]) / (realGGlist[[3]] + realGGlist[[4]])


# ###plot TF imputation RNA##

# myTFs <- c('STAT5A','STAT4','EBF2','MGA','PPARG','TCF7L2')



# tfGGlist <- lapply(myTFs,function(x) { 
#   plotMarker2D(umap.d,
#                RNAmat.s,
#                markers = x,
#                #minCutoff = 'q0.2',
#                maxCutoff = "q0.99",
#                colorPalette = "brewer_purple",
#                legend.position = 'right'
#   ) + ggtitle(paste0(x," RNA pairing ATAC")) + theme(plot.title=element_text(size = 18))
# })

# options(repr.plot.width = 15, repr.plot.height = 15)
# (tfGGlist[[1]] + tfGGlist[[2]]) / (tfGGlist[[3]] + tfGGlist[[4]]) /
# (tfGGlist[[5]] + tfGGlist[[6]])

###arranged to plot 

# options(repr.plot.width = 15, repr.plot.height = 21.0)
# (dorcGGlist[[1]] + rnaGGlist[[1]] + realGGlist[[1]]) / 
# (dorcGGlist[[2]] + rnaGGlist[[2]] + realGGlist[[2]]) / 
# (dorcGGlist[[3]] + rnaGGlist[[3]] + realGGlist[[3]]) / 
# (dorcGGlist[[4]] + rnaGGlist[[4]] + realGGlist[[4]]) 

##3 x 4 DORCs
options(repr.plot.width = 15, repr.plot.height = 5*length(dorcGGlist))
(dorcGGlist[[1]] + rnaGGlist[[1]] + realGGlist[[1]]) / 
(dorcGGlist[[2]] + rnaGGlist[[2]] + realGGlist[[2]]) /
(dorcGGlist[[3]] + rnaGGlist[[3]] + realGGlist[[3]]) /
(dorcGGlist[[4]] + rnaGGlist[[4]] + realGGlist[[4]])

ggsave('plot_network/tf-target.genes.module1.4_x_3.pdf',width = 15, height = 20)
ggsave('plot_network/tf-target.genes.module2.4_x_3.pdf',width = 15, height = 20)



# ##3TFs + 1 DORC
# options(repr.plot.width = 10, repr.plot.height = 10)
# (rnaGGlist[[1]] + rnaGGlist[[2]]) / (rnaGGlist[[3]] + rnaGGlist[[4]] )
# ggsave('plot_network/dorc.genes.TF_with_three_target.imputateRNA.pdf',width = 10, height = 10)


# options(repr.plot.width = 10, repr.plot.height = 10)
# (realGGlist[[1]] + realGGlist[[2]]) / (realGGlist[[3]] + realGGlist[[4]] )
# ggsave('plot_network/dorc.genes.TF_with_three_target.expression.pdf',width = 10, height = 10)





####################Network view################
library(networkD3)

res.network <- plotfigRNetwork(figR.d,
              score.cut = 1.5,
              weight.edges = TRUE) %>%
saveNetwork('test.score.cut1.5.html')


plotfigRNetwork(figR.d,
              score.cut = 1.4,
              weight.edges = TRUE) %>%
saveNetwork('test.score.cut1.4.html')

plotfigRNetwork(figR.d,
              score.cut = 1.2,
              weight.edges = TRUE) %>%
saveNetwork('test.score.cut1.2.html')


plotfigRNetwork(subset(figR.d,Motif == 'ESRRG' & abs(Score) > 0.2),
              score.cut = 0.2,
              weight.edges = TRUE)


# score.cut = 1.4
# score.cut = 1.2
# weight.edges = TRUE
# DORCs = NULL
# TFs = NULL

# plotfigRNetwork <- function (figR.d, score.cut = 1.2, DORCs = NULL, TFs = NULL, weight.edges = FALSE) 
# {
#     net.dat <- figR.d %>% filter(abs(Score) >= score.cut)
#     if (!is.null(DORCs)) 
#         net.dat <- net.dat %>% filter(DORC %in% DORCs)
#     if (!is.null(TFs)) 
#         net.dat <- net.dat %>% filter(Motif %in% TFs)
#     net.dat$Motif <- paste0(net.dat$Motif, ".")
#     net.dat$DORC <- paste0(net.dat$DORC)
#     dorcs <- data.frame(name = unique(net.dat$DORC), group = "DORC", 
#         size = 8)
#     tfs <- data.frame(name = unique(net.dat$Motif), group = "TF", 
#         size = 3)
#     nodes <- rbind(dorcs, tfs)
#     edges <- as.data.frame(net.dat)
#     links <- data.frame(source = unlist(lapply(edges$Motif, function(x) {
#         which(nodes$name == x) - 1
#     })), target = unlist(lapply(edges$DORC, function(x) {
#         which(nodes$name == x) - 1
#     })), corr = edges$Corr, enrichment = edges$Enrichment.P)
#     links$Value <- scales::rescale(edges$Score) * 20
#     colors <- c("Red", "Orange", "Yellow", "Green", "Blue", "Purple", 
#         "Tomato", "Forest Green", "Sky Blue", "Gray", "Steelblue3", 
#         "Firebrick2", "Brown")
#     nodeColorMap <- data.frame(color = colors, hex = gplots::col2hex(colors))
#     getColors <- function(tfColor, dorcColor = NULL) {
#         temp <- c(as.character(nodeColorMap[nodeColorMap$color == 
#             tfColor, ]$hex), as.character(nodeColorMap[nodeColorMap$color == 
#             dorcColor, ]$hex))
#         if (is.null(dorcColor)) {
#             temp <- temp[1]
#         }
#         colors <- paste(temp, collapse = "\", \"")
#         colorJS <- paste("d3.scaleOrdinal([\"", colors, "\"])")
#         colorJS
#     }
#     res.network <- networkD3::forceNetwork(Links = links, Nodes = nodes, Source = "target", 
#         Target = "source", NodeID = "name", Group = "group", 
#         Value = "Value", Nodesize = "size", radiusCalculation = "Math.sqrt(d.nodesize)*2", 
#         arrows = FALSE, opacityNoHover = 0.6, opacity = 1, zoom = TRUE, 
#         bounded = TRUE, charge = -15, fontSize = 13, legend = TRUE, 
#         colourScale = getColors(tfColor = "Tomato", dorcColor = "Sky Blue"), 
#         linkColour = ifelse(links$corr > 0, as.character(nodeColorMap[nodeColorMap$color == 
#             "Forest Green", ]$hex), as.character(nodeColorMap[nodeColorMap$color == 
#             "Purple", ]$hex)))
    
#     ##output TF-target table for cytoscape
#     stopifnot(all.equal(links$corr, edges$Corr))
#     stopifnot(all.equal(links$enrichment, edges$Enrichment.P))
#     write.table(edges[,c('Motif','DORC','Corr','Enrichment.P')],file = 'cytoscape/dorc_tf_enrichment.table.txt',sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)
    
    
    
# }

saveNetwork(network = res.network,'test.score.cut1.4.new.html')




####################quick and simple igraph tf-target network (code from quick_igraph_regulon.r)###########
###################extract node and edge table for cytoscape ?############


mapping_pair = list(
    #id2:rna      id1:ata
    '11' = '9',
    '6' = '5',    
    '4' = '6',
    '1' = '4',
    '3' = '1',
    '2' = '3',
    '7' = '7',
    '5' = '2',
    '8' = '8'
    
# #rna      atac
# '8' = '7',
# '5' = '5',
# '10' = '2',
# '1' = '1',
# '4' = '4',
# '6' = '6',
# '3' = '3',
# '9' = '8',
# '2' = '2'
)

figR.d.bk <- figR.d


getNetworkDataList <- function(figR.d, score.cut = 1.2, DORCs = NULL, TFs = NULL,save=NULL,data.rna = NULL, cid.rna = NULL,data.dorc=NULL, cid.dorc = NULL, module = NULL){
    ##copy from plotfigRNetwork
    #net.dat <- figR.d %>% filter(abs(Score) >= score.cut) #both positive and negative regulation

    net.dat <- figR.d %>% filter(Score >= score.cut) #only positive regulation
    
    if (!is.null(TFs))  #filter all tf target, not only in tf-target heatmap row
        net.dat <- net.dat %>% filter(Motif %in% TFs)  
#     if (!is.null(DORCs)) 
#         net.dat <- net.dat %>% filter(DORC %in% DORCs)

#     if (!is.null(TFs) & !is.null(DORCs)) 
#         net.dat <- net.dat %>% filter(Motif %in% TFs & DORC %in% DORCs) 
    
 
#     if (!is.null(DORCs) & !is.null(TFs)) 
#         net.dat <- net.dat %>% filter(DORC %in% DORCs | Motif %in% TFs)
    
    #net.dat$Motif <- paste0(net.dat$Motif, ".")
    net.dat$Motif <- paste0(net.dat$Motif)
    net.dat$DORC <- paste0(net.dat$DORC)
    
    tflist.tab <- split(x = net.dat, f = net.dat$Motif)
    
    stopifnot(sum(table(names(tflist.tab) %in% TFs)) == length(tflist.tab))
    
    
    ##1 for plorTFgraphList
    tflist <- lapply(tflist.tab, FUN = function(x){ x$DORC  } )
    
    
    #2 for cytoscape (simple but pretty network with two layers)??
    
    ##node table
    dorcs <- data.frame(name = unique(net.dat$DORC), group = "DORC", 
        size = 3)
    tfs <- data.frame(name = unique(net.dat$Motif), group = "TF", 
        size = 8)
    nodes <- rbind(dorcs, tfs)

    ##add node data (tf and dorc gene expr 
    stopifnot(cid.rna %in% colnames(data.rna))
    stopifnot(cid.dorc %in% colnames(data.dorc))
        
    value.rna <- data.rna[as.character(nodes$name),cid.rna] #color?
    value.dorc <- data.dorc[as.character(nodes$name),cid.dorc] #size?
    
    stopifnot(length(value.rna) == nrow(nodes) )
    stopifnot(length(value.dorc) == nrow(nodes) )
    
    nodes$value.rna <- value.rna
    nodes$value.dorc <- value.dorc
    
    
    ##is in tf-target heatmap rowid?
    
    nodes$isinHeatmap <- ifelse(nodes$name %in% DORCs,'yes','no')
    
    
    #rename dorc and tf id, if one is both dorc and tf
    #sum(duplicated(nodes$name))
    nodes$name <- paste0(nodes$name)
    
    idx <- which(duplicated(nodes$name))
    nodes$name[idx] <- paste0( nodes$name[idx],'_',nodes$group[idx]  )
    
    stopifnot( sum(duplicated(nodes$name) ) == 0 )
    
    #nodes[duplicated(nodes$name),]
    #nodes[which(nodes$name == "ZNF467"),]
    
    ###edge table
    edges <- as.data.frame(net.dat)
    
    links <- data.frame(
        source = unlist(lapply(edges$Motif, function(x) {
                        which(nodes$name == x) })),
        target = unlist(lapply(edges$DORC, function(x) {
                        which(nodes$name == x) })),
        
        corr = edges$Corr, 
        enrichment = edges$Enrichment.P
    )
    links$Value <- scales::rescale(edges$Score) * 20
    
    links$source_name <- as.character(nodes$name[links$source])
    links$target_name <- as.character(nodes$name[links$target])
    
    stopifnot(all.equal(links$source_name,edges$Motif)) #TRUE
    stopifnot(all.equal(links$target_name,edges$DORC)) #TRUE
    stopifnot(all.equal(links$corr,edges$Corr))#TRUE
    stopifnot(all.equal(links$enrichment,edges$Enrichment.P)) #TRUE
    
    
    if(save){    
       write.table(x = nodes,file = paste0('plot_network/nodes.table.',module,'-cidrna',cid.rna,'-ciddorc',cid.dorc,'.txt'),sep = '\t',quote = FALSE, row.names = FALSE, col.names = TRUE)
       write.table(x = edges,file = paste0('plot_network/edges.table.',module,'-cidrna',cid.rna,'-ciddorc',cid.dorc,'.txt'),sep = '\t',quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    }
    
    
   return(tflist)   
    
}



#module slice
#DORCs = module1$row_slice  
#TFs = module1$col_slice
#TFs <- TFs[!TFs %in% c('JUND','BACH1')]

# save=TRUE
# data.rna = exprMat.z.aggre
# cid.rna = '5'
# data.dorc= dorcMat.aggre
# cid.dorc = '5'

# module = 'module1'                          


# tflist_module1_module2 <- getNetworkDataList(figR.d = figR.d, 
#                    score.cut = 1.5,#1.2, 
#                    DORCs = c(module1$row_slice,module2$row_slice), 
#                    TFs = c(module1$col_slice,module2$col_slice),#[!module1$col_slice %in% c('JUND','BACH1')],
#                    save=TRUE,
#                    data.rna = exprMat.z.aggre, 
#                    cid.rna = '3',
#                    data.dorc= dorcMat.aggre , 
#                    cid.dorc = '1', 
#                    module = 'module1_module2'
#                   )


tflist_module1 <- getNetworkDataList(figR.d = figR.d, 
                   score.cut = 1.2, #we used 1.5 to identify tf-target heatmap, but use 1.2 to collect more tf target for exploring
                   DORCs = module1$row_slice, 
                   TFs = module1$col_slice,#[!module1$col_slice %in% c('JUND','BACH1')],
                   save=TRUE,
                   data.rna = exprMat.z.aggre, 
                   cid.rna = '2',
                   data.dorc= dorcMat.aggre , 
                   cid.dorc = '3', 
                   module = 'module1'
                  )

tflist_module2 <- getNetworkDataList(figR.d = figR.d, 
                   score.cut = 1.2, 
                   DORCs = module2$row_slice, 
                   TFs = module2$col_slice,#[!module2$col_slice %in% c('JUND','BACH1')],
                   save=TRUE,
                   data.rna = exprMat.z.aggre, 
                   cid.rna = '3',
                   data.dorc= dorcMat.aggre , 
                   cid.dorc = '1', 
                   module = 'module2'
                  )

tflist_module3 <- getNetworkDataList(figR.d = figR.d, 
                   score.cut = 1.2, 
                   DORCs = module3$row_slice, 
                   TFs = module3$col_slice,#[!module2$col_slice %in% c('JUND','BACH1')],
                   save=TRUE,
                   data.rna = exprMat.z.aggre, 
                   cid.rna = '3',
                   data.dorc= dorcMat.aggre , 
                   cid.dorc = '2', 
                   module = 'module3'
                  )


tflist_module4 <- getNetworkDataList(figR.d = figR.d, 
                   score.cut = 1.2, 
                   DORCs = module4$row_slice, 
                   TFs = module4$col_slice,#[!module2$col_slice %in% c('JUND','BACH1')],
                   save=TRUE,
                   data.rna = exprMat.z.aggre, 
                   cid.rna = '11',
                   data.dorc= dorcMat.aggre , 
                   cid.dorc = '9', 
                   module = 'module4'
                  )

tflist_module5 <- getNetworkDataList(figR.d = figR.d, 
                   score.cut = 1.0,#to keep FLT1 
                   DORCs = module5$row_slice, 
                   TFs = c(module5$col_slice,'GCM1'),#[!module2$col_slice %in% c('JUND','BACH1')], #
                   save=TRUE,
                   data.rna = exprMat.z.aggre, 
                   cid.rna = '5',
                   data.dorc= dorcMat.aggre , 
                   cid.dorc = '2', 
                   module = 'module5'
                  )


tflist_module6 <- getNetworkDataList(figR.d = figR.d, 
                   score.cut = 1.2, 
                   DORCs = module6$row_slice, 
                   TFs = c(module6$col_slice,"CREB1"),#[!module2$col_slice %in% c('JUND','BACH1')],
                   save=TRUE,
                   data.rna = exprMat.z.aggre, 
                   cid.rna = '7',
                   data.dorc= dorcMat.aggre , 
                   cid.dorc = '7', 
                   module = 'module6'
                  )

tflist_module7 <- getNetworkDataList(figR.d = figR.d, 
                   score.cut = 1.2, 
                   DORCs = module7$row_slice, 
                   TFs = module7$col_slice,#[!module2$col_slice %in% c('JUND','BACH1')],
                   save=TRUE,
                   data.rna = exprMat.z.aggre, 
                   cid.rna = '6',
                   data.dorc= dorcMat.aggre , 
                   cid.dorc = '5', 
                   module = 'module7'
                  )





saveRDS(exprMat.z.aggre,'plot_network/exprMat.z.aggre.rds')
saveRDS(dorcMat.aggre ,'plot_network/dorcMat.aggre.rds')


##############use igraph#######                           
                               
plotTFGraph_list <- function(tflist = NULL,layout = NULL, seed = NULL, downsample = NULL, geneHi = NULL,node_label = NULL, save = NULL,width= NULL, height = NULL){
    #need igraph, ggraph
    #input TF-target list (or nestted list), plot igraph by ggraph with kk layout
    
      set.seed(seed)
      if(save){
        tflist.df <- list2df(tflist)
        write.table(tflist.df, file =  paste('plot_TF_regulon_network/',paste( names(tflist),collapse = '_'  ),'-tf-target.txt',sep='') ,col.names = TRUE, row.names = FALSE, sep = '\t',quote = FALSE)  #for cytoscape input
      }
    
    
     #downsampling but keep geneHi  
     if(is.numeric(downsample) ){
         tflist_downsample <- lapply(tflist,function(x){  
                            lengthx <- length(x)
                            if(lengthx >= 200){
                              geneHi_hit <- x[x %in% geneHi]
                              x <- x[!x %in% geneHi]
                              return(c(sample(x,size = downsample*lengthx,replace = FALSE),geneHi_hit) )
                            }else{return(x)}
                          } 
                         )
      }else{ tflist_downsample <- tflist  }
    
    
      ##create graph
      g <- list2graph(tflist_downsample)
    
      #quick plot
      ##options(repr.plot.width= width, repr.plot.height = height)
      ##set.seed(10);plot(g)
      
      ##set vertex attribute
      vtx <- igraph::V(g) #vertex
      n_vtx <- length(vtx)
      id_vtx <- names(vtx)
      
      geneHi_idx <- match(geneHi, id_vtx)
      geneHi_idx <- geneHi_idx[!is.na(geneHi_idx)]
      #id_vtx[geneHi_idx]
    
      n_pathvtx <- length(tflist)
      id_pathvtx <- names(tflist) 
    
      n_genevtx <- n_vtx - n_pathvtx
      id_genevtx <- id_vtx[-(1:n_pathvtx)]
    
      path_size <- sapply(tflist, length)
      igraph::V(g)$size <- min(path_size)/2 #the default dot size
      
      igraph::V(g)$size[1:n_pathvtx] <- path_size
      
      ##also set each gene node size as foldchange
      #id_gene <- id_vtx[(n_pathvtx+1):n_vtx]
      
    
      V(g)$color <- "#B3B3B3"
      V(g)$color[1:n_pathvtx] <- "grey90"
      V(g)$color[geneHi_idx] <- 'red'
    
      #set edge attribute
      #edge <- igraph::E(g)
      E(g)$width <- 0.5
      E(g)$color <- 'darkgrey'
      #E(g)$betweenness <- edge.betweenness(g)
    
      ##quick plot 
      ##set.seed(10)
      ##options(repr.plot.width= width, repr.plot.height = height)
      ##plot(g)
    
     ##use pretty ggraph with layout kk
     n <- length(tflist)
    
     #node size
     cex_category = 3
     cex_gene = 1.5
     
     #label size
     node_label_size = NULL
     
     label_gene = 1
     cex_label_gene = 3
    
     label_category =1 
     cex_label_category = 4
    
     layout =  layout#'kk' #'circlepack'# 
     circular = FALSE
     colorEdge = FALSE
    
    node_scales <- c(rep(cex_category, n), rep(cex_gene, (length(V(g)) - n)))
    
    if (colorEdge) {
        E(g)$category <- rep(names(tflist), sapply(tflist, length))
        edge_layer <- geom_edge_link(aes_(color = ~category), alpha=.8) #change geom_edge to geom_edge_link
    } else {
        edge_layer <- geom_edge_link(alpha=.8, colour='darkgrey') #change geom_edge to geom_edge_link
    }
    
    show_legend <- c(FALSE, FALSE)
    names(show_legend) <- c("color", "size")
    
    set.seed(10)
    p <- ggraph(g, layout=layout, circular = circular)
    p <- p +
        edge_layer +
        #geom_node_point()
        geom_node_point(aes_(color=~I("#E5C494"), size=~size),
            data = p$data[1:n, ], show.legend = show_legend) +
        scale_size(range=c(3, 8) * cex_category) +
        ##ggnewscale::new_scale("size") +
    
        ##ggnewscale::new_scale_color() +
        geom_node_point(aes_(color=~color, size=~size),
            data = p$data[-(1:n), ], show.legend = show_legend) +
        ##geom_label(aes(x = x, y = y, label = name), nudge_y = 0.15, label.size = NA,fill = NA,size = 5) +
        scale_size(range=c(3, 3) * cex_gene) +
#         scale_colour_gradient2(name = "fold change", low = "blue",
#                                mid = "white", high = "red") +
        ggtitle(paste('TF-network ',paste( names(tflist),collapse = '_'  ),' layout:',layout,' downsample: ',downsample,sep='') ) +
        theme_void()+
        theme(legend.position = 'none')

    ##add node (category + gene) label with repelling
    if (node_label == "category") {
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,],
                size = label_category * cex_label_category, bg.color = "white")
        } else {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,],
                size = label_category * cex_label_category)
        }
    } else if (node_label == "gene") {
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                repel=TRUE, size = label_gene * cex_label_gene, bg.color = "white")
        } else {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
            repel=TRUE, size = label_gene * cex_label_gene)
        }
    } else if (node_label == "all") {
        if (utils::packageVersion("ggrepel") >= "0.9.0") {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                    repel=TRUE, size = label_gene * cex_label_gene, bg.color = "white") + 
                geom_node_text(aes_(label=~name), repel=TRUE,
                    size = label_category * cex_label_category, bg.color = "white", data = p$data[1:n,]) 
        } else {
            p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                    repel=TRUE, size = label_gene * cex_label_gene) + 
                geom_node_text(aes_(label=~name), data = p$data[1:n,],
                    repel=TRUE, size = label_category * cex_label_category)
        }

    }
    
    
    if(save){  
       pdf(file = paste('plot_network/',paste( names(tflist),collapse = '_'  ),'-tf-target.',layout,'.pdf',sep=''),width = width, height = height)
       print(p) 
       dev.off()
    }else{
      options(repr.plot.width= width, repr.plot.height = height)
      print(p)
    }
    return( list(g = g, p = p) )

}



list2graph <- function(inputList) {
    x <- list2df(inputList)
    g <- graph.data.frame(x, directed=FALSE)
    return(g)
}


list2df <- function(inputList) {
    # ldf <- lapply(1:length(inputList), function(i) {
    ldf <- lapply(seq_len(length(inputList)), function(i) {
        data.frame(categoryID=rep(names(inputList[i]),
                                  length(inputList[[i]])),
                   Gene=inputList[[i]])
    })

    do.call('rbind', ldf)
}


###start to plot simple igraph network

#geneHi <- c('FLT1','FLT4','PAPPA','PAPPA2','CSHL1','CSH2','LEP','FSTL4')
#geneHi <- c('FLT1','FLT4','PAPPA','PAPPA2','CSHL1','CSH2','LEP','FSTL4','ERVFRD-1','REL','ESRRG','PPARD','GLUL','FAR2','ANXA1','LIMD1','')

#geneHi <- c('LAMA3','ADCY5','SLC25A25','XDH')                               
#geneHi <- c('PAPPA', 'DTNB', 'C1ATNF6', 'OSMR', 'PDLIM5', 'XDH', 'LAMA3', 'LINC01483', 'RRAS2' ,'AR')

geneHi <- c('PAPPA','XDH')  #module1$row_slice
res.tftarget.kk.module1 <- plotTFGraph_list(tflist = tflist_module1, downsample = FALSE,layout = 'kk',seed = 123, geneHi = geneHi ,node_label = 'all', save = FALSE,width= 7.5, height = 7.5)

geneHi <- c('LAMA3','ADAMTSL1')  #module2$row_slice
res.tftarget.kk.module2 <- plotTFGraph_list(tflist = tflist_module2, downsample = FALSE,layout = 'kk',seed = 123, geneHi =  geneHi,node_label = 'all', save = FALSE,width= 7.5, height = 7.5)

geneHi <- c('LAMA3','ADCY5','SLC25A25','XDH')  
res.tftarget.kk.module3 <- plotTFGraph_list(tflist = tflist_module3, downsample = FALSE,layout = 'kk',seed = 123, geneHi = geneHi ,node_label = 'all', save = FALSE,width= 7.5, height = 7.5)

res.tftarget.kk.module4 <- plotTFGraph_list(tflist = tflist_module4, downsample = FALSE,layout = 'kk',seed = 123, geneHi = geneHi ,node_label = 'all', save = FALSE,width= 7.5, height = 7.5)

res.tftarget.kk.module5 <- plotTFGraph_list(tflist = tflist_module5, downsample = FALSE,layout = 'kk',seed = 123, geneHi = geneHi ,node_label = 'all', save = FALSE,width= 7.5, height = 7.5)

res.tftarget.kk.module6 <- plotTFGraph_list(tflist = tflist_module6, downsample = FALSE,layout = 'kk',seed = 123, geneHi = geneHi ,node_label = 'all', save = FALSE,width= 7.5, height = 7.5)


#res.tftarget.stress <- plotTFGraph_list(tflist = regulon.list, downsample = 0.5,layout = 'stress',seed = 123, geneHi = geneHi ,node_label = 'all', save = FALSE,width= 15, height = 10)


res.tftarget.nicely <- plotTFGraph_list(tflist = tflist_module1, downsample = FALSE,layout = 'nicely',seed = 123, geneHi = geneHi ,node_label = 'all', save = FALSE,width= 7.5, height = 7.5)




#res.tftarget.nodownsample <- plotTFGraph_list(tflist = tflist, downsample = FALSE,layout = 'kk',seed = 123, geneHi = geneHi ,node_label = 'all', save = FALSE,width= 7.5, height = 7.5)
#plot(layout_with_sparse_stress(res.tftarget.nodownsample$g, pivots = 100) )



#######script end??####










############get detail for given TF-DORC regulatory unit #########

##find TFs that regulate given dorc gene (often marker gene)

getRegulator <- function(figR.d = NULL, gene = NULL, corr.cutoff = 0.5, enr.cutoff = 1, save = NULL){

    figR.d.sel <- subset(figR.d,DORC==gene)
    figR.d.sel.hl1 <- subset(figR.d.sel, figR.d.sel$Corr.log10P > corr.cutoff & figR.d.sel$Enrichment.log10P > enr.cutoff) #highlight, use positive corr and enr only
    stopifnot(nrow(figR.d.sel) != 0)
    
    if(save){pdf( file = paste('tf-regulator.',gene,'.cor_',corr.cutoff,'.enr_',enr.cutoff,'.pdf',sep=''), width = 7.5, height = 7.5, useDingbats = FALSE ) }
    
    options(repr.plot.height = 7.5, repr.plot.width = 7.5)
    plot(figR.d.sel$Corr.log10P,figR.d.sel$Enrichment.log10P,pch=19,cex=0.1,col='grey', xlab='DORC accessibility vs TF expression, log10(pvalue)', ylab = 'TF enrichment in knn-DORC', main = gene)
    #abline(v=0,h=0,lty=2)
    #abline(v=0,h=0,lty=2)
    abline(v=c(-1*corr.cutoff,corr.cutoff),h=c(-1*enr.cutoff, enr.cutoff), lty=2, col = 'grey'  )
    points(figR.d.sel.hl1$Corr.log10P,figR.d.sel.hl1$Enrichment.log10P,pch=19,cex=0.1,col='red')
    text(figR.d.sel.hl1$Corr.log10P,figR.d.sel.hl1$Enrichment.log10P,labels = paste0(figR.d.sel.hl1$Motif,'-',figR.d.sel.hl1$DORC), pos = 4, cex = 0.5 )

    if(save){dev.off()}
    
    return(paste0(figR.d.sel.hl1$Motif,'-',figR.d.sel.hl1$DORC))
}



tfcandidator.FLT1 <- getRegulator(figR.d = figR.d, gene = 'FLT1', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE)
'BACH1-FLT1''CXXC5-FLT1''FOSL2-FLT1''GMEB2-FLT1''HLF-FLT1''JUNB-FLT1''JUND-FLT1''ZBTB43-FLT1'

tfcandidator.FLT4 <- getRegulator(figR.d = figR.d, gene = 'FLT4', corr.cutoff = 0.5, enr.cutoff = 1, save = TRUE)
'ATF3-FLT4''ATF6-FLT4''CREB3L1-FLT4''FOXN2-FLT4''GATA3-FLT4''GRHL1-FLT4''JUNB-FLT4''LHX2-FLT4''TWIST1-FLT4''XBP1-FLT4'

tfcandidator.PAPPA <- getRegulator(figR.d = figR.d, gene = 'PAPPA', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE)
'ATF3-PAPPA''BCL6-PAPPA''GRHL1-PAPPA''HBP1-PAPPA''MITF-PAPPA''MTF1-PAPPA''SREBF1-PAPPA''STAT3-PAPPA''STAT5A-PAPPA''STAT5B-PAPPA''TFE3-PAPPA''USF2-PAPPA''XBP1-PAPPA''ZSCAN9-PAPPA'

tfcandidator.PAPPA2 <- getRegulator(figR.d = figR.d, gene = 'PAPPA2', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE)#no hit

tfcandidator.CSHL1 <- getRegulator(figR.d = figR.d, gene = 'CSHL1', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE)

tfcandidator.STAT5A <- getRegulator(figR.d = figR.d, gene = 'STAT5A', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE)  #no STAT5A DORC

tfcandidator.STAT4 <- getRegulator(figR.d = figR.d, gene = 'STAT4', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE)
'ATF3-STAT4''GRHL1-STAT4''MITF-STAT4''MTF1-STAT4''SREBF1-STAT4''STAT3-STAT4''STAT5A-STAT4''STAT5B-STAT4''TFE3-STAT4''XBP1-STAT4'


tfcandidator.ERVFRD1 <- getRegulator(figR.d = figR.d, gene = 'ERVFRD-1', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE)
'CUX1-ERVFRD-1''FOXM1-ERVFRD-1''FOXP1-ERVFRD-1''TP53-ERVFRD-1''TP63-ERVFRD-1''ZIM2-ERVFRD-1''ZNF124-ERVFRD-1''ZNF16-ERVFRD-1''ZNF594-ERVFRD-1''ZNF680-ERVFRD-1''ZNF708-ERVFRD-1''ZNF71-ERVFRD-1''ZNF771-ERVFRD-1''ZNF823-ERVFRD-1'


tfcandidator.INHBA <- getRegulator(figR.d = figR.d, gene = 'INHBA', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE) 
'FOSL2-INHBA''HMG20B-INHBA''JUNB-INHBA''ZNF550-INHBA' #FLT1 like


tfcandidator.TGFB1 <- getRegulator(figR.d = figR.d, gene = 'TGFB1', corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE) 

##for given gene set to find comman TF?
goid_extract_genes.list <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/cisTarget/i-cisTarget/DARs_lift2hg19_fulldb_auc0.001/goid_extract_hypoxia_angiogenesis.rds')

#goid_extract_genes <- unique(unlist(goid_extract_genes.list) )
#goid_extract_genes <- c('LEP','EGLN3','LIMD1','HDAC2','ANGPTL4','ATG7','BNIP3L','NDRG1','CARD16','HK2','SLC2A1','PRKCE','NPEPPS','PAM','EGLN1','SERPINE1','FLT1','SASH1','ITGA5','ENG','GATA2','ADAM12','JUP','ALOX5')

#list(
#'response to oxygen levels' = c( 'LEP' 'EGLN3' 'LIMD1' 'HDAC2' 'ANGPTL4' 'ATG7' 'BNIP3L' 'NDRG1' 'CARD16' 'HK2', 'SLC2A1', 'PRKCE', 'NPEPPS', 'PAM', 'EGLN1' )
#positive regulation of angiogenesis :  SERPINE1 FLT1 ANGPTL4 SASH1 ITGA5 HK2 ENG GATA2 ADAM12 JUP 
#response to hypoxia :  LEP EGLN3 LIMD1 ANGPTL4 BNIP3L NDRG1 CARD16 HK2 SLC2A1 PRKCE NPEPPS PAM EGLN1 
#response to decreased oxygen levels :  LEP EGLN3 LIMD1 ANGPTL4 BNIP3L NDRG1 CARD16 HK2 SLC2A1 PRKCE NPEPPS PAM EGLN1 
#regulation of angiogenesis :  LEP SERPINE1 FLT1 ANGPTL4 SASH1 ALOX5 ITGA5 HK2 ENG GATA2 ADAM12 JUP EGLN1 

#)

genes_hypoxia <- unique(unlist(goid_extract_genes.list[grep('oxygen|hypoxia',names(goid_extract_genes.list) )]) )
#'LEP''EGLN3''LIMD1''HDAC2''ANGPTL4''ATG7''BNIP3L''NDRG1''CARD16''HK2''SLC2A1''PRKCE''NPEPPS''PAM''EGLN1'
for(i in genes_hypoxia){ cat(i,'\n') }

genes_angio <- unique(unlist(goid_extract_genes.list[grep('angiogene',names(goid_extract_genes.list) )]) )
#'SERPINE1''FLT1''ANGPTL4''SASH1''ITGA5''HK2''ENG''GATA2''ADAM12''JUP''LEP''ALOX5''EGLN1'

for(i in genes_angio){ cat(i,'\n') }


res.tfcandidator.list <- list()
#for(geneid in genes_hypoxia){
for(geneid in genes_angio){
    cat('do tfcandidator collection for gene ',geneid,'\n',sep='')
    if(geneid %in% figR.d$DORC){}else{  cat('skip '); next}
   res.tfcandidator.list[[geneid]] <- getRegulator(figR.d = figR.d, gene = geneid, corr.cutoff = 0.5, enr.cutoff = 1, save = FALSE) 
    
}



##find target dorc gene for given tf

exprMat.ave.z <- readRDS("/sda/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/chromVAR_TF_specific/full12526/TF_dev-vs-TF_expression/exprMat.ave.z.full11206.rds")

exprMat.ave.z.te <-readRDS("/sda/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/chromVAR_TF_specific/full12526/TF_dev-vs-TF_expression/exprMat.ave.z.rds")


colnames(exprMat.ave.z.te) <- paste0('c',colnames(exprMat.ave.z.te))


# score.cutoff <- 0.2

# tfid <- 'ESRRG'
# tfid <- 'MITF'
# tfid <- 'STAT5A'

# tfid <- 'STAT4'

# tfid <- 'BACH1'
# tfid <- 'FOSL2'
# tfid <- 'TWIST1'

getTFtarget <- function(tfid = NULL, score.cutoff = NULL,width = NULL, height = NULL){
    #only filtering by regulation score, no additional expression correlation of tf gene with target gene candidates
    
    ##tftarget <- as.character(subset(figR.d,Motif == tfid & Score > score.cutoff)$DORC)
    ##tftarget <- as.character(subset(figR.d,Motif == tfid & abs(Score) > score.cutoff)$DORC)
    #tftarget <- as.character(subset(figR.d,Motif == tfid & Score < score.cutoff)$DORC)
    
    tftarget.df <- subset(figR.d,Motif == tfid & abs(Score) > score.cutoff)
    
    stopifnot(nrow(tftarget.df) != 0 )
    
    res.stat <- table(tftarget.df$DORC %in% rownames(exprMat.ave.z))

    tftarget.df.pos <- subset(tftarget.df, Score > 0)
    tftarget.df.neg <- subset(tftarget.df, Score < 0)
    
    tftarget.df.pos$DORC <- as.character(tftarget.df.pos$DORC)
    tftarget.df.neg$DORC <- as.character(tftarget.df.neg$DORC)
    
    #for positive regulated target
    if( length(unique(tftarget.df.pos$DORC)) == 0 ){ cat('empty positive regulation\n')  }else{
        options(repr.plot.width = width, repr.plot.height = height)
        hp1 = Heatmap(as.matrix(exprMat.ave.z)[c(unique(tftarget.df.pos$DORC),tfid),], name = "tf-target expr", 
                         cluster_rows = TRUE, 
                         cluster_columns = TRUE, 
                         show_row_names = TRUE,
                         use_raster = FALSE,
                         ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
                         column_title = paste('Positive tf-target gene expression\n under regulation of ',tfid,' regulation score cutoff:',score.cutoff,sep=''),
                         clustering_method_row = "ward.D", ##ward.D,complete
                         clustering_distance_rows  = "pearson",
                         row_names_gp = gpar(fontsize = 8),
                         column_names_gp = gpar(fontsize = 20),
                         column_dend_height = unit(2, "cm"), 
                         row_dend_width = unit(2, "cm"),
                         #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
                         #heatmap_legend_param = list(color_bar = "continuous")
                         ##right_annotation = ha#,heatmap_width=unit(8, "cm")
            ) 

        ComplexHeatmap::draw(hp1)
    }
    
    #for negative regulated target
    if( length(unique(tftarget.df.neg$DORC)) == 0 ){ cat('empty negative regulation\n') }else{
        options(repr.plot.width = width, repr.plot.height = height)
        hp2 = Heatmap(as.matrix(exprMat.ave.z)[c(unique(tftarget.df.neg$DORC),tfid),], name = "tf-target expr", 
                         cluster_rows = TRUE, 
                         cluster_columns = TRUE, 
                         show_row_names = TRUE,
                         use_raster = FALSE,
                         ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
                         column_title = paste('Negative tf-target gene expression\n under regulation of ',tfid,' regulation score cutoff:',-1*score.cutoff,sep=''),
                         clustering_method_row = "ward.D", ##ward.D,complete
                         clustering_distance_rows  = "pearson",
                         row_names_gp = gpar(fontsize = 8),
                         column_names_gp = gpar(fontsize = 20),
                         column_dend_height = unit(2, "cm"), 
                         row_dend_width = unit(2, "cm"),
                         #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
                         #heatmap_legend_param = list(color_bar = "continuous")
                         ##right_annotation = ha#,heatmap_width=unit(8, "cm")
            ) 

        ComplexHeatmap::draw(hp2)
    }
    
    return(list(pos=tftarget.df.pos,neg=tftarget.df.neg))
}



tfid <- 'ESRRG'
score.cutoff <- 0.2
tftarget.ESRRG <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)


tfid <- 'MITF'
score.cutoff <- 0.5
tftarget.MITF <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)


tfid <- 'STAT5A'
score.cutoff <- 0.2
tftarget.STAT5A <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)

grep('PAPPA',tftarget.STAT5A[['pos']]$DORC,value=TRUE)
grep('MKLN1',tftarget.STAT5A[['pos']]$DORC,value=TRUE)
grep('GH',tftarget.STAT5A[['pos']]$DORC,value=TRUE)
grep('CSH',tftarget.STAT5A[['pos']]$DORC,value=TRUE)



tfid <- 'STAT4'
score.cutoff <- 0.5
tftarget.STAT4 <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)


tfid <- 'TWIST1'
score.cutoff <- 1
tftarget.TWIST1 <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)


tfid <- 'BACH1'
score.cutoff <- 1
tftarget.BACH1 <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)


tfid <- 'FOSL2'
score.cutoff <- 0.8
tftarget.FOSL2 <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)

tfid <- 'JUND'
score.cutoff <- 1
tftarget.JUND <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)

tfid <- 'JUNB'
score.cutoff <- 0.5
tftarget.JUNB <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)


tfid <- 'HLF'
score.cutoff <- 0.5
tftarget.HLF <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 12)

tfid <- 'MSX1'
score.cutoff <- 0.3
tftarget.MSX1 <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 15)

tfid <- 'MSX2'
score.cutoff <- 0.5
tftarget.MSX2 <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 15)


tfid <- 'TEAD4'
score.cutoff <- 0.5
tftarget.TEAD4 <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 15)

# tfid <- 'HIF1A'
# score.cutoff <- 0.5
# tftarget.HIF1A <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 15)
#no hit

tfid <- 'ARNT'
score.cutoff <- 0.2
tftarget.ARNT <- getTFtarget(tfid , score.cutoff,width = 7.5, height = 15)





####get DORC (gene) groups (dimension reduction for aggregated dorcMat with NMF??)

                               
                               
                               
                               



####get TF-DORC-target table###


dorcGenes#156  # 481 dorc genes (number of p2g >= 5)

enrichTab <- figR.d

dorcTab.split <- split(cisCorr.filt, f = cisCorr.filt$Gene ) #10359 genes with at least one peak-to-gene links

dorc.len<- sort(sapply(X = dorcTab.split, FUN = function(x){nrow(x)} ),decreasing = TRUE)
dorc.len.df <- data.frame(Gene=names(dorc.len),n=dorc.len, stringsAsFactors = FALSE)

table( dorcGenes %in% names(dorcTab.split) )

# dorc.len['PAPPA']
# dorc.len['FLT4'] #5
# dorc.len['STAT4'] #6

# idx <- grep(pattern = 'FLT4',x=names(dorc.len))
# dorc.len[idx]
# dorc.len['FLT4']

# sum(duplicated(names(dorc.len) )) #0


# dorc.len.new <- (table(cisCorr.filt$Gene))
# #dorc.len.new['FLT4']
# dorc.len.new.df <- data.frame(sort(dorc.len.new,decreasing = TRUE), stringsAsFactors = FALSE)
# colnames(dorc.len.new.df) <- c('Gene','n')
# dorc.len.new.df$Gene <- as.character(dorc.len.new.df$Gene)


numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))


all.equal(dorc.len.df,numDorcs, check.attributes = FALSE) #TRUE
#all.equal(dorc.len.df,dorc.len.new.df, check.attributes = FALSE) #TRUE

numDorcs[grep('FLT4',numDorcs$Gene ),]

numDorcs[grep('LAMA3',numDorcs$Gene ),]


table(dorcGenes %in% names(dorcTab.split ))
TRUE 
 156

TRUE 
 481

dorcTab.split[['PAPPA']]
dorcTab.split[['PAPPA2']]

dorcTab.split[['CSHL1']]
dorcTab.split[['CSH1']]

dorcTab.split[['FLT1']]
dorcTab.split[['INHBA']]
dorcTab.split[['FLT4']]


#motif_ix <- figR.d.new$motif_ix
#motif_pos <- figR.d.newnew$motif_pos

##############output motif occuring peaks bed/and motif_pos bed for each tf motif#############

dim(assay(motif_ix) )
127983 x 860

#299537 x 918

for(i in colnames(assay(motif_ix)) ){
    cat('output tf peak ',i,'\n',sep='')
    peaks_withmotif <- rownames(motif_ix)[as.matrix(assay(motif_ix[,i]))[,1]]
    peaks_withmotif.df <- as.data.frame(t(as.data.frame(strsplit(peaks_withmotif,split = ':|-') )))
    rownames(peaks_withmotif.df) <- NULL
    peaks_withmotif.df$id <- paste0(i,"_peak_with_motif","_",1:nrow(peaks_withmotif.df))
    write.table(peaks_withmotif.df, 
         #file=paste("beds_motif_match/",tfid,'.',i,'.',pos.motif[i],'.bed',sep=''),
         file=paste("beds_motif_match/",i,'.peak_with_motif.bed',sep=''),
         sep = '\t',
         quote=FALSE,
         row.names = FALSE,
         col.names = FALSE
    )
    
}




##output motif match cordinates as bed
as.data.frame(motif_pos[['STAT5A']])

writeMotif_position <- function(                       
                               tfid = tfid,
                               matches.pos = matches.pos,
                               length.r = length.r
                              ){
    ##write all subfamilies of TF-motif positions
    #pos.motif <- grep(paste("_",tfid,"_",sep=""),names(matches.pos),value = TRUE )
    stopifnot( tfid %in% names(matches.pos)  )
    #for (i in 1:length(pos.motif) ){
      df.out <- as.data.frame(matches.pos[[tfid]])
      #df.out <- transform(df.out,width= paste0("motif_",tfid,"_",pos.motif[i],"_",1:nrow(df.out)))
      #df.out <- transform(df.out,width= paste0("motif_",tfid,"_subfamily",i,"_",1:nrow(df.out)))
      df.out <- transform(df.out,width= paste0("motif_",tfid,"_",1:nrow(df.out)))
      write.table(df.out, 
                 #file=paste("beds_motif_match/",tfid,'.',i,'.',pos.motif[i],'.bed',sep=''),
                 file=paste("beds_motif_match/",tfid,'.motif_pos.bed',sep=''),
                 sep = '\t',
                 quote=FALSE,
                 row.names = FALSE,
                 col.names = FALSE
                )
    #}
    ##write all extended subfamilies of TF-motif positions
    #for (i in 1:length(pos.motif) ){
      pos.motif.gr <- matches.pos[[tfid]]
      pos.motif.region <- resize(pos.motif.gr,fix = 'center',width = length.r)
      df.out <- as.data.frame(pos.motif.region)
      df.out <- transform(df.out,width= paste0("motif_extend",length.r,"_",tfid,"_",1:nrow(df.out)))
      write.table(df.out, 
                 file=paste("beds_motif_match/",tfid,'.extend', length.r,'.motif_pos.bed',sep=''),
                 sep = '\t',
                 quote=FALSE,
                 row.names = FALSE,
                 col.names = FALSE
                )
    #}
    return('done')
}



length.r = 450

#dir.create('beds_motif_match')


stopifnot(all.equal(colnames(motif_ix),names(motif_pos) ))
stopifnot( sum(duplicated( names(motif_pos) ) ) == 0 )

motifs.fullist <- names(motif_pos)

for (i in motifs.fullist  ){
  
  writeMotif_position(                       
                     tfid = i,
                     matches.pos = motif_pos, #matches.pos,
                     length.r = length.r
                    )
}




###output dorc link bed file and cisCorr.filt table full peak_to_gene junction bed#########

dorcGenes
dorcTab.split

table( dorcGenes %in% names(dorcTab.split) )
#481 TRUE in total 10359 genes

##get TSS of gene from 

tss <- read.table('/home/mjwang/dataex/cellranger_genomes/refdata-cellranger-atac-GRCh38-1.1.0/regions/tss.bed',header = FALSE, sep = '\t', stringsAsFactors = FALSE) #slightly different with cellranger-atac.1.2

colnames(tss) <- c('chr','start','end','geneid','score','strand','type')


##output dorc link bed file

write.table(x = 'track name=junctions description="DORC (481 genes) peak to gene links"',
            file = 'dorc_peak_to_gene_link/dorc.dorcGenes.junction.bed',
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t',
            quote=FALSE
           )


for(i in dorcGenes){
    #cat('dorc gene peak to gene links for ',i,'\n',sep='')
    dorc.df <- dorcTab.split[[i]]
    stopifnot(nrow(dorc.df) > 0)
    
    ##define tss region
    tss.gene <- subset(tss,geneid == i)
    #stopifnot(nrow(tss.gene) > 0)
    
    if(nrow(tss.gene) == 0){
      cat('can not found tss for gene ',i,'\n',sep='')
      next
    }
    
    if(nrow(tss.gene) == 1){ 
      dorc.df.ext <- data.frame(dorc.df, tss=paste0(tss.gene$chr,':',tss.gene$start,'-',tss.gene$end) ,stringsAsFactors = FALSE)
    
    }else if(nrow(tss.gene) > 1){
        #choose the most upstream one?
        strand <- unique(tss.gene$strand)
        stopifnot(length(strand) == 1)
        if(strand == "+"){ 
          tss.gene <- tss.gene [which.min(tss.gene$start),]
            
        }else if(strand == '-'){
            tss.gene <- tss.gene [which.max(tss.gene$end),]
            
        }
        stopifnot(nrow(tss.gene) == 1)
        dorc.df.ext <- data.frame(dorc.df, tss=paste0(tss.gene$chr,':',tss.gene$start,'-',tss.gene$end) ,stringsAsFactors = FALSE)
        
        
    }
    
    
    ##igv junction format
    #track name=junctions description="Cicero_Granja connections.filter3.0 links"
    #chr1    634199  778977  link1   0.305989577883302       +
    #chr1    634199  845876  link2   0.435166666059608       +
    #chr1    634199  877093  link3   0.370885687189481       +
    dorc.df.ext.expand <- dorc.df.ext %>% 
      tidyr::separate(PeakRanges,c('chr_peak','start_peak','end_peak')) %>% 
      tidyr::separate(tss,c('chr_tss','start_tss','end_tss'))
    
    stopifnot( all(dorc.df.ext.expand$chr_peak == dorc.df.ext.expand$chr_tss) )
    dorc.table <- data.frame(chr=dorc.df.ext.expand$chr_peak,
               start =  ceiling(0.5*(as.numeric(dorc.df.ext.expand$start_peak)
                                     +as.numeric(dorc.df.ext.expand$end_peak)
                                    )
                               ),
               end =  ceiling(0.5*(as.numeric(dorc.df.ext.expand$start_tss)
                                     +as.numeric(dorc.df.ext.expand$end_tss)
                                    )
                               ),
               id = paste('dorc_',i,'_link',1:nrow(dorc.df.ext.expand),sep=''),
               #score = -1*log10(dorc.df.ext.expand$rObs ),
               score = dorc.df.ext.expand$rObs,
               strand = '+'
            )
    write.table(x = dorc.table,
            file = 'dorc_peak_to_gene_link/dorc.dorcGenes.junction.bed',
            row.names = FALSE,
            col.names = FALSE,
            sep='\t',
            quote=FALSE,
            append = TRUE
           )
    
}


can not found tss for gene LINC01136
can not found tss for gene LINC01095
can not found tss for gene PROSER2-AS1
can not found tss for gene ARHGAP26-AS1
can not found tss for gene FER1L6-AS2
can not found tss for gene ISM1
can not found tss for gene LINC00474
can not found tss for gene LINC00882
can not found tss for gene ABHD17C
can not found tss for gene CALML3-AS1
can not found tss for gene CFLAR-AS1
can not found tss for gene CSE1L-AS1
can not found tss for gene GRK3
can not found tss for gene GSN-AS1
can not found tss for gene IQCJ-SCHIP1-AS1
can not found tss for gene LINC01119
can not found tss for gene LINC01267
can not found tss for gene LINC01270
can not found tss for gene LINC01307
can not found tss for gene LINC01411
can not found tss for gene LINC01509
can not found tss for gene LZTS1-AS1
can not found tss for gene MAFA-AS1
can not found tss for gene N4BP3
can not found tss for gene OGFRP1




#########summary tf-peak-gene table for cytoscape, along with node data####



getTable <- function(cid = NULL, tf_select = NULL, figR.d = NULL, score.cutoff.list = NULL, dorcTab.split = NULL, motif_ix = NULL, save = NULL){    
    ####get TF-dorc first from enrich table with filtering, iterative for all dorc genes (number of p2g >= 5)
    ####also add edge attributions (tf expr - dorc score correlation;  peak to gene link correlation)###
    ##A  B  edge_attr1 edge_attr2
    ##B  C  edge_attr1 edge_attr2
    
    ##edge type
    ##TF  docr (tf-to-dorc(peaks))
    ##peak  gene (peak-to-gene)

    ##get node table but fill node data later
    
    
#     edge.table <- data.frame()
    
#     corr.cutoff = 0.5
#     enr.cutoff = 1
#     #for(gene in dorcGenes){
#     for(gene in c('FLT1','FLT4','PAPPA','PAPPA2','INHBA') ){
#       figR.d.sel <- subset(figR.d,DORC==gene)
#       if(nrow(figR.d.sel) == 0){cat('empty dorc of gene ',gene,'\n');next}
#       figR.d.sel.hl <- subset(figR.d.sel, figR.d.sel$Corr.log10P > corr.cutoff & figR.d.sel$Enrichment.log10P > enr.cutoff) #highlight
#       if(nrow(figR.d.sel.hl) == 0){cat('empty dorc sel of gene', gene,'\n');next}
        
#       tf_dorc.df <- data.frame(figR.d.sel.hl[,c('DORC','Motif')],type = 'TF-peak', stringsAsFactors = FALSE)
#       colnames(tf_dorc.df) <- c('V1','V2','V3')
#       #get dorc peak-to-gene links for this dorc gene (target gene)
#       if(gene %in% names(dorcTab)){
#           p2g <- data.frame(dorcTab[[gene]][,c('Gene','PeakRanges')],type='peak-gene',stringsAsFactors = FALSE)
#           colnames(p2g) <- c('V1','V2','V3')
#           edge.table <- rbind(edge.table,tf_dorc.df,p2g)
#       }else{stop('no dorc peaks for gene ',gene)}
    
#     }    
    
    edge.table <- data.frame()
    node.table <- data.frame()
    for(tfid in tf_select[[cid]] ){
    #for(tfid in c('STAT5A','STAT4','MITF') ){
        
        cat('##collect tf-peak-gene table for tf: ',tfid,'\n')
        
        node.table <- rbind.data.frame(node.table,data.frame(node=tfid,type='TF',stringsAsFactors = FALSE))
        
        score.cutoff <- score.cutoff.list[[cid]][[tfid]]
        
        tftarget.df <- subset(figR.d,Motif == tfid & abs(Score) > score.cutoff)
        stopifnot(nrow(tftarget.df) != 0 )

        #res.stat <- table(tftarget.df$DORC %in% rownames(exprMat.ave.z))
        #stopifnot(res.stat['TRUE'] == nrow(tftarget.df) )
        
        tftarget.df.pos <- subset(tftarget.df, Score > 0)
        tftarget.df.neg <- subset(tftarget.df, Score < 0)

        stopifnot(nrow(tftarget.df.pos) != 0 )
        
        tftarget.df.pos$DORC <- as.character(tftarget.df.pos$DORC)
        tftarget.df.neg$DORC <- as.character(tftarget.df.neg$DORC)

        stopifnot( sum(duplicated(tftarget.df.pos$DORC) ) == 0 )
        
        for(gene in unique(tftarget.df.pos$DORC)) {
            cat ('  >>collect peaks for dorc gene: ',gene,'\n')
           if(gene %in% names(dorcTab.split)){
              p2g <- data.frame(dorcTab.split[[gene]][,c('PeakRanges','Gene','rObs')],type='peak-gene',stringsAsFactors = FALSE)
               
               ##filter and keep peaks with this tf-motif
               
               peaks_withmotif <- rownames(motif_ix)[as.matrix(assay(motif_ix[,tfid]))[,1]]
               #length(peaks_withmotif)
               #9249

               #table(subset(cisCorr.filt,Gene == 'ADCY7')$PeakRanges %in% peaks_withmotif)
               #FALSE  TRUE 
               #    8     1

               #peak_withmotif_idx <- which(subset(cisCorr.filt,Gene == 'ADCY7')$PeakRanges %in% peaks_withmotif)
               peak_withmotif_idx <- which(p2g$PeakRanges %in% peaks_withmotif)

               if(length(peak_withmotif_idx) == 0){cat('  no tf motif occuring in cisCorr.filt table\n');next}
               
               #subset(cisCorr.filt,Gene == 'ADCY7')[peak_withmotif_idx,]
               p2g.filter <- p2g[peak_withmotif_idx,]

               p2g.filter <- p2g.filter[!duplicated(p2g.filter$PeakRanges),]
               
               colnames(p2g.filter) <- c('V1','V2','V3','V4')
               
               node.table <- rbind.data.frame(node.table,data.frame(node=p2g.filter$V1,type='region',stringsAsFactors = FALSE))
               node.table <- rbind.data.frame(node.table,data.frame(node=p2g.filter$V2,type='gene',stringsAsFactors = FALSE))
               
               
               tftarget.df.pos.sel <- subset(tftarget.df.pos,DORC == gene)
               stopifnot(nrow(tftarget.df.pos.sel) == 1)
               
               tf2peak <- data.frame( TF=tfid, peak =  p2g.filter$V1, Corr = tftarget.df.pos.sel$Corr ,type = 'tf-peak', stringsAsFactors = FALSE )
               colnames(tf2peak) <- c('V1','V2','V3','V4')
               

               edge.table <- rbind(edge.table,tf2peak,p2g.filter)
          }else{stop('no dorc peaks for dorc gene ',gene)}

        }
        
    }
    
    edge.table <- edge.table[!duplicated(edge.table[,c(1,2)]),]
    colnames(edge.table) <- c('source','target','weight','type')
    
    node.table <- node.table[!duplicated(node.table),]
    
    
    if(save){
        write.table(edge.table,paste('dorc_peak_to_gene_link/',cid,'.edge.table.txt',sep=''),sep = '\t',row.names=FALSE, col.names = TRUE,quote=FALSE) 
        #edge table
        write.table(node.table,paste('dorc_peak_to_gene_link/',cid,'.node.table.txt',sep='') ,sep = '\t',row.names=FALSE, col.names = TRUE,quote=FALSE) 
        #node table

    }
    
    
    return(list(edge.table=edge.table,node.table=node.table))
    
    
}



tf_select <- list( #select from tf-dorc heatmap (tf enrichment in dorc peaks of selected dorc genes)
    'c2' = c('STAT5A','STAT4','MITF'),
    'c1' = c('FOSL2','MYCN','POU2F3')
)

score.cutoff.list <- list( #regulatory score in figR.d
    'c2' = c('STAT5A' = 0.5,'STAT4' = 0.5,'MITF' = 0.5),
    'c1' = c('FOSL2' = 0.5,'MYCN' = 0.5,'POU2F3' = 0.5)

)


res.network.list <- list()

res.network.list[['c2']] <- getTable(cid = 'c2',tf_select =tf_select, figR.d = figR.d, score.cutoff.list = score.cutoff.list, dorcTab.split = dorcTab.split, motif_ix = motif_ix, save = TRUE)


res.network.list[['c1']] <- getTable(cid = 'c1',tf_select =tf_select, figR.d = figR.d, score.cutoff.list = score.cutoff.list, dorcTab.split = dorcTab.split, motif_ix = motif_ix, save = TRUE)



#all.equal(res.network.list.bk[['c2']]$edge.table[,c(1,2)], res.network.list[['c2']]$edge.table[,c(1,2)])
#TRUE

#test1 <- res.network.list.bk[['c2']]$node.table
#test2 <- res.network.list[['c2']]$node.table

#all.equal(test1[order(test1$node),1],test2[order(test2$node),1],check.attributes = FALSE) #TRUE
#all.equal(test1[order(test1$node),2],test2[order(test2$node),2],check.attributes = FALSE)




#####






#############add data  for node: tf (expression <and FC ??>), peak (accessibility) , gene (expression and percentage)################



##readin data
exprMat.perc <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/chromVAR_TF_specific/full12526/TF_dev-vs-TF_expression/exprMat.perc.rds')
#243070 x 5

exprMat.perc$id <- paste0('c',exprMat.perc$id)
table(exprMat.perc$id)
  c1   c10    c2    c3    c4    c5    c6    c7    c8    c9 
24307 24307 24307 24307 24307 24307 24307 24307 24307 24307

exprMat.perc$features.plot <- as.character(exprMat.perc$features.plot)



subset(exprMat.perc, features.plot == 'PAPPA')
exprMat.ave.z.te['PAPPA',]

cor(subset(exprMat.perc, features.plot == 'PAPPA')$avg.exp.scaled, unlist(exprMat.ave.z.te['PAPPA',]))
#0.982409310602705

cor(subset(exprMat.perc, features.plot == 'PAPPA')$avg.exp, unlist(exprMat.ave.z.te['PAPPA',]))
#0.982409310602705

cor(subset(exprMat.perc, features.plot == 'PAPPA')$pct.exp, unlist(exprMat.ave.z.te['PAPPA',]))
#0.967340836587857

sum(duplicated(rownames(exprMat.ave.z.te))) #0
exprMat.ave.z.te$geneid <- rownames(exprMat.ave.z.te) #to avoid rowid duplication when extract data



peakMat.aggre <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/cicero_Granja/peakMat.aggre.rds')

colnames(peakMat.aggre) <- paste( 'c', colnames(peakMat.aggre),sep='')

sum(duplicated(rownames(peakMat.aggre))) #0

peakMat.aggre$peakid <- rownames(peakMat.aggre) #to avoid rowid duplication 



#add data columns (with NA if not available) for nodes and edges table  of res.network.list

#cid <- 'c2'
#cid <- 'c1'



#colnames(tfmotif_nes_align)
#'c6''c3''c9''c5''c1''c2''c4'

#colnames(tfmotif_expr_align)
#'c9''c5''c10''c6''c7''c3''c2'

map_cluster <- list( #atac cid to rna cid
       'c6' = 'c9',
       'c3' = 'c5',
       'c9' = 'c10',
       'c5' = 'c6',
       'c1' = 'c7',
       'c2' = 'c3',
       'c4' = 'c2' 
)


#cid_map <- map_cluster[[cid]]


addNodeData <- function(res.network.list = NULL, map_cluster= NULL, exprMat.perc = NULL, exprMat.ave.z.te = NULL,peakMat.aggre =NULL, save = NULL){

    for(cid in names(res.network.list)){
    #for(cid in cid){
        stopifnot(!is.null(map_cluster[[cid]]))
        cid_map <- map_cluster[[cid]]
        cat ('add data for node and edge of atac cluster ', cid, ' with matched rna cluster ',cid_map,'\n')
        #tf_peak_gene.df.sel <- subset(tf_peak_gene.df, dar == cid)

        network.node <- res.network.list[[cid]][['node.table']]

        network.edge <- res.network.list[[cid]][['edge.table']]
        res.network.data.list[[cid]][['edge.table.data']] <-  network.edge #already has data


        ##add data for node table
    #     node.df <- rbind.data.frame(
    #         unique(data.frame(name = tf_peak_gene.df.sel$tf,type = 'TF',stringsAsFactors = FALSE) ),
    #         unique(data.frame(name = tf_peak_gene.df.sel$peak,type = 'region',stringsAsFactors = FALSE) ),
    #         unique(data.frame(name = tf_peak_gene.df.sel$target,type = 'gene',stringsAsFactors = FALSE) ),
    #         #unique(data.frame(name = tf_peak_gene.df.sel$target_p2g_combined,type = 'gene',stringsAsFactors = FALSE) ),
    #         stringsAsFactors = FALSE
    #     )

        ##1 add tf gene expression, and percentage for node table

        ##sum(duplicated(network.node$node)) != 0 #target gene or peak may duplicated , but is ok
        #network.node[which(duplicated(network.node$node)),]
        #network.node[grep('STAT4',network.node$node),]

        exprMat.perc.sel <- subset(exprMat.perc, id == cid_map)
        stopifnot( sum(duplicated(exprMat.perc.sel$features.plot)) == 0 )
        rownames(exprMat.perc.sel) <- exprMat.perc.sel$features.plot

        expr_perc <- exprMat.perc.sel[network.node$node,c('features.plot','pct.exp')]

        #expr_value <-  exprMat.ave.z[node.df$name,cid,drop = FALSE]
        #expr_value$geneid <- rownames(expr_value) #will add .1 if duplicated
        expr_value <-  exprMat.ave.z.te[network.node$node,c('geneid',cid_map)]
        #expr_value$geneid <- rownames(expr_value)
        colnames(expr_value) <- c('geneid','exp.value')
        #expr_value <- expr_value[,c(2,1)]

        network.node.data <- cbind.data.frame(network.node,expr_perc,expr_value)

          stopifnot(all.equal(network.node.data[!is.na(network.node.data$pct.exp),'node'],network.node.data[!is.na(network.node.data$pct.exp),'features.plot']) ) #TRUE
        stopifnot(all.equal(network.node.data[!is.na(network.node.data$value),'node'],network.node.data[!is.na(network.node.data$value),'geneid']) )#TRUE

        #network.node.data <- network.node.data[,c('node','class','features.plot','pct.exp','exp.value')]

        ##2 add peak accessibility

        #acc_value <- peakMat.aggre[node.df$name,cid,drop = FALSE]
        #acc_value$peakid <- rownames(acc_value)
        acc_value <- peakMat.aggre[network.node$node,c('peakid',cid)]
        colnames(acc_value) <- c('peakid','acc_value')
        #acc_value <- acc_value[,c(2,1)]

        network.node.data <- cbind.data.frame(network.node.data,acc_value)

        stopifnot(all.equal(network.node.data[!is.na(network.node.data$acc_value),'node'],network.node.data[!is.na(network.node.data$acc_value),'peakid']) )

        network.node.data <- network.node.data[,c('node','type','features.plot','pct.exp','exp.value','acc_value')]
        colnames(network.node.data) <- c('name','type',	'label','pct.exp','expr.value','acc_value')

        network.node.data$label <- ifelse(is.na(network.node.data$label),'',network.node.data$label)

    #     #NA to 0 for value column, not necessary for cytoscape!
    #     network.node.data$pct.exp <-  ifelse(is.na(network.node.data$pct.exp), 0  , network.node.data$pct.exp)
    #     network.node.data$expr.value <-  ifelse(is.na(network.node.data$pct.exp), 0  , network.node.data$pct.exp)
    #     network.node.data$acc_value <-  ifelse(is.na(network.node.data$pct.exp), 0  , network.node.data$pct.exp)

        res.network.data.list[[cid]][['node.table.data']] <-  network.node.data

    #     ##add data for edge table

    #     ##1, add tf-dorc correlation 

    #     network.edge

    #     figR.d
    #     cisCorr.filt
        

    }


    if(save){
        saveRDS(res.network.data.list,'dorc_peak_to_gene_link/res.network.data.list.rds')


        #write.table(res.network.data.list[['c2']][['node.table.data']],file='dorc_peak_to_gene_link/c2.node.table.data.txt',sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)

        #write.table(res.network.data.list[['c2']][['edge.table.data']],file='dorc_peak_to_gene_link/c2.edge.table.data.txt',sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)


        write.table(res.network.data.list[[cid]][['node.table.data']],file=paste('dorc_peak_to_gene_link/',cid,'.node.table.data.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)

        write.table(res.network.data.list[[cid]][['edge.table.data']],file=paste('dorc_peak_to_gene_link/',cid,'.edge.table.data.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)

    }
    
    return(res.network.data.list)


}



res.network.data.list <- list() 

#res.network.data.list.bk <- res.network.data.list



res.network.data.list <- addNodeData(res.network.list, map_cluster, exprMat.perc, exprMat.ave.z.te ,peakMat.aggre, save = TRUE)








