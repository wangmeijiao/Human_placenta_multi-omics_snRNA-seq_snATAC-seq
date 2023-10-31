

library(ArchR)


#library(SnapATAC)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)



corCutOff = 0.75
cellsToUse = NULL
k = 100
knnIteration = 500
overlapCutoff = 0.8 
maxDist = 250000
scaleTo = 10^4
log2Norm = TRUE
#predictionCutoff = 0.4
addEmpiricalPval = TRUE#FALSE
seed = 1
threads = 1 #max(floor(getArchRThreads() / 2), 1)
#verbose = TRUE




files <- list(
  "GeneScoreMatrix" = 'data/geneScoreMatrix.rds',
  "MotifMatrix" = 'data/MotifMatrix.full.rds', #zscored TF deviation score
  "GeneIntegrateMatrix" = 'data/geneIntegrateMatrix.rds',
  'PeakMatrix' = 'data/PeakMatrix.rds'
    
)


color_peak <- 'solarExtra'
color_tfdev <- 'solarExtra'                     
color_ga <- 'horizonExtra'
color_imputeRNA <- 'blueYellow'

#color_rna <- 'greenBlue'
color_rna <- 'whiteBlue'
#color_rna <- 'whitePurple'
#color_rna <- 'greyMagma'



map_cellcolor_atac <- c( #sip picked color (use this)
    '5' = '#4F1C47',
    '3' = '#8D1541',
    '4' = '#C31240',
    '8' = '#F95944',
    '6' = '#FB8D3C',
    '2' = '#FCC140',
    #'#ECDB68',
    #'#AFD16A',
    #'#60C589',
    '9' = 'darkgreen',#'#2BB9AD',
    '7' = '#4776B2',#'#337C99',
    '1' = '#3D3F69'
)



barplot(1:length(map_cellcolor_atac),col = map_cellcolor_atac,main ='map_cellcolor_atac' ,cex.main=2,names.arg = names(map_cellcolor_atac))

map_cellcolor <- map_cellcolor_atac

# reds <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf')
# blues <- c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')

# map_cellcolor <- list(
#     '1'=reds[3],#'STB-2',
#     '2'=reds[2],#'STB-3', 
#     '3'=blues[4],#'CTB-2',
#     '4'=reds[1],#'STB-4',
#     '5'=reds[6],#'naive STB',
#     '6'=blues[3],#'CTB-1',
#     '7'=reds[4],#'STB-1',
#     '8'=reds[5],#'STB-5',
#     '9'='darkgreen'#'Fusion Component'
    
# #     '1'=reds[5], #STB1
# #     '2'=reds[6],  #STB-naive
# #     '3'=reds[3],  #STB3
# #     '4'=reds[2], #STB4
# #     '5'=reds[4], #STB2
# #     '6'=reds[1], #Syncytial knot
# #     '7'='#694d9f', #STB-new
# #     '8'='darkgreen' #CTB
#  )


set.seed(seed)



#################################


ArchRProj = readRDS('/home/mjwang/pwdex/placenta_10X_combine/02.archR_harmony/PLA-8w-ATAC/proj5.rds')

# reducedDims = 'Harmony'#"IterativeLSI",
# useMatrix = "GeneIntegrateMatrix"
# dimsToUse = 1:30
# scaleDims = NULL
# corCutOff = 0.75
# cellsToUse = NULL
# k = 100
# knnIteration = 500
# overlapCutoff = 0.8 
# maxDist = 250000
# scaleTo = 10^4
# log2Norm = TRUE
# predictionCutoff = 0.4
# addEmpiricalPval = FALSE
# seed = 1
# threads = max(floor(getArchRThreads() / 2), 1)
# verbose = TRUE
# #logFile = createLogFile("addPeak2GeneLinks")


##AvailableMatrices <- getAvailableMatrices(ArchRProj)
#'GeneIntegrationMatrix''GeneScoreMatrix''MotifMatrix''PeakMatrix''TileMatrix'

#ArrowFiles <- getArrowFiles(ArchRProj)

#tstart <- Sys.time()






# dfAll <- ArchR:::.safelapply(seq_along(ArrowFiles), function(x){
#     cNx <- paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], paste0(useMatrix, "/Info/CellNames")))
#     pSx <- tryCatch({
#       h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
#     }, error = function(e){
#       if(getArchRVerbose()) message("No predictionScore found. Continuing without predictionScore!")
#       rep(9999999, length(cNx))
#     })
#     DataFrame(
#       cellNames = cNx,
#       predictionScore = pSx
#     )
# }, threads = threads) %>% Reduce("rbind", .)


# keep <- sum(dfAll[,2] >= predictionCutoff) / nrow(dfAll)
# dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),] #8925


#Get Peak Set
##peakSet <- getPeakSet(ArchRProj)
#.logThis(peakSet, "peakSet", logFile = logFile)
#peakSet <- readRDS('data/peaks.gr.rds')



# # getFeatureDF <- function(featureFile = NULL){

# #     ##featureDF <- readRDS('../atac_activity_old/gene.coords.rds') #gene coordinate as the same in gene activity score calculation, gene coordinates from Ensembl v86, total 63970
# #     #featureDF <- readRDS('../atac_activity_old/gene.coords.filter.rds') #filtered 36795

# #     ##featureDF <- readRDS('../atac_activity_old/gene.coords.cellrangerRNA.Env93.rds') #33538 total, cellranger RNA v93 gene gtf, with filter for protein coding, antisense, lncRNA gene etc

# #     #seqnames(featureDF)@values
# #     #seqnames(featureDF)@lengths

# #     featureDF <- readRDS(featureFile)
# #     idx <- vector()
# #     for(i in seqnames(featureDF)@lengths){idx <- c(idx,1:i)}

# #     featureDF <- DataFrame( seqnames=seqnames(featureDF),
# #                            start=start(featureDF),
# #                            end=end(featureDF),
# #                            strand=strand(featureDF),
# #                            #name = mcols(featureDF)$gene_name, 
# #                            idx = idx  )

# #     stopifnot(length(idx) == nrow(featureDF) )#TRUE
# #     cat('read in feature table: ',nrow(featureDF),'\n')
# #     return(featureDF)
# # }


# # peakSet <- getFeatureDF(featureFile ='data/peaks.gr.rds')




####ChIPseeker annotate peaks (not based on TSS, but based on genes) (use this?)

peakSet <- readRDS('data/peaks.gr.rds') #274189 #299537
peaks.gr.anno <- ChIPseeker::annotatePeak(
                                            peakSet,
                                            TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                                            annoDb='org.Hs.eg.db'
                                         )
peaks.gr.anno.df <- as.data.frame(peaks.gr.anno) #274189 #299537
#peaks.gr.anno.df[peaks.gr.anno.df$SYMBOL == "ERVFRD-1",]


##peaks.gr.anno.df.1 <- readRDS('data/peaks.gr.anno.rds')  #filter NA
##all.equal(peaks.gr.anno.df,peaks.gr.anno.df.1)


#Get the sequences and compute the GC content for each peak
#https://www.biostars.org/p/478444/
freqs <- Biostrings::alphabetFrequency(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, peakSet))
gc <- (freqs[,'C'] + freqs[,'G'])/rowSums(freqs) #274189 #299537


peaks.gr.anno.df$GC <- gc



idx <- vector()
for(i in seqnames(peakSet)@lengths){idx <- c(idx,1:i)}

all.equal(as.character(seqnames(peakSet) ), as.character(peaks.gr.anno.df$seqnames),check.attributes = FALSE )#TRUE
all.equal(start(peakSet), peaks.gr.anno.df$start,check.attributes = FALSE )
#TRUE

peaks.gr.anno.df$idx <- idx

##add to peak.gr mcol

mcols(peakSet) <- peaks.gr.anno.df[,c('distanceToTSS','SYMBOL','annotation','GC','idx')]
#distToGeneStart nearestGene    peakType distToTSS nearestTSS        GC       idx         N


saveRDS(peakSet,'gene_peak_links/peakSet.rds')

#Gene Info
#geneSet <- ArchR:::.getFeatureDF(ArrowFiles, useMatrix, threads = threads) #if - strand, swap start for end
# DataFrame with 17838 rows and 6 columns
#       seqnames     start       end  strand      name     idx



geneSet <- readRDS('correlate_trajectory/GeneIntegrateMatrix/featureDF.genes.trajectory1.1.rds')
#15294
#17185

all.equal(readRDS('correlate_trajectory/GeneIntegrateMatrix/featureDF.genes.trajectory1.1.rds'),readRDS('correlate_trajectory/GeneIntegrateMatrix/featureDF.genes.trajectory1.2.rds')) #TRUE

#reindex idx
idx <- vector()
for(i in geneSet$seqnames@lengths ){idx <- c(idx,1:i)}

geneSet$idx <- idx



#substitute +/- to 1/2
#swith end as start if 2

strand <- geneSet$strand
strand_number <- ifelse(test = strand == '+', 1,2)

flag <- ifelse(test = strand == '-', TRUE,FALSE)
table(flag)
FALSE  TRUE 
 7758  7536

#FALSE  TRUE 
# 8695  8490

#if (flag ){ temp <- geneSet$start; geneSet$start <- geneSet$end;  geneSet$end <- temp }

temp <- geneSet[flag,'start']
geneSet[flag,'start'] <- geneSet[flag,'end']
geneSet[flag,'end'] <- temp

geneSet$strand <- Rle(strand_number) #not rle in base

#saveRDS(geneSet,'gene_peak_links/geneSet.rds')

geneSet <- readRDS('result_good/gene_peak_links/geneSet.rds')


geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
#.logThis(geneStart, "geneStart", logFile = logFile)

#saveRDS(geneStart,'gene_peak_links/geneStart.rds')
##geneStart <- readRDS('result_good/gene_peak_links/geneStart.rds')

all.equal(geneStart , readRDS('result_good/gene_peak_links/geneStart.rds')) #TRUE


##########1 get aggregate grouplist (with dimension reduction matrix ) with ArchR method################


##########a) ArchR method, get aggreGroup from snapatac smat #########
#Get Reduced Dims


reducedDims = "Harmony"
dimsToUse = 1:30


scaleDims = NULL


rD_archr <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
# if(!is.null(cellsToUse)){
# rD <- rD[cellsToUse, ,drop=FALSE]
# }
#10121 x 30, reducedDims='Harmony'


rD_smat <- readRDS('../../../02.snapATAC_harmony/cicero_Granja/dimred.use_dr_smat.rds')
#22785 x 18 # x.sp@smat@dmat, from snapatac2 X_spectral_harmony_cstb.txt

#11093 x 16
#x.sp@smat@dmat

table(rownames(rD_archr) %in% rownames(rD_smat) )
# TRUE 
# 10121


#rD <- rD_archr
rD <- rD_smat


##rename rowid
cellid <- rownames(rD)
sum(duplicated(cellid)) #0

cellid <- sapply(stringr::str_split(cellid,'early|:|-',n=4), function(x){paste0( 'placenta.atac_',x[3],'-',x[2] ) } )

rownames(rD) <- cellid

#Subsample

idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration) #random choose 500 start point

#KNN Matrix
#.logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
knnObj <- ArchR:::.computeKNN(data = rD, query = rD[idx,], k = k) #a matrix, similar to cicero aggregroup
##knnIteration x k: 500 x 100


#Determin Overlap
#.logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k)) #vector
table(keepKnn)

  0 
500 

-1   0 
  1 499 

#  -1   0 
#   2 498

#Keep Above Cutoff
knnObj <- knnObj[keepKnn==0,]
#.logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)


#Convert To Names List
knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
rownames(rD)[knnObj[x, ]]
}) %>% SimpleList

length(unique(unlist(knnObj) ) ) #cover 8646 with rD_archr (lsi then harmony)
 #17553 #cover 10407 with rD_smat (snapatac)


lens_aggre <- unlist(lapply(knnObj,FUN = length))

table(lens_aggre)
100 
500 

# 100 
# 499 

saveRDS(knnObj,paste('gene_peak_links/','knnObj.rds',sep='') ) #reproducible if seed =1

########



#####b) get aggregate group by reusing cicero aggrGroup, as pseudotime groupList (use this?)####

aggrGroup <- readRDS('data/save-cicero-KNN-Groupings-cds.rds') #no pseudotime ordering, based on umap knn
#aggrGroup <- readRDS('matching_group_te/save-cicero-KNN-Groupings-cds.rds') #no pseudotime ordering
#4133 × 50

#2982 x 50

sum(duplicated(rownames(aggrGroup)) ) #0, no reduntance

####rename cellid
aggrGroup_new <- data.frame(row.names = rownames(aggrGroup) , stringsAsFactors =  FALSE )
for(i in 1:nrow(aggrGroup) ){
    cellid <- as.character(unlist(aggrGroup[i,] ) )
    stopifnot(sum(duplicated(cellid )) == 0 )
    
    cellid <- sapply(stringr::str_split(cellid,'early|:|-',n=4), function(x){paste0( 'placenta.atac_',x[3],'-',x[2] ) } )
    
#     idx1 <- grep(pattern = 'donor1',x = cellid )
#     idx2 <- grep(pattern = 'donor2',x = cellid )
#     stopifnot(length(c(idx1,idx2)) == length(cellid) )
    
#     cellid[idx1] <- gsub(pattern = '_donor1#',replacement = '.atac_',x = cellid[idx1])

#     cellid[idx2] <- gsub(pattern = '-1$',replacement = '-2',x = cellid[idx2])
#     cellid[idx2] <- gsub(pattern = '_donor2#',replacement = '.atac_',x = cellid[idx2])

    aggrGroup_new <- rbind.data.frame (aggrGroup_new,cellid, stringsAsFactors = FALSE )

    
}

colnames(aggrGroup_new) <- colnames(aggrGroup)

aggrGroup <- aggrGroup_new #4133 x 50 #2982 x 50



##data.frame to SimpleList
groupList <- lapply( seq_len(nrow(aggrGroup)),FUN = function(x){  as.character(unlist(aggrGroup[x,])) }   )  
names(groupList) <- paste0('aggr',1:length(groupList) )
groupList <- SimpleList( groupList )



# all.equal( 
#     unlist( groupList[['placenta_donor2#GCAGCTGTCAATCGTG-1']] ), as.character(unlist(aggrGroup['placenta_donor2#GCAGCTGTCAATCGTG-1',])), 
#     check.attributes = FALSE 
# ) #TRUE

##all.equal(names(groupList),rownames(aggrGroup) ) #TRUE

#groupList <- as.list(as.data.frame(t(aggrGroup)) ) #the same as apply method
#all.equal(groupList, SimpleList(groupList.1), check.attributes = FALSE) #TRUE

length(unique(unlist(aggrGroup) ) ) #22782 #11091 of total 11093
length(unique(unlist(groupList) ) ) #22782 #11091

lens_aggre <- unlist(lapply(groupList,FUN = length))
table(lens_aggre)
  50 
4133 

 50 
2982 

#saveRDS(groupList,paste('gene_peak_links/','GroupList.rds',sep='') )
groupList <- readRDS(paste('gene_peak_links/','GroupList.rds',sep='') )



######2 start to do peak mat and imputate (expr) mat correlation #######

#Check Chromosomes
chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
chrij <- intersect(chri, chrj)

#get Features
geneDF <- mcols(geneStart)
peakDF <- mcols(peakSet)
geneDF$seqnames <- seqnames(geneStart) #15294 #17185
peakDF$seqnames <- seqnames(peakSet) #274189 #299537

peakDF$name <- paste0(seqnames(peakSet),":",start(peakSet),'-',end(peakSet))



#saveRDS(geneDF,paste('gene_peak_links/','geneDF.rds',sep='') )
#saveRDS(peakDF,paste('gene_peak_links/','peakDF.rds',sep='') )



####get aggreGroup mat data

aggreGroupMatrix <- function(groupList = groupList,type = 'impute', featureDF = featureDF, useMatrix = useMatrix, files = files){

#     outdir <- paste('correlate_trajectory/',useMatrix,sep="")
#     if(dir.exists(outdir) ){}else{dir.create(outdir)}

    cat('use matrix:',useMatrix,' use file:',files[[useMatrix]],'\n')
    maty <- readRDS(files[[useMatrix]]) #18316 x 22013 #54662 x 11074
    #23032 x 11074, gene imputation RNA matrix, GeneIntegrateMatrix

    sum(duplicated(featureDF$name) )#0
    #24 #7327

    ##dupid <- featureDF$name[duplicated(featureDF$name)]
    #'RGS5''TBCE''PDE11A''LINC01238''PRSS50''CYB561D2''ATXN7''TXNRD3NB''CCDC39''MATR3''SOD2''POLR2J3''ABCF2''TMSB15B''PINX1''LINC01505''IGF2''HSPA14''EMG1''DIABLO''LINC02203''COG8''SCO2''H2BFS'

    ##table(featureDF$name[duplicated(featureDF$name)] %in% rownames(maty))

    # FALSE  TRUE 
    #   233  7094 

    # FALSE  TRUE 
    #    10    14

    ##remove duplication
    featureDF <- featureDF[!duplicated(featureDF$name) ,]
    stopifnot(sum(duplicated(featureDF$name)) == 0 ) #0

    rownames(featureDF) <- featureDF$name

    stopifnot(sum(duplicated(rownames(maty)) )  == 0)#0

    table(rownames(maty) %in% featureDF$name )
    #FALSE  TRUE 
    #3022 15294
    
    #FALSE  TRUE 
    #5847 17185
    
    maty <- maty[rownames(maty) %in% featureDF$name,]
    #15294 x 22013
    
    #match(dupid, rownames(mat)) #match will reture the first of the multiple positions
    #countMatches(dupid, rownames(mat))

    featureDF.sel <- featureDF[rownames(maty),]

    stopifnot(all.equal(rownames(maty),rownames(featureDF.sel)) ) #TRUE
    featureDF <- featureDF.sel


    saveRDS(featureDF,paste('gene_peak_links','/featureDF.',type,'.rds',sep='') )
    ##saveRDS(featureDF,'correlate_trajectory/featureDF.genes.RNAv93.rds')
    #saveRDS(maty,'correlate_trajectory/maty.geneExpressionMatrix.rds')
    saveRDS(maty,paste('gene_peak_links','/maty.',type,'.rds',sep='') )
    
    table (unique (unlist(groupList) ) %in% colnames(maty) )
    #FALSE  TRUE 
    #772 22010
    
    #FALSE  TRUE 
    #19 11072
    
    ####aggregate mat by pseudotime bin###
    matAggr <- matrix(0, nrow = nrow(featureDF), ncol = length(groupList))
    ##matAggr <- matrix(0, nrow = nrow(maty), ncol = length(groupList)) #for motifmat
    ##matAggr <- matrix(0, nrow = ncol(maty), ncol = length(groupList)) #for peakmat
    
    colnames(matAggr) <- names(groupList)
    rownames(matAggr) <- rownames(featureDF)
    #rownames(matAggr) <- rownames(maty) #for motifmat
    #rownames(matAggr) <- colnames(maty) #for peakmat
    #15294 × 4133
    
    pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
    for (z in seq_along(groupList)) {
      setTxtProgressBar(pb,round(z*100/length(groupList),0))
      cellsGroupz <- groupList[[z]]
      idx <- BiocGenerics::which(colnames(maty) %in%
      ##idx <- BiocGenerics::which(rownames(maty) %in% #for peakmat
        cellsGroupz)
      if (length(idx) > 0) {
        matAggr[, z] <- matAggr[, z] + Matrix::rowSums(maty[,idx
          , drop = FALSE])
      }else{message('wrong')}
    }

    #matAggr <- matAggr[rownames(featureDF), , drop = FALSE]
    
    stopifnot(all.equal(rownames(matAggr),rownames(maty)) )#TRUE
    #stopifnot(all.equal(rownames(matAggr),colnames(maty)) ) #for peakmat
    return(matAggr)
    
}



aggreGroupMatrix_reverse <- function(groupList = groupList,type = 'impute', featureDF = featureDF, useMatrix = useMatrix, files = files){

#     outdir <- paste('correlate_trajectory/',useMatrix,sep="")
#     if(dir.exists(outdir) ){}else{dir.create(outdir)}

    cat('use matrix:',useMatrix,' use file:',files[[useMatrix]],'\n')
    maty <- readRDS(files[[useMatrix]]) #54662 x 11074, 
    #11093 x 299537,peak reversed matrix, PeakMatrix

    sum(duplicated(featureDF$name) )
    #24 #7327

    ##dupid <- featureDF$name[duplicated(featureDF$name)]
    #'RGS5''TBCE''PDE11A''LINC01238''PRSS50''CYB561D2''ATXN7''TXNRD3NB''CCDC39''MATR3''SOD2''POLR2J3''ABCF2''TMSB15B''PINX1''LINC01505''IGF2''HSPA14''EMG1''DIABLO''LINC02203''COG8''SCO2''H2BFS'

    ##table(featureDF$name[duplicated(featureDF$name)] %in% rownames(maty))

    # FALSE  TRUE 
    #   233  7094 

    # FALSE  TRUE 
    #    10    14

    ##remove duplication
    featureDF <- featureDF[!duplicated(featureDF$name) ,]
    stopifnot(sum(duplicated(featureDF$name)) == 0 ) #0

    rownames(featureDF) <- featureDF$name

    stopifnot(sum(duplicated(colnames(maty)) )  == 0)#0

    #table(colnames(maty) %in% featureDF$name )
    #TRUE 
    #299537
    
    all.equal(colnames(maty) , featureDF$name) #TRUE
    
    maty <- maty[,colnames(maty) %in% featureDF$name]

    #match(dupid, rownames(mat)) #match will reture the first of the multiple positions
    #countMatches(dupid, rownames(mat))

    featureDF.sel <- featureDF[colnames(maty),]

    stopifnot(all.equal(colnames(maty),rownames(featureDF.sel)) ) #TRUE
    featureDF <- featureDF.sel


    saveRDS(featureDF,paste('gene_peak_links','/featureDF.',type,'.rds',sep='') )
    ##saveRDS(featureDF,'correlate_trajectory/featureDF.genes.RNAv93.rds')
    #saveRDS(maty,'correlate_trajectory/maty.geneExpressionMatrix.rds')
    saveRDS(maty,paste('gene_peak_links','/maty.',type,'.rds',sep='') )
    
    table (unique (unlist(groupList) ) %in% rownames(maty) )
    #TRUE 
    #11091 
    
    ####aggregate mat by pseudotime bin###
    #matAggr <- matrix(0, nrow = nrow(featureDF), ncol = length(groupList))
    ##matAggr <- matrix(0, nrow = nrow(maty), ncol = length(groupList)) #for motifmat
    matAggr <- matrix(0, nrow = ncol(maty), ncol = length(groupList)) #for peakmat
    
    colnames(matAggr) <- names(groupList)
    #rownames(matAggr) <- rownames(featureDF)
    #rownames(matAggr) <- rownames(maty) #for motifmat
    rownames(matAggr) <- colnames(maty) #for peakmat
    
    pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
    for (z in seq_along(groupList)) {
      setTxtProgressBar(pb,round(z*100/length(groupList),0))
      cellsGroupz <- groupList[[z]]
      #idx <- BiocGenerics::which(colnames(maty) %in%
      idx <- BiocGenerics::which(rownames(maty) %in% #for peakmat
        cellsGroupz)
      if (length(idx) > 0) {
        #matAggr[, z] <- matAggr[, z] + Matrix::rowSums(maty[,idx
        matAggr[, z] <- matAggr[, z] + Matrix::colSums(maty[idx,
          , drop = FALSE])
      }else{message('wrong')}
    }

    #matAggr <- matAggr[rownames(featureDF), , drop = FALSE]
    
    #stopifnot(all.equal(rownames(matAggr),rownames(maty)) )#TRUE
    stopifnot(all.equal(rownames(matAggr),colnames(maty)) ) #for peakmat
    return(matAggr)
    
}


##get grouped mat (use cicero method)
featureDF <- geneDF
useMatrix = "GeneIntegrateMatrix"
groupMatRNA <- aggreGroupMatrix(groupList = groupList,type = 'impute', featureDF = featureDF, useMatrix = useMatrix, files = files)
rawMatRNA <- groupMatRNA
15294 x 4133

#17185 x 2982

saveRDS(groupMatRNA,'gene_peak_links/groupMatRNA.rds')


featureDF <- peakDF
useMatrix <- 'PeakMatrix'
groupMatATAC <- aggreGroupMatrix_reverse(groupList = groupList,type = 'ciselement', featureDF = featureDF, useMatrix = useMatrix, files = files)
rawMatATAC <- groupMatATAC
274189 x 4133

#299537 x 2982

saveRDS(groupMatATAC,'gene_peak_links/groupMatATAC.rds')


all.equal (colnames(groupMatRNA),colnames(groupMatATAC) ) #TRUE

c('PAPPA','FLT1') %in% rownames(groupMatRNA) #TRUE



# #Check Chromosomes
# chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
# chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
# chrij <- intersect(chri, chrj)

# #Features
# geneDF <- mcols(geneStart)
# peakDF <- mcols(peakSet)
# geneDF$seqnames <- seqnames(geneStart)
# peakDF$seqnames <- seqnames(peakSet)



# #Group Matrix RNA
# #.logDiffTime(main="Getting Group RNA Matrix", t1=tstart, verbose=verbose, logFile=logFile)
# groupMatRNA <- ArchR:::.getGroupMatrix(
#     ArrowFiles = getArrowFiles(ArchRProj), 
#     featureDF = geneDF, 
#     groupList = knnObj, 
#     useMatrix = useMatrix, #GeneIntegrationMatrix
#     threads = 12,#threads,
#     verbose = FALSE
# )
# rawMatRNA <- groupMatRNA
# #.logThis(groupMatRNA, "groupMatRNA", logFile = logFile)

# #Group Matrix ATAC
# #.logDiffTime(main="Getting Group ATAC Matrix", t1=tstart, verbose=verbose, logFile=logFile)
# groupMatATAC <- ArchR:::.getGroupMatrix(
#     ArrowFiles = getArrowFiles(ArchRProj), 
#     featureDF = peakDF, 
#     groupList = knnObj, 
#     useMatrix = "PeakMatrix",
#     threads = 12,#threads,
#     verbose = TRUE#FALSE
# )#252545 × 497
# rawMatATAC <- groupMatATAC
# #.logThis(groupMatATAC, "groupMatATAC", logFile = logFile)

# #.logDiffTime(main="Normalizing Group Matrices", t1=tstart, verbose=verbose, logFile=logFile)


groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo
groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo

if(log2Norm){
  groupMatRNA  <- log2(groupMatRNA + 1)
  groupMatATAC <- log2(groupMatATAC + 1)    
}

dim(groupMatRNA)
15294 x 4133

#17185 x 2982
#17838 x 497

dim(groupMatATAC )
274189 x 4133

#299537 x 2982
#252545 x 497



#create se object for atac and rna
names(geneStart) <- NULL

seRNA <- SummarizedExperiment(
  assays = SimpleList(RNA = groupMatRNA, RawRNA = rawMatRNA), 
  rowRanges = geneStart
)
#metadata(seRNA)$KNNList <- knnObj
metadata(seRNA)$KNNList <- groupList
#.logThis(seRNA, "seRNA", logFile = logFile)

names(peakSet) <- NULL

seATAC <- SummarizedExperiment(
  assays = SimpleList(ATAC = groupMatATAC, RawATAC = rawMatATAC), 
  rowRanges = peakSet
)
#metadata(seATAC)$KNNList <- knnObj
metadata(seATAC)$KNNList <- groupList
#.logThis(seATAC, "seATAC", logFile = logFile)

rm(groupMatRNA, groupMatATAC)
gc()


#Save Group Matrices

saveRDS(seATAC,'gene_peak_links/seATAC-Group-KNN.rds')
saveRDS(seRNA,'gene_peak_links/seRNA-Group-KNN.rds')




#Overlaps 
#.logDiffTime(main="Finding Peak Gene Pairings", t1=tstart, verbose=verbose, logFile=logFile)
o <- DataFrame(
    findOverlaps( #function in rangedSummarizedExperiment obj
      resize(seRNA, 2 * maxDist + 1, "center"), #query, the same if use rowRanges(seRNA)
      #ArchR:::.suppressAll(resize(seRNA, 2 * maxDist + 1, "center")), 
      resize(rowRanges(seATAC), 1, "center"), #subject
      ignore.strand = TRUE
    )
)

#Get Distance from Fixed point A B 
o$distance <- GenomicRanges::distance(rowRanges(seRNA)[o[,1]] , rowRanges(seATAC)[o[,2]] ) #function in summarizedExperiment
colnames(o) <- c("B", "A", "distance")

#Null Correlations, Background Correlations
if(addEmpiricalPval){
  #.logDiffTime(main="Computing Background Correlations", t1=tstart, verbose=verbose, logFile=logFile)
  nullCor <- ArchR:::.getNullCorrelations(seATAC, seRNA, o, 1000)
}

##Computing Correlations
#.logDiffTime(main="Computing Correlations", t1=tstart, verbose=verbose, logFile=logFile)
o$Correlation <- ArchR:::rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(seATAC)-2))) #T-statistic P-value
o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC) - 2)
o$FDR <- p.adjust(o$Pval, method = "fdr")
out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")  
mcols(peakSet) <- NULL
names(peakSet) <- NULL
metadata(out)$peakSet <- peakSet
metadata(out)$geneSet <- geneStart
#1141795
#1333854

if(addEmpiricalPval){
  out$EmpPval <- 2*pnorm(-abs(((out$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
  out$EmpFDR <- p.adjust(out$EmpPval, method = "fdr")
}



# dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
# outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
# .safeSaveRDS(seATAC, outATAC, compress = FALSE)
# outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")
# .safeSaveRDS(seRNA, outRNA, compress = FALSE)
# metadata(out)$seATAC <- outATAC
# metadata(out)$seRNA <- outRNA

metadata(out)$seATAC <- 'gene_peak_links/seATAC-Group-KNN.rds'
metadata(out)$seRNA <- 'gene_peak_links/seRNA-Group-KNN.rds'

#saveRDS(out,'gene_peak_links/out.Peak2GeneLinks.rds') #all gene to peak link

saveRDS(out,'gene_peak_links/out.p2g.raw.rds')
out <- readRDS('gene_peak_links/out.p2g.raw.rds')
#out <- readRDS('result_good/gene_peak_links/out.p2g.raw.rds')

#metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out



#saveRDS(ArchRProj,'results_step_by_step/gene_peak_links/ArchRProj.rds')

#.logDiffTime(main="Completed Peak2Gene Correlations!", t1=tstart, verbose=verbose, logFile=logFile)
#.endLogging(logFile = logFile)

#ArchRProj



#####filter and output gene-peak link result
corCutOff = 0.45
FDRCutOff = 0.0001
varCutOffATAC = 0.25
varCutOffRNA = 0.25
resolution = 1
returnLoops = TRUE


#filter by cor >= 0.45, FDR <= 1e-4, var >= 0.25 of both RNA and ATAC


#p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
p2g <- out #1333854
p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
#112511
#96350

if(!is.null(varCutOffATAC)){
  p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
}

if(!is.null(varCutOffRNA)){
  p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
}

p2g #100682 #93053 x 8 #167108 x 8

#saveRDS(p2g,'gene_peak_links/out.Peak2GeneLinks.filter.rds')
saveRDS(p2g,'gene_peak_links/out.p2g.filter.rds')

p2g <- readRDS('gene_peak_links/out.p2g.filter.rds')

gene_p2g <- table(p2g$gene )
gene_p2g['PAPPA'] #65 #57
gene_p2g['FLT1'] #61 #63

gene_p2g['ERVFRD-1'] #3


mean(gene_p2g) #10 average p2g per gene



p2g.sel <- p2g[p2g$gene == 'PAPPA',] #65 row
p2g.sel <- p2g[p2g$gene == 'ERVFRD-1',] 
peakSummits.sel <- resize(metadata(p2g.sel)$peakSet, 1, "center")
geneStarts.sel <- resize(metadata(p2g.sel)$geneSet, 1, "start")

if(!is.null(resolution)){
  summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
  geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
}else{
  summitTiles <- start(peakSummits)
  geneTiles <- start(geneTiles)
}



if(returnLoops){ #for plot track 

  peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
  geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")

  if(!is.null(resolution)){
    summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
    geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
  }else{
    summitTiles <- start(peakSummits)
    geneTiles <- start(geneTiles)
  }

  loops <- ArchR:::.constructGR(
    seqnames = seqnames(peakSummits[p2g$idxATAC]),
    start = summitTiles[p2g$idxATAC],
    end = geneTiles[p2g$idxRNA]
  )
  mcols(loops)$value <- p2g$Correlation
  mcols(loops)$FDR <- p2g$FDR

  loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
  loops <- unique(loops)
  loops <- loops[width(loops) > 0]
  loops <- sort(sortSeqlevels(loops))

  loops <- SimpleList(Peak2GeneLinks = loops)

}

loop_width <- width(loops$Peak2GeneLinks)
hist(loop_width)

saveRDS(loops,'gene_peak_links/loops.rds')




###::::::::::::::plot peak-gene link heatmap (peak-cell mat) with rna heatmap:::::::::###
#ArchRProj = NULL, 
corCutOff = 0.45
FDRCutOff = 0.0001
varCutOffATAC = 0.25
varCutOffRNA = 0.25
k = 25
nPlot = 25000
limitsATAC = c(-2, 2)
limitsRNA = c(-2, 2)
groupBy = 'cluster'#"Clusters"
palGroup = NULL
palATAC = paletteContinuous("solarExtra")
palRNA = paletteContinuous("blueYellow")
verbose = TRUE
returnMatrices = TRUE #FALSE
seed = 1
#logFile = createLogFile("plotPeak2GeneHeatmap")


#########################################
# Get Inputs
#########################################
#ccd <- getCellColData(ArchRProj, select = groupBy) #cluster

# cluster.df.add <- readRDS('../cluster.df.add.te.rds')
# ccd <- DataFrame(cluster = cluster.df.add[,'cluster'] )
# rownames(ccd ) <- rownames(cluster.df.add)


# #p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
# p2g <- readRDS('gene_peak_links/out.Peak2GeneLinks.rds') #reload the raw peak2gene link,1333584
# p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
# #96350 #201762

# if(!is.null(varCutOffATAC)){
#   p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
# }

# if(!is.null(varCutOffRNA)){
#   p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
# }

# if(nrow(p2g) == 0){
#   stop("No peak2genelinks found with your cutoffs!")
# }
# #93053 #167108 

# all.equal(p2g, readRDS('gene_peak_links/out.Peak2GeneLinks.filter.rds')) #TRUE

#p2g <- readRDS('gene_peak_links/out.Peak2GeneLinks.filter.rds')
p2g <- readRDS('gene_peak_links/out.p2g.filter.rds') #100682
#p2g <- readRDS('result_good/gene_peak_links/out.p2g.filter.rds') #93053 #filter by cor >= 0.45, FDR <= 1e-4, var >= 0.25 of both RNA and ATAC, but with peak to more than one gene

# if(!file.exists(metadata(p2g)$seATAC)){
#  stop("seATAC does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
# }
# if(!file.exists(metadata(p2g)$seRNA)){
#   stop("seRNA does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
# }


##filter mATAC and mRNA for filtered p2g links
mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ] #100682 x 4133 #93053 x 2982  
mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ] #100682 x 4133 #93053 x 2982 

#mATAC <- seATAC[p2g$idxATAC, ]
#mRNA <- seRNA[p2g$idxRNA, ]

p2g$peak <- paste0(rowRanges(mATAC))
p2g$gene <- rowData(mRNA)$name
gc()

mATAC <- assay(mATAC) #100682 x 4133 #93053 x 2982
mRNA <- assay(mRNA) #100682 x 4133 #93053 x 2982


#########################################
# Determine Groups from KNN
#########################################

cluster.df.add <-readRDS('../../../02.snapATAC_harmony/cluster.df.add.snapatac2.rds')
#cluster.df.add <-readRDS('../../../02.snapATAC_harmony/cluster.df.add.final.final.rds')

#cluster.df.add.atac <- subset(x = cluster.df.add, subset = type == 'atac')
#table(cluster.df.add.atac$cluster)

table(cluster.df.add$cluster)
1    7    9    2    6    4    5    8    3 
3937 1636  814 3594 3713 2493 2230 1447 2921 

#    1    2    3    4    5    6    7    8    9 
# 1832 1783 1090 1546 1475 1509 1019  444  395

plot(cluster.df.add[,c('UMAP_1','UMAP_2')],cex = 0.1)



##rename cluster.df.add from snapatac to get ccd for heatmap column annotation
cellid <- rownames(cluster.df.add)
sum(duplicated(cellid )) #0

cellid <- sapply(stringr::str_split(cellid,'early|:|-',n=4), function(x){paste0( 'placenta.atac_',x[3],'-',x[2] ) } )


# idx1 <- grep(pattern = 'donor1',x = cellid )
# idx2 <- grep(pattern = 'donor2',x = cellid )
# #unique(idx1+idx2)#1
# length(c(idx1,idx2)) #11093
    
# cellid[idx1] <- gsub(pattern = '_donor1#',replacement = '.atac_',x = cellid[idx1])

# cellid[idx2] <- gsub(pattern = '-1$',replacement = '-2',x = cellid[idx2])
# cellid[idx2] <- gsub(pattern = '_donor2#',replacement = '.atac_',x = cellid[idx2])

row.names(cluster.df.add) <- cellid

##get ccd
ccd <- DataFrame(cluster = cluster.df.add[,'cluster'] )
rownames(ccd ) <- rownames(cluster.df.add)




#.logDiffTime(main="Determining KNN Groups!", t1=tstart, verbose=verbose, logFile=logFile)
KNNList <- as(metadata(readRDS(metadata(p2g)$seRNA))$KNNList, "list") #rna aggregate group

all.equal(SimpleList(KNNList), groupList,check.attributes = FALSE) #TRUE


KNNGroups <- lapply(seq_along(KNNList), function(x){
  KNNx <- KNNList[[x]]
  names(sort(table(ccd[KNNx, 1, drop = TRUE]), decreasing = TRUE))[1]
}) %>% unlist #reprensentive cluster in each aggre-group

cD <- DataFrame(row.names=paste0("K_", seq_len(ncol(mATAC))), groupBy = KNNGroups) #raw colname annotation
##pal <- ArchR:::paletteDiscrete(values=gtools::mixedsort(unique(ccd[,1])))

pal <- unlist(map_cellcolor)

# if(!is.null(palGroup)){
#   pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% names(pal)]
# }
colorMap <- list(groupBy = pal)
attr(colorMap[[1]], "discrete") <- TRUE

#########################################
# Organize Matrices
#########################################
mATAC <- ArchR:::.rowZscores(mATAC)
mRNA <- ArchR:::.rowZscores(mRNA)

all.equal(rownames(mATAC),p2g$peak ) #TRUE
all.equal(rownames(mRNA),p2g$gene ) #TRUE

rownames(mATAC) <- NULL
rownames(mRNA) <- NULL
colnames(mATAC) <- paste0("K_", seq_len(ncol(mATAC)))
colnames(mRNA) <- paste0("K_", seq_len(ncol(mRNA)))
rownames(mATAC) <- paste0("P2G_", seq_len(nrow(mATAC)))
rownames(mRNA) <- paste0("P2G_", seq_len(nrow(mRNA)))
rownames(p2g) <- paste0("P2G_", seq_len(nrow(p2g)))

all.equal(rownames(mATAC),rownames(mRNA)) #TRUE
all.equal(rownames(mATAC), rownames(p2g)) #TRUE


##Ordering Peak2Gene Links (kmean method)
#.logDiffTime(main="Ordering Peak2Gene Links!", t1=tstart, verbose=verbose, logFile=logFile)
if(!is.null(seed)){
  set.seed(seed)
}

k1 <- kmeans(mATAC, k) #atac as reference, k= 25 #slow!!

if(nrow(mATAC) > nPlot){ #random select each k-group by propertion up to total 25000
    nPK <- nPlot * table(k1$cluster) / length(k1$cluster) #sample ratio
    splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
    kDF <- lapply(seq_along(splitK), function(x){
      idx <- sample(splitK[[x]], floor(nPK[x]))
      k <- rep(x, length(idx))
      DataFrame(k = k, idx = idx)
    }) %>% Reduce("rbind", .)
}else{
   kDF <- DataFrame(k = k1$cluster, idx = seq_len(nrow(mATAC)))
 }

saveRDS(k1,'gene_peak_links/k1.rds')

##determin the row k_group and column order (with hclust)
bS <- ArchR:::.binarySort(t(ArchR:::.groupMeans(t(mATAC[kDF[,2],]), kDF[,1])),  clusterCols = TRUE, cutOff = 1) #will do hclust in this .binarySort
rowOrder <- rownames(bS[[1]])
colOrder <- colnames(bS[[1]])
kDF[,3] <- as.integer(ArchR::mapLabels(paste0(kDF[,1]), newLabels = paste0(seq_along(rowOrder)), oldLabels = rowOrder)) #to sort with roworder

saveRDS(bS,'gene_peak_links/bS.rds')
saveRDS(colOrder,'gene_peak_links/colOrder.rds')


saveRDS(kDF,'gene_peak_links/kDF.rds')


###ordering peak2gene links (similar to trajectory plot)


if(returnMatrices){
    
    outMatch <- SimpleList(
        ATAC = SimpleList(
          matrix = mATAC[kDF[,2],colOrder],
          kmeansId = kDF[,3],
          colData = cD[colOrder,,drop=FALSE]
        ),
        RNA = SimpleList(
          matrix = mRNA[kDF[,2],colOrder],
          kmeansId = kDF[,3],
          colData = cD[colOrder,,drop=FALSE]
        ),
        Peak2GeneLinks = p2g[kDF[,2],]
      )

}


saveRDS(outMatch,'gene_peak_links/outMatch.gene-peak-link-paired.rds') #24988

#out <- readRDS('gene_peak_links/gene-peak-link-paired.mat.rds')

write.table(outMatch$Peak2GeneLinks,file = 'gene_peak_links/outMatch.gene-peak-link-paired.table.txt',sep='\t',quote = FALSE,row.names = TRUE, col.names = TRUE)

# linksSig <- peakLinks[which(peakLinks$sigCorrelation)]

# outMatch <- list(
#   seA = seA[unique(mcols(linksSig)$peakName),], 
#   seB = seB[unique(mcols(linksSig)$gene_name),], 
#   linksSig = linksSig,
#   linksAll = peakLinks
#   )



#########################################
# Plot Heatmaps
#########################################
#.logDiffTime(main="Constructing ATAC Heatmap!", t1=tstart, verbose=verbose, logFile=logFile)

all.equal (rownames(mATAC),rownames(mRNA) ) #TRUE

htATAC <- ArchR:::.ArchRHeatmap(
    mat = mATAC[kDF[,2],colOrder],
    scale = FALSE,
    limits = limitsATAC,
    color = palATAC, 
    colData = cD[colOrder,,drop=FALSE],
    colorMap = colorMap,
    clusterCols = FALSE,
    clusterRows = FALSE,
    split = kDF[,3],
    labelRows = FALSE,
    labelCols = FALSE,
    draw = FALSE,
    name = paste0("ATAC Z-Scores\n", nrow(mATAC), " P2GLinks")
)

#.logDiffTime(main = "Constructing RNA Heatmap!", t1 = tstart, verbose = verbose, logFile = logFile)
htRNA <- ArchR:::.ArchRHeatmap(
    mat = mRNA[kDF[,2],colOrder], 
    scale = FALSE,
    limits = limitsRNA,
    color = palRNA, 
    colData = cD[colOrder,,drop=FALSE],
    colorMap = colorMap,
    clusterCols = FALSE,
    clusterRows = FALSE,
    split = kDF[,3],
    labelRows = FALSE,
    labelCols = FALSE,
    draw = FALSE,
    name = paste0("RNA Z-Scores\n", nrow(mRNA), " P2GLinks")
)






#.endLogging(logFile = logFile)

pdf('gene_peak_links/p2g_heatmap/gene-peak-link-paired.mat.pdf',width = 6.5,height = 6.5,useDingbats = FALSE)
set.seed(1);ComplexHeatmap::draw(htATAC + htRNA) #slow
dev.off()


saveRDS(mATAC,'gene_peak_links/p2g_heatmap/mATAC.rds')
saveRDS(mRNA,'gene_peak_links/p2g_heatmap/mRNA.rds')

saveRDS(kDF,'gene_peak_links/p2g_heatmap/kDF.rds') #rowOrder and split group
#saveRDS(colOrder,'gene_peak_links/p2g_heatmap/colOrder.rds')
saveRDS(colOrder,'gene_peak_links/p2g_heatmap/colOrder.usesnapatac.rds')
saveRDS(cD,'gene_peak_links/p2g_heatmap/cD.usesnapatac.rds')
saveRDS(colorMap,'gene_peak_links/p2g_heatmap/colorMap.usesnapatac.rds')

saveRDS(p2g, 'gene_peak_links/p2g_heatmap/p2g.rds')




###look into details?
set.seed(1);res.row.order <- row_order(htATAC) 

res.row.order.ids <- names(res.row.order) #with splicing of row
res.row.order.flat <- unlist(res.row.order) 
#all.equal(res.row.order.combine,res.row.order.flat) #TRUE
names(res.row.order.flat) <- NULL
set.seed(1); res.column.order <- column_order(htATAC) #no need to set seed = 1, no kmeans method applied

rowid <- rownames(mATAC[kDF[,2],colOrder])[res.row.order.flat]
colid <- colnames(mATAC)[res.column.order]

##res.row.order.flat <- head(res.row.order.flat,n=120)





#################plot gene-peak link and gene expression with trajectory ?#####################

###reload for replot

mATAC <- readRDS('gene_peak_links/p2g_heatmap/mATAC.rds')
mRNA <- readRDS('gene_peak_links/p2g_heatmap/mRNA.rds')

kDF <- readRDS('gene_peak_links/p2g_heatmap/kDF.rds') #rowOrder and split group
colOrder<- readRDS('gene_peak_links/p2g_heatmap/colOrder.rds') #use liger cluster
#colOrder<- readRDS('gene_peak_links/p2g_heatmap/colOrder.usesnapatac.rds')

p2g <- readRDS( 'gene_peak_links/p2g_heatmap/p2g.rds')




matATAC <- mATAC
matRNA <- mRNA

all.equal (rownames(p2g), rownames(matATAC) ) #TRUE
all.equal (rownames(p2g), rownames(matRNA) ) #TRUE

p2g #100682 #93053
max(p2g$Correlation) #0.97 #0.95
min(p2g$Correlation) #0.45 #0.45

min(p2g$FDR) #0 #0
max(p2g$FDR) #1.8327 #1.0138


##how many genes?
length(table(p2g$gene) ) #9785 #9922
hist(table(p2g$gene),breaks = 50)

##unique peaks ?
sum(duplicated(p2g$peak) ) 
#53940
#46495


res.stat <- table(p2g$gene,p2g$peak)

#9922 x 46558 unique 

peak_link_cnt <- Matrix::colSums(res.stat)

peak_link_cnt[peak_link_cnt > 1 ]

table(peak_link_cnt > 1) 

FALSE  TRUE 
22304 24438 

#FALSE  TRUE 
#23489 23069 





##############output p2g junction file (both raw and filter for the most correlation p2g link)##############



p2g.filter <- readRDS('gene_peak_links/out.p2g.filter.rds')#100682
#p2g.filter <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/liger/peak_gene_link/result_good/gene_peak_links/out.p2g.filter.rds') #93053
p2g.filter <- as.data.frame(p2g.filter)

all.equal(p2g.filter[,1:8],as.data.frame(p2g[,1:8]),check.attributes = FALSE) #TRUE

all.equal(p2g.filter, as.data.frame(p2g),check.attributes = FALSE) #FALSE

sum(duplicated(p2g.filter$peak))
53940
#46495

p2g.filter$rowid <- rownames(p2g.filter) 

p2g.filter.uniq <- p2g.filter %>% dplyr::group_by(peak) %>% dplyr::slice(which.max(Correlation))
#46742
#46558

p2g.filter.uniq <- as.data.frame(p2g.filter.uniq)

sum(duplicated(p2g.filter.uniq$peak)) #0

#p2g.filter.uniq <- p2g.filter.uniq[sort(p2g.filter.uniq$rowid),]

length(unique(p2g.filter.uniq$gene) )
#7325
#7678

rownames(p2g.filter.uniq) <- p2g.filter.uniq$peak

options(repr.plot.width = 4.5, repr.plot.height = 4.5)
hist(p2g.filter.uniq$Correlation,breaks = 30)

saveRDS(p2g.filter.uniq,'gene_peak_links/out.p2g.filter.uniq.rds')



p2gTab.split <- split(p2g.filter, f = p2g.filter$gene ) #9785 #9922 #10359 genes with at least one peak-to-gene links


p2gTab.split <- split(p2g.filter.uniq, f = p2g.filter.uniq$gene ) #7325 #7678 genes with one to one peak-to-gene links

p2g.len<- sort(sapply(X = p2gTab.split, FUN = function(x){nrow(x)} ),decreasing = TRUE)
p2g.len.df <- data.frame(Gene=names(p2g.len),n=p2g.len, stringsAsFactors = FALSE)

options(repr.plot.width = 7.5, repr.plot.height = 7.5)
hist(p2g.len.df$n,breaks = 100)


p2g.len.df[grep('PAPPA',p2g.len.df$Gene),]

PAPPA2	53
PAPPA	21

# PAPPA	50
# PAPPA2	62

p2g.len.df[grep('FLT',p2g.len.df$Gene),]
FLT1	55
FLT4	10

FLT1	57
FLT4	10

p2g.len.df[grep('CSH',p2g.len.df$Gene),]
CSH2	8
CSHL1	3
PRKCSH	2


# CSH1	3
# CSH2	1
# CSHL1	1

p2g.len.df[grep('GH[2|R]',p2g.len.df$Gene),]
GH2	GH2	14
GHR	GHR	9
GHRL	GHRL	2

# GH2	17
# GHR	15


p2g.len.df[grep('ERV',p2g.len.df$Gene),]
ERVV-1	4
ERVFRD-1	3
ERVV-2	3
ERVW-1	1

# ERVV-1	8
# ERVFRD-1	4
# ERVV-2	3
# ERVW-1	2

p2g.len.df[grep('STAT4',p2g.len.df$Gene),]
STAT4	14

#STAT4	19

p2g.len.df[grep('STAT5A',p2g.len.df$Gene),]
STAT5A	7

#STAT5A	5

p2g.len.df[grep('MYCN',p2g.len.df$Gene),]
MYCNUT	15
MYCN	1

#MYCN	1



qs <- quantile(p2g.len.df$n,probs = seq(from = 0,to = 1,by = 0.1))
0% 1 10% 1 20% 1 30% 2 40% 3 50% 4 60% 5 70% 7 80% 10 90% 15 100% 92

#0% 1 10% 1 20% 1 30% 2 40% 3 50% 4 60% 5 70% 7 80% 10 90% 14 100% 97

#0% 1 10% 2 20% 3 30% 4 40% 6 50% 8 60% 9 70% 12 80% 15 90% 19 100% 104

p2gGenes <- p2g.len.df[p2g.len.df$n >= qs['80%'],'Gene']
length(p2gGenes)
1633

#1555 : cutoff = 10, 80%, with p2g.filter.uniq

#4987 : cutoff = 8 50%
#2021 : cutoff = 15 80%




p2gGenes #1633
p2gTab.split #7325


table( p2gGenes %in% names(p2gTab.split) )
#1633 TRUE
#1555 TRUE in total 7678 genes

##table( p2gGenes %in% names(p2gTab.split) )
#481 TRUE in total 10359 genes

##get TSS of gene from 

tss <- read.table('/home/mjwang/dataex/cellranger_genomes/refdata-cellranger-atac-GRCh38-1.1.0/regions/tss.bed',header = FALSE, sep = '\t', stringsAsFactors = FALSE) #slightly different with cellranger-atac.1.2

##tss to identify p2g use featureDF from GeneIntegrationMat plotting featureDF, from gff of snapatac2


colnames(tss) <- c('chr','start','end','geneid','score','strand','type')

tss <- geneStart #used to get p2g (15294, from GeneIntegrationMat featureDF)

geneStart.df <- as.data.frame(geneStart, stringsAsFactors = FALSE)
colnames(geneStart.df) <- c('chr','start','end','width','strand','geneid','idx')
geneStart.df <- geneStart.df[,c('chr','start','end','geneid')]

#subset(geneStart.df,geneid == 'FLT1')
sum(duplicated(geneStart.df$geneid )) #0

##output p2g link bed file

write.table(
            #x = 'track name=junctions description="p2g (1555 genes) peak to gene links"',
            #file = 'gene_peak_links/p2g.p2gGenes.junction.bed',
            #x = 'track name=junctions description="p2g (all 7678 genes) peak to gene links"',
            x = 'track name=junctions description="p2g (all 7325 genes) peak to gene links"',
            file = 'gene_peak_links/p2g.all.junction.bed',
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t',
            quote=FALSE
           )

p2g.table.all <- data.frame()

#for(i in p2gGenes){
for(i in names(p2gTab.split)){
    #cat('p2g gene peak to gene links for ',i,'\n',sep='')
    p2g.df <- p2gTab.split[[i]]
    stopifnot(nrow(p2g.df) > 0)
    
    ##define tss region
    tss.gene <- subset(tss,geneid == i)
    #stopifnot(nrow(tss.gene) > 0)
    
    if(nrow(tss.gene) == 0){
      cat('can not found tss for gene ',i,'\n',sep='')
      next
    }
    
    if(nrow(tss.gene) == 1){ 
      p2g.df.ext <- data.frame(p2g.df, tss=paste0(tss.gene$chr,':',tss.gene$start,'-',tss.gene$end) ,stringsAsFactors = FALSE)
    
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
        p2g.df.ext <- data.frame(p2g.df, tss=paste0(tss.gene$chr,':',tss.gene$start,'-',tss.gene$end) ,stringsAsFactors = FALSE)
        
        
    }
    
    
    ##igv junction format
    #track name=junctions description="Cicero_Granja connections.filter3.0 links"
    #chr1    634199  778977  link1   0.305989577883302       +
    #chr1    634199  845876  link2   0.435166666059608       +
    #chr1    634199  877093  link3   0.370885687189481       +
    p2g.df.ext.expand <- p2g.df.ext %>% 
      tidyr::separate(peak,c('chr_peak','start_peak','end_peak')) %>% 
      tidyr::separate(tss,c('chr_tss','start_tss','end_tss'))
    
    #stopifnot( all(p2g.df.ext.expand$chr_peak == p2g.df.ext.expand$chr_tss) )
    
    if(!all(p2g.df.ext.expand$chr_peak == p2g.df.ext.expand$chr_tss)){
        cat('chr diffs: ',p2g.df.ext.expand$chr_peak[1],'  ',p2g.df.ext.expand$chr_tss[1],' ',p2g.df.ext.expand$gene[1],'\n',sep='')
        next
    }
    
    p2g.table <- data.frame(chr=p2g.df.ext.expand$chr_peak,
               start =  ceiling(0.5*(as.numeric(p2g.df.ext.expand$start_peak)
                                     +as.numeric(p2g.df.ext.expand$end_peak)
                                    )
                               ),
               end =  ceiling(0.5*(as.numeric(p2g.df.ext.expand$start_tss)
                                     +as.numeric(p2g.df.ext.expand$end_tss)
                                    )
                               ),
               id = paste('p2g_',i,'_link',1:nrow(p2g.df.ext.expand),sep=''),
               #score = -1*log10(p2g.df.ext.expand$rObs ),
               score = p2g.df.ext.expand$Correlation,
               strand = '+'
            )
    
    p2g.table.all <- rbind.data.frame(p2g.table.all,p2g.table)
    
    write.table(x = p2g.table,
            file = 'gene_peak_links//p2g.all.junction.bed',
            #file = 'gene_peak_links//p2g.p2gGenes.junction.bed',
            row.names = FALSE,
            col.names = FALSE,
            sep='\t',
            quote=FALSE,
            append = TRUE
           )
    
}


###plot p2g peak to tss distance
quantile(abs(p2g.table.all$end-p2g.table.all$start),probs = seq(0,1,0.1) )
#0% 0 10% 10154 20% 26624 30% 45690.6 40% 67990.4 50% 93231.5 60% 120762 70% 150045.3 80% 182929.4 90% 216893.8 100% 1158946
hist(abs(p2g.table.all$end-p2g.table.all$start),breaks = 100)


can not found tss for gene AADACL2-AS1
can not found tss for gene ABHD15-AS1
can not found tss for gene ABHD17C
can not found tss for gene ABRAXAS1
can not found tss for gene ABRAXAS2
can not found tss for gene ADAMTS19-AS1
can not found tss for gene AFDN
can not found tss for gene AFDN-DT
can not found tss for gene AFG1L
can not found tss for gene ALOX12-AS1
can not found tss for gene ARHGAP26-AS1
can not found tss for gene ARHGAP45
can not found tss for gene ARHGEF3-AS1
can not found tss for gene ARMH4
can not found tss for gene ATP5F1A
can not found tss for gene ATP5F1B
can not found tss for gene ATP5F1C
can not found tss for gene ATP5IF1
can not found tss for gene ATP5MC2
can not found tss for gene ATP5MC3
can not found tss for gene ATP5ME
can not found tss for gene ATP5PF
can not found tss for gene BASP1-AS1
can not found tss for gene BBOX1-AS1
can not found tss for gene BCLAF3
can not found tss for gene BDNF-AS
can not found tss for gene BFSP2-AS1
can not found tss for gene BHLHE40-AS1
can not found tss for gene BICDL2
can not found tss for gene BICRA-AS1
can not found tss for gene BMT2
can not found tss for gene C10orf143
can not found tss for gene C19orf48
can not found tss for gene C21orf91-OT1
can not found tss for gene C2CD6
can not found tss for gene CAPN10-DT
can not found tss for gene CARMIL1
can not found tss for gene CASC9
can not found tss for gene CASTOR1
can not found tss for gene CASTOR2
can not found tss for gene CATSPERE
can not found tss for gene CAVIN1
can not found tss for gene CAVIN4
can not found tss for gene CCAR2
can not found tss for gene CCDC144NL-AS1
can not found tss for gene CEMIP2
can not found tss for gene CENPX
can not found tss for gene CEP128
can not found tss for gene CEP83
can not found tss for gene CFAP298
can not found tss for gene CFAP299
can not found tss for gene CHODL-AS1
can not found tss for gene CHRM3-AS2
can not found tss for gene CLDN10-AS1
can not found tss for gene CLRN1-AS1
can not found tss for gene COA4
can not found tss for gene COLGALT1
can not found tss for gene COP1
can not found tss for gene COQ8B
can not found tss for gene CPLANE1
can not found tss for gene CRACR2B
can not found tss for gene CXXC4-AS1
can not found tss for gene CYP1B1-AS1
can not found tss for gene DDX11-AS1
can not found tss for gene DELE1
can not found tss for gene DEPP1
can not found tss for gene DIAPH2-AS1
can not found tss for gene DIRC3-AS1
can not found tss for gene DLGAP1-AS1
can not found tss for gene DLGAP1-AS2
can not found tss for gene DLGAP4-AS1
can not found tss for gene DM1-AS
can not found tss for gene DMAC1
can not found tss for gene DPEP2NB
can not found tss for gene DPH6-DT
can not found tss for gene DUBR
can not found tss for gene EGFR-AS1
can not found tss for gene EIPR1
can not found tss for gene ELMO1-AS1
can not found tss for gene ELP1
can not found tss for gene ENTPD3-AS1
can not found tss for gene EPHA1-AS1
can not found tss for gene EPN2-AS1
can not found tss for gene ESS2
can not found tss for gene EWSAT1
can not found tss for gene F11-AS1
can not found tss for gene FAM214A
can not found tss for gene FAM222A-AS1
can not found tss for gene FAM53B-AS1
can not found tss for gene FBH1
can not found tss for gene FLG-AS1
can not found tss for gene FRG1-DT
can not found tss for gene FRMD6-AS2
can not found tss for gene FSIP2-AS1
can not found tss for gene FYB2
can not found tss for gene GAS1RR
can not found tss for gene GATA3-AS1
can not found tss for gene GATD1
can not found tss for gene GMDS-DT
can not found tss for gene GPAT4
can not found tss for gene GPRC5D-AS1
can not found tss for gene GRAMD2A
can not found tss for gene GRAMD2B
can not found tss for gene GRK3
can not found tss for gene GSDME
can not found tss for gene GSEC
can not found tss for gene GVQW3
can not found tss for gene HDAC2-AS2
can not found tss for gene HIKESHI
can not found tss for gene HLX-AS1
can not found tss for gene HMGN3-AS1
can not found tss for gene HPF1
can not found tss for gene HYKK
can not found tss for gene IDH2-DT
can not found tss for gene IDI2-AS1
can not found tss for gene IFT43
can not found tss for gene IGF2-AS
can not found tss for gene INAVA
can not found tss for gene INHBA-AS1
can not found tss for gene IQCH-AS1
can not found tss for gene IQCM
can not found tss for gene ISM1
can not found tss for gene ITGB1-DT
can not found tss for gene ITGB5-AS1
can not found tss for gene JCAD
can not found tss for gene JHY
can not found tss for gene JPT1
can not found tss for gene JPT2
can not found tss for gene KANTR
can not found tss for gene KAT14
can not found tss for gene KCNMB2-AS1
can not found tss for gene KDM4A-AS1
can not found tss for gene KIF26B-AS1
can not found tss for gene KIF9-AS1
can not found tss for gene KNL1
can not found tss for gene KRTAP5-AS1
can not found tss for gene KYAT3
can not found tss for gene LARGE1
can not found tss for gene LARGE2
can not found tss for gene LCMT1-AS2
can not found tss for gene LENG8-AS1
can not found tss for gene LGALS8-AS1
can not found tss for gene LHFPL6
can not found tss for gene LINC00221
can not found tss for gene LINC00222
can not found tss for gene LINC00240
can not found tss for gene LINC00265
can not found tss for gene LINC00276
can not found tss for gene LINC00348
can not found tss for gene LINC00393
can not found tss for gene LINC00456
can not found tss for gene LINC00474
can not found tss for gene LINC00486
can not found tss for gene LINC00571
can not found tss for gene LINC00578
can not found tss for gene LINC00601
can not found tss for gene LINC00607
can not found tss for gene LINC00621
can not found tss for gene LINC00623
can not found tss for gene LINC00630
can not found tss for gene LINC00642
can not found tss for gene LINC00656
can not found tss for gene LINC00665
can not found tss for gene LINC00685
can not found tss for gene LINC00842
can not found tss for gene LINC00845
can not found tss for gene LINC00847
can not found tss for gene LINC00862
can not found tss for gene LINC00877
can not found tss for gene LINC00882
can not found tss for gene LINC00886
can not found tss for gene LINC00937
can not found tss for gene LINC00939
can not found tss for gene LINC00964
can not found tss for gene LINC00967
can not found tss for gene LINC01004
can not found tss for gene LINC01088
can not found tss for gene LINC01095
can not found tss for gene LINC01098
can not found tss for gene LINC01118
can not found tss for gene LINC01119
can not found tss for gene LINC01162
can not found tss for gene LINC01194
can not found tss for gene LINC01203
can not found tss for gene LINC01205
can not found tss for gene LINC01210
can not found tss for gene LINC01284
can not found tss for gene LINC01291
can not found tss for gene LINC01320
can not found tss for gene LINC01322
can not found tss for gene LINC01355
can not found tss for gene LINC01356
can not found tss for gene LINC01376
can not found tss for gene LINC01411
can not found tss for gene LINC01416
can not found tss for gene LINC01446
can not found tss for gene LINC01456
can not found tss for gene LINC01483
can not found tss for gene LINC01484
can not found tss for gene LINC01485
can not found tss for gene LINC01501
can not found tss for gene LINC01503
can not found tss for gene LINC01505
can not found tss for gene LINC01508
can not found tss for gene LINC01509
can not found tss for gene LINC01554
can not found tss for gene LINC01572
can not found tss for gene LINC01579
can not found tss for gene LINC01591
can not found tss for gene LINC01605
can not found tss for gene LINC01612
can not found tss for gene LINC01618
can not found tss for gene LINC01687
can not found tss for gene LINC01719
can not found tss for gene LINC01725
can not found tss for gene LINC01731
can not found tss for gene LINC01748
can not found tss for gene LINC01754
can not found tss for gene LINC01792
can not found tss for gene LINC01794
can not found tss for gene LINC01800
can not found tss for gene LINC01804
can not found tss for gene LINC01807
can not found tss for gene LINC01811
can not found tss for gene LINC01814
can not found tss for gene LINC01820
can not found tss for gene LINC01876
can not found tss for gene LINC01889
can not found tss for gene LINC01891
can not found tss for gene LINC01895
can not found tss for gene LINC01909
can not found tss for gene LINC01926
can not found tss for gene LINC01948
can not found tss for gene LINC01949
can not found tss for gene LINC01950
can not found tss for gene LINC01968
can not found tss for gene LINC01980
can not found tss for gene LINC01983
can not found tss for gene LINC02008
can not found tss for gene LINC02009
can not found tss for gene LINC02018
can not found tss for gene LINC02021
can not found tss for gene LINC02035
can not found tss for gene LINC02036
can not found tss for gene LINC02050
can not found tss for gene LINC02082
can not found tss for gene LINC02099
can not found tss for gene LINC02100
can not found tss for gene LINC02109
can not found tss for gene LINC02175
can not found tss for gene LINC02196
can not found tss for gene LINC02201
can not found tss for gene LINC02211
can not found tss for gene LINC02245
can not found tss for gene LINC02250
can not found tss for gene LINC02253
can not found tss for gene LINC02256
can not found tss for gene LINC02267
can not found tss for gene LINC02288
can not found tss for gene LINC02291
can not found tss for gene LINC02315
can not found tss for gene LINC02343
can not found tss for gene LINC02365
can not found tss for gene LINC02378
can not found tss for gene LINC02432
can not found tss for gene LINC02435
can not found tss for gene LINC02463
can not found tss for gene LINC02482
can not found tss for gene LINC02484
can not found tss for gene LINC02494
can not found tss for gene LINC02518
can not found tss for gene LINC02533
can not found tss for gene LINC02570
can not found tss for gene LIX1-AS1
can not found tss for gene LMCD1-AS1
can not found tss for gene LMO7-AS1
can not found tss for gene LRMDA
can not found tss for gene LRP4-AS1
can not found tss for gene LUCAT1
can not found tss for gene LY6E-DT
can not found tss for gene MAMDC2-AS1
can not found tss for gene MARF1
can not found tss for gene MCM3AP-AS1
can not found tss for gene MCPH1-AS1
can not found tss for gene MCRIP1
can not found tss for gene MCUB
can not found tss for gene MEF2C-AS1
can not found tss for gene METTL26
can not found tss for gene MIATNB
can not found tss for gene MINDY2
can not found tss for gene MINDY3
can not found tss for gene MIR2052HG
can not found tss for gene MIR222HG
can not found tss for gene MIR29B2CHG
can not found tss for gene MIR4300HG
can not found tss for gene MIR4435-2HG
can not found tss for gene MIR4713HG
can not found tss for gene MIR503HG
can not found tss for gene MIR5689HG
can not found tss for gene MIR646HG
can not found tss for gene MISP3
can not found tss for gene MME-AS1
can not found tss for gene MMP24OS
can not found tss for gene MRE11
can not found tss for gene MRM2
can not found tss for gene MRM3
can not found tss for gene MTERF2
can not found tss for gene MYCNUT
can not found tss for gene MYLK-AS1
can not found tss for gene MYO16-AS1
can not found tss for gene MYORG
can not found tss for gene N4BP3
can not found tss for gene NALCN-AS1
can not found tss for gene NCF4-AS1
can not found tss for gene NCK1-DT
can not found tss for gene NDUFAF8
can not found tss for gene NDUFB2-AS1
can not found tss for gene NDUFV2-AS1
can not found tss for gene NEPRO
can not found tss for gene NEURL1-AS1
can not found tss for gene NIPSNAP2
can not found tss for gene NKAPD1
can not found tss for gene NNT-AS1
can not found tss for gene NR2F2-AS1
can not found tss for gene NRSN2-AS1
can not found tss for gene NSG2
can not found tss for gene NUP50-DT
can not found tss for gene NUTM2A-AS1
can not found tss for gene ODAPH
can not found tss for gene OGA
can not found tss for gene OSBPL10-AS1
can not found tss for gene P3H2-AS1
can not found tss for gene PARPBP
can not found tss for gene PAXX
can not found tss for gene PCCA-AS1
can not found tss for gene PINX1
can not found tss for gene PIP4K2B
can not found tss for gene PIP4P1
can not found tss for gene PIP4P2
can not found tss for gene PJVK
can not found tss for gene PKIA-AS1
can not found tss for gene PLCE1-AS1
can not found tss for gene PLPBP
can not found tss for gene PPP1R26-AS1
can not found tss for gene PPP4R3A
can not found tss for gene PRC1-AS1
can not found tss for gene PRKCA-AS1
can not found tss for gene PRKN
can not found tss for gene PROSER2-AS1
can not found tss for gene PSMD6-AS2
can not found tss for gene PTPA
can not found tss for gene PXN-AS1
can not found tss for gene RAB11B-AS1
can not found tss for gene RAB5IF
can not found tss for gene RAP2C-AS1
can not found tss for gene RBMS3-AS2
can not found tss for gene RBMS3-AS3
can not found tss for gene RCC1L
can not found tss for gene RETREG1
can not found tss for gene RETREG3
can not found tss for gene REXO5
can not found tss for gene RFLNA
can not found tss for gene RIPOR3
can not found tss for gene RMC1
can not found tss for gene RMDN2-AS1
can not found tss for gene RNASEH2B-AS1
can not found tss for gene RORA-AS1
can not found tss for gene RTF2
can not found tss for gene RTL8A
can not found tss for gene RUBCNL
can not found tss for gene RUNDC3A-AS1
can not found tss for gene SAP30L-AS1
can not found tss for gene SAXO2
can not found tss for gene SBF2-AS1
can not found tss for gene SCOC-AS1
can not found tss for gene SCUBE3
can not found tss for gene SDCBP2-AS1
can not found tss for gene SEC24B-AS1
can not found tss for gene SELENOI
can not found tss for gene SELENOM
can not found tss for gene SELENON
can not found tss for gene SELENOP
can not found tss for gene SELENOT
can not found tss for gene SEMA6A-AS1
can not found tss for gene SEMA6A-AS2
can not found tss for gene SGMS1-AS1
can not found tss for gene SH3TC2-DT
can not found tss for gene SHLD1
can not found tss for gene SIAH2-AS1
can not found tss for gene SINHCAF
can not found tss for gene SLC12A9-AS1
can not found tss for gene SLC16A1-AS1
can not found tss for gene SLC16A12-AS1
can not found tss for gene SLC8A1-AS1
can not found tss for gene SLIRP
can not found tss for gene SNHG19
can not found tss for gene SNHG22
can not found tss for gene SNHG25
can not found tss for gene SNHG27
can not found tss for gene SPART
can not found tss for gene SPINDOC
can not found tss for gene SPON1-AS1
can not found tss for gene SPOUT1
chr diffs: chrX  chrY SPRY3
can not found tss for gene SPTY2D1OS
can not found tss for gene ST7-AS2
can not found tss for gene STARD4-AS1
can not found tss for gene STMP1
can not found tss for gene STN1
can not found tss for gene STX18-AS1
can not found tss for gene STXBP5-AS1
can not found tss for gene SYNPR-AS1
can not found tss for gene TAPT1-AS1
can not found tss for gene TENT2
can not found tss for gene TENT4A
can not found tss for gene TENT4B
can not found tss for gene TENT5A
can not found tss for gene THSD4-AS1
can not found tss for gene THUMPD3-AS1
can not found tss for gene TIMM29
can not found tss for gene TMC3-AS1
can not found tss for gene TMEM131L
can not found tss for gene TMEM202-AS1
can not found tss for gene TMEM268
can not found tss for gene TMEM273
can not found tss for gene TMEM44-AS1
can not found tss for gene TMEM72-AS1
can not found tss for gene TMEM94
can not found tss for gene TRAF3IP2-AS1
can not found tss for gene TRAM2-AS1
can not found tss for gene TRIR
can not found tss for gene TRMT9B
can not found tss for gene TSPEAR-AS1
can not found tss for gene TSPOAP1-AS1
can not found tss for gene TTI2
can not found tss for gene TTN-AS1
can not found tss for gene TTTY14
can not found tss for gene TUT4
can not found tss for gene TUT7
can not found tss for gene UBE2R2-AS1
can not found tss for gene UFD1
can not found tss for gene VASH1-AS1
can not found tss for gene VLDLR-AS1
can not found tss for gene VPS26C
can not found tss for gene VPS35L
can not found tss for gene WHRN
can not found tss for gene WWC3-AS1
can not found tss for gene XACT
can not found tss for gene XXYLT1-AS2
can not found tss for gene ZBED5-AS1
can not found tss for gene ZNF341-AS1
can not found tss for gene ZNF350-AS1
can not found tss for gene ZNF451-AS1
can not found tss for gene ZNF503-AS1
can not found tss for gene ZNF723
can not found tss for gene ZNF724
can not found tss for gene ZNRF3-AS1
can not found tss for gene ZSCAN16-AS1
can not found tss for gene ZUP1



# ##output only p2gGenes
# can not found tss for gene RRS1-AS1
# can not found tss for gene LINC00967
# can not found tss for gene LINC01119
# can not found tss for gene LINC01483
# can not found tss for gene RAP2C-AS1
# can not found tss for gene LINC01273
# can not found tss for gene MIR646HG
# can not found tss for gene DSG2-AS1
# can not found tss for gene LINC00589
# can not found tss for gene LINC01618
# can not found tss for gene MME-AS1
# can not found tss for gene LINC01411
# can not found tss for gene LINC01484
# can not found tss for gene DLGAP1-AS2
# can not found tss for gene LINC00882
# can not found tss for gene LINC01003
# can not found tss for gene DIRC3-AS1
# can not found tss for gene MRPL23-AS1
# can not found tss for gene MYLK-AS1
# can not found tss for gene SLC16A1-AS1
# can not found tss for gene AC087762.1
# can not found tss for gene AFDN
# can not found tss for gene EPHA1-AS1
# can not found tss for gene GAS1RR
# can not found tss for gene PINX1
# can not found tss for gene ABHD17C
# can not found tss for gene ARHGEF3-AS1
# can not found tss for gene DNMBP-AS1
# can not found tss for gene LINC01485
# can not found tss for gene MISP3
# can not found tss for gene BBOX1-AS1
# can not found tss for gene CARMIL1
# can not found tss for gene CNEP1R1
# can not found tss for gene COLGALT1
# can not found tss for gene GRK3
# can not found tss for gene HMGN3-AS1
# can not found tss for gene LINC00500
# can not found tss for gene LINC01210
# can not found tss for gene LINC01501
# can not found tss for gene LINC01569
# can not found tss for gene STX18-AS1
# can not found tss for gene FAM13A-AS1
# can not found tss for gene GVQW2
# can not found tss for gene LCMT1-AS2
# can not found tss for gene LINC00456
# can not found tss for gene N4BP3
# can not found tss for gene CD81-AS1
# can not found tss for gene LINC01301
# can not found tss for gene LINC01356
# can not found tss for gene NALCN-AS1
# can not found tss for gene RBMS3-AS3
# can not found tss for gene XXYLT1-AS2
# can not found tss for gene C8orf37-AS1
# can not found tss for gene HDAC11-AS1
# can not found tss for gene ISM1
# can not found tss for gene LARGE2
# can not found tss for gene LIN28B-AS1
# can not found tss for gene LINC01137
# can not found tss for gene LINC01509
# can not found tss for gene LUCAT1
# can not found tss for gene LY86-AS1
# can not found tss for gene MYCNUT
# can not found tss for gene PSMD6-AS2
# can not found tss for gene ST7-AS2
# can not found tss for gene TTTY14


# ##output full p2g gene
# can not found tss for gene A2ML1-AS1
# can not found tss for gene ABHD15-AS1
# can not found tss for gene ABHD17C
# can not found tss for gene AC006449.2
# can not found tss for gene AC008074.1
# can not found tss for gene AC008522.1
# can not found tss for gene AC009060.1
# can not found tss for gene AC012363.1
# can not found tss for gene AC018512.1
# can not found tss for gene AC087762.1
# can not found tss for gene AC138028.1
# can not found tss for gene ADM5
# can not found tss for gene AFDN
# can not found tss for gene AL358113.1
# can not found tss for gene AL591806.1
# can not found tss for gene ARHGAP26-AS1
# can not found tss for gene ARHGAP45
# can not found tss for gene ARHGEF3-AS1
# can not found tss for gene ASB16-AS1
# can not found tss for gene ATP6V0E2-AS1
# can not found tss for gene AZIN1-AS1
# can not found tss for gene BAALC-AS1
# can not found tss for gene BBOX1-AS1
# can not found tss for gene BDNF-AS
# can not found tss for gene BFSP2-AS1
# can not found tss for gene BICDL1
# can not found tss for gene BMT2
# can not found tss for gene BRWD1-AS1
# can not found tss for gene C15orf53
# can not found tss for gene C19orf48
# can not found tss for gene C21orf91-OT1
# can not found tss for gene C2orf83
# can not found tss for gene C6orf99
# can not found tss for gene C8orf37-AS1
# can not found tss for gene C9orf41-AS1
# can not found tss for gene CARMIL1
# can not found tss for gene CARMN
# can not found tss for gene CASC15
# can not found tss for gene CASC9
# can not found tss for gene CATIP-AS2
# can not found tss for gene CCDC144NL-AS1
# can not found tss for gene CCDC163
# can not found tss for gene CCDC183-AS1
# can not found tss for gene CD81-AS1
# can not found tss for gene CEP128
# can not found tss for gene CEP83
# can not found tss for gene CFLAR-AS1
# can not found tss for gene CHODL-AS1
# can not found tss for gene CLDN10-AS1
# can not found tss for gene CLRN1-AS1
# can not found tss for gene CNEP1R1
# can not found tss for gene CNOT9
# can not found tss for gene COA4
# can not found tss for gene COLGALT1
# can not found tss for gene COQ8B
# can not found tss for gene CPB2-AS1
# can not found tss for gene CRACR2A
# can not found tss for gene CRTC3-AS1
# can not found tss for gene CYP1B1-AS1
# can not found tss for gene CYP4A22-AS1
# can not found tss for gene CYYR1-AS1
# can not found tss for gene DCAF1
# can not found tss for gene DDX11-AS1
# can not found tss for gene DENND5B-AS1
# can not found tss for gene DHFR2
# can not found tss for gene DIAPH2-AS1
# can not found tss for gene DICER1-AS1
# can not found tss for gene DIO2-AS1
# can not found tss for gene DIRC3-AS1
# can not found tss for gene DLEU7-AS1
# can not found tss for gene DLGAP1-AS1
# can not found tss for gene DLGAP1-AS2
# can not found tss for gene DLGAP4-AS1
# can not found tss for gene DNMBP-AS1
# can not found tss for gene DSG2-AS1
# can not found tss for gene EBLN3P
# can not found tss for gene EGLN3-AS1
# can not found tss for gene ELMO1-AS1
# can not found tss for gene EMC6
# can not found tss for gene ENTPD1-AS1
# can not found tss for gene EP300-AS1
# can not found tss for gene EPHA1-AS1
# can not found tss for gene ERVH48-1
# can not found tss for gene ERVMER61-1
# can not found tss for gene FALEC
# can not found tss for gene FAM13A-AS1
# can not found tss for gene FAM214A
# can not found tss for gene FAM95C
# can not found tss for gene FER1L6-AS2
# can not found tss for gene FLG-AS1
# can not found tss for gene FRMD6-AS2
# can not found tss for gene FSIP2-AS1
# can not found tss for gene GAS1RR
# can not found tss for gene GATA3-AS1
# can not found tss for gene GON7
# can not found tss for gene GPAT4
# can not found tss for gene GRK2
# can not found tss for gene GRK3
# can not found tss for gene GRTP1-AS1
# can not found tss for gene GVQW2
# can not found tss for gene H1FX-AS1
# can not found tss for gene HDAC11-AS1
# can not found tss for gene HIF1A-AS2
# can not found tss for gene HIKESHI
# can not found tss for gene HM13-AS1
# can not found tss for gene HMGN3-AS1
# can not found tss for gene HORMAD2-AS1
# can not found tss for gene HPF1
# can not found tss for gene HYKK
# can not found tss for gene IFT43
# can not found tss for gene IGFBP7-AS1
# can not found tss for gene IL21-AS1
# can not found tss for gene INHBA-AS1
# can not found tss for gene IPO9-AS1
# can not found tss for gene IQCH-AS1
# can not found tss for gene ISM1
# can not found tss for gene JAZF1-AS1
# can not found tss for gene KANTR
# can not found tss for gene KAT14
# can not found tss for gene KCNMB2-AS1
# can not found tss for gene KCTD21-AS1
# can not found tss for gene KIF9-AS1
# can not found tss for gene LANCL1-AS1
# can not found tss for gene LARGE2
# can not found tss for gene LATS2-AS1
# can not found tss for gene LCMT1-AS2
# can not found tss for gene LIN28B-AS1
# can not found tss for gene LINC00159
# can not found tss for gene LINC00222
# can not found tss for gene LINC00240
# can not found tss for gene LINC00276
# can not found tss for gene LINC00323
# can not found tss for gene LINC00348
# can not found tss for gene LINC00370
# can not found tss for gene LINC00381
# can not found tss for gene LINC00456
# can not found tss for gene LINC00474
# can not found tss for gene LINC00500
# can not found tss for gene LINC00533
# can not found tss for gene LINC00571
# can not found tss for gene LINC00578
# can not found tss for gene LINC00589
# can not found tss for gene LINC00601
# can not found tss for gene LINC00607
# can not found tss for gene LINC00623
# can not found tss for gene LINC00630
# can not found tss for gene LINC00642
# can not found tss for gene LINC00656
# can not found tss for gene LINC00665
# can not found tss for gene LINC00667
# can not found tss for gene LINC00854
# can not found tss for gene LINC00862
# can not found tss for gene LINC00877
# can not found tss for gene LINC00879
# can not found tss for gene LINC00882
# can not found tss for gene LINC00884
# can not found tss for gene LINC00885
# can not found tss for gene LINC00909
# can not found tss for gene LINC00937
# can not found tss for gene LINC00958
# can not found tss for gene LINC00964
# can not found tss for gene LINC00967
# can not found tss for gene LINC00992
# can not found tss for gene LINC01003
# can not found tss for gene LINC01004
# can not found tss for gene LINC01023
# can not found tss for gene LINC01031
# can not found tss for gene LINC01088
# can not found tss for gene LINC01091
# can not found tss for gene LINC01095
# can not found tss for gene LINC01098
# can not found tss for gene LINC01099
# can not found tss for gene LINC01106
# can not found tss for gene LINC01119
# can not found tss for gene LINC01122
# can not found tss for gene LINC01134
# can not found tss for gene LINC01135
# can not found tss for gene LINC01137
# can not found tss for gene LINC01162
# can not found tss for gene LINC01184
# can not found tss for gene LINC01192
# can not found tss for gene LINC01194
# can not found tss for gene LINC01205
# can not found tss for gene LINC01210
# can not found tss for gene LINC01222
# can not found tss for gene LINC01234
# can not found tss for gene LINC01270
# can not found tss for gene LINC01273
# can not found tss for gene LINC01284
# can not found tss for gene LINC01291
# can not found tss for gene LINC01301
# can not found tss for gene LINC01307
# can not found tss for gene LINC01315
# can not found tss for gene LINC01317
# can not found tss for gene LINC01322
# can not found tss for gene LINC01336
# can not found tss for gene LINC01355
# can not found tss for gene LINC01356
# can not found tss for gene LINC01372
# can not found tss for gene LINC01411
# can not found tss for gene LINC01446
# can not found tss for gene LINC01456
# can not found tss for gene LINC01483
# can not found tss for gene LINC01484
# can not found tss for gene LINC01485
# can not found tss for gene LINC01495
# can not found tss for gene LINC01501
# can not found tss for gene LINC01503
# can not found tss for gene LINC01504
# can not found tss for gene LINC01505
# can not found tss for gene LINC01508
# can not found tss for gene LINC01509
# can not found tss for gene LINC01518
# can not found tss for gene LINC01527
# can not found tss for gene LINC01554
# can not found tss for gene LINC01569
# can not found tss for gene LINC01572
# can not found tss for gene LINC01588
# can not found tss for gene LINC01591
# chr diffs: chr12  chr20 LINC01598
# can not found tss for gene LINC01605
# can not found tss for gene LINC01607
# can not found tss for gene LINC01612
# can not found tss for gene LINC01618
# can not found tss for gene LRRC29
# can not found tss for gene LUCAT1
# can not found tss for gene LY86-AS1
# can not found tss for gene MAFTRR
# can not found tss for gene MAP3K14-AS1
# can not found tss for gene MCPH1-AS1
# can not found tss for gene MCRIP1
# can not found tss for gene MCUB
# can not found tss for gene MEF2C-AS1
# can not found tss for gene METTL26
# can not found tss for gene MIATNB
# can not found tss for gene MIR194-2HG
# can not found tss for gene MIR2052HG
# can not found tss for gene MIR222HG
# can not found tss for gene MIR3681HG
# can not found tss for gene MIR4300HG
# can not found tss for gene MIR4435-2HG
# can not found tss for gene MIR503HG
# can not found tss for gene MIR5689HG
# can not found tss for gene MIR646HG
# can not found tss for gene MIR9-3HG
# can not found tss for gene MISP3
# can not found tss for gene MME-AS1
# can not found tss for gene MMP25-AS1
# can not found tss for gene MPRIP-AS1
# can not found tss for gene MRM2
# can not found tss for gene MRM3
# can not found tss for gene MRPL23-AS1
# can not found tss for gene MTERF2
# can not found tss for gene MYCNUT
# can not found tss for gene MYLK-AS1
# can not found tss for gene MYO16-AS1
# can not found tss for gene N4BP3
# can not found tss for gene NALCN-AS1
# can not found tss for gene NCBP2-AS2
# can not found tss for gene NDUFV2-AS1
# can not found tss for gene NEPRO
# can not found tss for gene NEURL1-AS1
# can not found tss for gene NNT-AS1
# can not found tss for gene NR2F2-AS1
# can not found tss for gene NUTM2A-AS1
# can not found tss for gene OSBPL10-AS1
# can not found tss for gene OSMR-AS1
# can not found tss for gene OTUD6B-AS1
# can not found tss for gene OTX2-AS1
# can not found tss for gene PANO1
# can not found tss for gene PAPPA-AS1
# can not found tss for gene PARPBP
# can not found tss for gene PAXBP1-AS1
# can not found tss for gene PAXIP1-AS2
# can not found tss for gene PCCA-AS1
# can not found tss for gene PINK1-AS
# can not found tss for gene PINX1
# can not found tss for gene PIP4K2B
# can not found tss for gene PITRM1-AS1
# can not found tss for gene PKIA-AS1
# can not found tss for gene PLBD1-AS1
# can not found tss for gene PLCE1-AS1
# can not found tss for gene PLCG1-AS1
# can not found tss for gene PRC1-AS1
# can not found tss for gene PRKCA-AS1
# can not found tss for gene PROSER2-AS1
# can not found tss for gene PSMD6-AS2
# can not found tss for gene PTCHD1-AS
# can not found tss for gene RAB11B-AS1
# can not found tss for gene RAB30-AS1
# can not found tss for gene RAP2C-AS1
# can not found tss for gene RASAL2-AS1
# can not found tss for gene RBMS3-AS2
# can not found tss for gene RBMS3-AS3
# can not found tss for gene RCC1L
# can not found tss for gene RDH10-AS1
# can not found tss for gene RFX3-AS1
# can not found tss for gene RNASEH1-AS1
# can not found tss for gene RNF139-AS1
# can not found tss for gene RORA-AS1
# can not found tss for gene RRS1-AS1
# can not found tss for gene RUNDC3A-AS1
# can not found tss for gene SAMSN1-AS1
# can not found tss for gene SAP30L-AS1
# can not found tss for gene SBF2-AS1
# can not found tss for gene SCOC-AS1
# can not found tss for gene SCUBE3
# can not found tss for gene SEC24B-AS1
# can not found tss for gene SEMA3B-AS1
# can not found tss for gene SEMA3F-AS1
# can not found tss for gene SEMA6A-AS1
# can not found tss for gene SEPT4-AS1
# can not found tss for gene SEPT7-AS1
# can not found tss for gene SERPINB9P1
# can not found tss for gene SGMS1-AS1
# can not found tss for gene SH3BP5-AS1
# can not found tss for gene SLC16A1-AS1
# can not found tss for gene SLC16A12-AS1
# can not found tss for gene SLC2A1-AS1
# can not found tss for gene SLC8A1-AS1
# can not found tss for gene SLIRP
# can not found tss for gene SMC5-AS1
# can not found tss for gene SNHG19
# can not found tss for gene SNHG22
# can not found tss for gene SNHG25
# can not found tss for gene SPOUT1
# chr diffs: chrX  chrY SPRY3
# can not found tss for gene SRP14-AS1
# can not found tss for gene SRRM2-AS1
# can not found tss for gene ST7-AS2
# can not found tss for gene ST8SIA6-AS1
# can not found tss for gene STARD4-AS1
# can not found tss for gene STX18-AS1
# can not found tss for gene STXBP5-AS1
# can not found tss for gene SYNPR-AS1
# can not found tss for gene TAF1A-AS1
# can not found tss for gene TAPT1-AS1
# can not found tss for gene TENM4
# can not found tss for gene TESMIN
# can not found tss for gene TET2-AS1
# can not found tss for gene TFAP2A-AS1
# can not found tss for gene TH2LCRR
# can not found tss for gene THSD4-AS1
# can not found tss for gene TMC3-AS1
# can not found tss for gene TMCC1-AS1
# can not found tss for gene TMEM147-AS1
# can not found tss for gene TMEM161B-AS1
# can not found tss for gene TMEM72-AS1
# can not found tss for gene TMEM94
# can not found tss for gene TRAF3IP2-AS1
# can not found tss for gene TRAM2-AS1
# can not found tss for gene TRERNA1
# can not found tss for gene TSPEAR-AS1
# can not found tss for gene TSPOAP1-AS1
# can not found tss for gene TTC39A-AS1
# can not found tss for gene TTI2
# can not found tss for gene TTN-AS1
# can not found tss for gene TTTY14
# can not found tss for gene UBA6-AS1
# can not found tss for gene UBE2R2-AS1
# can not found tss for gene UBL7-AS1
# can not found tss for gene UCHL1-AS1
# can not found tss for gene UGDH-AS1
# can not found tss for gene URB1-AS1
# can not found tss for gene USP46-AS1
# chr diffs: chrX  chrY VAMP7
# can not found tss for gene VLDLR-AS1
# can not found tss for gene WAC-AS1
# can not found tss for gene WDCP
# can not found tss for gene WEE2-AS1
# can not found tss for gene XACT
# can not found tss for gene XXYLT1-AS2
# can not found tss for gene YTHDF3-AS1
# can not found tss for gene ZBED5-AS1
# can not found tss for gene ZBTB11-AS1
# can not found tss for gene ZIM2-AS1
# can not found tss for gene ZMYM4-AS1
# can not found tss for gene ZNF337-AS1
# can not found tss for gene ZNF341-AS1
# can not found tss for gene ZNF503-AS1
# can not found tss for gene ZNF529-AS1
# can not found tss for gene ZNF630-AS1
# can not found tss for gene ZRANB2-AS1
# can not found tss for gene ZSCAN16-AS1












