#for each atac cell, find the closest rna cell in the common liger space, cluster by cluster
#output a cell-cell match by cluster, with equal cell number 

##modified: pair within sample, use get.knnx(.. k = 3) to reduce duplication (unlike granja's simple method nearest distance)



#library('RANN') #2.6.1, nn2 only
library(FNN) #1.1.3, get.knn, get.knnx, knn.index, knnx.index

#in fact , RANN and FNN use the same algrithm (kd tree) to search nearset neighbor, using the C++ ANN library http://www.cs.umd.edu/~mount/ANN/

library(igraph)

library(ggplot2)
library(Rcpp) #1.0.8.2



color_good <- c("#E7D654", "#6F1482" ,"#DC7035", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", "#63AC4E", "#D181B0" ,
                "#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,"#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,
                "#CA362E" ,"#2B3918","#1E1E1E" )

sample <- 'placenta 8 weeks RNA integrate with ATAC'


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



sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX) > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY) > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }

    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);
    
    const int nx = X.ncol();
    const int ny = Y.ncol();

    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }

    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }

    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }

    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }

    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);
    } 

    return(cor);

  }'
)

scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}


# ###a test for RANN nn2, find the clostest point for query points####
# #https://stackoverflow.com/questions/24694707/get-closest-point-based-on-coordinates

# # your data
# df <- data.frame(X=c(1,2,2,2,3,4),Y=c(1,2,2,3,3,4))
# pts <- data.frame(X=c(2.1, 2.5), Y=c(2.3, 2.5))
# #library(RANN)

# plot(df)
# points(pts,col='red')


# # for each point in pts, find the nearest neighbor from df
# closest <- RANN::nn2(data = df, query = pts, k = 1)
# # argument k sets the number of nearest neighbours, here 1 (the closest)
# closest
# # $nn.idx
# # [,1]
# # [1,]    3
# # [2,]    5
# # 
# # $nn.dists
# # [,1]
# # [1,] 0.3162278
# # [2,] 0.7071068
# # Get coordinates of nearest neighbor
# pts$X_snap <- df[closest$nn.idx, "X"]
# pts$Y_snap <- df[closest$nn.idx, "Y"]
# pts
# #     X   Y X_snap Y_snap
# # 1 2.1 2.3      2      2
# # 2 2.5 2.5      3      3



# #######################a little test for FNN###############
# set.seed(123)
# x1 = cbind(runif(10),runif(10)) ##reference xy coordinates, can extend to n dims?
# #x1.bk <- x1
# #row.names(x1) = paste0('x1_',1:nrow(x1)) #no use to add rownames
# len1 = nrow(x1)
# row.names(x1) <- paste0('ref',1:nrow(x1))

# x2 = cbind(runif(12),runif(12)) ##query xy coordinates, can extend to n dims?
# #x2.bk <- x2
# #x2 = cbind(runif(2),runif(2))
# len2 = nrow(x2)
# #row.names(x2) = paste0('x2_',1:nrow(x2))
# row.names(x2) <- paste0('que',1:nrow(x2))

# ###plot in xy coordinate to view dots nearest directly####
# par(mfrow=c(1,3))
# options(repr.plot.width=15,repr.plot.height=5)
# plot(x1,cex=5,pch  = 19,col='black',main='reference',cex.main=3)
# text(x1,labels = 1:len1,col='grey',cex = 2 )
# plot(x2,cex=5,pch  = 19,col='grey',main='query',cex.main=3)
# text(x2,labels = 1:len2,col='red',cex = 2 )
# xall <- rbind(x1,x2)
# lenall <- nrow(xall)

# #options(repr.plot.width=10,repr.plot.height=8)
# plot(xall,cex=5,pch  = 19,col=c(rep('black',len1),rep('grey',len2)),main='merge',cex.main=3 )#  c('black','grey'),each=10 ) )
# text(xall,labels = c(1:len1,1:len2),col=c(rep('grey',len1),rep('red',len2)),cex = 2 )#rep(c('grey','red'),each=10 ) )


# #######search and plot with igraph########
# # k = 5
# # ###plot x1 graph and x2 graph
# # ##graph1 (reference)
# # par(mfrow=c(1,1))
# # options(repr.plot.width=5,repr.plot.height=5)
# # nn1 = get.knn(x1,k)
# # test.df1 = data.frame(from = rep(1:nrow(nn1$nn.index), k), 
# #                     to = as.vector(nn1$nn.index), 
# #                     weight = 1/(1 + as.vector(nn1$nn.dist))
# #                    )
# # test.df.nw1 = igraph::graph_from_data_frame(test.df1, directed = TRUE)
# # plot(test.df.nw1,main='graph of reference',cex.main=3)


# # ##graph2 (query)
# # par(mfrow=c(1,1))
# # options(repr.plot.width=5,repr.plot.height=5)
# # nn2 = get.knn(x2,k)
# # test.df2 = data.frame(from = rep(1:nrow(nn2$nn.index), k), 
# #                     to = as.vector(nn2$nn.index), 
# #                     weight = 1/(1 + as.vector(nn2$nn.dist))
# #                    )
# # test.df.nw2 = igraph::graph_from_data_frame(test.df2, directed = TRUE)
# # plot(test.df.nw2,main='graph of query',cex.main=3)


# ##query x2 from x1
# #set.seed(100) #graph varies even set seed, but topology is similar
# k = 3
# ##k = 1
# par(mfrow=c(1,1))
# options(repr.plot.width=6,repr.plot.height=6)
# nn = get.knnx(data = x1,query = x2,k = k) #return query length of reference point


# ##construct a pairwise result table, use the closest point

# #x1: ref
# #x2 query

# idx <- as.vector(nn$nn.index) #k = 1, one column mat
# dist <- as.vector(nn$nn.dist)

# if(length(idx) == nrow(x2) ){ #hit_in_ref length must eq query
#  result_pair <- data.frame(query_x = x2[,1], query_y = x2[,2], hit_in_ref_x = x1[idx, 1], hit_in_ref_y = x1[idx, 2]  )
#  result_pair_id <- data.frame(query = row.names(x2), hits_in_ref = row.names(x1)[idx] )
    
# }else{cat('can not create result pair table')}


# test.df = data.frame(from = paste0('q',rep(1:nrow(nn$nn.index), k)), 
#                     to = as.vector(nn$nn.index),
#                     #weight =  as.vector(nn$nn.dist)
#                     weight = 1/(1 + as.vector(nn$nn.dist))
#                    )
# test.df.nw = igraph::graph_from_data_frame(test.df, directed = TRUE)
# plot(test.df.nw,
#      main='graph of qeury to reference',
#      cex.main=3
#      )
#      #vertex.size=15,
#      #edge.width=E(test.df.nw)$weight)

# ##modify the network attributes
# get.data.frame(test.df.nw,what='vertices')

# #see vertex order
# V(test.df.nw)

# #set vertex color
# V(test.df.nw)$color <- c(rep('red',12),rep('orange',10))

# #plot the new graph
# plot(test.df.nw,
#      main='graph of qeury to reference',
#      cex.main=3
#      )

# # length(edge.betweenness(test.df.nw) )
# # length(E(test.df.nw)$weight )

# # #edge.betweenness(test.df.nw) = E(test.df.nw)$weight * 5 #not settable
# # edge_attr(test.df.nw) #see edge attribes

# # test.df.nw <- test.df.nw %>% set_edge_attr('color',value='red') 
# # plot(test.df.nw,main='graph of qeury to reference',cex.main=3)



# #get all nearest from reference
# unique(as.vector(nn$nn.index) )

# #filter eucliean dist?
# boxplot(as.vector(nn$nn.dist))
# filter <- quantile(as.vector(nn$nn.dist),probs = c(0.1,0.25,0.5,0.75,0.9,1) )[3]
# unique(as.vector(nn$nn.index) [as.vector(nn$nn.dist <= filter)])

# #####test end########






##########method 1, use simple smallest euclidean distance, Granja method#########
#Nearest Neighbor differential
findNN <- function(query, reference, method = "euclidean"){
    findClosest <- function(x, m, method = "euclidean"){
        if(method=="euclidean"){
            which.min(sqrt(colSums((t(m) - x) * (t(m) - x))))#one to all minimal euclidean distance
        }else if(method=="pearson"){
            which.max(cor(t(m),x,method = method)[,1])
        }else if(method=="spearman"){
            which.max(cor(t(m),x,method = method)[,1])
        }
    }

    pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
    mat <- data.frame(matrix(ncol = 4, nrow = nrow(query)))
    colnames(mat) <- c("x", "i", "y", "j")
    for(i in seq_len(nrow(query))){
        setTxtProgressBar(pb,round(i*100/nrow(query),0))
      j <- findClosest(query[i,], reference, method)
      mat[i,] <- c(x = rownames(query)[i], i = i, y = rownames(reference)[j], j = j)
    }
    return(mat)
}



    
findNN_knnx <- function(query, reference, k = 5){ 
    
    #search nearest dot with FNN::getknnx(query,reference, k = 5)
    #iteratively choose non-redundant paired id within k = 5
    
    #pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
    mat <- data.frame(matrix(ncol = 4, nrow = nrow(query)))
    colnames(mat) <- c("x", "i", "y", "j")
    
    match_df = FNN::get.knnx(reference, query, k = k)$nn.index #query:atac, reference:rna; return query row
    
    stopifnot(nrow(match_df) == nrow(query) )
    stopifnot(max(match_df) <= nrow(reference))
    
    rownames(match_df) <- rownames(query)
    
    mat$x <- rownames(query)
    mat$i <- 1:nrow(query)
    
    hits <- match_df[,1]
    
    for(i in 2:k){
        
        idx <- which(duplicated(hits))
        hits[idx] <- match_df[idx,i]
        
        
    }
    
    cat('before reduce redundancy: ',sum(duplicated(match_df[,1])), '\n' ,sep='') #1781 of 3594 duplicated
    cat('after reduce redundancy: ',sum(duplicated(hits)), '\n' ,sep='')   #925 left 
    
    mat$j <- hits
    mat$y <- rownames(reference)[hits]
    
    return(mat)

}






a.placenta.te <- readRDS('../a.placenta.te.rds')

cluster.df.add.te <- readRDS('../cluster.df.add.te.rds')

all.equal(rownames(cluster.df.add.te),rownames(a.placenta.te@H.norm) ) #TRUE


###intersection gene id to generate merged data slot###
# Create Seurat object from liger object, keeping liger highly variable genes
##seurat_obj = ligerToSeurat(a.placenta.te, use.liger.genes = T) #do not use this!!, merged row.data slot changed at certain genes!
##saveRDS(seurat_obj,'../seurat_obj_fromliger.te.rds')

seurat_obj <- readRDS('../seurat_obj_fromliger.te.rds')
#DefaultAssay(seurat_obj) <- "RNA"


#plot(seurat_obj@reductions$tsne@cell.embeddings,pch=19,cex=0.1)
library(Seurat)
Seurat::DimPlot(seurat_obj, reduction = "tsne", #UMAP in fact
              #group.by = "Eday", 
              #cols = color_set3_ext17,
              label.size = 10,
              label = TRUE,
              #repel = TRUE, 
              pt.size = 0.1) +

              guides(colour = guide_legend(override.aes = list(size=8)))
        

mat.combine <- seurat_obj@assays$RNA@data #62681 x 21268  #renormed combined (with gene combine, not intersect) 
#mat.combine <- seurat_obj@assays$RNA@counts #62681 x 21268

table(seurat_obj@meta.data$orig.ident)    
atac1 atac2  rna1  rna2 
 4549  6525  5978  4216 

table(a.placenta.te@cell.data$dataset)
 atac1 atac2  rna1  rna2 
 4549  6525  5978  4216    


# #########################do not use liger2seurat !!, merged row.data slot changed at certain genes!
# genes.diff <- readRDS('../peak_gene_link/genes.diff.rds')    


# counts <- seurat_obj@assays$RNA@counts 
# cellid <- colnames(counts)
# cellid <- gsub(pattern = '^atac1_|^atac2_|^rna1_|^rna2_',replacement = '',x=cellid)
# colnames(counts) <- cellid 

    
# shared_cellid <- intersect(colnames(counts),c(colnames(a.placenta.te@raw.data[['atac1']]),colnames(a.placenta.te@raw.data[['atac2']]) ) )
# temp <- cbind(a.placenta.te@raw.data[['atac1']][genes.diff,],a.placenta.te@raw.data[['atac2']][genes.diff,])

# all.equal(temp[genes.diff,shared_cellid],counts[genes.diff,shared_cellid]) #FALSE

# matGS <- readRDS('../../seuratv3CCA/mat.atac.modified.rds')

# genes.diff <- genes.diff[genes.diff %in% rownames(matGS)]
# shared_cellid <- intersect(colnames(temp),colnames(matGS) )
# all.equal(temp[genes.diff,shared_cellid],matGS[genes.diff,shared_cellid] ) #TRUE

# shared_cellid <- intersect(colnames(counts),colnames(matGS) )
# all.equal(counts[genes.diff,shared_cellid],matGS[genes.diff,shared_cellid] ) #FALSE
    

# shared_geneid <- intersect(rownames(counts), rownames(matGS) )
# shared_geneid <- shared_geneid[!shared_geneid %in% genes.diff ]

# all.equal(counts[shared_geneid,shared_cellid],matGS[shared_geneid,shared_cellid] )


# # #############liger2seurat not use#########################

# allGenes <- unique(unlist(lapply(a.placenta.te@raw.data, rownames))) #62681
# shared_geneid <- Reduce(intersect,lapply(a.placenta.te@raw.data, rownames) ) #17201 (17512 in seurat cca)

# all.equal(allGenes, rownames(mat.combine)) #193 diff, name change with underscore

# allGenes[which(allGenes != rownames(mat.combine) ) ]#193 diff
# rownames(mat.combine) [which(allGenes != rownames(mat.combine) ) ]

# mat.shared <- mat.combine[shared_geneid,] #17201 x 21268 (still use intersect mat)

# #shared_rowid_CCA <- readRDS("../../greenleaf_lsi_RNA_ATAC_seuratV2CCA/results/shared_rowid.rds")
# #share_share <- intersect(shared_geneid,shared_rowid_CCA) #16868

# cellid <- colnames(mat.shared)
# cellid <- gsub(pattern = '^atac1_|^atac2_|^rna1_|^rna2_',replacement = '',x=cellid)
# colnames(mat.shared) <- cellid 
# mat.data <- mat.shared

# matGS <- readRDS('../../seuratv3CCA/mat.atac.modified.rds')
# shared_geneid <- intersect(rownames(mat.data ), rownames(matGS) ) #16936
# shared_cellid <- intersect(colnames(mat.data),colnames(matGS) )

# all.equal(mat.data[shared_geneid,shared_cellid],matGS[shared_geneid,shared_cellid] )


# saveRDS(mat.shared,'matching_group_te/mat.shared.rds')    
# #saveRDS(mat.shared,'mat.shared.raw.rds') 



####use 20 dims (use this?)
matINMF <- a.placenta.te@H.norm #43207 x 35 #21268 x 20
cellid_atac <- grep(pattern = "atac",x = rownames(matINMF) , value = TRUE) #22013 #11074
cellid_rna <- grep(pattern = "rna",x = rownames(matINMF) , value = TRUE) #21194 #10194
intersect(cellid_atac,cellid_rna) #0

#group by sample

cellid_atac.df <- data.frame(cellid_atac = cellid_atac, sample= paste0('D',sapply(stringr::str_split(cellid_atac,'-',n=2) , function(x){x[2]}  )),stringsAsFactors = FALSE )

cellid_atac.list <-  split(cellid_atac.df, f = cellid_atac.df$sample)

lapply(cellid_atac.list,function(x){nrow(x)} )
D1 3594
D2 5105
D3 3842
D4 3513
D5 1908
D6 4051

cellid_rna.df <- data.frame(cellid_rna = cellid_rna, sample= paste0('D',sapply(stringr::str_split(cellid_rna,'-',n=2) , function(x){x[2]}  )),stringsAsFactors = FALSE )

cellid_rna.list <-  split(cellid_rna.df, f = cellid_rna.df$sample)

lapply(cellid_rna.list,function(x){nrow(x)} )
D1 4375
D2 3197
D3 3812
D4 3063
D5 2589 
D6 4158


table(cluster.df.add.te$anno)
rna_D1 atac_D1  rna_D2 atac_D2  rna_D3 atac_D3  rna_D4 atac_D4  rna_D5 atac_D5 
   4375    3594    3197    5105    3812    3842    3063    3513    2589    1908 
 rna_D6 atac_D6 
   4158    4051 




saveRDS(matINMF,'matching_group_te/matINMF.rds')    
saveRDS(cellid_atac.list,'matching_group_te/cellid_atac.list.rds') 
saveRDS(cellid_rna.list,'matching_group_te/cellid_rna.list.rds') 

# ##use UMAP 2-dim space
# matUMAP  <- as.matrix(cluster.df.add.te[,c('UMAP_1','UMAP_2')] )
# cellid_atac <- grep(pattern = "atac",x = rownames(matUMAP) , value = TRUE) #11074
# cellid_rna <- grep(pattern = "rna",x = rownames(matUMAP) , value = TRUE) #10194
# intersect(cellid_atac,cellid_rna) #0

# all.equal(a.placenta.te@tsne.coords,matUMAP,check.attributes = FALSE ) #TRUE

# saveRDS(matUMAP,'matching_group_te/matUMAP.rds')


#KNN Search
#Alternatively for speed FNN::get.knnx(query, reference, k = 1)
#We just used a simple function
matchedCells <- findNN(
  query = matINMF[cellid_atac,], #22013 x 35 #11074 x 20
  reference = matINMF[cellid_rna,], #21194 x 35 #10194 x 20
  #query = matUMAP[cellid_atac,], #11074 x 2
  #reference = matUMAP[cellid_rna,], #10194 x 2
  method = "euclidean")

##use get.knnx() with k , choose non-redundancy paired id iteratively
matchedCells.bk <- matchedCells    
matchedCells.list <- lapply( names(cellid_atac.list), function(x){
    cellid_atac_in <- cellid_atac.list[[x]]$cellid_atac
    cellid_rna_in <- cellid_rna.list[[x]]$cellid_rna
#     res <- findNN(
#       query = matINMF[cellid_atac_in,], 
#       reference = matINMF[cellid_rna_in,],
#       method = "euclidean"
#     )
    
    res <- findNN_knnx(
      query = matINMF[cellid_atac_in,], 
      reference = matINMF[cellid_rna_in,],
      k = 3
    )
        
    return(res)
})

names(matchedCells.list) <- names(cellid_atac.list)

lapply(matchedCells.list,function(x){ sum(duplicated(x$y))  } )

$D1 1166 $D2 2777 $D3 1589 $D4 1755 $D5 693 $D6 1843 #k = 3
$D1 925 $D2 2566 $D3 1336 $D4 1569 $D5 581 $D6 1625 #k = 5

$D1 1781 #graja's simple method
$D2 3308 
$D3 2131 
$D4 2183 
$D5 983 
$D6 2394

matchedCells <- Reduce(rbind.data.frame,matchedCells.list)

#compare with original method (granja's simple smallest euclidean distance method)
table(matchedCells$x == matchedCells.bk$x) #TRUE
table(matchedCells$y == matchedCells.bk$y)
FALSE  TRUE 
17695  4318 

length(intersect(unique(matchedCells$y),unique(matchedCells.bk$y))) #7507 #6344
length(unique(matchedCells$y)) #12190 #9233
length(unique(matchedCells.bk$y)) #9360

##RESULT: match atac with rna, all-to-all largely similar with match atac with rna, sample-to-sample; both have many duplications;

sum(duplicated(matchedCells[,'y'])) #9823 #8602 of 22013 dup #12653 of 22013 dup #5892 of 11074 #5395 dup #5472 of 10023    
length(unique(matchedCells[,'y']) ) #12190 #13411 uniq #9360 uniq #5182 unique #5679 unique  #4551 unique

all.equal(matchedCells$x,cellid_atac ) #TRUE


#filter duplicated in rna cellid (no need to do this)
##matchedCells.filter <- matchedCells[!duplicated(matchedCells[,'y']),] #5182

##readin merged mat
##mat.data <- readRDS('data/liger.mat.rds') #17201 x 21268
##mat.data <- readRDS('data/liger.norm.mat.rds') #17203 x 21268
mat.data <- readRDS('data/liger.norm.toseurat.mat.rds')#15261 x 43207  #17201 x 21268, use this
#liger ga + rna raw.data to seurat then normalize by seurat

##mat.data <- readRDS('data/cca.norm.mat.rds') #17512 x 20221
#mat.data.1 <- readRDS('../../greenleaf_lsi_RNA_ATAC_seuratV2CCA/results/CCA.data.rds') #20221, the same with above
    
all.equal(rownames(matINMF),colnames(mat.data) )   #TRUE
table(matchedCells$x %in% colnames(mat.data) ) #TRUE
table(matchedCells$y %in% colnames(mat.data) ) #TRUE

    
# shared_cellid <- intersect(rownames(matINMF),colnames(mat.data))    
    
# matINMF <- matINMF[shared_cellid,]
# mat.data <- mat.data[,shared_cellid]
# all.equal(rownames(matINMF),colnames(mat.data) ) 

# matchedCells <- matchedCells[matchedCells$x %in%  colnames(mat.data),]

#readin var.gene
var.genes <- readRDS('data/var.genes.liger.rds')  #2121 #1584 #use this
##var.genes <- readRDS('data/var.genes.liger.toseurat.rds') #do not use this
##var.genes <- readRDS('data/var.genes.cca.rds')

table(var.genes %in% rownames(mat.data))
    TRUE 
    2121 
    
FALSE  TRUE 
    2  1582

# FALSE  TRUE 
#   104  3908 

var.genes <- var.genes[var.genes %in% rownames(mat.data)]

saveRDS(var.genes,'matching_group_te/var.genes.rds')
saveRDS(mat.data,'matching_group_te/mat.data.rds')

###add cor column
matchedCells$corCCA <- rowCorCpp(
  match(matchedCells$x, colnames(mat.data)), 
  match(matchedCells$y, colnames(mat.data)),
  matINMF,matINMF)
  #alignedCCA, alignedCCA)

matchedCells$corVarRNA <- rowCorCpp( #correlated with column data (all gene data), seurat normalized tn5 ga vs rna
  match(matchedCells$x, colnames(mat.data)), 
  match(matchedCells$y, colnames(mat.data)), 
  t(as.matrix(mat.data[var.genes,])), 
  t(as.matrix(mat.data[var.genes,])))

#matchx <- match(matchedCells$x, colnames(mat.data))
#matchy <- match(matchedCells$y, colnames(mat.data))
#mat <- as.matrix(mat.data[var.genes,]) #1582 × 21268

saveRDS(matchedCells,'matching_group_te/matchedCells.rds') #22013 × 6 #11074 x 6
saveRDS(matchedCells,'matching_group_te/matchedCells.pair_sample.get_knnx_k5.rds') #22013 × 6
saveRDS(matchedCells,'matching_group_te/matchedCells.pair_sample.get_knnx_k3.rds') #22013 × 6

##duplication rate if not filter

y.count <- table(matchedCells$y)
options(repr.plot.width=4.5,repr.plot.height=4.5)
hist(y.count,breaks =20,xlim=c(1,20))

table (duplicated(matchedCells$y) )

FALSE  TRUE #k = 3
12190  9823 

FALSE  TRUE 
 9360 12653

FALSE  TRUE 
 5182  5892

#   FALSE  TRUE 
#  4551  5472

#plot distribution of correlation ratio within CCA@dr$cca.aligned and CCA@data    

table(matchedCells$corCCA >= 0.45)

FALSE  TRUE #k = 3
   13 22000 
    
FALSE  TRUE #k = 5
   18 21995 
    
FALSE  TRUE 
    1 22012
    
TRUE 
11074 
    
#     FALSE  TRUE 
#    28  9995
  
hist(matchedCells$corVarRNA)
abline(v=0.3)

table(matchedCells$corVarRNA >= 0.25) #?
FALSE  TRUE #k = 3
 4826 17187

FALSE  TRUE #k = 5
 4906 17107

table(matchedCells$corVarRNA >= 0.3)


FALSE  TRUE #k = 3
11532 10481

FALSE  TRUE # k = 5
11661 10352
    
    FALSE  TRUE 
12146  9867
    

FALSE  TRUE 
 1035 10039
    
# FALSE  TRUE 
#  5870  5204

# FALSE  TRUE 
#  1044 10030

# FALSE  TRUE 
#  1409  8614
    
table(matchedCells$corVarRNA >= 0.35)
    
FALSE  TRUE #k = 3
18108  3905

FALSE  TRUE #k =5
18233  3780
    
    
FALSE  TRUE 
 2408  8666
    
# FALSE  TRUE 
#  7149  3925 

# FALSE  TRUE 
#  2425  8649 


matchedCells.filter0.3 <- matchedCells[matchedCells$corVarRNA >= 0.3,] #10481 #10352 #9867 #10039 #8614
matchedCells.filter0.25 <- matchedCells[matchedCells$corVarRNA >= 0.25,] #17107

matchedCells.filter <- matchedCells.filter0.3 #10481, use k = 3, filter with corVarRNA >= 0.3

##add UMAP and cluster to matchedCell df
cluster.df.add.te <- readRDS("../cluster.df.add.te.rds") #21268
cluster.all.coord <- readRDS('../cluster.all.coord.rds') #23732

table (rownames(cluster.df.add.te) %in% rownames(cluster.all.coord) )
TRUE 
43207
    
# TRUE 
# 21268 

table(cluster.df.add.te$type)
atac   rna 
22013 21194
    
#  atac   rna 
# 11074 10194 

all.equal(rownames(cluster.df.add.te),rownames(matINMF)) #TRUE
#all.equal(rownames(cluster.df.add.te),rownames(matUMAP)) #TRUE
all.equal(matchedCells$x,rownames(subset( cluster.df.add.te, type == 'atac') )) #TRUE


cluster.df.add_pair <- data.frame(cellid_atac = matchedCells.filter$x,
    cluster.df.add.te[matchedCells.filter$x,c('cluster_lib','UMAP_1','UMAP_2','type','anno')],
                                  cellid_rna = matchedCells.filter$y,
    cluster.df.add.te[matchedCells.filter$y,c('cluster_lib','UMAP_1','UMAP_2','type','anno')],
                                    stringsAsFactors = FALSE
)

all.equal(rownames(cluster.df.add_pair),cluster.df.add_pair$cellid_atac,check.attributes = FALSE) #TRUE
row.names(cluster.df.add_pair) <- NULL
colnames(cluster.df.add_pair) <- c('cellid_atac','cluster_atac','UMAP_1_atac','UMAP_2_atac','type_atac','anno_atac','cellid_rna','cluster_rna','UMAP_1_rna','UMAP_2_rna','type_rna','anno_rna')
    
#cluster.df.add_pair.filter <- cluster.df.add_pair[!duplicated(cluster.df.add_pair[,'cellid_rna']),] #5182
cluster.df.add_pair.filter <- cluster.df.add_pair    
    
saveRDS(cluster.df.add_pair,'matching_group_te/cluster.df.add_pair.rds')
saveRDS(cluster.df.add_pair,'matching_group_te/cluster.df.add_pair.use_getknnx_k3.filterRNACor0.3.rds') #use this?
#saveRDS(cluster.df.add_pair.filter,'cluster.df.add_pair.filter.rds')

#cluster.df.add <- cluster.df.add_pair
#cluster.df.add <- cluster.df.add_pair.filter

table(cluster.df.add_pair.filter$cluster_atac)
     1    7    9    2    6    4    5    8    3   10   11 #0.3
1666  448  328 2005 1908 1226  911  625 1364    0    0
    
      1    7    9    2    6    4    5    8    3   10   11 #0.25
3769 1568  794 3483 3575 2423 2178 1380 2843    0    0
    
    
  1    2    3    4    5    6    7    8    9   10 
1827 1782 1086 1545 1474 1506 1017  443  394    0 

#   1   2   3   4   5   6   7   8   9  10 
# 842 874 467 745 714 692 427 244 177   0

#    1    2    3    4    5    6    7    8    9 
# 1105  853  682  638  641  587  532  126   18 

table(cluster.df.add_pair.filter$cluster_rna)
     1    7    9    2    6    4    5    8    3   10   11 #0.3
1288  296  949 1030  899 2179  802  750  570 1475  243
    
     1    7    9    2    6    4    5    8    3   10   11  #0.25
2083  659 2007 1730 2767 4245 1802 1529 1331 3288  572
    
1    2    3    4    5    6    7    8    9   10 
1836 1557 1438 1059 1475  976 1124  457  787  365 
    
#   1   2   3   4   5   6   7   8   9  10 
# 890 756 676 529 641 459 484 246 344 157 
    
#     1    2    3    4    5    6    7    8    9 
# 1173  805  740  669  479  617  531  125   43

dist <- vector()    

options(repr.plot.width=12,repr.plot.height=12)
plot(cluster.df.add_pair.filter[,c('UMAP_1_atac','UMAP_2_atac')],pch=19,cex = 0.1, col = 'red')
points(cluster.df.add_pair.filter[,c('UMAP_1_rna','UMAP_2_rna')],pch=19,cex = 0.1, col = 'blue')
#for(i in seq_len(nrow(cluster.df.add_pair.filter))){
for(i in sample(x = seq_len(nrow(cluster.df.add_pair.filter)),size = 500,replace = FALSE ) ){
  x0 <- cluster.df.add_pair.filter[i,'UMAP_1_atac']
  y0 <- cluster.df.add_pair.filter[i,'UMAP_2_atac']
  x1 <- cluster.df.add_pair.filter[i,'UMAP_1_rna']
  y1 <- cluster.df.add_pair.filter[i,'UMAP_2_rna']

  dist <- c(dist,sqrt((x1-x0)^2 + (y1-y0)^2) )
  segments(x0,y0,x1,y1,col = 'black',lwd = 2)
}

x0 <- cluster.df.add_pair.filter[,'UMAP_1_atac']
y0 <- cluster.df.add_pair.filter[,'UMAP_2_atac']
x1 <- cluster.df.add_pair.filter[,'UMAP_1_rna']
y1 <- cluster.df.add_pair.filter[,'UMAP_2_rna']

dist <- c(dist,sqrt((x1-x0)^2 + (y1-y0)^2) )

hist(dist,breaks = 100)
    
#how close are they ?    
    

##confusion matrix
cM <- as.matrix(ArchR::confusionMatrix(cluster.df.add_pair.filter$cluster_atac, cluster.df.add_pair.filter$cluster_rna))


preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

pheatmap::pheatmap(as.data.frame(cM))

options(repr.plot.width=5.5,repr.plot.height=5)
#pheatmap::pheatmap(as.data.frame(log(cM+1))
pheatmap::pheatmap(as.data.frame(cM)[c('1','7','9','2','6','8','3','4','5'),c('9','6','11','8','4','1','3','2','10')],cluster_rows = FALSE,cluster_cols = FALSE,color = color_tfdev,border_color = 'white',fontsize = 24,scale = 'column',breaks = seq(-0.5,2.5,3/length(color_tfdev)) )

    
atac:
1-7-9-2-6-8-3-4-5

rna:
9-6-11-8-4-1-3-2-10


    
    
####pairing atac and rna nuclei done


    

    
##have a quick  (slow to pheatmap dcgmatrix) look on gene activity and expression correlation degree with filtered matched cell pair

matRNA <- readRDS('../../seuratv3CCA/mat.rna.modified.rds') #24307 x 10198
colnames(matRNA) <- gsub(pattern = '__',replacement = "_", x = colnames(matRNA))


matGS <- readRDS('data/geneScoreMatrix.norm.rds')#37956 x 10041


matRNA <- readRDS('data/geneExpressionMatrix.norm.rds')

    
table(cluster.df.add_pair.filter$cellid_atac %in% colnames(matGS) )
    
FALSE  TRUE 
 1049 10025
    
# FALSE  TRUE 
#   481  4701

table(cluster.df.add_pair.filter$cellid_rna %in% colnames(matRNA) )    
TRUE 
11074
    
# TRUE 
# 5182  

cluster.df.add_pair.filter.filter <- cluster.df.add_pair.filter[cluster.df.add_pair.filter$cellid_atac %in% colnames(matGS) ,] #10025 #4701

#Gene Universe
geneUniverse <- intersect(rownames(matGS),rownames(matRNA)) #17512


    
#transform data
matGS <- edgeR::cpm(y = matGS,log = TRUE, prior.count = 1,)
matRNA <- edgeR::cpm(y = matRNA,log = TRUE, prior.count = 1)





gene.list <- var.genes
    
gene.test<- c("FOS", "JUNB", "BACH1", "SMARCC1", "TP63", "TOPORS")

gene.test <- c("FOS", "JUNB",'ATF3', "STAT4", "STAT5A", "STAT5B", "STAT3")

gene.test <- gene.list[1:500]    
#
#gene.test  <- gene.list   
table(gene.test %in% geneUniverse) #TRUE


# data_GS <- matGS[gene.test,cluster.df.add_pair.filter.filter$cellid_atac]
# dim(data_GS) #500 x 10025  #7 x 10025 #7 x 4701

# options(repr.plot.width=12,repr.plot.height=5)
# res.p_GS <- pheatmap::pheatmap(data_GS,show_colnames = FALSE)
# roworder <- res.p_GS$tree_row$order
# colorder <- res.p_GS$tree_col$order

# options(repr.plot.width=12,repr.plot.height=5)
# pheatmap::pheatmap(data_GS[roworder,colorder],show_colnames = FALSE,cluster_rows = FALSE, cluster_cols = FALSE,scale = 'column')


# data_RNA <- matRNA[gene.test,cluster.df.add_pair.filter.filter$cellid_rna]    
# dim(data_RNA)#500 x 10025  #7 x 10025#7 x 4701
    
# pheatmap::pheatmap(data_RNA[roworder,colorder],show_colnames = FALSE,cluster_rows = FALSE, cluster_cols = FALSE,scale = 'column' )




mat.data <- matGS[gene.test,] #500 x 10041
mat.data.z <- scale_rows(mat.data) # #zscore by row

table(cluster.df.add_pair.filter.filter$cellid_atac %in% colnames(mat.data.z)     )
 TRUE 
10025 

#plot with hclust and color breaks
options(repr.plot.width=8,repr.plot.height=10)

breaks <- seq(-2,2,by=4/(length(color_peak)-1))
res.p <- pheatmap::pheatmap(mat.data.z[,cluster.df.add_pair.filter.filter$cellid_atac],show_colnames = FALSE,show_rownames = TRUE,cluster_rows = TRUE, cluster_cols = TRUE,treeheight_col = 0,treeheight_row = 0,color = color_peak,breaks=breaks)
    
roworder <- res.p$tree_row$order
colorder <- res.p$tree_col$order


mat.data <- matRNA[gene.test,] #500 x 10198
mat.data.z <- scale_rows(mat.data) # #zscore by row

table(cluster.df.add_pair.filter.filter$cellid_rna[colorder] %in% colnames(matRNA))
    
breaks <- seq(-2,2,by=4/(length(color_tfdev)-1))
pheatmap::pheatmap(mat.data.z[roworder,cluster.df.add_pair.filter.filter$cellid_rna[colorder]],show_colnames = FALSE,show_rownames = TRUE,cluster_rows = FALSE, cluster_cols = FALSE,color = color_tfdev, breaks =  breaks)



##order by pre-calculated trajectory
dfD.1.1 <- readRDS('../../../02.snapATAC_harmony/pseudotime_ordering_Granja/dfD.1.1.order.rds')
#7798 
dfD.1.2 <- readRDS('../../../02.snapATAC_harmony/pseudotime_ordering_Granja/dfD.1.2.order.rds')
#7320

mat.atac <- matGS #37956 x 10041
mat.rna <- matRNA  #24307 x 10198



table(row.names(dfD.1.1) %in%  colnames(mat.atac)  )#no
table(row.names(dfD.1.2) %in%  colnames(mat.atac)  )#no
    
cellid <- row.names(dfD.1.1)
cellid <- row.names(dfD.1.2)
sum(duplicated(cellid )) #0
    
idx1 <- grep(pattern = 'donor1',x = cellid )
idx2 <- grep(pattern = 'donor2',x = cellid )
#unique(idx1+idx2)#1
length(c(idx1,idx2)) #7320 #7798
    
cellid[idx1] <- gsub(pattern = '_donor1#',replacement = '.atac_',x = cellid[idx1])

cellid[idx2] <- gsub(pattern = '-1$',replacement = '-2',x = cellid[idx2])
cellid[idx2] <- gsub(pattern = '_donor2#',replacement = '.atac_',x = cellid[idx2])

table(cellid %in%  colnames(mat.atac) )

FALSE  TRUE 
  760  7038

FALSE  TRUE 
  706  6614 

row.names(dfD.1.1) <- cellid
row.names(dfD.1.2) <- cellid
    

##read in cicero aggregate (low overlap) groups
cicero_KNN_groups <- readRDS('data/save-cicero-KNN-Groupings-cds.rds')    
length(unique(unlist(cicero_KNN_groups) ) ) #22782 #11091




#plot trajectory correlation heatmap

#gene.test <- var.genes[1:500] 

#gene.test<- c("FOS", "JUNB", "BACH1", "SMARCC1", "TP63", "TOPORS")

gene.test <- c("FOS", "JUNB",'ATF3', "STAT4", "STAT5A", "STAT5B", "STAT3")


##for dfD.1.1
mat.data <- mat.atac[gene.test,]
mat.data.z <- scale_rows(mat.data) # #zscore by row

table(colnames(mat.data.z) %in% rownames(dfD.1.1) )
FALSE  TRUE 
 3003  7038

mat.data.z.sel <- mat.data.z[,colnames(mat.data.z) %in% rownames(dfD.1.1)]
mat.data.z.sel.order <- mat.data.z.sel[,order(dfD.1.1[colnames(mat.data.z.sel),'pseudotime'])]
#pseudotime.1.1 <- dfD.1.1[colnames(mat.data.z),'pseudotime']
#mat.data.z.traj.1.1 <- mat.data.z[,match()]

    
#plot with hclust and color breaks
options(repr.plot.width=8,repr.plot.height=10)
options(repr.plot.width=4,repr.plot.height=5)

breaks <- seq(-2,2,by=4/(length(color_peak)-1))
res.p <- pheatmap::pheatmap(mat.data.z.sel.order,show_colnames = FALSE,show_rownames = TRUE,cluster_rows = TRUE, cluster_cols = FALSE,treeheight_col = 0,treeheight_row = 0,color = color_peak,breaks=breaks)
    
roworder <- res.p$tree_row$order
#colorder <- res.p$tree_col$order


mat.data <- mat.rna[gene.test,]
mat.data.z <- scale_rows(mat.data) # #zscore by row

rownames(matchedCells) <- matchedCells$x
matched.rna.traj <- matchedCells[colnames(mat.data.z.sel.order),'y']
#matched.rna.traj <- matchedCells.filter0.3[colnames(mat.data.z.sel.order),'y']
matched.rna.traj <- matched.rna.traj[!is.na(matched.rna.traj)]
mat.data.z.order <- mat.data.z[gene.test[roworder],matched.rna.traj]

breaks <- seq(-2,2,by=4/(length(color_tfdev)-1))
pheatmap::pheatmap(mat.data.z.order,show_colnames = FALSE,show_rownames = TRUE,cluster_rows = FALSE, cluster_cols = FALSE,color = color_tfdev, breaks =  breaks)


##for dfD.1.2
mat.data <- mat.atac[gene.test,]
mat.data.z <- scale_rows(mat.data) # #zscore by row

table(colnames(mat.data.z) %in% rownames(dfD.1.2) )
# FALSE  TRUE 
#  3427  6614 
mat.data.z.sel <- mat.data.z[,colnames(mat.data.z) %in% rownames(dfD.1.2)]
mat.data.z.sel.order <- mat.data.z.sel[gene.test[roworder],order(dfD.1.2[colnames(mat.data.z.sel),'pseudotime'])]
    

#plot with hclust and color breaks
options(repr.plot.width=8,repr.plot.height=10)
options(repr.plot.width=4,repr.plot.height=5)
    
breaks <- seq(-2,2,by=4/(length(color_peak)-1))
res.p <- pheatmap::pheatmap(mat.data.z.sel.order,show_colnames = FALSE,show_rownames = TRUE,cluster_rows = FALSE, cluster_cols = FALSE,treeheight_col = 0,treeheight_row = 0,color = color_peak,breaks=breaks)
    
# roworder <- res.p$tree_row$order
# #colorder <- res.p$tree_col$order


mat.data <- mat.rna[gene.test,]
mat.data.z <- scale_rows(mat.data) # #zscore by row

rownames(matchedCells) <- matchedCells$x
matched.rna.traj <- matchedCells[colnames(mat.data.z.sel.order),'y']
#matched.rna.traj <- matchedCells.filter0.3[colnames(mat.data.z.sel.order),'y']
matched.rna.traj <- matched.rna.traj[!is.na(matched.rna.traj)]
mat.data.z.order <- mat.data.z[gene.test[roworder],matched.rna.traj]

breaks <- seq(-2,2,by=4/(length(color_tfdev)-1))
pheatmap::pheatmap(mat.data.z.order,show_colnames = FALSE,show_rownames = TRUE,cluster_rows = FALSE, cluster_cols = FALSE,color = color_tfdev, breaks =  breaks)


table(dfD.1.1$Group) #PAPPA+
#   6    3    9    5    2    4 
# 1509 1090  395 1475 1783 1546 

table(dfD.1.2$Group) #FLT1+

#    6    3    9    5    7    1 
# 1509 1090  395 1475 1019 1832 
    
    
    


    
    
    
    
    
    
    
    
    
    
    
    

    


    
    
    
    
    
    
    
    
    
    


# ############start to customized plot################

# UMAP = cluster.df.add[,c('UMAP_1','UMAP_2')]

# left <- 1.3*min(UMAP[,1])
# right <- 0.7*max(UMAP[,1])
# top <- 1*max(UMAP[,2])
# bottom <- 1*min(UMAP[,2])


# #######visualize integrated library#########
# options(repr.plot.width=7.5,repr.plot.height=7.5)
# #ggplot(data=cluster.all.coord,aes(x=UMAP1,y=UMAP2,col=type) ) +
# ggplot(data=cluster.df.add,aes(x=UMAP_1_,y=UMAP_2,col=anno) ) +
# #ggplot(data=cluster.all.coord,aes(x=UMAP1,y=UMAP2,col=anno) ) +
#   geom_point(size=0.2,alpha=0.9) +
#   scale_color_manual(values = c('rna_D1'='#94C6DD', 'rna_D2'='#1273AE'  , 'atac_D1'='pink','atac_D2'='#C80927'),name='8w library') +
#   xlim(left,right) + ylim(bottom,top) +
#   #theme_classic() +
#   theme(legend.position = 'top',
#         axis.text=element_blank(), 
#         axis.title = element_text(size = 18, face = "bold"),
#         axis.ticks = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(color="black", fill = NA,size=1),
#         plot.title = element_text(size = 15, face = "bold"),
#         #complete = TRUE
#         plot.margin = unit(c(1,1,1,1), "lines"),
#         legend.text = element_text(size  = 15),
#         legend.title = element_text(size  = 15)
#        ) +
#   guides(colour=guide_legend(override.aes=list(size=8)))
# ggsave(filename = "pdfs_te/PLA-8w-RNA-ATAC-liger.source_of_donor.pdf",height=7.5,width=7.5,useDingbats=FALSE)

# # ggplot(data=cluster.all.coord,aes(x=dim1,y=dim2,col=clusters_int) ) +
# #   geom_point(size=0.5) +
# #   scale_color_manual(values = color_good) +
# #   theme_classic() +
# #   guides(colour=guide_legend(override.aes=list(size=6)))


# # ggplot(data=cluster.all.coord,aes(x=dim1,y=dim2,col=clusters_int) ) +
# #   geom_point(size=0.5) +
# #   scale_color_manual(values = color_good) +
# #   theme_classic() +
# #   guides(colour=guide_legend(override.aes=list(size=6)))


# ######plot the integrated clusters
# options(repr.plot.width=7.5,repr.plot.height=7.5)
# ggplot(data=cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster) ) +
# #ggplot(data=cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster_lib) ) +
# #ggplot(data=cluster.all.coord,aes(x=UMAP1,y=UMAP2,col=clusters_int) ) +
#   geom_point(size=0.2) +
#   scale_color_manual(values = color_good) +
#   xlim(left,right) + ylim(bottom,top) +
#   #theme_classic() +
#   theme(legend.position = 'none',
#         axis.text=element_blank(), 
#         axis.title = element_text(size = 25, face = "bold"),
#         axis.ticks = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(color="black", fill = NA,size=1),
#         plot.title = element_text(size = 15, face = "bold"),
#         #complete = TRUE
#         plot.margin = unit(c(1,1,1,1), "lines"),
#         legend.text = element_text(size  = 15),
#         legend.title = element_text(size  = 15)
#        ) +
#   guides(colour=guide_legend(override.aes=list(size=8)))
# ggsave(filename = "pdfs_te/PLA-8w-RNA-ATAC-liger.umap.nolabel.pdf",height=7.5,width=7.5,useDingbats=FALSE)




# ############customized way to plot umap-cluster with text halo###########
# centers <- cluster.df.add %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
#         y = median(x = UMAP_2))

# centers_shift = data.frame()
# ##add little shift for text x and y, to plot text halo, borrow from snapATAC
# theta= seq(0, 2*pi, length.out=50)
# r=0.1
# strwidth = 0.5 
# strheight = 0.5
# xo <- r*strwidth # r*strwidth('A')
# yo <- r*strheight #r*strheight('A')
# for (i in seq_len(nrow(centers))){
#   for (j in theta) {
#         centers_shift = rbind(centers_shift,
#                               data.frame(
#                                   cluster=as.character(unlist(centers[i,'cluster'])),
#                                   x=centers[i,'x'] + cos(j)*xo, 
#                                   y=centers[i,'y'] + sin(j)*yo
#                                  )
#                        )
#       }
# }

# ####the UMAP plot with annotation
# #label right
# options(repr.plot.height=5,repr.plot.width=5)
# ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
#   geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
#   scale_colour_manual(values = color_good)  +
#   #scale_colour_manual(values = color_snap_mod1)  +
#   #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
#   #theme_classic() +
#   #theme_bw() +
#   xlim(left,right) + ylim(bottom,top) +
#   theme(
#         legend.position = 'right',
#         axis.text=element_blank(), 
#         axis.title = element_text(size = 15, face = "bold"),
#         axis.ticks = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(color="black", fill = NA,size=1),
# #         panel.background = element_rect(fill = "white", colour = "white", 
# #                 size = rel(1)),
#         #panel.border = element_blank(),
#         plot.title = element_text(size = 15, face = "bold"),
#         #complete = TRUE
#         plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
#        ) +
#  #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
#   ggtitle(paste(sample, "\ntotal cells:",nrow(cluster.df.add),  sep=" ") ) +
# #   geom_text(data = centers_shift, #the halo
# #             mapping = aes(x=x,y=y,label = cluster), 
# #             colour = "white", 
# #             size = 6.5) +
# #   geom_text(data = centers, 
# #             mapping = aes(x=x,y=y,label = cluster), 
# #             colour = "black", 
# #             size = 6) +
#   guides(col = guide_legend(override.aes = list(size = 3))) +  ##no effect ??
#   #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
#   #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
#   labs(x = "UMAP1", y = "UMAP2")
#   #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

# ggsave(filename = "pdfs_te/PLA-8w-RNA-ATAC-liger.UMAP.labelright.pdf",height=5,width=6,useDingbats=FALSE)

# ##label on cluster
