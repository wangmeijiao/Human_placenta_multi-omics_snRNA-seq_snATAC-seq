

##find conserved marker gene among samples (integrated)
##find differential expressed gene among samples

##visualize marker gene, TF, hormone gene expression in vivo and in vitro by heatmap (with split and annotation)


library('Seurat')
packageVersion('Seurat') #4.1.1


library('ggplot2')
library('pheatmap')
library('ComplexHeatmap')

library(Rcpp) #1.0.9

library(ggtern)

##install.packages('ggtern')
#also installing the dependencies ‘lifecycle’, ‘rlang’, ‘vctrs’, ‘tensorA’, ‘bayesm’, ‘ggplot2’, ‘compositions’, ‘latex2exp’, ‘proto’





#BiocManager::install("sva")


#install metap for FindConservedMarkers
#install.packages('qqconf') #need R >= 4.0.0
#BiocManager::install('multtest')
#install.packages('metap')


color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E",'grey' )




####read in color palette list###
color_set_yellowbrick.flat <- readRDS('/home/mjwang/progs/misc-tools/r/color_set_yellowbrick.flat.rds')
color_list_archr <- readRDS('/home/mjwang/progs/misc-tools/r/ArchR.color_list.rds')


options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_list_archr)){ 
  len = length(color_list_archr[[name]])
  color = color_list_archr[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}

options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_set_yellowbrick.flat)){ 
  len = length(color_set_yellowbrick.flat[[name]])
  color = color_set_yellowbrick.flat[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}


color_set_BuenColors <- readRDS('/home/mjwang/progs/misc-tools/r/color_set_BuenColors.rds')


options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_set_BuenColors)){ 
  len = length(color_set_BuenColors[[name]])
  color = color_set_BuenColors[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}


color_use <- color_set_yellowbrick.flat[['OrRd.6']]
color_use <- color_set_yellowbrick.flat[['RdPu.6']]
color_use <- color_set_yellowbrick.flat[['Greens.6']]
color_use <- color_set_yellowbrick.flat[['Purples.6']]
color_use <- color_set_yellowbrick.flat[['Blues.6']]
color_use <- color_set_yellowbrick.flat[['Reds.6']]
color_use <- rev(color_set_yellowbrick.flat[['RdYlBu.10']])

color_use <- colorRampPalette(BuenColors::jdb_palette('brewer_red'))(10)#(256)
color_use1 <- c('grey95',colorRampPalette(BuenColors::jdb_palette('brewer_red'))(256) )
color_use2 <- c('grey95',colorRampPalette(BuenColors::jdb_palette('brewer_blue'))(256) )


color_use <- color_set_BuenColors[['color_set1']]
color_discrete <- color_list_archr[['summerNight']] #seven discrete color, use this


barplot(1:length(color_use),col = color_use,cex.main=2)
barplot(1:length(color_use1),col = color_use1,cex.main=2)
barplot(1:length(color_use2),col = color_use2,cex.main=2)



pretty_color <- list(
       '0' = '#00441B',
        '1' = '#46A040',
        '2' = '#00AF99',
         '3' ='#FFC179',
        '4' = '#98D9E9',
        '5' = '#F6313E',
        '6' = '#FFA300',
        '7' = '#C390D4',
        '8' = '#FF5A00',
        '9' = '#AFAFAF',
        '10' = '#7D7D7D',
        '11' = '#4B4B4B',
        '12' = '#8F1336',
        '13' = '#0081C9',
        '14' = '#001588',
        '15' = '#490C65',
        '16' = '#BA7FD0',
        '17' = '#E7D654',
       '18' =  '#6F1482',
       '19' =  'navy'
    
    
)


pretty_color <- unlist(pretty_color)
names(pretty_color) <- NULL


color_set <- color_use



sample <- 'integration-invivo-vitroSTB'



##Granja 's Rcpp correlation
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




###########

placenta.integrated <-  readRDS('placenta.integrated.early6-BL-CT30.final.final.reclustering.rds')
##placenta.integrated <-  readRDS('placenta.integrated.early6-BL-CT30.final.final.rds')#load seurat 3.6.3 CCA integrated placenta villus and CT30 ,BL celline data

33478 x 31369 

gc()

#cluster.df.add <- readRDS('cluster.df.add.final.final.rds')
cluster.df.add <- readRDS('cluster.df.add.final.final.reclustering.rds')

placenta.integrated[['cluster.sample.new']] <- cluster.df.add$cluster.sample.new
table(placenta.integrated[['cluster.sample.new']])


  STB.BL_1   STB.BL_2   STB.BL_3   STB.BL_4   STB.BL_5   STB.BL_6 STB.CT30_1 
      1266        874        648        603        492         47       1215 
STB.CT30_2 STB.CT30_3 STB.CT30_4 STB.CT30_5 STB.CT30_6 STB.CT30_7 STB.CT30_8 
       967        725        205        318        675        230        322 
STB.CT30_9   villus_1  villus_10  villus_11   villus_2   villus_3   villus_4 
       228       2407       2649        508       2188       2081       4941 
  villus_5   villus_6   villus_7   villus_8   villus_9 
      1810       1733       1654       1375       1208 




# DefaultAssay(placenta.integrated) <- "integrated"


# ###fill scale.data from placenta.integrated.full


# expr.scale.data <- readRDS('exprMat.integrate.scale.data.rds')
# 25150 x 33827

# table (  colnames(expr.scale.data) %in% colnames(placenta.integrated) )

# FALSE  TRUE 
#  2458 31369

# table ( colnames(placenta.integrated)   %in%  colnames(expr.scale.data)   )
#  TRUE 
# 31369

# table ( rownames(placenta.integrated)   %in%  rownames(expr.scale.data)   )
# TRUE 
# 2000 


# all.equal (colnames(placenta.integrated)  ,  colnames(placenta.integrated@assays$integrated@scale.data) )
# #TRUE

# placenta.integrated@assays$integrated@scale.data  <-  expr.scale.data[,colnames(placenta.integrated)]





# ###fill data from placenta.integrated.full

# expr.data <- readRDS('exprMat.integrated.data.rds')
# 25150 x 33827


# all.equal(dimnames(expr.scale.data ),dimnames(expr.data) )
# #TRUE

# placenta.integrated@assays$integrated@data  <-  expr.data[,colnames(placenta.integrated)]







##check embedding
options(repr.plot.height=5.5, repr.plot.width = 15)
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'sample_meta',split.by = 'sample_meta') + NoLegend()

options(repr.plot.height=5.5, repr.plot.width = 5.5)
DimPlot(placenta.integrated, reduction = "umap_rotate",label=TRUE,cols=color_good,label.size = 8,pt.size = .1) 

#placenta.integrated <- RenameIdents(object = placenta.integrated, '0' = '1', '1' = '2', '2'='3','3'='4' ,  '4'='5',   '5'='6','6'='7','7'='8','8'='9','9'='10'   )
options(repr.plot.height=5.5, repr.plot.width = 5.5)
DimPlot(placenta.integrated, reduction = "umap_rotate",label=TRUE,cols=color_good,label.size = 8,pt.size = .1) 


all.equal(placenta.integrated[['cluster']][,1],Idents(placenta.integrated), check.attributes = FALSE ) #TRUE
#cluster is the new integration cluster

all.equal(placenta.integrated[['cluster']][,1],cluster.df.add$cluster, check.attributes = FALSE ) 


#look at new cluster of integration with splitting
options(repr.plot.height=5.5, repr.plot.width = 15)
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'cluster',split.by = 'sample_meta') + NoLegend()

#look at cluster.sample with splitting
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'cluster.sample',split.by = 'sample_meta') #+ NoLegend()

#look at cluster.sample with splitting
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'cluster.sample.new',split.by = 'sample_meta') #+ NoLegend()




##:::::::::::::on map differential expressed gene test:::::::::::###### 


######identify conserved marker gene in the integrated space
#https://satijalab.org/seurat/articles/integration_introduction.html
# For performing differential expression after integration, we switch back to the original
# data

#DefaultAssay(placenta.integrated) <- "RNA" #full gene set, instead of 2000 anchor.feature gene

DefaultAssay(placenta.integrated) <- "integrated" #fill scale.data for rerun IntegrateData with all gene


table(Idents(placenta.integrated))
 5    4    2    8   11   10    1    7    3    9    6 
2407 2924 3550 1361 1179 1523 5032 2249 4928 3860 2356

 1    2    3    4    5    6    7    8    9   10 
6264 5375 4729 3677 2776 2428 2079 1524 1357 1160 

#  1    2    3    4    5    6    7    8    9   10 
# 6273 5418 4750 3721 2846 2447 2088 1629 1366 1160 

# 0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 4180 3948 3804 3748 3730 2449 2440 2075 2060 1645 1104 1049  576  420


table(placenta.integrated[['sample_meta']])
 STB.BL STB.CT30   villus 
    3930     4885    22554

#    STB.BL STB.CT30   villus 
#     3930     5214    22554 

table(placenta.integrated[['sample_meta']][,1],placenta.integrated[['cluster']][,1])
              5    4    2    8   11   10    1    7    3    9    6
  STB.BL     84  321  298  481   11 1365  394  158  568  202   48
  STB.CT30  700 1928  953  381    9  156  192  159  228  163   16
  villus   1623  675 2299  499 1159    2 4446 1932 4132 3495 2292


#              1    2    3    4    5    6    7    8    9   10
#   STB.BL    597  464  174  312  298   94  136 1367  479    9
#   STB.CT30  339  207   76 1032 1840  708  138  155  381    9
#   villus   5328 4704 4479 2333  638 1626 1805    2  497 1142

#              1    2    3    4    5    6    7    8    9   10
#   STB.BL    597  464  174  312  298   94  136 1367  479    9
#   STB.CT30  348  250   97 1076 1910  727  147  260  390    9
#   villus   5328 4704 4479 2333  638 1626 1805    2  497 1142


# ##merge CT30 and BL?
# sample_meta <- placenta.integrated[['sample_meta']]

# table(sample_meta)
#  STB.BL STB.CT30   villus 
#   3388     5651    24189

##or filter and keep only placenta.villus and CT30 ?
placenta.integrated.villus_CT30 <- subset(placenta.integrated, subset =  sample_meta %in% c('STB.CT30', 'villus')  )

table(placenta.integrated.villus_CT30[['sample_meta']])

STB.CT30   villus 
    4885    22554

# STB.CT30   villus 
#     5214    22554

# STB.CT30   villus 
#     5651    24189



##conserved marker gene among placenta.villus CT30 BL
#share.markers.villus_vs_CT30_vs_BL.df <- FindConservedMarkers(placenta.integrated, ident.1 = '1',grouping.var = "sample_meta", verbose = TRUE)

#c9 (default assay use RNA and slot to data)
share.markers.villus_vs_CT30_vs_BL.df.c9 <- FindConservedMarkers(placenta.integrated, ident.1 = '9',grouping.var = "sample_meta", verbose = TRUE)
head(share.markers.villus_vs_CT30_vs_BL.df.c9)

rownames(share.markers.villus_vs_CT30_vs_BL.df.c9)[1:10]
'PRKCE''PAPPA''FNDC3A''PAG1''ARHGAP26''KIAA1211''ICA1''FER1L5''PTPRJ''SLCO2A1'

#'ADAMTS6''VGLL3''NECTIN3''ELMO1''GULP1''SCIN''PRKCE''CSGALNACT1''RYBP''CHSY1'

saveRDS(share.markers.villus_vs_CT30_vs_BL.df.c9,'DEGs_conserved/share.markers.villus_vs_CT30_vs_BL.df.c9.rds')



#c3
share.markers.villus_vs_CT30_vs_BL.df.c3 <- FindConservedMarkers(placenta.integrated, ident.1 = '3',grouping.var = "sample_meta", verbose = TRUE)

head(share.markers.villus_vs_CT30_vs_BL.df.c3)
rownames(share.markers.villus_vs_CT30_vs_BL.df.c3)[1:10]

'ADAMTS6''VGLL3''GULP1''SCIN''PRKCE''RYBP''FNDC3A''MAN1A2''SH3TC2-DT''FBLN1'

saveRDS(share.markers.villus_vs_CT30_vs_BL.df.c3,'DEGs_conserved/share.markers.villus_vs_CT30_vs_BL.df.c3.rds')


#c7
share.markers.villus_vs_CT30_vs_BL.df.c7 <- FindConservedMarkers(placenta.integrated, ident.1 = '7',grouping.var = "sample_meta", verbose = TRUE)

head(share.markers.villus_vs_CT30_vs_BL.df.c7)
rownames(share.markers.villus_vs_CT30_vs_BL.df.c7)[1:10]
'PFKP''FSTL3''SH3BP5''MYCNUT''MYCN''MYCNOS''PLCL2''LIMD1''MTSS2''GDPD4'

#grep('FLT1',rownames(share.markers.villus_vs_CT30_vs_BL.df.c7))


saveRDS(share.markers.villus_vs_CT30_vs_BL.df.c7 ,'DEGs_conserved/share.markers.villus_vs_CT30_vs_BL.df.c7.rds')


##c1
share.markers.villus_vs_CT30_vs_BL.df.c1 <- FindConservedMarkers(placenta.integrated, ident.1 = '1',grouping.var = "sample_meta", verbose = TRUE)

head(share.markers.villus_vs_CT30_vs_BL.df.c1)
rownames(share.markers.villus_vs_CT30_vs_BL.df.c1)[1:10]
'CYP19A1''LIMCH1''SASH1''ADAM12''ANXA4''CBLB''LYN''ENG''HTRA4''PLIN2'


#grep('FLT1',rownames(share.markers.villus_vs_CT30_vs_BL.df.c1))

saveRDS(share.markers.villus_vs_CT30_vs_BL.df.c1,'DEGs_conserved/share.markers.villus_vs_CT30_vs_BL.df.c1.rds')


##conserved marker gene among placenta.villus CT30

share.markers.villus_vs_CT30.c9.df <- FindConservedMarkers(placenta.integrated.villus_CT30, ident.1 = '9', grouping.var = "sample_meta", verbose = TRUE)
#share.markers.villus_vs_CT30.c1.df <- FindConservedMarkers(placenta.integrated.villus_CT30, ident.1 = '1', grouping.var = "sample_meta", verbose = TRUE)
head(share.markers.villus_vs_CT30.c9.df)


share.markers.villus_vs_CT30.c7.df <- FindConservedMarkers(placenta.integrated.villus_CT30, ident.1 = '7', grouping.var = "sample_meta", verbose = TRUE)
head(share.markers.villus_vs_CT30.c7.df)



#share.markers.villus_vs_CT30.df <- share.markers.villus_vs_CT30.c8.df
#share.markers.villus_vs_CT30 <- rownames(share.markers.villus_vs_CT30.df)



#saveRDS(share.markers.villus_vs_CT30.df,'share.markers.villus_vs_CT30_vs_BL.rds')
#459 × 17 #c6 villus_vs_CT30_vs_BL
#557 x 17 #c8 villus_vs_CT30





###manually plot conserved marker gene by assay

#DefaultAssay(placenta.integrated) <- "integrated"
DefaultAssay(placenta.integrated) <- "RNA"


#for(i in share.markers.villus_vs_CT30[1:10]){
#for(i in share.markers.villus_vs_CT30.tf){
#for(i in rownames(share.markers.villus_vs_CT30_vs_BL.df)[1:10]){
##for(i in sample(marker.gene.df$gene,size = 5)  ){

#for(i in c('PAPPA','FLT1','SH3TC2','ERVFRD-1','TEAD4')) {
#for(i in c('AR','ESRRG','PAPPA','DNMT1','PSG8','ADAMTSL1')) {

##important gene
for(i in c('STAT5A','MITF','STAT5B','CEBPB','FOSL2','PAPPA','FLT1')) {
  options(repr.plot.height=5.5,repr.plot.width=15)
  res.p <- FeaturePlot(placenta.integrated, features = i, reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ scale_color_gradientn(colours = c('grey',color_use)) 
  print(res.p)

}



##split plot with cutoff and colorset
gene = 'PAPPA'
gene = 'FLT1'

#options(repr.plot.height=12.5,repr.plot.width=12.5)
options(repr.plot.height=5.5,repr.plot.width=5.5)
FeaturePlot(subset(placenta.integrated,  subset = sample_meta == "STB.BL"), features = gene, reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q90',split.by = 'sample_meta') + NoAxes()

options(repr.plot.height=5.5,repr.plot.width=5.5)
FeaturePlot(subset(placenta.integrated,  subset = sample_meta == "STB.CT30"), features = gene, reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q30',max.cutoff = 'q99',split.by = 'sample_meta') + NoAxes()

options(repr.plot.height=5.5,repr.plot.width=5.5)
FeaturePlot(subset(placenta.integrated,  subset = sample_meta == "villus"), features = gene, reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q30',max.cutoff = 'q99',split.by = 'sample_meta') + NoAxes()



##similar gene (find conserved marker list?)

options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'ESRRG', reduction = "umap_rotate",pt.size = .1,slot = 'data', min.cutoff = 'q50',max.cutoff = 'q99',split.by = 'sample_meta') 

FeaturePlot(placenta.integrated, features = 'PAPPA', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 

FeaturePlot(placenta.integrated, features = 'ADAMTSL1', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 

FeaturePlot(placenta.integrated, features = 'ADAMTS6', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q30',max.cutoff = 'q99',split.by = 'sample_meta') 


FeaturePlot(placenta.integrated, features = 'PRKCE', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q30',max.cutoff = 'q99',split.by = 'sample_meta') 


FeaturePlot(placenta.integrated, features = 'FNDC3A', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q30',max.cutoff = 'q99',split.by = 'sample_meta')  

FeaturePlot(placenta.integrated, features = 'PAG1', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q30',max.cutoff = 'q99',split.by = 'sample_meta') 







FeaturePlot(placenta.integrated, features = 'FLT1', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 

FeaturePlot(placenta.integrated, features = 'LEP', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 


FeaturePlot(placenta.integrated, features = 'MYCN', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 


FeaturePlot(placenta.integrated, features = 'CYP19A1', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q10',max.cutoff = 'q99',split.by = 'sample_meta')

FeaturePlot(placenta.integrated, features = 'FSTL3', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q10',max.cutoff = 'q99',split.by = 'sample_meta')



FeaturePlot(placenta.integrated, features = 'LIMD1', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 


FeaturePlot(placenta.integrated, features = 'ENG', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 



##not so similar gene

options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'DNMT1', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 

FeaturePlot(placenta.integrated, features = 'PSG8', reduction = "umap_rotate",pt.size = .5,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 


FeaturePlot(placenta.integrated, features = 'INHBA', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q10',max.cutoff = 'q99',split.by = 'sample_meta') 


##diffent gene

options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'CSH2', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 
FeaturePlot(placenta.integrated, features = 'CSH1', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 
FeaturePlot(placenta.integrated, features = 'CSHL1', reduction = "umap_rotate",pt.size = 1,slot = 'data', min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') 



# ##tf gene
# gene.tf <- c('CEBPB','MYCN','BACH1','MAFF')
# gene.tf <- c('AHR','ZNF117','RUNX1')

# for(i in gene.tf){
#   options(repr.plot.height=5.5,repr.plot.width=15)
#   res.p <- FeaturePlot(placenta.integrated, features = i, reduction = "umap",pt.size = .8,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') + scale_color_gradientn(colours = c('grey',color_use)) 
#   print(res.p)

# }








######compare different genes in vivo and in vitro by heatmap#################



##marker gene
marker.gene <- list(
    #'Quality control' = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.ribo', 'percent.hb', 'percent.chrY', 'percent.xist','XIST'),
    'Trophoblast' = c('KRT7', 'GATA3', 'TFAP2A','TFAP2C'), #KRT7 = CK7
    'CTB proliferation' = c('TEAD4','BRCA1','MKI67','TOP2A'),
    'CTB' = c('DNMT1', 'CDH1','EGFR','PPARG', 'TEAD3', 'TP53','TP63', 'TP73', 'BCAM','ITGA6'),
    'CTB fusion' = c('ERVFRD-1', 'GCM1', 'OVOL1','ANXA1','PPARD'), 
    'STB nascent' = c('SH3TC2', 'BACE2', 'ESRRG','ESR1'), 
    'STB general' = c('PSG8', 'CGA', 'PSG2', 'PSG5', 'LEP'), 
    'STB Mature 1' = c('PAPPA', 'ADAMTSL1','ADAMTS6', 'GH2', 'GHR', 'JAK2', 'LAMA3', 'AR', 'VDR', 'CSHL1', 'CSH1', 'CSH2', 'STAT5A', 'STAT5B', 'STAT4','AHR','MITF'), 
    'STB Mature 2' = c('FLT1', 'ENG', 'ANGPTL4', 'FSTL3', 'INHBA', 'INHA','MYCN', 'POU2F3', 'LVRN', 'TGFB1', 'FOSL2', 'JUND','FOS', 'FOSB', 'JUNB', 'JUN','HDAC2',"CEBPB"),
    #'STB apoptosis' = c('DDX60', 'DDX58', 'MAP4K4', 'SPATA5', 'GDF15', 'CROT', 'CDKN1A', 'ADCY5'),
    #'EVT' = c('HLA-G', 'LAIR2', 'PLAC8', 'MMP2'), 
    #'STR general' = c('VIM', 'DLK1'), 
    #'Vascular Endothelial Cell' = c('PECAM1'), 
    #'STR' = c('HIVEP3', 'HLA-A', 'HLA-DPA1', 'HLA-DPB1'), 
    #'Mesenchymal STR'= c('THY1'), 
    #'Hofbauer Cell'= c('CD68', 'CD14'), 
    #'Red blood' = c('HBA1', 'HBZ')
    'conservedMarker' = c('PRKCE','FNDC3A','PAG1','ARHGAP26', 'PFKP','SH3BP5','PLCL2','LIMD1')
)


marker.gene.df <- Reduce(f = rbind.data.frame,lapply(names(marker.gene), function(x){ data.frame(type=x, gene = marker.gene[[x]] )     }  ) )

sum(duplicated(marker.gene.df$gene)) #0

rownames(marker.gene.df) <- marker.gene.df$gene




###hormone gene 
hormone_list = c("INSIG1","STC2","ACTN4","ACTN1","PGF","ANGPT4","GAL","ANGPT1","ADM","FSTL3","FSTL1","IGF2","ANGPTL4","CSH1","CSH2","POMC","TAC3","FST","INSIG2","ANG","INHA","PSG2","PSG5","PSG3","CTGF","STC1","CGA","INSL4","LEP","PSG9","PSG4","CGB7","CGB8","CGB5","CGB2","CGB1","GH2","CSHL1","INHBA","PSG11","PSG1","PSG6","ANGPT2","LHB","KISS1","CRH","PSG8","CCK","PSG10P","PSG7","IGF1","RETN","GAST","GH1","PRL","REN","ANGPTL2","ANGPTL1","RLN2","EPO","INSL3","RLN3","EDN1","EDN3","GCG","BDNG","GNRH1","ADM2","INSL6","NMU","RLN1","APLN","AGT","FSTL4","NPPA","UTS2","NTS","UCN","SST","GNRH2","ANGPTL3","PPY","NPPB","ANGPTL16","UCN2","AMH","OXT","GIP","NPPC","CALCB","INHBE","SCT","PTHLH","ANGPTL7","ANGPTL5","PNOC","VIP","TAC1","CORT","ADIPOQ")#,'PAPPA','FLT1','GCM1')
#100

hormone.peptide.gene <- list( #see ppt peptide hormone, from Yawei 
   'Metabolism' = c('ADIPOQ','CALCA','CALCB','CRH','CSH1','CSH2','CSHL1','GAL','GCG','GH1','GH2','GHRH','GHRL','GIP','IGF1','INS','INSL3','INSL4','INSL5','INSL6','LEP','PPY','PRL','PTH','PTHLH','RETN','RETNLB','SOST','SST','STC1','STC2','TRH','TSHB','UCN','UCN2','VIP'),
    'VASCULAR' = c('ADM','ADM2','AGT','ANGPTL3','APLN','CTGF','EDN1','EDN3','NPPA','NPPB','NPPC','NTS2'),
    'BLOOD_CELL' = c('EPO','THPO'),
    'GUT_PEPTIDE' = c('CCK','GAST','MLN','SCT','TAC1'),
    'NEUROPEPTIDE' = c('BDNF','CROT','NMU','NTS','PMCH','PNOC'),
    'REPRODUCTION' = c('AMH','CGA','CGB7','FSHB','FST','GNRH2','GPHA2','INHA','INHBA','INHBB','INHBC','INHBE','KISS1','LHB','OXT','RLN1','RLN2')


)
#78


length(intersect(hormone_list, unlist(hormone.peptide.gene) ))
#58


#data.table::rbindlist(hormone.peptide.gene)
hormone.peptide.gene.df <- Reduce(f = rbind.data.frame,lapply(names(hormone.peptide.gene), function(x){ data.frame(type=x, gene = hormone.peptide.gene[[x]] )     }  ) )

sum(duplicated(hormone.peptide.gene.df)) #0

rownames(hormone.peptide.gene.df) <- hormone.peptide.gene.df$gene




##identify and plot TF heatmap pairwise heatmap with CT30 vs placental.villus, along cluster.sample of placenta.villus?###############


tflist_rowid <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.seurat_harmony/tflist_row_order.rds') #194

sum(duplicated(tflist_rowid)) #0


tflist <- read.table('/home/mjwang/dataex/cisTargetDb_v10/tf_lists/allTFs_hg38.txt',header = FALSE, stringsAsFactors = FALSE)$V1 #1892

sum(duplicated(tflist)) #0



length(intersect(tflist_rowid,tflist)) #190

tflist_rowid[!tflist_rowid %in% tflist]
'ARID2''GHR''NPAS3''LIN54' #not in tflist





##hypoxia gene
hypoxia.gene.1 <- c('HDAC2', 'EGLN1', 'EGLN3','ANGPTL4','SLC11A2', 'SLC2A1', 'LIMD1', 'HK2', 'NDRG1','HILPDA','ATG7' ,'CARD16')

hypoxia.gene.2 <- c('FLT1', 'VEGFA', 'LDHA', 'HIF1A', 'VHL', 'EDN1', 'EP300', 'TAF1', 'RPA1','MDM2', 'IGFBP3', 'FHL2', 'CDKN1A')


hypoxia.gene <- c(hypoxia.gene.1,hypoxia.gene.2)

sum(duplicated(hypoxia.gene)) #0


##FindConservedMarker result


intersect(rownames(share.markers.villus_vs_CT30_vs_BL.df.c9), tflist)
'AHR''PDLIM5''TCF7L2''TEAD1''FHL2''NFE2L3''PHF20''THRB''LRRFIP1''NCOA1''ZNF431''ZFHX3''JUND'


share.markers.villus_vs_CT30.tf <- intersect(share.markers.villus_vs_CT30,tflist)
'MYCN''MSI2''GATA2''MYLK''NCOA3''CREB3L2''DAB2''SCMH1''HMG20B''ZFHX3''THRB''NR4A3''PHF20''RFX7''TEAD4''NCOA1''CEBPB''SP3''CUX1''MAFF''TBX3''LRRFIP1''ZNF43''TBX20''HMGB1''RB1''ZNF652''RBFOX2''MEF2A''MEIS2''KMT2A''HUNK''ZNF431''E2F3''RBPJ''BACH1' #c6  placenta.villus CT30 BL


#'ZNF117''AHR''ARID3A''PDLIM5''RUNX1''DNMT3A''MYLK''CREB5''DAB2''TCF7L2''RBMS1''TEAD1''DNMT1''NR2F2''ZNF331''FHL2''MSI2''TGIF1''UGP2''ZNF440''PHF20''AEBP2''TBX20''TP63''PARP1''THRB''LRRFIP1''HUNK''KMT2A''RB1''NCOA1''RFX7''MEF2A''RBFOX2''GLIS3''ELF1''ZNF43''ZNF69''JUND''CUX1''GATA2''NAP1L1''FOS''ZNF626''FOXP1''MEIS2''TET1''ZFHX3' #c8 placenta.villus CT30




###################get average exprMat  ################
#(use slot data ave of assay RNA  )


#villus
placenta.integrated.villus <- subset(placenta.integrated, subset =  sample_meta %in% c('villus')  )

table(placenta.integrated.villus[['sample_meta']])
villus 
 22554 

# villus 
#  24189


table(Idents(placenta.integrated.villus))
  5    4    2    8   11   10    1    7    3    9    6 
1623  675 2299  499 1159    2 4446 1932 4132 3495 2292

#  1    2    3    4    5    6    7    8    9   10 
# 5328 4704 4479 2333  638 1626 1805    2  497 1142


DefaultAssay(placenta.integrated.villus) #<- "RNA"

#dim(placenta.integrated.villus@assays$integrated@data)
#25150 x 22554

#dim(placenta.integrated.villus@assays$integrated@scale.data)
#25150 x 22554



##filter out low expression percentage gene

gene.expr.perc.villus <- rowMeans(GetAssayData( placenta.integrated.villus,slot="counts", assay="RNA" ) > 0) *100

gene.expr.perc.villus['PAPPA'] #72.9
gene.expr.perc.villus['FLT1'] #98.6
gene.expr.perc.villus['CSH1'] #72.04
gene.expr.perc.villus['CSH2'] #76.07




# 0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 3878 3781  118 2979 3246 1416 2162  895 1936 1558  482 1039  281  418 

#DefaultAssay(placenta.integrated.villus) <- "RNA" #for use of all raw count
#DefaultAssay(placenta.integrated.villus) <- "integrated" 

#exprMat.ave.z.villus <- AverageExpression(placenta.integrated.villus,slot = "scale.data",assays = 'integrated')[[1]] #average by cluster of integrated, 25150 x 10 2000 x 14
#exprMat.ave.z.villus <- AverageExpression(placenta.integrated.villus,slot = "data",assays = 'integrated')[[1]]

#exprMat.ave.z.villus <- AverageExpression(placenta.integrated.villus,slot = "scale.data",assays = 'RNA')[[1]] #use this?
exprMat.ave.z.villus <- AverageExpression(placenta.integrated.villus,slot = "data",assays = 'RNA')[[1]]#for convinent, still use 'exprMat.ave.z.villus'
#exprMat.ave.z.villus <- AverageExpression(placenta.integrated.villus,slot = "counts",assays = 'RNA')[[1]]#for convinent, still use 'exprMat.ave.z.villus'


exprMat.villus <- GetAssayData(placenta.integrated.villus,slot = "data")
31478 x 22554

##use Magic package to get imputation (smooth)
library(Rmagic) #v2.0.3.99 from git pip install magic-impute==3.0.0 #v1.4.0 need python magic-imputate, install magic-imputate package in r-reticular conda env:pip install magic-impute==1.4.0
exprMat.villus.magic <- Rmagic::magic(t(exprMat.villus),
                                 genes = NULL,#marker.genes,
                                 n.jobs=1, #broken with > 1!!
                                ) #need gene as column mat #batch effect?


timetag <- Sys.time() #quick

timetag


exprMat.villus.magic <- t(exprMat.villus.magic$result)

saveRDS(exprMat.villus.magic,'exprMat.villus.magic.rds')


# DefaultAssay(placenta.integrated.villus) #<- "RNA"
# exprMat.villus <- GetAssayData(placenta.integrated.villus,slot = "counts")

# exprMat.cpm.villus <- edgeR::cpm(exprMat.villus) ##use this, normalize to colsums,lib.size <- colSums(y)
# exprMat.log2cpm.villus <- edgeR::cpm(exprMat.villus,prior.count = 1, log=TRUE) #use this, log is log2

gc()

aggregateClusters <- function(idents,ids,data){
   data.res <- data.frame(row.names = row.names(data))
   if(ids == "all"){ids = levels(idents)}
   for(id in ids){
     barcodes.sel <- names(idents[idents == id])
     #data.sel <- Matrix::rowSums(data[,barcodes.sel])
     data.sel <- Matrix::rowMeans(data[,barcodes.sel])
     data.res <- cbind(data.res,data.sel)

   }
   colnames(data.res) <- ids
   return(data.res)

}

# idents <- Idents(placenta.integrated.villus)
# idents <- droplevels(idents)
# all.equal(colnames(exprMat.cpm.villus),names(idents)) #TRUE
# exprMat.ave.z.villus <- aggregateClusters (idents,'all',exprMat.cpm.villus) #31478 x 11

idents <- Idents(placenta.integrated.villus)
idents <- droplevels(idents)
all.equal(colnames(exprMat.villus.magic),names(idents)) #TRUE
exprMat.ave.z.villus <- aggregateClusters (idents,'all',exprMat.villus.magic) #31478 x 11





##BL
placenta.integrated.BL <- subset(placenta.integrated, subset =  sample_meta %in% c('STB.BL')  )

table(placenta.integrated.BL[['sample_meta']])


STB.BL 
  3930

# STB.BL 
#   3388

table(Idents(placenta.integrated.BL))
5    4    2    8   11   10    1    7    3    9    6 
  84  321  298  481   11 1365  394  158  568  202   48 

#  1    2    3    4    5    6    7    8    9   10 
#  597  464  174  312  298   94  136 1367  479    9 

# 0    1    2    3    4    5    6    7    8    9   10   11   12 
#   77   53 1695  165  152  411   56  364   16   39  198    2  160


gene.expr.perc.BL <- rowMeans(GetAssayData( placenta.integrated.BL,slot="counts", assay="RNA" ) > 0) *100

gene.expr.perc.BL['PAPPA'] #3.4
gene.expr.perc.BL['FLT1'] #1.8
gene.expr.perc.BL['CSH1'] #1.01
gene.expr.perc.BL['CSH2'] #0.45

#> 1.5

#exprMat.ave.z.BL <- AverageExpression(placenta.integrated.BL,slot = "scale.data",assays = 'integrated')[[1]] #average by cluster of integrated 25150 x 10

#exprMat.ave.z.BL <- AverageExpression(placenta.integrated.BL,slot = "data",assays = 'integrated')[[1]] 


#exprMat.ave.z.BL <- AverageExpression(placenta.integrated.BL,slot = "scale.data",assays = 'RNA')[[1]]#use this?
exprMat.ave.z.BL <- AverageExpression(placenta.integrated.BL,slot = "data",assays = 'RNA')[[1]]
#exprMat.ave.z.BL <- AverageExpression(placenta.integrated.BL,slot = "counts",assays = 'RNA')[[1]]

#exprMat.ave.z.BL #31478 × 11 



exprMat.BL <- GetAssayData(placenta.integrated.BL,slot = "data")
31478 x 3930

##use Magic package to get imputation (smooth)
#library(Rmagic) #v2.0.3.99 from git pip install magic-impute==3.0.0 #v1.4.0 need python magic-imputate, install magic-imputate package in r-reticular conda env:pip install magic-impute==1.4.0
exprMat.BL.magic <- Rmagic::magic(t(exprMat.BL),
                                 genes = NULL,#marker.genes,
                                 n.jobs=1, #broken with > 1!!
                                ) #need gene as column mat #batch effect?


timetag <- Sys.time() #quick

timetag


exprMat.BL.magic <- t(exprMat.BL.magic$result)

saveRDS(exprMat.BL.magic,'exprMat.BL.magic.rds')


gc()



idents <- Idents(placenta.integrated.BL)
idents <- droplevels(idents)
all.equal(colnames(exprMat.BL.magic),names(idents)) #TRUE
exprMat.ave.z.BL <- aggregateClusters (idents,'all',exprMat.BL.magic) #31478 x 11





##CT30
placenta.integrated.CT30 <- subset(placenta.integrated, subset =  sample_meta %in% c('STB.CT30')  )

table(placenta.integrated.CT30[['sample_meta']])

STB.CT30 
    4885

# STB.CT30 
#     5651

table(Idents(placenta.integrated.CT30))
  5    4    2    8   11   10    1    7    3    9    6 
 700 1928  953  381    9  156  192  159  228  163   16

# 1    2    3    4    5    6    7    8    9   10 
#  339  207   76 1032 1840  708  138  155  381    9

#   0    1    2    3    4    5    6    7    8    9   10   11   12   13 
#  225  114 1991  604  332  622  222  816  108   48  424    8  135    2 


gene.expr.perc.CT30 <- rowMeans(GetAssayData( placenta.integrated.CT30,slot="counts", assay="RNA" ) > 0) *100

gene.expr.perc.CT30['PAPPA'] #48
gene.expr.perc.CT30['FLT1'] #23
gene.expr.perc.CT30['CSH1'] #5.13
gene.expr.perc.CT30['CSH2'] #0.26

quantile(x = gene.expr.perc.CT30, probs = seq(0,1,0.1))
0% 0 10% 0 20% 0.0409416581371546 30% 0.122824974411464 40% 0.327533265097236 50% 0.921187308085977 60% 2.92732855680655 70% 7.98362333674514 80% 17.1873080859775 90% 36.6837256908905 100%100


#exprMat.ave.z.CT30 <- AverageExpression(placenta.integrated.CT30,slot = "scale.data",assays = 'integrated')[[1]] #average by cluster of integrated 25150 x 10
#exprMat.ave.z.CT30 <- AverageExpression(placenta.integrated.CT30,slot = "data",assays = 'integrated')[[1]]


#exprMat.ave.z.CT30 <- AverageExpression(placenta.integrated.CT30,slot = "scale.data",assays = 'RNA')[[1]]#use this?
exprMat.ave.z.CT30 <- AverageExpression(placenta.integrated.CT30,slot = "data",assays = 'RNA')[[1]]
#exprMat.ave.z.CT30 <- AverageExpression(placenta.integrated.CT30,slot = "counts",assays = 'RNA')[[1]]



exprMat.CT30 <- GetAssayData(placenta.integrated.CT30,slot = "data")
31478 x 4885

##use Magic package to get imputation (smooth)
#library(Rmagic) #v2.0.3.99 from git pip install magic-impute==3.0.0 #v1.4.0 need python magic-imputate, install magic-imputate package in r-reticular conda env:pip install magic-impute==1.4.0
exprMat.CT30.magic <- Rmagic::magic(t(exprMat.CT30),
                                 genes = NULL,#marker.genes,
                                 n.jobs=1, #broken with > 1!!
                                ) #need gene as column mat #batch effect?


timetag <- Sys.time() #quick

timetag


exprMat.CT30.magic <- t(exprMat.CT30.magic$result)

saveRDS(exprMat.CT30.magic,'exprMat.CT30.magic.rds')


gc()

idents <- Idents(placenta.integrated.CT30)
idents <- droplevels(idents)
all.equal(colnames(exprMat.CT30.magic),names(idents)) #TRUE
exprMat.ave.z.CT30 <- aggregateClusters (idents,'all',exprMat.CT30.magic) #31478 x 11






###get gene of interest##


gene.use <- marker.gene.df$gene #75 #65

gene.use <- hormone_list #100
#gene.use <- hormone.peptide.gene.df$gene #78

gene.use <- tflist #1892
gene.use <- tflist_rowid #194

gene.use <- hypoxia.gene #25

sum(duplicated(gene.use )) #0

length(gene.use )



table(gene.use %in% rownames(exprMat.ave.z.villus) )


FALSE  TRUE #tf long list
  279  1613

# FALSE  TRUE  #tf long list
#   122  1770 



FALSE  TRUE# hormone gene
   20    80

# FALSE  TRUE  #hormone gene
#     6    94

TRUE 
  75

# TRUE #marker gene
#   67

# TRUE  #marker gene
#   65

# TRUE 
#   25

# TRUE 
#   65


# TRUE 
#   59


# FALSE  TRUE 
#     6    94

# FALSE  TRUE 
#     6    97

# FALSE  TRUE 
#     8    70

# FALSE  TRUE 
#   116  1776

# # TRUE 
# #  194

# FALSE  TRUE 
#     8    70


table(gene.use %in% rownames(exprMat.ave.z.CT30) )
#the same

table(gene.use %in% rownames(exprMat.ave.z.BL) )
#the same



gene.use.sel <- gene.use[gene.use %in% rownames(exprMat.ave.z.villus)]


exprMat.ave.z.villus.sel <- exprMat.ave.z.villus[gene.use.sel,c( '5','4','2','8','11','10','1','7','3','9','6' )]
all.equal(rownames(exprMat.ave.z.villus.sel) ,gene.use.sel) #TRUE

colnames(exprMat.ave.z.villus.sel) <- paste0('villus_',colnames(exprMat.ave.z.villus.sel))



exprMat.ave.z.CT30.sel <- exprMat.ave.z.CT30[gene.use.sel,c('5','4','2','8','11','10','1','7','3','9','6')]
all.equal(rownames(exprMat.ave.z.CT30.sel) ,gene.use.sel) #TRUE

colnames(exprMat.ave.z.CT30.sel) <- paste0('CT30_',colnames(exprMat.ave.z.CT30.sel))


#idx <- which(rownames(exprMat.ave.z.CT30.sel) != c(deg.c10,deg.c3))
#rownames(exprMat.ave.z.CT30.sel)[idx] <- c(deg.c10,deg.c3)[idx]
#colnames(exprMat.ave.z.CT30.sel) <- c('CT30_Mature1','CT30_Mature2')


exprMat.ave.z.BL.sel <- exprMat.ave.z.BL[gene.use.sel,c('5','4','2','8','11','10','1','7','3','9','6')]
all.equal(rownames(exprMat.ave.z.BL.sel) ,gene.use.sel) #TRUE
#idx <- which(rownames(exprMat.ave.z.BL.sel) != c(deg.c10,deg.c3))
colnames(exprMat.ave.z.BL.sel) <- paste0('BL_',colnames(exprMat.ave.z.BL.sel))
#colnames(exprMat.ave.z.BL.sel) <- c('BL_Mature1?','BL_Mature2')






all.equal(rownames(exprMat.ave.z.villus.sel) ,gene.use.sel) #TRUE
all.equal(rownames(exprMat.ave.z.CT30.sel) ,gene.use.sel) #TRUE
all.equal(rownames(exprMat.ave.z.BL.sel) ,gene.use.sel) #TRUE



##get normalize factor, with CTB as reference
gene <- 'STAT5A'
gene <- 'MITF'
gene <- 'CEBPB'
gene <- 'FOSL2'


gene <- 'PAPPA'
gene <- 'FLT1'

gene <- 'ADAMTS6'


gene <- 'LAMA3'

data.villus <- unlist(exprMat.ave.z.villus.sel[gene,])
data.CT30 <- unlist(exprMat.ave.z.CT30.sel[gene,])
data.BL <- unlist(exprMat.ave.z.BL.sel[gene,])



# ##get normalize factor from CTB data point (CEBPB)
# ave.villus <- mean(data.villus[c('villus_5','villus_4','villus_2')])
# ave.ct30 <- mean(data.CT30[c('CT30_5','CT30_4','CT30_2')])
# ave.bl <- mean(data.BL[c('BL_5','BL_4','BL_2')])

# norm.f.ct30 <- ave.ct30/ave.villus
# norm.f.bl <- ave.bl/ave.villus

# data.CT30 <- data.CT30/norm.f.ct30
# data.BL <- data.BL/norm.f.bl


##zscore
data.villus <- (data.villus-mean(data.villus))/sd(data.villus)
data.CT30 <- (data.CT30-mean(data.CT30))/sd(data.CT30)
data.BL <- (data.BL-mean(data.BL))/sd(data.BL)


##get ymin ymax
ymax.villus <- max(data.villus)
ymax.CT30 <- max(data.CT30)
ymax.BL <- max(data.BL)
ymax <- Reduce(max,list(ymax.villus, ymax.CT30 ,ymax.BL) )

ymin.villus <- min(data.villus)
ymin.CT30 <- min(data.CT30)
ymin.BL <- min(data.BL)
ymin <- Reduce(min,list(ymin.villus, ymin.CT30 ,ymin.BL) )



pdf( paste0("pdfs/heatmap_invivo_vitro_compare/compare_villus_CT30_BL.",gene,".pdf"),width=5.5, height=5.5)

options(repr.plot.width=5.5, repr.plot.height=5.5)
plot(data.villus,lty=1,type='b',col='red',ylab = 'expression', xlab = 'cell type',main = paste0(gene, ' expression along differention', sep ='' ),ylim =c(ymin,ymax  ) , xaxt  = 'n' )
points(data.CT30, lty=1,type='b',col='blue')
points(data.BL, lty=1,type='b',col='green')
axis(side = '1',at = 1:length(data.villus), labels = gsub(pattern = "villus_",replacement = "cluster_",x=names(data.villus)), las = 2 )
legend("topleft",legend = c('early villus','hTSC-CT30','hTSC-BL'),fill = c('red','blue','green') )

dev.off()



# exprMat.ave.z.CT30.sel <- exprMat.ave.z.CT30.sel/norm.f.ct30
# exprMat.ave.z.BL.sel <- exprMat.ave.z.BL.sel/norm.f.bl


# exprMat.ave.z.CT30.sel <- exprMat.ave.z.CT30.sel/norm.f.ct30
# exprMat.ave.z.BL.sel <- exprMat.ave.z.BL.sel/norm.f.bl




# ##filter by row correlation?

# exprMat.ave.z.villus.sel #gene x group


# mat.data1 <- exprMat.ave.z.villus.sel
# rownames(mat.data1) <- paste0('villus_',rownames(mat.data1) )

# mat.data2 <- exprMat.ave.z.CT30.sel
# rownames(mat.data2) <- paste0('CT30_',rownames(mat.data2) )

# #all(rownames(mat.data1) == rownames(mat.data2))


# matchedGenes <- data.frame(x = rownames(mat.data1), y = rownames(mat.data2) )


# matchedGenes$corVar <- rowCorCpp( #correlated with column data (all gene data)
#   match(matchedGenes$x, rownames(mat.data1)), 
#   match(matchedGenes$y, rownames(mat.data2)), 
#   as.matrix(mat.data1), 
#   as.matrix(mat.data2)
# )


# matchedGenes.sel <- matchedGenes.sel[rowSums(is.na(matchedGenes.sel)) ==0,]

# matchedGenes.sel <- matchedGenes[matchedGenes$corVar > 0.8,] #214

# matchedGenes.sel$x_gene <- sapply(stringr::str_split(matchedGenes.sel$x,pattern = "_",n=2), function(x){x[2]} )
# matchedGenes.sel$y_gene <- sapply(stringr::str_split(matchedGenes.sel$y,pattern = "_",n=2), function(x){x[2]} )


# all.equal(matchedGenes.sel$x_gene,matchedGenes.sel$y_gene)
# #TRUE


# gene.cor0.8 <- matchedGenes.sel$x_gene

# #gene.cor0.8.tf <- gene.cor0.8[gene.cor0.8 %in% gene.use] #205




########get combined exprMat.ave.z
# exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel)
# title <- 'villus only'

# exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel,exprMat.ave.z.CT30.sel)#,exprMat.ave.z.BL.sel)

# title <- 'villus vs CT30'

# exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel,exprMat.ave.z.BL.sel)
# title <- 'villus vs BL'




#use this
exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel,exprMat.ave.z.CT30.sel,exprMat.ave.z.BL.sel)#,exprMat.ave.z.BL.sel)

#exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel,exprMat.ave.z.CT30.sel,exprMat.ave.z.BL.sel*2)#,exprMat.ave.z.BL.sel)


#exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel,exprMat.ave.z.CT30.sel,exprMat.ave.z.BL.sel)



exprMat.ave.z.villus.sel.z <- scale_rows(exprMat.ave.z.villus.sel)
exprMat.ave.z.CT30.sel.z <- scale_rows(exprMat.ave.z.CT30.sel)
exprMat.ave.z.BL.sel.z <- scale_rows(exprMat.ave.z.BL.sel)

exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel.z,exprMat.ave.z.CT30.sel.z,exprMat.ave.z.BL.sel.z)



#####
#exprMat.ave.z.sel.combine <- exprMat.ave.z.sel.combine[gene.cor0.8,]



##filter gene with low expression and NA

exprMat.ave.z.sel.combine.filter <- exprMat.ave.z.sel.combine[!rowSums(exprMat.ave.z.sel.combine == 0) > 28 ,]


exprMat.ave.z.sel.combine.filter <- exprMat.ave.z.sel.combine.filter[!rowSums(is.na(exprMat.ave.z.sel.combine.filter)) > 1,]

sum (matrixStats::rowSds( as.matrix(exprMat.ave.z.sel.combine.filter) ) < 0.001)

exprMat.ave.z.sel.combine.filter <- exprMat.ave.z.sel.combine.filter[!matrixStats::rowSds( as.matrix(exprMat.ave.z.sel.combine.filter) ) < 0.001,]



geneid <- rownames(exprMat.ave.z.sel.combine.filter)

flag.villus <- gene.expr.perc.villus[geneid] > 1.5
flag.BL <- gene.expr.perc.BL[geneid] >1.5
flag.CT30 <- gene.expr.perc.CT30[geneid] >1.5

flag.combine <- flag.villus & flag.BL & flag.CT30

geneid.filter <- geneid[flag.combine]

grep("STAT5A",geneid.filter)
grep("CEBPB",geneid.filter)
grep("FOSL2",geneid.filter)
grep("MITF",geneid.filter)
grep("PSG8",geneid.filter)

#geneid.filter <- geneid.filter[-35]

exprMat.ave.data <- exprMat.ave.z.sel.combine.filter[geneid.filter,]


#exprMat.ave.data <- exprMat.ave.z.sel.combine

#title <- 'villus vs CT30 vs BL (raw counts no ComBat)'
#title <- 'villus vs CT30 vs BL (raw counts with ComBat)'

title <- 'villus vs CT30 vs BL (rna assay data slot\n MAGIC smooth, zscore then combine))'
#title <- 'villus vs CT30 vs BL (rna assay data slot\n MAGIC smooth, BL*2 combine, zscore ))'
#title <- 'villus vs CT30 vs BL (rna assay data slot\n MAGIC smooth, norm_by_CEBPB combine, zscore ))'
#title <- 'villus vs CT30 vs BL (rna assay data slot\n combine then zscore by row))'
#title <- 'villus vs CT30 vs BL (rna assay data slot\n combine then zscore by column))'
#title <- 'villus vs CT30 vs BL (rna assay data slot\n zscore then combine))'

#title <- 'villus vs CT30 vs BL (rna assay count slot\n log2cpm ave then combine))'
#title <- 'villus vs CT30 vs BL (rna assay count slot\n cpm ave then combine zscore))'
#title <- 'villus vs CT30 vs BL (rna assay count slot\n cpm ave then zscore combine))'

##combat normalize?

col_split <- factor(sapply(stringr::str_split(string = colnames(exprMat.ave.data),pattern = '_',n=2), function(x){x[1]} ),levels = c('villus','CT30','BL'))

#exprMat.ave.data.combat <- sva::ComBat(dat = exprMat.ave.data,batch = col_split)

#exprMat.ave.data <- exprMat.ave.data.combat


###zscore row data
scale_rows = function(x){ #pheatmap code
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

# scale_cols = function(x){ #pheatmap code
#   m = apply(x, 2, mean, na.rm = T)
#   s = apply(x, 2, sd, na.rm = T)
#   return((x - m) / s)
# }

#zscore by row
exprMat.ave.data.z <- scale_rows(exprMat.ave.data)

#zscore by column
#exprMat.ave.data.z <- t(scale_rows(t(exprMat.ave.data)))


#already zscored before combine by column
exprMat.ave.data.z <- exprMat.ave.data

##begin to plot



gene_type = 'marker gene'
#gene_type = 'hormone gene'
#gene_type = 'tf gene'
#gene_type = 'Hypoxia gene'

#title <- 'villus vs CT30 vs BL '

#title <- 'hTSCs CT30 vs villus'

#color_discrete <- color_discrete
#names(color_discrete) <- unique(hormone.peptide.gene.df[rownames(exprMat.ave.data.z ),'type',drop = TRUE])


#names(color_discrete) <- unique(marker.gene.df[rownames(exprMat.ave.data.z ),'type',drop = TRUE])



# ha_left <- rowAnnotation(
#                        df = hormone.peptide.gene.df[rownames(exprMat.ave.data.z),'type',drop = FALSE],
#                        col = list( 'type' = color_discrete)
    
    
#                      )

# ha_left <- rowAnnotation(
#                        df = marker.gene.df[rownames(exprMat.ave.data.z),'type',drop = FALSE],
#                        col = list( 'type' = color_discrete)
    
    
#                      )


ha_top <- columnAnnotation(
                       df = data.frame(dataset=col_split),
                       col = list( 'dataset' = c('villus'='#a64027','CT30'='#022336', 'BL' = '#2a7185') )
    
    
                     )



ha_top_text = HeatmapAnnotation(
        foo = anno_block(gp = gpar(fill = c('villus'='#a64027','CT30'='#022336', 'BL' = '#2a7185')),
                         labels = c('villus', 'CT30', 'BL'), 
                         labels_gp = gpar(col = "white", fontsize = 15)
                        )
)



map_cellname_cca <- list(
"villus_6" = "villus CTB proliferation",
"villus_5" = "villus CTB-1",
"villus_4" = "villus CTB-2",
"villus_9" = "villus CTB Fusion",
"villus_10" = "villus STB Nascent (villus)",
"villus_8" = "villus STB Nascent (celline)",
"villus_2" = "villus STB Premature 2",
"villus_7" = "villus STB Mature 2",
"villus_1" = "villus STB Mature 1",
"villus_3" = "villus STB Mixed",

"CT30_6" = "CT30 CTB proliferation",
"CT30_5" = "CT30 CTB-1",
"CT30_4" = "CT30 CTB-2",
"CT30_9" = "CT30 CTB Fusion",
"CT30_10" = "CT30 STB Nascent (villus)",
"CT30_8" = "CT30 STB Nascent (celline)",
"CT30_2" = "CT30 STB Premature 2",
"CT30_7" = "CT30 STB Mature 2",
"CT30_1" = "CT30 STB Mature 1",
"CT30_3" = "CT30 STB Mixed",

"BL_6" = "BL CTB proliferation",
"BL_5" = "BL CTB-1",
"BL_4" = "BL CTB-2",
"BL_9" = "BL CTB Fusion",
"BL_10" = "BL STB Nascent (BL)",
"BL_8" = "BL STB Nascent (celline)",
"BL_2" = "BL STB Premature 2",
"BL_7" = "BL STB Mature 2",
"BL_1" = "BL STB Mature 1",
"BL_3" = "BL STB Mixed"


)


all.equal(colnames(exprMat.ave.data.z), names(map_cellname_cca) ) #TRUE
colnames(exprMat.ave.data.z) <- sapply(map_cellname_cca[colnames(exprMat.ave.data.z)],function(x){x} )
colnames(exprMat.ave.data.z) <- gsub("\\.1",replacement = "",x = colnames(exprMat.ave.data.z))

n_km = 10
#n_km = 5
#n_km = 1
#n_km = 3
#n_km = 10 #for marker gene
#n_km = 10 #for hypoxia gene
#n_km = 10 #for tf gene

#options(repr.plot.height=30,repr.plot.width=9) #for TF correlation selection 

options(repr.plot.height=15,repr.plot.width=9) #for marker gene, hormone gene, hypoxia gene
#options(repr.plot.height=65,repr.plot.width=6) 
hp = ComplexHeatmap::Heatmap(exprMat.ave.data.z, 
                             name = 'zscore', 
             cluster_rows = TRUE, 
             cluster_columns = FALSE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             col = circlize::colorRamp2(seq(-0.5,2,by=2.5/(length(color_use)-1)), color_use),
             ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             #row_title = "hclust: ward.D",
             clustering_method_row = "ward.D", ##ward.D,complete
             show_row_dend = FALSE,
             row_km = n_km,
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 12),
             column_title = paste0(title,' \n',gene_type," expression"),
             #row_title = 'kmeans groups',
             column_title_gp = gpar(fontsize = 15),
             column_names_gp = gpar(fontsize = 12),
             ##row_split = hormone.peptide.gene.df[rownames(exprMat.ave.z.sel.combine.z),'type'],
             column_split = col_split,
             #title = paste0(title,' ',gene_type," gene expression")
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
             #left_annotation = ha_left,
             top_annotation = ha_top_text
) 
#width = max(grobWidth(textGrob(labels))))

##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
set.seed(1);hp = draw(hp, heatmap_legend_side = "left",gap = unit(0.1, "cm"),row_title = 'Kmeans group') #8.895833 8.843750


############get kmeans groupid from hp, then add left annotation and redraw################
set.seed(1);row_groupid_list <- row_order(hp)

#exprMat.ave.data.z[c(41,35,33),]


kmeans.group.df <- Reduce(f = rbind.data.frame,lapply(names(row_groupid_list), function(x){ data.frame(kmeans=x, geneidx = row_groupid_list[[x]] )     }  ) )

kmeans.group.df$gene <- rownames(exprMat.ave.data.z)[kmeans.group.df$geneidx]

groupid <- unique(kmeans.group.df[,'kmeans',drop = TRUE])

color_discrete_use <- c(color_discrete,c('orange','grey','navy'))
#color_discrete_use <- color_discrete[1:length(groupid)]

names(color_discrete_use) <- groupid

ha_left <- rowAnnotation(
                       df = kmeans.group.df[,'kmeans',drop = FALSE],
                       col = list( 'kmeans' = color_discrete_use[])
    
    
                     )

row_split <- kmeans.group.df[,'kmeans',drop = TRUE]

row_split <- factor(row_split, levels = groupid)

###redraw without km and row clustering


options(repr.plot.height=15,repr.plot.width=9)
#options(repr.plot.height=30,repr.plot.width=9) #for tf gene
hp = ComplexHeatmap::Heatmap(#exprMat.ave.data.z[kmeans.group.df[,'geneidx',drop = TRUE],], 
                             data_out,
                             name = 'zscore', 
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             show_row_names = FALSE,#TRUE,
             use_raster = FALSE,
             col = circlize::colorRamp2(seq(-0.5,2,by=2.5/(length(color_use)-1)), color_use),
             ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             #row_title = "hclust: ward.D",
             #clustering_method_row = "ward.D", ##ward.D,complete
             #show_row_dend = FALSE,
             #row_km = 10,
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 12),
             column_title = paste0(title,' \n',gene_type," expression"),
             #row_title = 'kmeans groups',
             column_title_gp = gpar(fontsize = 15),
             column_names_gp = gpar(fontsize = 12),
             row_split = row_split,
             column_split = col_split,
             #title = paste0(title,' ',gene_type," gene expression")
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
             left_annotation = ha_left,
             top_annotation = ha_top_text,
             right_annotation = ha_right
) 
#width = max(grobWidth(textGrob(labels))))

##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
set.seed(1);draw(hp, heatmap_legend_side = "left",gap = unit(0.1, "cm"),row_title = 'Kmeans group') #8.895833 8.843750


##save to pdf
#outf <- paste0( 'pdfs/',gsub(pattern = " ",replacement = "_",x=gene_type),'.',gsub(pattern = " |\\(|\\)",replacement = "_",x=title),'.pdf'  )
#outf <- "pdfs/heatmap_invivo_vitro_compare/marker_gene.villus_vs_CT30_vs_BL_rna_assay_data_slot_magic_smooth_zscore_then_combine.pdf"
outf <- "pdfs/heatmap_invivo_vitro_compare/marker_gene.villus_vs_CT30_vs_BL_rna_assay_data_slot_magic_smooth_zscore_then_combine.reannotation.pdf"
#outf <- "pdfs/marker_gene.villus_vs_CT30_vs_BL__rna_assay_data_slot_zscore_then_combine.pdf"
pdf(file = outf,height = 15,width=9,useDingbats = FALSE)
#pdf(file = outf,height = 30,width=9,useDingbats = FALSE) #for tf gene
set.seed(1);draw(hp, heatmap_legend_side = "left",gap = unit(0.1, "cm"),row_title = 'Kmeans group')
dev.off()




###output data mat for tenary plot??####


##reannotate rowid and colid
data_out <- exprMat.ave.data.z[kmeans.group.df[,'geneidx',drop = TRUE],]


colnames(data_out)

#'villus_5''villus_4''villus_2''villus_8''villus_11''villus_10''villus_1''villus_7''villus_3''villus_9''villus_6''CT30_5''CT30_4''CT30_2''CT30_8''CT30_11''CT30_10''CT30_1''CT30_7''CT30_3''CT30_9''CT30_6''BL_5''BL_4''BL_2''BL_8''BL_11''BL_10''BL_1''BL_7''BL_3''BL_9''BL_6'


col.list <- list('villus_5' = 'villus_CTB',
                  'villus_4' = 'villus_CTB',
                  'villus_2' = 'villus_CTB',
                  'villus_8' = 'villus_CTB_fusion',
                  'villus_11' = 'villus_STB',
                  'villus_10' = 'villus_STB',
                  'villus_1' = 'villus_STB',
                  'villus_7' = 'villus_STB',
                  'villus_3' = 'villus_STB',
                  'villus_9' = 'villus_STB',
                  'villus_6' = 'villus_STB',
                  'CT30_5' = 'CT30_CTB',
                  'CT30_4' = 'CT30_CTB',
                  'CT30_2' = 'CT30_CTB',
                  'CT30_8' = 'CT30_CTB_fusion',
                  'CT30_11' = 'CT30_STB',
                  'CT30_10' = 'CT30_STB',
                  'CT30_1' = 'CT30_STB',
                  'CT30_7' = 'CT30_STB',
                  'CT30_3' = 'CT30_STB',
                  'CT30_9' = 'CT30_STB',
                  'CT30_6' = 'CT30_STB',
                  'BL_5' = 'BL_CTB',
                  'BL_4' = 'BL_CTB',
                  'BL_2' = 'BL_CTB',
                  'BL_8' = 'BL_CTB_fusion',
                  'BL_11' = 'BL_STB',
                  'BL_10' = 'BL_STB',
                  'BL_1' = 'BL_STB',
                  'BL_7' = 'BL_STB',
                  'BL_3' = 'BL_STB',
                  'BL_9' = 'BL_STB',
                  'BL_6' = 'BL_STB'
                 )


colnames(data_out) <- unlist(col.list[colnames(data_out)])

##redraw with above code?


##reannotation rowid with highlight

# Set stylings for row names and make our selected rows unique

geneid_similar <- c('FLT1','ESRRG','ENG','CEBPB','PAPPA')
geneid_different <- c('FOSL2','PSG8','PSG2','LAMA3','MITF','STAT5A')

#similar
row_idx <- which(rownames(data_out) %in% geneid_similar)

#fontsizes <- rep(10, nrow(data_out))
#fontsizes[row_idx] <- 18

fontcolors <- rep('black', nrow(data_out))
fontcolors[row_idx] <- 'black'
fontfaces <- rep('plain',nrow(data_out))
fontfaces[row_idx] <- 'bold'

#diff
row_idx <- which(rownames(data_out) %in% geneid_different)

#fontsizes <- rep(10, nrow(data_out))
#fontsizes[row_idx] <- 18

#fontcolors <- rep('black', nrow(data_out))
fontcolors[row_idx] <- 'red'
#fontfaces <- rep('plain',nrow(data_out))
fontfaces[row_idx] <- 'bold'



# Create text annotation object for displaying row names
ha_right <- rowAnnotation(rows = anno_text(rownames(data_out), gp = gpar(fontsize = 12, fontface = fontfaces, col = fontcolors)))



###save data_out with grouping



saveRDS(data_out, 'heatmap_invivo_vitro_compare/data_out.rds')

data_out <- readRDS( 'heatmap_invivo_vitro_compare/data_out.rds')



openxlsx::write.xlsx( data_out, 
            file = 'heatmap_invivo_vitro_compare/data_out.xlsx',
            rowNames = TRUE
          )




colnames(data_out)
'villus_CTB' 1
'villus_CTB' 2
'villus_CTB' 3
'villus_CTB_fusion' 4
'villus_STB' 5
'villus_STB' 6
'villus_STB' 7
'villus_STB' 8
'villus_STB' 9
'villus_STB' 10
'villus_STB' 11
'CT30_CTB' 12
'CT30_CTB' 13
'CT30_CTB' 14
'CT30_CTB_fusion' 15
'CT30_STB' 16
'CT30_STB' 17
'CT30_STB' 18
'CT30_STB' 19
'CT30_STB' 20
'CT30_STB' 21
'CT30_STB' 22
'BL_CTB' 23
'BL_CTB' 24
'BL_CTB' 25
'BL_CTB_fusion' 26
'BL_STB' 27
'BL_STB' 28
'BL_STB' 29
'BL_STB' 30
'BL_STB' 31
'BL_STB' 32
'BL_STB' 33


data_mat_CTB_ave <- data.frame( "villus_CTB_ave" = rowMeans(data_out[,c(1,2,3)]),"CT30_CTB_ave" =  rowMeans(data_out[,c(11,12,13)]), "BL_CTB_ave" =  rowMeans(data_out[,c(23,24,25)]) )

data_mat_CTB_fusion_ave <- data.frame( "villus_CTB_fusion_ave" = data_out[,4],"CT30_CTB_fusion_ave" = data_out[,15], "BL_CTB_fusion_ave" = data_out[,26] )
rownames(data_mat_CTB_fusion_ave) <- rownames(data_out)

data_mat_STB_ave <- data.frame( "villus_STB_ave" = rowMeans(data_out[,5:11]),"CT30_STB_ave" =  rowMeans(data_out[,16:22]), "BL_STB_ave" =  rowMeans(data_out[,27:33]) )

write.table(data_mat_CTB_ave,file = 'heatmap_invivo_vitro_compare/data_mat_CTB_ave.txt',sep='\t',quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(data_mat_CTB_fusion_ave,file = 'heatmap_invivo_vitro_compare/data_mat_CTB_fusion_ave.txt',sep='\t',quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(data_mat_STB_ave,file = 'heatmap_invivo_vitro_compare/data_mat_STB_ave.txt',sep='\t',quote = FALSE, row.names = TRUE, col.names = TRUE)



