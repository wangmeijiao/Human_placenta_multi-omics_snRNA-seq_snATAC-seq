
##extended script (based on work_seurat_harmony_finalize.r) for identify DEG genes for single cell/nuclei RNA-seq
##find more de gene with more relax parameters
##find negative foldchange genes
##do in situ on UMAP for any given two clusters (also for intetegrated UMAP if applicable)

##do nearest knn bulk within cluster for different conditions (stage? disease vs normal) ? (Granja's method)

#wang meijiao @ 2022.1.24




##https://satijalab.org/seurat/reference/findallmarkers

library(Seurat)

library(edgeR)

library(magrittr)
library(dplyr)

library(ggplot2)
library(ggrepel)

library(patchwork)


outdir = 'DEGs_extend'
sample = "PLA-early_combined-RNA"

color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" )


##add six new color at the begining

##seaborn_pair <- c('c8_neg' = '#145b7d', 'c8_pos' ='#11264f','c5_neg' = '#78cdd1', 'c5_pos' ='#008792','c9_neg' = '#2a5caa','c9_pos' ='#102b6a','c10_neg' = '#b2df8a', 'c10_pos' ='#33a02c', 'c6_neg' = '#c85d44', 'c6_pos' ='#f15a22','c1_neg' = '#f58f98','c1_pos' = '#f15b6c', 'c3_neg' ='#fb9a99','c3_pos' = '#e31a1c', 'c2_neg' ='#b69968', 'c2_pos' ='#decb00', 'c4_neg' ='#fdbf6f','c4_pos' = '#ff7f00', 'c7_neg' ='#ffff99', 'c7_pos' ='#b15928' ) #modified for celltype aggreement
seaborn_pair <- c('c8_neg' = '#145b7d', 'c8_pos' ='#11264f','c5_neg' = '#78cdd1', 'c5_pos' ='#008792','c9_neg' = '#2a5caa','c9_pos' ='#102b6a','c11_neg' = '#2a5caa','c11_pos' ='#102b6a','c10_neg' = '#b2df8a', 'c10_pos' ='#33a02c', 'c6_neg' = '#c85d44', 'c6_pos' ='#f15a22','c1_neg' = '#f58f98','c1_pos' = '#f15b6c', 'c3_neg' ='#fb9a99','c3_pos' = '#e31a1c', 'c2_neg' ='#cab2d6', 'c2_pos' ='#6a3d9a', 'c4_neg' ='#fdbf6f','c4_pos' = '#ff7f00', 'c7_neg' ='#ffff99', 'c7_pos' ='#b15928' ) #original colorset

##sns.color_palette("Paired").as_hex()
barplot(1:length(seaborn_pair),col=seaborn_pair,names.arg = names(seaborn_pair), las =2 )




reds <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf')
blues <- c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','darkblue')
blues_darker <- darker(blues, 0.3)

oranges <- c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6')
purples <- c('#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')

colorset <- blues
colorset <- blues_darker
barplot(1:length(colorset),col=colorset,names.arg = 1:length(colorset), las =2 )

rna_map_cellcolor <- list(
    
    '1'=color_good[1],
    '2'=color_good[2],
    '3'=color_good[3],
    '4'=color_good[4],
    '5'=color_good[5],
    '6'=color_good[6],
    '7'=color_good[7],
    '8'=color_good[8],
    '9'=color_good[9],
    '10'=color_good[10],
    '11'= color_good[11]

    
##ok1   
#     '1'=purples[2],#'STB1',
#     '2'=purples[3],#'STB5', 
#     '3'=purples[4],#'STB2',
#     '4'=blues[2],#'STB3',
#     '5'=blues[5],#'CTB-2',
#     '8'=purples[1],#naive STB',
#     '7'= 'darkgreen',#'STB4',
#     '6'=blues[3],#'lightblue',#'CTB-1',
#     '9'=blues[4],#'CTB-3',
#     '11'='darkorange',#'Fusion component'
#     '10'= purples[5]

    
#     '1'=purples[2],#'STB1',
#     '2'=purples[3],#'STB5', 
#     '3'=purples[5],#'STB2',
#     '4'=purples[2],#'STB3',
#     '5'=blues[5],#'CTB-2',
#     '6'=purples[1],#naive STB',
#     '7'=purples[4],#'STB4',
#     '8'=blues[3],#'lightblue',#'CTB-1',
#     '9'=blues[4],#'CTB-3',
#     '10'='darkgreen'#'Fusion component'

    
#     '1'=purples[4], #STB4
#     '2'=purples[5],  #STB5
#     '3'=purples[2],  #STB2
#     '4'=purples[1], #STB1
#     '5'=purples[3], #STB3
#     '6'='#8B0000', #Syncytial knot
#     '7'='#d8daeb', #naive STB
#     '8'='#7f3b08', #STB-new
#     '9'='darkgreen', #CTB
#     '10'=''
 )


cellcolor <- unlist(rna_map_cellcolor )#[c('10','6','4','7','1','3','2')]

barplot(1:length(cellcolor),col=cellcolor,names.arg = names(cellcolor), las =2 )



cellcolor_darker <- darker(cellcolor, 0.3)
names(cellcolor_darker) <- names(cellcolor)

barplot(1:length(cellcolor),col=cellcolor,names.arg = names(cellcolor), las =2 )
barplot(1:length(cellcolor_darker),col=cellcolor_darker,names.arg = names(cellcolor_darker), las =2 )




darker <- function(colorset = NULL, ratio = NULL){
    ##darker a given color by ratio (lighter if negative)
    #https://graphicdesign.stackexchange.com/questions/75417/how-to-make-a-given-color-a-bit-darker-or-lighter
    ##color code to rgb -> get dist from 255 -> get increased or decreased 
    ##ratio can not > 1
    
    res.color <- vector()
    
    stopifnot(ratio <=1)
    
    for (color in colorset){
      rgb.df <- as.data.frame(col2rgb(color))
      rgb.df.do <- round( rgb.df * (1-ratio)   )
      #rgb.df.do <- round(rgb.df - abs(rgb.df - 255) * ratio)
      rgb <- unlist(rgb.df.do)
      rgb.hex <- paste("#",paste(as.hexmode(rgb),collapse=""),sep="")
      #rgb.df.do.hex <- apply(rgb.df.do, 2, as.hexmode )#  as.hexmode(rgb.df.do)
      res.color <- c(res.color, rgb.hex)
    }
    return(res.color)
    
    
    
}


quickDimPlot_labelon <- function(data = cluster.df.add, feature = 'cluster', color_use = NULL, title= '8w all', xlim = NULL, ylim = NULL ,shrink.x = 1.0, shrink.y = 1.0){

  options(repr.plot.height=7.5,repr.plot.width=8.5)
  centers <- data %>% dplyr::group_by_at(feature) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))   
  
  colnames(centers) <- c('cluster','x','y')
  
  centers_shift = data.frame()
  ##add little shift for text x and y, to plot text halo, borrow from snapATAC
  theta= seq(0, 2*pi, length.out=50)
  r=0.1
  strwidth = 0.5 
  strheight = 0.5
  xo <- r*strwidth # r*strwidth('A')
  yo <- r*strheight #r*strheight('A')
  for (i in seq_len(nrow(centers))){
    for (j in theta) {
          centers_shift = rbind(centers_shift,
                                data.frame(
                                    cluster=as.character(unlist(centers[i,'cluster'])),
                                    x=centers[i,'x'] + cos(j)*xo, 
                                    y=centers[i,'y'] + sin(j)*yo
                                   )
                         )
        }
  }

  
  if( is.null(xlim) ){ 
                       xmin = min(data[,'UMAP_1']); 
                       xmax = max(data[,'UMAP_1'])

                     }
    else{
        xmin <- xlim[1]
        xmax <- xlim[2]
        
    }
  if( is.null(xlim) ){ 
                       ymin = min(data[,'UMAP_2']); 
                       ymax = max(data[,'UMAP_2'])

                     }
    else{
        ymin <- ylim[1]
        ymax <- ylim[2]
        
    }

  p <- ggplot(data=data, aes_string(x = "UMAP_1", y = "UMAP_2", col = `feature`)) +
  #ggplot(data=cluster.df.add_TM_ok, aes(x = UMAP_1, y = UMAP_2, col = library)) +
  geom_point(size = .1) +
#   geom_text(data = centers, 
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "black", 
#             size = 10) +
  
  scale_colour_manual(values = color_use)  +
  guides(col = guide_legend(override.aes = list(size = 6))) + 
  #ggtitle(paste(title,', total nuclei: ',nrow(data),sep='') ) +
  #theme_bw() +
  theme(
        legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "six donors, total cells:",nrow(data),  sep=" ") ) +
  geom_text(data = centers_shift, #the halo
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "white", 
            size = 6.5) +
  geom_text(data = centers, 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            size = 6) +
  ##guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  ylim(shrink.y*ymin,shrink.y*ymax) + xlim(shrink.x*xmin,shrink.x*xmax) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
    
    
#   theme(legend.position = 'right',
#         axis.text=element_blank(), 
#         axis.title = element_text(size = 15, face = "bold"),
#         axis.ticks = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(color="black", fill = NA,size=1),
#         plot.title = element_text(size = 15, face = "bold"),
#         #complete = TRUE
#         plot.margin = unit(c(1,1,1,1), "lines")
#        )
  
  print(p)
  return('done')
}




##########readin data#####
# placenta <- readRDS('placenta.final.final.rds')
# cluster.df.add <- readRDS('cluster.df.add.final.rds')


placenta <- readRDS('snapshot_rds/placenta.rds') #placenta.cstb.rds in fact
cluster.df.add <- readRDS('snapshot_rds/cluster.df.add.rds') #cluster.df.add in fact


options(repr.plot.width = 7.5, repr.plot.height=7.5)
quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', color_use = rna_map_cellcolor, title= 'early combined_RNA', shrink.x = 2.2, shrink.y = 1.05)

DimPlot(object = placenta, label = TRUE,cols= color_good,pt.size = 0.5,label.size = 10,reduction = "umap_rotate") + NoLegend() 



all.equal(Idents(placenta), cluster.df.add$cluster ,check.attributes = FALSE)
#TRUE


placenta.stb <- subset(placenta,idents = c('8','4','2','10','5','1','3'))

Idents(placenta.stb) <- factor(Idents(placenta.stb), levels = c('8','4','2','10','5','1','3') )


######do findAllMarkers() and visualization######

#check the clusters
options(repr.plot.height = 7.5, repr.plot.width = 7.5)
DimPlot(placenta,label = TRUE,label.size = 15,reduction = 'umap_rotate') + NoLegend()

DimPlot(placenta.stb,label = TRUE,label.size = 15,reduction = 'umap_rotate') + NoLegend()

table(Idents(placenta))
1    2    3    4    5    6    7    8    9   10   11 
2504 2238 2229 5336 1857 1740 1659 1501 1275 2838  525

#   1    2    3    4    5    6    7    8    9   10 
# 1720 1428 1225 1093 1083  962  892  852  655  288 

table(Idents(placenta.stb))
 8    4    2   10    5    1    3 
1501 5336 2238 2838 1857 2504 2229

#  6    1    3    2    4    7 
#  962 1720 1225 1428 1093  892

#   1    2    3    4    6    7 
# 1720 1428 1225 1093  962  892

plot(cluster.df.add[,c('UMAP_1','UMAP_2')])



# #get average matrix by seurat function##
# #exprMat.ave <- AverageExpression(placenta, slot = 'data')[[1]]

# #zscore normalize
# #exprMat.ave.z <-  AverageExpression(placenta,slot = "scale.data")[[1]]#,return.seurat = TRUE)


# saveRDS(object = exprMat.ave,'exprMat.ave.rds')
# saveRDS(object = exprMat.ave.z,'exprMat.ave.z.rds')

# exprMat.ave <- readRDS('exprMat.ave.rds')
# exprMat.ave.z <- readRDS('exprMat.ave.z.rds') #24307 x 10 #use this (from TF expr vs TF dev dir)


##get aggregated mat
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

exprMat.z <- GetAssayData(placenta, slot = 'scale.data')
#all.equal(exprMat.z, readRDS('exprMat.scaledata.rds') ) #TRUE


all.equal(rownames(cluster.df.add),colnames(exprMat.z) ) #TRUE
cluster <- cluster.df.add$cluster
names(cluster) <- rownames(cluster.df.add)

cluster <- factor(cluster,levels = c('7','9','6','11','8','4','2','10','5','1','3' )  )

exprMat.z.aggre <- aggregateClusters(cluster,'all',exprMat.z) #already zscored
saveRDS(exprMat.z.aggre,'exprMat.z.aggre.rds')



all.equal(dimnames(exprMat.z.aggre), dimnames(exprMat.ave.z[,c('7','9','6','11','8','4','2','10','5','1','3' ) ]) ) #TRUE

all.equal(exprMat.z.aggre, exprMat.ave.z[,c('7','9','6','11','8','4','2','10','5','1','3' ) ]) 
#TRUE !


###zscore row data
scale_rows = function(x){ #pheatmap code
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

##as background geneset???
rownames(exprMat.z)



##do all DEG identification (will find negative genes by default)


##i very loosen method
#placenta.markers <- FindAllMarkers(placenta, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
placenta.markers.veryloose <- FindAllMarkers(placenta, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1) #very loose
#35830 x 7 
#10222 x 7, use this, will keep STAT4

placenta.markers.veryloose[grep('STAT5A',placenta.markers.veryloose$gene),] #hit in c10
placenta.markers.veryloose[grep('STAT4',placenta.markers.veryloose$gene),] #hit in c10

#saveRDS(placenta.markers.veryloose,'DEGs_extend/placenta.markers.veryloose.rds')
placenta.markers.veryloose <- readRDS('DEGs_extend/placenta.markers.veryloose.rds')



##ii loose method
placenta.markers.loose <- FindAllMarkers(placenta, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25) #loose
#12924 x 7
#saveRDS(placenta.markers.loose,'DEGs_extend/placenta.markers.loose.rds')
placenta.markers.loose <- readRDS('DEGs_extend/de.markers.all.loose.rds')

#timetag <- Sys.time()
#"2023-02-03 05:08:32 CST"

placenta.markers.loose[grep('STAT5A',placenta.markers.loose$gene),] #no hit
placenta.markers.loose[grep('STAT4',placenta.markers.loose$gene),] #hit in c10


###iii strict method
placenta.markers.strict <- FindAllMarkers(placenta, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25) #the same as orinigal de gene identification but output negative ones, will lost STAT4
#11676 x 7

placenta.markers.strict[grep('STAT5A',placenta.markers.strict$gene),] #no hit
placenta.markers.strict[grep('STAT4',placenta.markers.strict$gene),] #no hit


#saveRDS(placenta.markers.strict,'DEGs_extend/placenta.markers.strict.rds')
placenta.markers.strict <- readRDS('DEGs_extend/placenta.markers.strict.rds')




#de.markers.all <- placenta.markers.strict #11676
#de.markers.all <- placenta.markers.loose #12924 
de.markers.all <- placenta.markers.veryloose #35830, with STAT5A in



table(de.markers.all$avg_logFC > 0 )

FALSE  TRUE  #loose
 5594  7330

FALSE  TRUE  #veryloose
16814 19016




de.markers.all$cluster <- factor(de.markers.all$cluster, c('7','9','6','11','8','4','2','10','5','1','3'))
table(de.markers.all$cluster)

7    9    6   11    8    4    2   10    5    1    3  #loose
2891 2275 2388 1489  479  488  553  553  546  567  695

7    9    6   11    8    4    2   10    5    1    3 #veryloose
6183 5399 5524 3410 1964 2088 2056 2119 2173 2343 2571



for(i in c('7','9','6','11','8','4','2','10','5','1','3') ){
    cat('cluster',i," FALSE  TRUE\n")
    cat(table( subset(de.markers.all, cluster == i)$avg_logFC > 0 ),'\n' )
    
    
}

#loose
cluster 7  FALSE  TRUE
774 2117 
cluster 9  FALSE  TRUE
721 1554 
cluster 6  FALSE  TRUE
705 1683 
cluster 11  FALSE  TRUE
675 814 
cluster 8  FALSE  TRUE
344 135 
cluster 4  FALSE  TRUE
422 66 
cluster 2  FALSE  TRUE
365 188 
cluster 10  FALSE  TRUE
311 242 
cluster 5  FALSE  TRUE
354 192 
cluster 1  FALSE  TRUE
450 117 
cluster 3  FALSE  TRUE
473 222

#veryloose
cluster 7  FALSE  TRUE
1398 4785 
cluster 9  FALSE  TRUE
1317 4082 
cluster 6  FALSE  TRUE
1299 4225 
cluster 11  FALSE  TRUE
1860 1550 
cluster 8  FALSE  TRUE
1263 701 
cluster 4  FALSE  TRUE
1666 422 
cluster 2  FALSE  TRUE
1460 596 
cluster 10  FALSE  TRUE
1329 790 
cluster 5  FALSE  TRUE
1508 665 
cluster 1  FALSE  TRUE
1876 467 
cluster 3  FALSE  TRUE
1838 733




# table(placenta.markers$avg_logFC > 0 )
# FALSE  TRUE 
# 16814 19016

# # FALSE  TRUE 
# #  4682  5540

# # table(placenta.markers.pct_25$avg_logFC > 0 )
# # FALSE  TRUE 
# #  4332  5317 

# table(placenta.markers.loose$avg_logFC > 0 )
# FALSE  TRUE 
#  5594  7330


# # table(placenta.markers$cluster)
# #  8    5    9   10    6    1    3    2    4    7 
# # 2264 1946 1326 1283  507  582  609  516  513  676 

# # # 1    2    3    4    5    6    7    8    9   10 
# # #  582  516  609  513 1946  507  676 2264 1326 1283 


# # table(placenta.markers.pct_25$cluster)
# #   8    5    9   10    6    1    3    2    4    7 
# # 2163 1855 1285 1223  467  531  554  476  481  614 

# # #  1    2    3    4    5    6    7    8    9   10 
# # #  531  476  554  481 1855  467  614 2163 1285 1223 



# placenta.markers$cluster <- factor(placenta.markers$cluster, c('7','9','6','11','8','4','2','10','5','1','3'))
# table(placenta.markers$cluster)

# 7    9    6   11    8    4    2   10    5    1    3 
# 6183 5399 5524 3410 1964 2088 2056 2119 2173 2343 2571

# #    8    5    9   10    6    1    3    2    4    7 
# # 2264 1946 1326 1283  507  582  609  516  513  676


# # placenta.markers.pct_25$cluster <- factor(placenta.markers.pct_25$cluster, c('8','5','9','10','6','1','3','2','4','7'))
# # table(placenta.markers.pct_25$cluster)
# #   8    5    9   10    6    1    3    2    4    7 
# # 2163 1855 1285 1223  467  531  554  476  481  614

# placenta.markers.logfc0.25$cluster <- factor(placenta.markers.logfc0.25$cluster, c('7','9','6','11','8','4','2','10','5','1','3'))
# table(placenta.markers.logfc0.25$cluster)
# 7    9    6   11    8    4    2   10    5    1    3 
# 2891 2275 2388 1489  479  488  553  553  546  567  695


# for(i in c('7','9','6','11','8','4','2','10','5','1','3') ){
#     cat('cluster',i," FALSE  TRUE\n")
#     cat(table( subset(placenta.markers, cluster == i)$avg_logFC > 0 ),'\n' )
    
    
# }

# cluster 7  FALSE  TRUE
# 1398 4785 
# cluster 9  FALSE  TRUE
# 1317 4082 
# cluster 6  FALSE  TRUE
# 1299 4225 
# cluster 11  FALSE  TRUE
# 1860 1550 
# cluster 8  FALSE  TRUE
# 1263 701 
# cluster 4  FALSE  TRUE
# 1666 422 
# cluster 2  FALSE  TRUE
# 1460 596 
# cluster 10  FALSE  TRUE
# 1329 790 
# cluster 5  FALSE  TRUE
# 1508 665 
# cluster 1  FALSE  TRUE
# 1876 467 
# cluster 3  FALSE  TRUE
# 1838 733

# # cluster 8  FALSE  TRUE
# # 647 1617 
# # cluster 5  FALSE  TRUE
# # 637 1309 
# # cluster 9  FALSE  TRUE
# # 432 894 
# # cluster 10  FALSE  TRUE
# # 589 694 
# # cluster 6  FALSE  TRUE
# # 370 137 
# # cluster 1  FALSE  TRUE
# # 437 145 
# # cluster 3  FALSE  TRUE
# # 367 242 
# # cluster 2  FALSE  TRUE
# # 383 133 
# # cluster 4  FALSE  TRUE
# # 390 123 
# # cluster 7  FALSE  TRUE
# # 430 246 


# # for(i in c('8','5','9','10','6','1','3','2','4','7') ){
# #     cat('cluster',i," FALSE  TRUE\n")
# #     cat(table( subset(placenta.markers.pct_25, cluster == i)$avg_logFC > 0 ),'\n' )

# # }

# # cluster 8  FALSE  TRUE
# # 620 1543 
# # cluster 5  FALSE  TRUE
# # 609 1246 
# # cluster 9  FALSE  TRUE
# # 422 863 
# # cluster 10  FALSE  TRUE
# # 561 662 
# # cluster 6  FALSE  TRUE
# # 330 137 
# # cluster 1  FALSE  TRUE
# # 386 145 
# # cluster 3  FALSE  TRUE
# # 323 231 
# # cluster 2  FALSE  TRUE
# # 344 132 
# # cluster 4  FALSE  TRUE
# # 360 121 
# # cluster 7  FALSE  TRUE
# # 377 237 



# for(i in c('7','9','6','11','8','4','2','10','5','1','3') ){
#     cat('cluster',i," FALSE  TRUE\n")
#     cat(table( subset(placenta.markers.logfc0.25, cluster == i)$avg_logFC > 0 ),'\n' )
    
    
# }

# cluster 7  FALSE  TRUE
# 774 2117 
# cluster 9  FALSE  TRUE
# 721 1554 
# cluster 6  FALSE  TRUE
# 705 1683 
# cluster 11  FALSE  TRUE
# 675 814 
# cluster 8  FALSE  TRUE
# 344 135 
# cluster 4  FALSE  TRUE
# 422 66 
# cluster 2  FALSE  TRUE
# 365 188 
# cluster 10  FALSE  TRUE
# 311 242 
# cluster 5  FALSE  TRUE
# 354 192 
# cluster 1  FALSE  TRUE
# 450 117 
# cluster 3  FALSE  TRUE
# 473 222 





##for STB only
placenta.markers.stb <- FindAllMarkers(placenta.stb, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
#within stb clusters, one cluster vs all other (average)

#saveRDS(placenta.markers.stb,'DEGs_extend/de.markers.all.stb.rds')
placenta.markers.stb <- readRDS('DEGs_extend/de.markers.all.stb.rds')


table(placenta.markers.stb$cluster)
  8   4   2  10   5   1   3 
341  17 267 228 150 156 265 

# 6   1   3   2   4   7 
# 283  72 199  75 171 178 


placenta.markers.stb.perc_25 <- FindAllMarkers(placenta.stb, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

#saveRDS(placenta.markers.stb.perc_25,'DEGs_extend/de.markers.all.stb.perc_25.rds')
placenta.markers.stb.perc_25 <- readRDS('DEGs_extend/de.markers.all.stb.perc_25.rds')

table(placenta.markers.stb.perc_25$cluster)
  8   4   2  10   5   1   3 
315  13 253 211 142 137 250 

#  6   1   3   2   4   7 
# 275  69 188  73 164 170 

# ##general method with both negative and positive de genes
# placenta.markers.addnegative <- FindAllMarkers(placenta, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

# placenta.markers.addnegative <- FindAllMarkers(placenta, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)

# ##min.pct, any of the compared two must has percentage expressed (becareful for different genes with scarce expression)


table(placenta.markers.stb$avg_logFC > 0 )
FALSE  TRUE 
  737   687 

# FALSE  TRUE 
#   506   472

table(placenta.markers.stb.perc_25$avg_logFC > 0 )
FALSE  TRUE 
  665   656

# FALSE  TRUE 
#   484   455 



for(i in c('8','4','2','10','5','1','3' ) ){
    cat('cluster',i," FALSE  TRUE\n")
    cat(table( subset(placenta.markers.stb, cluster == i)$avg_logFC > 0 ),'\n' )
    
}

cluster 8  FALSE  TRUE
201 140 
cluster 4  FALSE  TRUE
17 
cluster 2  FALSE  TRUE
154 113 
cluster 10  FALSE  TRUE
60 168 
cluster 5  FALSE  TRUE
64 86 
cluster 1  FALSE  TRUE
114 42 
cluster 3  FALSE  TRUE
127 138

# cluster 6  FALSE  TRUE
# 185 98 
# cluster 1  FALSE  TRUE
# 47 25 
# cluster 3  FALSE  TRUE
# 54 145 
# cluster 4  FALSE  TRUE
# 124 47 
# cluster 7  FALSE  TRUE
# 46 132 
# cluster 2  FALSE  TRUE
# 50 25 

for(i in c('8','4','2','10','5','1','3' ) ){
    cat('cluster',i," FALSE  TRUE\n")
    cat(table( subset(placenta.markers.stb.perc_25, cluster == i)$avg_logFC > 0 ),'\n' )
    
}

cluster 8  FALSE  TRUE
176 139 
cluster 4  FALSE  TRUE
13 
cluster 2  FALSE  TRUE
144 109 
cluster 10  FALSE  TRUE
57 154 
cluster 5  FALSE  TRUE
60 82 
cluster 1  FALSE  TRUE
95 42 
cluster 3  FALSE  TRUE
120 130

# cluster 6  FALSE  TRUE
# 178 97 
# cluster 1  FALSE  TRUE
# 44 25 
# cluster 3  FALSE  TRUE
# 52 136 
# cluster 4  FALSE  TRUE
# 118 46 
# cluster 7  FALSE  TRUE
# 44 126 
# cluster 2  FALSE  TRUE
# 48 25 



# saveRDS(placenta.markers,'DEGs_extend/de.markers.all.rds')
# saveRDS(placenta.markers.pct_25,'DEGs_extend/de.markers.all.pct_25.rds')

# saveRDS(placenta.markers.stb,'DEGs_extend/de.markers.all.stb.rds')
# saveRDS(placenta.markers.stb.perc_25,'DEGs_extend/de.markers.all.stb.perc_25.rds')


#placenta.markers <- readRDS('DEGs_extend/de.markers.all.rds') #loose (perc: 0.1) all
#placenta.markers.pct_25 <- readRDS('DEGs_extend/de.markers.all.pct_25.rds') #strict all

#placenta.markers.stb <- readRDS('DEGs_extend/de.markers.all.stb.rds') #loose stb
#placenta.markers.stb.perc_25 <- readRDS('DEGs_extend/de.markers.all.stb.perc_25.rds') #strict stb


##########compare positive deg of full method deg with positive deg with positive only########


list2DF_keepref <- function(datalist = NULL, refid = 'de.PAPPA', plot = TRUE,title = 'compare pappa with previous', height = 8.5, width = 3.5){ 
    
    len.list <- length(datalist)
    
    ##get all shared id and length
    allid.shared <- Reduce(f = intersect, x = datalist)
    allid.combined <- unique(unlist(datalist))
    cat('shared ids length is: ',length(allid.shared),'\n',
        'shared ids are: ',paste(allid.shared,collapse = ','),'\n\n',
        'combined id length is ',length(allid.combined),
        sep = ''
       
       )

    #check
    #intersect(datalist[[1]],datalist[[2]])
    #intersect(datalist[[2]],datalist[[1]])
    
    rowid <- vector()
    ##get combined id (order by ref if set, if not set, use the sortedd allid.combined as row id)
    if( is.null(refid) ){
#       ref <- datalist[[1]]
        rowid <- sort(allid.combined)
    }else{ 
           ref = datalist[[refid]]
           notref <- allid.combined[!allid.combined %in% ref]
           if( length(notref) == 0 ){rowid = ref}else{
             rowid <- c(ref,notref)
           }
         }
    
    ##the count mat
    nes.df <- data.frame(matrix(data = 0, nrow = length(rowid), ncol = len.list ),row.names = rowid )
    colnames(nes.df) <- names(datalist)
    
    for(i in names(datalist) ){
        for(j in datalist[[i]] ){
        #for(j in names(datalist[[i]]) ){ 
            nes.df[j,i] <- nes.df[j,i]+1 #datalist[[i]][[j]]
        }
        
        
    }
    
    
    if(plot){
        options(repr.plot.width = width, repr.plot.height = height)
        res.p <- pheatmap::pheatmap(
                               nes.df, 
                               cluster_cols = FALSE,
                               cluster_rows = FALSE, 
                               treeheight_col = 2.5,
                               treeheight_row = 2.5,
                               #cellwidth = 20, 
                               #cellheight = 2.8,
                               na_col = 'grey', 
                               color = c('white','red','red','red'),#color_use,
                               breaks= c(-1,0,1,2),
                               border =TRUE,
                               border_color = 'black',
                               #labels_col = colid, 
                #                labels_row = make_bold_names(marker.mat.z,
                #                                             rownames, 
                #                                             rowid_hi
                #                                            ),
                               show_rownames = TRUE,

                               angle_col = 270,
                               #display_numbers = data.df.text,
                               #number_color = 'white',
                               #fontsize_number = 10,
                               fontsize = 15,
                               fontsize_col = 13,
                               fontsize_row = 10,
                               main = paste(' deg list compare\n',title,sep=''),
                              # main = paste(title,' gene expression zscore',sep=''),
                               silent = FALSE,
                               #scale = 'row',
                               scale = 'none',
                               #legend_breaks = c(2,4,6,8,10),
                               #legend_labels = c(2,4,6,8,10),
                               legend = FALSE
                               #legend = TRUE
                #                        height = height,
                #                        width = width,
                #                        filename = paste('result_do_GO_Pathway_quickonestep/ck.GO_BP.pvalue_',
                #                                         pvalue,'.qvalue_', qvalue,'.heatmap.pdf' ,sep = '' 
                #                                        )

                              )



         print(res.p)

    }
    
    return(nes.df)
    
}




deg_positive_only <-  marker.genes.de <- readRDS('DEGs/marker.genes.de.top25.sorted.rds') #positive only, very loose, top25, no c4

deg_positive_full <- subset(de.markers.all,avg_logFC > 0)


###compare deg with previous deg result

#for (i in levels(deg_positive_only$cluster)  ){
for (i in c('7','9','6','11','8','2','10','5','1','3') ){
    de.list <- list(
      #'deg_positive_only' = subset(deg_positive_only,cluster == '10')[1:25,'gene',drop=TRUE],
       #'deg_positive_full' = subset(deg_positive_full ,cluster == '10')[1:25,'gene',drop=TRUE]

       'deg_positive_only' = subset(deg_positive_only,cluster == i)[1:25,'gene',drop=TRUE],
       'deg_positive_full' = subset(deg_positive_full ,cluster == i)[1:25,'gene',drop=TRUE]


       #   'deg_positive_only' = subset(deg_positive_only,cluster == '3')[1:25,'gene',drop=TRUE],
       #'deg_positive_full' = subset(deg_positive_full ,cluster == '3')[1:25,'gene',drop=TRUE]


    )

    df.compare.top25 <- list2DF_keepref(datalist = de.list, refid = 'deg_positive_only',plot =TRUE , title = paste0('c',i,' top25 \n deg_positive_only\n vs deg_positive_full'), height = 10.5, width = 3.5 )
    ##nearly the same (include order), seurat findAllMarker positive only and full method
}






#######plot de table with grid volcano plot, labeling top 10 genes of both up/down expressed#########
##output de gene table for top50 top100 and top1000

res.de.list <- list()
##options(repr.plot.height = 10, repr.plot.width = 12)
##par(mfrow = c(2,3))
#for(i in levels(placenta.markers$cluster) ){
for(i in levels(de.markers.all$cluster) ){
#for(i in levels(placenta.markers.stb$cluster) ){

    
    cat('de gene for cluster ',i,'\n',sep= '')
    
    #sample <- paste('loose','use_all',sep="_")
    sample <- paste('veryloose','use_all',sep="_")
    #sample <- paste('strict','use_all',sep="_")
    #sample <- paste('loose','stb_only',sep="_")
    #sample <- paste('strict','stb_only',sep="_")

    if(!dir.exists( paste('DEGs_extend/de_gene_split/',sample,sep='') ) ){ dir.create(paste('DEGs_extend/de_gene_split/',sample,sep=''),recursive = TRUE) }
    if(!dir.exists( paste('DEGs_extend/pdfs/',sample,sep='') ) ){ dir.create(paste('DEGs_extend/pdfs/',sample,sep=''),recursive = TRUE) }
    
    
    stb.markers <- subset(de.markers.all, cluster == i)
    #stb.markers <- subset(placenta.markers, cluster == i) #use global de table 
    #stb.markers <- subset(placenta.markers.pct_25, cluster == i) #use global de table 
    #stb.markers <- subset(placenta.markers.stb, cluster == i)
    #stb.markers <- subset(placenta.markers.stb.perc_25, cluster == i)

    
    stb.markers.pos <- subset(stb.markers, avg_logFC >= 0)
    stb.markers.neg <- subset(stb.markers, avg_logFC < 0)
    
    
    x.pos <-  stb.markers.pos[,c('avg_logFC')]
    y.pos <-  -log10(stb.markers.pos[,c('p_val_adj')])
    
    x.neg <-  stb.markers.neg[,c('avg_logFC')]
    y.neg <-  -log10(stb.markers.neg[,c('p_val_adj')])
    
    
    stb.markers.pos$type <- "positive"
    stb.markers.neg$type <- "negative"
    stb.markers.combine <- rbind.data.frame(stb.markers.pos,stb.markers.neg)
    stb.markers.combine$type <- factor(stb.markers.combine$type, levels = c("positive","negative"))
    
    stb.topn.positive <- stb.markers.pos %>% dplyr::top_n(n = 50, wt = avg_logFC)
    stb.topn.negative <- stb.markers.neg %>% dplyr::top_n(n = -50, wt = avg_logFC)

    stb.topn.positive <- head(stb.topn.positive,n=10) #to label
    stb.topn.negative <- head(stb.topn.negative,n=10)
    
    ##output for GO analysis
    stb.top100.gene.positive <- (stb.markers.pos %>% dplyr::top_n(n = 100, wt = avg_logFC))$gene
    stb.top100.gene.negative <- (stb.markers.neg %>% dplyr::top_n(n = -100, wt = avg_logFC))$gene

    stb.top1000.gene.positive <- (stb.markers.pos %>% dplyr::top_n(n = 1000, wt = avg_logFC))$gene
    stb.top1000.gene.negative <- (stb.markers.neg %>% dplyr::top_n(n = -1000, wt = avg_logFC))$gene

     #top 100 gene
     write.table(stb.top100.gene.positive,
             file=paste("DEGs_extend/de_gene_split/",
                        sample,"/cluster.",
                        i,
                        ".de.genes.top100.positive.txt"
                        ,sep=""),
             sep = '\t',
             row.names = FALSE,
             col.names = FALSE,
             quote = FALSE
            )    
    write.table(stb.top100.gene.negative,
             file=paste("DEGs_extend/de_gene_split/",
                        sample,"/cluster.",
                        i,
                        ".de.genes.top100.negative.txt"
                        ,sep=""),
             sep = '\t',
             row.names = FALSE,
             col.names = FALSE,
             quote = FALSE
            )    
    
    #top 1000 gene
    write.table(stb.top1000.gene.positive,
             file=paste("DEGs_extend/de_gene_split/",
                        sample,"/cluster.",
                        i,
                        ".de.genes.top1000.positive.txt"
                        ,sep=""),
             sep = '\t',
             row.names = FALSE,
             col.names = FALSE,
             quote = FALSE
            )    
    write.table(stb.top1000.gene.negative,
             file=paste("DEGs_extend/de_gene_split/",
                        sample,"/cluster.",
                        i,
                        ".de.genes.top1000.negative.txt"
                        ,sep=""),
             sep = '\t',
             row.names = FALSE,
             col.names = FALSE,
             quote = FALSE
            )    
    
    ##output genelist with foldchange for GSEA?
    stb.top100.df.positive <- (stb.markers.pos %>% dplyr::top_n(n = 100, wt = avg_logFC))
    stb.top100.df.negative <- (stb.markers.neg %>% dplyr::top_n(n = -100, wt = avg_logFC))

    stb.top1000.df.positive <- (stb.markers.pos %>% dplyr::top_n(n = 1000, wt = avg_logFC))
    stb.top1000.df.negative <- (stb.markers.neg %>% dplyr::top_n(n = -1000, wt = avg_logFC))

    #top 100 gene df
    write.table(stb.top100.df.positive,
         file=paste("DEGs_extend/de_gene_split/",
                    sample,"/cluster.",
                    i,
                    ".de.genes.top100.df.positive.txt"
                    ,sep=""),
         sep = '\t',
         row.names = FALSE,
         col.names = TRUE,
         quote = FALSE
        )    
    write.table(stb.top100.df.negative,
         file=paste("DEGs_extend/de_gene_split/",
                    sample,"/cluster.",
                    i,
                    ".de.genes.top100.df.negative.txt"
                    ,sep=""),
         sep = '\t',
         row.names = FALSE,
         col.names = TRUE,
         quote = FALSE
        ) 
    
    #top 1000 gene df
    write.table(stb.top1000.df.positive,
         file=paste("DEGs_extend/de_gene_split/",
                    sample,"/cluster.",
                    i,
                    ".de.genes.top1000.df.positive.txt"
                    ,sep=""),
         sep = '\t',
         row.names = FALSE,
         col.names = TRUE,
         quote = FALSE
        )    
    write.table(stb.top1000.df.negative,
         file=paste("DEGs_extend/de_gene_split/",
                    sample,"/cluster.",
                    i,
                    ".de.genes.top1000.df.negative.txt"
                    ,sep=""),
         sep = '\t',
         row.names = FALSE,
         col.names = TRUE,
         quote = FALSE
        ) 
    
    
    #stb.topn.positive <- stb.markers.pos %>% dplyr::slice_max(avg_logFC,n = 10)
    #stb.topn.negative <- stb.markers.neg %>% dplyr::slice_min(avg_logFC,n = 10)
    
    
    color_pos <- seaborn_pair[paste('c',i,'_pos',sep='')]
    color_neg <- seaborn_pair[paste('c',i,'_neg',sep='')]

    names(color_pos) <- 'positive'
    names(color_neg) <- 'negative'
    
    color_use <- c(color_pos, color_neg)
    
    ##ggrepel method
    options(repr.plot.height = 5.5, repr.plot.width = 4.5)
    res.p <- ggplot( stb.markers.combine, aes(x = avg_logFC, y = -log10(p_val_adj), col= type ) ) +
              geom_point(size = 0.8) +
              geom_point(data = stb.topn.positive,size = 0.9, col = 'black', shape =1) +
              geom_point(data = stb.topn.negative,size = 0.9, col = 'black', shape = 1) +
              scale_color_manual( values = color_use  ) +
              geom_text_repel(data=stb.topn.positive,
                              aes(x = avg_logFC, y = -log10(p_val_adj), label = gene, col = type),
                              size=2.5,color = 'black',max.overlaps = 25) +
              geom_text_repel(data=stb.topn.negative,
                              aes(x = avg_logFC, y = -log10(p_val_adj), label = gene, col = type),
                              size=2.5, color = 'black',max.overlaps = 25
                             ) +
              ggtitle(paste("DEG genes of cluster ",i,sep='')) +
              xlab('logFC') +
              ylab('-log10(FDR)') +
              xlim(c(-3,3)) +
              #ylim(c(0,300)) +
              ylim(c(0,350)) +
              theme(legend.position = 'none',
                    #axis.text=element_blank(), 
                    #axis.title = element_text(size = 15, face = "bold"),
                    #axis.ticks = element_blank(),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(color="black", fill = NA,size=.25),
                    plot.title = element_text(size = 15, face = "bold"),
                    #complete = TRUE
                    plot.margin = unit(c(1,1,1,1), "lines")
                    #strip.text = element_text(size = 20), #for facet strip label
                    #strip.background = element_rect(fill='white') #the strip label fill
                   )
    print(res.p)
    res.de.list[[i]] <- res.p
    ggsave(paste("DEGs_extend/pdfs/",sample,"/cluster.",i,".de.genes.volcanoplot.pdf",sep=""), width = 4.5, height = 5.5, useDingbats = FALSE )
    
    
    ##ggsave(paste("DEGs_extend/pdfs/cluster.",i,".de.genes.volcanoplot.pdf",sep=""), width = 4.5, height = 5.5, useDingbats = FALSE )
    #ggsave(paste("DEGs_extend/pdfs/cluster.",i,".de.genes.volcanoplot.pct_25.pdf",sep=""), width = 4.5, height = 5.5, useDingbats = FALSE )
    ##ggsave(paste("DEGs_extend/pdfs/cluster.",i,".de.genes.volcanoplot.stb.pdf",sep=""), width = 4.5, height = 5.5, useDingbats = FALSE )


#     ##simple method
#     plot(stb.markers[,c('avg_logFC')],-log10(stb.markers[,c('p_val_adj')]), xlim= c(-3,3), ylim= c(0,300), main = paste("DEG genes of cluster ",i,sep=''), cex = 0.8,pch = 19, col = 'grey', xlab = "logFC", ylab = "-log10(FDR)", cex.lab = 1.5) #col = cellcolor_darker[i] )
#     points(x.pos,y.pos, cex = 0.8,pch = 19, col = seaborn_pair[paste('c',i,'_pos',sep='')] )
#     points(x.neg,y.neg, cex = 0.8,pch = 19, col = seaborn_pair[paste('c',i,'_neg',sep='')] )
    
#     points(stb.topn.positive[,c('avg_logFC')],
#            -log10(stb.topn.positive[,c('p_val_adj')]),
#            col="black",pch=1,cex=1)
#     points(stb.topn.negative[,c('avg_logFC')],
#            -log10(stb.topn.negative[,c('p_val_adj')]),
#            col="black",pch=1,cex=1)
    
#   ##add text info
#   ##text(x=14,y=15.8,pos=4,srt=45,labels = paste("FC =",format(2^threshold.FC,digits = 2),sep = " ")) 
#   ##text(x=15.3,y=14.3,pos=4,srt=45,labels = paste("FC =",format(-1*2^threshold.FC,digits = 2),sep = " ")) 

#   for(i in 1:nrow(stb.topn.positive)){
#     lines(c(stb.topn.positive[i,"avg_logFC"],
#             stb.topn.positive[i,"avg_logFC"]+.2),
#           c(-log10(stb.topn.positive[i,"p_val_adj"]),
#             -log10(stb.topn.positive[i,"p_val_adj"])+20),
#           lwd=1,col='grey')
#     text(x=stb.topn.positive[i,"avg_logFC"]+.1,
#          y=-log10(stb.topn.positive[i,"p_val_adj"])+20,
#          pos=4,labels=stb.topn.positive$gene[i],cex = .8 )
#   }
#   for(i in 1:nrow(stb.topn.negative)){
#     lines(c(stb.topn.negative[i,"avg_logFC"],
#             stb.topn.negative[i,"avg_logFC"]-.5),
#           c(-log10(stb.topn.negative[i,"p_val_adj"]),
#             -log10(stb.topn.negative[i,"p_val_adj"])+20),
#           lwd=1,col='grey')
#     text(x=stb.topn.negative[i,"avg_logFC"]-.5,
#          y=-log10(stb.topn.negative[i,"p_val_adj"])+20,
#          pos=2,labels=stb.topn.negative$gene[i],cex = .8 )
#   }
    

    
}


##combine to pdf
options(repr.plot.width=12,repr.plot.height=10)
res.de.list[['11']] + res.de.list[['8']]  + res.de.list[['2']]  + res.de.list[['10']] + res.de.list[['1']] + res.de.list[['3']] + plot_layout(ncol=3,nrow=2 ) + 
# res.de.list[[1]] + res.de.list[[2]]  + res.de.list[[3]]  + res.de.list[[4]] + res.de.list[[5]] + res.de.list[[6]] + plot_layout(ncol=3,nrow=2 ) + 
plot_annotation(title = paste('de gene volcanoplot( ',sample,")", sep =''),
                tag_levels = 'A' , 
                theme=theme(plot.title = element_text(size = 25, hjust = 0.5),
                            plot.tag = element_text(size = 8, face = 2)
                           ) 
               )

ggsave(paste("DEGs_extend/pdfs/",sample,"/cluster.all.de.genes.volcanoplot.part1.pdf",sep=""),height=10,width=12,useDingbats = FALSE)


options(repr.plot.width=12,repr.plot.height=10)
res.de.list[['4']] + res.de.list[['5']] + res.de.list[['7']]  + res.de.list[['9']]  + res.de.list[['6']] + plot_layout(ncol=3,nrow=2 ) + 

#res.de.list[[7]] + res.de.list[[8]]  + res.de.list[[9]]  + res.de.list[[10]] + plot_layout(ncol=3,nrow=2 ) + 
plot_annotation(title = paste('de gene volcanoplot( ',sample,")", sep =''),
                tag_levels = 'A' , 
                theme=theme(plot.title = element_text(size = 25, hjust = 0.5),
                            plot.tag = element_text(size = 8, face = 2)
                           ) 
               )

ggsave(paste("DEGs_extend/pdfs/",sample,"/cluster.all.de.genes.volcanoplot.part2.pdf",sep=""),height=10,width=12,useDingbats = FALSE)





####doDEG script done####





###########iteratively FindMarkers for all-to-all pairwise STB (only) in-depth de genes identification########
##or all use naive as control (ident.2)??
##include both positive and negative ones
##visualize with vocano plot arranged in grid


idents <-  c('8','5','9','10','6','1','3','2','4','7')
#idents.stb <-  c('6','1','3','4','7','2')

pairs_vs_control <- list(
     '1' = '6',
     '3' = '6',
     '4' = '6',
     '7' = '6',
     '2' = '6',
     '6' = 'all' #all other as background
)


##start to do iterative de test and visualization
res.pairs_vs_control <- list()
for(i in c('3') ){
#for(i in names(pairs_vs_control) ){
    cat (i,'_vs_',pairs_vs_control[[i]],'\n',sep="")
    id1 <- i
    id2 <- pairs_vs_control[[i]]
    stb.markers <- FindMarkers(placenta, ident.1 = id1, ident.2 = id2, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0) #slow if logfc.threshold = 0
    stb.markers <- FindMarkers(placenta, ident.1 = id1, ident.2 = id2, only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
    res.pairs_vs_control[[i]] <- stb.markers


    data.mat <- exprMat.z.aggre[,c(id1,id2)]
    #stopifnot(colnames(data.mat))
    
    
    ##scatter plot and output de gene table
    sample <- paste(i,'_vs_',pairs_vs_control[[i]],sep='')
    
    #log2thresh_exp <- 1 # log2TPM >= 1, TPM >= 2
    #log2thresh_exp <- 2.32 # log2TPM >= 2.32, TPM >= 5
    #logFC_thresh <- 0.678 # FC >= 1.6
    #logFC_thresh <- 0.5849 # FC >= 1.5
    ##logFC_thresh <- 1 # FC >= 2
    #logFC_thresh <- 0.848 #FC 1.8

    log2thresh_exp <- 0 # log2CPM >= 0, CPM >= 1
    logFC_thresh <- 0.25 # FC >= 1.18

#     select.up <- ( data.combine[,1] >= log2thresh_exp & data.combine[,2] >= log2thresh_exp  & 
#                      data.combine[,2]-data.combine[,1] >= logFC_thresh) #log2TPM >=1, log2FC >= 0.5849, 1.5 fold change
#     select.down <- (  data.combine[,1] >= log2thresh_exp & data.combine[,2] >= log2thresh_exp & 
#                         data.combine[,2]-data.combine[,1] <= -1*logFC_thresh) ##log2TPM >=1, log2FC >= 0.5849, 1.5 fold change

    ##stb.top50.positive <- stb.markers %>% dplyr::top_n(n = 50, wt = avg_logFC) #will keep original order
    ##stb.top50.negative <- stb.markers %>% dplyr::top_n(n = 50, wt = -avg_logFC)
    
    stb.top50.positive <- stb.markers %>% dplyr::slice_max(avg_logFC,n = 100)
    stb.top50.negative <- stb.markers %>% dplyr::slice_min(avg_logFC,n = 100)
    
#     stb.top50.positive.my <- head(stb.markers[order(stb.markers$avg_logFC,decreasing = TRUE),], n = 50 )
#     stb.top50.positive.slicemax <- stb.markers %>% dplyr::slice_max(avg_logFC,n = 50)
    
#     stb.top50.negative.my <- head(stb.markers[order(stb.markers$avg_logFC,decreasing = FALSE),], n = 50 )
#     stb.top50.negative.slicemin <- stb.markers %>% dplyr::slice_min(avg_logFC,n = 50)
    
#     all.equal(stb.top50.positive.my,stb.top50.positive.slicemax) #TRUE
#     all.equal( sort(rownames(stb.top50.positive)), sort(rownames(stb.top50.positive.slicemax)) ) #TRUE
    
#     all.equal(stb.top50.negative.my,stb.top50.negative.slicemin) #TRUE
#     all.equal( sort(rownames(stb.top50.negative)), sort(rownames(stb.top50.negative.slicemin)) ) #TRUE
    

    select.up <- rownames(data.mat) %in%  rownames(stb.top50.positive) 
    select.down <- rownames(data.mat) %in%  rownames(stb.top50.negative) 
    
    plot(data.mat,pch = 19, cex = 0.8,col='grey')
    points(data.mat[select.up,],pch = 19, cex = 0.8,col ='red')
    points(data.mat[select.down,],pch = 19, cex = 0.8,col ='blue')
    
    outfile = paste("DEGs_extend/pdfs/select.up.down.mat_FC1.5/out.pairwise.scatter.all.FC1.5.nolabel.",sample,".pdf",sep='')
    pdf(file = outfile,height = 7.5,width = 7.5,useDingbats = FALSE  )
    drawSimpleDE(data.mat,sample,select.up = select.up,select.down = select.down, 
               threshold.exp = log2thresh_exp, threshold.FC = logFC_thresh)

    dev.off()  


    ##volcano plot and output de gene table
    outfile = paste("select.up.down.mat_FC1.5/out.pairwise.scatter.all.FC1.5.nolabel.",sample,".pdf",sep='')
    pdf(file = outfile,height = 7.5,width = 7.5,useDingbats = FALSE  )
    drawVolcanoDE(res.table.new,sample,select.up,select.down,FC.cutoff,padj.cutoff )
    dev.off()  
    
    select.up <- rownames(stb.markers) %in%  rownames(stb.top50.positive) 
    select.down <- rownames(stb.markers) %in%  rownames(stb.top50.negative) 
    
    
    plot(stb.markers[,c('avg_logFC')],-log10(stb.markers[,c('p_val_adj')]), xlim= c(-2,2), main = sample, cex = 0.8, col = 'grey',pch = 19 )
    abline(h=30,v=c(-0.5,0.5),lty=2)
    points(stb.markers[select.up,c('avg_logFC')],-log10(stb.markers[select.up,c('p_val_adj')]),pch = 19, cex = 0.8,col ='red')
    points(stb.markers[select.down,c('avg_logFC')],-log10(stb.markers[select.down,c('p_val_adj')]),pch = 19, cex = 0.8,col ='blue')
    

    
    
}




##########################do the simple pairwise fold change visulization and output step
#height=4.166667,width = 4.177083
drawSimpleDE <- function(data  = NULL,sample = NULL,select.up = NULL, select.down = NULL, threshold.exp = NULL, threshold.FC = NULL){
  r <- format(cor(data[,1],data[,2],method = 'spearman'),digits = 3)
  #plot(data,pch=20,cex=0.65,xlim=c(-1,1000),ylim=c(-1,1000),col='lightgrey',
      #xlab=paste("WT",sample,"RPM",sep=" "),ylab=paste("KO",sample,"RPM",sep=" "))
  plot(data,pch=20,cex=0.65,xlim=c(-1,1),ylim=c(-1,1),col='lightgrey',
       xlab=paste("WT",sample,"log2CPM",sep=" "),ylab=paste("KO",sample,"log2CPM",sep=" "))
  points(data[select.up,],pch=20,cex=0.65,col="red" )
  points(data[select.down,],pch=20,cex=0.65,col="blue" )
  text(x = -1, y = 12,labels = paste("Up =",sum(select.up),sep = ""), col='black', pos = 4) #align left end
  text(x = 5.5, y = 1,labels = paste("Down =",sum(select.down),sep = ""), col='black', pos = 4) 
  text(x = -1, y = 17.5,labels = paste("R=",r,sep = ""), col='black', pos = 4)
  text(x = -1, y = 19,labels = paste("N =",nrow(data),sep = ""), col='black', pos = 4)
  text(x = -1, y = 16,labels = paste("Log2CPM >=",threshold.exp, "FC >=",format(2^threshold.FC,digits = 2),sep = " "),col='black', pos = 4)
  #text(x = 0, y = 500,labels = paste("n=",sum(select.up),sep = ""), col='black')
  #text(x = 500, y = 1,labels = paste("n=",sum(select.down),sep = ""), col='black')
  #text(x = -1, y = 18,labels = paste("R=",r,sep = ""), col='black')
  #abline(a=0.5849,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=-0.5849,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
  #abline(a=-1,b=1,col='darkgrey',lty=2,lwd=1.2,untf = FALSE)
  abline(a=threshold.FC,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  abline(a=-1*threshold.FC,b=1,col='black',lty=2,lwd=1.2,untf = FALSE)
  ##abline(a=0,b=1,col='darkgrey',lty=1,lwd=1.2)
  #abline(h=0,v=0)
  
  ##add some gene point annotation
  
  ##expected.up.new <- read.table("expected.up.all",stringsAsFactors = F)$V1
  ##expected.down.new <- read.table("expected.down.all",stringsAsFactors = F)$V1
  
  ##expected.up.remove <- read.table("expected.up.remove",stringsAsFactors = F)$V1
  ##expected.down.remove <- read.table("expected.down.remove",stringsAsFactors = F)$V1
 
  ##expected.up.new <- expected.up.new[!(expected.up.new %in% expected.up.remove  )]
  ##expected.down.new <- expected.down.new[!(expected.down.new %in% expected.down.remove  )]
  
  #expected.up.new <- read.table("expected.up.new",stringsAsFactors = F)$V1
  #expected.down.new <- read.table("expected.down.new",stringsAsFactors = F)$V1
  
  ##expected.up.new.data  <- data[match (expected.up.new,rownames(data)), ]
  ##expected.down.new.data  <- data[match (expected.down.new,rownames(data)), ]
  
  ##points(expected.up.new.data,col="yellow",pch=1,cex=0.65)
  ##points(expected.down.new.data,col="darkgreen",pch=1,cex=0.65)
  
  ##add text info
  text(x=14,y=15.8,pos=4,srt=45,labels = paste("FC =",format(2^threshold.FC,digits = 2),sep = " ")) 
  text(x=15.3,y=14.3,pos=4,srt=45,labels = paste("FC =",format(-1*2^threshold.FC,digits = 2),sep = " ")) 
#   for(i in 1:nrow(expected.up.new.data)){
#     lines(c(expected.up.new.data[i,1],expected.up.new.data[i,1]+3),c(expected.up.new.data[i,2],expected.up.new.data[i,2]+1),lwd=1,col='black')
#     text(x=expected.up.new.data[i,1]+3,y=expected.up.new.data[i,2]+1,pos=4,labels=row.names(expected.up.new.data)[i],cex = 0.8 )
#   }
#   for(i in 1:nrow(expected.down.new.data)){
#     lines(c(expected.down.new.data[i,1],expected.down.new.data[i,1]-1),c(expected.down.new.data[i,2],expected.down.new.data[i,2]+3),lwd=1,col='black')
#     text(x=expected.down.new.data[i,1]-0.5,y=expected.down.new.data[i,2]+3,pos=2,labels=row.names(expected.down.new.data)[i],cex = 0.8 )
#   }
  
  return(1)
}
###########




##############valcano plot DGE results

drawVolcanoDE <- function(data  = NULL,sample = NULL,select.up = NULL, select.down = NULL, 
                          FC.cutoff = NULL, padj.cutoff = NULL){
  #input DGE result table
  plot(data$logFC,-log10(data$padj),pch=20,cex = 0.5,xlim = c(-10,10),col='black',
       main=paste(sample,"\n(FDR < ",padj.cutoff,", FC > ",FC.cutoff,")",sep=""),xlab = expression(Log[2]~(Fold~Change)), ylab=expression(-Log[10]~(FDR)),
       cex.lab=1.5)
       #cex.lab=1.5,family='Calibri Light')
  points(data[select.up,"logFC"],-log10(data[select.up,"padj"]),col='red' ,pch=20,cex = 0.5)
  points(data[select.down,"logFC"],-log10(data[select.down,"padj"]),col='navy' ,pch=20,cex = 0.5)
  abline(v=c(-log2(FC.cutoff),log2(FC.cutoff)),lty=2,col='grey')
  abline(h=-log10(padj.cutoff),lty=2,col='grey')
  y_max <- max(-log10(data$padj))
  text(x = 11, y = y_max,labels = paste("Up =",sum(select.up),sep = ""), col='red', pos = 2) #align right end
  text(x = -11, y = y_max,labels = paste("Down =",sum(select.down),sep = ""), col='navy', pos = 4) #align left end
  text(x=log2(FC.cutoff)+1.2,y=7-0.1,pos=2,srt=90,labels = paste("FC = ",FC.cutoff,sep = " ")) 
  text(x=-log2(FC.cutoff),y=7,pos=2,srt=90,labels = paste("FC = ",-FC.cutoff,sep = " "))
  
  ##add some gene point annotation
  
  # expected.up.new <- read.table("expected.up.all",stringsAsFactors = F)$V1
  # expected.down.new <- read.table("expected.down.all",stringsAsFactors = F)$V1
  # 
  # expected.up.remove <- read.table("expected.up.remove",stringsAsFactors = F)$V1
  # expected.down.remove <- read.table("expected.down.remove",stringsAsFactors = F)$V1
  # 
  # expected.up.new <- expected.up.new[!(expected.up.new %in% expected.up.remove  )]
  # expected.down.new <- expected.down.new[!(expected.down.new %in% expected.down.remove  )]
  # 
  expected.up.new <- read.table("expected.up.new.new",stringsAsFactors = F)$V1
  expected.down.new <- read.table("expected.down.new.new",stringsAsFactors = F)$V1
  
  expected.up.new.data  <- data[match (expected.up.new,rownames(data)), ]
  expected.down.new.data  <- data[match (expected.down.new,rownames(data)), ]
  
  points(expected.up.new.data$logFC,-log10(expected.up.new.data$padj),col="yellow",pch=1,cex=0.65)
  points(expected.down.new.data$logFC,-log10(expected.down.new.data$padj),col="darkgreen",pch=1,cex=0.65)
  
  ##add text info
  ###text(x=14,y=15.8,pos=4,srt=45,labels = paste("FC =",format(2^threshold.FC,digits = 2),sep = " ")) 
  ###text(x=15.3,y=14.3,pos=4,srt=45,labels = paste("FC =",format(-1*2^threshold.FC,digits = 2),sep = " ")) 

  # for(i in 1:nrow(expected.down.new.data)){
  #   lines(c(expected.down.new.data[i,"logFC"],expected.down.new.data[i,"logFC"]+3),c(-log10(expected.down.new.data[i,"padj"]),-log10(expected.down.new.data[i,"padj"])+1),lwd=1,col='grey')
  #   text(x=expected.down.new.data[i,"logFC"]+3,y=-log10(expected.down.new.data[i,"padj"])+1,pos=4,labels=row.names(expected.down.new.data)[i],cex = 1.5 )
  # }
  # for(i in 1:nrow(expected.up.new.data)){
  #   lines(c(expected.up.new.data[i,"logFC"],expected.up.new.data[i,"logFC"]-3),c(-log10(expected.up.new.data[i,"padj"]),-log10(expected.up.new.data[i,"padj"])+3),lwd=1,col='grey')
  #   text(x=expected.up.new.data[i,"logFC"]-2.5,y=-log10(expected.up.new.data[i,"padj"])+3,pos=2,labels=row.names(expected.up.new.data)[i],cex = 1.5 )
  # }

  return(1)
}
















