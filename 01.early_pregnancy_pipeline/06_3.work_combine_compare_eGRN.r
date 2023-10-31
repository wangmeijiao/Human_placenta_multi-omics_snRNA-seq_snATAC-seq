

##because i-cisTarget and pycistarget use different motif ranking database (feather file, with motif as row and predefine meta-region set as column), perhaps it is reasonable to combine two results?


library(pheatmap)
library(ComplexHeatmap)

library(ggplot2)


#########set colorset
#color_gradient_biasedBR  <- readRDS('../../../../../02.snapATAC_harmony/chromVAR_TF_specific/full12526/color_gradient_biasedBR.rds')
color_gradient_biasedBR  <- readRDS('/home/mjwang/progs/misc-tools/r/color_gradient_biasedBR.rds')
#plot
barplot(1:length(color_gradient_biasedBR),col=rev(color_gradient_biasedBR) )


#color_set_yellowbrick.flat <- readRDS('../../../../../03.snRNA_snATAC/liger/color_set_yellowbrick.flat.rds') 
#saveRDS(color_set_yellowbrick.flat,'complexHeatmap_interactive/color_set_yellowbrick.flat.rds')
color_set_yellowbrick.flat <- readRDS('/home/mjwang/progs/misc-tools/r/color_set_yellowbrick.flat.rds')

colorset1 <- color_set_yellowbrick.flat[['RdYlGn.11']]
colorset2 <- rev(color_set_yellowbrick.flat[['YlGnBu.9']]) #or this
colorset3 <- rev(color_set_yellowbrick.flat[['RdPu.9']]) #use this !
colorset4 <- color_set_yellowbrick.flat[['Spectral.8']]
colorset5 <- color_set_yellowbrick.flat[['RdGy.11']]

#colorset <- colorset1

colorset <- c('white',color_set_yellowbrick.flat[['Reds.6']])#for freq heatmap

barplot(1:length(colorset),col=colorset )

colorset_go <- rev(color_set_yellowbrick.flat[['YlGnBu.9']])
colorset_pathway <- rev(color_set_yellowbrick.flat[['RdPu.9']])

barplot(1:length(colorset_go),col=colorset_go )
barplot(1:length(colorset_pathway),col=colorset_pathway )

##color set done


map_cluster <- list( #atac cid to rna cid
    'c1' = 'c9',
    'c7' = 'c6',
    'c9' = 'c11',
    'c2' = 'c8',
    'c8' = 'c1',
    'c3' = 'c3',
    'c6' = 'c4',
    'c4' = 'c2',
    'c5' = 'c10'
    #- ' = '5
    
#        'c6' = 'c9',
#        'c3' = 'c5',
#        'c9' = 'c10',
#        'c5' = 'c6',
#        'c1' = 'c7',
#        'c2' = 'c3',
#        'c4' = 'c2' 
)




#########1 start from tfmotif.list.filter: motif was choosen for max nes, include peak, target gene, peak_combined, targetgene_combined

#i-cisTarget dar_filter_p2g_auc0.005

tfmotif.list.filter.icistarget <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/cisTarget/i-cisTarget/dar_filter_p2g_auc0.005/tfmotif.list.filter.rds')

#pycistarget result_pycisTarget_auc0.005, with dar full 
tfmotif.list.filter.pycistarget <- readRDS('../pycisTarget/result_pycisTarget_auc0.005/tfmotif.list.filter.rds')

##2 compare and combine record of each tf in each cluster (only for share tf record?)

#use p2g overlap dar
lapply(tfmotif.list.filter.icistarget, function(x){length(x)} )
$c1 58 $c7 66 $c9 61 $c2 94 $c8 154 $c3 125 $c6 85 $c4 119 $c5 151


#use dar full
lapply(tfmotif.list.filter.pycistarget, function(x){length(x)} )
$c1 114 $c7 115 $c9 99 $c2 55 $c8 56 $c3 96 $c6 93 $c4 77 $c5 130

all.equal(names(tfmotif.list.filter.icistarget),names(tfmotif.list.filter.pycistarget) ) #TRUE


tfmotif.list.filter.combine <- list()

tf_highconf <- list()

for(cid in c('c1','c7','c9','c2','c8','c3','c6','c4','c5')){
    cat('process cid ',cid,'\n',sep='')
    tf_share <- intersect(names(tfmotif.list.filter.icistarget[[cid]]),names(tfmotif.list.filter.pycistarget[[cid]]) )
    tf_highconf[[cid]] <- tf_share
    tf_combine <- unique(c(names(tfmotif.list.filter.icistarget[[cid]]),names(tfmotif.list.filter.pycistarget[[cid]])))
    
    cat('shared tf: ',paste(tf_share,collapse = ','),'\n',sep='')
    cat('combined tf: ',paste(tf_combine,collapse = ','),'\n',sep='')
    
    for(tfid in tf_combine){
        cat('  for tf ',tfid,'\n',sep='')
        tfmotif.list <- list()
        
        if(tfid %in% tf_share){# merge the two list for all keys except heatmap (but fill heatmap with 'NA')
            tf.icistarget <- tfmotif.list.filter.icistarget[[cid]][[tfid]]
            tf.pycistarget <- tfmotif.list.filter.pycistarget[[cid]][[tfid]]
            
            stopifnot(all.equal( names(tf.icistarget), names(tf.pycistarget) ) )    #'cid''tf''NES''motif''peaks''nearGenes''peaks_combined''nearGenes_combined''targetGenes_p2g''targetGenes_p2g_combined''shared_targets''heatmap'
            
            #check
            stopifnot(all.equal( tf.icistarget[['cid']], tf.pycistarget[['cid']] ) ) 
            stopifnot(all.equal( tf.icistarget[['tf']], tf.pycistarget[['tf']] ) ) 
            tfmotif.list[['cid']] <- tf.icistarget[['cid']]
            tfmotif.list[['tf']] <- tf.icistarget[['tf']]
            
            #max NES (or use pycistarget NES only?)
            tfmotif.list[['NES']] <- max(as.numeric(tf.icistarget[['NES']]),as.numeric(tf.pycistarget[['NES']]))
            
            #combine motif id
            tfmotif.list[['motif']] <- paste(c(tf.icistarget[['motif']], tf.pycistarget[['motif']] ),collapse = ',')
            
            #check ,combine peaks and remove duplicate
            checkCombine <- function(keyid = NULL, check_ratio = FALSE){
                
                stopifnot( keyid %in% names(tf.icistarget) )
                stopifnot( keyid %in% names(tf.pycistarget) )
                
                peaks_share <- intersect(tf.icistarget[[keyid]], tf.pycistarget[[keyid]]) 
                #stopifnot( length(peaks_share) != 0 )
                
                flag <- TRUE
                
                if(length(peaks_share) == 0 ){ 
                    cat('    Error: 0 overlap of icistarget and pycistarget result of key ',keyid, ' (still combine)\n',sep=''); 
                    peaks <- unique(c(tf.icistarget[[keyid]], tf.pycistarget[[keyid]]) ) 

                    if(check_ratio != FALSE){
                                             return(list('peaks'=peaks,'flag'=flag))
                                            }else{return(peaks)}
                }
                
                
                if(check_ratio != FALSE){
                      
                  if( length(peaks_share)/length(tf.icistarget[[keyid]]) < check_ratio   ){ cat('    Warnning: less than ',check_ratio,'(',round(length(peaks_share)/length(tf.icistarget[[keyid]]),2),')',' overlap in icistarget of key ',keyid, '(still combine)\n',sep='') ;flag <- FALSE}
                  if( length(peaks_share)/length(tf.pycistarget[[keyid]]) < check_ratio   ){ cat('    Warnning: less than ',check_ratio,'(',round(length(peaks_share)/length(tf.pycistarget[[keyid]]),2),')',' overlap in pycistarget of key ',keyid, '(still combine)\n',sep='') ;flag <- FALSE}
                }
                
                peaks <- unique(c(tf.icistarget[[keyid]], tf.pycistarget[[keyid]]) ) #order lost
                
                if(check_ratio != FALSE){
                                         return(list('peaks'=peaks,'flag'=flag))
                                        }else{return(peaks)}
                
                
            }
            
            
            #combine but ordre lost
            tfmotif.list[['peaks']] <- checkCombine(keyid = 'peaks', check_ratio = FALSE)
            tfmotif.list[['nearGenes']] <- checkCombine(keyid = 'nearGenes', check_ratio = FALSE)
            tfmotif.list[['peaks_combined']] <- checkCombine(keyid = 'peaks_combined', check_ratio = FALSE)
            tfmotif.list[['nearGenes_combined']] <- checkCombine(keyid = 'nearGenes_combined', check_ratio = FALSE)
            
            temp1 <- checkCombine(keyid = 'targetGenes_p2g', check_ratio = 0.6)
            tfmotif.list[['targetGenes_p2g']] <- temp1[['peaks']]
            
            temp2 <- checkCombine(keyid = 'targetGenes_p2g_combined', check_ratio = 0.6)
            tfmotif.list[['targetGenes_p2g_combined']] <- temp2[['peaks']]
            if(temp2[['flag']]){ tf_highconf[[cid]] <- c(tf_highconf[[cid]], tfid)  } 
            
            
            tfmotif.list[['shared_targets']] <- checkCombine(keyid = 'shared_targets', check_ratio = FALSE)
            tfmotif.list[['heatmap']] <- NA
            
            
            tfmotif.list.filter.combine[[cid]][[tfid]] <- tfmotif.list
            
            
        }else{#just add in
            
            tf.icistarget <- tfmotif.list.filter.icistarget[[cid]][[tfid]]
            tf.pycistarget <- tfmotif.list.filter.pycistarget[[cid]][[tfid]]
            
            if(tfid %in% names(tfmotif.list.filter.icistarget[[cid]]) ){
                tfmotif.list.filter.combine[[cid]][[tfid]] <- tfmotif.list.filter.icistarget[[cid]][[tfid]]
                
            }else{ tfmotif.list.filter.combine[[cid]][[tfid]] <- tfmotif.list.filter.pycistarget[[cid]][[tfid]] }
        }
        
    }
    
    
}

saveRDS(tfmotif.list.filter.combine,'tfmotif.list.filter.combine.rds')
saveRDS(tf_highconf,'tf_highconf.rds')

str(tfmotif.list.filter.combine,max.level = 3)



##########2 format tf-nes-cluster mat among all cluster and compared to tf-expr mat######

tfmotif.list.filter <- tfmotif.list.filter.combine #for code compatiblility

tfmotif.nes.list <- list()
for(cid in names(tfmotif.list.filter) ){
    tfmotif.nes.list[[cid]] <- lapply( tfmotif.list.filter[[cid]], function(x){ 
          #grep('DUSP8',x[['shared_targets']],value=TRUE)
          x[['NES']]
    }   )
    
    tfmotif.nes.list[[cid]] <- unlist(tfmotif.nes.list[[cid]])
}



list2DF <- function(datalist = NULL){
    len.list <- length(datalist)
    allid.list <- sapply( datalist, function(x){ names(x) }  )
    allid <- Reduce(f = c, x = allid.list)
    allid <- unique(allid)
    
    nes.df <- data.frame(matrix(data = 0, nrow = length(allid), ncol = len.list ),row.names = allid )
    colnames(nes.df) <- names(datalist)
    
    for(i in names(datalist) ){
        for(j in names(datalist[[i]]) ){  
            nes.df[j,i] <- datalist[[i]][[j]]
        
        }
        
        
    }
    return(nes.df)
    
}




tfmotif.nes.df <- list2DF(datalist = tfmotif.nes.list)
rowid <- rownames(tfmotif.nes.df)
tfmotif.nes.df <- apply(tfmotif.nes.df,2,as.numeric)
rownames(tfmotif.nes.df) <- rowid
564 x 9 #combined i-cistarget and pycistarget tf
#343 x 9 
#604 x 9
#486 x 9



options(repr.plot.height = 55, repr.plot.width = 5.5)
#options(repr.plot.width = width, repr.plot.height = height)
nes.p <- pheatmap(
               tfmotif.nes.df, 
               cluster_cols = TRUE,
               cluster_rows = TRUE, 
               treeheight_col = 2.5,
               treeheight_row = 2.5,
               #cellwidth = 20, 
               #cellheight = 2.8,
               na_col = 'grey', 
               color = colorset,
               breaks =  seq( 0, max(tfmotif.nes.df,na.rm = TRUE), by = max(tfmotif.nes.df,na.rm = TRUE)/length(colorset)  )  ,
               border =TRUE,
               border_color = 'black',
               #labels_col = colid, 
#                labels_row = make_bold_names(tftarget.genesets.freqDF,
#                                             rownames, 
#                                             rowid_hi
#                                            ),

               ##angle_col = 315,
               ##number_format = "%.3f",
               ##display_numbers = TRUE,
               #display_numbers = data.df.text,
               #number_color = 'white',
               #fontsize_number = 10,
               fontsize_col = 20,
               fontsize_row = 8,
               main = paste('raw ATAC_TF_enrichment_in_dar (',nrow(tfmotif.nes.df),')',sep=''),
               silent = FALSE
               #legend_breaks = c(2,4,6,8,10),
               #legend_labels = c(2,4,6,8,10),
               #legend = TRUE
#                        height = height,
#                        width = width,
#                        filename = paste('result_do_GO_Pathway_quickonestep/ck.GO_BP.pvalue_',
#                                         pvalue,'.qvalue_', qvalue,'.heatmap.pdf' ,sep = '' 
#                                        )

              )




####readin expression mat to add TF gene expression ###############


exprMat.ave.z <-readRDS('/home/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.seurat_harmony/exprMat.ave.z.rds')


table(rownames(tfmotif.nes.df) %in% rownames(exprMat.ave.z))
FALSE  TRUE 
   46   518 

# FALSE  TRUE 
#    28   315

# FALSE  TRUE 
#    30   393

# FALSE  TRUE 
#    55   549 

# FALSE  TRUE 
#    95   686 

# FALSE  TRUE 
#    50   436


tfmotif.exp.df <-  exprMat.ave.z[rownames(tfmotif.nes.df)[rownames(tfmotif.nes.df) %in% rownames(exprMat.ave.z)] ,]

colnames(tfmotif.exp.df) <- paste('c',colnames(tfmotif.exp.df),sep='')
#518 x 11 

colnames(exprMat.ave.z) <- paste0('c',colnames(exprMat.ave.z))



rowsds <- matrixStats::rowSds(as.matrix(tfmotif.exp.df)) 
quantile(rowsds)
0% 0 25% 0.0341119380404436 50% 0.105664705399811 75% 0.20960621349174 100% 0.862695880097478

#0% 0 25% 0.0363122441095136 50% 0.116759898798046 75% 0.202081333125182 100% 0.862695880097478
#0% 9.09696665366478e-19 25% 0.035196195352777 50% 0.106537710274705 75% 0.22504298929468 100% 0.862695880097478
#0% 0 25% 0.0338533934349017 50% 0.0928136970060022 75% 0.173794584630954 100% 0.745274926860947
#0 %0 25% 0.0340040628757603 50% 0.092181068482616 75% 0.179327464130576 100% 0.745274926860947

#0% 0 25% 0.0343054017574775 50% 0.0964485621056574 75% 0.18938655327171 100% 0.745274926860947

tfmotif.exp.df.sel <- tfmotif.exp.df[rowsds >= 0.209,]
#tfmotif.exp.df.sel <- tfmotif.exp.df[rowsds >= 0.20,]
#tfmotif.exp.df.sel <- tfmotif.exp.df[rowsds >= 0.225,]
#tfmotif.exp.df.sel <- tfmotif.exp.df[rowsds >= 0.17,]
#tfmotif.exp.df.sel <- tfmotif.exp.df[rowsds >= 0.18,]
#115

options(repr.plot.height=55,repr.plot.width=6)
hp.tfexpr.all = Heatmap(tfmotif.exp.df, name = "tfmotif.exp.df", 
#options(repr.plot.height=20,repr.plot.width=6)
#hp.tfexpr.all.sel = Heatmap(tfmotif.exp.df.sel, name = "tfmotif.exp.df.sel", 
             cluster_rows = TRUE, 
             cluster_columns = TRUE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             row_title = "TFmotif expression (sd >= 0.209,hclust: ward.D)",
             clustering_method_row = "ward.D", ##ward.D,complete
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 8),
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
#width = max(grobWidth(textGrob(labels))))

##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
set.seed(1);draw(hp.tfexpr.all, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750



###

table(rownames(tfmotif.exp.df) %in% rownames(tfmotif.nes.df) )
TRUE 
 518 

# TRUE
# 315

tfmotif_nes_order <- tfmotif.nes.df[rownames(tfmotif.exp.df),c('c1','c7','c9','c2','c8','c3','c6','c4','c5')] #393 x 9 #549 x 7 #686 x 7#486 x 9
tfmotif_expr_order <-  tfmotif.exp.df[, c('c9','c6','c11','c8','c1','c3','c4','c2','c10')] #393 x 9 #549 x 7 #686 x 7 #436 x 9


all.equal(rownames(tfmotif_nes_order),rownames(tfmotif_expr_order))
#TRUE

ncol(tfmotif_nes_order) == ncol(tfmotif_expr_order)
#TRUE





############3 try to align tfmotif.nes and tfmotif.expr#########################


###1) simply filter nes mat by sd only and output aligned mat


rowsds <- matrixStats::rowSds(as.matrix(tfmotif_nes_order)) 
quantile(rowsds)
0 % 1.001534864 25% 1.3776429035599 50% 1.79490451117222 75% 2.70589030896802 100% 8.99409851497449

#0% 1.00024370764536 25% 1.56135735179872 50% 1.85810023684149 75% 2.73892451682364 100% 8.91745465096806
#0% 1.001534864 25% 1.17049302433333 50% 1.56724911366667 75% 2.29690767344758 100% 5.67053096288452
#0% 0 25% 1.24365301086292 50% 1.69424550107429 75% 2.04032305229728 100% 4.4030093758897
#0 %0 25% 1.32590404210684 50% 1.72133111753708 75% 2.1376933095494 100% 4.50416386545977
#0 %0 25% 1.23685615927494 50% 1.69746221464998 75% 2.05908753706225 100% 3.16393933673023
#0% 0.802997574237999 25% 1.14129335983333 50% 1.58681988662351 75% 2.01129544418732 100% 3.12067000400633

sd.cutoff <- 1.3776#1.5613 #1.56 #1.69 #1.72#1.69
tfmotif_nes_order.sel <- tfmotif_nes_order[rowsds >= sd.cutoff,]
388
#241
#198
#275
#353
#219
#183

options(repr.plot.height = 45, repr.plot.width = 5.5)
#options(repr.plot.width = width, repr.plot.height = height)
nes.sel.p <- pheatmap(
               tfmotif_nes_order.sel, 
               cluster_cols = TRUE,
               cluster_rows = TRUE, 
               treeheight_col = 2.5,
               treeheight_row = 2.5,
               #cellwidth = 20, 
               #cellheight = 2.8,
               na_col = 'white', 
               color = colorset,#rev(colorset_pathway),
               breaks =  seq( 0, max(tfmotif_nes_order.sel,na.rm = TRUE), by = max(tfmotif_nes_order.sel,na.rm = TRUE)/length(colorset)  )  ,
               border =TRUE,
               border_color = 'black',
               #labels_col = colid, 
#                labels_row = make_bold_names(tftarget.genesets.freqDF,
#                                             rownames, 
#                                             rowid_hi
#                                            ),

               ##angle_col = 315,
               #display_numbers = data.df.text,
               #number_color = 'white',
               #fontsize_number = 10,
               fontsize_col = 20,
               fontsize_row = 8,
               main = paste('filter sd ATAC_TF_enrichment_in_dar ( 20% sd cutoff:',sd.cutoff,' )',sep=''),
               silent = FALSE
               #legend_breaks = c(2,4,6,8,10),
               #legend_labels = c(2,4,6,8,10),
               #legend = TRUE
#                        height = height,
#                        width = width,
#                        filename = paste('result_do_GO_Pathway_quickonestep/ck.GO_BP.pvalue_',
#                                         pvalue,'.qvalue_', qvalue,'.heatmap.pdf' ,sep = '' 
#                                        )

              )


rowid <- rownames(tfmotif_nes_order.sel)[nes.sel.p$tree_row$order]
colid <- colnames(tfmotif_nes_order.sel)[nes.sel.p$tree_col$order]

table(rowid %in% rownames(tfmotif_expr_order))
#388 TRUE
#219 TRUE

tfmotif_expr_order.sel <- tfmotif_expr_order[rowid,unlist(map_cluster[colid])]
#388 x 9


options(repr.plot.height=45,repr.plot.width=5.5)
exrp.sel.p = Heatmap(tfmotif_expr_order.sel, name = 'zscore', 
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             row_title = paste('TF_expression ( aligned to TF_NES )',sep=''),
             clustering_method_row = "ward.D", ##ward.D,complete
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 8),
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
set.seed(1);draw(exrp.sel.p, heatmap_legend_side = "left",gap = unit(0.1, "cm")) 




##2) filter by correlation of tf nes mat with TF expr mat and output aligned mat###################


tfmotif_nes_order <- tfmotif.nes.df[,c('c1','c7','c9','c2','c8','c3','c6','c4','c5')] #423 x 9 #604 x 9 #486 x 9
tfmotif_expr_order <-  tfmotif.exp.df[, c('c9','c6','c11','c8','c1','c3','c4','c2','c10')] #393 x 9 #549 x 7 #436 x 9

all.equal (nrow(tfmotif_nes_order), nrow(tfmotif_expr_order) )#no 

length(intersect(rownames(tfmotif_nes_order),rownames(tfmotif_expr_order)))
#518

#315
#393
#549
#686
#436

i <- 'STAT4'
i <- 'STAT5A'
i <- 'MITF'
i <- 'CEBPB'

unlist(tfmotif_expr_order[i,])
unlist(tfmotif_nes_order[i,])

cor(unlist(tfmotif_expr_order[i,]), unlist(tfmotif_nes_order[i,]))

#0.833021795964144 STAT4
#0.813862083337557 STAT5A
#0.0128766126069388 MITF
#0.260419445911582 CEBPB


#0.833021795964144 STAT4
#0.813862083337557 STAT5A

#0.949077518228907 STAT4
#0.693993551910119


#0.819235091282893 STAT4
#0.75978 STAT5A


res.tfmotif.cor <- vector()
for(i in rownames(tfmotif_expr_order)){
  res.tfmotif.cor[i] <- cor(unlist(tfmotif_expr_order[i,]), unlist(tfmotif_nes_order[i,]))
  
}

res.tfmotif.cor <- res.tfmotif.cor[!is.na(res.tfmotif.cor)]

options(repr.plot.width = 7.5, repr.plot.height = 7.5)
plot( 1:length(res.tfmotif.cor), res.tfmotif.cor[order(res.tfmotif.cor)], pch = 19, cex = 0.2, main = 'Rank of correlation coefficency between \ntf nes and tf expr',xlab = 'TF Rank', ylab = 'Pearson correlation coefficency'  )
abline(h=0.25,lty=2,col='red')

res.tfmotif.cor.rank <- res.tfmotif.cor[order(res.tfmotif.cor,decreasing = TRUE)]
res.tfmotif.cor.rank.0.25 <- res.tfmotif.cor.rank[res.tfmotif.cor.rank > 0.25]

names(res.tfmotif.cor.rank.0.25) #top cor tf
'XRCC4''TFDP1''E2F3''SUPT20H''BRF1''ZNF32''TEAD4''RORB''NR5A1''ZNF205''RBBP5''EZH2''FIGLA''ZCCHC14''DLX5''TEAD1''TCF7L1''TP63''SMARCA4''ZFP14''CTNNB1''E2F2''ZMIZ1''MCTP2''STAT4''HIC1''STAT5A''TP53''SETBP1''ZEB1''ZSCAN29''HNF4G''TAL1''ZNF554''STAT6''IRX4''CEBPD''IRX2''TRIM21''BHLHE41''SP100''ZNF71''NELFE''TCF7''LCORL''IRX6''SALL2''USF1''NR4A1''TCF3''ESRRA''TFAP2B''TFEB''GATA3''UBB''E2F4''CEBPG''PAXIP1''NR1H4''MEF2D''SP4''MAZ''E2F1''SUZ12''NFATC1''GRHL3''PPARG''TFAP2C''ZBTB14''EGR1''MXD4''MYB''STAT5B''AKR1A1''ZBTB21''POLR2A''ATF3''GATA2''KLF10''RNF114''GLIS2''TP73''ZNF799''RELA''SMAD6''RORA''BACH1''HNF1B''HAND1''STAT3''FOSL2''PAX3''SMAD7''XBP1''ZNF284''TCF4''ZFP2''MAFK''BCL3''INSM2''TCF21''ZNF770''FOXL2''ZNF611''BNC1''NFE2''ARID5B''CNOT3''MAFG''SP3''NR2C1''FOXJ3''HES1''BACH2''FOSL1''GCM1''ZBED6''ZNF217''GCM2''MAFF''TFAP2E''HOXB2''CREB5''PPARA''MAX''SMARCB1''PGAM2''NR3C1''TCF12''DLX3''KLF4''RXRA''RAB7A''HES2''ZNF83''MYOG''ZNF664''WT1''CEBPA''CEBPB''EIF5A2'


##filter by cor
res.tfmotif.cor.sel.cor0.6 <- res.tfmotif.cor[res.tfmotif.cor >= 0.6]
res.tfmotif.cor.sel.cor0.5 <- res.tfmotif.cor[res.tfmotif.cor >= 0.5]
res.tfmotif.cor.sel.cor0.4 <- res.tfmotif.cor[res.tfmotif.cor >= 0.4]
res.tfmotif.cor.sel.cor0.25 <- res.tfmotif.cor[res.tfmotif.cor >= 0.25]

tf.cor.sel.cor0.5 <- names(res.tfmotif.cor.sel.cor0.5)
tf.cor.sel.cor0.4 <- names(res.tfmotif.cor.sel.cor0.4)
tf.cor.sel.cor0.25 <- names(res.tfmotif.cor.sel.cor0.25)

tf.cor.sel <- tf.cor.sel.cor0.25
#tf.cor.sel <- tf.cor.sel.cor0.5

'TP63''TP53''TP73''MYB''TEAD4''TEAD1''ZBTB14''NELFE''IRX4''EZH2''ESRRA''POLR2A''MCTP2''BRF1''SMARCA4''RBBP5''E2F2''E2F1''ZCCHC14''ZMIZ1''CEBPB''CNOT3''MAZ''XRCC4''ZSCAN29''SP4''PGAM2''SUPT20H''TFDP1''AKR1A1''EGR1''E2F3''SUZ12''TCF12''WT1''ZFP14''NR3C1''HIC1''DLX5''ZNF611''TFAP2B''TFAP2E''TFAP2C''PPARG''RXRA''NR1H4''HNF4G''PPARA''NR2C1''ZNF205''ZNF217''CTNNB1''TCF7''TCF7L1''BHLHE41''USF1''TFEB''TCF3''MAX''HAND1''TCF4''MXD4''TAL1''MYOG''LCORL''PAXIP1''ZNF32''IRX2''ZEB1''DLX3''STAT3''NR4A1''SP3''E2F4''SP100''ARID5B''SETBP1''UBB''ZNF71''RORA''FOSL1''ATF3''GRHL3''FOSL2''NFATC1''GATA2''GATA3''RELA''BACH2''NFE2''BCL3''BACH1''RNF114''INSM2''MAFK''MAFF''MAFG''CEBPD''CEBPG''CEBPA''HES2''RORB''SMARCB1''PAX3''FOXJ3''NR5A1''EIF5A2''GCM2''GCM1''ZNF554''RAB7A''IRX6''FOXL2''GLIS2''KLF4''BNC1''HOXB2''TRIM21''CREB5''ZNF770''XBP1''ZBTB21''STAT5A''STAT4''STAT6''STAT5B''ZBED6''FIGLA''ZNF799''ZNF83''HNF1B''TCF21''ZNF664''ZNF284''MEF2D''KLF10''ZFP2''SMAD7''SMAD6''SALL2''HES1'

all.equal(sort(names(res.tfmotif.cor.rank.0.25)),sort(tf.cor.sel) ) #TRUE


#'TP63''TP53''TP73''MYB''IRX4''TEAD1''TEAD4''ZFP14''ZBTB14''NR3C1''HIC1''DLX5''ZNF611''TFAP2B''TFAP2E''TFAP2C''PPARG''RXRA''HNF4G''THRB''RXRB''PPARA''NR1H2''NR2C1''ZNF205''TFCP2''UBP1''TFCP2L1''CTNNB1''TCF7''TCF7L2''TCF7L1''BHLHE41''USF1''TFEB''TCF3''MYC''HAND1''TCF4''MXD4''TAL1''LCORL''NR2F2''PAXIP1''ZNF32''IRX2''ZEB1''DLX3''STAT3''NR4A1''ESRRA''RARA''RARG''NR1H3''ZNF71''RORA''FOSL1''INSM2''ATF3''FOSL2''BACH2''NFE2''MAFK''BACH1''MAFF''MAFG''ZNF655''ESR2''CEBPD''CEBPG''CEBPA''CEBPB''HES2''RORB''GCM2''GCM1''HOXB2''XBP1''CREB5''ZBTB21''STAT5A''STAT4''STAT6''STAT5B''SP3''FIGLA''ZNF799''ZNF83''TCF21''GLIS2''KLF7''KLF10''KLF5''ZFP2''SMAD7''SMAD6''HES1' #0.25


#'TP63''TP73''MYB''TEAD4''TEAD1''ZBTB14''NELFE''IRX4''EZH2''EP300''ESRRA''POLR2A''MCTP2''BRF1''SMARCA4''RBBP5''E2F2''E2F1''ZCCHC14''ZMIZ1''CNOT3''MAZ''XRCC4''ZSCAN29''SP4''PGAM2''SUPT20H''TFDP1''AKR1A1''MECOM''EGR1''E2F3''SUZ12''TCF12''WT1''SP3''E2F4''SP100''ARID5B''SETBP1''UBB''GRHL1''ZNF217''FOSL1''JUNB''ATF3''GRHL3''FOSL2''NFATC1''GATA2''STAT3''GATA3''RELA''HNF4A''NFE2''BCL3''BACH1''RNF114''MAFK''CEBPD''SMARCB1''NFE2L1''MAFG''PAX3''CEBPG''PURA''FOXJ3''NR5A1''EIF5A2''MAX''RAB7A''GCM2''IRX6''FOXL2''GLIS3''KLF9''KLF4''BNC1''TRIM21''ZNF554''BATF3''ZNF770''NR1H4''ZBED6''ZNF664''ZNF284''STAT6''STAT4''STAT5B''STAT5A''MEF2D''SALL2''MITF''HNF1B''MYF6''MYOG' #0.25

#'TP63''TEAD4''TEAD1''NELFE''IRX4''EZH2''ESRRA''POLR2A''MCTP2''BRF1''SMARCA4''RBBP5''E2F2''E2F1''ZCCHC14''ZMIZ1''MAZ''XRCC4''ZSCAN29''SP4''SUPT20H''TFDP1''AKR1A1''MECOM''EGR1''E2F3''SUZ12''TCF12''E2F4''SP100''SETBP1''UBB''ZNF217''ATF3''GRHL3''FOSL2''NFATC1''GATA2''GATA3''RELA''BACH1''RNF114''CEBPD''NR5A1''GCM2''IRX6''GLIS3''TRIM21''ZNF554''STAT6''STAT4''STAT5A''MEF2D''SALL2''MITF''HNF1B''MYF6''MYOG' #0.5






##manually add these tf in case of filtered by low cor


c('CEBPB','ESRRG','PPARD','REL','MYCN','FOSL2','TBX3','BCL6','GCM1','OVOL1','POU2F3','MITF') %in% tf.cor.sel

tf.cor.sel <- c(tf.cor.sel,'CEBPB','ESRRG','PPARD','REL','MYCN','FOSL2','TBX3','BCL6','GCM1','OVOL1','POU2F3','MITF')

tf.cor.sel[duplicated(tf.cor.sel)]
'CEBPB''FOSL2''GCM1'

#FOSL2
#no dup
tf.cor.sel <- tf.cor.sel[!duplicated(tf.cor.sel)]

tf.cor.sel[!tf.cor.sel %in% rownames(tfmotif_nes_order) | !tf.cor.sel %in% rownames(tfmotif_expr_order)]
'ESRRG''REL''TBX3''POU2F3'

#'ESRRG''REL''TBX3''OVOL1''POU2F3'
#'ESRRG''PPARD''REL''TBX3''POU2F3'
#TBX3, POU2F3

tf.cor.sel <- tf.cor.sel[tf.cor.sel %in% rownames(tfmotif_nes_order) & tf.cor.sel %in% rownames(tfmotif_expr_order)]


#for(i in names(res.tfmotif.cor.sel)){
#for(i in names(res.tfmotif.cor.sel.cor0.5)){
for(i in tf.cor.sel){
  options(repr.plot.height=5.5,repr.plot.width=5.5)
  plot(unlist(tfmotif_expr_order[i,]), unlist(tfmotif_nes_order[i,]) , xlab = 'TF gene expression (normalized)',ylab = 'TF NES', pch = 19, cex = 0.5, main = paste('TF gene ',i) )
}


#options(repr.plot.height=5.5,repr.plot.width=5)
options(repr.plot.height=15.5,repr.plot.width=5)
#hp.expr = Heatmap(tfmotif_expr_order[tf.cor.sel.cor0.5,], name = "tfmotif.exp.high_corr", 
hp.expr = Heatmap(tfmotif_expr_order[tf.cor.sel,], name = "tfmotif.exp.high_corr",
             cluster_rows = TRUE, 
             cluster_columns = FALSE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             #col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis::viridis(n = 11,option = "C")),
             row_title = "TFmotif expression (aligned,hclust: ward.D)",
             clustering_method_row = "ward.D", ##ward.D,complete
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 8),
             rect_gp = gpar(col = 'black', lwd=0.5),
             ##border_gp = gpar(col = 'black', lwd = 2),
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
#width = max(grobWidth(textGrob(labels))))

##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
set.seed(1);draw(hp.expr, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750



set.seed(1);res.row.order <- row_order(hp.expr) 
#rowid.order <- rownames(tfmotif_expr_order[tf.cor.sel.cor0.5,])[res.row.order]
rowid.order <- rownames(tfmotif_expr_order[tf.cor.sel,])[res.row.order]

#options(repr.plot.height=5.5,repr.plot.width=5)
options(repr.plot.height=15.5,repr.plot.width=5)
hp.nes = Heatmap(tfmotif_nes_order[rowid.order ,], name = "tfmotif.NES.high_corr", 
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             #col = c('white',colorset),
             row_title = "TFmotif NES (aligned,hclust: ward.D)",
             clustering_method_row = "ward.D", ##ward.D,complete
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 8),
             rect_gp = gpar(col = 'black', lwd=0.5),
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
#width = max(grobWidth(textGrob(labels))))

##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
set.seed(1);draw(hp.nes, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750





##filter aligned tf gene nes and expr mat for low expression tf gene


tfmotif_nes_align <- hp.nes@matrix
tfmotif_expr_align <- hp.expr@matrix
#rowname not equal

all.equal( sort(rownames(tfmotif_expr_align)),sort(rownames(tfmotif_nes_align)) ) #TRUE

tfmotif_expr_align <- tfmotif_expr_align[rownames(tfmotif_nes_align),]

all.equal( rownames(tfmotif_expr_align),rownames(tfmotif_nes_align) ) #TRUE


colnames(tfmotif_nes_align)
'c1''c7''c9''c2''c8''c3''c6''c4''c5'
#'c6''c3''c9''c5''c1''c2''c4'

colnames(tfmotif_expr_align)
'c9''c6''c11''c8''c1''c3''c4''c2''c10'
#'c9''c5''c10''c6''c7''c3''c2'

##exclude low expressed tf

idx1 <- grep ('^MAX$',rownames(tfmotif_expr_align) )
idx2 <-  grep ('^NFATC1$',rownames(tfmotif_expr_align) )

rowid_exclude <- c(rownames(tfmotif_expr_align)[idx1:idx2])

rowid_exclude <- c(rowid_exclude,'MAZ','ZEB1','E2F2','TCF4')


##

tfmotif_expr_align <- tfmotif_expr_align[!rownames(tfmotif_expr_align) %in% rowid_exclude,]
tfmotif_nes_align <- tfmotif_nes_align[!rownames(tfmotif_nes_align) %in% rowid_exclude,]

#tfmotif_expr_align <- tfmotif_expr_align[rowid_mannual,]
#tfmotif_nes_align <- tfmotif_nes_align[rowid_mannual,]

all.equal( rownames(tfmotif_expr_align),rownames(tfmotif_nes_align) )
#TRUE

#options(repr.plot.height=5.6,repr.plot.width=5)
options(repr.plot.height=10.5,repr.plot.width=5)
hp.expr.manual = Heatmap(tfmotif_expr_align, name = "tfmotif.exp.high_corr", 
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             col = circlize::colorRamp2(seq(-0.8,0.8,length.out = 255), colorRampPalette(c('blue','white','red'))(255) ),
             row_title = "TF gene expression (aligned)",
             clustering_method_row = "ward.D", ##ward.D,complete
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 8),
             rect_gp = gpar(col = 'black', lwd=0.5),
             ##border_gp = gpar(col = 'black', lwd = 2),
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
#width = max(grobWidth(textGrob(labels))))

##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
set.seed(1);draw(hp.expr.manual, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750

zScore <- function (m = NULL, min = -2, max = 2, limit = FALSE) 
{
    z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    if (limit) {
        z[z > max] <- max
        z[z < min] <- min
    }
    return(z)
}#ArchR


tfmotif_nes_align.z <- zScore(tfmotif_nes_align,limit = TRUE)


#options(repr.plot.height=5.5,repr.plot.width=5)
options(repr.plot.height=10.5,repr.plot.width=5)
hp.nes.manual = Heatmap(tfmotif_nes_align.z , name = "tfmotif.NES.high_corr", 
#hp.nes.manual = Heatmap(tfmotif_nes_align , name = "tfmotif.NES.high_corr", 
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             #col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             #col = colorset,
             col = circlize::colorRamp2(seq(-0.5,2.5,length.out=length(colorset)), colorset),
             row_title = "TF motif NES (aligned)",
             clustering_method_row = "ward.D", ##ward.D,complete
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 8),
             rect_gp = gpar(col = 'black', lwd=0.5),
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
#width = max(grobWidth(textGrob(labels))))

##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
set.seed(1);draw(hp.nes.manual, heatmap_legend_side = "right",gap = unit(0.1, "cm")) #8.895833 8.843750




######4 quick query taret gene for upstream TF regulators#######

tgene <- 'CSH2' #MAFK? NFE2?
tgene <- 'PAPPA'
tgene <- 'OSMR'
tgene <- 'DTNB'
tgene <- 'GH2'
tgene <- 'LAMA3' #BACH1? #ZBTB21? NFE2?
#tgene <- 'ERV'
tgene <- 'ERVFRD-1' #TP63 TP53 

tgene <- 'SH3TC2$'
tgene <- 'PSG8'


tgene <- 'FLT1' #CEBPB, FOSL2
tgene <- 'ENG' #CEBPB, FOSL2
tgene <- 'INHBA'
tgene <- 'FSTL3' #empty (MYCN according to dorc method)
tgene <- 'LVRN'
tgene <- 'EGLN3' #Egl-9 Family Hypoxia Inducible Factor 3 #CEBPB?

tgene <- 'MYCN'
tgene <- 'POU2F3'

#tgene <- 'ESRRG'

tgene <- 'STAT5A'
tgene <- 'STAT5B'
tgene <- 'STAT4'


tgene <- 'DNMT1' #TP53 TP63
tgene <- 'CDH1'


for(cid in names(tfmotif.list.filter) ){
    tfid_hit <- vector()
    for(tfid in names(tfmotif.list.filter[[cid]]) ){
        peak_tgene <-  tfmotif.list.filter[[cid]][[tfid]][['targetGenes_p2g_combined']]
        #peak_tgene <-  tfmotif.list.filter[[cid]][[tfid]][['nearGenes_combined']]
        ix <- grep (tgene,peak_tgene,value=TRUE)
        if(length(ix) != 0){ cat('Found tgene: ', paste(ix,collapse = ','), ' for gene ',tgene, ' in cid: ',cid, ' in tf: ',tfid,'\n',sep=''); tfid_hit <- c(tfid_hit,tfid)  }
        
    }
    tfid_hit <- sort(unique(tfid_hit))
    cat('cid:',cid,' ',paste(tfid_hit,collapse = ','),'\n',sep='')

    ##plot heatmap
    
    #table(tfid_hit %in% rownames(exprMat.ave.z))
    
    tfid_hit <- c(tgene,tfid_hit) #include tgene as a reference gene
    
    tfid_hit <- tfid_hit[tfid_hit %in% rownames(exprMat.ave.z)] 
    
    tfid_hit.exp.df <-  exprMat.ave.z[tfid_hit ,c('c7','c9','c6','c11','c8','c1','c3','c4','c2','c10','c5')]

    if(nrow(tfid_hit.exp.df) >= 3){ #include tgene
        options(repr.plot.height=12.5,repr.plot.width=4.5)
        hp.hitexpr = Heatmap(as.matrix(tfid_hit.exp.df), name = "zscore",
                     cluster_rows = TRUE, 
                     cluster_columns = FALSE, 
                     show_row_names = TRUE,
                     use_raster = FALSE,
                     #col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis::viridis(n = 11,option = "C")),
                     row_title = paste0("tf expression of given target gene: ",tgene,' in cid: ',cid),
                     clustering_method_row = "ward.D", ##ward.D,complete
                     #clustering_distance_rows  = "pearson",
                     row_names_gp = gpar(fontsize = 8),
                     rect_gp = gpar(col = 'black', lwd=0.5),
                     ##border_gp = gpar(col = 'black', lwd = 2),
                     #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
                     #heatmap_legend_param = list(color_bar = "continuous")
                     ##right_annotation = ha#,heatmap_width=unit(8, "cm")
        ) 
        #width = max(grobWidth(textGrob(labels))))

        ##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
        set.seed(1);draw(hp.hitexpr, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750
        
        
        ##plot the related tf nes df heatmap
        
        set.seed(1);res.row.order <- row_order(hp.hitexpr) 
        #rowid.order <- rownames(tfmotif_expr_order[tf.cor.sel.cor0.5,])[res.row.order]
        rowid.order <- rownames(tfid_hit.exp.df)[res.row.order]
        rowid.order <- rowid.order[!rowid.order %in% tgene]

        #heatmap type 1
        #options(repr.plot.height=5.5,repr.plot.width=5)
        options(repr.plot.height=12.5,repr.plot.width=4.5)
        hp.hitnes = Heatmap(tfmotif.nes.df[rowid.order ,], name = "NES", 
                     cluster_rows = TRUE, 
                     cluster_columns = FALSE, 
                     show_row_names = TRUE,
                     use_raster = FALSE,
                     ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
                     #col = c('white',colorset),
                     row_title = paste0("TFmotif NES  of given target gene: ",tgene,' in cid: ',cid),
                     clustering_method_row = "ward.D", ##ward.D,complete
                     #clustering_distance_rows  = "pearson",
                     row_names_gp = gpar(fontsize = 8),
                     rect_gp = gpar(col = 'black', lwd=0.5),
                     #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
                     #heatmap_legend_param = list(color_bar = "continuous")
                     ##right_annotation = ha#,heatmap_width=unit(8, "cm")
        ) 
        #width = max(grobWidth(textGrob(labels))))

        ##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
        set.seed(1);draw(hp.hitnes, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750

#         #heatmap type 2
#         options(repr.plot.height = 12.5, repr.plot.width = 4.5)
        
#         h.hitnes <- pheatmap(
#                tfmotif.nes.df[rowid.order ,c('c1','c7','c9','c2','c8','c3','c6','c4','c5')], 
#                cluster_cols = FALSE,
#                cluster_rows = FALSE, 
#                treeheight_col = 2.5,
#                treeheight_row = 2.5,
#                #cellwidth = 20, 
#                #cellheight = 2.8,
#                na_col = 'white', 
#                color = colorset,#rev(colorset_pathway),
#                breaks =  seq( 0, max(tfmotif.nes.df[rowid.order ,],na.rm = TRUE), by = max(tfmotif.nes.df[rowid.order ,],na.rm = TRUE)/length(colorset)  )  ,
#                border =TRUE,
#                border_color = 'black',
#                #labels_col = colid, 
# #                labels_row = make_bold_names(tftarget.genesets.freqDF,
# #                                             rownames, 
# #                                             rowid_hi
# #                                            ),

#                ##angle_col = 315,
#                #display_numbers = data.df.text,
#                #number_color = 'white',
#                #fontsize_number = 10,
#                fontsize_col = 20,
#                fontsize_row = 8,
#                main = paste0("TFmotif NES  of given target gene: ",tgene,' in cid: ',cid),
#                silent = FALSE
#                #legend_breaks = c(2,4,6,8,10),
#                #legend_labels = c(2,4,6,8,10),
#                #legend = TRUE
# #                        height = height,
# #                        width = width,
# #                        filename = paste('result_do_GO_Pathway_quickonestep/ck.GO_BP.pvalue_',
# #                                         pvalue,'.qvalue_', qvalue,'.heatmap.pdf' ,sep = '' 
# #                                        )

#               )

        
    }
    
    
}

##get peaks linked to this gene that enriched this tf


##ERVFRD-1: TP63, TP53
grep("ERVFRD-1",tfmotif.list.filter[['c9']][['TP63']][['targetGenes_p2g_combined']],value=TRUE)
#'chr6:11129830-11130560|ERVFRD-1'

grep("ERVFRD-1",tfmotif.list.filter[['c9']][['TP53']][['targetGenes_p2g_combined']],value=TRUE)
#chr6:11129830-11130560|ERVFRD-1


##LAMA3: STAT5A, BACH1?
grep("LAMA3",tfmotif.list.filter[['c5']][['STAT5A']][['targetGenes_p2g_combined']],value=TRUE)
'chr18:23675729-23676608|LAMA3'

grep("LAMA3",tfmotif.list.filter[['c5']][['BACH1']][['targetGenes_p2g_combined']],value=TRUE)
'chr18:23446537-23447232|LAMA3''chr18:23716593-23716862|LAMA3'


##CSH2: MITF
grep("CSH2",tfmotif.list.filter[['c5']][['MITF']][['targetGenes_p2g_combined']],value=TRUE)
#'chr17:63896549-63896943|CSH2'


#FLT1
grep("FLT1",tfmotif.list.filter[['c3']][['CEBPB']][['targetGenes_p2g_combined']],value=TRUE)
grep("FLT1",tfmotif.list.filter[['c3']][['CEBPG']][['targetGenes_p2g_combined']],value=TRUE)
#both tf the same
'chr13:28468163-28469121|FLT1''chr13:28573903-28574421|FLT1'

#ENG
grep("ENG",tfmotif.list.filter[['c3']][['CEBPB']][['targetGenes_p2g_combined']],value=TRUE)
'chr9:127864834-127865205|ENG''chr9:127919368-127919792|ENG''chr9:127850294-127850743|ENG'
#chr9:127850294-127850743|ENG

grep("ENG",tfmotif.list.filter[['c3']][['CEBPG']][['targetGenes_p2g_combined']],value=TRUE)
'chr9:127919368-127919792|ENG'



##hypoxia gene?
hypoxia_gene1 <- c('HDAC2', 'EGLN1', 'EGLN3','ANGPTL4','SLC11A2', 'SLC2A1', 'LIMD1', 'HK2', 'NDRG1','HILPDA','ATG7' ,'CARD16')

hypoxia_gene2 <- c('FLT1', 'VEGFA', 'LDHA', 'HIF1A', 'VHL', 'EDN1', 'EP300', 'TAF1', 'RPA1','MDM2', 'IGFBP3', 'FHL2', 'CDKN1A')


hit_gene1 <- vector()
for(i in hypoxia_gene1){
  hit_gene1 <- c(hit_gene1,grep(i,tfmotif.list.filter[['c3']][['CEBPB']][['targetGenes_p2g_combined']],value = TRUE) )

}
'chr14:34280352-34280894|EGLN3''chr2:74739656-74741665|HK2''chr2:74747587-74747996|HK2''chr2:74874747-74875895|HK2''chr8:133216547-133216905|NDRG1'


hit_gene2 <- vector()
for(i in hypoxia_gene2){
  hit_gene2 <- c(hit_gene2,grep(i,tfmotif.list.filter[['c3']][['CEBPB']][['targetGenes_p2g_combined']],value = TRUE) )

}

'chr13:28468163-28469121|FLT1''chr13:28573903-28574421|FLT1'



##5 select tf as tf-nes-expr heatmap and output eGRN with comnined tfmotif.list.filter.combine







##############output to long format for dotplot with grid#################

wide2long <- function(df1 = NULL,df2=NULL){
    ##turn a wide dataframe to long format, with row_idx, col_idx and row_id, col_id and value 
    row_name <- rownames(df1)
    col_name <- colnames(df1)
    
    #make sure that df2 has the same rowname with df1, the same ncol with df1 (but no need to have the same colname, but be paired)
    
    stopifnot( all.equal(rownames(df1), rownames(df2) ) )
    stopifnot( ncol(df1) == ncol(df2) )
    
    df.long <- data.frame()
    for(i in 1:nrow(df1)){
        for(j in 1:ncol(df1) ){
            row_idx <- i
            col_idx <- j
            row_id <- row_name[i]
            col_id <- col_name[j]
            value1 <- df1[i,j]
            value2 <- df2[i,j]
            df.long <- rbind.data.frame(df.long, data.frame( row_idx = row_idx, col_idx = col_idx, row_id = row_id, col_id= col_id,  value1=value1, value2=value2  ) )
            
        }
        
    }
    
    return(df.long)
}




#tf_keep <- unique(unlist(tf_select))
#tf_keep <- c('STAT5A','STAT4','AR','MITF',)


all.equal(rownames(tfmotif_nes_align.z),rownames(tfmotif_expr_align) )


##reorder rowid before wide2long

##tf before ordered

tf_order <- c('ARID5B','MCTP2','CTNNB1','PPARG','TP63','TEAD1','TFEB','TEAD4','EGR1','TFAP2C','RXRA','BHLHE41','IRX2','SUZ12','PPARA','SP100','TFDP1','SETBP1','TCF7L1','E2F3','ESRRA','POLR2A','NR4A1','GRHL3','TFAP2E','TCF7','UBB','SP4','DLX5','ZNF32','BRF1','TP53','PAXIP1','SUPT20H','E2F4','ZCCHC14','SMARCA4','NR2C1','NR3C1','LCORL','ZSCAN29','TFAP2B','ZNF71','ZNF205','ZFP14','XRCC4','MXD4','AKR1A1','RORA','ZMIZ1','E2F1','TCF3','EZH2','PPARD','GATA3','FOXJ3','ATF3','ZNF217','CREB5','OVOL1','DLX3','CNOT3','FOSL1','SMAD7','TCF12','RELA','STAT3','KLF4','FOSL2','ZNF770','ZNF664','ZBED6','TRIM21','CEBPD','BCL3','NFE2','CEBPG','ZNF799','ZNF611','ZFP2','MEF2D','RNF114','GLIS2','HOXB2','HES1','STAT6','GATA2','MAFG','CEBPA','SP3','RAB7A','CEBPB','MYCN','GCM1','ZNF554','XBP1','ZBTB21','SMARCB1','MITF','STAT4','STAT5A','SMAD6','ZNF83','MAFK','HES2','MAFF','BACH1','BCL6','STAT5B')
#tf_order <- c('CTNNB1','PPARG','TP63','NR2F2','TCF7L2','TEAD1','TFEB','TFCP2L1','THRB','TEAD4','PPARD','ATF3','SMAD7','ZNF655','FOSL1','DLX5','NR1H3','ZNF32','PAXIP1','UBP1','TP53','TFCP2','LCORL','NR2C1','KLF5','ZFP14','TFAP2B','ZNF71','ZNF205','MXD4','BHLHE41','IRX2','TFAP2C','RXRA','PPARA','TCF7L1','GCM1','CEBPB','MAFG','CEBPA','SP3','STAT3','FOSL2','HES1','ZFP2','GLIS2','CEBPD','CEBPG','NFE2','STAT6','MYCN','STAT4','XBP1','ZBTB21','MITF','SMAD6','ZNF83','STAT5A','MAFK','HES2','MAFF','BACH1','BCL6','STAT5B')

all.equal(rownames(tfmotif_nes_align.z),tf_order) #TRUE


#tf_order <- c('ARID5B','MCTP2','TP63','TEAD1','TEAD4','ESRRA','POLR2A','EGR1','SUPT20H','E2F4','ZCCHC14','SMARCA4','BRF1','SUZ12','SP100','TFDP1','SETBP1','E2F3','AKR1A1','MAZ','XRCC4','ZMIZ1','E2F1','EZH2','OVOL1','ATF3','ZNF217','FOXJ3','TCF12','CNOT3','GRHL3','EP300','FOSL1','RELA','NFE2L1','MAX','UBB','SP4','IRX4','STAT3','KLF4','FOSL2','ZNF664','MEF2D','JUNB','TRIM21','CEBPD','BCL3','NFE2','CEBPG','STAT6','GATA3','GRHL1','GCM1','GATA2','MAFG','RAB7A','CEBPB','SP3','MYCN','ZNF554','MITF','STAT4','STAT5A','RNF114','PURA','SMARCB1','MAFK','BACH1','GLIS3','BCL6','STAT5B')

#tf_order <- c('ARID5B','MCTP2','TP63','TEAD1','TEAD4','ESRRA','POLR2A','EGR1','SUPT20H','E2F4','ZCCHC14','SMARCA4','BRF1','SUZ12','SP100','TFDP1','SETBP1','E2F3','GATA3','FOXJ3','OVOL1','ATF3','ZNF217','GCM1','GRHL1','GATA2','MAFG','RAB7A','ZNF554','MITF','STAT4','STAT5A','MAFK','BACH1','GLIS3','BCL6','STAT5B','TCF12','CNOT3','GRHL3','EP300','FOSL1','RELA','NFE2L1','MAX','UBB','SP4','IRX4','AKR1A1','MAZ','XRCC4','ZMIZ1','E2F1','EZH2','MYCN','STAT3','KLF4','FOSL2','ZNF664','MEF2D','JUNB','TRIM21','CEBPD','BCL3','NFE2','CEBPG','STAT6','RNF114','PURA','SMARCB1','SP3')

#tf_order <- c('MCTP2','TP63','PPARG','NR2F2','TEAD4','TEAD1','ESRRG','PPARD','REL','DLX3','POLR2A','SP4','SP110','ZNF16','RAD21','CPSF4','MZF1','RBBP5','VEZF1','HNF4G','TCF7','TEF','ZNF263','ZMAT2','UBB','EGR1','TFAP2B','POLE4','NFYA','PATZ1','KLF13','SUPT20H','APEX1','NFYC','BRF1','ZNF808','ZNF148','E2F4','GRHPR','ASCC1','TP53','ZCCHC14','ZNF43','PPARA','ANXA11','CTBP2','TCF7L1','GATA2','NR4A3','AFF4','FOSL2','NR2C2','JUNB','ZDHHC5','ZBTB7B','OVOL1','ZNF217','EGR3','SP2','ZNF134','ZNF460','KLF16','STAT3','RELA','TRIM21','CEBPG','ZFP62','HOXB2','RNF114','KLF15','BATF','CEBPD','BCL3','VDR','NFE2','GLIS2','STAT6','ATF3','BACH1','MXD1','RAB7A','CEBPB','GCM1','POU2F1','ZNF91','MAFK','STAT4','AR','ELK1','ARNT2','RUNX1','STAT5A','XBP1','JUND','ETV5','SMARCB1','ETV4','MITF','HP1BP3','SMAD6','MYCN','GLIS3','TFAP2A','MXI1','SREBF1','BCL6','MYLK','TBL1XR1','STAT5B')


rowid.nes <- rownames(tfmotif_nes_align.z)
rowid.expr <- rownames(tfmotif_expr_align)

all.equal(rowid.nes,rowid.expr) #TRUE
all.equal(tf_order, rowid.nes) #TRUE




########calculate correlation plot and value for selected TFs  below

i <- 'STAT4'
i <- 'STAT5A'
i <- 'MITF'
i <- 'STAT5B'
i <- 'CEBPB'
#i <- 'KLF5'
i <- 'CEBPG'
i <- 'FOSL2'
i <- 'TP53'
i <- 'TP63'

unlist(tfmotif_expr_order[i,])
unlist(tfmotif_nes_order[i,])
cor(unlist(tfmotif_expr_order[i,]), unlist(tfmotif_nes_order[i,]))
#0.833021795964144 STAT4
#0.813862083337557 STAT5A
#0.0128766126069388 MITF
#0.543131476214119 STAT5B
#0.260419445911582 CEBPB
#0.614995318526313 CEBPG
#0.468366405808652 FOSL2
#0.799246617537414 TP53
#0.909793214810372 TP63



#0.94907751822890 STAT4
#0.693993551910119 STAT5A
#0.566398891258084 MITF
#0.439910975667492 STAT5B
#-0.119836779170299 CEBPB

#0.819235091282893 STAT4
#0.75978 STAT5A


res.tfmotif.cor.new <- vector()
for(i in order_list){
  res.tfmotif.cor.new[i] <- cor(unlist(tfmotif_expr_align[i,]), unlist(tfmotif_nes_align.z[i,]))
  
}



color_set = list(
    'c1' = '#3D3F69',
    'c7' = '#4776B2',
    'c9' = 'darkgreen',
    'c2' = '#FCC140',
    'c8' = '#F95944',
    'c3' = '#8D1541',
    'c6' = '#FB8D3C',
    'c4' = '#C31240',
    'c5' = '#4F1C47'
#            'c6'='#74add1',
#            'c3'='#4575b4',
#            'c9'='darkgreen',
#            'c5'='#ffffbf',
#            'c7'='#fdae61',
#            'c1'='#f46d43',
#            'c8'='#fee090',
#            'c2'='#d73027',
#            'c4'='#a50026'
)



#par(mfrow = c(4,5))
par(mfrow = c(3,3))
options(repr.plot.height = 12, repr.plot.width = 12)


all.equal(rownames(tfmotif_expr_align),rownames(tfmotif_nes_align.z)) #TRUE

options(repr.plot.height=15,repr.plot.width=15)
par(mfrow = c(3,3))
for(i in order_list){
#for(i in rowid_ctb_tf){
  #options(repr.plot.height=5.5,repr.plot.width=5.5)
  plot(unlist(tfmotif_expr_align[i,]), unlist(tfmotif_nes_align.z[i,]) , 
       xlab = 'TF gene expression (normalized)',ylab = 'TF NES', pch = 19, cex = 2, 
       main = paste('TF gene ',i,' correlation: \n r = ',round(res.tfmotif.cor.new[i],digits = 3),sep=''),
       col = unlist(color_set[colnames(tfmotif_nes_align.z)]),
       cex.main = 2
       #cex.xlab = 3
      )
  abline(h=0.3,lty=2)
  text(unlist(tfmotif_expr_align[i,]), unlist(tfmotif_nes_align.z[i,]), 
       labels = colnames(tfmotif_nes_align.z),
       pos = 1,
       col = unlist(color_set[colnames(tfmotif_nes_align.z)]),
       cex = 2
      )
  
}


########


##tf rowid show in tf-nes-expr heatmpa


order_list <-  c(  'TP63','TEAD4','E2F3','TP53', #CTB
                   'RELA','EP300','GRHL3','ZNF217','FOSL1','GCM1', 'MAFK','MAFG',#Fusion and Nascent (EP300 will be lost)
                   'GLIS3','CEBPB','CEBPG','FOSL2','SP3', #Mature2 #GLIS3 will lost
                   'STAT5B','STAT4','STAT5A','ZNF554','MITF','SMAD6','GLIS2','BACH1','ZBTB21', 'NFE2' #Premature1 and Mature1
                )

order_list_pycistarget <-  c('TP63','TEAD4','TP53','GCM1','MAFK','MAFG','CEBPB','CEBPG','FOSL2','SP3','STAT5B','STAT4','STAT5A','MITF','SMAD6','GLIS2','BACH1','ZBTB21', 'NFE2') #pycistarget heatmap rowid



order_list_icistarget <- c('TP63','TEAD4','E2F3','RELA','ZNF217','EP300','GRHL3','FOSL1','MAFK','GCM1','GLIS3','MAFG','CEBPG','CEBPB','FOSL2','SP3','MEF2D','ZNF554','MITF','STAT5A','STAT5B','STAT4') #i-cisTarget heatmap rowid

# sum(duplicated(order_list)) #0
# sum(duplicated(order_list_pycistarget)) #0
# sum(duplicated(order_list_icistarget)) #0


# length(order_list) #27
# length((order_list_pycistarget)) #19
# length(order_list_icistarget) #22


# length(intersect( order_list, c(order_list_pycistarget,order_list_icistarget) )) #27

# all.equal(sort(order_list), sort(intersect( order_list, c(order_list_pycistarget,order_list_icistarget) )))
# #TRUE  order_list in order_list_pycistarget + order_list_icistarget


##check
sort(order_list)
'BACH1''CEBPB''CEBPG''E2F3''EP300''FOSL1''FOSL2''GCM1''GLIS2''GLIS3''GRHL3''MAFG''MAFK''MITF''NFE2''RELA''SMAD6''SP3''STAT4''STAT5A''STAT5B''TEAD4''TP53''TP63''ZBTB21''ZNF217''ZNF554'

sort(unique(c(order_list_pycistarget, order_list_icistarget)) )
'BACH1''CEBPB''CEBPG''E2F3''EP300''FOSL1''FOSL2''GCM1''GLIS2''GLIS3''GRHL3''MAFG''MAFK''MEF2D''MITF''NFE2''RELA''SMAD6''SP3''STAT4''STAT5A''STAT5B''TEAD4''TP53''TP63''ZBTB21''ZNF217''ZNF554'

#'MEF2D'



sum(duplicated(order_list)) #0
table(order_list %in% tf_order )
FALSE  TRUE 
    2    25

# TRUE 
#   19 

# TRUE 
#   20

# TRUE 
#   17

# TRUE
# 16

# TRUE
# 22

# TRUE
# 20

# TRUE 
#   19

# TRUE 
#   20

order_list[!order_list %in% tf_order]
'EP300''GLIS3'

order_list <- order_list[order_list %in% tf_order]


idx <- which(tf_order %in% order_list)
#for(i in idx){cat(i,',',sep='')}
#1,5,10,48,49,51,53,57,64,66,74,87,88,89,91,92,98,102,104,106

tf_order[idx]

'TP63''TEAD4''E2F3''GRHL3''TP53''ZNF217''FOSL1''RELA''FOSL2''NFE2''CEBPG''GLIS2''MAFG''SP3''CEBPB''GCM1''ZNF554''ZBTB21''MITF''STAT4''STAT5A''SMAD6''MAFK''BACH1''STAT5B'

#'TP63''TEAD4''TP53''GCM1''CEBPB''MAFG''SP3''FOSL2''GLIS2''CEBPG''NFE2''STAT4''ZBTB21''MITF''SMAD6''STAT5A''MAFK''BACH1''STAT5B'

#'TP63''TEAD4''TP53''KLF5''GCM1''CEBPB''MAFG''SP3''FOSL2''GLIS2''CEBPG''NFE2''STAT4''ZBTB21''MITF''SMAD6''STAT5A''MAFK''BACH1''STAT5B'

#'TP63''TEAD4''TP53''KLF5''GCM1''CEBPB''MAFG''SP3''FOSL2''GLIS2''CEBPG''NFE2''STAT4''ZBTB21''MITF''SMAD6''STAT5A''MAFK''BACH1''STAT5B'

#'TP63''TEAD4''TP53''KLF5''GCM1''CEBPB''MAFG''SP3''FOSL2''GLIS2''CEBPG''STAT4''MITF''SMAD6''STAT5A''MAFK''STAT5B'

#'TP63''TP53''TEAD4''KLF5''GCM1''MAFK''MAFG''CEBPB''CEBPG''FOSL2''SP3''GLIS2''STAT5A''STAT5B''STAT4''MITF''SMAD6'

#'TP63''TEAD4''TP53''GCM1''CEBPB''MAFG''SP3''FOSL2''GLIS2''CEBPG''STAT4''MITF''SMAD6''STAT5A''MAFK''STAT5B'

#'TP63''TEAD4''E2F3''ZNF217''GRHL3''EP300''FOSL1''RELA''FOSL2''MEF2D''CEBPG''GCM1''MAFG''CEBPB''SP3''ZNF554''MITF''STAT4''STAT5A''MAFK''GLIS3''STAT5B'

#'TP63''TEAD4''E2F3''EP300''FOSL1''RELA''FOSL2''MEF2D''CEBPG''GCM1''MAFG''CEBPB''SP3''ZNF554''MITF''STAT4''STAT5A''MAFK''GLIS3''STAT5B'

#'TP63''TEAD4''E2F3''GCM1''MAFG''ZNF554''MITF''STAT4''STAT5A''MAFK''GLIS3''STAT5B''EP300''FOSL1''RELA''FOSL2''MEF2D''CEBPG''SP3'

#'MCTP2''TEAD4''DLX3''GATA2''NR4A3''FOSL2''JUNB''ZNF217''RELA''CEBPG''VDR''STAT4''AR''ELK1''RUNX1''STAT5A''MITF''GLIS3''MXI1''BCL6'

all.equal(sort(order_list),sort(tf_order[idx])  )
#TRUE

tf_order[idx] <- order_list


sum(duplicated(tf_order)) #0


tfmotif_align_long <- wide2long(df1 = tfmotif_nes_align.z[tf_order,], df2 = tfmotif_expr_align[tf_order,])


#tfmotif_align_long <- wide2long(df1 = tfmotif_nes_align, df2 = tfmotif_expr_align)
#tfmotif_align_long <- wide2long(df1 = tfmotif_nes_align.z, df2 = tfmotif_expr_align)



colnames(tfmotif_align_long) <- c('y_idx','x_idx','y_id','x_id','NES','expression')


#tf_hi <- c('STAT4','STAT5A','JUNB','CEBPG','JUND','RELA','MAFG','TP63','PPARG','NR2F2','AR','MYCN','MITF','VDR')
#tf_hi <- c('STAT4','STAT5A','JUNB','CEBPG','JUND','ZNF66','ZNF134','RELA','MAFG','TP63','PPARG','NR2F2')
##chosen for network construction

tf_hi <- c('STAT4','STAT5A','MITF','FOSL2','CEBPG','CEBPB','TEAD4','TP53')
#tf_hi <- c('STAT4','STAT5A','MITF','FOSL2','CEBPG','CEBPB','RELA','GRHL3','TEAD4','TP63')

#tf_hi <- c('STAT4','STAT5A','AR','MITF','FOSL2','JUNB','CEBPG','JUND','RELA','ZNF217','GRHL1','TEAD4','TP63','PPARG')


tfmotif_align_long$edge <- ifelse(tfmotif_align_long$y_id %in% tf_hi, 'black','grey' )





#####################ggplot dot plot with grid or use external python script####################

# ggplot(data = tfmotif_align_long, aes(x= x_id, y = y_id, size = NES, col = expression) ) +
#   geom_point()


temp <- tfmotif_align_long$expression

options(repr.plot.height=7.5,repr.plot.width=7.5)
hist(temp,breaks = 100)

#cutoff by q45
#round(quantile(x = temp,probs = seq(0,1,0.1)),3)
#0% -1.174 10% -0.218 20% -0.15 30% -0.103 40% -0.07 50% -0.044 60% 0 70% 0.053 80% 0.137 90% 0.284 100% 1.571

qs <- round(quantile(x = temp,probs = seq(0,1,0.05)),3)
qs['40%']
40%: -0.07
45% -0.058

#temp[temp < qs['45%']] <- qs['45%']
#min(temp) == qs['45%'] #TRUE


temp[temp < qs['40%']] <- qs['40%']
min(temp) == qs['40%'] #TRUE

#transform to positive number for dot size
temp <- temp + abs(min(temp))
min(temp) == 0 #TRUE

#temp[temp < 0] <- 0


##scale to min max
temp.min <- min(temp)
temp.max <- max(temp)

# (temp - temp.min)/(temp.max-temp) = x/(1-x)
# x*(temp.max-temp) = (temp - temp.min)*1-  (temp - temp.min)* x

# x*temp.max-x*temp = temp - temp.min - x*temp + x* temp.min

# x*temp.max = temp - temp.min + x* temp.min

# x*temp.max - x* temp.min = temp - temp.min
# x * (temp.max -temp.min) = temp - temp.min

x = (temp - temp.min)/(temp.max -temp.min)

tfmotif_align_long$expression.transform <- x

# options(repr.plot.height=5.5,repr.plot.width=8) #too urgly, use python script
# ggplot(data = tfmotif_align_long, aes(x= x_id, y = y_id, size = expression.transform, col = NES) ) +
#   geom_point() +
#   scale_y_discrete(lim=rev)


#tfmotif_align_long <- tfmotif_align_long[!tfmotif_align_long$y_id %in% c('EBF1','GLI2'),]



tfmotif_align_long.keep <- subset(tfmotif_align_long, y_id %in% order_list)


outdir <- 'result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005'
saveRDS(tfmotif_align_long,paste0(outdir,'/tfmotif_align_long.rds') )
saveRDS(tfmotif_align_long,paste0(outdir,'/tfmotif_align_long.cutoffq40.rds') )


#write.table(tfmotif_align_long,file = 'tfmotif_align_long.txt',sep = '\t',col.names = TRUE,row.names = FALSE)


write.table(tfmotif_align_long.keep,file = paste0(outdir,'/tfmotif_align_long.txt'),sep = '\t',col.names = TRUE,row.names = FALSE)
write.table(tfmotif_align_long.keep,file = paste0(outdir,'/tfmotif_align_long.cutoff.q40.txt'),sep = '\t',col.names = TRUE,row.names = FALSE)




###################extract data from tfmotif.list.filter and output tf-peak-targetgene table###############


##tf_hi <- c('STAT4','STAT5A','JUNB','CEBPG','JUND','ZNF66','ZNF134','RELA','MAFG','TP63','PPARG','NR2F2')
tfmotif.list #TF-motif with multiple
tfmotif.list.filter #TF-motif-choose one and get heatmap, targetgene, use this
#tfmotif.list.filter.withdup


tf_select_full <- list( #rectangle area in tf nes-expr heatmap
    'c1' = c('TP63','TEAD4','E2F3','TP53'),
    'c7' = c('TP63','TEAD4','E2F3','TP53'),
    'c9' = c('TP63','TP53','GCM1','RELA','GRHL3','ZNF217','FOSL1','MAFK','MAFG'),
    'c2' = c('MAFK','MAFG','GCM1'),
    'c8' = c('GCM1'),
    'c3' = c('CEBPG','CEBPB','FOSL2','SP3'),
    #'c6' = c(),
    'c4' = c('STAT5B','STAT4','STAT5A','ZNF554','MITF','SMAD6','GLIS2'),
    'c5' = c('STAT5B','STAT4','STAT5A','ZNF554','MITF','BACH1','ZBTB21','NFE2')
)




#check 
for(cid in names(tf_select_full) ){
    tfs <- tf_select_full[[cid]]
    tfs <- tfs[tfs %in% names(tfmotif.list.filter[[cid]])]
    #tfs <- tfs[tfs %in% names(tfmotif.list.filter.withdup[[cid]])]
    cat(cid,' tfs found: ',paste(tfs,collapse = ','),'\n',sep='' )
    
    
}

c1 tfs found: TP63,TEAD4,E2F3,TP53
c7 tfs found: TP63,TEAD4,E2F3,TP53
c9 tfs found: TP63,TP53,RELA,GRHL3,ZNF217,FOSL1,MAFK,MAFG
c2 tfs found: MAFK,MAFG,GCM1
c8 tfs found: GCM1
c3 tfs found: CEBPG,CEBPB,FOSL2,SP3
c4 tfs found: STAT5B,STAT4,STAT5A,ZNF554,MITF,SMAD6,GLIS2
c5 tfs found: STAT5B,STAT4,STAT5A,ZNF554,MITF,BACH1,ZBTB21,NFE2


# c1 tfs found: TP63,TEAD4,TP53
# c7 tfs found: TP63,TEAD4,TP53
# c9 tfs found: TP63,MAFK,MAFG
# c2 tfs found: MAFK,MAFG,GCM1
# c8 tfs found: GCM1
# c3 tfs found: CEBPG,CEBPB,FOSL2,SP3
# c4 tfs found: STAT5B,STAT4,STAT5A,MITF,SMAD6,GLIS2
# c5 tfs found: STAT5B,STAT4,STAT5A,BACH1,ZBTB21,NFE2

# c1 tfs found: TP63,TEAD4,TP53
# c7 tfs found: TP63,TEAD4,TP53
# c9 tfs found: MAFK,MAFG
# c2 tfs found: MAFK,GCM1
# c8 tfs found: GCM1
# c3 tfs found: CEBPG,CEBPB,FOSL2,SP3
# c4 tfs found: GLIS2,STAT5A,STAT5B,STAT4,MITF
# c5 tfs found: STAT5A,STAT5B,STAT4


# c1 tfs found: TP63,TEAD4,E2F3
# c7 tfs found: TP63,TEAD4,E2F3
# c9 tfs found: RELA,ZNF217,EP300,GRHL3,FOSL1
# c2 tfs found: MAFK
# c8 tfs found: GLIS3,GCM1
# c3 tfs found: MAFG,CEBPG,CEBPB,FOSL2,SP3
# c4 tfs found: MEF2D,ZNF554,STAT5A,STAT5B,STAT4
# c5 tfs found: ZNF554,MITF,STAT5A,STAT5B,STAT4

# c6 tfs found: TP63,PPARG,TEAD4
# c3 tfs found: TP63,PPARG,NR2F2
# c9 tfs found: RELA,ZNF217,TP63,GRHL1,EGR3
# c5 tfs found: BCL6,TEAD3
# c1 tfs found: JUNB,JUND,CEBPG,FOSL2,PPARD,GRHL1
# c2 tfs found: STAT4,STAT5A,AR,MITF,ESRRG

##after filtering (full tf set but remove lost)
tf_select <- list(
    
    'c1' = c('TP63','TEAD4','E2F3','TP53'),
    'c7' = c('TP63','TEAD4','E2F3','TP53'),
    'c9' = c('TP63','TP53','RELA','GRHL3','ZNF217','FOSL1','MAFK','MAFG'), #no GCM1 enrichment but expression
    'c2' = c('MAFK','MAFG','GCM1'),
    'c8' = c('GCM1'),
    'c3' = c('CEBPG','CEBPB','FOSL2','SP3'),
    #'c6' = c(),
    'c4' = c('STAT5B','STAT4','STAT5A','ZNF554','MITF','SMAD6','GLIS2'),
    'c5' = c('STAT5B','STAT4','STAT5A','ZNF554','MITF','BACH1','ZBTB21','NFE2')
    
#         'c1' = c('TP63','TEAD4','TP53'),
#     'c7' = c('TP63','TEAD4','TP53'),
#     'c9' = c('TP63','MAFK','MAFG'),
#     'c2' = c('MAFK','MAFG','GCM1'),
#     'c8' = c('GCM1'),
#     'c3' = c('CEBPG','CEBPB','FOSL2','SP3'),
#     #'c6' = c(),
#     'c4' = c('STAT5B','STAT4','STAT5A','MITF','SMAD6','GLIS2'),
#     'c5' = c('STAT5B','STAT4','STAT5A','BACH1','ZBTB21','NFE2')
    
#     'c1' = c('TP63','TEAD4','E2F3'),
#     'c7' = c('TP63','TEAD4','E2F3'),
#     'c9' = c('RELA','ZNF217','EP300','GRHL3','FOSL1'),
#     'c2' = c('MAFK'),
#     'c8' = c('GLIS3','GCM1'),
#     'c3' = c('MAFG','CEBPG','CEBPB','FOSL2','SP3'),
#     #'c6' = c(),
#     'c4' = c('MEF2D','ZNF554','STAT5A','STAT5B','STAT4'),
#     'c5' = c('ZNF554','MITF','STAT5A','STAT5B','STAT4')
    
    
#    'c6' = c('TP63','PPARG','TEAD4'),#,'NR2F2'), #NR2F2 expressed (highly), but has NES = 0 (filled value, no found actually)
#     'c3' = c('TP63','PPARG','NR2F2'),
#     'c9' = c('RELA','ZNF217','TP63','GRHL1','EGR3'),
#     'c5' = c('BCL6','TEAD3'),
#     'c1' = c('JUNB','JUND','CEBPG','FOSL2','PPARD','GRHL1'),
#     'c2' = c('STAT4', 'STAT5A' ,'AR','MITF','ESRRG' )

)


##TF bundle 1 (reduced tf set for clear eGRN construction, no more than 4 tf for each cluster, select most unique enriched and expressed tf)

tf_select_1 <- list(
    
   'c1' = c('TP63','TEAD4','TP53'),
    'c7' = c('TP63','TEAD4','TP53'),
    'c9' = c('TP53','RELA','GRHL3','ZNF217'), #no GCM1 enrichment but expression
    'c2' = c('MAFK','MAFG','GCM1'),
    'c8' = c('GCM1'),
    'c3' = c('CEBPB','FOSL2','SP3'),
    #'c6' = c(),
    'c4' = c('STAT5B','STAT5A','MITF','SMAD6'),
    'c5' = c('STAT5B','STAT5A','STAT4','MITF')
    
    
    
    
#     'c1' = c('TP63','TEAD4','TP53'),
#     'c7' = c('TP63','TEAD4','TP53'),
#     'c9' = c('MAFK','MAFG'),
#     'c2' = c('MAFK','GCM1'),
#     'c8' = c('GCM1'),
#     'c3' = c('CEBPG','CEBPB','FOSL2','SP3'),
#     #'c6' = c(),
#     'c4' = c('GLIS2','STAT5A','STAT5B','STAT4','MITF'),
#     'c5' = c('STAT5A','STAT5B','STAT4')
    
#     'c1' = c('TP63','TEAD4','E2F3'),
#     'c7' = c('TP63','TEAD4','E2F3'),
#     'c9' = c('RELA','ZNF217','EP300','GRHL3','FOSL1'),
#     'c2' = c('MAFK'),
#     'c8' = c('GLIS3','GCM1'),
#     'c3' = c('MAFG','CEBPG','CEBPB','FOSL2','SP3'),
#     #'c6' = c(),
#     'c4' = c('MEF2D','ZNF554','STAT5A','STAT5B','STAT4'),
#     'c5' = c('ZNF554','MITF','STAT5A','STAT5B','STAT4')
    
    
    
#    'c6' = c('TP63','PPARG','TEAD4'),#,'NR2F2'), #NR2F2 expressed (highly), but has NES = 0 (filled value, no found actually)
#     'c3' = c('TP63','PPARG','NR2F2'),
#     'c9' = c('RELA','ZNF217','TP63','GRHL1'),
#     'c5' = c('BCL6','TEAD3'),
#     'c1' = c('JUNB','CEBPG','FOSL2','PPARD','GRHL1'),
#     'c2' = c('STAT4', 'STAT5A' ,'AR','MITF','ESRRG' )

)


tf_select_2 <- list( #c3 add CEBPG, c9 add MAFK del TP53?
    
    
   'c1' = c('TP63','TEAD4','TP53'),
    'c7' = c('TP63','TEAD4','TP53'),
    'c9' = c('RELA','GRHL3','ZNF217','MAFK'), #no GCM1 enrichment but expression
    'c2' = c('MAFK','MAFG','GCM1'),
    'c8' = c('GCM1'),
    'c3' = c('CEBPB','CEBPG','FOSL2','SP3'),
    #'c6' = c(),
    'c4' = c('STAT5B','STAT5A','MITF','SMAD6'),
    'c5' = c('STAT5B','STAT5A','STAT4','MITF')
    
#    'c6' = c('TP63','PPARG','TEAD4'),#,'NR2F2'), #NR2F2 expressed (highly), but has NES = 0 (filled value, no found actually)
#     'c3' = c('TP63','PPARG','NR2F2'),
#     'c9' = c('RELA','ZNF217','TP63','GRHL1'),
#     'c5' = c('BCL6','TEAD3'),
#     'c1' = c('JUNB','CEBPG','FOSL2'),
#     'c2' = c('STAT4', 'STAT5A' ,'AR','MITF')

)


##TF bundle 3 #add ERVFRD-1 regulatory TF TEAD4

tf_select_3 <- list( #to reduce CTB fusion tf number
    
       
   'c1' = c('TP63','TEAD4','TP53'),
    'c7' = c('TP63','TEAD4','TP53'),
    'c9' = c('RELA','MAFK'), #no GCM1 enrichment but expression
    'c2' = c('MAFK','MAFG','GCM1'),
    'c8' = c('GCM1'),
    'c3' = c('CEBPB','CEBPG','FOSL2','SP3'),
    #'c6' = c(),
    'c4' = c('STAT5B','STAT5A','MITF','SMAD6'),
    'c5' = c('STAT5B','STAT5A','STAT4','MITF')
    
#    'c6' = c('TP63','PPARG','TEAD4'),#,'NR2F2'), #NR2F2 expressed (highly), but has NES = 0 (filled value, no found actually)
#     'c3' = c('TP63','PPARG','NR2F2'),
#     'c9' = c('RELA','ZNF217','GRHL1','TEAD4'),
#     'c5' = c('BCL6','TEAD3'),
#     'c1' = c('JUNB','CEBPG','FOSL2'),
#     'c2' = c('STAT4', 'STAT5A' ,'AR','MITF' )

)


#tf_select <- tf_select

#tf_select <- tf_select_1

#tf_select <- tf_select_2

tf_select <- tf_select_3

##outdir

#outdir <- 'tf_peak_gene-network_useuniq/'
#outdir <- 'tf_peak_gene-network_useuniq/tf_select_test'

#use uniq p2g to remove self looping

#outdir <- 'result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005/tf_peak_gene-network_useuniq/tf_select'
#outdir <- 'result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005/tf_peak_gene-network_useuniq/tf_select_1'
#outdir <- 'result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005/tf_peak_gene-network_useuniq/tf_select_2'
outdir <- 'result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005/tf_peak_gene-network_useuniq/tf_select_3'
#outdir <- 'tf_peak_gene-network_useuniq/tf_select_2'
#outdir <- 'tf_peak_gene-network_useuniq/tf_select_3'

if(!dir.exists(outdir) ){dir.create(outdir,recursive = TRUE)}


awk '($6 > 0.6 && $6 != "Cor"){print $6}' network.table.df.c5.txt |wcl
159

awk '($6 > 0.5 && $6 != "Cor"){print $6}' network.table.df.c5.txt |wcl
298


#awk '($6 > 0.6 && $6 != "Cor"){print $6}' network.table.df.c6.txt |wcl
#796

# awk '($6 > 0.6 && $6 != "Cor"){print $6}' network.table.df.c3.txt |wcl
# #654

# awk '($6 > 0.5 && $6 != "Cor"){print $6}' network.table.df.c1.txt |wcl
# 160 of total 218


cutoff.cor <- list(
    
    
    'c1' = 0,
    'c7' = 0,
    'c9' = 0,
    'c2' = 0,
    'c8' = 0,
    'c3' = 0,
    'c6' = 0,
    'c4' = 0,
    'c5' = 0

    
    
#     'c6' = 0.82,
#     'c3' = 0.82,
#     'c9' = 0,
#     'c5' = 0,
#     'c1' = 0,
#     'c2' = 0


)



#::::::::::::::::::::##################start to build network ##############::::::::::::::::::::::#


#save tf_select related rds

if(!dir.exists( paste0(outdir,'/data_rds/') ) ){dir.create(paste0(outdir,'/data_rds/'),recursive = TRUE)}

#tf_select
saveRDS(tf_select,paste0(outdir,'/data_rds/tf_select.rds') )


# ##exprPerc mat
# #exprMat.perc <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/chromVAR_TF_specific/full12526/TF_dev-vs-TF_expression/exprMat.perc.rds')
# #243070 x 5

# exprMat.perc <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.seurat_harmony/exprMat.perc.rds')
# 320452 × 5

# exprMat.perc$id <- paste0('c',exprMat.perc$id)
# table(exprMat.perc$id)
# c1   c10   c11    c2    c3    c4    c5    c6    c7    c8    c9 
# 29132 29132 29132 29132 29132 29132 29132 29132 29132 29132 29132 

# #   c1   c10    c2    c3    c4    c5    c6    c7    c8    c9 
# # 24307 24307 24307 24307 24307 24307 24307 24307 24307 24307

# exprMat.perc$features.plot <- as.character(exprMat.perc$features.plot)


# #exprMat

# colnames(exprMat.ave.z) <- paste0('c',colnames(exprMat.ave.z))

# subset(exprMat.perc, features.plot == 'PAPPA')
# exprMat.ave.z['PAPPA',]

# cor(subset(exprMat.perc, features.plot == 'PAPPA')$avg.exp.scaled, unlist(exprMat.ave.z['PAPPA',]))
# 0.987113377671874
# #0.982409310602705

# cor(subset(exprMat.perc, features.plot == 'PAPPA')$avg.exp, unlist(exprMat.ave.z['PAPPA',]))
# 0.987113377671874
# #0.982409310602705

# cor(subset(exprMat.perc, features.plot == 'PAPPA')$pct.exp, unlist(exprMat.ave.z['PAPPA',]))
# 0.984788632221227
# #0.967340836587857

# sum(duplicated(rownames(exprMat.ave.z)))  #0
# exprMat.ave.z$geneid <- rownames(exprMat.ave.z) #to avoid rowid duplication when extract data


# ##peak accessibility mat
# #peakMat.aggre <- readRDS('/sda/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/cicero_Granja/peakMat.aggre.rds')

# peakMat.aggre <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/cicero_Granja/peakMat.aggre.rds')

# colnames(peakMat.aggre) <- paste( 'c', colnames(peakMat.aggre),sep='')

# rowid <- sapply(stringr::str_split(rownames(peakMat.aggre),pattern = "_",n=3),function(x){ paste0( x[1],":",x[2],"-",x[3]  )  })

# sum(duplicated(rowid) ) #0

# rownames(peakMat.aggre) <- rowid

# sum(duplicated(rownames(peakMat.aggre))) #0

# peakMat.aggre$peakid <- rownames(peakMat.aggre) #to avoid rowid duplication 


#dir.create('data_rds')

##read general data #link from i-cisTarget result dir

#tfmotif.list.filter
saveRDS(tfmotif.list.filter,'data_rds/tfmotif.list.filter.rds') #combined of two ranking db

#map_cluster
saveRDS(map_cluster,'data_rds/map_cluster.rds')


p2g.filter.uniq.modify <- readRDS('data_rds/p2g.filter.uniq.modify.rds')  #46742 × 11
exprMat.perc <- readRDS('data_rds/exprMat.perc.rds') #320452 × 5
exprMat.ave.z <- readRDS('data_rds/exprMat.ave.z.rds') #29132 × 12
peakMat.aggre <- readRDS('data_rds/peakMat.aggre.rds') #274189 × 10
deg.list.up <- readRDS('data_rds/deg.list.up.rds')


#saveRDS(exprMat.perc,'data_rds/exprMat.perc.rds')
#saveRDS(exprMat.ave.z,'data_rds/exprMat.ave.z.rds')
#saveRDS(peakMat.aggre,'data_rds/peakMat.aggre.rds')
#saveRDS(deg.list.up,'data_rds/deg.list.up.rds')


##deg labeling

length(deg.list.up[['c2']]) #188
length(deg.list.up[['c10']]) #242
length(deg.list.up[['c1']]) #117
length(deg.list.up[['c3']]) #222
length(deg.list.up[['c11']]) #814
length(deg.list.up[['c8']]) #135
length(deg.list.up[['c9']]) #1000
length(deg.list.up[['c6']]) #1000



# length(deg.list.up[['c7']]) #246 total
# length(deg.list.up[['c10']]) #694 total
# length(deg.list.up[['c3']]) #242


lapply(deg.list.up,FUN = function(x){ length(x)  } )
$c9 1000 $c6 1000 $c11 814 $c8 135 $c1 117 $c3 222 $c4 66 $c2 188 $c10 242

# $c8 1000 $c5 1000 $c9 894 
# $c10 694 
# $c6 137 
# $c1 145 $c3 242 
# $c2 133 $c4 123 $c7 246



#######automatic function#######

# ####read in presaved rds######
# ##p2g filter uniq modify table
# p2g.filter.uniq.modify <- readRDS('data_rds/p2g.filter.uniq.modify.rds')

# #map_cluster
# map_cluster <- readRDS('data_rds/map_cluster.rds')

# ##data rds
# exprMat.perc <- readRDS('data_rds/exprMat.perc.rds')
# exprMat.ave.z <- readRDS('data_rds/exprMat.ave.z.rds')
# peakMat.aggre <- readRDS('data_rds/peakMat.aggre.rds')

# deg.list.up <- readRDS('data_rds/deg.list.up.rds')


source('eGRN_build.r') #from modified version of pycistarget


eGRN_build(tf_select = tf_select, #modified for pycistarget with DAR as input
            outdir = outdir,
            p2g.filter.uniq.modify = p2g.filter.uniq.modify,
            map_cluster = map_cluster,
            exprMat.perc  = exprMat.perc,
            exprMat.ave.z  = exprMat.ave.z,
            peakMat.aggre = peakMat.aggre,
            deg.list.up = deg.list.up,
            cutoff.cor = cutoff.cor
          )

#tf_select_3: for ctb fusion network
1. get tf-peak-targetgene table
cluster  c1 
 tfid  TP63 
 tfid  TEAD4 
 tfid  TP53 
cluster  c7 
 tfid  TP63 
 tfid  TEAD4 
 tfid  TP53 
cluster  c9 
 tfid  RELA 
 tfid  MAFK 
cluster  c2 
 tfid  MAFK 
 tfid  MAFG 
 tfid  GCM1 
cluster  c8 
 tfid  GCM1 
cluster  c3 
 tfid  CEBPB 
 tfid  CEBPG 
 tfid  FOSL2 
 tfid  SP3 
cluster  c4 
 tfid  STAT5B 
 tfid  STAT5A 
 tfid  MITF 
 tfid  SMAD6 
cluster  c5 
 tfid  STAT5B 
 tfid  STAT5A 
 tfid  STAT4 
 tfid  MITF 
2. add data for node table
add data for node of cluster  c1  with matched rna cluster  c9 
add data for node of cluster  c7  with matched rna cluster  c6 
add data for node of cluster  c9  with matched rna cluster  c11 
add data for node of cluster  c2  with matched rna cluster  c8 
add data for node of cluster  c8  with matched rna cluster  c1 
add data for node of cluster  c3  with matched rna cluster  c3 
add data for node of cluster  c4  with matched rna cluster  c2 
add data for node of cluster  c5  with matched rna cluster  c10 
'ok to outdir: result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005/tf_peak_gene-network_useuniq/tf_select_3\n'



#tf_select_2: for mature2 network
1. get tf-peak-targetgene table
cluster  c1 
 tfid  TP63 
 tfid  TEAD4 
 tfid  TP53 
cluster  c7 
 tfid  TP63 
 tfid  TEAD4 
 tfid  TP53 
cluster  c9 
 tfid  RELA 
 tfid  GRHL3 
 tfid  ZNF217 
 tfid  MAFK 
cluster  c2 
 tfid  MAFK 
 tfid  MAFG 
 tfid  GCM1 
cluster  c8 
 tfid  GCM1 
cluster  c3 
 tfid  CEBPB 
 tfid  CEBPG 
 tfid  FOSL2 
 tfid  SP3 
cluster  c4 
 tfid  STAT5B 
 tfid  STAT5A 
 tfid  MITF 
 tfid  SMAD6 
cluster  c5 
 tfid  STAT5B 
 tfid  STAT5A 
 tfid  STAT4 
 tfid  MITF 
2. add data for node table
add data for node of cluster  c1  with matched rna cluster  c9 
add data for node of cluster  c7  with matched rna cluster  c6 
add data for node of cluster  c9  with matched rna cluster  c11 
add data for node of cluster  c2  with matched rna cluster  c8 
add data for node of cluster  c8  with matched rna cluster  c1 
add data for node of cluster  c3  with matched rna cluster  c3 
add data for node of cluster  c4  with matched rna cluster  c2 
add data for node of cluster  c5  with matched rna cluster  c10 
'ok to outdir: result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005/tf_peak_gene-network_useuniq/tf_select_2\n'


#tf_select_1: for mature 1 network
1. get tf-peak-targetgene table
cluster  c1 
 tfid  TP63 
 tfid  TEAD4 
 tfid  TP53 
cluster  c7 
 tfid  TP63 
 tfid  TEAD4 
 tfid  TP53 
cluster  c9 
 tfid  TP53 
 tfid  RELA 
 tfid  GRHL3 
 tfid  ZNF217 
cluster  c2 
 tfid  MAFK 
 tfid  MAFG 
 tfid  GCM1 
cluster  c8 
 tfid  GCM1 
cluster  c3 
 tfid  CEBPB 
 tfid  FOSL2 
 tfid  SP3 
cluster  c4 
 tfid  STAT5B 
 tfid  STAT5A 
 tfid  MITF 
 tfid  SMAD6 
cluster  c5 
 tfid  STAT5B 
 tfid  STAT5A 
 tfid  STAT4 
 tfid  MITF 
2. add data for node table
add data for node of cluster  c1  with matched rna cluster  c9 
add data for node of cluster  c7  with matched rna cluster  c6 
add data for node of cluster  c9  with matched rna cluster  c11 
add data for node of cluster  c2  with matched rna cluster  c8 
add data for node of cluster  c8  with matched rna cluster  c1 
add data for node of cluster  c3  with matched rna cluster  c3 
add data for node of cluster  c4  with matched rna cluster  c2 
add data for node of cluster  c5  with matched rna cluster  c10 
'ok to outdir: result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005/tf_peak_gene-network_useuniq/tf_select_1\n'



##tf_select
1. get tf-peak-targetgene table
cluster  c1 
 tfid  TP63 
 tfid  TEAD4 
 tfid  E2F3 
 tfid  TP53 
cluster  c7 
 tfid  TP63 
 tfid  TEAD4 
 tfid  E2F3 
 tfid  TP53 
cluster  c9 
 tfid  TP63 
 tfid  TP53 
 tfid  RELA 
 tfid  GRHL3 
 tfid  ZNF217 
 tfid  FOSL1 
 tfid  MAFK 
 tfid  MAFG 
cluster  c2 
 tfid  MAFK 
 tfid  MAFG 
 tfid  GCM1 
cluster  c8 
 tfid  GCM1 
cluster  c3 
 tfid  CEBPG 
 tfid  CEBPB 
 tfid  FOSL2 
 tfid  SP3 
cluster  c4 
 tfid  STAT5B 
 tfid  STAT4 
 tfid  STAT5A 
 tfid  ZNF554 
 tfid  MITF 
 tfid  SMAD6 
 tfid  GLIS2 
cluster  c5 
 tfid  STAT5B 
 tfid  STAT4 
 tfid  STAT5A 
 tfid  ZNF554 
 tfid  MITF 
 tfid  BACH1 
 tfid  ZBTB21 
 tfid  NFE2 
2. add data for node table
add data for node of cluster  c1  with matched rna cluster  c9 
add data for node of cluster  c7  with matched rna cluster  c6 
add data for node of cluster  c9  with matched rna cluster  c11 
add data for node of cluster  c2  with matched rna cluster  c8 
add data for node of cluster  c8  with matched rna cluster  c1 
add data for node of cluster  c3  with matched rna cluster  c3 
add data for node of cluster  c4  with matched rna cluster  c2 
add data for node of cluster  c5  with matched rna cluster  c10 
'ok to outdir: result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005/tf_peak_gene-network_useuniq/tf_select\n'




