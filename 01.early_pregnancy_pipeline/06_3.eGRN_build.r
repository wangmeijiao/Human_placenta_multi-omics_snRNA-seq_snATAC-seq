

eGRN_build <- function( tf_select = NULL,
                        outdir = NULL,
                        p2g.filter.uniq.modify = NULL,
                        map_cluster = NULL,
                        exprMat.perc  = NULL,
                        exprMat.ave.z  = NULL,
                        peakMat.aggre = NULL,
                        deg.list.up = NULL,
                        cutoff.cor = NULL ){
    ##construction of  tf-cisElement-targetGene enhancer gene regulatory network

    cat('1. get tf-peak-targetgene table\n')

    ##get tf-peak-target_gene table
    tf_peak_gene.df <- data.frame()
    for(cid in names(tf_select) ){
        cat('cluster ',cid,'\n')

        for(tfid in tf_select[[cid]]){
          cat(' tfid ', tfid, '\n')

          stopifnot(!is.null(tfmotif.list.filter[[cid]][[tfid]]))

          hp <- tfmotif.list.filter[[cid]][[tfid]][['heatmap']]
          motifid <- tfmotif.list.filter[[cid]][[tfid]][['motif']]
          nes <-  tfmotif.list.filter[[cid]][[tfid]][['NES']]

          peaks_combined <- tfmotif.list.filter[[cid]][[tfid]][['peaks_combined']]

          targetGenes_p2g_combined <- tfmotif.list.filter[[cid]][[tfid]][['targetGenes_p2g_combined']]
          targetGenes_dorc_combined <- tfmotif.list.filter[[cid]][[tfid]][['targetGenes_dorc_combined']]
          nearGene_combined <- tfmotif.list.filter[[cid]][[tfid]][['nearGenes_combined']]

          shared_targets <- tfmotif.list.filter[[cid]][[tfid]][['shared_targets']]

          targetGenes_all <- unique(targetGenes_p2g_combined)  #use all motif targets hit p2g_filter_uniq
          #targetGenes_all <- unique(targetGenes_dorc_combined) #p2g has many one-to-many peak2gene links
          #targetGenes_all <- unique(c(targetGenes_p2g_combined,targetGenes_dorc_combined,nearGene_combined))
          targetGenes_all.peaks <- sapply(strsplit(targetGenes_all,split = '\\|'), function(x){x[1]}  )
          targetGenes_all.genes <- sapply(strsplit(targetGenes_all,split = '\\|'), function(x){x[2]}  )
          targetGenes_all.peaks <- unique(targetGenes_all.peaks)
          targetGenes_all.genes <- unique(targetGenes_all.genes) 

            
          ##filter peak_combined to targetGenes_all.peaks  (for pycistarget with DAR as input)
          peaks_combined <- peaks_combined[peaks_combined %in% targetGenes_all.peaks]    
        
            
          ##for dar full then use p2g for annotation, will comment below (no need to comment out)
          stopifnot( isTRUE(all.equal(sort(targetGenes_all.peaks),sort(peaks_combined)) ) )#TRUE
        
          ##stopifnot(all.equal(sort(targetGenes_all.genes),sort(shared_targets)) )#TRUE

          #peaks_combined and shared_targets are the complete set

          ##choose targetGenes_all (if combined equal to peaks_combined and shared_targets, because use dar_overlap_p2g)

          if( isTRUE(all.equal(sort(targetGenes_all.peaks) ,sort(peaks_combined)) ) ){ 
              #for(i in targetGenes_all){    
                  targetGenes_all.df <- as.data.frame(Reduce(rbind,strsplit(x = targetGenes_all, split = "\\|")),stringsAsFactors = FALSE)
                  rownames(targetGenes_all.df) <- NULL
                  colnames(targetGenes_all.df) <- c('peak','target')
                  targetGenes_all.df$tf <- tfid
                  targetGenes_all.df$dar <- cid
                  targetGenes_all.df$NES <- nes

                  #get peak2gene correlation from p2g.filter.uniq table
                  p2gid <- paste0(targetGenes_all.df$peak,'|',targetGenes_all.df$target)
                  stopifnot(all.equal(targetGenes_all,p2gid))
                  stopifnot(sum(p2gid %in% paste0(p2g.filter.uniq.modify$peak,'|',p2g.filter.uniq.modify$gene) ) == length(p2gid) ) 

                  cor.df <- p2g.filter.uniq.modify[p2gid,c('peak','Correlation')]
                  stopifnot(all.equal(targetGenes_all.df$peak, cor.df$peak))
                  cor.value <- cor.df$Correlation
                  stopifnot(sum(is.na(cor.value)) == 0)
                  targetGenes_all.df$cor <- cor.value

                  stopifnot(sum(cor.value <= 0) == 0)
              
                  ##add additional filter layers for network with high density
                  if(cutoff.cor[[cid]] == 'auto'){
                      qs <- quantile(x = targetGenes_all.df$cor, probs = seq(from = 0, to =1, by = 0.1) )
                      cutoff <- qs['25%']
                      targetGenes_all.df <- subset(targetGenes_all.df,cor >= cutoff)
                      
                  }else{
                    len.filtered <- sum(targetGenes_all.df$cor < cutoff.cor[[cid]])
                    len.total <- length(targetGenes_all.df$cor)
                    if(len.filtered != 0){
                      cat('   filter peak-to-gene with cor cutoff: ',cutoff.cor[[cid]],', ',len.filtered,' of total ',len.total,' filtered \n',sep='')
                    }else{  }
                      
                    targetGenes_all.df <- subset(targetGenes_all.df,cor >= cutoff.cor[[cid]])
                  }
                  stopifnot(nrow(targetGenes_all.df) != 0 )
              
                  #edge table
                  targetGenes_all.df <- targetGenes_all.df[,c('dar','tf','peak','target','NES','cor')]
              #}
          }else{
              stop("targetGenes_p2g_combined not eq peak_combined, will stop \n")
              
              #use peaks_combined as the whole set, add target gene as 'NA'
              ##cat("targetGenes_p2g_combined not eq peak_combined, use peak_combined and fill target gene as NA \n")
              
              #use targetGenes_p2g_combined to get edge table and discard peaks without p2g overlap?
              
#               peaks.df <- data.frame(peak = peaks_combined, stringsAsFactors = FALSE)
#               temp.df <- as.data.frame(Reduce(rbind,strsplit(x = targetGenes_all, split = "\\|")),stringsAsFactors = FALSE)
#               rownames( temp.df ) <- NULL
#               colnames( temp.df ) <- c('peak','target')
#                temp.df$tf <- tfid
#                temp.df$dar <- cid
#               stopifnot(  sum(temp.df$peak %in% peaks_combined) == nrow(temp.df) )
#               #merge
#               targetGenes_all.df <- base::merge(x = peaks.df, y = temp.df, by = 'peak', all.x = TRUE)

#               targetGenes_all.df$tf <- tfid
#               targetGenes_all.df$dar <- cid
#               #targetGenes_all.df$NES <- nes

#               targetGenes_all.df <- targetGenes_all.df[,c('dar','tf','peak','target','NES','cor')]
          }
            
            
          tf_peak_gene.df <- rbind.data.frame(tf_peak_gene.df,targetGenes_all.df)


        }
    }

    saveRDS(tf_peak_gene.df,paste(outdir,'/tf_peak_gene.df.motifcombined.p2guniq.rds',sep='') )
    write.table(tf_peak_gene.df,file = paste(outdir,'/tf_peak_gene.df.motifcombined.p2guniq.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)



    ##split into pairwise format (edge table format?) as tf-peak (cistrome), peak-targetgene (regulon)

    cistrome.df <- tf_peak_gene.df[,c('dar','tf','peak','NES')]
    cistrome.df$type <- 'cistrome'
    colnames(cistrome.df) <- c('dar','name1','name2','score','type')


    regulon.df <- tf_peak_gene.df[,c('dar','peak','target','cor')]
    #regulon.df <- tf_peak_gene.df[,c('dar','peak','target_p2g_combined')]
    regulon.df$type <- 'regulon'
    colnames(regulon.df) <- c('dar','name1','name2','score','type')


    network.table.df <- rbind.data.frame(cistrome.df,regulon.df)


    #########dedup#########

    flag1 <- duplicated(paste(network.table.df$dar,network.table.df$name1,network.table.df$name2,sep='|') )
    #sum(flag1) #1081, only filter within dar duplicated edges

    flag2 <- duplicated(paste(network.table.df$name1,network.table.df$name2,sep='|') )
    #sum(flag2) #2355, filter all duplicated edges 

    flag3 <- flag1 == FALSE & flag2 == TRUE
    #sum(flag3) #1274 duplicated between dars (keep, c1 c2 very few!)

    #sum(flag1)+sum(flag3) == sum(flag2)

    flag <- flag1 #only remove within dar duplicated edge, keep among dar edge duplications
    #network.table.df[flag,]


    network.table.df.dedup <- network.table.df[!flag,]
 

    ##split score to NES and Cor, fill with NA (cytoscape will fail to read in ,use 0)

    nes_value <- ifelse(network.table.df.dedup$type =='cistrome', network.table.df.dedup$score ,0)
    cor_value <- ifelse(network.table.df.dedup$type =='regulon', network.table.df.dedup$score ,0)

    network.table.df.dedup$NES <- sprintf("%.5f",as.numeric(nes_value))
    network.table.df.dedup$Cor <- sprintf("%.5f",as.numeric(cor_value)) 

    network.table.df.dedup$score <- NULL
    #network.table.df.dedup$test <- 1


    ##save edge table to file


    write.table(network.table.df.dedup,file = paste(outdir,'/network.table.df.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)

#     write.table(subset(network.table.df.dedup,dar=='c2'),file = paste(outdir,'/network.table.df.c2.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)

#     write.table(subset(network.table.df.dedup,dar=='c1'),file = paste(outdir,'/network.table.df.c1.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)

#     write.table(subset(network.table.df.dedup,dar=='c9'),file = paste(outdir,'/network.table.df.c9.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)

#     write.table(subset(network.table.df.dedup,dar=='c5'),file = paste(outdir,'/network.table.df.c5.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)

#     write.table(subset(network.table.df.dedup,dar=='c6'),file = paste(outdir,'/network.table.df.c6.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)

#     write.table(subset(network.table.df.dedup,dar=='c3'),file = paste(outdir,'/network.table.df.c3.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)

    for(i in unique(network.table.df.dedup$dar) ){
        
        write.table(subset(network.table.df.dedup,dar==i),file = paste(outdir,'/network.table.df.',i,'.txt',sep=''),col.names = TRUE, row.names = FALSE, sep='\t',quote=FALSE)
        
    }
    
    
    
    cat('2. add data for node table\n')

    #############add data (of this cluster) for node: tf (expression <and FC >), peak (accessibility) , gene (expression and percentage)################



    #create node attribution table
    ##start to fill node data
    tf_peak_gene.node.list <- list()

    for(cid in names(tf_select)){

        stopifnot(!is.null(map_cluster[[cid]]))
        cid_map <- map_cluster[[cid]]

        cat ('add data for node of cluster ', cid, ' with matched rna cluster ',cid_map,'\n')
        tf_peak_gene.df.sel <- subset(tf_peak_gene.df, dar == cid)

        node.df <- rbind.data.frame(
            unique(data.frame(name = tf_peak_gene.df.sel$tf,type = 'TF',stringsAsFactors = FALSE) ),
            unique(data.frame(name = tf_peak_gene.df.sel$peak,type = 'region',stringsAsFactors = FALSE) ),
            unique(data.frame(name = tf_peak_gene.df.sel$target,type = 'gene',stringsAsFactors = FALSE) ),
            #unique(data.frame(name = tf_peak_gene.df.sel$target_p2g_combined,type = 'gene',stringsAsFactors = FALSE) ),
            stringsAsFactors = FALSE
        )

        ##add tf gene expression, and percentage

        sum(duplicated(node.df$name)) != 0 #target gene or peak may duplicated 
        #node.df[which(duplicated(node.df$name)),]

        ##isDEG or isShow ?

        deg <- deg.list.up[[cid_map]]
        deg.sel <- deg[node.df$name ]

        deg.label <- names(deg.sel)

        deg.label[is.na(deg.label)] <- ''
        deg.label[which(node.df$type == 'TF')] <- node.df$name[which(node.df$type == 'TF')]

        deg.log2FC <- deg.sel
        isDEG <- ifelse(deg.label == '','No','Yes')
        names(deg.log2FC) <- NULL

        deg.sel.df <- data.frame('isDEG'=isDEG,'deg_label'=deg.label,'log2FC'=deg.log2FC, stringsAsFactors = FALSE)

        ##exprPerc
       exprMat.perc.sel <- subset(exprMat.perc, id == cid_map)
        stopifnot( sum(duplicated(exprMat.perc.sel$features.plot)) == 0 )
        rownames(exprMat.perc.sel) <- exprMat.perc.sel$features.plot

        expr_perc <- exprMat.perc.sel[node.df$name,c('features.plot','pct.exp')]

        ##exprValue

        #expr_value <-  exprMat.ave.z[node.df$name,cid,drop = FALSE]
        #expr_value$geneid <- rownames(expr_value) #will add .1 if duplicated
        expr_value <-  exprMat.ave.z[node.df$name,c('geneid',cid_map)]
        #expr_value$geneid <- rownames(expr_value)
        colnames(expr_value) <- c('geneid','value')
        #expr_value <- expr_value[,c(2,1)]

        #node.df.data <- cbind.data.frame(node.df,expr_perc,expr_value)
        node.df.data <- cbind.data.frame(node.df,deg.sel.df,expr_perc,expr_value)

        stopifnot(all.equal(node.df.data[!is.na(node.df.data$log2FC),'name'],node.df.data[!is.na(node.df.data$log2FC),'deg_label']) )
        stopifnot(all.equal(node.df.data[!is.na(node.df.data$pct.exp),'name'],node.df.data[!is.na(node.df.data$pct.exp),'features.plot']) ) #TRUE
        stopifnot(all.equal(node.df.data[!is.na(node.df.data$value),'name'],node.df.data[!is.na(node.df.data$value),'geneid']) )#TRUE

        ##add peak accessibility

        #acc_value <- peakMat.aggre[node.df$name,cid,drop = FALSE]
        #acc_value$peakid <- rownames(acc_value)
        acc_value <- peakMat.aggre[node.df$name,c('peakid',cid)]
        colnames(acc_value) <- c('peak','acc_value')
        #acc_value <- acc_value[,c(2,1)]

        node.df.data <- cbind.data.frame(node.df.data,acc_value)

        stopifnot(all.equal(node.df.data[!is.na(node.df.data$acc_value),'name'],node.df.data[!is.na(node.df.data$acc_value),'peak']) )

        colnames(node.df.data) <- c('name',	'type', 'isDEG', 'deg_label','log2FC','label',	'pct.exp',	'geneid', 'expr.value',	'peakid','acc_value')

    #     #NA to 0 for value column
    #     node.df.data$pct.exp <-  ifelse(is.na(node.df.data$pct.exp), 0  , node.df.data$pct.exp)
    #     node.df.data$expr.value <-  ifelse(is.na(node.df.data$pct.exp), 0  , node.df.data$pct.exp)
    #     node.df.data$acc_value <-  ifelse(is.na(node.df.data$pct.exp), 0  , node.df.data$pct.exp)

        tf_peak_gene.node.list[[cid]] <- node.df.data
    }

 
    saveRDS(tf_peak_gene.node.list,paste(outdir,'/tf_peak_gene.node.list.rds',sep='') )


    for(i in names(tf_peak_gene.node.list) ){
        write.table(tf_peak_gene.node.list[[i]],file=paste(outdir,'/tf_peak_gene.node.',i,'.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)

    }
    
  
#     write.table(tf_peak_gene.node.list[['c2']],file=paste(outdir,'/tf_peak_gene.node.c2.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)

#     write.table(tf_peak_gene.node.list[['c1']],file=paste(outdir,'/tf_peak_gene.node.c1.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)



#     write.table(tf_peak_gene.node.list[['c9']],file=paste(outdir,'/tf_peak_gene.node.c9.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)



#     write.table(tf_peak_gene.node.list[['c5']],file=paste(outdir,'/tf_peak_gene.node.c5.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)



#     write.table(tf_peak_gene.node.list[['c6']],file=paste(outdir,'/tf_peak_gene.node.c6.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)


#     write.table(tf_peak_gene.node.list[['c3']],file=paste(outdir,'/tf_peak_gene.node.c3.txt',sep=''),sep='\t',col.names = TRUE, row.names = FALSE, quote = FALSE)

  return(paste('ok to outdir: ',outdir,'\n',sep=''))


}
