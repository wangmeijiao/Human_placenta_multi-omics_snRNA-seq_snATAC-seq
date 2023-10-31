#https://github.com/r3fang/SnapATAC/blob/master/examples/10X_snATAC/README.md#viz

####install SnapATAC
#install.packages("plot3D")
#library(devtools)
#install_github("r3fang/SnapATAC")


#rhdf5, GenomicRanges, IRanges, edgeR, Rhdf5lib, BiocGenerics, S4Vectors, GenomeInfoDb, XVector, limma, GenomeInfoDbData, zlibbioc
library('SnapATAC') #1.0.0, snaptools(py3.7) 1.4.8
library(ggplot2)
library('GenomicRanges')
library(magrittr)

library(leiden)

library(viridis)

library(hrbrthemes)

library(tictoc)

library(patchwork)
library(grid)
library(gridExtra)

#compare plot many cluster results
library(clustree)

##for Silhoutte coefficent
library(bluster)


##for Dunn index
library(fpc)
library(NbClust)

##for SnapATAC chromVAR

library('JASPAR2016') #The Jaspar motif database 2016
library('JASPAR2020') #The Jaspar motif database 2020
#devtools::install_github("GreenleafLab/chromVARmotifs")

#https://github.com/GreenleafLab/chromVARmotifs
library('chromVARmotifs') #Greenleaf lab motifs, 
data("human_pwms_v1") #1764 TFs

##get version 2, the nonredundant, 870 TFs

sout <- sapply(strsplit(names(human_pwms_v1), split = "_"), function(s) {c(s[3])}  )
human_pwms_v2 <- human_pwms_v1[match(unique(sout), sout)]

#save(human_pwms_v2, file = 'human_pwms_v2.rda')




############Buen colors########
library(BuenColors)
color_set0 <- jdb_color_maps #17 different colors
names(color_set0) <- NULL
#plot(1:17,1:17,pch = 19, cex = 5,col=jdb_color_maps)

#discrete colors
color_set1 <- jdb_palette("solar_extra") #9 discrete but gradient colors
color_set2 <- jdb_palette("brewer_spectra") #9 discrete but gradient colors
color_set3 <- jdb_palette("flame_light") #9 discrete but gradient colors, good!

color_set3_ext12 <- colorRampPalette(colors = as.character(color_set3))(12)
color_set3_ext17 <- colorRampPalette(colors = as.character(color_set3))(17)

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


#########customized colors########
color_snap = c('1'='grey','2'='#E31A1C','3'='#FFD700','4'='#771122','5'='#777711','6'='#1F78B4','7'='#68228B','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')
#names(color_snap) <- NULL

#modified for CTB with dark red colors
color_snap_mod1 = c('1'='#777711','2'='#E31A1C','3'='#68228B','4'='#771122','5'='grey','6'='#1F78B4','7'='#FFD700','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')

color_signac = c(
'0'='#E6D55E','1'='#792B8A','2'='#DA703D','3'='#9DC8E5','4'='#BA273C','5'='#C2C184','6'='#7F8084','7'='#65AB53','8'='#D082AF','9'='#496EAB','10'='#DE896D','11'='#491F8B','12'='#E1AD49','13'='#8E1B85','14'='#E7EE77','15'='#7D1A1D','16'='#96B355')
names(color_signac) <- NULL

color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" )

color_good <- c("#E7D654", "#6F1482" ,"#DC7035", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,"#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,
 "#CA362E" ,"#2B3918","#1E1E1E" )




color_gradient_my <- c(
    rgb(5,48,97,maxColorValue = 255),
    rgb(42,113,178,maxColorValue = 255),
    rgb(147,198,222,maxColorValue = 255),
    rgb(239,243,245,maxColorValue = 255),
    rgb(253,219,199,maxColorValue = 255),
    rgb(214,96,77,maxColorValue = 255),
    rgb(121,5,34,maxColorValue = 255)

)

color_cellranger <-c('#820610','#C50F1E','#F42428','#F86D30','#FBB33D','#FCFB4E','#C0FB61','#87FB8F','#41FEF9','#2BAED7','#155CB1','#08238D')





reds <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf')
blues <- c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')

map_cellcolor_atac <- list(
    '4'=reds[5], #STB1
    '2'=reds[6],  #STB-naive
    '6'=reds[3],  #STB3
    '1'=reds[2], #STB4
    '5'=reds[4], #STB2
    '3'=reds[1], #Syncytial knot
    '7'='#694d9f', #STB-new
    '8'='darkgreen', #CTB
    '9' = '#92c5de'
 )


map_cellcolor_atac <- list( #similar to rna color palette, use this?
  '7' = '#67001f',
 '3'= '#b2182b',
 '8'=  'black',#'grey35',#'darkgreen',#'black',
 '6'= '#f4a582',
 '4'= '#FB8D3C',#'#e27f36',#'#c97130',#'#b0632a',#'#e27f36',#'#FB8D3C',
 '1'= '#d6604d',
 '2'= 'darkred',
 '9'= '#92c5de',
 '5'= '#fedcbd'#'#4393c3',
 #'11'= 'darkgreen',#'#2166ac',
 #'9'= '#053061'
)


map_cellcolor_atac <- list( #equal to color_good
    
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

 )



cluster.df.add.flipy <- cluster.df.add

cluster.df.add.flipy$UMAP_2 <- -1* cluster.df.add.flipy$UMAP_2

quickDimPlot_labelon(data = cluster.df.add.flipy, feature = 'cluster', color_use = map_cellcolor_atac, title= 'late combined-ATAC', shrink.x = .3, shrink.y = .1,shuffle = FALSE,pt.size = .3,height = 7.5, width = 8.5)

ggsave('pdfs/PLA-late-combine.cstb.hotcolor.pdf',height=7.5,width=8.5,useDingbats=FALSE)



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

darker (colorset = c('#FB8D3C'), ratio = 0.1) #'#e27f36'#'#b0632a'




###select one global color set###
color <- color_good


sample <- "PLA-late-combine-ATAC"



###########step 1, create snap object ##########

#ls ../../placenta_10X_late{1..6}/02.snapATAC/snap_from_10xbam/*.snap 


snap.files = c(
    
'../../placenta_10X_late1/02.snapATAC/snap_from_10xbam/placenta-late1-atac.snap',
'../../placenta_10X_late2/02.snapATAC/snap_from_10xbam/placenta-late2-atac.snap',
'../../placenta_10X_late3/02.snapATAC/snap_from_10xbam/placenta-late3-atac.snap',
#'../../placenta_10X_late4/02.snapATAC/snap_from_10xbam/placenta-late4-atac.snap',
'../../placenta_10X_late5/02.snapATAC/snap_from_10xbam/placenta-late5-atac.snap',
'../../placenta_10X_late6/02.snapATAC/snap_from_10xbam/placenta-late6-atac.snap',
'../../placenta_10X_late9/02.snapATAC/snap_from_10xbam/placenta-late9-atac.snap'
    #"placenta-term-atac-1.snap", 
    #"placenta-term-atac-2.snap"
  )



#snap.file = "snap_from_10xbam/placenta-term-atac-1.snap"

sample.names = c(
    "placenta_donor1",
    "placenta_donor2",
    "placenta_donor3",
    #"placenta_donor4",
    "placenta_donor5",
    "placenta_donor6",
    "placenta_donor9"
  );
#prefix = c('donor1','donor2')

#sample.name = 'placenta-term-atac-1'

#ls ../../placenta_10X_late{1..6}/01.data_cellranger_atac/singlecell.csv

barcode.files = c(
    #"/home/mjwang/pwdex/placenta_10X_term1/01.data_cellranger_atac/ATAC_term1/singlecell.csv",
    #"/home/mjwang/pwdex/placenta_10X_term2/01.data_cellranger_atac/ATAC_term2/singlecell.csv",
    '../../placenta_10X_late1/01.data_cellranger_atac/singlecell.csv',
    '../../placenta_10X_late2/01.data_cellranger_atac/singlecell.csv',
    '../../placenta_10X_late3/01.data_cellranger_atac/promotor_ratio_singlecell/singlecell.add_promoter_frag_count.csv',
    #'../../placenta_10X_late4/01.data_cellranger_atac/singlecell.csv',
    '../../placenta_10X_late5/01.data_cellranger_atac/promotor_ratio_singlecell/singlecell.add_promoter_frag_count.csv',
    '../../placenta_10X_late6/01.data_cellranger_atac/promotor_ratio_singlecell/singlecell.add_promoter_frag_count.csv',
    '../../placenta_10X_late9/01.data_cellranger_atac/promotor_ratio_singlecell/singlecell.add_promoter_frag_count.csv'
  );

#barcode.file = '../../01.data_cellranger_atac/PLA-term-ATAC-1/singlecell.csv'


x.sp.ls = lapply(seq(snap.files), function(i){
    createSnap(
        file=snap.files[i],
        sample=sample.names[i]
    );
  })

names(x.sp.ls) = sample.names

#x.sp = createSnap(file=snap.file,sample=sample.name) ##20000 barcodes from snap file without filtering
#names(x.sp) = sample.name

barcode.ls = lapply(seq(snap.files), function(i){
    barcodes = read.csv(
        barcode.files[i], 
        head=TRUE,
        stringsAsFactors = FALSE
    )
    # remove NO BAROCDE line
    barcodes = barcodes[2:nrow(barcodes),];
    barcodes$logUMI = log10(barcodes$passed_filters + 1);
    barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
    barcodes$peak_ratio = (barcodes$peak_region_fragments+1) / (barcodes$passed_filters + 1)
    #barcodes$barcode <- gsub(pattern = '[0-9]+$',replacement = i, x = barcodes$barcode), no need to do this, since the sample slot can do 
    barcodes
  })

names(barcode.ls) = sample.names



# barcode.ls[['placenta_donor1']][,'sample'] <- 'D1'
# barcode.ls[['placenta_donor2']][,'sample'] <- 'D2'
# barcode.ls[['placenta_donor3']][,'sample'] <- 'D3'
# barcode.ls[['placenta_donor4']][,'sample'] <- 'D4'
# barcode.ls[['placenta_donor5']][,'sample'] <- 'D5'
# barcode.ls[['placenta_donor6']][,'sample'] <- 'D6'


#barcode.ls[['placenta_donor2']][,'barcode'] <- gsub(pattern = "-1$",replacement = '-2',x = barcode.ls[['placenta_donor2']][,'barcode'])


# barcodes <- rbind.data.frame(barcode.ls[['placenta_donor1']],barcode.ls[['placenta_donor2']])

# barcodes <- Reduce(f = rbind.data.frame, x = barcode.ls)

# barcodes$sample <- factor(barcodes$sample,levels=c('D1','D2'))

# table(barcodes$sample)


#look for barcode file length

lapply(seq(barcode.files), function(i){ nrow(barcode.ls[[i]])   }  )
563137
548114
430703
433944
403010
408249

# 563137
# 548114
# 430703
# 397189
# 433944
# 403010

# 657285
# 572249
# 425371
# 430167
# 415602
# 397870


barcode.ls[['placenta_donor1']][,'sample'] <- 'D1' #22 column
barcode.ls[['placenta_donor2']][,'sample'] <- 'D2' #22 column
barcode.ls[['placenta_donor3']][,'sample'] <- 'D3' #23 column
#barcode.ls[['placenta_donor4']][,'sample'] <- 'D4' #23 column
barcode.ls[['placenta_donor5']][,'sample'] <- 'D5' #23 column
barcode.ls[['placenta_donor6']][,'sample'] <- 'D6' #23 column
barcode.ls[['placenta_donor9']][,'sample'] <- 'D9' #23 column

barcode.ls[['placenta_donor1']][,'library'] <- 'placenta_10X_late1'
barcode.ls[['placenta_donor2']][,'library'] <- 'placenta_10X_late2'
barcode.ls[['placenta_donor3']][,'library'] <- 'placenta_10X_late3'
#barcode.ls[['placenta_donor4']][,'library'] <- 'placenta_10X_late4'
barcode.ls[['placenta_donor5']][,'library'] <- 'placenta_10X_late5'
barcode.ls[['placenta_donor6']][,'library'] <- 'placenta_10X_late6'
barcode.ls[['placenta_donor9']][,'library'] <- 'placenta_10X_late9'



##modify barcode id
barcode.ls <- lapply(barcode.ls, function(x){ x$barcode.add = paste(x$library,x$barcode,sep=':') ;x  }
      
      )



#barcode.ls[['placenta_donor1']]$promoter_region_fragments
#barcode.ls[['placenta_donor3']]$promoter_region_fragments


##get shared column name

colid1 <- colnames(barcode.ls[['placenta_donor1']])
colid2 <- colnames(barcode.ls[['placenta_donor2']])
colid3 <- colnames(barcode.ls[['placenta_donor3']])
#colid4 <- colnames(barcode.ls[['placenta_donor4']])
colid5 <- colnames(barcode.ls[['placenta_donor5']])
colid6 <- colnames(barcode.ls[['placenta_donor6']])
colid9 <- colnames(barcode.ls[['placenta_donor9']])


all.equal(colid1,colid2) #TRUE
#all.equal(colid3,colid4) #TRUE
all.equal(colid5,colid6) #TRUE
all.equal(colid3,colid5) #TRUE
all.equal(colid6,colid9) #TRUE

all.equal(colid1,colid3) #one diff: 'excluded_reason' 'nonprimary'
'barcode''total''duplicate''chimeric''unmapped''lowmapq''mitochondrial''passed_filters''cell_id''is__cell_barcode''TSS_fragments''DNase_sensitive_region_fragments''enhancer_region_fragments''promoter_region_fragments''on_target_fragments''blacklist_region_fragments''peak_region_fragments''peak_region_cutsites''logUMI''promoter_ratio''peak_ratio''sample''library''barcode.add'

'barcode''total''duplicate''chimeric''unmapped''lowmapq''mitochondrial''nonprimary''passed_filters''is__cell_barcode''excluded_reason''TSS_fragments''DNase_sensitive_region_fragments''enhancer_region_fragments''promoter_region_fragments''on_target_fragments''blacklist_region_fragments''peak_region_fragments''peak_region_cutsites''logUMI''promoter_ratio''peak_ratio''sample''library''barcode.add'


# 'barcode''total''duplicate''chimeric''unmapped''lowmapq''mitochondrial''passed_filters''cell_id''is__cell_barcode''TSS_fragments''DNase_sensitive_region_fragments''enhancer_region_fragments''promoter_region_fragments''on_target_fragments''blacklist_region_fragments''peak_region_fragments''peak_region_cutsites''logUMI''promoter_ratio''peak_ratio''sample'

# 'barcode''total''duplicate''chimeric''unmapped''lowmapq''mitochondrial''nonprimary''passed_filters''is__cell_barcode''excluded_reason''TSS_fragments''DNase_sensitive_region_fragments''enhancer_region_fragments''promoter_region_fragments''on_target_fragments''blacklist_region_fragments''peak_region_fragments''peak_region_cutsites''logUMI''promoter_ratio''peak_ratio''sample'

colid <- intersect(colid1,colid3)

'barcode''total''duplicate''chimeric''unmapped''lowmapq''mitochondrial''passed_filters''is__cell_barcode''TSS_fragments''DNase_sensitive_region_fragments''enhancer_region_fragments''promoter_region_fragments''on_target_fragments''blacklist_region_fragments''peak_region_fragments''peak_region_cutsites''logUMI''promoter_ratio''peak_ratio''sample''library''barcode.add'

# 'barcode''total''duplicate''chimeric''unmapped''lowmapq''mitochondrial''passed_filters''is__cell_barcode''TSS_fragments''DNase_sensitive_region_fragments''enhancer_region_fragments''promoter_region_fragments''on_target_fragments''blacklist_region_fragments''peak_region_fragments''peak_region_cutsites''logUMI''promoter_ratio''peak_ratio''sample'

barcodes <- rbind.data.frame(barcode.ls[['placenta_donor1']][,colid],
                             barcode.ls[['placenta_donor2']][,colid],
                             barcode.ls[['placenta_donor3']][,colid],
                             #barcode.ls[['placenta_donor4']][,colid],
                             barcode.ls[['placenta_donor5']][,colid],
                             barcode.ls[['placenta_donor6']][,colid],
                             barcode.ls[['placenta_donor9']][,colid]
                            )

2787157 × 23
#2776097 × 23


barcodes$sample <- factor(barcodes$sample,levels=c('D1','D2','D3','D5','D6','D9'))

table(barcodes$sample)
    D1     D2     D3     D5     D6     D9 
563137 548114 430703 433944 403010 408249

#     D1     D2     D3     D4     D5     D6 
# 563137 548114 430703 397189 433944 403010 

#  D1     D2     D3     D4     D5     D6 
# 657285 572249 425371 430167 415602 397870 


saveRDS(barcodes,'barcodes.rds')
barcodes <- readRDS('barcodes.rds')



##boxplot barcode table column



par(mfrow = c(2,3))
boxplot(subset(barcodes, sample == 'D1' )$passed_filters )
boxplot(subset(barcodes, sample == 'D2' )$passed_filters )
boxplot(subset(barcodes, sample == 'D3' )$passed_filters )
#boxplot(subset(barcodes, sample == 'D4' )$passed_filters )
boxplot(subset(barcodes, sample == 'D5' )$passed_filters )
boxplot(subset(barcodes, sample == 'D6' )$passed_filters )
boxplot(subset(barcodes, sample == 'D9' )$passed_filters )


par(mfrow = c(2,3))
boxplot(subset(barcodes, sample == 'D1' )$promoter_region_fragments ) 
boxplot(subset(barcodes, sample == 'D2' )$promoter_region_fragments )
boxplot(subset(barcodes, sample == 'D3' )$promoter_region_fragments ) #manually added #0, because cellranger-atac-2 use arc no promoter.bed
#boxplot(subset(barcodes, sample == 'D4' )$promoter_region_fragments ) #manually added#0
boxplot(subset(barcodes, sample == 'D5' )$promoter_region_fragments ) #manually added#0
boxplot(subset(barcodes, sample == 'D6' )$promoter_region_fragments ) #manually added#0
boxplot(subset(barcodes, sample == 'D9' )$promoter_region_fragments )#manually added


par(mfrow = c(2,3))
boxplot(subset(barcodes, sample == 'D1' )$TSS_fragments )
boxplot(subset(barcodes, sample == 'D2' )$TSS_fragments )
boxplot(subset(barcodes, sample == 'D3' )$TSS_fragments )
#boxplot(subset(barcodes, sample == 'D4' )$TSS_fragments )
boxplot(subset(barcodes, sample == 'D5' )$TSS_fragments )
boxplot(subset(barcodes, sample == 'D6' )$TSS_fragments )
boxplot(subset(barcodes, sample == 'D9' )$TSS_fragments )


par(mfrow = c(2,3))
boxplot(subset(barcodes, sample == 'D1' )$logUMI )
boxplot(subset(barcodes, sample == 'D2' )$logUMI )
boxplot(subset(barcodes, sample == 'D3' )$logUMI )
#boxplot(subset(barcodes, sample == 'D4' )$logUMI )
boxplot(subset(barcodes, sample == 'D5' )$logUMI )
boxplot(subset(barcodes, sample == 'D6' )$logUMI )
boxplot(subset(barcodes, sample == 'D9' )$logUMI )

par(mfrow = c(2,3))
boxplot(subset(barcodes, sample == 'D1' )$promoter_ratio ) #
boxplot(subset(barcodes, sample == 'D2' )$promoter_ratio )
boxplot(subset(barcodes, sample == 'D3' )$promoter_ratio )
#boxplot(subset(barcodes, sample == 'D4' )$promoter_ratio )
boxplot(subset(barcodes, sample == 'D5' )$promoter_ratio )
boxplot(subset(barcodes, sample == 'D6' )$promoter_ratio )
boxplot(subset(barcodes, sample == 'D9' )$promoter_ratio )

par(mfrow = c(2,3))
boxplot(subset(barcodes, sample == 'D1' )$peak_ratio ) #
boxplot(subset(barcodes, sample == 'D2' )$peak_ratio )
boxplot(subset(barcodes, sample == 'D3' )$peak_ratio )
#boxplot(subset(barcodes, sample == 'D4' )$peak_ratio )
boxplot(subset(barcodes, sample == 'D5' )$peak_ratio )
boxplot(subset(barcodes, sample == 'D6' )$peak_ratio )
boxplot(subset(barcodes, sample == 'D9' )$peak_ratio )




######step 2. select barcode in one step#############
sample.name <- 'donor1 to donor6'


barcodes

ggplot(barcodes, aes(x=logUMI, y=promoter_ratio)) + 
    geom_point(size=0.1, col="grey") +
    theme_classic()	+
    ggtitle(sample.name) +
    ylim(0, 1) + xlim(0, 6) + 
    labs(x = "log10(UMI)", y="promoter ratio")
ggplot(barcodes, aes(x=logUMI, y=promoter_ratio,col=sample)) + 
    geom_point(size=0.1) +
    theme_classic()	+
    ggtitle(sample.name) +
    ylim(0, 1) + xlim(0, 6) + 
    labs(x = "log10(UMI)", y="promoter ratio")

ggplot(barcodes, aes(x=logUMI, y=peak_ratio)) + 
    geom_point(size=0.1, col="grey") +
    theme_classic()	+
    ggtitle(sample.name) +
    ylim(0, 1) + xlim(0, 6) + 
    labs(x = "log10(UMI)", y="peak ratio")    
ggplot(barcodes, aes(x=logUMI, y=peak_ratio,col=sample)) + 
    geom_point(size=0.1,alpha=0.2) +
    theme_classic()	+
    ggtitle(sample.name) +
    ylim(0, 1) + xlim(0, 6) + 
    labs(x = "log10(UMI)", y="peak ratio")
     

cutoff.logUMI.low = 4 # 10k
#cutoff.logUMI.low = 3.5 #3162
cutoff.logUMI.high = 5 #100k  (will filter >50,000 later)
cutoff.FRIP.low = 0.2 #0.15 #promoter ratio min
cutoff.FRIP.high = 0.5 #promoter ratio max
#cutoff.FRIP.low = 0.2
#cutoff.FRIP.high = 1

#cutoff.FRIP.low = 0.65#peak ratio min
#cutoff.FRIP.high = 0.9 #peak ratio max


idx = which(
    barcodes$logUMI >= cutoff.logUMI.low & 
    barcodes$logUMI <= cutoff.logUMI.high & 
    barcodes$promoter_ratio >= cutoff.FRIP.low &
    barcodes$promoter_ratio <= cutoff.FRIP.high
    #barcodes$peak_ratio >= cutoff.FRIP.low &
    #barcodes$peak_ratio <= cutoff.FRIP.high
) #51787 #10951 #20938 #10951 #8107 #6618

#peak ratio
#3.3-5 (1995 - 100000), 0.25 - 0.9: 62596
#3.3-5 (1995 - 100000), 0.5 - 0.9: 57908
#3.5-5 (3162 - 100000), 0.6 - 0.9: 57908
#3.5-5 (3162 - 100000), 0.6 - 0.9: 40570


#promoter ratio (no use)
#3.5-5 (3162 - 100000), 0.2-0.5: 20938
#3.3-5 (1995 - 100000), 0.2-0.5: 28386 

ggplot(barcodes, aes(x=logUMI, y=promoter_ratio)) + 
    geom_point(size=0.1, col="grey") +
    geom_point(data=barcodes[idx,], aes(x=logUMI, y=promoter_ratio),size=0.3,col='red'  ) +
    theme_classic()	+
    ggtitle( paste0(sample.name,' promoter_ratio plot, cutoff.logUMI: ',cutoff.logUMI.low,'-',cutoff.logUMI.high, '\ncutoff.FRIP: ',cutoff.FRIP.low, '-',cutoff.FRIP.high, ' total selected:',length(idx))) +
    ylim(0, 1) + xlim(0, 6) + 
    labs(x = "log10(UMI)", y="promoter ratio")

ggplot(barcodes, aes(x=logUMI, y=peak_ratio)) + 
    geom_point(size=0.1, col="grey") +
    geom_point(data=barcodes[idx,], aes(x=logUMI, y=peak_ratio),size=0.3,col='red'  ) +
    theme_classic()	+
    ggtitle(paste0(sample.name,' peak_ratio plot, cutoff.logUMI: ',cutoff.logUMI.low,'-',cutoff.logUMI.high, '\ncutoff.FRIP: ',cutoff.FRIP.low, '-',cutoff.FRIP.high, ' total selected:',length(idx))) +
    ylim(0, 1) + xlim(0, 6) + 
    labs(x = "log10(UMI)", y="peak ratio")
        


barcodes.bk <- barcodes
barcodes = barcodes[idx,]
#barcodes$barcode <- droplevels(barcodes$barcode)

table(barcodes$sample)

 D1    D2    D3    D5    D6    D9  #3.5 - 5; 0.2 - 0.5 51787
11377  9561  6047  7089  9237  8476 


  D1   D2   D3   D5   D6   D9  #4.0 - 5.0; 0.2 - 0.5 31684, use this?
6508 4443 5079 4105 5785 5764

#  D1    D2    D3    D4    D5    D6 
# 10668  8643  4464  4013  5480  7302

#original term
6508 4443


########filter by snapatac2 barcode

#barcodes <- readRDS('barcodes.rds')

#barcodes_peakratio <- readRDS('barcode.filtered.peak_ratio_logUMI3.5.rds')

# barcodes_snapatac2 = read.table('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_late_combine/02.snapatac2_harmony/cluster_df_add.txt',sep='\t',header=TRUE,comment.char = "",stringsAsFactors = FALSE,row.names = 1)
# 42579 × 13

# table(rownames(barcodes_snapatac2) %in% barcodes$barcode.add )
# FALSE  TRUE 
# 18582 23997

# TRUE 
# 42579

# # FALSE  TRUE 
# #  9071 33508

# table(barcodes$barcode.add %in% rownames(barcodes_snapatac2)  )
# FALSE  TRUE 
#  7687 23997

#   FALSE    TRUE 
# 2733518   42579 

# # FALSE  TRUE 
# #  7062 33508 


# #barcode_share <- intersect(rownames(barcodes_snapatac2),barcodes_peakratio$barcode.add)
# #33508

# idx <- barcodes$barcode.add %in% rownames(barcodes_snapatac2)

# #idx <- barcodes$barcode.add %in% barcode_share

# table(barcodes[idx,'sample'])

#   D1   D2   D3   D4   D5   D6 
# 8742 7321 5754 5217 7417 8128


# ##data.obs['sample'].value_counts().sort_index()
# PLA-atac-late1    8742
# PLA-atac-late2    7321
# PLA-atac-late3    5754
# PLA-atac-late4    5217
# PLA-atac-late5    7417
# PLA-atac-late6    8128

# #   D1   D2   D3   D4   D5   D6 
# # 7817 6418 4037 3775 5082 6379 

# #  D1   D2   D3   D4   D5   D6 
# # 8742 7321 5754 5217 7417 8128

# #  D1   D2   D3   D4   D5   D6 
# # 7817 6418 4037 3775 5082 6379

# ggplot(barcodes, aes(x=logUMI, y=peak_ratio)) + 
#     geom_point(size=0.1, col="grey") +
#     geom_point(data=barcodes[idx,], aes(x=logUMI, y=peak_ratio),size=0.3,col='red'  ) +
#     theme_classic()	+
#     #ggtitle(paste0(sample.name,' peak_ratio plot, cutoff.logUMI: ',cutoff.logUMI.low,'-',cutoff.logUMI.high, '\ncutoff.FRIP: ',cutoff.FRIP.low, '-',cutoff.FRIP.high, ' total selected:',sum(idx))) +
#     ggtitle(paste0(sample.name,' use snapatac2 barcode')) +
#     ylim(0, 1) + xlim(0, 6) + 
#     labs(x = "log10(UMI)", y="peak ratio")



# barcodes.bk <- barcodes
# barcodes = barcodes[idx,]






######

saveRDS(object = barcodes,'barcode.filtered.promoter_ratio_0.2-0.5.logUMI4-5.rds')
barcodes <- readRDS('barcode.filtered.promoter_ratio_0.2-0.5.logUMI4-5.rds')

#saveRDS(object = barcodes,'barcode.filtered.peak_ratio_logUMI3.5.rds')
#saveRDS(object = barcodes,'barcode.filtered.logUMI4.rds')

#saveRDS(object = barcodes,'barcode.filtered.peak_ratio_logUMI3.5.overlap.snapatac2.rds')

#saveRDS(barcodes, 'barcode.overlap_snapatac2.rds')


barcode.list <- split(barcodes,f=barcodes$sample)
saveRDS(barcode.list,'barcode.list.rds')

barcode.list <- readRDS('barcode.list.rds')


all.equal(barcodes ,Reduce(rbind,barcode.list)) #TRUE



# ##############step 2 select barcode (with calculated snapATAC result)#######
# ##read in filtered cell barcode from both snapATAC dir

# # barcode.sel1 <- read.table("../../placenta_10X/02.snapATAC/PLA-term-ATAC-1/EVT_rescue/barcode.selected",stringsAsFactors = FALSE)$V1
# # #6539

# # barcode.sel2 <- read.table("../../placenta_10X_new/03.snapATAC/barcode.selected",stringsAsFactors = FALSE)$V1
# # #8107

# barcode.file.list =c("../../placenta_10X/02.snapATAC/PLA-term-ATAC-1/EVT_rescue/barcode.selected","../../placenta_10X_new/03.snapATAC/barcode.selected")


# barcode.list = lapply(barcode.file.list, function(file){
#     read.table(file)[,1];
#   })

# names(barcode.list) <- sample.names

names(x.sp.ls)
'placenta_donor1''placenta_donor2''placenta_donor3''placenta_donor5''placenta_donor6''placenta_donor9'
#'placenta_donor1''placenta_donor2''placenta_donor3''placenta_donor4''placenta_donor5''placenta_donor6'

#'placenta_donor1' 'placenta_donor2'
names(barcode.list)
'D1''D2''D3''D5''D6''D9'

#'D1''D2''D3''D4''D5''D6'

#'D1' 'D2'

lapply(barcode.list,function(x){nrow(x)})
$D1 6508 $D2 4443 $D3 5079 $D5 4105$ D6 5785 $D9 5764

$D1 8742
$D2 7321
$D3 5754
$D4 5217
$D5 7417
$D6 8128



x.sp.list = lapply(seq(x.sp.ls), function(i){
    x.sp = x.sp.ls[[i]];
    x.sp  = x.sp[x.sp@barcode %in% barcode.list[[i]][,'barcode'],];
  })

names(x.sp.list) = sample.names

$placenta_donor1
number of barcodes: 6508
number of bins: 0
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor2
number of barcodes: 4443
number of bins: 0
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor3
number of barcodes: 5079
number of bins: 0
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor5
number of barcodes: 4105
number of bins: 0
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor6
number of barcodes: 5785
number of bins: 0
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor9
number of barcodes: 5764
number of bins: 0
number of genes: 0
number of peaks: 0
number of motifs: 0

# $placenta_donor1
# number of barcodes: 8742
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 7321
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor3
# number of barcodes: 5754
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor4
# number of barcodes: 5217
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor5
# number of barcodes: 7417
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor6
# number of barcodes: 8128
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


# $placenta_donor1
# number of barcodes: 7817
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 6418
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor3
# number of barcodes: 4037
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor4
# number of barcodes: 3775
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor5
# number of barcodes: 5082
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor6
# number of barcodes: 6379
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# x.sp.list
# $placenta_donor1
# number of barcodes: 6508
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 4443
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor1
# number of barcodes: 11377
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 9561
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor1
# number of barcodes: 6249
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 8107
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


# placenta_donor1
# number of barcodes: 6539
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# placenta_donor2
# number of barcodes: 8107
# number of bins: 0
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


##



######step 3 add cell-by-bin matrix
x.sp.list = lapply(seq(x.sp.list), function(i){
    x.sp = addBmatToSnap(x.sp.list[[i]], bin.size=5000);
    x.sp
  })


names(x.sp.list) = sample.names

x.sp.list

$placenta_donor1
number of barcodes: 6508
number of bins: 620094
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor2
number of barcodes: 4443
number of bins: 620094
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor3
number of barcodes: 5079
number of bins: 620059
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor5
number of barcodes: 4105
number of bins: 620059
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor6
number of barcodes: 5785
number of bins: 620059
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor9
number of barcodes: 5764
number of bins: 620059
number of genes: 0
number of peaks: 0
number of motifs: 0

# $placenta_donor1
# number of barcodes: 8742
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 7321
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor3
# number of barcodes: 5754
# number of bins: 620059
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor4
# number of barcodes: 5217
# number of bins: 620059
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor5
# number of barcodes: 7417
# number of bins: 620059
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor6
# number of barcodes: 8128
# number of bins: 620059
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor1
# number of barcodes: 7817
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 6418
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor3
# number of barcodes: 4037
# number of bins: 620059
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor4
# number of barcodes: 3775
# number of bins: 620059
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor5
# number of barcodes: 5082
# number of bins: 620059
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor6
# number of barcodes: 6379
# number of bins: 620059
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor1
# number of barcodes: 6508
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 4443
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor1
# number of barcodes: 11377
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 9561
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


# $placenta_donor1
# number of barcodes: 6249
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 8107
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


# placenta_donor1
# number of barcodes: 6539
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# placenta_donor2
# number of barcodes: 8107
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0



#######step 4 combine two snap object############

#all.equal (x.sp.list[[1]]@feature$name,x.sp.list[[2]]@feature$name ) #TRUE


#no need to intersect bin

##if not, To combine these two snap objects, common bins are selected
# bin.shared = Reduce(intersect, lapply(x.sp.list, function(x.sp) {x.sp@feature$name}  ))
# x.sp.list <- lapply(x.sp.list, function(x.sp){
#     idy = match(bin.shared, x.sp@feature$name);
#     x.sp[,idy, mat="bmat"];
#   })


all.equal (x.sp.list[[1]]@feature$name,x.sp.list[[2]]@feature$name ) #TRUE
all.equal (x.sp.list[[3]]@feature$name,x.sp.list[[4]]@feature$name ) #TRUE
all.equal (x.sp.list[[5]]@feature$name,x.sp.list[[6]]@feature$name ) #TRUE
all.equal (x.sp.list[[3]]@feature$name,x.sp.list[[5]]@feature$name ) #TRUE


all.equal (x.sp.list[[1]]@feature$name,x.sp.list[[3]]@feature$name ) #FALSE different at patch chromosomes

bin.shared.1and3 = intersect(x.sp.list[[1]]@feature$name,x.sp.list[[3]]@feature$name)

bin.diff1 <- x.sp.list[[1]]@feature$name[!(x.sp.list[[1]]@feature$name %in% bin.shared.1and3)]
unique(sort(sapply( X =  stringr::str_split(string = bin.diff1, pattern = ":", n=2), function(x){ x[1]  } )) )
'chr1_KI270706v1_random''chr1_KI270707v1_random''chr1_KI270708v1_random''chr1_KI270709v1_random''chr1_KI270710v1_random''chr1_KI270711v1_random''chr1_KI270712v1_random''chr1_KI270713v1_random''chr1_KI270714v1_random''chr11_KI270721v1_random''chr14_GL000009v2_random''chr14_GL000194v1_random' ... 'chrUn_KI270748v1''chrUn_KI270749v1''chrUn_KI270750v1''chrUn_KI270751v1''chrUn_KI270752v1''chrUn_KI270753v1''chrUn_KI270754v1''chrUn_KI270755v1''chrUn_KI270756v1''chrUn_KI270757v1''chrY_KI270740v1_random'


bin.diff3 <- x.sp.list[[3]]@feature$name[!(x.sp.list[[3]]@feature$name %in% bin.shared.1and3)]
unique(sort(sapply( X =  stringr::str_split(string = bin.diff3, pattern = ":", n=2), function(x){ x[1]  } )) )
'GL000008.2''GL000009.2''GL000194.1''GL000195.1''GL000205.2''GL000208.1''GL000213.1''GL000214.1''GL000216.2''GL000218.1''GL000219.1''GL000220.1''GL000221.1''GL000224.1''GL000225.1''GL000226.1''KI270302.1''KI270303.1''KI270304.1''KI270305.1''KI270310.1''KI270311.1''KI270312.1''KI270315.1''KI270316.1''KI270317.1''KI270320.1''KI270322.1''KI270329.1''KI270330.1''KI270333.1''KI270334.1' ... 'KI270740.1''KI270741.1''KI270742.1''KI270743.1''KI270744.1''KI270745.1''KI270746.1''KI270747.1''KI270748.1''KI270749.1''KI270750.1''KI270751.1''KI270752.1''KI270753.1''KI270754.1''KI270755.1''KI270756.1''KI270757.1'


lapply(seq(x.sp.list), function(x){ length(x.sp.list[[x]]@feature$name)   }   )
620094
620094
620059
620059
620059
620059

# 620094
# 620094
# 620059
# 620059
# 620059
# 620059


#need to intersect bin

#if not, To combine these two snap objects, common bins are selected
bin.shared = Reduce(intersect, lapply(x.sp.list, function(x.sp) {x.sp@feature$name}  )) #617669 shared

all.equal(bin.shared.1and3,bin.shared) #TRUE

all.equal(x.sp.list[[1]]@feature$name[match(bin.shared, x.sp.list[[1]]@feature$name)],bin.shared )#TRUE

x.sp.list <- lapply(x.sp.list, function(x.sp){
    idy = match(bin.shared, x.sp@feature$name);
    x.sp[,idy, mat="bmat"];
  })


x.sp.list

$placenta_donor1
number of barcodes: 6508
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor2
number of barcodes: 4443
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor3
number of barcodes: 5079
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor5
number of barcodes: 4105
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor6
number of barcodes: 5785
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

$placenta_donor9
number of barcodes: 5764
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

# $placenta_donor1
# number of barcodes: 8742
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 7321
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor3
# number of barcodes: 5754
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor4
# number of barcodes: 5217
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor5
# number of barcodes: 7417
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor6
# number of barcodes: 8128
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


# $placenta_donor1
# number of barcodes: 7817
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor2
# number of barcodes: 6418
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor3
# number of barcodes: 4037
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor4
# number of barcodes: 3775
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor5
# number of barcodes: 5082
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# $placenta_donor6
# number of barcodes: 6379
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0



all.equal (x.sp.list[[1]]@feature$name,x.sp.list[[3]]@feature$name ) #TRUE


##merge bin mat with common bin id

x.sp = Reduce(snapRbind, x.sp.list)

x.sp

number of barcodes: 31684
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

number of barcodes: 42579
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

# number of barcodes: 33508
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 10951
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 20938 (11377 + 9561)
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 14356 (6249 + 8107)
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 14646  (6539 + 8107)
# number of bins: 620094
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


table(x.sp@sample)

placenta_donor1 placenta_donor2 placenta_donor3 placenta_donor5 placenta_donor6 
           6508            4443            5079            4105            5785 
placenta_donor9 
           5764 

# placenta_donor1 placenta_donor2 placenta_donor3 placenta_donor4 placenta_donor5 
#            8742            7321            5754            5217            7417 
# placenta_donor6 
#            8128


# placenta_donor1 placenta_donor2 placenta_donor3 placenta_donor4 placenta_donor5 
#            7817            6418            4037            3775            5082 
# placenta_donor6 
#            6379 

# placenta_donor1 placenta_donor2 
#            6508            4443 

# placenta_donor1 placenta_donor2 
#           11377            9561

# placenta_donor1 placenta_donor2 
#            6249            8107 


#placenta_donor1 placenta_donor2 
#           6539            8107 

rm(x.sp.list)
rm(x.sp.ls)
gc()


saveRDS(x.sp,'x.sp.stage0.rds') #one can restart from here ?
#x.sp = readRDS('x.sp.stage0.rds')



####################step 5 matrix binarization########
x.sp = makeBinary(x.sp, mat="bmat") #bmat 20938 x 620094 #bmat: 14646 x 620094

number of barcodes: 31684
number of bins: 617669
number of genes: 0
number of peaks: 0
number of motifs: 0

# number of barcodes: 42579
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 33508
# number of bins: 617669
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


# #fill the feature slot (instead of bin, use peaks)
# feature_file = "../../01.data_cellranger_atac/PLA-term-ATAC-1/filtered_peak_bc_matrix/peaks.bed"
# features <- readr::read_tsv(feature_file, col_names = F)  
# features = as.data.frame(features)
# features.ranges = GRanges(
#     seqnames = features[,1], 
#     ranges = IRanges(features[,2], features[,3]),
	
#   );

# x.sp@feature = features.ranges



#######read in late1 late2 bin feature and use (and skip following steps), not use###

feature.gr <- readRDS('../02.snapATAC_harmony_tuning/02.snapATAC_harmony_late1_late2/feature.gr.rds')


idy = queryHits(
    findOverlaps(x.sp@feature, feature.gr)
  );


if(length(idy) > 0){
    x.sp = x.sp[,idy, mat="bmat"];
  };

x.sp

number of barcodes: 31684
number of bins: 538607
number of genes: 0
number of peaks: 0
number of motifs: 0


# number of barcodes: 16030
# number of bins: 539783
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 16715
# number of bins: 539783
# number of genes: 0
# number of peaks: 0
# number of motifs: 0





###################step 6 bin filtering###############
library(GenomicRanges)
black_list = read.table("blacklist.bed")
black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );

idy = queryHits(
    findOverlaps(x.sp@feature, black_list.gr)
  );

##remove bins hitting blacklist
if(length(idy) > 0){
    x.sp = x.sp[,-idy, mat="bmat"];
  };

x.sp

number of barcodes: 31684
number of bins: 617643
number of genes: 0
number of peaks: 0
number of motifs: 0

number of barcodes: 42579
number of bins: 617643
number of genes: 0
number of peaks: 0
number of motifs: 0

# number of barcodes: 33508
# number of bins: 617643
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 10951
# number of bins: 620068
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 20938
# number of bins: 620068
# number of genes: 0
# number of peaks: 0
# number of motifs: 0


#number of barcodes: 14356
# number of bins: 620068
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 14646
# number of bins: 620068
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

##remove unwanted chromosomes
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM|chrUn|chrEBV", seqlevels(x.sp@feature))];

#chr.exclude.1  = grep("random|chrM", seqlevels(x.sp@feature),value = TRUE)
#all.equal(chr.exclude.1,chr.exclude) #TRUE

idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){
    x.sp = x.sp[,-idy, mat="bmat"]
  };

x.sp

number of barcodes: 31684
number of bins: 617639
number of genes: 0
number of peaks: 0
number of motifs: 0

number of barcodes: 42579
number of bins: 617639
number of genes: 0
number of peaks: 0
number of motifs: 0

# number of barcodes: 33508
# number of bins: 617639
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 10951
# number of bins: 617639
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 20938
# number of bins: 617639
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 14356
# number of bins: 617639
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 14646
# number of bins: 617639
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

##remove the top 5% bins that overlap with invariant features such as the house keeping gene promoters.
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(bin Cov)", 
    col="lightblue", 
    xlim=c(0, 5)
  );


bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95) #remove 0 then quantile
95%: 3.62304243424638

95%: 3.6612

#95%: 3.58927922123597

#95%: 3.20112389720738
#95%: 3.32674537956532
##95%: 3.37
##95%:  3.383

#95%: 3.18525876529659
###95%: 3.43
#########95%: 2.97726621242729 in example
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);

545479

546157

#544977

x.sp = x.sp[, idy, mat="bmat"];

x.sp

number of barcodes: 31684
number of bins: 545479
number of genes: 0
number of peaks: 0
number of motifs: 0

number of barcodes: 42579
number of bins: 546157
number of genes: 0
number of peaks: 0
number of motifs: 0

# number of barcodes: 33508
# number of bins: 544977
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 10951
# number of bins: 539783
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 20938
# number of bins: 540667
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

#number of barcodes: 14356
# number of bins: 544717
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 14646
# number of bins: 545204
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

## remove any cells of bin coverage less than 1,000,The rational behind this is that some cells may have high number of unique fragments but end up with low bin coverage after filtering. This step is optional but highly

idx = which(Matrix::rowSums(x.sp@bmat) > 1000); #31684 #33508 barcode,no need to filter #10951 #20935 #14356 #14646, no need to filter
#x.sp = x.sp[idx,];


# x.sp 
# number of barcodes: 20935
# number of bins: 540667
# number of genes: 0
# number of peaks: 0
# number of motifs: 0

#saveRDS(x.sp,"x.sp.stage1.rds")
#x.sp <- readRDS("x.sp.stage1.rds")


############step 7 dimentionality reduction (diffusion map with jaccard similarity index)#############

# ########for large cell population####
# #diffusion maps algorithm, a nonlinear dimensionality reduction technique that discovers low-dimension manifold by performing random walk on the data and is highly robust to noise and perturbation.
# row.covs.dens <- density(
#     x = x.sp@metaData[,"logUMI"], 
#     bw = 'nrd', adjust = 1
#   );

# #density-based sampling approach
# sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps);
# set.seed(1);
# #idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = 6820, prob = sampling_prob)); #no sampling
# idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = 3000, prob = sampling_prob)); #half sampling
# #idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = 1000, prob = sampling_prob)); #less sampling


# ##Split the x.sp into landmark (x.landmark.sp) and query (x.query.sp) cells.
# x.landmark.sp = x.sp[idx.landmark.ds,]; #3000 x 544771
# x.query.sp = x.sp[-idx.landmark.ds,]; #3820 x 544771


# #Run diffusion maps on the landmark cells.

# x.landmark.sp = runDiffusionMaps(
#     obj= x.landmark.sp,
#     input.mat="bmat", 
#     num.eigs=50
#   );

# x.landmark.sp@metaData$landmark = 1;

# #Project query cells to landmarks

# x.query.sp = runDiffusionMapsExtension(
#     obj1=x.landmark.sp, 
#     obj2=x.query.sp,
#     input.mat="bmat"
#   );
# x.query.sp@metaData$landmark = 0;

# ##Combine landmark and query cells.
# x.sp = snapRbind(x.landmark.sp, x.query.sp);
# x.sp = x.sp[order(x.sp@metaData[,"sample"])]; #IMPORTANT

# x.sp
# #number of barcodes: 6820
# #number of bins: 544771
# #number of genes: 0
# #number of peaks: 0
# #number of motifs: 0


#####simply run dimensionality reduction with diffusion map (small cell population)
# x.sp = runDiffusionMaps( #ok for 14356 cells #ok for 14646 cells
#     obj=x.sp,
#     input.mat="bmat", 
#     num.eigs=50
#   );

#####jmat -> nmat #normalized jaccard similarity index



##run without normalization, 
# source("runDiffusionMaps_new.r")
# x.sp = runDiffusionMaps_new( #without normalization for jmat
#     obj=x.sp,
#     input.mat="bmat", 
#     num.eigs=50
#   );
#all.equal(x.sp@jmat@jmat,x.sp@jmat@nmat)#TRUE


##saveRDS(x.sp,"x.sp.onestep_filter.rds")
#x.sp.1 = readRDS("x.sp.stage2.use10xmat.rds")
#all.equal(x.sp,x.sp.1)

###x.sp <- readRDS("x.sp.stage2.rds")

##x.sp = x.landmark.sp

###########step 8 determine significant components#####################
# #simply looking at a pairwise plot and select the number of eigen vectors that the scatter plot starts looking like a blob
# plotDimReductPW(  #for smat, cell x 50
#     obj=x.sp, 
#     eigs.dims=1:50,
#     point.size=0.3,
#     point.color="grey",
#     point.shape=19,
#     point.alpha=0.6,
#     down.sample=5000,
#     pdf.file.name=NULL, 
#     pdf.height=7, 
#     pdf.width=7
#   );
##choose k = 29 #half sampling
##choose k = 15 #1000 sampling
##choose k = 12 #no sampling ?

##use eigs = 16

# ######step 9 remove batch effect (will read and rewrite x.sp@smat@dmat etc.)####
# library(harmony)
# x.after.sp = runHarmony(
#     obj=x.sp, 
#     eigs.dim=1:16, 
#     #eigs.dim=1:50, 
#     meta_data=x.sp@sample # sample slot
#   );

# #converged after 10 iterations

#############step 10 graph-based clustering###########
# x.sp = runKNN(
#     obj=x.sp,
#     eigs.dims=1:16,
#     #k=29 
# 	k=15
# 	#k=12
#     #k=30 #for not norm smat
#   );

# library(leiden); #an extension of louvain

# x.sp=runCluster( #from knn graph, cell x cell
#     obj=x.sp,
#     tmp.folder=tempdir(),
#     louvain.lib="R-igraph", #"R-igraph" "leiden"
#     seed.use=10,
#     resolution=1 #0.7, 1.5
#   );


# ###get neightbor graph
# x.after.sp = runKNN(
#     obj=x.after.sp,
#     eigs.dims=1:16,
#     #k=29 
# 	k=15
# 	#k=12
#     #k=30 #for not norm smat
#   );


# ###clustering 
# library(leiden); #an extension of louvain
# x.after.sp=runCluster( #from knn graph, cell x cell
#     obj=x.after.sp,
#     tmp.folder=tempdir(),
#     louvain.lib="R-igraph", #"R-igraph" "leiden"
#     seed.use=10,
#     resolution=1 #0.7, 1.5
#   );



################################step 11 visulization######################



#library('umap');
##umap
# x.sp = runViz(
#     obj=x.sp, 
#     tmp.folder=tempdir(),
#     dims=2,
#     #eigs.dims=1:20, 
# 	#eigs.dims=2:30,
# 	eigs.dims=1:16,
#     method="umap",
#     seed.use=10
#   );

# x.after.sp = runViz(
#     obj=x.after.sp, 
#     tmp.folder=tempdir(),
#     dims=2,
#     #eigs.dims=1:20, 
# 	#eigs.dims=2:30,
# 	eigs.dims=1:16,
#     method="umap",
#     seed.use=10 #
#   );




##save bin feature

saveRDS(x.after.sp@feature,'feature.gr.rds') #545479

#feature <- readRDS('feature.gr.rds') #d1 d2 bin feature


#feature <- readRDS('feature.rds')
#538747

feature.df <- as.data.frame(x.after.sp@feature)

options(scipen=999)
feature.df$start <- feature.df$start - 1
feature.df$name <- paste0( feature.df$seqnames,":",feature.df$start,'-',  feature.df$end)


write.table(x=feature.df[,c(1,2,3,6)],file = "feature.gr.df.bed",col.names = FALSE,row.names = FALSE, sep = '\t',quote = FALSE)






####step 7 dimentionality reduction (diffusion map with jaccard similarity index)#
###simply run dimensionality reduction with diffusion map (small cell population)
x.sp = runDiffusionMaps( #slow but ok for 20935 cells #ok for 14356 cells #ok for 14646 cells
    obj=x.sp,
    input.mat="bmat", 
    num.eigs=50
  );

##run without normalization, 
# source("runDiffusionMaps_new.r")
# x.sp = runDiffusionMaps_new( #without normalization for jmat
#     obj=x.sp,
#     input.mat="bmat", 
#     num.eigs=50
#   );
#all.equal(x.sp@jmat@jmat,x.sp@jmat@nmat)#TRUE



###step 8 determine significant components#######
#simply looking at a pairwise plot and select the number of eigen vectors that the scatter plot starts looking like a blob
res.p <- plotDimReductPW(  #for smat, cell x 50
    obj=x.sp, 
    eigs.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
  );
##choose k = 29 #half sampling
##choose k = 15 #1000 sampling
##choose k = 12 #no sampling ?

##use eigs = 16
use eigs = 25

##save before harmony 
saveRDS(x.sp,"x.sp.before_harmony.rds")

x.sp <- readRDS("x.sp.before_harmony.rds")

######step 9 remove batch effect (will read and rewrite x.sp@smat@dmat etc.)####
library(harmony)

table(x.sp@sample)

placenta_donor1 placenta_donor2 placenta_donor3 placenta_donor5 placenta_donor6 
           6508            4443            5079            4105            5785 
placenta_donor9 
           5764 

placenta_donor1 placenta_donor2 placenta_donor3 placenta_donor4 placenta_donor5 
           8742            7321            5754            5217            7417 
placenta_donor6 
           8128 


# placenta_donor1 placenta_donor2 placenta_donor3 placenta_donor4 placenta_donor5 
#            7817            6418            4037            3775            5082 
# placenta_donor6 
#            6379 

source('runHarmony.r')

x.after.sp = runHarmony(
    obj=x.sp, 
    eigs.dim=1:25,
    #eigs.dim=1:16, 
    #eigs.dim=1:50, 
    meta_data=x.sp@sample, # sample slot
    max.iter.harmony = 20
  );

Harmony converged after 7 iterations

Harmony converged after 6 iterations

#Harmony converged after 3 iterations
#Harmony converged after 9 iterations

#converged after 10 iterations



#######step 10 graph-based clustering and umap#######
###get neightbor graph from x.after.sp@smat@dmat(after harmony) cell x cell
x.after.sp = runKNN(
  obj=x.after.sp,
  #eigs.dims=1:16,
  eigs.dims=1:25,
  #k=29 
  #k=15 #a higher k corresponding to a lower resolution
  #k=12
  k=30
);


##get graph-based clustering from knn graph, 
x.after.sp=runCluster( 
    obj=x.after.sp,
    tmp.folder=tempdir(),
    louvain.lib='leiden', #"R-igraph" "leiden"  #leiden is a py package, which also used by signac
    seed.use=10, #keep static? seed = 10
    resolution=0.9 #0.7,1, 1.5
  );


#source('runViz_new.r')
#x.after.sp = runViz_new( 
x.after.sp = runViz( 
    obj=x.after.sp, 
    tmp.folder=tempdir(),
    dims=2,
    #eigs.dims=1:20, 
    #eigs.dims=2:30,
    eigs.dims=1:25,
    #eigs.dims=1:16,
    method="umap",#default use umap::umap, modified to use  uwot::umap
    seed.use=123 #
  );

source('plotViz.new.r')

options(repr.plot.height=7.5,repr.plot.width=7.5)
plotViz.new( #with text label halo, point size = 1
#plotViz( #with text label halo, point size = 1
    obj=x.after.sp,
    method="umap", 
    main=paste("snapATAC harmony dims=1:25 K=30"," UMAP\n seed:",123,' res:',0.9,' repeat:',1,' cluster umap:umap',sep=''),
    #main=paste("snapATAC harmony dims=1:16 K=30"," UMAP\n seed:",123,' res:',0.9,' repeat:',1,' cluster umap:umap',sep=''),
    #main=paste("snapATAC harmony dims=1:16 K=",k," UMAP\n seed:",seed,' res:',res,' cluster:',cluster_method,sep=''),
    point.color=x.after.sp@cluster, 
    point.size=0.3, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    #down.sample=10000,
    legend.add=FALSE,
    xlim = c(-10,10),
    ylim = c(-10,15)
    );




###restart tuning from here
#saveRDS(x.after.sp,"x.after.sp.rds")
x.after.sp <- readRDS('x.after.sp.rds')



##



#################step 11  umap embedding with tunning (better before cluster decision)###################



#seed = 123 ##choose one from the iterative umap tunning
#can also tuning umap parameters, m_dist spread, but this seems good at this moment


##first quick look at the original umap
source("runViz_new.r")

##get umap from @smat@dmat with uwot::umap or umap::umap 
# x.after.sp = runViz_new( 
#     obj=x.after.sp, 
#     tmp.folder=tempdir(),
#     dims=2,
#     #eigs.dims=1:20, 
#     #eigs.dims=2:30,
#     eigs.dims=1:25,
#     #eigs.dims=1:16,
#     method="umap",#default use umap::umap, modified to use  uwot::umap
#     seed.use=seed #
#   );



#plot original umap with cluster
options(repr.plot.height=7.5,repr.plot.width=7.5)
p=plotViz( #with text label halo, point size = 1
    obj=x.after.sp,
    method="umap", 
    main=paste("snapATAC harmony dims=1:25 K=",30," UMAP\n seed:",seed,' res:',res,' cluster:',cluster_method,sep=''),
    #main=paste("snapATAC harmony dims=1:16 K=",30," UMAP\n seed:",seed,' res:',res,' cluster:',cluster_method,sep=''),
    point.color=x.after.sp@cluster, 
    point.size=0.3, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    #down.sample=10000,
    legend.add=FALSE
    );



##directly tuning umap by uwot::umap
spread = 1
min_dist = 0.01
seed = 123
eigs.dims = 1:25

data.use = x.after.sp@smat@dmat[, eigs.dims]
#31684 × 25

rownames(data.use) <- rownames(x.after.sp@metaData)
colnames(data.use) <- paste0('X',1:ncol(data.use))

saveRDS(data.use,'smat_dmat.rds')


cluster <- x.after.sp@cluster

res.p.umap <- list()

#par(mfrow=c(2,2))
for (min_dist in c(0.1,0.2,0.3,0.4,0.5)){
#for (mdist in c(0.3)){
    cat('min_dist  ',min_dist,'\n')
    for (spread in c(0.5,1.0,1.5,2.0)){
#    for (spread in c(1.0)){    
        cat('  spread  ',spread,'\n')
        umap.df = uwot::umap(data.use,min_dist = min_dist,spread = spread )
        options(repr.plot.height = 7.5, repr.plot.width = 7.5)
        #plot(umap.df, cex=0.2, pch = 19, main = paste( 'seed: ',seed, 'eigs.dims: ', eigs.dims, 'min_dist: ', min_dist,' spread: ',spread,sep=''  )  )
        
        colnames(umap.df) <- c('UMAP_1','UMAP_2')
        cluster.df  <- data.frame(cluster=cluster, umap.df)
        res.p <- ggplot(cluster.df,aes (x=UMAP_1, y=UMAP_2, col = cluster) ) +
          geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
          #scale_colour_manual(values = color_good)  +
          scale_colour_manual(values = color_snap_mod1)  +

          theme(
                legend.position = 'right',
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
        ggtitle(paste0( 'min_dist:',min_dist,'-spread:',spread )  ) +
        guides(col = guide_legend(override.aes = list(size = 3))) +  ##no effect ??
       #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
       #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
       labs(x = "UMAP1", y = "UMAP2")+
       geom_point(size=0.1) +
       theme()
        
       res.p.umap[[paste0( 'min_dist:',min_dist,'-spread:',spread ) ]] <- res.p
        
        
    }
    
    
}

names(res.p.umap)
for(i in names(res.p.umap) ){ print(res.p.umap[[i]])   }


saveRDS(res.p.umap,'res.p.umap.rds')


###====try to fix parameters for umap embedding========##

min_dist = 0.4
spread = 1
res.p <- res.p.umap[[paste0( 'min_dist:',min_dist,'-spread:',spread ) ]]
plot(res.p$data[,2:3],pch=19,cex=0.1)


############rotate the umap xy coordinate if needed#############
#test <- data.frame(x=c(10,2),y=c(25,45))
#plot(test,xlim = c(-50,50),ylim = c(-50,50),cex=2 );abline(h = 0);abline(v=0)
#test.rotate90 <- data.frame(x=test$x*cos(90)-test$y*sin(90),y=test$y*cos(90)+test$x*sin(90))
#points(test.rotate90,cex=2,pch=19,col='red')

umap.df <- res.p$data[,2:3]

#umap.df <- data.frame(UMAP_1=x.after.sp@umap[,1],UMAP_2=x.after.sp@umap[,2])

rownames(umap.df) <- rownames(x.after.sp@metaData)

##iterative rotate the umap direction
#for(degree in seq(from = 50,100,10) ) { 
for(degree in c(50) ) { 
    rotate_d <- degree*(pi/180) #degree to radian
    cat('degree,rotate_d:',degree,' ',rotate_d)
    plot(umap.df,cex=0.2,col='grey',main= paste('degree:',degree,sep=''))#,xlim = c(-8,5),ylim = c(-10,10));
    abline(h = 0);
    abline(v=0) 
    UMAP.rotate <- data.frame(
        UMAP_1=umap.df$UMAP_1*cos(rotate_d)-umap.df$UMAP_2*sin(rotate_d),
        UMAP_2=umap.df$UMAP_1*sin(rotate_d)+umap.df$UMAP_2*cos(rotate_d)
    )
    #UMAP.rotate <- data.frame(UMAP_1=-1*umap.df$UMAP_2,UMAP_2=umap.df$UMAP_1)
    points(UMAP.rotate,pch=19,col='red',cex=0.2)
}

fix degree = 50

##flip y
UMAP.rotate[,1] <- -1 * UMAP.rotate[,1]

#rewrite the umap slot
#rownames(UMAP.rotate) <- rownames(x.after.sp@metaData)
UMAP.bk = x.after.sp@umap 
x.after.sp@umap = as.matrix(UMAP.rotate)

#saveRDS(object=UMAP.rotate,file = 'uwot.umap.seed177.rotate150.rds')


##look at the tuned, rotated umap
options(repr.plot.height=7.5,repr.plot.width=7.5)
plotViz.new( #with text label halo, point size = 1
#plotViz( #with text label halo, point size = 1
    obj=x.after.sp,
    method="umap", 
    main=paste("snapATAC harmony dims=1:25 K=30"," uwot::umap\n seed:",123,' res:',0.9,' min_dist:',min_dist, ' spread: ',spread,' rotate:',degree,sep=''),
    point.color=x.after.sp@cluster, 
    point.size=.5, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=TRUE,
    text.size=2,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    #down.sample=10000,
    legend.add=FALSE,
    xlim = c(-10,15),
    ylim = c(-10,8)
    );



#################





#saveRDS(object = x.after.sp@umap,'uwot.umap.tuned.rds')

#####umap evaluation??######
##unlike tunning the tsne embedding with  Kullback-Leibler divergence reported by tsne with certain Plexity, tuning of the umap embedding, can alternate the final umap embedding by tunning up/down spread and min_dist
##https://jlmelville.github.io/uwot/abparams.html

##not doing this tuning here, use the default spread and min_dist in uwot::umap 



############################



###########do some iterative tuning for clustering and go back to fix the clustering parameters###
##compare cluster assignment by plot as http://bioconductor.org/books/release/OSCA/clustering.html#cluster-separation-redux
##better do umap iterative first

#k = 15 #30  #a higher k corresponding to a lower resolution
#res = 0.9 #0.8
cluster_method = 'leiden' #"R-igraph" "leiden"
k = 30
seed = 123

source('plotViz.new.r')

res.clusters = list()
rep <- 0
# for(k in c(10,15,20,25,30)){ #use k = 30
#     cat('k: ',k,'\n')
# for(cluster_method in c("R-igraph", "leiden")){ #use leiden
#     cat('cluster_method: ',cluster_method,'\n')   
for(r in seq(0.1,1,0.1) ){ #use 
#for(r in c(0.8,0.82,0.85,0.88,0.9,0.9,0.92,0.95) ){ #use 0.9?
#for(r in c(0.82,0.82,0.82,0.82,0.82,0.82,0.82,0.82) ){ #use 0.82, and repeatively until the final clustering
    cat('res: ',r,'\n') 
    rep <- rep + 1
        
#     tic('runKNN')
#     ###get neightbor graph from x.after.sp@smat@dmat(after harmony) cell x cell
#     x.after.sp = runKNN(
#       obj=x.after.sp,
#       eigs.dims=1:16,
#       k=k
#     );
#     toc()
    
    tic('runCluster')
    ##get graph-based clustering from knn graph, 
    x.after.sp=runCluster( 
        obj=x.after.sp,
        tmp.folder=tempdir(),
        #louvain.lib='leiden',
        louvain.lib=cluster_method,
        seed.use=seed, #keep static? seed = 10
        #resolution=res
        resolution=r
      );
    toc()
    
    #plot umap with cluster
    options(repr.plot.height=7.5,repr.plot.width=7.5)
    
    p=plotViz.new( #with text label halo, point size = 1
        obj=x.after.sp,
        method="umap", 
       # main=paste("snapATAC harmony dims=1:16 K=",k," UMAP\n seed:",seed,' res:',r,' repeat:',rep,' cluster:',cluster_method,sep=''),
        main=paste("snapATAC harmony dims=1:25 K=",k," UMAP\n seed:",seed,' res:',r,' repeat:',rep,' cluster:',cluster_method,sep=''),
        point.color=x.after.sp@cluster, 
        point.size=0.3, 
        point.shape=19, 
        point.alpha=0.8, 
        text.add=TRUE,
        text.size=1,
        text.color="black",
        text.halo.add=TRUE,
        text.halo.color="white",
        text.halo.width=0.2,
        #down.sample=10000,
        legend.add=FALSE,
        xlim = c(-8,10),
        ylim = c(-8,10)
        );
    
    
    #res.clusters[[paste0('k',k)]] = x.after.sp@cluster
    #res.clusters[[cluster_method]] = x.after.sp@cluster
    #res.clusters[[paste0('r',r)]] = x.after.sp@cluster
    res.clusters[[paste0('r',rep)]] = x.after.sp@cluster
    
}

saveRDS(res.clusters,'res.clusters.rds')

names(res.clusters) <- paste0('r',seq(0,1,0.1))





# ##pairwise confusion map by pheatmap for two clustering method
# tab <- table(louvain=res.clusters[["R-igraph"]], leiden=res.clusters[["leiden"]])
# tab <- tab/rowSums(tab)
# pheatmap::pheatmap(tab, color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)


########multiple plot by clustree in one  and try to choose the most stable (shared) clusters#########
#names(res.clusters) <- c('m1','m2')
options(repr.plot.height=7.5,repr.plot.width=15)
res.clusters.df <- as.data.frame(do.call(cbind,res.clusters) )
#combined <- cbind(k.50=clust.50, k.10=clust, k.5=clust.5)
clustree(res.clusters.df, prefix="r", edge_arrow=FALSE,node_text_size = 5)
#clustree(res.clusters.df, prefix="m", edge_arrow=FALSE)
#clustree(res.clusters.df, prefix="k", edge_arrow=FALSE)
#use k =30


###################clustering evaluation############
##https://www.datanovia.com/en/lessons/cluster-validation-statistics-must-know-methods/

##by Dunn index (Preissl et al 2018) and compared to shuffed cell set to verify independent of confounding factors, need to be highest###
res.dunn <- cluster.stats(dist(x.after.sp@smat@dmat),as.numeric(as.character(x.after.sp@cluster)))
res.dunn$dunn
res.dunn$dunn2 
##Preissl'paper choose the clustering profile with the highest dunn index


## get all Dunn index for all res.clusters repeats

d <- dist(x.after.sp@smat@dmat)
501922086


res.clusters[['r0']] <- NULL

#very slow, skip
res.dunn <- sapply(names(res.clusters),FUN = function(x){
  res.dunn <- cluster.stats(d,as.numeric(as.character(res.clusters[[x]])))
  cat(paste('dunn:',round(res.dunn$dunn,digits = 2),'\t','dunn2:',round(res.dunn$dunn2,digits = 2),'\n',sep=''   ) )
   c(round(res.dunn$dunn,digits = 2),round(res.dunn$dunn2,digits = 2))
} )



# r0.8'dunn:0.027292708162731 dunn2:0.495886348494452'
# r0.85'dunn:0.0298403011000761 dunn2:0.610401844162175'
# r0.88'dunn:0.0323365714363431 dunn2:0.605350233520075'
# r0.9'dunn:0.0271804499307676 dunn2:0.538559905826147'
# r0.92'dunn:0.0284082086046624 dunn2:0.566460677539086'
# r0.95'dunn:0.0258828999975614 dunn2:0.536055764743704'
# r0.82'dunn:0.0137925786064152 dunn2:0.53698310301756'

#0.0335325929609384
#0.61000495031456

##0.028 #0.0265 in online example
#0.61

# ##res = 0.9, repeat1-8
# dunn:0.03	dunn2:0.58
# dunn:0.04	dunn2:0.58
# dunn:0.03	dunn2:0.62 #use repeat 3
# dunn:0.03	dunn2:0.62
# dunn:0.01	dunn2:0.54
# dunn:0.03	dunn2:0.62
# dunn:0.02	dunn2:0.6
# dunn:0.02	dunn2:0.59


#saveRDS(res.clusters,'res.clusters.tunning.rds')
#saveRDS(res.clusters,'res.clusters.tunning_res0.9.rds')
#saveRDS(res.clusters,'res.clusters.tunning_res0.9.rep1-8.rds')

#res.choose = 'r0.88'
#res.choose = 'r0.85'
res.choose = 'r0.9'

##rep.choose = 'r3' #choose repeat 3 of res = 0.9
#rep.choose = 'r7' #choose repeat 7 of res = 0.9
#rep.choose = 'r1' #choose repeat 1 of res = 0.9

######by Performing the calculations of Silhouette coefficient by approxSilhouette on the @smat@dmat  (copy from OSCA book)
## based on distance to centroid?
#sil.approx <- approxSilhouette(x.after.sp@smat@dmat, clusters=res.clusters[[rep.choose]] )
sil.approx <- approxSilhouette(x.after.sp@smat@dmat, clusters=res.clusters[[res.choose]] )

sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, res.clusters[[res.choose]], sil.data$other))
sil.data$cluster <- factor(res.clusters[[res.choose]])
#sil.data$closest <- factor(ifelse(sil.data$width > 0, res.clusters[[rep.choose]], sil.data$other))
#sil.data$cluster <- factor(res.clusters[[rep.choose]])

options(repr.plot.height=5,repr.plot.width=5)
ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley",size=0.01) +
    scale_color_manual(values = color_snap) +
    ggtitle(paste('res.choose:',res.choose)) + #,'dunn:','0.03')) +
    theme_bw()+
    theme(legend.position = 'none')
ggsave('cluser_evaluation_silhouette.pdf',height=5,width=5)


#look into c1 negative Silhouette value for mix with which?
table(sil.data$width > 0)

FALSE  TRUE 
 9603 22081 

FALSE  TRUE 
 2890 11466  #79.8% Silhouette coefficent above 0

table(sil.data$width > 0.15)

FALSE  TRUE 
27445  4239

FALSE  TRUE 
 8237  6119 

# table(subset(sil.data,cluster == '1' & width < 0)[,'closest'])
#   1   2   3   4   5   6   7   8    9  10  11  12  13  14  15 
#   0  95   4  195 217  2  266 169   0   0   0   0   0   0   0 

# table(subset(sil.data,cluster == '7' & width < 0)[,'closest'])
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
#   6  75  12   9 141   1   0  54   0   0   0   0   0   0   0 

##plot all Silhouette coefficent
for(rep in paste0('r',1:8)   ){
    cat("Silhouette coefficient plot for repeat ",rep,'\n')
    sil.approx <- approxSilhouette(x.after.sp@smat@dmat, clusters=res.clusters[[rep]] )
#sil.approx <- approxSilhouette(x.after.sp@smat@dmat, clusters=res.clusters[[res.choose]] )

    sil.data <- as.data.frame(sil.approx)
    sil.data$closest <- factor(ifelse(sil.data$width > 0, res.clusters[[rep]], sil.data$other))
    sil.data$cluster <- factor(res.clusters[[rep]])

    options(repr.plot.height=5,repr.plot.width=5)
    p <- ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
        ggbeeswarm::geom_quasirandom(method="smiley",size=0.01) +
        scale_color_manual(values = color_good) +
        ggtitle(paste('rep.choose:',rep,'dunn:',res.dunn[1,rep])) +
        theme_bw()+
        theme(legend.position = 'none')

    print(p)    
    
}





###plot all dotDistri of each res.clusters repeat
dotDistri_new = function (cluster = NULL, id = NULL,rep=NULL){
    colids = colnames(cluster) #must 'cluster','dim1','dim2'
    cluster.sel = cluster[ cluster[,1] == id,]
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    color = 'red'
    #color = ifelse('cluster_sg' == id,'red','navy')
    plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,' res=',rep,sep=''),cex.main = 2.25,xaxt = 'n' ) 
    points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
    return(paste("cluster ",id," ok",sep='') )
}

for(rep in names(res.clusters)){
    cat('rep:',rep,'\n')
    par(mfrow=c(4,4))
    options(repr.plot.height=15,repr.plot.width=15)
    for(i in levels(res.clusters[[rep]]) ){
        
      df <- data.frame(cluster=res.clusters[[rep]],x.after.sp@umap)
      colnames(df) <- c('cluster','UMAP_1','UMAP_2')
      dotDistri_new(cluster = df, id = i,rep=rep)
    }
}


###--choose cluster result---###
#x.after.sp@cluster <- res.clusters[["m2"]] #leiden use k=30 (seed default as 10)
#x.after.sp@cluster <- res.clusters[["r0.92"]] #leiden use k=30 (seed default as 10), res = 0.92
#x.after.sp@cluster <- res.clusters[["r0.88"]] #leiden use k=30 (seed default as 10), res = 0.88
#x.after.sp@cluster <- res.clusters[["r0.85"]] #leiden use k=30 (seed default as 10), res = 0.85
##x.after.sp@cluster <- res.clusters[["r0.9"]] #leiden use k=30 (seed default as 10), res = 0.9
x.after.sp@cluster <- res.clusters[["r3"]] #leiden use k=30 (seed default as 10), res = 0.9 repeat 3
##x.after.sp@cluster <- res.clusters[["r7"]] #leiden use k=30 (seed default as 10), res = 0.9 repeat 7
##x.after.sp@cluster <- res.clusters[["r1"]] #leiden use k=30 (seed default as 10), res = 0.9 repeat 1
#saveRDS(object = x.after.sp@cluster,'leiden.clusters.r0.92.rds')
#saveRDS(object = x.after.sp@cluster,'leiden.clusters.r0.88.rds')
#saveRDS(object = x.after.sp@cluster,'leiden.clusters.r0.85.rds')
saveRDS(object = x.after.sp@cluster,'leiden.clusters.r0.9.repeat3.rds')
##saveRDS(object = x.after.sp@cluster,'leiden.clusters.r0.9.repeat7.rds')
##saveRDS(object = x.after.sp@cluster,'leiden.clusters.r0.9.repeat1.rds')
#saveRDS(object = x.after.sp@cluster,'leiden.clusters.rds')

#res = 0.88
#res = 0.85
#res = 0.82

res = 0.9


x.after.sp@cluster <- res.clusters[['r0.9']]


###plot the fixed umap and cluster
options(repr.plot.height=7.5,repr.plot.width=7.5)

plotViz.new( #with text label halo, point size = 1
    obj=x.after.sp,
    method="umap", 
   # main=paste("snapATAC harmony dims=1:16 K=",k," UMAP\n seed:",seed,' res:',r,' repeat:',rep,' cluster:',cluster_method,sep=''),
    main=paste("snapATAC harmony dims=1:25 K=",k," UMAP\n seed:",seed,' res:',res,' cluster:',cluster_method,sep=''),
    point.color=x.after.sp@cluster, 
    point.size=0.3, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    #down.sample=10000,
    legend.add=FALSE,
    xlim = c(-8,10),
    ylim = c(-8,10)
    );








##############add metadata to x.after.sp from barcodes df (this need do only once)############


##get snap object metadata
metadata <- x.after.sp@metaData[,1:6] 

saveRDS(metadata,'metadata.ori.rds')


#length(levels(metadata$barcode) )
metadata$barcode <- as.character(metadata$barcode)

##bmat == metadata barcode == x.after.sp@barcode
all.equal(  rownames(x.after.sp@bmat), metadata$barcode, check.attributes = FALSE  )#TRUE
all.equal(x.after.sp@barcode, metadata$barcode) #TRUE


##compare to barcodes df and barcode.list
all.equal(barcodes,Reduce(rbind,barcode.list)) #TRUE, barcode.list come from barcodes df
all.equal(x.after.sp@barcode, barcodes$barcode) #TRUE, x.after.sp qc by barcodes df



metadata$sample <- x.after.sp@sample
metadata$sample <- factor(metadata$sample,levels=sample.names)

rownames(metadata) <- paste0(metadata$sample,"#", metadata$barcode )

sum(duplicated(rownames(metadata))) #0

barcode.df <- barcodes #make a new df
rownames(barcode.df) <- rownames(metadata)





##check before combine two dataframe
all.equal(barcode.df$barcode, metadata$barcode) #TRUE
all.equal(rownames(metadata), rownames(barcode.df) )


sum(is.na(barcode.df)) #0
sum(is.na(metadata)) #0


metadata.add <- cbind.data.frame(metadata,barcode.df)


#remove duplicated barcode column after check equal
all.equal(metadata.add[,1], metadata.add[,8] )#TRUE
metadata.add[,8]  <- NULL 

sum(duplicated(colnames(metadata.add))) #1

which(colnames(metadata.add) == 'sample')
#7, 27

colnames(metadata.add)[27] <- 'sampleid'



##rename sample colname as library after check
#colnames(metadata.add)[length( colnames(metadata.add))] <- "library"
# table(metadata.add$sample == 'placenta_donor1' & metadata.add$sampleid == 'D1' |
#       metadata.add$sample == 'placenta_donor2' & metadata.add$sampleid == 'D2')

table(metadata.add$sample,metadata.add$sampleid)
                 
                    D1   D2   D3   D5   D6   D9
  placenta_donor1 6508    0    0    0    0    0
  placenta_donor2    0 4443    0    0    0    0
  placenta_donor3    0    0 5079    0    0    0
  placenta_donor5    0    0    0 4105    0    0
  placenta_donor6    0    0    0    0 5785    0
  placenta_donor9    0    0    0    0    0 5764

#aligned

all.equal(metadata.add$barcode , as.character(x.after.sp@metaData$barcode)) #TRUE

##replace x.after.sp@metaData with new df

x.after.sp@metaData <- metadata.add


#######################





###########################save clusters and do customized visualization############################



# umap <- x.after.sp@umap
# cluster <- x.after.sp@cluster

# cluster.df <- data.frame(cluster,umap)
# plot(cluster.df[,2:3],pch = 19, cex = 0.1)



all.equal (x.after.sp@barcode,row.names(x.after.sp@bmat) ) #TRUE

#row.names(x.after.sp@umap) = row.names(x.after.sp@bmat)
#row.names(x.after.sp@tsne) = row.names(x.after.sp@bmat)

cluster.df = data.frame(cluster=x.after.sp@cluster,UMAP_1=x.after.sp@umap[,1],UMAP_2=x.after.sp@umap[,2]) 
rownames(cluster.df) = rownames(x.after.sp@metaData)

plot(cluster.df[,2:3],pch = 19, cex = .5)

#all.equal (rownames(cluster.df),rownames(x.after.sp@metaData) )#TRUE
cluster.df.add <- cbind(cluster.df,x.after.sp@metaData)


all.equal( as.character(cluster.df.add$barcode) ,barcodes$barcode,check.attributes = FALSE) #TRUE

all.equal( as.character(cluster.df.add$barcode), sapply(stringr::str_split(string = barcodes$barcode.add, pattern = ':',n = 2 ), function(x){x[2]} ) ) #TRUE


##add cluster_merge (from snapatac2) to cluster_df_add_cstb
all.equal(rownames(cluster.df.add), rownames(cluster.df.add.cstb.snapatac2 )  ) #TRUE
all.equal(rownames(cluster.df.add), rownames(cluster.df.add.snapatac2 )  )

cluster.df.add[,c('leiden','cluster_merge','cell_type') ] <- cluster.df.add.cstb.snapatac2[,c('leiden','cluster_merge','cell_type') ]

cluster.df.add[,c('leiden') ] <- cluster.df.add.snapatac2[,c('leiden') ]

cluster.df.add$leiden <- factor( as.character(cluster.df.add$leiden ) , levels =  as.character(sort(unique(cluster.df.add$leiden )))  )



#cluster.df.add <- cbind.data.frame(cluster.df.add,barcodes)

#sum(duplicated(cluster.df.add$barcode.add)) #0

#rownames(cluster.df.add) <- cluster.df.add$barcode.add



saveRDS(cluster.df.add,"cluster.df.add.rds")
write.table(x = cluster.df.add,file = 'cluster.df.add.txt',quote = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)


#saveRDS(cluster.df.add,"cluster.df.add.cstb.rds")
cluster.df.add <- readRDS("cluster.df.add.cstb.rds")
#write.table(x = cluster.df.add,file = 'cluster.df.add.cstb.txt',quote = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)


saveRDS(cluster.df.add,"cluster.df.add.final.rds")
write.table(x = cluster.df.add,file = 'cluster.df.add.final.txt',quote = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)




##annotate by previous d1 d2 cluster
cluster.df.add.late1_late2_filtered <- read.table('../../placenta_10X_late1/02.snapatac2/cluster.df.add.late1_late2_filtered.txt',header=TRUE,sep='\t',stringsAsFactors = FALSE,comment.char = "")
#8334 x 34

#placenta_donor2#
#placenta_10X_late6:

cellid <- rownames(cluster.df.add.late1_late2_filtered)

cellid <- gsub(pattern = "placenta_donor",replacement = "placenta_10X_late",x=cellid )
cellid <- gsub(pattern = "#",replacement = ":",x=cellid )

cluster.df.add.late1_late2_filtered$barcode.add <- cellid

table(rownames(cluster.df.add.late1_late2_filtered) %in% rownames(cluster.df.add) )
FALSE  TRUE 
  679  7655

TRUE 
8334 
table(cluster.df.add.late1_late2_filtered$barcode.add %in% cluster.df.add$barcode.add )

FALSE  TRUE 
  679  7655

TRUE 
8334

# table(cluster.df.add.late1_late2_filtered$barcode_add %in% rownames(cluster.df.add) )
# FALSE  TRUE 
#   594  7740

cluster.df.add.late1_late2_filtered.sel <- cluster.df.add.late1_late2_filtered[rownames(cluster.df.add.late1_late2_filtered) %in% rownames(cluster.df.add), ]

#rownames(cluster.df.add.late1_late2_filtered.sel) <- cluster.df.add.late1_late2_filtered.sel$barcode_add

# cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == 6) )
# cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == 7) )
# cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == 2) )
# cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == 4) )
# cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == 3) )
# cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == 8) )

# cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == 5) )
# cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == 1) )


options(repr.plot.width = 7.5,repr.plot.height = 15)
par(mfrow = c(4,2))
for(i in c(8,2,5,3,1,4,7,6)){
  #cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered,cluster == i) )
  cellid_sel <- rownames(subset(cluster.df.add.late1_late2_filtered.sel,cluster == i) )
  plot(cluster.df.add[,c('UMAP_1','UMAP_2')],pch=16,cex=0.1,xlim = c(-8,10),ylim = c(-8,10), main=i )
  points( cluster.df.add[cellid_sel,c('UMAP_1','UMAP_2')] ,pch=16,cex=0.1,xlim = c(-10,10),ylim = c(-30,10),col='red'   )

}

options(repr.plot.width = 7.5,repr.plot.height = 7.5)
#plot(cluster.df.add.late1_late2_filtered[,c('UMAP_1','UMAP_2')],pch=16)
plot(cluster.df.add.late1_late2_filtered.sel[,c('UMAP_1','UMAP_2')],pch=16)





############customized way to plot umap-cluster with text halo###########

##define cell name and color see below 

# centers <- cluster.df.add %>% dplyr::group_by(cellname) %>% dplyr::summarize(x = median(x = UMAP_1), 
#         y = median(x = UMAP_2))

centers <- cluster.df.add %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
       y = median(x = UMAP_2))

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
                                  #cluster=as.character(unlist(centers[i,'cellname'])),
                                  cluster=as.character(unlist(centers[i,'cluster'])),
                                  x=centers[i,'x'] + cos(j)*xo, 
                                  y=centers[i,'y'] + sin(j)*yo
                                 )
                       )
      }
}



xmin = min(cluster.df.add[,'UMAP_1'])
xmax = max(cluster.df.add[,'UMAP_1'])
ymin = min(cluster.df.add[,'UMAP_2'])
ymax = max(cluster.df.add[,'UMAP_2'])

shrink.x = 0.3
shrink.y = 0.1

cluster.df.add.shuffle <- cluster.df.add[sample(1:nrow(cluster.df.add),size=nrow(cluster.df.add),replace=FALSE ),]


####the UMAP plot with annotation
#label right
options(repr.plot.height=5,repr.plot.width=6)
ggplot(cluster.df.add.shuffle,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = color_good)  +
  scale_colour_manual(values = color_snap_mod1)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  
  theme(
        legend.position = 'right',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ylim( ymin-shrink.y*abs(ymin-ymax),ymax+shrink.y*abs(ymin-ymax) ) + 
  xlim( xmin-shrink.x*abs(xmin-xmax),xmax+shrink.x*abs(xmin-xmax) ) +
  ggtitle(paste(sample, "six donors, \ntotal cells:",nrow(cluster.df.add),  sep=" ") ) +
#   geom_text(data = centers_shift, #the halo
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "white", 
#             size = 6.5) +
#   geom_text(data = centers, 
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "black", 
#             size = 6) +
  guides(col = guide_legend(override.aes = list(size = 3))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/UMAP/PLA-late-combine-ATAC-UMAP.pdf",height=5,width=6,useDingbats=FALSE)
ggsave(filename = "pdfs/UMAP/PLA-late-combine-ATAC-UMAP.merge.pdf",height=5,width=6,useDingbats=FALSE)

##label on cluster
options(repr.plot.height=5,repr.plot.width=5.5,repr.plot.res = 150)
ggplot(cluster.df.add.shuffle,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cellcolor  )) +
  geom_point(size = 0.2,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
  #scale_colour_manual(values = unlist(map_cellcolor) )  +
  ##scale_colour_manual(values = color_snap_mod1)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  theme(
        legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "six donors, \ntotal cells:",nrow(cluster.df.add),  sep=" ") ) +
  geom_text(data = centers_shift, #the halo
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "white", 
            size = 4.5) +
            #size = 6.5) +
  geom_text(data = centers, 
            #mapping = aes(x=x,y=y,label = cellname), 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            size = 4.5) +
            #size = 6.5) +
  ##guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  #ylim(-5,5) + xlim(-2,6) +  ##will overwrite scale_x/y_continuous
  ylim( ymin-shrink.y*abs(ymin-ymax),ymax+shrink.y*abs(ymin-ymax) ) + 
  xlim( xmin-shrink.x*abs(xmin-xmax),xmax+shrink.x*abs(xmin-xmax) ) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/UMAP/PLA-late-combined-ATAC-UMAP.labelon.pdf",height=5,width=5.5)
ggsave(filename = "pdfs/UMAP/PLA-late-combined-ATAC-UMAP.labelon.merge.pdf",height=5,width=5.5)

###############


##by donors
options(repr.plot.height=7.5,repr.plot.width=7.5)
ggplot(cluster.df.add.shuffle,aes(x=UMAP_1,y=UMAP_2,col= sample  )) +
  geom_point(size = .1,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_manual(values = color_good)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Source of donor ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(-10,10) + xlim(-10,10) +  ##will overwrite scale_x/y_continuous
  ylim( ymin-shrink.y*abs(ymin-ymax),ymax+shrink.y*abs(ymin-ymax) ) + 
  xlim( xmin-shrink.x*abs(xmin-xmax),xmax+shrink.x*abs(xmin-xmax) ) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/qc/PLA-late-combine-ATAC-source-of-donor.pdf",height=15,width=15)
####


###by depth
#options(repr.plot.height=5.5,repr.plot.width=5)
options(repr.plot.height=5.5,repr.plot.width=5)
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= passed_filters)) +
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= log(UQ) )) +
ggplot(cluster.df.add.shuffle,aes(x=UMAP_1,y=UMAP_2,col= logUMI )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = rev(color_cellranger))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'D'))  + 
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Sequence depth \n(log fragment count)",  sep=" ") ) +
  #ggtitle(paste(sample, "Sequence depth ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  ylim( ymin-shrink.y*abs(ymin-ymax),ymax+shrink.y*abs(ymin-ymax) ) + 
  xlim( xmin-shrink.x*abs(xmin-xmax),xmax+shrink.x*abs(xmin-xmax) ) +
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/qc/PLA-late-combine-ATAC-depth.pdf",height=5.5,width=5)


##by TSS score
##see following 



##by FRiP (promoter percentage)
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add.shuffle,aes(x=UMAP_1,y=UMAP_2,col= promoter_ratio )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = rev(color_cellranger))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'D'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Promoter ratio",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  ylim( ymin-shrink.y*abs(ymin-ymax),ymax+shrink.y*abs(ymin-ymax) ) + 
  xlim( xmin-shrink.x*abs(xmin-xmax),xmax+shrink.x*abs(xmin-xmax) ) +
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/qc/PLA-late-combine-ATAC-promotorRatio.pdf",height=5.5,width=5)


########################





#########split plot distribution of clusters (to evaluate the separation trend of clusters)####

##plot distribution of each clusters
dotDistri = function (cluster = NULL, id = NULL){
    colids = colnames(cluster) #must 'cluster','dim1','dim2'
    cluster.sel = cluster[ cluster[,1] == id,]
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    color = 'red'
    #color = ifelse('cluster_sg' == id,'red','navy')
    plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
    points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
    return(paste("cluster ",id," ok",sep='') )
}


#pdf( "pdfs/PLA-term-ATAC-cluster.distri1.pdf",height=15,width=15,useDingbats = FALSE)
par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('1','2','3','4','5','6','7','8','9') ){
  dotDistri(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i)
}
#dev.off()

#pdf( "pdfs/PLA-term-ATAC-cluster.distri2.pdf",height=15,width=15,useDingbats = FALSE)
par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
#for(i in c('12','10','13','15','11','14') ){
#for(i in c('9','11','10','14','12','15') ){
#for(i in c('11','10','13','15','14','12') ){
#for(i in c('12','10','13','15','14','11') ){
#for(i in c('13','11','9','14','12','15') ){
#for(i in c('12','13','10','15','11','14') ){
for(i in c('10','11','12','13','14','15','16') ){
  dotDistri(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i)
}
#dev.off()


#saveRDS(x.after.sp,"x.after.sp.tuning.rds")



######plot cluster count barplot, use SnapATAC object colors

# colPanel = SnapATAC:::createColorPanel(length(unique(x.after.sp@cluster)))

# color = scales::alpha(colPanel[factor(unique(x.after.sp@cluster))], 0.8)

table(cluster.df.add$cluster)
1    2    3    4    5    6    7    8  #cstb merge
3457 3583 6282 1449 2983 5465 1050  423

1    2    3    4    5    6    7    8    9   10 #cstb
3457 3583 3186 3096 2983 2814 2651 1449 1050  423

 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
3639 3595 3353 3259 3139 2962 2791 2528 1755 1526 1106  959  447  234  213  178 

#    1    2    3    4    5    6    7    8    9   10   11   12 
# 1801 1749 1539 1044 1048  613  441  346  375  194  127   59 

#  1    2    3    4    5    6    7    8    9   10   11   12   13 
# 1801 1749 1539 1044 1048    0  613  441  346  375  194  127   59 

#    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 1880 1829 1624 1309 1209  737  684  539  385  375  194  127   59 

#    1    2    3    4    5    6    7    8    9   10   11   12 
# 5720 4228 3024 2844 2011  747  691  566  526  281  192  105

#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1832 1783 1090 1546 1475 1509 1019  444  395  352  290  253  235  171  132 

#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 2025 1826 1756 1669 1570 1544 1280  684  428  396  293  291  262  186  146 

#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 2103 1955 1786 1768 1456 1217 1098  758  659  395  286  282  263  186  144 

#color = color_good
color = color_snap_mod1

options(repr.plot.width=2,repr.plot.height=7)
barplot(table(cluster.df.add$cluster),horiz = TRUE,col = color )

stb_sum <- sum(table(cluster.df.add$cluster)[c('9','5','7','1','2','8','4')]) #8494 #9482 #8835
stb_sum/sum(table(cluster.df.add$cluster)) #12526 #14356
8494/12526 ~= 67.8%
#9482/14356 ~= 66%
#8835/14356 ~= 61.5%


# stb_sum <- sum(table(cluster.df.add$cluster)[c(18,3,10,4,6,7,13,19)]) #8371
# stb_sum/sum(table(cluster.df.add$cluster)) #14148
# 8471/14148 ~= 60%


##############################################




####save object before marker gene 


#saveRDS(x.after.sp,"x.after.sp.onestep_filter.rds")

#saveRDS(x.after.sp,"x.after.sp.rds")

#x.after.sp <- readRDS("x.after.sp.rds")



############annotate clusters with signac cisTopic episcanpy cicero, cellranger-atac#############
# #read in meta and clusters
# combined = readRDS("../../04.scRNA_scATAC/liger/all.atac.methods.cluster.tab.withmeta.rds") #7416
# which (is.na.data.frame(combined) ) #will return row idx with NA in any columns

# id.shared = intersect(rownames(cluster),rownames(combined)) #6772 of 6820 , 48 cells less

# ##combine with snapATAC clusters
# combined.addsn = cbind(cluster[id.shared,],combined[id.shared,]) #6772 of 6820
# #all.equal( rownames(combined.addsn), rownames(cluster)  ) #diffs
# which (is.na.data.frame(combined.addsn) ) #has na #will return row idx with NA in any columns
# combined.addsn$cluster.ct = factor(as.numeric(combined.addsn$cluster.ct)  ) #will -1

# colnames(combined.addsn)[1:5] = c("cluster.sn",'UMAP_1.sn',"UMAP_2.sn",'tSNE_1.sn','tSNE_2.sn')


# ####plot meta data#####

# cellDecorate = function(data=NULL,coord1=NULL,coord2=NULL,value=NULL,title=NULL){
#   p=ggplot(data=data,aes_string(x=coord1,y=coord2,col=value ) ) +
#     geom_point(size=0.5,alpha=0.8) +
#     #scale_color_manual(values = color_good) +
#     scale_color_gradientn(colours = rainbow(n = 7,rev = TRUE)) +
#     labs(title = title, subtitle="",caption = '',tag="",label="") +
#     xlab (coord1) +
#     ylab (coord2) +
#     #scale_x_continuous(limits = c(-10,10)) +
#     #scale_y_continuous(limits = c(-10,10)) +
#     #ggtitle(title) +
#     theme_classic() #+
#     #guides( col = guide_legend(override.aes = list(size = 5)) )
#  p
# }  


# ###sequencing depth
# cellDecorate(data=combined.addsn,coord1="UMAP_1.sn",coord2="UMAP_2.sn",value="log10(passed_filters+1)",title='snapATAC UMAP sequencing depth')

# ###nucleosome signal
# cellDecorate(data=combined.addsn,coord1="UMAP_1.sn",coord2="UMAP_2.sn",value="nucleosome_signal",title='snapATAC UMAP nucleosome_signal')


# ##TSS_enrichment score
# cellDecorate(data=combined.addsn,coord1="UMAP_1.sn",coord2="UMAP_2.sn",value="TSS.enrichment",title='snapATAC UMAP TSS.enrichment')


# ################cross cluster comparison############

# plotClusterDist = function(combined.add=combined.add){

# ##wrapper for plot distribution of each clusters
#   dotDistri = function (cluster = NULL, id = NULL){
#     colids = colnames(cluster) #must 'cluster','dim1','dim2'
#     cluster.sel = cluster[ cluster[,1] == id,]
#     n_sel = nrow(cluster.sel)
#     #cat ('select for cluster ',id,' n = ',n_sel," \n")
#     color = 'red'
#     #color = ifelse('cluster_sg' == id,'red','navy')
#     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25 ) 
#     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
#     return(paste("cluster ",id," ok",sep='') )
#   }

# #signac
#   par(mfrow=c(3,3))
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '1')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '10')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '3')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '8')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '2')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '0')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '4')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '5')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '6')

#   par(mfrow=c(3,3))
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '7')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '14')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '9')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '11')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '12')
#   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '16')
# #   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '14')
# #   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '9')
# #   dotDistri(cluster = combined.add[,c('cluster.sg','UMAP_1.sn','UMAP_2.sn')], id = '12')

#   par(mfrow=c(3,3))
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '1')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '10')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '3')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '8')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '2')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '0')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '4')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '5')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '6')

#   par(mfrow=c(3,3))
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '7')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '14')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '9')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '11')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '12')
#   dotDistri(cluster = combined.add[,c('cluster.sg','tSNE_1.sn','tSNE_2.sn')], id = '16')
    
    
# 	#cistopic
#   par(mfrow=c(3,3))
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '5')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '16')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '21')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '18')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '8')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '2')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '6')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '1')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '10')

#   par(mfrow=c(3,3))
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '4')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '13')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '19')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '20')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '3')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '15')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '14')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '9')
#   dotDistri(cluster = combined.add[,c('cluster.ct','UMAP_1.sn','UMAP_2.sn')], id = '12')

# 	#cellranger-atac
# par(mfrow=c(3,3))
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '0')
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '11')
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '2')
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '8')
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '9')
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '5')
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '3')
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '1')
#   dotDistri(cluster = combined.add[,c('cluster.cr','UMAP_1.sn','UMAP_2.sn')], id = '7')
	
# 	return(" ")

# }

# ##by snapatac cluster
# ggplot(data=combined.addsn,aes(x=UMAP_1.sn,y=UMAP_2.sn,col=cluster.sg) ) +
#   geom_point(size=1) +
#   scale_color_manual(values = color_good) +
#   #xlim(-10,10) +
#   #ylim(-10,10) +
#   theme_classic() +
#   guides( col = guide_legend(override.aes = list(size = 10)) )


# ##plot each clusters
# plotClusterDist(combined.add=combined.addsn)

# # #################

# saveRDS(x.sp,"x.sp.nonorm_jmat.rds")


# x.sp <- readRDS("x.sp.nonorm_jmat.rds")


# ########################################jump these steps in this script#######################

# ##################step 9 scRNA based annotation (CCA)##################
# library(Seurat);
# placenta.rna = readRDS("../../03.seurate/PLA-term-RNASEQ-1/placenta.PLA-term-RNA-1.rds");
# placenta.rna$tech = "rna";
# variable.genes = VariableFeatures(object = placenta.rna); #2000
# genes.df = read.table("genes.gtf.symbol.bed");
# genes.gr = GRanges(genes.df[,1], 
# 				   IRanges(genes.df[,2], genes.df[,3]), 
# 				   name=genes.df[,4]);
# genes.sel.gr = genes.gr[which(genes.gr$name %in% variable.genes)]; #1609

# ## reload the bmat, this is optional but highly recommanded
# x.sp = addBmatToSnap(x.sp);
# saveRDS(x.sp,"x.sp.stage4.use10xmat.rds")

# x.sp = createGmatFromMat(
#     obj=x.sp, 
#     input.mat="bmat",
#     genes=genes.sel.gr,
#     do.par=TRUE,
#     num.cores=10
#   );

# ###convert the snap object to Seurat object in preparation of integration.
# placenta.atac <- snapToSeurat(
#     obj=x.sp, 
#     eigs.dims=1:20, 
#     norm=TRUE,
#     scale=TRUE
#   );
# transfer.anchors <- FindTransferAnchors(
#     reference = placenta.rna, 
#     query = placenta.atac, 
#     features = variable.genes, 
#     reference.assay = "RNA", 
#     query.assay = "ACTIVITY", 
#     reduction = "cca"
#   );
# celltype.predictions <- TransferData(
#     anchorset = transfer.anchors, 
#     refdata = placenta.rna$seurat_clusters,
#     weight.reduction = placenta.atac[["SnapATAC"]],
#     dims = 1:20
#   );
# x.sp@metaData$predicted.id = celltype.predictions$predicted.id;
# x.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max);
# x.sp@cluster = as.factor(x.sp@metaData$predicted.id);


# ############step 10 create psudo multiomics cells#######
# refdata <- GetAssayData(
#     object = placenta.rna, 
#     assay = "RNA", 
#     slot = "data"
#   );
# imputation <- TransferData(
#     anchorset = transfer.anchors, 
#     refdata = refdata, 
#     weight.reduction = placenta.atac[["SnapATAC"]], 
#     dims = 1:20
#   );
# x.sp@gmat = t(imputation@data);
# #rm(imputation); # free memory
# #rm(refdata);    # free memory
# #rm(placenta.rna);   # free memory
# #rm(placenta.atac); # free memory

# ######step 11, remove cells of low prediction score#####
# hist(
#     x.sp@metaData$predict.max.score, 
#     xlab="prediction score", 
#     col="lightblue", 
#     xlim=c(0, 1),
#     main="placenta 10X"
#   );
# abline(v=0.3, col="red", lwd=2, lty=2);
# table(x.sp@metaData$predict.max.score > 0.3);
# #FALSE  TRUE 
# # 1256  6160

# #> 0.5
# #FALSE  TRUE 
# # 4046  3370 

# x.sp = x.sp[x.sp@metaData$predict.max.score > 0.3,];
# x.sp
# saveRDS(x.sp,"x.sp.stage5.use10xmat.rds")

# plotViz(
#     obj=x.sp,
#     method="umap", 
#     main="placenta scRNA to scATAC",
#     point.color=x.sp@metaData[,"predicted.id"], 
#     point.size=0.8, 
#     point.shape=19, 
#     text.add=TRUE,
#     text.size=2,
#     text.color="black",
#     down.sample=10000,
#     legend.add=FALSE
#   );


#   dotDistri = function (cluster = NULL, id = NULL){
#     colids = colnames(cluster) #must 'cluster','dim1','dim2'
#     cluster.sel = cluster[ cluster[,1] == id,]
#     n_sel = nrow(cluster.sel)
#     #cat ('select for cluster ',id,' n = ',n_sel," \n")
#     color = 'red'
#     #color = ifelse('cluster_sg' == id,'red','navy')
#     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25 ) 
#     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
#     return(paste("cluster ",id," ok",sep='') )
#   }
  
#  combined.cca = cbind('cluster'=x.sp@metaData[,"predicted.id"],x.sp@umap)
 
# ##scRNA
#   par(mfrow=c(3,3))
#   dotDistri(cluster = combined.cca, id = '7')
#   dotDistri(cluster = combined.cca, id = '8')
#   dotDistri(cluster = combined.cca, id = '6')
#   dotDistri(cluster = combined.cca, id = '10')
#   dotDistri(cluster = combined.cca, id = '2')
#   dotDistri(cluster = combined.cca, id = '1')
#   dotDistri(cluster = combined.cca, id = '9')
#   dotDistri(cluster = combined.cca, id = '3')
#   dotDistri(cluster = combined.cca, id = '4')
	
#   par(mfrow=c(2,2))
#   dotDistri(cluster = combined.cca, id = '0')
#   dotDistri(cluster = combined.cca, id = '5')
#   dotDistri(cluster = combined.cca, id = '11')
#   dotDistri(cluster = combined.cca, id = '12')
		

# ##step 12 gene expression projecte onto umap


# marker.genes = c(
#     "CGA", "PSG1",
#     "DNMT1", "PLAC8", "FLT1", 
#     "VIM", "COL1A1", "FN1", 
#     "HLA-G", "PAPPA2", "CD14" 
# #     "AIF1", "DKK1", "IGFBP1", 
# #      "CNN1", "MYH11", 
# #     "HBB", "HBG1", 
# #     "GZMA"
#   );
# #marker.genes[which(!marker.genes %in% colnames(x.sp@gmat) )]
# par(mfrow = c(3, 3));
# for(i in 1:9){
#     j = which(colnames(x.sp@gmat) == marker.genes[i])
#     plotFeatureSingle(
#         obj=x.sp,
#         feature.value=x.sp@gmat[,j],
#         method="umap", 
#         main=marker.genes[i],
#         point.size=0.1, 
#         point.shape=19, 
#         down.sample=10000,
#         quantiles=c(0.01, 0.99)
#  )};
# par(mfrow = c(3, 3));
# for(i in 10:19){
#     j = which(colnames(x.sp@gmat) == marker.genes[i])
#     plotFeatureSingle(
#         obj=x.sp,
#         feature.value=x.sp@gmat[,j],
#         method="umap", 
#         main=marker.genes[i],
#         point.size=0.1, 
#         point.shape=19, 
#         down.sample=10000,
#         quantiles=c(0.01, 0.99)
#  )};


# #################################jump these steps in this script################################



###step 12 calculate the gmat and annotation cluster with marker gene activity####
genes = read.table("genes.gtf.symbol.bed") #cellranger-atac gtf, gencode v28, ensembl v92
genes <- genes[!grepl(pattern = '^MT-',x = genes[,4]) ,]
genes.gr = GRanges(genes[,1], 
    IRanges(genes[,2], genes[,3]), name=genes[,4]
  ) #without promoter 





#genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)]

########re-add bmat to object (to use the raw bin count instead of binary count)


x.sp.stage0 <- readRDS("x.sp.stage0.rds")

#x.after.sp = addBmatToSnap(x.after.sp,bin.size = 5000) #will cause error, snap file bin total size diff

##add bmat raw from x.sp.stage0

all.equal( x.after.sp@barcode,x.sp.stage0@barcode ) #TRUE

table(mcols(x.after.sp@feature)$name %in% mcols(x.sp.stage0@feature)$name)
TRUE 
545479

idy <- match(mcols(x.after.sp@feature)$name , mcols(x.sp.stage0@feature)$name)
x.sp.stage0.sel <- x.sp.stage0[, idy, mat="bmat"]

x.sp.stage0.sel
number of barcodes: 31684
number of bins: 545479
number of genes: 0
number of peaks: 0
number of motifs: 0


all.equal(x.after.sp@feature , x.sp.stage0.sel@feature) #TRUE


all.equal( x.after.sp@barcode,x.sp.stage0.sel@barcode )

x.after.sp@bmat <- x.sp.stage0.sel@bmat
max(x.after.sp@bmat) #14


saveRDS(x.after.sp@graph,'graph.rds')
saveRDS(x.after.sp@smat,'smat.rds')
saveRDS(x.after.sp@metaData, 'metaData.rds')

##add raw bmat done




#14356 58344
x.after.sp = createGmatFromMat(
    obj=x.after.sp, 
    input.mat="bmat",
    #genes=genes.sel.gr,
    genes=genes.gr, #calculate for all genes, very slow
    do.par=TRUE,
    num.cores=20
  )

# normalize the cell-by-gene matrix
x.after.sp = scaleCountMatrix( #very slow for ~20000, quick for 9336
    obj=x.after.sp, 
    cov=x.after.sp@metaData$passed_filters + 1,
    mat="gmat",
    method = "RPM"
  )
# smooth the cell-by-gene matrix
x.after.sp = runMagic( #will take 1-2 days!!
    obj=x.after.sp,
    input.mat="gmat",
    step.size=3
  )


##saveRDS(x.after.sp,"x.after.sp.addgmat.onestep_filter.rds") ##11G

##save gmat of SnapATAC
rownames(x.after.sp@gmat) <- rownames(x.after.sp@metaData)


saveRDS(x.after.sp@gmat,"gmat.rds") #big and slow, 16G


rm(x.sp.stage0.sel)
x.after.sp@gmat <- readRDS('empty.rds')#as(matrix(0),'dgCMatrix')

gc()



#restart from here
#saveRDS(x.after.sp,'x.after.sp.rm_gmat.rds')

x.after.sp <- readRDS('x.after.sp.rm_gmat.rds')
number of barcodes: 31684
number of bins: 545479
number of genes: 1
number of peaks: 0
number of motifs: 0




# # ############calculate gmat for given gene only (still cost much time!, run full gene again in stand alone script)###########

# # x.sp.stage0 <- readRDS("x.sp.stage0.rds")

# # number of barcodes: 31684
# # number of bins: 617669
# # number of genes: 0
# # number of peaks: 0
# # number of motifs: 0


# # all.equal( x.after.sp@barcode,x.sp.stage0@barcode ) #TRUE

# # x.sp.stage0@metaData <- x.after.sp@metaData
# # x.sp.stage0@graph <- x.after.sp@graph

# # genes.gr.sel <- genes.gr[which(mcols(genes.gr)$name %in% c('PAPPA','FLT1','CSH1','CSH2','CSHL1','GH2','MAFK','JUND','CD14')),]

# # x.sp.stage0 = createGmatFromMat(
# #     obj=x.sp.stage0, 
# #     input.mat="bmat",
# #     genes=genes.gr.sel, #calculate for given gene
# #     do.par=FALSE,
# #     num.cores=1
# #   )

# # # normalize the cell-by-gene matrix
# # x.sp.stage0 = scaleCountMatrix( 
# #     obj=x.sp.stage0, 
# #     cov=x.sp.stage0@metaData$passed_filters + 1,
# #     mat="gmat",
# #     method = "RPM"
# #   )
# # # smooth the cell-by-gene matrix
# # x.sp.stage0 = runMagic( #will take 
# #     obj=x.sp.stage0,
# #     input.mat="gmat",
# #     step.size=3
# #   )


# # rownames(x.sp.stage0@gmat) <- rownames(x.sp.stage0@metaData)
# # saveRDS(x.sp.stage0@gmat,'gmat.addgene.rds')


# # rm(x.sp.stage0)
# # gc()

# # ######





# marker.genes = c(
#     "DNMT1", "CDH1", "MKI67",
#     "FLT1", "CSHL1", "PSG8", 
#     "ERVFRD-1", "LAIR2", "PLAC8",
#     'VIM','PECAM1','CD14'
#   );

# marker.genes = c( #try to distinguish STB terminals 
#     "FLT1", 'LEP','INSIG2', #FLT1-enriched early stage
#     "CSHL1",'CSH1','PAPPA',  #PAPPA-enriched late stage
#     'PSG1','CGA','GCM1' #general

#   );

# ##STR markers
# marker.genes = c(
#     "VIM", "PECAM1", "THY1",
#     "NR2F1", "DLK1", "HIVEP3", 
#     "CD68", "CD14", "HLA-A",
#     'HLA-DPA1','HLA-DPB1','MKI67'
#   );

# marker.genes %in% colnames(x.after.sp@gmat)


# ################ageing marker gene##########
# ##https://genomics.senescence.info/genes/microarray.php?show=4&sort=1&page=1
# ageing.OE.genes <- read.table('ageing.gene.OE.txt',header=FALSE,stringsAsFactors = FALSE,sep='\t')
# #overexpressed 56
# ageing.UE.genes <- read.table('ageing.gene.UE.txt',header=TRUE,stringsAsFactors = FALSE,sep='\t')
# #underexpressed  17

# ageing.genes <- ageing.OE.genes[,2]
# ageing.genes <- ageing.UE.genes[,2]

# ageing.genes = c(

#     "CASP3", "CDKN1A", "CDKN2A", #the CDKNxx p16, p21
#     "ANXA3", "ANXA5", "HIST1H1C", 
#     "SERPING1", "TMED10", "TXNIP",
#     'MPEG1','SPP1','LAPTM5'
#   );

# ageing.genes = c('CASP3','HIST1H1C',"TMED10",'TP53')


# ageing.genes %in% colnames(x.after.sp@gmat)
# marker.genes <- ageing.genes[ageing.genes %in% colnames(x.after.sp@gmat)]


# ##Syncytial knot

# marker.genes <- c('SPATA18','CROT','PTCHD4','CDKN1A')



# ##antiapoptotic
# #BCL1 = CCND1
# #CFLAR = FLIP
# #BORIS = CTCFL
# #IAP1 = BIRC3

# marker.genes <- c('BAG1','CCND1','BCL2','MCL1','A1','BIRC3','IAP2','IAP','CTCF','CTCFL','DNMT3A','CFLAR','BACH2','NFIL3')
# marker.genes <- marker.genes[marker.genes %in% colnames(x.after.sp@gmat)]


# #########



# #########plot marker gene scatter plot by SnapATAC #######
# #source('plotFeatureSingle_new.r')
# par(mfrow = c(1, 1));
# options(repr.plot.width=10,repr.plot.height=10)
# # for(i in 1:6){
# #     plotFeatureSingle(
# #     #plotFeatureSingle_new(
# #         obj=x.after.sp,
# #         feature.value=x.after.sp@gmat[, marker.genes[i]],
# #         method="umap", 
# #         main=marker.genes[i],
# #         point.size=0.2, 
# #         point.shape=19, 
# #         down.sample=10000,
# #         quantiles=c(0, 1)
# #   )}

# # for(i in 7:11){
# #     plotFeatureSingle(
# #     #plotFeatureSingle_new(
# #         obj=x.after.sp,
# #         feature.value=x.after.sp@gmat[, marker.genes[i]],
# #         method="umap", 
# #         main=marker.genes[i],
# #         point.size=0.2, 
# #         point.shape=19, 
# #         down.sample=10000,
# #         quantiles=c(0, 1)
# #   )}


# ##use snapATAC plot style###
# for (gene in marker.genes){
#     options(repr.plot.height=7,repr.plot.width=7.5)
#     #options(repr.plot.height=15,repr.plot.width=15.5)
    
#     p<- plotFeatureSingle( #use plot3D, can not append ggplot theme
#     #plotFeatureSingle_new(
#         obj=x.after.sp,
#         #feature.value=x.after.sp@gmat[, gene],
#         feature.value=gmat[, gene],
#         method="umap", 
#         main=gene,
#         point.size=0.2, 
#         point.shape=19, 
#         #down.sample=10000,
#         quantiles=c(0.1, 1)##make sharper
#       )
#     #ggsave(filename = paste('marker_gene.red_blue',gene,'.pdf',sep=''),height=5,width=5.5,useDingbats=FALSE  )

# }



# ######customized plot marker gene scatter plot#########
# ##marker.gmat <- as.data.frame(as(gmat[, marker.genes],'matrix'))
# marker.gmat <- as.data.frame(as(x.after.sp@gmat[, marker.genes],'matrix'))
# rownames(marker.gmat) <- rownames(x.after.sp@metaData)
# all.equal (rownames(cluster.df.add),rownames(marker.gmat) ) #TRUE
# marker.gmat.df <- cbind(cluster.df.add[,c('cluster','UMAP_1','UMAP_2')],marker.gmat )


# ##cutoff quantiles
# #low.q <- 0.1
# #high.q <- 0.9

# low.q <- c('DNMT1'=0.1,'CDH1'=0.1,'PAGE4'=0.1,'FLT1'=0,'CSHL1'=0.1,'PSG8'=0,'ERVFRD-1'=0,'LAIR2'=0,'PLAC8'=0.1,'VIM'=0,'PECAM1'=0,'MKI67'=0,'CD14'=0)
# high.q <- c('DNMT1'=1,'CDH1'=1,'PAGE4'=1,'FLT1'=1,'CSHL1'=1,'PSG8'=1,'ERVFRD-1'=1,'LAIR2'=1,'PLAC8'=1,'VIM'=1,'PECAM1'=1,'MKI67'=1,'CD14'=1)

# low.q <- c('FLT1'=0,'LEP'=0,'INSIG2'=0,'CSHL1'=0,'CSH1'=0,'PAPPA'=0,'PSG1'=0,'CGA'=0,'GCM1'=0)
# high.q <- c('FLT1'=1,'LEP'=1,'INSIG2'=1,'CSHL1'=1,'CSH1'=1,'PAPPA'=1,'PSG1'=1,'CGA'=1,'GCM1'=1)

# low.q <- c('CASP3'=0.1,'CDKN1A'=0.1,'CDKN2A'=0,'ANXA3'=0.1,'ANXA5'=0,'HIST1H1C'=0,'SERPING1'=0,'TMED10'=0,'TXNIP'=0,'MPEG1'=0,'SPP1'=0,'LAPTM5'=0)
# high.q <- c('CASP3'=0.9,'CDKN1A'=0.9,'CDKN2A'=1,'ANXA3'=1,'ANXA5'=1,'HIST1H1C'=1,'SERPING1'=1,'TMED10'=1,'TXNIP'=1,'MPEG1'=1,'SPP1'=1,'LAPTM5'=1)



# low.q <- c('DNMT1'=0.1,'CDH1'=0.1,'PAGE4'=0.1,'FLT1'=0,'CSHL1'=0.1,'PSG8'=0,'ERVFRD-1'=0,'LAIR2'=0,'PLAC8'=0.1,'VIM'=0,'PECAM1'=0,'MKI67'=0,'CD14'=0)
# high.q <- c('DNMT1'=1,'CDH1'=1,'PAGE4'=1,'FLT1'=1,'CSHL1'=1,'PSG8'=1,'ERVFRD-1'=1,'LAIR2'=1,'PLAC8'=1,'VIM'=1,'PECAM1'=1,'MKI67'=1,'CD14'=1)



# low.q <- c("VIM"=0, "PECAM1"=0, "THY1"=0,"NR2F1"=0, "DLK1"=0, "HIVEP3"=0, "CD68"=0, "CD14"=0, "HLA-A"=0,'HLA-DPA1'=0,'HLA-DPB1'=0,'MKI67'=0)
# high.q <- c("VIM"=1, "PECAM1"=1, "THY1"=1,"NR2F1"=1, "DLK1"=1, "HIVEP3"=1, "CD68"=1, "CD14"=1, "HLA-A"=1,'HLA-DPA1'=1,'HLA-DPB1'=1,'MKI67'=1)


# marker.gmat.df.cutoff <- marker.gmat.df[,1:3]
# for(gene in marker.genes){
#     #cat('gene ',gene)
#     feature.value <- marker.gmat.df[,gene]
#     #cutoff <- quantile(feature.value,probs = c(0,1) )
#     cutoff <- quantile(feature.value,probs = c(low.q[gene],high.q[gene])) 
#     cutoff.low <- cutoff[1]
#     cutoff.high <- cutoff[2]
#     feature.value[feature.value < cutoff.low] <- cutoff.low
#     feature.value[feature.value > cutoff.high] <- cutoff.high
#     marker.gmat.df.cutoff <- cbind(marker.gmat.df.cutoff,feature.value)
# }
# colnames(marker.gmat.df.cutoff)[4:ncol(marker.gmat.df.cutoff)] <- marker.genes


# ######
# #options(repr.plot.height=5.8,repr.plot.width=5)
# res.marker <- list()
# for(gene in marker.genes){
#     p <- ggplot(marker.gmat.df.cutoff,aes_string(x="UMAP_1",y="UMAP_2",col= paste0("`",gene,"`") )) +
#     #p <- ggplot(marker.gmat.df,aes_string(x="UMAP_1",y="UMAP_2",col= paste0("`",gene,"`") )) +
#       geom_point(size = .2,show.legend = TRUE,alpha= .5 ) +
#       #scale_colour_manual(values = c('red','navy'))  +
#       #cale_colour_gradientn(colors = viridis(6,option = 'C'))  +
#       scale_colour_gradientn(colors = color_gradient_my )  +
#       ##scale_colour_gradientn(colors = color_tfdev )  +
#       #theme_classic() +
#       #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
#       theme(
#             legend.position = 'none',
#             #legend.position = 'top',
#             axis.text=element_blank(), 
#             axis.title = element_text(size = 15, face = "bold"),
#             axis.ticks = element_blank(),
#             panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             panel.background = element_rect(color="black", fill = NA,size=1),
#             plot.title = element_text(size = 15, face = "bold"),
#             #complete = TRUE
#             plot.margin = unit(c(1,1,1,1), "lines")
#            )+
#       ggtitle(gene) +
#       #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
#       #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
#       #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
#       #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
#       labs(x = "UMAP1", y = "UMAP2")
#       #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
#     #ggsave(filename = "PLA-term-RNA-depth.pdf",height=5.5,width=5)
#     res.marker[[gene]] <- p
#     #print(p)
# }

# ##arrange plot by patchwork
# options(repr.plot.height=15,repr.plot.width=15)
# res.marker[['DNMT1']] + res.marker[['CDH1']] + res.marker[['ERVFRD-1']] + res.marker[['FLT1']] +
# res.marker[['CSHL1']] + res.marker[['PSG8']] + res.marker[['VIM']] + res.marker[['PECAM1']] +
# res.marker[['CD14']] #+ res.marker[['MKI67']] + res.marker[['PLAC8']] + res.marker[['CD14']] +
# plot_layout(ncol=3,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ggsave(filename = 'marker.gene.umap.pdf',height=13.5,width=18,useDingbats=FALSE) #can save grid multiple plots

# ##arrange plot by patchwork
# options(repr.plot.height=15,repr.plot.width=15)
# res.marker[['FLT1']] + res.marker[['LEP']] + res.marker[['INSIG2']] + res.marker[['CSHL1']] +
# res.marker[['CSH1']] + res.marker[['PAPPA']] + res.marker[['PSG1']] + res.marker[['CGA']] +
# res.marker[['GCM1']] 
# plot_layout(ncol=3,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ggsave(filename = 'marker.gene.hormone.bias.umap.pdf',height=13.5,width=16,useDingbats=FALSE) #can save grid multiple plots


# ##arrange plot by patchwork
# options(repr.plot.height=13.5,repr.plot.width=18)
# res.marker[[1]] + res.marker[[2]] + res.marker[[3]] + res.marker[[4]] +
# res.marker[[5]] + res.marker[[6]] + res.marker[[7]] + res.marker[[8]] +
# res.marker[[9]] + res.marker[[10]] + res.marker[[11]] + res.marker[[12]] +
# plot_layout(ncol=4,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ggsave(filename = 'marker.gene.ageing.umap.pdf',height=13.5,width=18,useDingbats=FALSE) #can save grid multiple plots


# ##arrange plot by patchwork
# options(repr.plot.height=13.5,repr.plot.width=18)
# res.marker[['VIM']] + res.marker[['PECAM1']] + res.marker[['THY1']] + res.marker[['NR2F1']] +
# res.marker[['DLK1']] + res.marker[['HIVEP3']] + res.marker[['CD68']] + res.marker[['CD14']] +
# res.marker[['HLA-A']] + res.marker[['HLA-DPA1']] + res.marker[['HLA-DPB1']] + res.marker[['MKI67']] +
# plot_layout(ncol=4,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ggsave(filename = 'marker.gene.STR.umap.pdf',height=13.5,width=18,useDingbats=FALSE) #can save grid multiple plots





## Heretical clustering ###
# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.after.sp@cluster)), x.after.sp@cluster), function(x){
	SnapATAC::colMeans(x.after.sp[x,], mat="bmat");
	})
# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");

pdf( "pdfs/hclust.cluster.pdf",height=5.5,width=4,useDingbats = FALSE)
options(repr.plot.height=5.5,repr.plot.width=4)
plot(hc, hang=-1, xlab="");
dev.off()


####bin ave coverage(depth) of cluster (bmat)
ensemble.ls.df <- do.call(cbind,ensemble.ls) #620094
ensemble.ls.df <- ensemble.ls.df[rowSums(ensemble.ls.df) != 0,] #575002
ensemble.ls.df

options(repr.plot.height=5,repr.plot.width=5)
boxplot(ensemble.ls.df,outline = FALSE,col = color_snap_mod1,las=2)


##bin sum coverage(depth) of cluster (bmat)
ensembleSum.ls = lapply(split(seq(length(x.after.sp@cluster)), x.after.sp@cluster), function(x){
	SnapATAC::colSums(x.after.sp[x,], mat="bmat");
	})
ensembleSum.ls.df <- do.call(cbind,ensembleSum.ls) #620094

ensembleSum.ls.df <- ensembleSum.ls.df[rowSums(ensembleSum.ls.df) != 0,] #574062
ensembleSum.ls.df

options(repr.plot.height=5,repr.plot.width=5)
boxplot(ensembleSum.ls.df,outline = FALSE,col = color_snap_mod1,las=2) #quick

##bin sum coverage(depth) of cluster (bmat) by donor
#all.equal(cluster.df.add[,-c(1,2,3)],x.after.sp@metaData) #TRUE
ensembleSum.ls = lapply(split(seq(length(x.after.sp@cluster)),  #group by donor and cluster
                              list(cluster=x.after.sp@cluster,
                                   lib=x.after.sp@metaData$library)
                             ), 
                        function(x){SnapATAC::colSums(x.after.sp[x,], mat="bmat")}
                       )
ensembleSum.ls.df <- do.call(cbind,ensembleSum.ls) #620094
saveRDS(ensembleSum.ls.df,'ensembleSum.ls.df.rds')


# ensembleSum.ls.df.d1 <- ensembleSum.ls.df[,grep(pattern = 'D1',x=colnames(ensembleSum.ls.df))]
# ensembleSum.ls.df.d2 <- ensembleSum.ls.df[,grep(pattern = 'D2',x=colnames(ensembleSum.ls.df))]



# ensembleSum.ls.df.d1 <- ensembleSum.ls.df.d1[rowSums(ensembleSum.ls.df.d1) != 0,] #572203
# ensembleSum.ls.df.d1
# ensembleSum.ls.df.d2 <- ensembleSum.ls.df.d2[rowSums(ensembleSum.ls.df.d2) != 0,] #568099
# ensembleSum.ls.df.d2


pdf( "pdfs/depth_by_donor.pdf",height=5,width=10,useDingbats = FALSE)
par(mfrow = c(1,2))
options(repr.plot.height=5,repr.plot.width=10)
boxplot(ensembleSum.ls.df.d1,outline = FALSE,col = color_snap_mod1,las=2) #quick
boxplot(ensembleSum.ls.df.d2,outline = FALSE,col = color_snap_mod1,las=2) #quick
dev.off()





######stat sample and cluster cell number correlation######
res.stat <- table(cluster.df.add$cluster,cluster.df.add$library) 

res.stat <- t(as.matrix(res.stat))
res.stat.perc <- 100*res.stat/rowSums(res.stat)

rowSums(res.stat.perc)
#100

rownames(res.stat.perc) <- gsub('placenta_10X_late','Sample',rownames(res.stat.perc))

res.stat.perc.df <- as.data.frame(res.stat.perc) #wide to long
colnames(res.stat.perc.df) <- c('sample','cluster','percentage')



##stacked barplot with ggplot
options(repr.plot.width = 6.5, repr.plot.height = 3.5)
ggplot(res.stat.perc.df, aes( fill = cluster, x = sample, y = percentage )  )+
  geom_hline(yintercept = c(0,25,50,75,100),linetype='solid',size = .3, col = 'black') +
  geom_hline(yintercept = c(0), linetype = 'solid', size = 1, col = 'black') +
  geom_bar(position=position_stack(reverse=TRUE),stat='identity', width = 0.5  ) +
  coord_flip()+
  scale_fill_manual(values = color_good, name = 'cluster',labels=paste0(1:8) ) +
  theme_ipsum(base_family = 'sans') +
  ylab('Percentage of cluster (%)')+
  xlab('Samples of late pregnancy')+
  ggtitle('Reproducibility of clusters')

ggsave('pdfs/qc/cluster_reproducible.pdf',width = 6.5, height = 3.5)



# placenta_10X_late1 placenta_10X_late2 placenta_10X_late3
#   1                 800                495                600
#   2                 667                570                451
#   3                 647                504                515
#   4                 611                471                519
#   5                 559                389                653
#   6                 576                410                464
#   7                 719                402                457
#   8                 541                349                551
#   9                 301                222                272
#   10                334                306                221
#   11                176                168                192
#   12                288                 81                105
#   13                169                 17                 19
#   14                 50                 21                 31
#   15                 30                 26                 15
#   16                 40                 12                 14
    
#      placenta_10X_late5 placenta_10X_late6 placenta_10X_late9
#   1                 472                616                656
#   2                 450                770                687
#   3                 465                631                591
#   4                 490                583                585
#   5                 446                554                538
#   6                 373                636                503
#   7                 349                423                441
#   8                 262                454                371
#   9                 379                242                339
#   10                160                268                237
#   11                161                201                208
#   12                 48                176                261
#   13                 11                146                 85
#   14                 15                 24                 93
#   15                 11                 30                101
#   16                 13                 31                 68

#        D1   D2
#   1  1068  733
#   2   925  824
#   3   992  547
#   4   659  385
#   5   597  451
#   6   321  292
#   7   243  198
#   8   181  165
#   9   295   80
#   10  170   24
#   11   81   46
#   12   36   23





#total sample 
#6249 8107 
res.cor <- cor(res.stat) #0.977 #0.983 #0.988
res.stat[,1] = -1 * res.stat[,1]

res.stat.df <- as.data.frame(res.stat)
#res.stat.df <- reshape2::melt(res.stat)
colnames(res.stat.df) <- c('cluster','sample','count')

##plot horizonal barplot  #  '#94C6DD', '#1273AE'  vs '#F3A585','#C80927'
options(repr.plot.height=5.5,repr.plot.width=4)
ggplot(res.stat.df, aes(fill=sample, y=count, x=cluster )) + 
    #geom_hline(yintercept = c(0,25,50,75),linetype='solid',size=.3,col='black') +
    #geom_hline(yintercept = c(0),linetype='solid',size=1,col='black') +
    #geom_bar(position="stack", stat="identity",alpha=1,width = 0.5) +  
    geom_bar(position=position_stack(reverse=TRUE), stat="identity",width = 0.5) +  
    #xlim(100,0) +
    #scale_x_continuous(breaks = seq(100,0,-25),labels=seq(100,0,-25)) +
    #scale_y_reverse() +
    #annotate('text',x = 16, y = 1000,label = paste('spearman cor=',cor(res.stat)[1,2]) ) +
    coord_flip() +
    #scale_fill_viridis(discrete = T,option = "E") +
    scale_fill_manual(values = c('#94C6DD', '#1273AE'),labels=c('D1','D2'),name='sample' ) +
    #ggtitle("cell number count") +
    labs(title = "cell number count", subtitle=paste('spearman cor=',round(res.cor[1,2],digits=2)) ) +
    theme_ipsum(base_family = 'sans') + #to avoid Arial narrow problem in ggsave
    ylab("count")

ggsave(filename="donor1.vs.donor2.cluster.cor.pdf",width = 4, height = 5.5)







###################compare with signac clusters, mapping and plot tss score from signac object########
metadata.signac <- readRDS('../02.signac_harmony/cluster.df.add.rds')

#modify the rownames 
maptab <- c('D1' = 'placenta_donor1','D2'='placenta_donor2')
rowid <- paste0( maptab[metadata.signac$sample],"#",rownames(metadata.signac) )
rowid <- gsub(pattern = '-2$',replacement = '-1',x = rowid)
rownames(metadata.signac) <- rowid

#rowid: 11391
#cluster.df.add: 14356

shareid <- intersect(rownames(cluster.df.add),rowid ) #10846 #10481

# table (rowid %in% rownames(cluster.df.add) )
# FALSE  TRUE 
#   910 10481
# table(rownames(cluster.df.add) %in%  rowid  )
# FALSE  TRUE 
#  3875 10481

cluster.df.add.addtss <- cbind(cluster.df.add[shareid,],metadata.signac[shareid,c('nucleosome_signal','nucleosome_percentile','TSS.enrichment','TSS.percentile','high.tss','nucleosome_group','cluster','UMAP_1','UMAP_2')])
#10481
colnames(cluster.df.add.addtss)[c(ncol(cluster.df.add.addtss)-2,ncol(cluster.df.add.addtss)-1,ncol(cluster.df.add.addtss)-0)] <- c('cluster_sg','UMAP_1_sg','UMAP_2_sg')


##by TSS score (add from signac pipeline)
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add.addtss,aes(x=UMAP_1,y=UMAP_2,col= TSS.enrichment )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = rev(color_cellranger))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'D'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "TSS enrichment",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "PLA-term-ATAC-tssScore.pdf",height=5.5,width=5)


##plot by Signac clusters
#label right
options(repr.plot.height=5,repr.plot.width=5.5)
ggplot(cluster.df.add.addtss,aes(x=UMAP_1_sg,y=UMAP_2_sg,col=cluster  )) +
#ggplot(cluster.df.add.addtss,aes(x=UMAP_1,y=UMAP_2,col=cluster_sg  )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = c(color_snap_mod1) ) +
  #scale_colour_manual(values = color_good)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  
  theme(
        legend.position = 'right',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  #ggtitle("SnapATAC umap with Signac clusters" ) +
  ggtitle("Signac umap with SnapATAC clusters" ) +
  ##ggtitle("Signac UMAP" ) +
#   geom_text(data = centers_shift, #the halo
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "white", 
#             size = 6.5) +
#   geom_text(data = centers, 
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "black", 
#             size = 6) +
  guides(col = guide_legend(override.aes = list(size = 3))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

#ggsave(filename = "PLA-term-ATAC-UMAP-with-signac.pdf",height=5,width=5.5,useDingbats=FALSE)




# ###============fix final cluster result before call peaks (must redraw above plots)====================###
# ##---filter(doublet and low depth cells) x.after.sp with kept cells before fix SnapATAC object----#


# # ###(1) detect and filter doublets with scrublet (py), use snATACutils method##
# # ##see standalone script snapATAC.rmDoublets.R
# # saveRDS(object = x.after.sp@bmat,"bmat.rds")
# # res.doublet <- readRDS('res.doublet.scrublet.rds')

# # ##plot doublet with umap

# # nrow(cluster.df.add) == length(res.doublet$doublet_scores)
# # res.doublet.df <- data.frame(cluster.df.add[,1:3],doublet=res.doublet$doublet_scores)

# # ###by doublet score
# # options(repr.plot.height=5.5,repr.plot.width=5)
# # #ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= passed_filters)) +
# # ggplot(res.doublet.df[!flag.doublet,],aes(x=UMAP_1,y=UMAP_2,col= doublet )) +
# #   geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
# #   #scale_colour_manual(values = c('red','navy'))  +
# #   scale_colour_gradientn(colors = rev(color_cellranger))  +
# #   #scale_colour_gradientn(colors = viridis(6,option = 'D'))  + 
# #   #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
# #   #theme_classic() +
# #   #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
# #   theme(legend.position = 'top',
# #         axis.text=element_blank(), 
# #         axis.title = element_text(size = 15, face = "bold"),
# #         axis.ticks = element_blank(),
# #         panel.grid.major = element_blank(), 
# #         panel.grid.minor = element_blank(),
# #         panel.background = element_rect(color="black", fill = NA,size=1),
# #         plot.title = element_text(size = 15, face = "bold"),
# #         #complete = TRUE
# #         plot.margin = unit(c(1,1,1,1), "lines")
# #        )+
# #   ggtitle(paste(sample, "doublet",  sep=" ") ) +
# #   #ggtitle(paste(sample, "Sequence depth ",  sep=" ") ) +
# #   #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
# #   #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
# #   #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
# #   #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
# #   labs(x = "UMAP1", y = "UMAP2")
# #   #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
# # ggsave(filename = "PLA-term-ATAC-doublet_usebmat.pdf",height=5.5,width=5)


# # ##hist and boxplot of doublet distribution


# # res.bx.doublet <- boxplot(res.doublet$doublet_scores);abline(h=res.bx.doublet$stats[4],lty=2,lwd=1,col='red')
# # hist(res.doublet$doublet_scores);abline(v=res.bx.doublet$stats[4],lty=2,lwd=1,col='red')
# # #use threshold = 0.188

# # flag.doublet <- ifelse(res.doublet$doublet_scores >0.2,TRUE,FALSE)
# # table(flag.doublet)

# # flag.doublet
# # FALSE  TRUE 
# # 11553  2803

# # #flag.doublet
# # #FALSE  TRUE 
# # #10755  3601 

# # # flag.doublet
# # # FALSE  TRUE 
# # # 13945   411 

# # ##simple mark the predicted doublet
# # options(repr.plot.height=15,repr.plot.width=15)
# # plot(res.doublet.df[,2:3],pch = 16, type='p',col='grey',cex=0.5,xlab='UMAP1',ylab='UMAP2',main='predicted doublet (>0.2)',cex.main = 2.25 ) 
# # points(res.doublet.df[flag.doublet,2:3],pch = 16, cex=0.5,col='red')
# # table(flag.doublet)
# # #flag.doublet
# # FALSE  TRUE 
# # 11553  2803

# ##(2)filter high depth?

# flag.hidepth <- ifelse(cluster.df.add$logUMI >4.698,TRUE,FALSE) #50,000
# table(flag.hidepth)

# #flag.hidepth
# FALSE  TRUE 
# 13492   864 

# # #flag.hidepth
# # FALSE  TRUE 
# # 13503   853 
# # #flag.hidepth
# # FALSE  TRUE 
# # 12424  1932 

# ##(3)filter c8  and c10 in repeat 7?

# # flag.c8 <- cluster.df.add$cluster == '8'
# # table(flag.c8)
# # #flag.c8
# # FALSE  TRUE 
# # 13672   684 

# flag.c6.c8 <- (cluster.df.add$cluster == '6') | (cluster.df.add$cluster == '8')
# table(flag.c6.c8)
# FALSE  TRUE 
#  9675  1276 

# flag.c6 <- (cluster.df.add$cluster == '6') 
# #possible doublet 
# table(flag.c6)
# FALSE  TRUE 
# 10214   737 

# ##(4) filter dots that too disperse from cluster centroid (c8, or all?), binomial model fitting of d distribusion ?

# # flag.c8 <- cluster.df.add$cluster == '8'
# # cluster.sel = cluster.df.add[ flag.c8,]
# # #n_sel = nrow(cluster.sel)
# # #cat ('select for cluster ',id,' n = ',n_sel," \n")
# # #color = 'red'
# # #color = ifelse('cluster_sg' == id,'red','navy')
# # #     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
# # #     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
# # #     points(center[center$cluster == id,2:3],pch = 16, cex=1,col='black')

# # dx <- cluster.sel[,2] - as.numeric(center[center$cluster == '8',2])
# # dy <- cluster.sel[,3] - as.numeric(center[center$cluster == '8',3])
# # d <- sqrt(dx**2 + dy**2)

# # options(repr.plot.height=5,repr.plot.width=5)
# # hist(d,breaks=100,main=paste("cluster ",'8',sep=''),cex.main=2);abline(v=3.8)

# # d.filter <- 3.8
# # flag.dist <- d  >= 3.8

# # flag.dist.c8 <- rownames(cluster.df.add)  %in% rownames(cluster.sel[flag.dist,])
# # table(flag.dist.c8)
# # #flag.dist
# # FALSE  TRUE 
# # 14144   212 


# dotDist = function (cluster = NULL, id = NULL,center = centers){
#     colids = colnames(cluster) #must 'cluster','dim1','dim2' and rowname
#     cluster.sel = cluster[ cluster[,1] == id,]
#     n_sel = nrow(cluster.sel)
#     #cat ('select for cluster ',id,' n = ',n_sel," \n")
#     #color = 'red'
#     #color = ifelse('cluster_sg' == id,'red','navy')
# #     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
# #     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
# #     points(center[center$cluster == id,2:3],pch = 16, cex=1,col='black')
       
#     dx <- cluster.sel[,2] - as.numeric(center[center$cluster == id,2])
#     dy <- cluster.sel[,3] - as.numeric(center[center$cluster == id,3])
#     d <- sqrt(dx**2 + dy**2)
    
#     #options(repr.plot.height=15,repr.plot.width=15)
#     hist(d,breaks=100,main=paste("cluster ",id,sep=''),cex.main=2)
#     d.quantile <- quantile(d,prob = seq(0,1,0.1))
#     abline(v=d.quantile,lty=2,lwd=1,col='red')
#     #flag <- d <= as.numeric(d.quantile['80%'])
#     #flag <- d <=2
#     return(paste("cluster ",id," ok",sep='') )
# }


# dotClean = function (cluster = NULL, id = NULL,center = centers,d.filter = d.filter){
#     #return a flag.dist for each cluster
#     colids = colnames(cluster) #must 'cluster','dim1','dim2'
#     cluster.sel = cluster[ cluster[,1] == id,]
#     d.filter.sel = d.filter[id]
#     n_sel = nrow(cluster.sel)
#     #cat ('select for cluster ',id,' n = ',n_sel," \n")
#     color = 'red'
#     cat(paste("cluster ",id," filtering dist",sep=''),'\n')
#     #color = ifelse('cluster_sg' == id,'red','navy')
# #     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
# #     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
# #     points(center[center$cluster == id,2:3],pch = 16, cex=1,col='black')
       
#     dx <- cluster.sel[,2] - as.numeric(center[center$cluster == id,2])
#     dy <- cluster.sel[,3] - as.numeric(center[center$cluster == id,3])
#     d <- sqrt(dx**2 + dy**2)
    
#     flag.dist <- d > d.filter.sel
#     flag.dist.cl <- rownames(cluster)  %in% rownames(cluster.sel[flag.dist,])
    
#     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
#     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col='red')
#     points(cluster.sel[!flag.dist,2],cluster.sel[!flag.dist,3],pch = 16, cex=0.5,col='blue') 
#     points(center[center$cluster == id,2:3],pch = 16, cex=2.5,col='black')
    
#     #return(paste("cluster ",id," ok",sep='') )
#     return(flag.dist.cl)
    
# }

# par(mfrow=c(3,3))
# options(repr.plot.height=15,repr.plot.width=15)
# for(i in c('1','2','3','4','5','6','7','8','9') ){
#   dotDist(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers)
# }



# # par(mfrow=c(2,3))
# # options(repr.plot.height=10,repr.plot.width=15)
# # for(i in c('12','9','13','15','11','14') ){
# #   dotDist(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers)
# # }

# d.filter <- c( '1'=2.5, '2'=2.5,'3'=2,'4'=2,'5'=2,'6'=2,
#               '7'=1.5,'8'=2,'9'=1.5,'10'=0,'11'=0,'12'=0,
#               '13'=0 )


# par(mfrow=c(3,3))
# res.flag.dist <- list()
# options(repr.plot.height=15,repr.plot.width=15)
# for(i in c(1:9) ){
#   res.flag.dist[[i]] <- dotClean(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i,
#            center = centers,d.filter = d.filter)
# }


# flag.dist.combine <- apply(do.call(cbind,res.flag.dist),1,any)
# table(flag.dist.combine)

# FALSE  TRUE 
#  9991   960 

# #flag.dist.combine
# # FALSE  TRUE 
# # 13722   634 


# ######



# # ##(5) overall contour filter by kde python, output flat_contour###
# # write.table(cluster.df.add[,1:3],file='cluster_input_for_kde.txt',sep='\t',quote = FALSE, row.names = TRUE, col.names = TRUE) (unwork)




# #########quick and simple plot cluster umap depth ... after filtering##

# cluster.df.add.bk <- cluster.df.add 
# #saveRDS(cluster.df.add.bk,"cluster.df.add.bk.rds")
# ##cluster.df.add <- cluster.df.add.bk #recover

# #length(flag.doublet) == length(flag.hidepth) #TRUE
# #length(flag.c8) == length(flag.r7.c10) #TRUE
# #length(flag.doublet) == length(flag.r7.c10) #TRUE
# #length(flag.dist.c8) == length(flag.r7.c10) #TRUE

# #length(flag.hidepth) == length(flag.c6) #TRUE
# #length(flag.dist.combine) == length(flag.r7.c10) #TRUE
# length(flag.dist.combine) == length(flag.c6) #TRUE

# #flag.filter <- flag.doublet | flag.hidepth | flag.c8 | flag.r7.c10
# #flag.filter <-  flag.c8 | flag.r7.c10 #1254
# #flag.filter <- flag.doublet | flag.c8 | flag.r7.c10 #3804
# #flag.filter <- flag.hidepth | flag.c8 | flag.r7.c10 #2961
# #flag.filter <- flag.hidepth | flag.dist.c8 | flag.r7.c10 
# #flag.filter <- flag.hidepth | flag.dist.combine | flag.r7.c10 
# flag.filter <- flag.dist.combine | flag.c6
# table(flag.filter)

# flag.filter
# FALSE  TRUE 
#  9336  1615 

# # #flag.filter
# # FALSE  TRUE 
# # 12526  1830 

# #flag.filter
# #FALSE  TRUE 
# #12873  1483 

# # flag.filter
# # FALSE  TRUE 
# # 12401  1955 




# ##quick look for filtered 
# options(repr.plot.height=15,repr.plot.width=15)
# plot(cluster.df.add[,2:3],
#      pch = 16, type='p',col='grey',
#      cex=0.5,xlab='UMAP1',ylab='UMAP2',main='filter ...',cex.main = 2.25 ) 
# #points(cluster.df.add[flag.c8,2:3],pch = 16, cex=0.5,col='red')
# #points(cluster.df.add[flag.r7.c10,2:3],pch = 16, cex=0.5,col='red')
# #points(cluster.df.add[flag.hidepth,2:3],pch = 16, cex=0.5,col='red')
# #points(cluster.df.add[flag.doublet,2:3],pch = 16, cex=0.5,col='red')
# points(cluster.df.add[flag.filter,2:3],pch = 16, cex=0.5,col='red')

# options(repr.plot.height=5.5,repr.plot.width=5)
# plot(cluster.df.add[!flag.filter,2:3],
#      pch = 16, type='p',col='grey',
#      cex=.2,xlab='UMAP1',ylab='UMAP2',main='filter ...',cex.main = 2.25 ) 



#######read in cluster_df_add after snapatac2 umap tuning renaming, cleaning for snapATAC filling###########

##full cell type snapatac2
cluster.df.add.snapatac2 <- read.table('cluster_df_add_final.snapatac2.txt',sep='\t',stringsAsFactors = FALSE,header=TRUE,row.names = 1,comment.char = "")
#26196 × 13, rm c8 c9(doublet)
colnames(cluster.df.add.snapatac2)[1:2] <- c('UMAP_1','UMAP_2')
cluster.df.add.snapatac2$cluster <- factor( as.character(cluster.df.add.snapatac2$cluster), levels = as.character(sort(unique(cluster.df.add.snapatac2$cluster))) )

##cstb snapatac2
cluster.df.add.cstb.snapatac2 <- read.table('cluster_df_add_cstb.snapatac2.txt',sep='\t',stringsAsFactors = FALSE,header=TRUE,row.names = 1,comment.char = "")
#24692 × 15 cstb only
colnames(cluster.df.add.cstb.snapatac2)[1:2] <- c('UMAP_1','UMAP_2')
cluster.df.add.cstb.snapatac2$cluster <- factor( as.character(cluster.df.add.cstb.snapatac2$cluster), levels = as.character(sort(unique(cluster.df.add.cstb.snapatac2$cluster))) )

cluster.df.add.cstb.snapatac2$cluster_merge <- factor( as.character(cluster.df.add.cstb.snapatac2$cluster_merge), levels = as.character(sort(unique(cluster.df.add.cstb.snapatac2$cluster_merge))) )


cluster.df.add <- read.table('cluster.df.add.txt',sep='\t',stringsAsFactors = FALSE,header=TRUE,row.names = 1,comment.char = "")
#31684 × 32
cluster.df.add$cluster <- factor( as.character(cluster.df.add$cluster), levels = as.character(sort(unique(cluster.df.add$cluster))) )




table(rownames(cluster.df.add.snapatac2) %in% rownames(cluster.df.add) )
TRUE 
26196
table(rownames(cluster.df.add.cstb.snapatac2) %in% rownames(cluster.df.add) )
TRUE 
24692

all.equal(rownames(x.after.sp@metaData),rownames(cluster.df.add)) #TRUE
31684 x 545479



source('quickDimPlot_labelon.r')
quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', color_use = color_good, title= 'late combined-ATAC', shrink.x = 0.5, shrink.y = 0.1,shuffle = TRUE,pt.size = .2)

quickDimPlot_labelon(data = cluster.df.add.snapatac2, feature = 'cluster', color_use = color_good, title= 'late combined-ATAC', shrink.x = .3, shrink.y = .1,shuffle = TRUE,pt.size = .2)
ggsave('pdfs/PLA-late-combine.final.pdf',height=7.5,width=8.5,useDingbats=FALSE)

quickDimPlot_labelon(data = cluster.df.add.cstb.snapatac2, feature = 'cluster', color_use = color_good, title= 'late combined-ATAC', shrink.x = 0.3, shrink.y = 0.1,shuffle = FALSE,pt.size = .5)
ggsave('pdfs/PLA-late-combine.cstb.pdf',height=7.5,width=8.5,useDingbats=FALSE)

quickDimPlot_labelon(data = cluster.df.add.cstb.snapatac2, feature = 'cluster_merge', color_use = color_good, title= 'late combined-ATAC', shrink.x =0.3, shrink.y = .1,shuffle = FALSE,pt.size = .3)
ggsave('pdfs/PLA-late-combine.cstb.merge.pdf',height=7.5,width=8.5,useDingbats=FALSE)




#############merge c4 and c6 as intermediate cluster?###########################3
#https://stackoverflow.com/questions/11737193/replace-given-value-in-vector
# map_cluster_list <- list( 
#     '1' = '1',
#     '2' = '2',
#     '3' = '3',
#     '4' = '4',
#     '5' = '5',
#     '6' = '4',
#     '7' = '6',
#     '8' = '7',
#     '9' = '8',
#     '10' = '9'
# )

table(cluster.df.add$cluster)
1    2    3    4    5    6    7    8    9   10 
3457 3583 3186 3096 2983 2814 2651 1449 1050  423

table(cluster.df.add$cluster_merge)
  1    2    3    4    5    6    7    8 
3457 3583 6282 1449 2983 5465 1050  423

quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', color_use = color_good, title= 'late combined-ATAC', shrink.x = .3, shrink.y = .1,shuffle = FALSE,pt.size = .3)

quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster_merge', color_use = color_good, title= 'late combined-ATAC', shrink.x = .3, shrink.y = .1,shuffle = FALSE,pt.size = .3)



cid <- as.character(cluster.df.add$cluster)

cid_mod <- plyr::revalue(cid,   
              c('1' = '1',
                '2' = '2',
                '3' = '3',
                '4' = '4',
                '5' = '5',
                '6' = '4',
                '7' = '6',
                '8' = '7',
                '9' = '8',
                '10' = '9'
                )
             
)

cid.df <- data.frame(cid=cid,cid_mod = cid_mod, stringsAsFactors = FALSE)

# cid.df <- data.frame(ori=cid,change = cid, stringsAsFactors = FALSE)
# cid.df %>% dplyr::case_when(
#     ori == '1' ~ '1',
#     ori == '2' ~ '2',
#     ori == '3' ~ '3',
#     ori == '4' ~ '4',
#     ori == '5' ~ '5',
#     ori == '6' ~ '4',
#     ori == '7' ~ '6',
#     ori == '8' ~ '7',
#     ori == '9' ~ '8',
#     ori == '10' ~ '9'
# )

table(cid_mod)
 1    2    3    4    5    6    7    8    9 
3457 3583 3186 5910 2983 2651 1449 1050  423


cluster.df.add$cluster_merge_new <- factor(cid_mod,levels = sort(unique(cid_mod)) ) 

quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster_merge_new', color_use = color_good, title= 'late combined-ATAC', shrink.x = .3, shrink.y = .1,shuffle = FALSE,pt.size = .3)
ggsave('pdfs/PLA-late-combine.cstb.merge.new.pdf',height=7.5,width=8.5,useDingbats=FALSE)


cluster.df.add$cluster_ori <- cluster.df.add$cluster
cluster.df.add$cluster <- cluster.df.add$cluster_merge_new


all.equal( rownames(x.after.sp@metaData), rownames(cluster.df.add)  ) #TRUE
quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', color_use = color_good, title= 'late combined-ATAC', shrink.x = .3, shrink.y = .1,shuffle = FALSE,pt.size = .3)

x.after.sp@cluster <- cluster.df.add$cluster



#save again at the end of this script



####################use above code to redraw cluster, umap , depth, marker gene###############
# cluster.df.add <- cluster.df.add[!flag.filter,] #10951 -> 9336
# #cluster.df.add$cluster <- droplevels(cluster.df.add$cluster)


# ###filter STR and cluster 7

# cluster.df.add.filter <- subset(cluster.df.add, ! cluster %in% c('7','9','10','11','12') ) #8140







#####################the final filtering object###################
idx1 <- which (rownames(x.after.sp@metaData) %in% rownames(cluster.df.add.snapatac2)     )
idx2 <- which (rownames(x.after.sp@metaData) %in% rownames(cluster.df.add.cstb.snapatac2)     )

table(rownames(x.after.sp@metaData) %in% rownames(cluster.df.add.snapatac2))
FALSE  TRUE 
 5488 26196 

# FALSE  TRUE 
#  1615  9336 

# FALSE  TRUE 
#  1830 12526
#FALSE  TRUE 
#  498 14148 




#x.after.sp.bk <- x.after.sp
x.after.sp.final <- x.after.sp[idx1,]

number of barcodes: 26196
number of bins: 545479
number of genes: 0
number of peaks: 0
number of motifs: 0

x.after.sp.cstb <- x.after.sp[idx2,]

number of barcodes: 24692
number of bins: 545479
number of genes: 0
number of peaks: 0
number of motifs: 0



# number of barcodes: 9336
# number of bins: 620094
# number of genes: 58344
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 12526
# number of bins: 620094
# number of genes: 58344
# number of peaks: 0
# number of motifs: 0

# number of barcodes: 14148
# number of bins: 620094
# number of genes: 11
# number of peaks: 317817
# number of motifs: 0

##check rowname 
all.equal(rownames(x.after.sp.final@metaData),rownames(cluster.df.add.snapatac2) ) #TRUE
table(x.after.sp.final@metaData$passed_filters == cluster.df.add.snapatac2$n_fragment ) 
FALSE  TRUE 
  641 25555


all.equal(rownames(x.after.sp.cstb@metaData),rownames(cluster.df.add.cstb.snapatac2) ) #TRUE
table(x.after.sp.cstb@metaData$passed_filters == cluster.df.add.cstb.snapatac2$n_fragment ) 
FALSE  TRUE 
  588 24104



##replace snapatac2 tuned umap and cluster to snapATAC obj

x.after.sp.final@umap <- as.matrix(cluster.df.add.snapatac2[,c('UMAP_1','UMAP_2')])
plot(x.after.sp.final@umap)

x.after.sp.final@cluster <- cluster.df.add.snapatac2$cluster

table(x.after.sp.final@cluster)
1    2    3    4    5    6    7    8    9   10   11   12   13   14 
3457 3583 3186 3096 2983 2814 2651 1449 1050  423  909  222  204  169



x.after.sp.cstb@umap <- as.matrix(cluster.df.add.cstb.snapatac2[,c('UMAP_1','UMAP_2')])
plot(x.after.sp.cstb@umap)

x.after.sp.cstb@cluster <- cluster.df.add.cstb.snapatac2$cluster

table(x.after.sp.cstb@cluster)
 1    2    3    4    5    6    7    8    9   10 
3457 3583 3186 3096 2983 2814 2651 1449 1050  423 


# #drop levels
# x.after.sp@cluster <- droplevels(x.after.sp@cluster)
#  1    2    3    4    5    7    8    9   10   11   12   13 
# 1801 1749 1539 1044 1048  613  441  346  375  194  127   59 

# #rename cluster id

# mapid <- list('1'='1','2'='2', '3'='3','4'='4','5'='5',
# '7'='6',
# '8'='7',
# '9'='8',
# '10'='9',
# '11'='10',
# '12'='11',
# '13'='12'
#  )

# names(mapid)

# cluster <- as.character(x.after.sp@cluster)

# for(i in 1:length(cluster)){ cluster[i] = mapid[[ cluster[i] ]] }

# cluster <- factor(cluster,levels=gtools::mixedsort( names(table(cluster)) ))
# x.after.sp@cluster <- cluster

# table(x.after.sp@cluster)
#   1    2    3    4    5    6    7    8    9   10   11   12 
# 1801 1749 1539 1044 1048  613  441  346  375  194  127   59 


# ##add cell name and cell color


# map_cellname <- list(
#     '1'='STB1',
#     '2'='STB naive', 
#     '3'='STB3',
#     '4'='STB4',
#     '5'='STB2',
#     '6'='Syncytial knot',
#     '7'='unknow1',
#     '8'='STB-new',
#     '9'='STR1',
#     '10'='CTB',
#     '11'='STR2',
#     '12'='STR3'
#  )



# #STB ()  navy gradient
# #Syncitial knot: 
# #CTB 
# #STR greens (3)

# reds <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf')
# blues <- c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')

# map_cellcolor <- list(
#     '1'=reds[5],
#     '2'=reds[6], 
#     '3'=reds[3],
#     '4'=reds[2],
#     '5'=reds[4],
#     '6'=reds[1],
#     '7'='grey', #unknow1
#     '8'='black', #unknow2
#     '9'=blues[1],
#     '10'='darkgreen', #CTB
#     '11'=blues[2],
#     '12'=blues[3]
#  )


# all.equal(x.after.sp@cluster,cluster.df.add$cluster)#TRUE
# cellname <- as.character(x.after.sp@cluster)
# cellcolor<- as.character(x.after.sp@cluster)
# for(i in 1:length(cellname)){ cellname[i] = map_cellname[[ cellname[i] ]] }
# for(i in 1:length(cellcolor)){ cellcolor[i] = map_cellcolor[[ cellcolor[i] ]] }

# cluster.df.add[,'cellname'] <- factor(cellname)
# cluster.df.add[,'cellcolor'] <- factor(cellcolor)



###save object final
#saveRDS(x.after.sp.final,"x.after.sp.final.rds")     
x.after.sp.final <- readRDS("x.after.sp.final.rds")
number of barcodes: 26196
number of bins: 545479
number of genes: 0
number of peaks: 0
number of motifs: 0


#saveRDS(x.after.sp.cstb,"x.after.sp.cstb.rds")     
x.after.sp.cstb <- readRDS("x.after.sp.cstb.rds")


# write.table(cluster.df,file='snapATAC.PLA-term-ATAC.umap.bin5k.final.txt',sep='\t',quote = FALSE, row.names = TRUE, col.names = TRUE)
# saveRDS(cluster.df.add,"cluster.df.add.final.rds")
# cluster.df.add <- readRDS("cluster.df.add.final.rds")

x.after.sp <- x.after.sp.cstb
#x.after.sp <- x.after.sp.final


rm(x.after.sp.cstb)
rm(x.after.sp.final)

gc()


##plot use above plot code again



###################






#################post-clustering steps start from step 13###########################
##########step 13 Identify peaks with macs2########

system('which snaptools',intern = TRUE)
#/home/mjwang/anaconda3/envs/myenv_ori/bin/snaptools

#'/home/mjwang/anaconda3/envs/myenv_ori/bin/snaptools'
#'/home/mjwang/.conda/envs/myenv/bin/snaptools'

system('which macs2',intern = TRUE)
#/home/mjwang/anaconda3/envs/myenv_ori/bin/macs2

#/home/mjwang/anaconda3/envs/myenv_ori/bin/macs2
#'/home/mjwang/.conda/envs/myenv/bin/macs2'
clusters.sel = names(table(x.after.sp@cluster))[which(table(x.after.sp@cluster) > 100)];
#'1''2''3''4''5''6''7''8''9'
#'1''2''3''4''5''7''8''9''10''11''12'




peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i]);
    peaks = runMACS(
        obj=x.after.sp[which(x.after.sp@cluster==clusters.sel[i]),], 
        output.prefix=paste0("MACS/placenta.cluster", gsub(" ", "_", clusters.sel)[i]),
        path.to.snaptools="/home/mjwang/anaconda3/envs/myenv_ori/bin/snaptools",
        path.to.macs="/home/mjwang/anaconda3/envs/myenv_ori/bin/macs2",
        gsize="hs", # mm, hs, etc
        buffer.size=500, 
        num.cores=1,
        macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
        tmp.folder=tempdir()
   );
	peaks
 	}, mc.cores=30);
    #}, mc.cores=15);


#mkdir MACS first
#dir.create('MACS')


peaks.names = system("ls MACS/* | grep narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
    peak.df = read.table(x)
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
  })

peak.gr = reduce(Reduce(c, peak.gr.ls)); #192048 cstb #173870 left #remove redundant, result in 299708

chrid <- as.character(seqnames(peak.gr))
#idx <- grep('^GL|^KI|^chrUn|random|unknown',chrid)
#chrid[idx]

idx <- grepl('^GL|^KI|^chrUn|random|unknown',chrid)
peak.gr.sel <- peak.gr[!idx,] #191816

peaks.df = as.data.frame(peak.gr.sel)[,1:3];
write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
		quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
		row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
		fileEncoding = "")


##merge peak (Granja's method?)


##stat peak 

hist(GenomicRanges::width(peak.gr),breaks = 100 ,xlim  = c(1,2000));abline(v=200)
#hist(GenomicRanges::width(x.sp@peak),breaks = 100 ,xlim  = c(1,2000));abline(v=200)




##use snaptools to snap-add-pmat for each snap files
# snaptools snap-del --snap-file placenta-term-atac-1.snap --session-name PM
# snaptools snap-add-pmat --snap-file placenta-term-atac-1.snap --peak-file peaks.combined.bed > log.snapaddpmat.donor1 2>&1 &

# snaptools snap-del --snap-file placenta-term-atac-new.snap --session-name PM
# snaptools snap-add-pmat --snap-file placenta-term-atac-1.snap --peak-file peaks.combined.bed > log.snapaddpmat.donor2 2>&1 &




######step 14 create a cell-by-peak matrix#####


#unique(gsub(pattern = "/mnt/disk1/mjwang/pwdex/placenta_10X_new/03.snapATAC/snap_from_10xbam/|/mnt/disk1/mjwang/pwdex/placenta_10X/02.snapATAC/PLA-term-ATAC-1/snap_from_10xbam/",replacement = "/home/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/",x = x.after.sp@file ) )
#x.after.sp@file = gsub(pattern = "/mnt/disk1/mjwang/pwdex/placenta_10X_new/03.snapATAC/snap_from_10xbam/|/mnt/disk1/mjwang/pwdex/placenta_10X/02.snapATAC/PLA-term-ATAC-1/snap_from_10xbam/",replacement = "/home/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/",x = x.after.sp@file )

x.after.sp = addPmatToSnap(x.after.sp,num.cores = 6) #12526 x 299708 #6820 x 272901 pmat #quick


x.after.sp
number of barcodes: 24692
number of bins: 545479
number of genes: 0
number of peaks: 191816
number of motifs: 0

# number of barcodes: 9336
# number of bins: 620094
# number of genes: 58344
# number of peaks: 173870
# number of motifs: 0

# number of barcodes: 12526
# number of bins: 620094
# number of genes: 58344
# number of peaks: 299708
# number of motifs: 0

#number of barcodes: 14646
# number of bins: 620094
# number of genes: 11
# number of peaks: 317817
# number of motifs: 0


#saveRDS(x.after.sp,"x.after.sp.onestep_filter.fixed.addpmat.rds") 
#x.after.sp <- readRDS("x.after.sp.addpmat.rds") 


max(x.after.sp@pmat) #255

pmat <- x.after.sp@pmat
dim(pmat)
24692 x 191816

all.equal(rownames(pmat),x.after.sp@metaData$barcode) #TRUE

rownames(pmat) <- rownames(x.after.sp@metaData)

colnames(pmat) <- mcols(x.after.sp@peak)$name

saveRDS(pmat,'pmat.rds')



###############step 15 identify differentially accessible peaks##########

################DARs for one cluster####################
#cluster.pos = 7
#cluster.pos = 8
cluster.pos = 9
DARs = findDAR( 
    obj=x.after.sp,
    input.mat="pmat",
    cluster.pos=cluster.pos,
    cluster.neg.method="knn",
    test.method="exactTest",
    bcv=0.1, #0.4 for human, 0.1 for mouse
    seed.use=10
  );
# test for each peaks

DARs$FDR = p.adjust(DARs$PValue, method="BH");
#idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
idy = which (DARs$PValue < 0.05 & DARs$logFC >0)

#bcv = 0.1, FDR < 0.05, logFC > 0
#cluster.pos=4 ,  1689 DAR peaks
#cluster.pos = 6,  23664 DAR peaks (bcv = 0.1, FDR < 0.05, logFC > 0)


options(repr.plot.width=15,repr.plot.height=15)
#par(mfrow = c(2, 1));
plot(DARs$logCPM, DARs$logFC, 
    pch=19, cex=0.1, col="grey", 
    ylab="logFC", xlab="logCPM",
    main=paste("cluster ",cluster.pos,sep="")
  );
points(DARs$logCPM[idy], 
    DARs$logFC[idy], 
    pch=19, 
    cex=0.5, 
    col="red"
  );
abline(h = 0, lwd=1, lty=2);
text(x = 10,y=5,labels = paste("n_DARs = ",length(idy),sep='') )

covs = Matrix::rowSums(x.after.sp@pmat);
vals = Matrix::rowSums(x.after.sp@pmat[,idy]) / covs; #percentage maximum? Greenleaf's method?
vals.zscore = (vals - mean(vals)) / sd(vals);
plotFeatureSingle(
    obj=x.after.sp,
    feature.value=vals.zscore,
    method="umap", 
    main=paste("Cluster ",cluster.pos,sep=""),
    point.size=0.3, 
    point.shape=19, 
    #down.sample=5000,
    quantiles=c(0.01, 0.99)
  );

mcols(x.after.sp@peak )['name'][idy,] 
#2835 for cluster 9


#6724 for cluster 6
#1592 for cluster 11

DARs.filterp <- DARs[DARs$PValue <= 0.01,] #653 of 272901
DARs.filterp <- DARs.filterp[order(DARs.filterp$PValue,decreasing = FALSE),]

####compare p.adjust####
#http://rcompanion.org/rcompanion/f_01.html

#https://www.researchgate.net/post/Can_anyone_explain_how_to_calculate_adjusted_p-values_q-values_following_Benjamini_Hochberg_correction
#The problem with Bonferroni, Hochberg and Holm is that they were developed for small n, e.g. n < 100. If your n is large, these methods are very conservative
DARs.filterp$FDR <- p.adjust(p = DARs.filterp$PValue,method = 'BH')
DARs.filterp$FDR.BH <- p.adjust(p = DARs.filterp$PValue,method = 'BH')
DARs.filterp$FDR.Bonferroni <- p.adjust(p = DARs.filterp$PValue,method = 'bonferroni')
DARs.filterp$FDR.Holm <- p.adjust(p = DARs.filterp$PValue,method = 'holm')
DARs.filterp$FDR.Hochberg <- p.adjust(p = DARs.filterp$PValue,method = 'hochberg')


Y=data.frame(
   X=DARs.filterp$PValue,
   BH=DARs.filterp$FDR.BH,  
   Bonferroni=DARs.filterp$FDR.Bonferroni,
   Holm=DARs.filterp$FDR.Holm,
   Hochberg=DARs.filterp$FDR.Hochberg 

)

Y.melt <- reshape2::melt(Y,id.vars = 'X')
ggplot(data=Y.melt,aes(x=X,y=value,group=variable,col=variable,shape=variable) ) +
  geom_point(size=5) +
  theme_bw()
#########


# matplot(X, Y,
#         xlab="Raw p-value",
#         ylab="Adjusted p-value",
#         type="l",
#         #asp=1,
#         col=1:4,
#         lty=1,
#         lwd=2
#        )


###############DARs for all clusters###################
#cluster.df.add$cluster <- droplevels(cluster.df.add$cluster)

all.equal(cluster.df.add$cluster,x.after.sp@cluster)#TRUE
#idy.ls = lapply(levels(cluster.df.add$cluster), function(cluster_i){ ##quick
idy.ls = lapply(levels(x.after.sp@cluster), function(cluster_i){ ##quick
    cat ('cluster ',cluster_i,sep='')
	DARs = findDAR(
		obj=x.after.sp,
		input.mat="pmat",
		cluster.pos=cluster_i,
		cluster.neg=NULL,
		cluster.neg.method="knn",
		bcv=0.1, #in fact, 0.4 for human, 0.1 for mouse
		test.method="exactTest",
		seed.use=10
		);
	DARs$FDR = p.adjust(DARs$PValue, method="BH");
	idy = which(DARs$PValue < 5e-2 & DARs$logFC > 0);
    #idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
	if((x=length(idy)) < 1L){
    #if((x=length(idy)) < 2000L){
			PValues = DARs$PValue;
			PValues[DARs$logFC < 0] = 1;
			idy = order(PValues, decreasing=FALSE)[1:2000];
			rm(PValues); # free memory
	}
	idy
  })
names(idy.ls) = levels(x.after.sp@cluster);
#names(idy.ls) = levels(cluster.df.add$cluster);

##Greenleaf's paper heatmap method?
covs = Matrix::rowSums(x.after.sp@pmat); #per cell
#saveRDS(covs,'covs.rds')

all.equal(covs,readRDS('covs.rds')) #TRUE

#pdf( "pdfs/PLA-term-ATAC-DARs.visualize2.pdf",height=15,width=15,useDingbats = FALSE)
options(repr.plot.width=15,repr.plot.height=15)
#par(mfrow = c(3, 3));
for(cluster_i in levels(x.after.sp@cluster) ){
#for(cluster_i in levels(cluster.df.add$cluster) ){
	print(cluster_i)
	idy = idy.ls[[cluster_i]];
	vals = Matrix::rowSums(x.after.sp@pmat[,idy]) / covs;
	vals.zscore = (vals - mean(vals)) / sd(vals);
    peaknames = as.data.frame(x.after.sp@peak[idy])
    write.table(peaknames,file = paste('DARs_exactTest/DARs.cluster',cluster_i,".exactTest.txt",sep=""),
               sep='\t',quote=FALSE,row.names=FALSE)
    
	plotFeatureSingle(
		obj=x.after.sp,
		feature.value=vals.zscore,
		method="umap", 
		main=paste("Cluster ",cluster_i,sep = ""),
		point.size=0.6, 
		point.shape=19 
		#down.sample=5000,
		#quantiles=c(0.01, 0.99)
		);
  }
#dev.off()

##the DAR size
temp <- sapply(idy.ls,length)
paste0(names(temp),':',temp,collapse = ',')
'1:1605,2:35969,3:4006,4:512,5:1603,6:1133,7:2143,8:5183,9:28677'

#1:808,2:24592,3:2163,4:621,5:761,6:1674,7:937,8:2140,9:22377,10:9344,11:6948,12:781

##certain cluster DARs
# saveRDS(DARs,"DARs_exactTest/DARs.fdr0.05.rds")
# saveRDS(DARs,"DARs_exactTest/DARs.pvalue0.05.rds")
# saveRDS(DARs,"DARs_exactTest/DARs.pvalue0.05.realout.rds") #not order 2000 then output if few peaks detected

##all DARs
saveRDS(idy.ls,"DARs_exactTest/idy.ls.rds")
idy.ls <- readRDS("DARs_exactTest/idy.ls.rds") 
#all.equal(idy.ls,idy.ls.1) #TRUE


#saveRDS(x.after.sp,"x.after.sp.addpmat_doDARs.rds") 
x.after.sp <- readRDS("x.after.sp.addpmat_doDARs.rds")

number of barcodes: 9336
number of bins: 620094
number of genes: 58344
number of peaks: 173870
number of motifs: 0
               
               
# ######step 16 Motif variable analysis (see more comprehensive analysis in dir chromVAR_TF_specific )#####

# library(chromVAR);
# library(motifmatchr);
# library(SummarizedExperiment);
# library("BSgenome.Hsapiens.UCSC.hg38")

# pmat <- x.after.sp@pmat #bk of pmat before binarize pmat
# #saveRDS(pmat,'pmat.raw.rds')
# pmat <- readRDS('pmat.raw.rds')

# #x.after.sp@pmat <- pmat

# max(x.after.sp@pmat) #255 #1

# x.after.sp = makeBinary(x.after.sp, "pmat");

# x.after.sp@mmat = runChromVAR( #cause problem
#     obj=x.after.sp,
#     input.mat="pmat",
#     genome=BSgenome.Hsapiens.UCSC.hg38,
#     min.count=10,
#     species="Homo sapiens"
#   );
# x.after.sp;

# #Error in loadFUN(x, seqname, ranges): trying to load regions beyond the boundaries of non-circular sequence "chrUn_KI270438v1"

# #"chrUn_KI270438v1" %in% names(seqlengths(BSgenome.Hsapiens.UCSC.hg38) ) #TRUE



# ##try step by step way instead

# data.use = x.after.sp@pmat
# peak.use = x.after.sp@peak #173870 #299708
# min.count = 10
# genome = BSgenome.Hsapiens.UCSC.hg38
# species = 'Homo sapiens'

# idy = which(Matrix::colSums(data.use) >= min.count) #170904 #294466
# data.use = data.use[, idy, dropping = TRUE]
# peak.use = peak.use[idy]

# #idy = which(  seqnames(peak.use) != 'chrUn_KI270438v1'  )
# #chrUn_KI270438v1
# #chrUn_KI270466v1

# idy = which(!grepl(pattern = "^chrUn",x = as.character(seqnames(peak.use)) ) )
# data.use = data.use[, idy, dropping = TRUE] #170852 #294404
# peak.use = peak.use[idy] 

# cat("Epoch: creating chromVAR object ... \n", file = stderr())
# rse <- SummarizedExperiment(assays = list(counts = t(data.use)), 
#     rowRanges = peak.use, colData = DataFrame(Cell_Type = 1:nrow(data.use), 
#         depth = Matrix::rowSums(data.use)))

# cat("Epoch: computing GC bias ... \n", file = stderr())
# rse <- addGCBias(rse, genome = genome)

# cat("Epoch: getting JASPAR motifs ... \n", file = stderr())
# motifs <- getJasparMotifs(collection = "CORE", species = species)

# cat("Epoch: scanning motifs in the peaks ... \n", file = stderr())
# motif_mm <- matchMotifs(motifs, rse, genome = genome)

# saveRDS(object = rse, 'rse.rds')
# saveRDS(object = motif_mm, 'motif_mm.rds')

# #all.equal(rse,readRDS('rse.rds') ) #TRUE
# #all.equal(motif_mm,readRDS('motif_mm.rds') ) #TRUE

# cat("Epoch: calculating motif variability between cells ... \n" )
# dev <- computeDeviations(object = rse, annotations = motif_mm) #very slow when load is high? #quick  in fact, if mem is enough
# dev_mat = t(assay(dev))

# saveRDS(file='dev_mat.SnapATAC.rds',object = dev_mat )
# cat("Epoch: Done ... \n", file = stderr())





#saveRDS(x.after.sp,"x.after.sp.stage7.rds")
#x.after.sp = readRDS("x.after.sp.stage7.rds")

#saveRDS(x.after.sp,"x.after.sp.nonorm_addpmat_doDARs_addmmat.rds") 

#####step 17  findMotifGenome homer search for known motifs#####

# DARs = findDAR(
#     obj=x.after.sp,
#     input.mat="pmat",
#     cluster.pos=4,
#     cluster.neg.method="knn",
#     test.method="exactTest",
#     bcv=0.4, #0.4 for human, 0.1 for mouse
#     seed.use=10
#   );
# DARs$FDR = p.adjust(DARs$PValue, method="BH");
# idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);

##homer findMotifGenome.pl
#/home/mjwang/.conda/envs/myenv/share/homer-4.10-0/bin/findMotifsGenome.pl ./homer/DoubleNegativeTcell/target_625245a041713 hg19 ./homer/DoubleNegativeTcell -len 10 -size 300 -S 2 -p 5 -cache 100 -fdr 5 -nomotif 

#add hg19 for homer 
#perl /home/mjwang/.conda/envs/myenv/share/homer-4.10-0/.//configureHomer.pl -install hg19

#add hg38 for homer
#perl /home/mjwang/.conda/envs/myenv/share/homer-4.10-0/.//configureHomer.pl -install hg38

#cluster.id = '4'

system("which findMotifsGenome.pl",intern = TRUE);
#/home/mjwang/anaconda3/envs/myenv_ori/bin/findMotifsGenome.pl
#/home/mjwang/.conda/envs/myenv/bin/findMotifsGenome.pl


cluster.id = '9'
idy = idy.ls[[cluster.id]] #use already identified DARs

result.dir = paste("./DARs_exactTest/homer/cluster",cluster.id,sep='')

if(dir.exists(result.dir)){ unlink(result.dir,recursive = TRUE)   } #need recursive = TRUE
               
motifs = runHomer(
	x.after.sp[,idy,"pmat"], 
	mat = "pmat",
	path.to.homer = "/home/mjwang/anaconda3/envs/myenv_ori/bin/findMotifsGenome.pl",
	result.dir = paste("./DARs_exactTest/homer/cluster",cluster.id,sep=''),
	num.cores=25,
	genome = 'hg38',
	motif.length = 10,
	scan.size = 300,
	optimize.count = 2,
	background = 'automatic',
	local.background = FALSE,
	only.known = FALSE,
    #only.known = TRUE,
	only.denovo = FALSE,
	fdr.num = 5,
	cache = 100,
	overwrite = TRUE,
	keep.minimal = FALSE
	);


##homer for all clusters (known motif only)
motifs.ls = lapply(levels(x.after.sp@cluster), function(cluster.id){
    
    idy = idy.ls[[cluster.id]] #use already identified DARs
    motifs = runHomer(
        x.after.sp[,idy,"pmat"], 
        mat = "pmat",
        path.to.homer = "/home/mjwang/anaconda3/envs/myenv_ori/bin/findMotifsGenome.pl",
        result.dir = paste("./DARs_exactTest/homer/cluster",cluster.id,sep=''),
        num.cores=25,
        genome = 'hg38',
        motif.length = 10,
        scan.size = 300,
        optimize.count = 2,
        background = 'automatic',
        local.background = FALSE,
        only.known = FALSE,
        #only.known = TRUE,
        only.denovo = FALSE,
        fdr.num = 20,
        cache = 100,
        overwrite = TRUE,
        keep.minimal = FALSE
        )  
    
   }
         
)

names(motifs.ls) <- levels(x.after.sp@cluster)
saveRDS(motifs.ls,file = 'DARs_exactTest/homer/motifs.rds')



##https://github.com/r3fang/SnapATAC/blob/master/examples/10X_brain_5k/README.md#diffusion_maps
##GREAT function annotation
cluster.pos='8'

library(rGREAT);
# DARs = findDAR(
#     obj=x.sp,
#     input.mat="pmat",
#     cluster.pos=13,
#     cluster.neg.method="knn",
#     test.method="exactTest",
#     bcv=0.1, #0.4 for human, 0.1 for mouse
#     seed.use=10
#   );
# DARs$FDR = p.adjust(DARs$PValue, method="BH");
# idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
idy = idy.ls[[cluster.pos]]
job = submitGreatJob(
    gr                    = x.after.sp@peak[idy],
    bg                    = NULL,
    species               = "hg38",
    includeCuratedRegDoms = TRUE,
    rule                  = "basalPlusExt",
    adv_upstream          = 5.0,
    adv_downstream        = 1.0,
    adv_span              = 1000.0,
    adv_twoDistance       = 1000.0,
    adv_oneDistance       = 1000.0,
    request_interval = 300,
    max_tries = 10,
    version = "default",
    base_url = "http://great.stanford.edu/public/cgi-bin"
  );
job





#step 18 predict gene-enhancer pairs (like cicero?)

TSS.loci = GRanges("chr12", IRanges(8219067, 8219068));
pairs = predictGenePeakPair(
	x.after.sp, 
	input.mat="pmat",
	gene.name="C3AR1", 
	gene.loci=resize(TSS.loci, width=500000, fix="center"),
	do.par=FALSE
	);
# convert the pair to genome browser arc plot format
pairs.df = as.data.frame(pairs);
pairs.df = data.frame(
	chr1=pairs.df[,"seqnames"],
	start1=pairs.df[,"start"],
	end1=pairs.df[,"end"],
	chr2="chr2",
	start2=8219067,
	end2=8219068,
	Pval=pairs.df[,"logPval"]
	);
head(pairs.df)
##    chr1  start1    end1 chr2  start2    end2       Pval
## 1 chr12 7984734 7985229 chr2 8219067 8219068 14.6075918
## 2 chr12 7987561 7988085 chr2 8219067 8219068  5.6718381
## 3 chr12 7989776 7990567 chr2 8219067 8219068 24.2564608
## 4 chr12 7996454 7996667 chr2 8219067 8219068  0.6411017
## 5 chr12 8000059 8000667 chr2 8219067 8219068  2.0324922
## 6 chr12 8012404 8013040 chr2 8219067 8219068  0.0000000

saveRDS(pairs.df,"pairs_gene_enhancer.rds")



######final save with bmat gmat pmat mmat snap obj#####
#saveRDS(x.after.sp,"x.after.sp.final.rds")
#saveRDS(x.after.sp,"x.after.sp.final.rds")
#saveRDS(x.after.sp,"x.after.sp.onestep_filter.final.bk.rds")

#write.table(cluster.df.add[,1:3],file='snapATAC.PLA-term-ATAC.umap.bin5k.final.txt',sep='\t',quote = FALSE, row.names = TRUE, col.names = TRUE)
##saveRDS(cluster.df.add,"cluster.df.add.final.rds")





##########reload#########################
x.after.sp <- readRDS("x.after.sp.final.rds")
cluster.df.add <- readRDS("cluster.df.add.final.rds")


x.after.sp <- readRDS("x.after.sp.cstb.rds")
cluster.df.add <- readRDS("cluster.df.add.cstb.rds")


all.equal(rownames(cluster.df.add),rownames(x.after.sp@metaData) ) #TRUE

options(repr.plot.width = 7.5, repr.plot.height = 7.5)
plot(cluster.df.add[,2:3],pch = 19, cex=0.1)



#gmat <- readRDS('gmat.rds')
#all.equal(gmat , x.after.sp@gmat)  #not eq, do not know why

#x.after.sp@gmat <- gmat
#all.equal(gmat , x.after.sp@gmat)


###rewrite pmat with raw pmat
#pmat <- readRDS('pmat.raw.rds')
#max(pmat) #255
#max(x.after.sp@pmat) #1

#all.equal(pmat , x.after.sp@pmat) #TRUE

#x.after.sp@pmat <- pmat


#dim(pmat) #9336 x 173870 #12526 x 299708
#dim(x.after.sp@pmat) #9336 x 173870 #12526 x 299708
#all.equal(rownames(pmat) ,rownames(x.after.sp@pmat) ) #TRUE
#all.equal(colnames(pmat) ,colnames(x.after.sp@pmat) ) #TRUE


# pmat.b <- pmat
# pmat.b[pmat.b > 1] <- 1   
# all.equal(pmat.b , x.after.sp@pmat) #

##
#x.after.sp@pmat <- pmat



#saveRDS(x.after.sp,"x.after.sp.final.rds")

#x.after.sp <- readRDS("x.after.sp.final.rds")

####save modified object rds again#########



saveRDS(x.after.sp,"x.after.sp.cstb.rds")

saveRDS(cluster.df.add,"cluster.df.add.cstb.rds")
write.table(x = cluster.df.add,file = 'cluster.df.add.cstb.txt',quote = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)





###save a slim version to reduce mem usage after loading
#saveRDS(x.after.sp@gmat,"gmat.rds") #not sparse, as big as 6.5G,  12526 x 58344

#x.after.sp@gmat <- readRDS('empty.rds')
#saveRDS(x.after.sp,"x.after.sp.final.slim.rds") #2.0G #pmat is raw



