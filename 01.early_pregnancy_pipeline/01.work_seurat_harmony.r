#library(Seurat,lib.loc = "/home/mjwang/.conda/envs/myenv/lib/R/library_opt/") #3.2.2

##upstream analysis of seurat

library(Seurat) #3.2.3
library(viridis)
library(magrittr)
library(dplyr)
library(tidyr)
library(Matrix)
library(ggplot2)
#library(matrixStats)
library(ComplexHeatmap)
#library(viridis)
#library('ggplot2')

#source("../../04.scRNA_scATAC/greenleaf_archR_doublets/ArchR/R/ColorPalettes.R")

library(DoubletFinder)

library(harmony) #Rcpp 1.0.8.3, harmony 1.0

library(hrbrthemes)

library(patchwork)

#library("writexl")
library('openxlsx')

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
names(color_snap) <- NULL

color_signac = c(
'0'='#E6D55E','1'='#792B8A','2'='#DA703D','3'='#9DC8E5','4'='#BA273C','5'='#C2C184','6'='#7F8084','7'='#65AB53','8'='#D082AF','9'='#496EAB','10'='#DE896D','11'='#491F8B','12'='#E1AD49','13'='#8E1B85','14'='#E7EE77','15'='#7D1A1D','16'='#96B355')
names(color_signac) <- NULL

color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" )

color_gradient_my <- c(
    rgb(5,48,97,maxColorValue = 255),
    rgb(42,113,178,maxColorValue = 255),
    rgb(147,198,222,maxColorValue = 255),
    rgb(239,243,245,maxColorValue = 255),
    rgb(253,219,199,maxColorValue = 255),
    rgb(214,96,77,maxColorValue = 255),
    rgb(121,5,34,maxColorValue = 255)

)

color_pyreds <- c(
"#fff5f0","#fef4ef","#fef3ee","#fef3ed","#fef2ec","#fef1eb","#fef1ea","#fef0e9","#feefe8","#feefe7","#feeee6","#feede5","#feede4","#feece3","#feebe2","#feebe1","#feeae0","#fee9e0","#fee9df","#fee8de"
)


#color_cellranger <- c("#820610",'#C60F1E','#F7592D','#FBB33D','#FDEC4A','#DFFB52','#B4FB68','#8BFB89','#70FBA7','#52FDCE','#3DEEF6','#33CEE5','#2BAED7','#2392CA','#1A70B9','#1354AB','#0A2E94','#061D88','#07177E')


color_cellranger <-c('#820610','#C50F1E','#F42428','#F86D30','#FBB33D','#FCFB4E','#C0FB61','#87FB8F','#41FEF9','#2BAED7','#155CB1','#08238D')




##color palette 1: warm (similar to ATAC), use this?
map_cellcolor_rna <- list(
  '10' = '#67001f',
 '3'= '#b2182b',
 '8'= '#FCC140',#'#FAC669',#'#EABFE1',#'#FBC04D',#'#d6604d',
 '1'= '#f4a582',
 '4'= '#FB8D3C',#'#FFA300',#'#fddbc7',
 '2'= '#d6604d',#''#AB855A',#'#7A5BA1',#'#f7f7f7',
 '5'='#d6604d', #'#d1e5f0',
 '7'= '#92c5de',
 '6'= '#4393c3',
 '11'= 'darkgreen',#'#2166ac',
 '9'= '#053061'
)

# ##temp
# rna_map_cellcolor <- list(
#   '10' = '#67001f',
#  '3'= '#b2182b',
#  '8'= '#EABFE1',#'#FBC04D',#'#d6604d',
#  '1'= '#f4a582',
#  '4'= '#FFA300',
#  '2'= '#d6604d',#''#AB855A',#'#7A5BA1',#'#f7f7f7',
#  '5'='#d6604d', #'#d1e5f0',
#  '7'= '#92c5de',
#  '6'= '#46A040',
#  '11'= '#FFA300',#'#2166ac',
#  '9'= '#00441B'
# )

options(repr.plot.height = 7.5, repr.plot.width = 7.5)
barplot(1:length(map_cellcolor_rna),col = unlist(map_cellcolor_rna),main ='rna_map_cellcolor' ,cex.main=2,names.arg = names(map_cellcolor_rna))

saveRDS(map_cellcolor_rna,'map_cellcolor_rna.hotcolor.rds')

#color palette 2: cold-purple (use this?)

###rna clusters
reds <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf')
blues <- c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','darkblue')
#blues_darker <- darker(blues, 0.3)

oranges <- c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6')
purples <- c('#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')

# colorset <- blues
# colorset <- blues_darker
# barplot(1:length(colorset),col=colorset,names.arg = 1:length(colorset), las =2 )


map_cellcolor_rna <- list(
           
    '1'=purples[2],#'STB1',
    '2'=purples[3],#'STB5', 
    '3'=purples[4],#'STB2',
    '4'='#C390D4',#blues[2],#'STB3',
    '5'=blues[5],#'CTB-2',
    '8'=purples[1],#naive STB',
    '7'= '#94C5DD',#'#4895C1',#'#92c5de',#'darkgreen',#'STB4',
    '6'= '#4977B0',#'#94C5DD', #'#46A040',#blues[3],#'lightblue',#'CTB-1',
    '9'= '#083160',#'#00441B',#blues[4],#'CTB-3',
    '11'='#116314',#'#FFA300',#'darkgreen',#'darkorange',#'Fusion component'
    '10'= purples[5]
    
    #yellow for fusion    
#     '1'=purples[2],#'STB1',
#     '2'=purples[3],#'STB5', 
#     '3'=purples[4],#'STB2',
#     '4'='#C390D4',#blues[2],#'STB3',
#     '5'=blues[5],#'CTB-2',
#     '8'=purples[1],#naive STB',
#     '7'= blues[4],#'#92c5de',#'darkgreen',#'STB4',
#     '6'='#46A040',#blues[3],#'lightblue',#'CTB-1',
#     '9'='#00441B',#blues[4],#'CTB-3',
#     '11'='#FFA300',#'darkgreen',#'darkorange',#'Fusion component'
#     '10'= purples[5]

)


options(repr.plot.height = 7.5, repr.plot.width = 7.5)
barplot(1:length(map_cellcolor_rna),col = unlist(rna_map_cellcolor),main ='rna_map_cellcolor' ,cex.main=2,names.arg = names(rna_map_cellcolor))

saveRDS(map_cellcolor_rna,'map_cellcolor_rna.coldpurple.rds')


#color palette 3: discrete


cellcolor <- list(
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


map_cellcolor_rna <- list(
    
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


options(repr.plot.height = 7.5, repr.plot.width = 7.5)
barplot(1:length(map_cellcolor_rna),col = unlist(rna_map_cellcolor),main ='rna_map_cellcolor' ,cex.main=2,names.arg = names(rna_map_cellcolor))




######quick look and save cluster.df.add with color
source('quickDimPlot_labelon.r')

options(repr.plot.width = 7.5, repr.plot.height=8.5)
quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', color_use = map_cellcolor_rna, title= 'early combined-RNA', shrink.x = 3, shrink.y = 1.5)

ggsave(filename = "pdfs/UMAP/PLA-early-RNA-UMAP.quickplot.pdf",height=7.5,width=8.5,useDingbats=FALSE)




##cluster name
map_cellname_rna <- list(
    
    '1'='STB Pre-mature 2',
    '2'='STB Pre-mature 1-1',
    '3'='STB Mature 2',
    '4'='STB Pre-mature 1-2',
    '5'='STB Mature 1-2',
    '6'= 'CTB-2',
    '7'= 'CTB Proliferation',
    '8'= 'STB Nascent',
    '9'= 'CTB-1',
    '10'='STB Mature 1-1',
    '11'= 'CTB Fusion'

)


########quantitative colors
###buencolors
#library(BuenColors)
palcolor <- 'brewer_red'
#palcolor <- 'brewer_heat' #BuenColors::jdb_palette
#palcolor <- 'brewer_purple'
#palcolor <- 'brewer_blue'
#palcolor <- 'brewer_green'
#palcolor <- 'brewer_violet'
#palcolor <- 'brewer_yes' 

#color_use <- BuenColors::jdb_palette(palcolor)

color_use <- colorRampPalette(BuenColors::jdb_palette(palcolor))(256)


###select one global color set###
color <- color_good



####
sample = "PLA-early_combined-RNA"


##eload object

#placenta <- readRDS( "PLA-8w-RNA.final.rds")




#####


file.list <- list(

   'PLA-rna-early1' = '../../placenta_10X_early1/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix.h5',
    'PLA-rna-early2' = '../../placenta_10X_early2/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix.h5',
    'PLA-rna-early3' = '../../placenta_10X_early3/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix.h5',
    'PLA-rna-early4' = '../../placenta_10X_early4/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix.h5',
    'PLA-rna-early5' = '../../placenta_10X_early5/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix.h5',
    'PLA-rna-early6' = '../../placenta_10X_early6/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix.h5'

)


file.rds <- list( #rds without any filtering

   'PLA-rna-early1' = '../../placenta_10X_early1/02.seurat/PLA-early1-RNA.rds',
    'PLA-rna-early2' = '../../placenta_10X_early2/02.seurat/PLA-early2-RNA.rds',
    'PLA-rna-early3' = '../../placenta_10X_early3/02.seurat/PLA-early3-RNA.rds',
    'PLA-rna-early4' = '../../placenta_10X_early4/02.seurat/PLA-early4-RNA.rds',
    'PLA-rna-early5' = '../../placenta_10X_early5/02.seurat/PLA-early5-RNA.rds',
    'PLA-rna-early6' = '../../placenta_10X_early6/02.seurat/PLA-early6-RNA.rds'

)



dir.list <- list(

    'PLA-rna-early1' = '../../placenta_10X_early1/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix',
    'PLA-rna-early2' = '../../placenta_10X_early2/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix',
    'PLA-rna-early3' = '../../placenta_10X_early3/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix',
    'PLA-rna-early4' = '../../placenta_10X_early4/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix',
    'PLA-rna-early5' = '../../placenta_10X_early5/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix',
    'PLA-rna-early6' = '../../placenta_10X_early6/01.data_cellranger_RNASEQ/filtered_feature_bc_matrix'

)

#for (i in names(dir.list)){   }


#data.dir = '../01.data_cellranger_RNASEQ/filtered_feature_bc_matrix/'
#data.dir = '../01.aggregation/PLA-8w-RNA-aggre_nonorm/filtered_feature_bc_matrix/'
#data.dir = '../01.aggregation/PLA-8w-RNA-aggre/filtered_feature_bc_matrix/'
#data.dir = "../../placenta_10X/01.data_cellranger_RNASEQ/PLA-8w-RNA-1/filtered_feature_bc_matrix/"
#data.dir = "../01.data_cellranger_RNASEQ/filtered_feature_bc_matrix/"
##placenta.count <- Read10X(data.dir = data.dir)
##placenta <- CreateSeuratObject(counts = placenta.count, project = sample, min.cells = 3, min.features = 200)
#26139 features across 7613
#24307 x 15407


##read seurat obj
placenta.list <- list()

for(i in names(file.rds)){
    placenta.list[[i]] <- readRDS(file.rds[[i]]) #use preprocessed rds (but no filtering)
    
}

gc()


#check colnames

for (i in names(placenta.list)){
    
    cat( dim(placenta.list[[i]])[1],'\t',dim(placenta.list[[i]])[2],'\n' )
    
}

#gene x cell
22473 	 9496 
22720 	 5911 
26139 	 7613 
24915 	 8104 
23639 	 6856 
24596 	 9285

#check umap
for (i in names(placenta.list)){
    
    cat( i,'\n' )
    options(repr.plot.height=7.5,repr.plot.width=7.5)
    res.p <- DimPlot(object = placenta.list[[i]], label = TRUE,cols=c(color_good),pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() + ggtitle (i)
    print(res.p)   
}

#check quality and marker gene
for (i in names(placenta.list)){
    cat( i,'\n' )
    options(repr.plot.height=15,repr.plot.width=15)
    res.p <- FeaturePlot(placenta.list[[i]], 
                         features = c('DNMT1','ERVFRD-1','ESRRG','SH3TC2','PAPPA','FLT1','STAT4','MYCN','CDKN1A','HBZ'), 
                         reduction = "umap",
                         pt.size = .1,
                         slot = 'scale.data',
                         min.cutoff = 'q10',
                         max.cutoff = 'q99'
                   ) #+ scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
    print(res.p)
    
#     for(j in c('DNMT1','ERVFRD-1','SH3TC2','PAPPA','FLT1','STAT4','MYCN')){
#         res.p <- FeaturePlot(placenta.list[[i]], features = j, reduction = "umap",pt.size = .1,
#                     slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99'
#                    ) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
#         print(res.p)

#     }
}







##merge obj, similar to scanpy concat? similar to cellranger aggre? will do intersection if rownames (gene) diff
placenta <- merge(  placenta.list[['PLA-rna-early1']], 
               y = c(placenta.list[['PLA-rna-early2']],
                     placenta.list[['PLA-rna-early3']],
                     placenta.list[['PLA-rna-early4']],
                     placenta.list[['PLA-rna-early5']],
                     placenta.list[['PLA-rna-early6']]
                    ), 
               add.cell.ids = c('PLA-rna-early1','PLA-rna-early2','PLA-rna-early3','PLA-rna-early4','PLA-rna-early5','PLA-rna-early6'), 
               project = "PLA-rna-early-combined",
               merge.data = TRUE #merge data slot and remove scale.data, recommend by Seurat::merge manual
            )

#29132 x 47265 (no any filtering)


##clean some columns in merged obj, leaving room for new column


placenta[['seurat_clusters']] <- NULL
placenta[['cluster']] <- NULL


str(placenta)
saveRDS(placenta,'placenta.raw.rds')

table(placenta[['orig.ident']])
PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
          9496           5911           7613           8104           6856 
PLA-early6-RNA 
          9285 


all.equal(placenta[['orig.ident']],placenta[['sample']],check.attributes = FALSE) #TRUE


table(placenta[['sample']][,1],placenta[['sex']][,1])
         
                 female male
  PLA-early1-RNA      0 9496
  PLA-early2-RNA   5911    0
  PLA-early3-RNA   7613    0
  PLA-early4-RNA      0 8104
  PLA-early5-RNA      0 6856
  PLA-early6-RNA   9285    0


sex.list <- list(
  'PLA-early1-RNA' = 'male',
  'PLA-early2-RNA' =  'female',
  'PLA-early3-RNA' =  'female',
  'PLA-early4-RNA' =  'male',
  'PLA-early5-RNA' =  'male',
  'PLA-early6-RNA' =  'female'

)


write.table(x = colnames(placenta),file = 'cellid.all.txt',sep = '\t',row.names = FALSE, col.names = FALSE,quote = FALSE)

##clean
rm(placenta.list)
gc()



##get qc percentage
#placenta[["percent.mt"]] = PercentageFeatureSet(placenta,pattern = "^MT-")
#placenta[["percent.ribo"]] = PercentageFeatureSet(placenta, "^RP[SL]")
#placenta[["percent.hb"]] <- PercentageFeatureSet(placenta, "^HB[^(P)]")
#placenta[["percent.xist"]] <- PercentageFeatureSet(placenta, "^XIST")

#grep(pattern = "^XIST", x = rownames(placenta), value = TRUE) #within PercentageFeatureSet


##get sex information , read in cellranger 7 gene.gtf chrY genes

# genes_chrY.df <- read.table("data/genes.chrY.bed") #111

# genes_chrY.split.df <- genes_chrY.df %>% separate(col = V4, sep = "\\|", into = c('geneid','gene_type') )

# colnames(genes_chrY.split.df) <- c('chr','start','end','geneid','gene_type','strand')

# genes_chrY.split.protein.df <- subset(genes_chrY.split.df, gene_type == 'protein_coding')
# #48

# chrY.genes <- genes_chrY.split.protein.df$geneid
# sum(duplicated(chrY.genes) ) #0



# chrY.genes.sel <- chrY.genes[ chrY.genes %in% rownames(placenta)  ]
# placenta[['percent.chrY']] = colSums(placenta@assays$RNA@counts[chrY.genes.sel,])/colSums(placenta@assays$RNA@counts)
    
    

# #####cell cycle scoring#####
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# placenta <- CellCycleScoring(placenta, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)



#####
options(repr.plot.height = 7.5, repr.plot.width=15)
VlnPlot(placenta,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1,group.by = 'orig.ident')

plot1 <- FeatureScatter(placenta, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'orig.ident',pt.size = .1)
plot2 <- FeatureScatter(placenta, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'orig.ident',pt.size = .1)
options(repr.plot.height=6,repr.plot.width=15)
CombinePlots(plots = list(plot1, plot2))

# ##filter now or later??
# #placenta <- subset(placenta, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 ) #use 5000 instead of 5500, 
# placenta <- subset(placenta, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 & nCount_RNA >3500 & nCount_RNA < 5e4 & percent.mt < 5) #the same condition in scanpy

# #after filter #24307 x 11560  #5646 #6471
# options(repr.plot.height=7.5,repr.plot.width=12.5)
# VlnPlot(placenta,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1 )



###filter apoptosis cluster 9

table(Idents(placenta))
  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
3455 3107 3089 2846 2829 2700 2084 1887 1577 1256  549  490  396  375  232  200 
  16   17 
 136   36

DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() 

placenta <- subset(placenta, idents = levels(Idents(placenta))[!levels(Idents(placenta)) %in% c('9')]  )


table(Idents(placenta))
  0    1    2    3    4    5    6    7    8   10   11   12   13   14   15   16 
3455 3107 3089 2846 2829 2700 2084 1887 1577  549  490  396  375  232  200  136 
  17 
  36

DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() 
#29132 x 25988 gene x cell

##normalize
#placenta <- NormalizeData(placenta, normalization.method = "LogNormalize", scale.factor = 10000)

##identification of highly variable features
placenta <- FindVariableFeatures(placenta, selection.method = "vst", nfeatures = 2000) #(reanalysis from here?)
placenta <- FindVariableFeatures(placenta, selection.method = "vst", nfeatures = 2500)

#this step may be the critical step. if cells were filtered, should rerun this step

mostVar_before <- VariableFeatures(placenta) #2000
mostVar_after <- VariableFeatures(placenta) #3500




mostVar_filter #3176

length(intersect(mostVar_after, mostVar_filter)) #1915




top30 = head(VariableFeatures(placenta),30  )
plot1 = VariableFeaturePlot(placenta   )
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
options(repr.plot.height=7.5,repr.plot.width=15)
CombinePlots(plots = list(plot1, plot2))

##scale the data
all.genes <- rownames(placenta)
placenta <- ScaleData(placenta, features = all.genes) #scale.data 29132 x 27244



##filtering mostVar genes for Mt ribo sex cell-cycling genes?


mostVar <- VariableFeatures(placenta) #2000


flag_mt <- grep(pattern = "^MT-", x = mostVar)#, value = TRUE)
flag_ribo <- grep(pattern = "^RP[SL]", x = mostVar)#, value = TRUE)


#get sex information , read in cellranger 7 gene.gtf chrY genes (not chrX!)
genes_chrY.df <- read.table("data/genes.chrY.bed") #111
genes_chrY.split.df <- genes_chrY.df %>% separate(col = V4, sep = "\\|", into = c('geneid','gene_type') )
colnames(genes_chrY.split.df) <- c('chr','start','end','geneid','gene_type','strand')
sex_genes <- c("XIST",genes_chrY.split.df$geneid) #all chrY genes include protein_coding
flag_sex <- which(mostVar %in% sex_genes)

##cell cycling gene
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


flag_cellcycle_1 <- which(mostVar %in% s.genes)
flag_cellcycle_2 <- which(mostVar %in% g2m.genes)

mostVar[as.numeric (c(flag_mt,flag_ribo,flag_sex,flag_cellcycle_1,flag_cellcycle_2) )]
#removed: 'RPS6KA2''RPL3L''PCDH11Y''AC010737.1''TTTY14''LINC00278''NLGN4Y''USP9Y''XIST''UTY''RRM2''BRIP1''DTL''ATAD2''CLSPN''UHRF1''POLA1''BLM''CENPE''TOP2A''CENPF''GTSE1''GAS2L3''KIF23''NUSAP1''CDCA2''CKS2''ECT2''ANLN''MKI67''CDC25C''KIF11''TPX2''HMMR''KIF2C''TACC3''SMC4'

mostVar_filter <- mostVar[-as.numeric (c(flag_mt,flag_ribo,flag_sex,flag_cellcycle_1,flag_cellcycle_2) )]
#'RPS6KA2''PCDH11Y''AC010737.1''TTTY14''LINC00278''NLGN4Y''RRM2''BRIP1''DTL''CENPE''TOP2A''CENPF''GTSE1''GAS2L3''KIF23''NUSAP1'
#3463 left
#1984 left

#linear dimensional reduction: PCA
#placenta <- RunPCA(placenta, features = VariableFeatures(object = placenta)) #total 50 PCs

mostVar_filter <- read.table('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.scanpy_harmony/highly_variable_genes.filter.txt',sep = '\t',stringsAsFactors = FALSE)$V1
#3176
#

placenta <- RunPCA(placenta, features = mostVar_filter)


#quick look the PCs
options(repr.plot.height=7.5,repr.plot.width=7.5)
VizDimLoadings(placenta, dims = 1:5, reduction = "pca")
DimPlot(placenta, reduction = "pca",group.by = 'orig.ident')
DimHeatmap(placenta, dims = 1:15, cells = 500, balanced = TRUE)

#choose PC number
placenta <- JackStraw(placenta, num.replicate = 100,dims = 35) #20 by default
placenta <- ScoreJackStraw(placenta, dims = 1:25)
#JackStrawPlot(placenta, dims = 1:15)
ElbowPlot(placenta,ndims = 40) #use n_pcs = 30

m

######## plot cluster size##
Idents(placenta)
all.equal(Idents(placenta),placenta@meta.data$seurat_clusters) #name id diff
cluster.stat = table(placenta@meta.data$seurat_clusters)
total <- sum(cluster.stat)  #15107 #9461
res.bp <- barplot(cluster.stat,col = color_good,las=2,names.arg = paste0("cluster",names(cluster.stat)),
                  ylab = "cell count", main='cluster cell number' )
text(x=res.bp[,1],y=cluster.stat+55,labels = cluster.stat,srt=90,xpd=TRUE)
text(x=16,y=1500,labels = paste("total = ",total,sep='') )
###

#system("which python")
#Sys.getenv()
#Sys.setenv("PATH"="/home/mjwang/bin:/usr/local/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin")
#system("which python") #use /home/mjwang/bin/python
#reticulate::py_config()

##non-linear dimensional reduction (UMAP/tSNE)
placenta = RunUMAP(placenta,dims=1:30,seed.use = 123) 
placenta =  RunTSNE(placenta,dims=1:30,seed.use = 123)


#saveRDS(placenta, file = "placenta.PLA-8w-RNA-new.rds",compress = TRUE)
#placenta <- readRDS(file = "placenta.PLA-8w-RNA-new.rds") #dgCMatrix

options(repr.plot.height=7.5,repr.plot.width=7.5)
#DimPlot(placenta, reduction = "tsne",label=TRUE,cols = color_good,label.size = 8,pt.size = 1)
DimPlot(placenta, reduction = "umap",label=TRUE,cols=c(color_good,color_good),label.size = 8,pt.size = 1)

options(repr.plot.height=7.5,repr.plot.width=15)
DimPlot(placenta, reduction = "umap",split.by = 'sample' ,label=FALSE,cols=c(color_good,color_good),label.size = 8,pt.size = 1)

options(repr.plot.height=7.5,repr.plot.width=12)
DimPlot(placenta, reduction = "umap",split.by = 'sex' ,label=FALSE,cols=c(color_good,color_good),label.size = 8,pt.size = 1)

# ####detect and mark doublet
# # define the expected number of doublet cellscells.
# nExp <- round(ncol(placenta) * 0.1)  #0.1,   #0.04  #301,  expect 4% doublets
# placenta <- doubletFinder_v3(placenta, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:30)#1:10)




# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(placenta@meta.data)[grepl("DF.classification", colnames(placenta@meta.data))]

table(placenta[[DF.name[1]]])
Doublet Singlet 
    761    6852

# Doublet Singlet 
#     301    7236

# table(placenta[[DF.name[2]]])
# Doublet Singlet 
#     754    6783


options(repr.plot.width= 12, repr.plot.height = 7.5)
res.p <- cowplot::plot_grid(ncol = 2, 
                   DimPlot(placenta, group.by = "seurat_clusters") + NoAxes(),
                   DimPlot(placenta, group.by = DF.name,pt.size = 2) + NoAxes()
                  )
print(res.p)

options(repr.plot.width= 7.5, repr.plot.height = 7.5)
DimPlot(placenta, group.by = DF.name,pt.size = 1) + NoAxes()

# res.p <- VlnPlot(placenta, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
# print(res.p)



##filter for singlet
##placenta = placenta[, placenta@meta.data[, DF.name] == "Singlet"]
#dim(placenta)
##placenta



##########recalculate with harmony remove batch effect###
V <- placenta@reductions$pca@cell.embeddings
27244 x 50 #reanalysis after filtering
47265 x 50
#11560 x 50
#15107 x 50

meta_data <- placenta@meta.data
27244 x 32
47265 x 29

#11560 x 6
#15107 x 6

#meta_data$sample = ifelse(grepl(pattern = "-1$", x=rownames(meta_data) ), 'D1', 'D2'  )
#meta_data$sample <- factor(meta_data$sample,levels=c('D1','D2'))
table(meta_data$sample)
PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
          5896           4079           4916           3992           3144 
PLA-early6-RNA 
          5217 

PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
          9496           5911           7613           8104           6856 
PLA-early6-RNA 
          9285

#  D1   D2 
#6996 4564 

#  D1   D2 
#9461 5646 

V_harmony <- HarmonyMatrix(V,meta_data, 'sample', do_pca = FALSE,max.iter.cluster = 20)

##should be max.iter.harmony = 20!

#Harmony converged after 10 iterations

#Harmony converged after 8 iterations
#47265 x 50

#Harmony converged after 6 iterations
#11560 x 50

#do not converged after 10 iterations
#15107 x 50

placenta@reductions$pca@cell.embeddings <- V_harmony #overwrite the PCA


#DefaultAssay(placenta) <- 'peaks'

#placenta <- RunUMAP(object = placenta, reduction = 'pca', dims = 1:30)
placenta = RunUMAP(placenta,dims=1:30,seed.use = 123,min.dist = 0.5,spread=1.2) 



#######tuning UMAP embedding#####

seed.use = 123
min.dist = 0.5
spread=1.2


# random_state = 0
# n_comps = 2
# min.dist = 0.
# spread = 1.0

umap.bk <- Embeddings(placenta,reduction = 'umap')

for (min.dist in c(0.1,0.2,0.3,0.4,0.5)){
#for mdist in [0.3]:
    for (spread in c(0.5,1.0,1.5,2.0)){
#    for spread in [1.0]:
      placenta = RunUMAP(placenta,dims=1:30,
                         seed.use = seed.use,
                         min.dist = min.dist,
                         spread=spread,
                         return.model=TRUE,
                         verbose=0
                        ) 
     options(repr.plot.height=7.5,repr.plot.width=7.5)
     res.p <- DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() + ggtitle(paste( c('umap tunning: mdist ',min.dist,' spread ',spread,' seed.use ',seed.use) ,collapse =' '))
     print(res.p)
   } 
    
    
}



##fix umap embedding

seed.use = 123
min.dist = 0.3
spread=0.5

placenta = RunUMAP(placenta,dims=1:30,
                 seed.use = seed.use,
                 min.dist = min.dist,
                 spread=spread,
                 return.model=TRUE,
                 verbose=0
                ) 


options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() #+ ggtitle(paste( c('umap tunning: mdist ',min.dist,' spread ',spread,' seed.use ',seed.use) ,collapse =' '))



############rotate the umap xy coordinate if needed#############
#test <- data.frame(x=c(10,2),y=c(25,45))
#plot(test,xlim = c(-50,50),ylim = c(-50,50),cex=2 );abline(h = 0);abline(v=0)
#test.rotate90 <- data.frame(x=test$x*cos(90)-test$y*sin(90),y=test$y*cos(90)+test$x*sin(90))
#points(test.rotate90,cex=2,pch=19,col='red')

umap.df <- data.frame(UMAP_1=Embeddings(placenta,reduction = 'umap')[,1],UMAP_2=Embeddings(placenta,reduction = 'umap')[,2])
#rownames(umap.df) <- colnames(placenta)

##iterative rotate the umap direction
#for(degree in seq(from = 50,90,5) ) { 
#for(degree in seq(from = 65,90,2) ) { 
for(degree in c(-90) ) { 
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

fix degree = -90
#fix degree = 175
#fix degree = 150
##fix degree = 175



#rewrite the umap slot
rownames(UMAP.rotate) <- rownames(Embeddings(placenta,reduction = 'umap'))
UMAP.bk = Embeddings(placenta,reduction = 'umap')
#Embeddings(placenta,reduction = 'umap_rotate') = as.matrix(UMAP.rotate) #error


placenta[['umap_rotate']] <- CreateDimReducObject(embeddings = as.matrix(UMAP.rotate), key = 'umaprotate_', assay = 'RNA')

#saveRDS(object=UMAP.rotate,file = 'uwot.umap.seed177.rotate175.rds')

options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap_rotate") + NoLegend() #+ ggtitle(paste( c('umap tunning: mdist ',min.dist,' spread ',spread,' seed.use ',seed.use) ,collapse =' '))



#################





###do clustering

placenta <- FindNeighbors(object = placenta, reduction = 'pca', dims = 1:30)
placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 1,resolution=0.9) #ori louvain
#placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 2,resolution=1.2) #enhanced louvain
##placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 4,resolution=1.2) #leiden (problem too large, why?)


#####tuning clustering###########

#louvain
#for (res in seq(0.3,0.7,0.1) ){
#for (res in seq(0.8,1.5,0.1) ){
#for (res in seq(1.5,2.5,0.1) ){

#leiden
#for (res in seq(0.3,0.7,0.1) ){
for (res in seq(0.7,2.5,0.1) ){
    #placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 1,resolution=res)
    placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 4,resolution=res) #slow
    
    options(repr.plot.height=7.5,repr.plot.width=7.5)
    res.p <- DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() + ggtitle(paste( c('resolution ',res) ,collapse =' '))
     print(res.p)
    
}

#res = 0.65
res = 1.5 #louvain
placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 1,resolution=res)

res = 1.4 #0.65 #leiden
placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 4,resolution=res)

table(Idents(placenta))
  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
3311 2556 2284 2275 2134 1895 1777 1776 1693 1532 1301 1119  536  420  409  316 
  16   17   18   19 
 203  184  136  131 

   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
4067 4032 3928 3705 3080 1809 1695 1361  535  420  404  361  203  183  136   47 
  16  #louvain
  22 


   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
2787 2347 2182 2001 2001 1928 1694 1670 1606 1492 1461 1014  807  665  546  420 
  17   18   19   20   21   22   23 
 406  317  203  184  136   97   24 

1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16  #leiden
5129 4526 3913 3336 3077 1993 1713  536  420  388  364  203  184  136   48   22 

options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() + ggtitle(paste( c('resolution ',res) ,collapse =' '))

##enlarge
options(repr.plot.height=15.5,repr.plot.width=15.5)
DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.2,label.size = 15,reduction = "umap") + NoLegend() + ggtitle(paste( c('resolution ',res) ,collapse =' '))


# ##plot by cluster (see below ggplot code)

# placenta <- AddMetaData(placenta,metadata = Idents(placenta), col.name = 'cluster')

# options(repr.plot.height=5.5,repr.plot.width=50)
# DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap",split.by='cluster') + NoLegend()


##quality check for each clusters
options(repr.plot.height=7.5,repr.plot.width=12.5)
VlnPlot(placenta,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1 ,group.by = 'orig.ident')

options(repr.plot.height=5.5, repr.plot.width = 25)
DimPlot(object = placenta, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap", group.by = 'sample',split.by = 'sample') + NoLegend()


options(repr.plot.height=5.5, repr.plot.width = 10)
DimPlot(object = placenta, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap", group.by = 'sex',split.by = 'sex') + NoLegend()



#############quick annotation############
# options(repr.plot.height=15,repr.plot.width=15)
# FeaturePlot(placenta, features = c("CGA", "PSG1","DNMT1","PAGE4",'MKI67','ERVFRD-1','LAIR2','PLAC8',"VIM"), reduction = "umap")

# FeaturePlot(placenta, features = c('CSH1','CSH2',"CSHL1",'PAPPA','GH2','FLT1','INSIG1', "LEP"), reduction = "umap")

# FeaturePlot(placenta, features = c('CSH1','CSH2',"CSHL1",'PAPPA','GH2','FLT1','INSIG1', "PRL"), reduction = "umap")
# FeaturePlot(placenta, features = c('nCount_RNA','nFeature_RNA','percent.mt'), reduction = "umap",cols = viridis(12,option = 'D') )#c('black','red') )



###plot qc values
options(repr.plot.height=7.5,repr.plot.width=7.5)
FeaturePlot(placenta, features = 'nFeature_RNA', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey90',color_use)) + NoAxes() ; 
FeaturePlot(placenta, features = 'nCount_RNA', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'percent.mt', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'percent.ribo', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'percent.hb', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'percent.chrY', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'percent.xist', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'XIST', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


##plot DoubletFinder 
nameid = colnames(placenta@meta.data)[grepl("pANN", colnames(placenta@meta.data))]
FeaturePlot(placenta, features = nameid, reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use))

nameid = colnames(placenta@meta.data)[grepl("DF.classification", colnames(placenta@meta.data))]
DimPlot(placenta, group.by = nameid,pt.size = .1) + NoAxes()


##sex related gene
sex_genes_sel <- sex_genes[sex_genes %in% row.names(placenta) ]

options(repr.plot.width=12.5, repr.plot.height = 12.5)
#options(repr.plot.width=7.5, repr.plot.height = 5.5)
#for(i in sex_genes_sel){
for(i in c('XIST','DDX3Y','USP9Y','RPS4Y1')){
  res.p <- FeaturePlot(placenta, features = i, reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
   print(res.p)
}


##sex gene by sex
options(repr.plot.width = 9, repr.plot.height = 4.5)
FeaturePlot(placenta, features = 'PCDH11Y', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes() #sex baised, because this gene has ortholog on chrX , will mismatch to chrY in female (with all chr for mapping)
FeaturePlot(placenta, features = 'USP9Y', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes() #few baise
FeaturePlot(placenta, features = 'RPS4Y1', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'DDX3Y', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'XIST', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

##selected gene by sex
options(repr.plot.width = 9, repr.plot.height = 4.5)
FeaturePlot(placenta, features = 'PAPPA', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'LAMA3', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'FLT1', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'ENG', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'CDKN1A', split.by = 'sex' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()



##sex gene by sample
options(repr.plot.width = 35, repr.plot.height = 6)
FeaturePlot(placenta, features = 'PCDH11Y', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1)# + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'USP9Y', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) #few sex biased 
FeaturePlot(placenta, features = 'RPS4Y1', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1)
FeaturePlot(placenta, features = 'DDX3Y', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q100',pt.size=.1)
FeaturePlot(placenta, features = 'XIST', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1)


##selected gene by sample 
options(repr.plot.width = 35, repr.plot.height = 6)
FeaturePlot(placenta, features = 'PAPPA', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q30', max.cutoff = 'q99',pt.size=.1)
FeaturePlot(placenta, features = 'LAMA3', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q30', max.cutoff = 'q99',pt.size=.1)
FeaturePlot(placenta, features = 'FLT1', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q50', max.cutoff = 'q99',pt.size=.1)
FeaturePlot(placenta, features = 'ENG', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q30', max.cutoff = 'q100',pt.size=.5)

FeaturePlot(placenta, features = 'CDKN1A', split.by = 'sample' ,reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1)


#########pixel style hexbin density scatterplot (grid)#######
##https://cran.r-project.org/web/packages/hexbin/vignettes/hexagon_binning.pdf
pdf( "pdfs_te/pixel_style.cluster.density.pdf",height=7.5,width=7.5,useDingbats = FALSE)

library(hexbin)
library(RColorBrewer)

options(repr.plot.height=7.5,repr.plot.width=7.5)
umap <- as.data.frame(Embeddings(placenta,reduction = 'umap'))
bin<-hexbin(umap$UMAP_1, umap$UMAP_2, xbins=40,xbnds = range(umap$UMAP_1)*1.3,ybnds = range(umap$UMAP_2)*1.3 )

my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))

hb <- plot(bin,main="Nuclei pileup density" , colramp=my_colors , legend=F,xlab='UMAP1',ylab='UMAP2')

dev.off()


##trophoblast
options(repr.plot.height=7.5,repr.plot.width=7.5)
FeaturePlot(placenta, features = 'KRT7', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use))
FeaturePlot(placenta, features = 'GATA3', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use))
FeaturePlot(placenta, features = 'TFAP2A', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use))


##CTB
FeaturePlot(placenta, features = 'DNMT1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'CDH1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99')+ scale_color_gradientn(colours = c('grey95',color_use))+ NoAxes()
FeaturePlot(placenta, features = 'PPARG', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


FeaturePlot(placenta, features = 'TEAD4', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'TEAD3', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


FeaturePlot(placenta, features = 'MKI67', reduction = "umap",pt.size = .5,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q100') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()#?


FeaturePlot(placenta, features = 'TP53', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'TP63', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes() #high
#FeaturePlot(placenta, features = 'TP73', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

#FeaturePlot(placenta, features = 'CDX2', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use))
FeaturePlot(placenta, features = 'BCAM', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()



##CTB fusion
options(repr.plot.height=7.5,repr.plot.width=7.5)
FeaturePlot(placenta, features = 'ERVFRD-1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q90') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'GCM1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'OVOL1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'REL', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()




##STB nascent
FeaturePlot(placenta, features = 'SH3TC2', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
#FeaturePlot(placenta, features = 'PDE4D', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes() #STB apoptosis?
FeaturePlot(placenta, features = 'BACE2', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'ESRRG', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'ESRRG', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


#FeaturePlot(placenta, features = 'PPARG', reduction = "umap",pt.size = 0.1,slot = 'data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

#FeaturePlot(placenta, features = 'PPARD', reduction = "umap",pt.size = 0.1,slot = 'data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


##STB general
options(repr.plot.height=7.5,repr.plot.width=7.5)

FeaturePlot(placenta, features = 'PSG8', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q100') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'CGA', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'PSG2', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'PSG5', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
#FeaturePlot(placenta, features = 'PPARG', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'LEP', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
#FeaturePlot(placenta, features = 'LIFR', reduction = "umap",pt.size = 2,slot = 'data',min.cutoff = 'q30',max.cutoff = 'q99')



##STB PAPPA
options(repr.plot.height=12.5,repr.plot.width=12.5)
FeaturePlot(placenta, features = 'PAPPA', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q100') + scale_color_gradientn(colours = c('grey90',color_use)) + NoAxes()
#FeaturePlot(placenta, features = 'ADAM18', reduction = "umap",pt.size = .1,slot = 'data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use))
FeaturePlot(placenta, features = 'ADAMTSL1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'ADAMTS6', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

#FeaturePlot(placenta, features = 'ANGPTL4', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'GH1', reduction = "umap_rotate",pt.size = .1,slot = 'data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'GH2', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'GHR', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey90',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'JAK1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'JAK2', reduction = "umap_rotate",pt.size = .2,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q100') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'LAMA3', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q100') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'AR', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'VDR', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'MITF', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q50',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'CSHL1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'CSH1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use))  + NoAxes()

FeaturePlot(placenta, features = 'CSH2', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use))  + NoAxes()#no



#for(i in grep('STAT',rownames(placenta),value=TRUE)){
for(i in grep('ADAM',rownames(placenta),value=TRUE)){
  res.p <- FeaturePlot(placenta, features = i, reduction = "umap",slot = 'scale.data',min.cutoff = 'q0', max.cutoff = 'q99',pt.size=.1) + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
   print(res.p)
}


FeaturePlot(placenta, features = 'STAT5A', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q60',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'STAT5B', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q50',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'STAT4', reduction = "umap",pt.size = .1,slot = 'data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'FOS', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'FOSB', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


FeaturePlot(placenta, features = 'JUNB', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'JUN', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()



##STB FLT1
options(repr.plot.height=12.5,repr.plot.width=12.5)
FeaturePlot(placenta, features = 'FLT1', reduction = "umap",pt.size = .2,slot = 'scale.data',min.cutoff = 'q50',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'ENG', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q20',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'ANGPTL4', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'FSTL3', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'INHBA', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'INHA', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'MYCN', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'POU2F3', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes() #?
FeaturePlot(placenta, features = 'LVRN', reduction = "umap",pt.size = .2,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q100') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'TGFB1', reduction = "umap",pt.size = .2,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

#for(i in grep('FOS',rownames(placenta),value=TRUE)){
#for(i in grep('JUN',rownames(placenta),value=TRUE)){
#for(i in grep('ATF',rownames(placenta),value=TRUE)){
#for(i in grep('SMAD',rownames(placenta),value=TRUE)){
for(i in grep('ANGP',rownames(placenta),value=TRUE)){
  res.p <- FeaturePlot(placenta, features = i, reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
  print(res.p)
}


FeaturePlot(placenta, features = 'FOSL2', reduction = "umap_rotate",pt.size = .2,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = 'JUND', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'JUNB', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()



FeaturePlot(placenta, features = 'FOSL1', reduction = "umap_rotate",pt.size = .2,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q100') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


#FeaturePlot(placenta, features = 'ATF3', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes() #?


##STB apoptosis
FeaturePlot(placenta, features = 'DDX60', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'DDX58', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
#FeaturePlot(placenta, features = 'MAPKAPK3', reduction = "umap",pt.size = 2,slot = 'data',min.cutoff = 'q0',max.cutoff = 'q99')
FeaturePlot(placenta, features = 'MAP4K4', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'SPATA5', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'GDF15', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
#FeaturePlot(placenta, features = 'CCDC30', reduction = "umap",pt.size = 2,slot = 'data',min.cutoff = 'q30',max.cutoff = 'q99')
#FeaturePlot(placenta, features = 'SASH1', reduction = "umap",pt.size = 2,slot = 'data',min.cutoff = 'q30',max.cutoff = 'q99')
#FeaturePlot(placenta, features = 'MDM2', reduction = "umap",pt.size = 2,slot = 'data',min.cutoff = 'q30',max.cutoff = 'q99')

FeaturePlot(placenta, features = 'CROT', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'CDKN1A', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'ADCY5', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

options(repr.plot.height=7.5, repr.plot.width = 7.5)
#for(i in marker.genes.de100[['3']][1:30]){
#for(i in marker.genes.de100[['2']][1:30]){
#for(i in  c('DDX60','DDX58','ZC3HAV1','GAS7','NMNAT2','CTIF'    ) ){
for(i in grep('BCL',rownames(placenta),value=TRUE)){
  res.p <- FeaturePlot(placenta, features = i, reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
  print(res.p)
}


##EVT

FeaturePlot(placenta, features = 'HLA-G', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'LAIR2', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'PLAC8', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'MMP2', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

FeaturePlot(placenta, features = '?', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()



##STR general
FeaturePlot(placenta, features = 'VIM', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'DLK1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

##Vascular Endothelial
FeaturePlot(placenta, features = 'PECAM1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

##STR
FeaturePlot(placenta, features = 'HIVEP3', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'HLA-A', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'HLA-DPA1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'HLA-DPB1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

##MesSTR
FeaturePlot(placenta, features = 'THY1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()

##Hofbauer cell (placental resident macrophage)
FeaturePlot(placenta, features = 'CD68', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'CD14', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


##red blood
FeaturePlot(placenta, features = 'HBA1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
FeaturePlot(placenta, features = 'HBZ', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()


for(i in grep('^HB',rownames(placenta),value=TRUE)){
  res.p <- FeaturePlot(placenta, features = i, reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99')  + scale_color_gradientn(colours = c('grey95',color_use)) + NoAxes()
  print(res.p)
}



##combined plot above genes with loop###

#get all genes above
#grep "^FeaturePlot" work_seurat.r|awk -vFS=" +" '{print $4}' |grep -v "nameid" |awk -vORS=" " '{print $1}'
marker.genes.full <- list(
    'Quality control' = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.ribo', 'percent.hb', 'percent.chrY', 'percent.xist','XIST'),
    'Trophoblast' = c('KRT7', 'GATA3', 'TFAP2A'),
    'CTB' = c('DNMT1', 'CDH1', 'PPARG', 'TEAD4', 'TEAD3', 'MKI67', 'TP53','TP63', 'TP73', 'BCAM'),
    'CTB fusion' = c('ERVFRD-1', 'GCM1', 'OVOL1','PPARD'), 
    'STB nascent' = c('SH3TC2', 'BACE2', 'ESRRG'), 
    'STB general' = c('PSG8', 'CGA', 'PSG2', 'PSG5', 'LEP'), 
    'STB PAPPA' = c('PAPPA', 'ADAMTSL1','ADAMTS6', 'GH2', 'GHR','JAK1', 'JAK2', 'LAMA3', 'AR', 'VDR', 'CSHL1', 'CSH1', 'CSH2', 'STAT5A', 'STAT5B', 'STAT4','FOS', 'FOSB', 'JUNB', 'JUN'), 
    'STB FLT1' = c('FLT1', 'ENG', 'ANGPTL4', 'FSTL3', 'INHBA', 'INHA','MYCN', 'POU2F3', 'LVRN', 'TGFB1', 'FOSL2', 'JUND'), 
    'STB apoptosis' = c('DDX60', 'DDX58', 'MAP4K4', 'SPATA5', 'GDF15', 'CROT', 'CDKN1A', 'ADCY5'),
    'EVT' = c('HLA-G', 'LAIR2', 'PLAC8', 'MMP2'), 
    'STR general' = c('VIM', 'DLK1'), 
    'Vascular Endothelial Cell' = c('PECAM1'), 
    'STR' = c('HIVEP3', 'HLA-A', 'HLA-DPA1', 'HLA-DPB1'), 
    'Mesenchymal STR'= c('THY1'), 
    'Hofbauer Cell'= c('CD68', 'CD14'), 
    'Red blood' = c('HBA1', 'HBZ')
)


hormone.peptide.gene <- list( #see ppt peptide hormone
   'Metabolism' = c('ADIPOQ','CALCA','CALCB','CRH','CSH1','CSH2','CSHL1','GAL','GCG','GH1','GH2','GHRH','GHRL','GIP','IGF1','INS','INSL3','INSL4','INSL5','INSL6','LEP','PPY','PRL','PTH','PTHLH','RETN','RETNLB','SOST','SST','STC1','STC2','TRH','TSHB','UCN','UCN2','VIP'),
    'VASCULAR' = c('ADM','ADM2','AGT','ANGPTL3','APLN','CTGF','EDN1','EDN3','NPPA','NPPB','NPPC','NTS2'),
    'BLOOD_CELL' = c('EPO','THPO'),
    'GUT_PEPTIDE' = c('CCK','GAST','MLN','SCT','TAC1'),
    'NEUROPEPTIDE' = c('BDNF','CROT','NMU','NTS','PMCH','PNOC'),
    'REPRODUCTION' = c('AMH','CGA','CGB7','FSHB','FST','GNRH2','GPHA2','INHA','INHBA','INHBB','INHBC','INHBE','KISS1','LHB','OXT','RLN1','RLN2')


)






#ERBB2: HER2
#Progestin , Progesterone gene???

#PR/PGR AR ESRRG



##plot all genes iteratively



marker.genes <- marker.genes.full
#marker.genes <- hormone.peptide.gene

outdir = 'pdfs'
#outdir = 'pdfs_final'
#outdir = 'pdfs_final_final'

if(dir.exists( paste0(outdir,'/marker_genes/') ) ){ 

 }else{
    dir.create(paste0(outdir,'/marker_genes/'),recursive=TRUE)
  }


# if(dir.exists( paste0(outdir,'/hormone_genes/') ) ){ 

#  }else{
#     dir.create(paste0(outdir,'/hormone_genes/'),recursive=TRUE)
#   }


for(id in names(marker.genes) ){
    
    cat('do for group ',id,'\n',sep='')
    
    id_simple <- gsub(pattern = " +",replacement = "_",x = id)
    

    if(dir.exists( paste0(outdir,'/marker_genes/',id_simple) ) ){  }else{dir.create(paste0(outdir,'/marker_genes/',id_simple))}

#     if(dir.exists( paste0(outdir,'/hormone_genes/',id_simple) ) ){  }else{dir.create(paste0(outdir,'/hormone_genes/',id_simple))}

    
    #res.marker <- list()
    if(id == 'Quality control'){marker.genes.sel = marker.genes[[id]]}else{
       marker.genes.sel <- marker.genes[[id]][marker.genes[[id]] %in% rownames(placenta)]
    }
    for (gene in marker.genes.sel){
        options(repr.plot.height=7.5,repr.plot.width=7.5)
        
        ##add short xaxis and yaxis with arrow
        umap.df <- Embeddings(placenta,reduction = 'umap_rotate')
        x.min <- min(umap.df[,1])
        y.min <- min(umap.df[,2])
        x.max <- max(umap.df[,1])
        y.max <- max(umap.df[,2])
        
        x.extend <- 0.2 * abs(x.max - x.min)
        y.extend <- 0.2 * abs(y.max - y.min)
        
        arrow <- data.frame(xstart=c(x.min,x.min),
                            xend=c(x.min+x.extend, x.min),
                            ystart=c(y.min,y.min),
                            yend=c(y.min,y.min+y.extend)
                           )
        
        res.p <- FeaturePlot(placenta, features = gene, reduction = "umap_rotate",
                             #pt.size = .5,
                             pt.size = 0.1,
                             slot = 'scale.data',
                             min.cutoff = 'q1',
                             max.cutoff = 'q99') +
                 scale_color_gradientn(colours = c('grey95',color_use)) + 
                 NoAxes() +
                 geom_segment(data = arrow, aes(x=xstart,y=ystart,xend=xend,yend=yend),col = 'grey50',size=.5, arrow = arrow(length = unit(0.02,'npc') )   ) +
                 ggtitle(paste0('marker gene of group ',id, ': ' ,gene) )+
                 theme(
                         legend.position = 'right'
#                         axis.text=element_blank(), 
#                         axis.title = element_text(size = 15, face = "bold"),
#                         axis.ticks = element_blank(),
#                         panel.grid.major = element_blank(), 
#                         panel.grid.minor = element_blank(),
#                         panel.background = element_rect(color="black", fill = NA,size=0.8),
#                 #         panel.background = element_rect(fill = "white", colour = "white", 
#                 #                 size = rel(1)),
#                         #panel.border = element_blank(),
#                         plot.title = element_text(size = 15, face = "bold"),
#                         #complete = TRUE
#                         plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
#                        )
#                   #labs(x = "UMAP1", y = "UMAP2")
                  )+
                  guides(fill = guide_legend(title="Expression"), guide_axis(angle=90))
        
        print(res.p)
        
        ggsave(
               filename = paste(paste0(outdir,'/marker_genes/',id_simple),
                                '/marker_gene.RNA.',gene,'.pdf',sep=''
                               ),
#               filename = paste(paste0(outdir,'/hormone_genes/',id_simple),
#                     '/hormone_gene.RNA.',gene,'.pdf',sep=''
#                    ),
               height=7.5,
               width=7.5,
               useDingbats=FALSE  
              )

  }
}  





#######add sample metadata###
placenta <- AddMetaData(placenta,metadata = Idents(placenta), col.name = 'cluster')
#placenta <- AddMetaData(placenta,metadata = rownames(placenta@meta.data), col.name = 'cellid') #for subset object by cell id

#placenta <- AddMetaData(placenta,metadata = 'female', col.name = 'sex')
#placenta <- AddMetaData(placenta,metadata = 'PLA-early3-RNA', col.name = 'sample')




#####get cluster.df.add####

cluster <- Idents(placenta)
umap <- placenta@reductions$umap@cell.embeddings
#umap <- placenta@reductions$umap_rotate@cell.embeddings

colnames(umap) <- c('UMAP_1','UMAP_2')

all.equal(names(cluster),rownames(umap)) #TRUE


cluster.df <- data.frame(cluster=cluster,umap)

metadata <- placenta@meta.data
all.equal(rownames(cluster.df),rownames(metadata))#TRUE

intersect(colnames(cluster.df), colnames(metadata) ) #cluster

grep('cluster',colnames(metadata) )
#49,50
#13,16

metadata[,c(49,50)] <- NULL
#metadata[,c(13,16)] <- NULL

all.equal(metadata,placenta@meta.data[,-c(49,50)])
#all.equal(metadata,placenta@meta.data[,-c(13,16)]) #TRUE


cluster.df.add <- cbind(cluster.df, metadata)

table(cluster.df.add$cluster)

   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
3311 2556 2284 2275 2134 1895 1777 1776 1693 1532 1301 1119  536  420  409  316 
  16   17   18   19 
 203  184  136  131


# cluster <- Idents(placenta)
# umap <- placenta@reductions$umap@cell.embeddings
# all.equal(names(cluster),rownames(umap)) #TRUE


# cluster.df <- data.frame(cluster=cluster,umap)

# metadata <- placenta@meta.data
# all.equal(rownames(cluster.df),rownames(metadata))#TRUE

# cluster.df.add <- cbind(cluster.df, metadata)


# ###filter low depth cells####
# table(cluster.df.add$cluster)
# cluster.df.add.filter <- subset(cluster.df.add,! cluster %in% c('10','17')  )
# cluster.df.add <- cluster.df.add.filter

# cluster.df.add$cluster <- droplevels(cluster.df.add$cluster)
# table(cluster.df.add$cluster)
# #  0    1    2    3    4    5    6    7    8    9   11   12   13   14   15   16 
# #1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74

# #11206 cells 



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


quickDimPlot_labelon <- function(data = cluster.df.add, feature = 'cluster', title= '8w all', xlim = NULL, ylim = NULL ,shrink.x = 1.0, shrink.y = 1.0){

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
  
  scale_colour_manual(values = color_good)  +
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




###quick looking and split looking

#quickDimPlot(data = cluster.df.add, feature = 'cluster', title= 'early combined')
quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', title= 'early combined')

par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
#for(i in c('8','7','0','5','2','1','3','15','16','4') ){
for(i in c('8','7','10','12','9','2','6','11','1','3','0','4','5','14','17','16','13','18','15') ){ #louvain 1.5
#for(i in c( '7','11','9','15','10','12','1','13','4','5','3','14','6','8'  )){ #leiden 1.4
#for(i in levels(cluster.df.add$cluster) ){
  dotDistri(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i)
}



saveRDS(cluster.df.add,'cluster.df.add.rds')
#write.table(file='cluster.df.add.txt',x=cluster.df.add,quote = FALSE,sep = '\t',row.names = TRUE,col.names = TRUE)


########save object before filtering###
saveRDS(placenta, "PLA-early-combine-RNA.rds")





####filter with customized qc threshold###

#modify doubletFinder columns
nameid = colnames(placenta@meta.data)[grepl("pANN", colnames(placenta@meta.data))]

##merge to one column
pANN.df <- placenta@meta.data[nameid] %>% mutate('pANN'= coalesce(pANN_0.25_0.09_950,pANN_0.25_0.09_591,pANN_0.25_0.09_761,pANN_0.25_0.09_810,pANN_0.25_0.09_686,pANN_0.25_0.09_928)) #%>% select(pANN_0.25_0.09_950,pANN)

placenta@meta.data[['pANN']] <- pANN.df$pANN
table(is.na(placenta@meta.data[['pANN']]))

FALSE 
47265


nameid = colnames(placenta@meta.data)[grepl("DF.classifications", colnames(placenta@meta.data))]
##merge to one column
DF.classifications.df <- placenta@meta.data[nameid] %>% mutate('DF.classifications'= coalesce(DF.classifications_0.25_0.09_950,DF.classifications_0.25_0.09_591,DF.classifications_0.25_0.09_761,DF.classifications_0.25_0.09_810,DF.classifications_0.25_0.09_686,DF.classifications_0.25_0.09_928)) #%>% select(pANN_0.25_0.09_950,pANN)


placenta@meta.data[['DF.classifications']] <- DF.classifications.df$DF.classifications
table(is.na(placenta@meta.data[['DF.classifications']]))
FALSE 
47265


options(repr.plot.height=12.5,repr.plot.width=12.5)
VlnPlot(placenta,features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo','percent.hb','pANN'),pt.size=0.1 ,group.by = 'orig.ident')


#placenta <- subset(placenta, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 ) #use 5000 instead of 5500, 
placenta.filter <- subset(placenta, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 & nCount_RNA >3500 & nCount_RNA < 5e4 & percent.mt < 5) #original filtering parameters, the same condition in scanpy

placenta.filter <- subset(placenta, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000 & nCount_RNA >3500 & nCount_RNA < 5e4 & percent.mt < 1.5 & percent.ribo < 2 & percent.hb < 2  &  pANN < 0.4 & DF.classifications == "Singlet") 


placenta.filter <- subset(placenta.filter, idents = levels(Idents(placenta.filter))[!levels(Idents(placenta.filter)) %in% c('13','5','14')]  )


##look at umap
options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(placenta.filter, reduction = "umap",label=TRUE,cols=color_good,label.size = 8,pt.size = 1)


#after filter 29132 x 27244  26139 x 4962    #24307 x 11560  #5646 #6471
options(repr.plot.height=7.5,repr.plot.width=12.5)
VlnPlot(placenta.filter,features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo','percent.hb','pANN'),pt.size=0.1,group.by = 'orig.ident' )





##rerun above pca->neighbor->umap codes


placenta <- placenta.filter




#####do we need to modify the placenta seurat object?yes####

#placenta@meta.data$seurat_clusters %in% levels(cluster.df.add$cluster)
#WhichCells(placenta,ident=levels(cluster.df.add$cluster) )


##rename and merge from cluster.df.add, louvain 1.5


cluster.raw <- as.character(cluster.df.add$cluster)
names(cluster.raw) <- rownames(cluster.df.add)

table(cluster.raw)
   0    1   10   11   12   13   14   15   16   17   18   19    2    3    4    5 
3311 2556 1301 1119  536  420  409  316  203  184  136  131 2284 2275 2134 1895 
   6    7    8    9 
1777 1776 1693 1532 

#    0    1   10   11   12   13   14   15   16   17   18   19    2    3    4    5 
# 3251 2556 1301 1097  536  420  409  316  203  184  136  131 2284 2275 2194 1895 
#    6    7    8    9 
# 1777 1776 1715 1532

map_cellname <- list(
    '8' = '8',
    '7' = '7',
    '10' = '10',
    '12' = '12',
    '9' = '9',
    '2' = '2',
    '6' = '11',
    '11' = '11',
    '1' = '1',
    '3' = '3',
    '0' = '4',
    '4' = '4',
    '5' = '5',
    '14' = '14',
    '17' = '17',
    '16' = '16',
    '13' = '13',
    '15' = '15',
    '18' = '18',
    '19' = '19'
)


for (i in seq_len(length(cluster.raw)) ){ 
    if(is.null(map_cellname[[cluster.raw[i]]])){message('error')}
    cluster.raw[i] <- map_cellname[[cluster.raw[i]]]
}

table(cluster.raw)
   1   10   11   12   13   14   15   16   17   18   19    2    3    4    5    7 
2556 1301 2896  536  420  409  316  203  184  136  131 2284 2275 5445 1895 1776 
   8    9 
1693 1532

   1   10   11   12   13   14   15   16   17   18   19    2    3    4    5    6 
2556 1301 1119  536  420  409  316  203  184  136  131 2284 2275 5445 1895 1777 
   7    8    9 
1776 1693 1532 


cluster.rename <- factor(cluster.raw,levels = as.character(sort(as.numeric(names(table(cluster.raw))))))


all.equal(names(cluster.rename),rownames(cluster.df.add)) #True

cluster.df.add$cluster_rename <- cluster.rename


quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster_rename', title= 'early combined')


cluster.df.add.filter<- subset(cluster.df.add,!(cluster_rename %in% c('15','19')) )

table(cluster.df.add.filter$cluster_rename)
  1    2    3    4    5    7    8    9   10   11   12   13   14   15   16   17 
2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409    0  203  184 
  18   19 
 136    0

#  1    2    3    4    5    6    7    8    9   10   12   13   14   15   16   17 
# 2556 2284 2275 5445 1895 1777 1776 2812 1532 1301  536  420  409    0  203  184 
#   18   19 
#  136    0 


cluster.df.add.filter$cluster_rename <- droplevels(cluster.df.add.filter$cluster_rename)

table(cluster.df.add.filter$cluster_rename)
   1    2    3    4    5    7    8    9   10   11   12   13   14   16   17   18 
2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409  203  184  136

  1    2    3    4    5    6    7    8    9   10   12   13   14   16   17   18 
2556 2284 2275 5445 1895 1777 1776 2812 1532 1301  536  420  409  203  184  136



cluster.raw <- as.character(cluster.df.add.filter$cluster_rename)
table(cluster.raw)
 1   10   11   12   13   14   16   17   18    2    3    4    5    7    8    9 
2556 1301 2896  536  420  409  203  184  136 2284 2275 5445 1895 1776 1693 1532 

#   1    2    3    4    5    6    7    8    9   10   12   13   14   16   17   18 
# 2556 2284 2275 5445 1895 1777 1776 2812 1532 1301  536  420  409  203  184  136 

as.character(sort(as.numeric(names(table(cluster.raw)))))

'1''2''3''4''5''7''8''9''10''11''12''13''14''16''17''18'

names(cluster.raw) <- rownames(cluster.df.add.filter)


##rename again by shifting and connecting
map_cellname <- list(
    '1' = '1',
    '2' = '2',
    '3' = '3',
    '4' = '4',
    '5' = '5',
    #'6' = '6',
    '7' = '6',
    '8' = '7',
    '9' = '8',
    '10' = '9',
    '11' = '10',
    '12' = '11',
    '13' = '12',
    '14' = '13',
    #'15' = '17',
    '16' = '14',
    '17' = '15',
    '18' = '16'
    #'19' = '18'
)


for (i in seq_len(length(cluster.raw)) ){ 
    if(is.null(map_cellname[[cluster.raw[i]]])){message('error')}
    cluster.raw[i] <- map_cellname[[cluster.raw[i]]]
}

cluster.rename <- factor(cluster.raw,levels = as.character(sort(as.numeric(names(table(cluster.raw))))))
table(cluster.rename)
  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409  203  184  136 


all.equal(names(cluster.rename), rownames(cluster.df.add.filter) ) #TRUE

cluster.df.add.filter$cluster_rename <- cluster.rename


quickDimPlot_labelon(data = cluster.df.add.filter, feature = 'cluster_rename', title= 'early combined')



cluster.df.add.bk <- cluster.df.add

cluster.df.add <- cluster.df.add.filter

table(cluster.df.add$cluster_rename)
 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409  203  184  136

 1    2    3    4    5    7    8    9   10   11   12   13   14   16   17   18 
2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409  203  184  136

quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster_rename', title= 'early combined')


#######modify seurat obj####

# ###rename ident from seurat obj
# placenta  <- RenameIdents(placenta,'0'='1', '1'='2','2'='3','3'='4','4'='5',
#            '5'='6','6'='7','7'='8','8'='9','9'='10',
#            '11'='11','12'='12','13'='13','14'='14','15'='15','16'='16' )


# ###


#placenta@meta.data

table(colnames(placenta) %in% rownames(cluster.df.add))
FALSE  TRUE 
  447 25541

# FALSE  TRUE 
#   354 11206 

#placenta_filter <- subset(placenta, subset = cellid %in% rownames(cluster.df.add) )
#11206

placenta_filter <- subset(placenta, cells = rownames(cluster.df.add) )
dim(placenta_filter)
#29132 x 25541

#24307 11206


all.equal (colnames(placenta_filter), rownames(cluster.df.add) )#TRUE
#all.equal (colnames(placenta_filter), colnames(placenta_filter1) )#TRUE

#all.equal(placenta_filter,placenta_filter1) #TRUE
#rm(placenta_filter1)


###rename

table(Idents(placenta_filter))

   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   16 
3311 2556 2284 2275 2134 1895 1777 1776 1693 1532 1301 1119  536  420  409  203 
  17   18 
 184  136 

#  0    1    2    3    4    5    6    7    8    9   11   12   13   14   15   16 
# 1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74



# placenta_filter  <- RenameIdents(placenta_filter,'0'='1', '1'='2','2'='3','3'='4','4'='5',
#            '5'='6','6'='7','7'='8','8'='9','9'='10',
#            '11'='11','12'='12','13'='13','14'='14','15'='15','16'='16' )


Idents(placenta_filter) <- cluster.df.add$cluster_rename


table(Idents(placenta_filter))
1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409  203  184  136

#  1    2    3    4    5    7    8    9   10   11   12   13   14   16   17   18 
# 2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409  203  184  136

#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74 

###

placenta <- placenta_filter


DimPlot(object = placenta, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() 


#####get cluster.df again and save
cluster <- Idents(placenta)
#umap <- placenta@reductions$umap@cell.embeddings
umap <- placenta@reductions$umap_rotate@cell.embeddings

colnames(umap) <- c('UMAP_1','UMAP_2')

all.equal(names(cluster),rownames(umap)) #TRUE


cluster.df <- data.frame(cluster=cluster,umap)

metadata <- placenta@meta.data
all.equal(rownames(cluster.df),rownames(metadata))#TRUE

intersect(colnames(cluster.df), colnames(metadata) ) #cluster

grep('cluster',colnames(metadata) )
#49,50
#13,16

metadata[,c(49,50)] <- NULL
#metadata[,c(13,16)] <- NULL

all.equal(metadata,placenta@meta.data[,-c(49,50)])
#all.equal(metadata,placenta@meta.data[,-c(13,16)]) #TRUE


cluster.df.add <- cbind(cluster.df, metadata)

table(cluster.df.add$cluster)

 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409  203  184  136 

#   1    2    3    4    5    7    8    9   10   11   12   13   14   16   17   18 
# 2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536  420  409  203  184  136 

# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 3311 2556 2284 2275 2134 1895 1777 1776 1693 1532 1301 1119  536  420  409  316 
#   16   17   18   19 
#  203  184  136  131

#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 4067 4032 3928 3705 3080 1809 1695 1361  535  420  404  361  203  183  136   47 
#   16 
#   22

#  0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 1133  738  645  543  481  411  382  204  166   89   45   44   43   38





##plot distribution of each clusters again

par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
#for(i in c('8','7','0','5','2','1','3','15','16','4') ){
for(i in c('8','7','10','12','9','2','11','1','3','4','5','14','17','13','16','18') ){ #louvain 1.5
#for(i in c( '7','11','9','15','10','12','1','13','4','5','3','14','6','8'  )){ #leiden 1.4
#for(i in levels(cluster.df.add$cluster) ){
  dotDistri(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i)
}




##filter for cstb obs###


cluster.df.add.cstb <- subset(cluster.df.add, cluster %in% c('7','9','6','11','8','2','10','5','4','1','3') )
#24189 x 57

cluster.df.add.cstb$cluster <- droplevels(cluster.df.add.cstb$cluster)

table(cluster.df.add.cstb$cluster)

   1    2    3    4    5    6    7    8    9   10   11 
2556 2284 2275 5445 1895 1776 1693 1532 1301 2896  536

##no need to rename

all.equal(colnames(placenta), rownames(cluster.df.add)) #TRUE


placenta.cstb <- subset(placenta, cells = rownames(cluster.df.add.cstb) )
all.equal(colnames(placenta.cstb), rownames(cluster.df.add.cstb))  #TRUE


#########filtering by distance q99 for cluster.df.add.cstb?##########


dotDist = function (cluster = NULL, id = NULL,center = centers, q = 0.98){
    colids = colnames(cluster) #must 'cluster','dim1','dim2' and rowname
    cluster.sel = cluster[ cluster[,1] == id,]
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    #color = 'red'
    #color = ifelse('cluster_sg' == id,'red','navy')
#     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
#     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
#     points(center[center$cluster == id,2:3],pch = 16, cex=1,col='black')
       
    dx <- cluster.sel[,2] - as.numeric(center[center$cluster == id,2])
    dy <- cluster.sel[,3] - as.numeric(center[center$cluster == id,3])
    d <- sqrt(dx**2 + dy**2)
    
    #options(repr.plot.height=15,repr.plot.width=15)
    res.h <- hist(d,breaks=100,main=paste("cluster ",id,sep=''),cex.main=2)
    #d.quantile <- quantile(d,prob = seq(0,1,0.1)) #min 0.1-0.9 max
    d.quantile <- quantile(d,prob = seq(0.9,1,0.01)) #0.9-1 total 11, 0.99 10
    d.qx <- quantile(d,prob = q) #select percentile
    abline(v=d.quantile,lty=2,lwd=1,col='red')
    text(x=d.quantile[10],y=0.9*(max(res.h$counts)),paste0('q99 is ',round(d.quantile[10],digits = 3),sep=''), pos =  4, adj=0.5,cex = 2,xpd = TRUE )
    text(x=d.qx,y=0.8*(max(res.h$counts)),paste0('q',100*q, ' is ',round(d.qx,digits = 3),sep=''), pos =  4, adj=0.5,cex = 2,xpd = TRUE )
    #flag <- d <= as.numeric(d.quantile['80%'])
    #flag <- d <=2
    return(d.qx)
    #return(paste("cluster ",id," ok",sep='') )
}


dotClean = function (cluster = NULL, id = NULL,center = centers,d.filter = d.filter){
    #return a flag.dist for each cluster
    colids = colnames(cluster) #must 'cluster','dim1','dim2'
    cluster.sel = cluster[ cluster[,1] == id,]
    d.filter.sel = d.filter[id]
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    color = 'red'
    cat(paste("cluster ",id," filtering dist",sep=''),'\n')
    #color = ifelse('cluster_sg' == id,'red','navy')
#     plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
#     points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col=color)
#     points(center[center$cluster == id,2:3],pch = 16, cex=1,col='black')
       
    dx <- cluster.sel[,2] - as.numeric(center[center$cluster == id,2])
    dy <- cluster.sel[,3] - as.numeric(center[center$cluster == id,3])
    d <- sqrt(dx**2 + dy**2)
    
    flag.dist <- d > d.filter.sel
    flag.dist.cl <- rownames(cluster)  %in% rownames(cluster.sel[flag.dist,])
    
    plot(cluster[,2],cluster[,3],pch = 16, type='p',col='grey',cex=0.5,xlab=colids[2],ylab=colids[3],main=paste(" cells cluster ",id," of ",colids[1],"\nn = ",n_sel,sep=''),cex.main = 2.25,xaxt = 'n' ) 
    points(cluster.sel[,2],cluster.sel[,3],pch = 16, cex=0.5,col='red')
    points(cluster.sel[!flag.dist,2],cluster.sel[!flag.dist,3],pch = 16, cex=0.5,col='blue') 
    points(center[center$cluster == id,2:3],pch = 16, cex=2.5,col='black')
    
    #return(paste("cluster ",id," ok",sep='') )
    return(flag.dist.cl)
    
}



##look for cluster dot distribution
#quickDimPlot(data = cluster.df.add, feature = 'cluster', title= 'early combined')
quickDimPlot_labelon(data = cluster.df.add.cstb, feature = 'cluster', title= 'early combined')



par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('7','9','6','11','8','4','2','10','5','1','3') ){
  dotDistri(cluster = cluster.df.add.cstb[,c('cluster','UMAP_1','UMAP_2')], id = i)
  
}


##look for distance distribution

centers <- cluster.df.add.cstb %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))


q.cutoff.list <- list()

par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('7','9','6','11','8','4','2','10','5','1','3') ){
  qx <- dotDist(cluster = cluster.df.add.cstb[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers, q = 0.98)
  q.cutoff.list[i] <- qx
}

# par(mfrow=c(2,3))
# options(repr.plot.height=10,repr.plot.width=15)
# for(i in c('12','9','13','15','11','14') ){
#   dotDist(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers)
# }


##get dist cutoff
d.filter <- unlist(q.cutoff.list)#c( '7'=,'9','6','11','8','4','2','10','5','1','3'   )

names(d.filter)
'7''9''6''11''8''4''2''10''5''1''3'

for(i in round(d.filter ,3)){cat(i,',',sep='')}
1.232,1.576,1.148,2.263,1.732,1.616,1.372,1.759,1.291,1.36,2.278 #q98
#1.28,1.686,1.228,2.378,2.433,1.715,1.477,1.873,1.395,1.461,3.01 #q99


# d.filter <- c( '1'=2.5, '2'=0,'3'=2,'4'=2,'5'=3,'6'=0,
#               '7'=2,'8'=2,'9'=0,'10'=0,'11'=0,'12'=0,
#               '13'=0,'14'=0,'15'=0 )


##do cleaning
par(mfrow=c(3,3))
res.flag.dist <- list()
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('7','9','6','11','8','4','2','10','5','1','3') ){
  res.flag.dist[[i]] <- dotClean(cluster = cluster.df.add.cstb[,c('cluster','UMAP_1','UMAP_2')], id = i,
           center = centers,d.filter = d.filter)
}


flag.dist.combine <- apply(do.call(cbind,res.flag.dist),1,any)
table(flag.dist.combine)
#flag.dist.combine
FALSE  TRUE 
23702   487

# FALSE  TRUE 
# 13722   634 

cluster.df.add.cstb.sel <- cluster.df.add.cstb[!flag.dist.combine,]
cluster.df.add.cstb.rm <- cluster.df.add.cstb[flag.dist.combine,]

###quick plot cluster distribution again
#quickDimPlot(data = cluster.df.add, feature = 'cluster', title= 'early combined')
quickDimPlot_labelon(data = cluster.df.add.cstb.sel, feature = 'cluster', title= 'early combined', shrink.x = 2, shrink.y = 1.05)
quickDimPlot_labelon(data = cluster.df.add.cstb.rm, feature = 'cluster', title= 'early combined')



par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('7','9','6','11','8','4','2','10','5','1','3') ){
  dotDistri(cluster = cluster.df.add.cstb.sel[,c('cluster','UMAP_1','UMAP_2')], id = i)
  
}





######filter obj##


all.equal(rownames(cluster.df.add.cstb),  colnames(placenta) ) #TRUE

cluster.df.add.cstb.bk <- cluster.df.add.cstb
cluster.df.add.cstb <- cluster.df.add.cstb.sel


placenta.sel <- subset(placenta,cells = rownames(cluster.df.add.cstb.sel))

options(repr.plot.width = 7.5, repr.plot.height=7.5)
DimPlot(object = placenta.sel, label = TRUE,cols=c(color_good),pt.size = 0.5,label.size = 10,reduction = "umap_rotate") + NoLegend() 



placenta <- placenta.sel

cluster.df.add <- cluster.df.add.cstb 


##modify placenta obj meta and Idents


all.equal(colnames(placenta),rownames(cluster.df.add)) #TRUE

placenta@meta.data$cluster <- cluster.df.add$cluster
Idents(placenta) <- cluster.df.add$cluster

#save rds in snapshot dir (see the bottom of this script)




########save object after filtering and fixing###

#final
saveRDS(cluster.df.add,'cluster.df.add.final.rds')
write.table(file='cluster.df.add.final.txt',x=cluster.df.add,quote = FALSE,sep = '\t',row.names = TRUE,col.names = TRUE)

cluster.df.add <- readRDS('cluster.df.add.final.rds')


##final final
saveRDS(cluster.df.add,'cluster.df.add.final.final.rds')
write.table(file='cluster.df.add.final.final.txt',x=cluster.df.add,quote = FALSE,sep = '\t',row.names = TRUE,col.names = TRUE)

cluster.df.add <- readRDS('cluster.df.add.final.final.rds')


##cstb
saveRDS(cluster.df.add.cstb,'cluster.df.add.cstb.rds')
write.table(file='cluster.df.add.cstb.txt',x=cluster.df.add.cstb,quote = FALSE,sep = '\t',row.names = TRUE,col.names = TRUE)



##final obj
saveRDS(placenta, "PLA-early-combine-RNA.final.rds")
placenta <- readRDS("PLA-early-combine-RNA.final.rds")

##final final obj
saveRDS(placenta, "PLA-early-combine-RNA.final.final.rds")
placenta <- readRDS("PLA-early-combine-RNA.final.final.rds")


##cstb obj
saveRDS(placenta.cstb, "PLA-early-combine-RNA.cstb.rds")
#placenta.cstb <- readRDS("PLA-early-combine-RNA.cstb.rds")



##for plotting convinency##
cluster.df.add <- cluster.df.add.cstb
placenta <- placenta.cstb


####find out which clusters to filter snRNA donor1 low quality cells

# cluster.donor1 <- read.table("../../placenta_10X/03.seurate/PLA-8w-RNASEQ-1/placenta.PLA-8w-RNA-1.umap.cl.txt",header = TRUE,sep = '\t',row.names = 1)

# cluster.donor1$cluster <- factor(cluster.donor1$cluster,levels=sort(unique(cluster.donor1$cluster)) )

# cluster.donor1.cents = cluster.donor1 %>% dplyr::group_by(cluster) %>%
#                    dplyr::summarise (meanx = mean(UMAP_1),
# 									 meany=mean(UMAP_2))


# options(repr.plot.width=7,repr.plot.height=7)
# ggplot(data=cluster.donor1,aes(x=UMAP_1,y=UMAP_2,col=cluster) ) +
#   geom_point(size=1) +
#   scale_color_manual(values = color_good) +
#   #xlim(-10,10) +
#   #ylim(-10,10) +
#   annotate(geom = 'text',x=cluster.donor1.cents$meanx,
# 		   y=cluster.donor1.cents$meany,
# 		  label=cluster.donor1.cents$cluster,
# 		  size =10,color='black') +
#   theme_classic() +
#   guides( col = guide_legend(override.aes = list(size = 6)) )

# ##filter cluster 1 

# id1 <- rownames(subset( cluster.donor1, cluster == '1'  ))

# cluster.df.add.filter <- cluster.df.add[!rownames(cluster.df.add) %in% id1,]
# #15107

# saveRDS(cluster.df.add.filter,'cluster.df.add.filter.rds')

# cluster.df.add<- cluster.df.add.filter
# rm(cluster.df.add.filter)




###############customized way to plot umap-cluster with text halo##########
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
                                  cluster=as.character(unlist(centers[i,'cluster'])),
                                  x=centers[i,'x'] + cos(j)*xo, 
                                  y=centers[i,'y'] + sin(j)*yo
                                 )
                       )
      }
}

##for x y shrink
xmin = min(cluster.df.add[,'UMAP_1']) 
xmax = max(cluster.df.add[,'UMAP_1'])
ymin = min(cluster.df.add[,'UMAP_2']) 
ymax = max(cluster.df.add[,'UMAP_2']) 

shrink.x = 2
shrink.y = 1.05

################the UMAP plot with annotation##############
#label right
options(repr.plot.height=7.5,repr.plot.width=7.5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = .5,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
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
  ggtitle(paste(sample, "six donors, total cells:",nrow(cluster.df.add),  sep=" ") ) +
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

ggsave(filename = "pdfs/PLA-early-combined-RNA-UMAP.pdf",height=5,width=5.5,useDingbats=FALSE)


##label on cluster
options(repr.plot.height=6.5,repr.plot.width=6.5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = map_cellcolor_rna) + #color_good)  +
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
        panel.background = element_blank(),#element_rect(color="black", fill = NA,size=0),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "six donors, total cells:",nrow(cluster.df.add),  sep=" ") ) +
  geom_text(data = centers_shift, #the halo
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "white", 
            size = 6.5) +
  geom_text(data = centers, 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            size = 6) +
  ##guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  ylim(shrink.y*ymin,shrink.y*ymax) +
  xlim(shrink.x*xmin,shrink.x*xmax) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

#print(res.p, vp=viewport(angle=-90)) #will not rotate text

ggsave(filename = "pdfs/UMAP/PLA-early-combined-RNA-UMAP.labelon.hotcolor.pdf",height=6.5,width=6.5,useDingbats=FALSE)
#ggsave(filename = "pdfs/UMAP/PLA-early-combined-RNA-UMAP.labelon.coldcolor.pdf",height=6.5,width=6.5,useDingbats=FALSE)

###############


color_donor <- list('PLA-early1-RNA'='#001C7F',
                    'PLA-early2-RNA'='#017517',
                    'PLA-early3-RNA'='#8C0900',
                    'PLA-early4-RNA'='#7600A1',
                    'PLA-early5-RNA'='#B8860B',
                    'PLA-early6-RNA'='#006374'
                   )

#PALETTES['sns_dark']
#'#001C7F', '#017517', '#8C0900', '#7600A1', '#B8860B', '#006374'


##by donors
#options(repr.plot.height=7.5,repr.plot.width=7)
options(repr.plot.height=7.5,repr.plot.width=25) #for facet_grid
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= sample  )) +
  geom_point(size = .1,show.legend = TRUE,alpha= 0.25 ) +
  scale_colour_manual(values = color_donor )+#c('red','navy','orange','green','purple','black'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  facet_grid(cols = vars(sample)) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Source of donor ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  ylim(shrink.y*ymin,shrink.y*ymax) +
  xlim(shrink.x*xmin,shrink.x*xmax) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

#ggsave(filename = "pdfs/QC/PLA-early-combined-RNA-source-of-donor.pdf",height=7.5,width=7)
ggsave(filename = "pdfs/QC/PLA-early-combined-RNA-source-of-donor.split.pdf",height=7.5,width=25)

####

##by sex
color_sex <- list('female'= 'darkred', 'male' = 'darkblue')

options(repr.plot.height=7.5,repr.plot.width=7)
#options(repr.plot.height=7.5,repr.plot.width=15) #for facet_grid
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= sex  )) +
  geom_point(size = .1,show.legend = TRUE,alpha= 0.3 ) +
  scale_colour_manual(values = color_sex)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ##facet_grid(cols = vars(sex)) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Source of sex ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  ylim(shrink.y*ymin,shrink.y*ymax) +
  xlim(shrink.x*xmin,shrink.x*xmax) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/QC/PLA-early-combined-RNA-source-or-sex.pdf",height=7.5,width=7)
#ggsave(filename = "pdfs/QC/PLA-early-combined-RNA-source-or-sex.split.pdf",height=7.5,width=15)

###

###by depth
options(repr.plot.height=7.5,repr.plot.width=7)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= log10(nCount_RNA) )) +
  geom_point(size = .2,show.legend = TRUE,alpha= .5 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  ##scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  scale_colour_gradientn(colors = rev(color_cellranger) )  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Sequence depth",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  ylim(shrink.y*ymin,shrink.y*ymax) +
  xlim(shrink.x*xmin,shrink.x*xmax) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/QC/PLA-early-combined-RNA-depth.pdf",height=7.5,width=7)


#by gene captured
options(repr.plot.height=7.5,repr.plot.width=7)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= nFeature_RNA )) +
  geom_point(size = 0.5,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Gene captured",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  ylim(shrink.y*ymin,shrink.y*ymax) +
  xlim(shrink.x*xmin,shrink.x*xmax) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/QC/PLA-early-combined-RNA-gene_captured.pdf",height=7.5,width=7)

##by cell cycling stage
options(repr.plot.height=7.5,repr.plot.width=7)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= Phase  )) +
  geom_point(size = 0.2,show.legend = TRUE,alpha= 0.5 ) +
  scale_colour_manual(values = c('#282F76','brown','darkgreen'))  +
  #scale_colour_manual(values = c('red','navy','orange'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "cell cycling phase",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  ylim(shrink.y*ymin,shrink.y*ymax) +
  xlim(shrink.x*xmin,shrink.x*xmax) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/QC/PLA-early-combined-RNA-cell_cycling.pdf",height=7.5,width=7)

##by cell cycling stage: plot layers of Phase score
cluster.df.add.G1 <- subset(cluster.df.add,Phase == 'G1') #6246
cluster.df.add.G2M <- subset(cluster.df.add,Phase == 'G2M') #1950
cluster.df.add.S <- subset(cluster.df.add,Phase == 'S') #3010

options(repr.plot.height=5,repr.plot.width=5)
ggplot(cluster.df.add.G1,aes(x=UMAP_1,y=UMAP_2,col= Phase  )) +
  #geom_point(size = 0.1,alpha= 0.5,col='#282F76' ) +
  ##geom_point(data=cluster.df.add.G1,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col='#282F76' ) +
  ##geom_point(data=cluster.df.add.G2M,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col='brown' ) +
  geom_point(data=cluster.df.add.S,aes(x=UMAP_1,y=UMAP_2,col= Phase),size = 0.1,alpha= 0.5,col='darkgreen' ) +
  #scale_colour_manual(values = c('#282F76','brown','darkgreen'))  +
  #scale_colour_manual(values = c('red','navy','orange'))  +
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
  ggtitle(paste(sample, "cell cycling phase S",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")





######plot cluster count barplot, use the same colors as UMAP
# colPanel = SnapATAC:::createColorPanel(length(unique(x.after.sp@cluster)))
# color = scales::alpha(colPanel[factor(unique(x.after.sp@cluster))], 0.8)

color = color_good

options(repr.plot.width=2,repr.plot.height=7)
barplot(table(cluster.df.add$cluster),horiz = TRUE,col = color )

##stb percentage
stb_sum <- sum(table(cluster.df.add$cluster)[c('1','2','3','4','6','7')]) #7322 #5777 #6779
stb_sum/sum(table(cluster.df.add$cluster)) #11206 #15107
#65.3% #51.6%  44.9%


######stat d1 d2 and cluster cell number correlation######
res.stat <- table(cluster.df.add$cluster,cluster.df.add$sample)

     PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA
  1             574            321            359            373            294
  2             486            340            377            357            339
  3             323            295            219            300            358
  4             990            687            811            674            765
  5             339            244            298            261            241
  6             499            379            546            266             34
  7             464            282            462            370             59
  8             319            287            231            260            189
  9             335            250            505            132             41
  10            496            435            465            369            446
  11            134             90            145             74             33
    
     PLA-early6-RNA
  1             583
  2             339
  3             734
  4            1409
  5             474
  6              16
  7              22
  8             215
  9              12
  10            627
  11             49

      D1  D2
  1  992 728
  2  827 601
  3  712 513
  4  669 424
  5  645 438
  6  529 434
  7  473 420
  8  544 309
  9  406 250
  10 279  86
  11 185 103
  12 111  78
  13 124  49
  14  77  46
  15  53  27
  16  32  42

#       D1  D2
#   0  992 728
#   1  827 601
#   2  712 513
#   3  669 424
#   4  645 438
#   5  529 434
#   6  473 420
#   7  544 309
#   8  406 250
#   9  279  86
#   11 185 103
#   12 111  78
#   13 124  49
#   14  77  46
#   15  53  27
#   16  32  42
res.cor <- cor(res.stat)#0.978
res.stat[,1] = -1 * res.stat[,1]

res.stat.df <- as.data.frame(res.stat)
#res.stat.df <- reshape2::melt(res.stat)
colnames(res.stat.df) <- c('cluster','sample','count')

res.stat.df$count <- log2(res.stat.df$count)

# ##plot horizonal barplot  #  '#94C6DD', '#1273AE'  vs '#F3A585','#C80927'
# options(repr.plot.height=5.5,repr.plot.width=4)
# ggplot(res.stat.df, aes(fill=sample, y=count, x=cluster )) + 
#     #geom_hline(yintercept = c(0,25,50,75),linetype='solid',size=.3,col='black') +
#     #geom_hline(yintercept = c(0),linetype='solid',size=1,col='black') +
#     #geom_bar(position="stack", stat="identity",alpha=1,width = 0.5) +  
#     geom_bar(position=position_stack(reverse=TRUE), stat="identity",width = 0.5) +  
#     #xlim(100,0) +
#     #scale_x_continuous(breaks = seq(100,0,-25),labels=seq(100,0,-25)) +
#     #scale_y_reverse() +
#     #annotate('text',x = 16, y = 1000,label = paste('spearman cor=',cor(res.stat)[1,2]) ) +
#     coord_flip() +
#     #scale_fill_viridis(discrete = T,option = "E") +
#     scale_fill_manual(values = c('#F3A585', '#C80927'),labels=c('D1','D2'),name='sample' ) +
#     #ggtitle("cell number count") +
#     labs(title = "cell number count", subtitle=paste('spearman cor=',round(res.cor[1,2],digits=2)) ) +
#     theme_ipsum(base_family = 'sans') + #to avoid Arial narrow problem in ggsave
#     ylab("count")

# ggsave(filename="donor1.vs.donor2.cluster.cor.pdf",width = 4, height = 5.5)


options(repr.plot.height=10.5,repr.plot.width=10)
ggplot(res.stat.df, aes(fill=count, y=cluster, x=sample )) + 
    #geom_hline(yintercept = c(0,25,50,75),linetype='solid',size=.3,col='black') +
    #geom_hline(yintercept = c(0),linetype='solid',size=1,col='black') +
    #geom_bar(position="stack", stat="identity",alpha=1,width = 0.5) +  
    #geom_bar(position=position_stack(reverse=TRUE), stat="identity",width = 0.5) +  
    geom_tile() +
    geom_text(aes(label=count),size = 6) +
    #xlim(100,0) +
    #scale_x_continuous(breaks = seq(100,0,-25),labels=seq(100,0,-25)) +
    #scale_y_reverse() +
    #annotate('text',x = 16, y = 1000,label = paste('spearman cor=',cor(res.stat)[1,2]) ) +
    coord_flip() +
    scale_fill_gradientn(colours = c('blue','white','red') ) +
    #scale_fill_viridis(discrete = T,option = "E") +
    #scale_fill_manual(values = c('#F3A585', '#C80927'),labels=c('D1','D2'),name='sample' ) +
    #ggtitle("cell number count") +
    ##labs(title = "cell number count", subtitle=paste('spearman cor=',round(res.cor[1,2],digits=2)) ) +
    theme_ipsum(base_family = 'sans') + #to avoid Arial narrow problem in ggsave
    ylab("count")+
    theme(legend.title = element_blank(),
            axis.text.x = element_text(angle=30,hjust=1,vjust=1.0),
            axis.text.y = element_text(size = 12))

ggsave(filename="pdfs/QC/sample.cluster.count.pdf",width = 10, height = 10.5)




###########marker gene plotting###########



marker.genes.villi = list(
    'CTB'=c("DNMT1", "CDH1", "PAGE4",'TP63','MKI67','PCNA'),
    'STB'=c('FLT1','LVRN' ,'LEP', 'ENG','PAPPA','CSHL1','CSH1','PSG8','CGA','LAMA3','BMP1'), #c("FLT1", "CSHL1", "PSG8", "PAPPA",'GCM1','CGA','CGB8'),
    'Fusion_CTB'=c("ERVFRD-1", "ANXA1",'ESRRG'),
    'EVT'=c("LAIR2", "PLAC8",'HLA-G','MMP2'),
    'STR'=c('VIM','HLA-A','HLA-DPA1','HLA-DPB1','DLK1','HIVEP3')
    #'MesSTR' = c('THY'),
    #'VascularEndo' = c('PECAM1'),
    #'Macro' = c('CD68','CD14')
  )



marker.genes <- marker.genes.villi

marker.genes <- marker.genes.villi_STB


marker.genes.full <- list( #from seurat code
    #'Quality control' = c('mitochondrial','passed_filters', 'logUMI', 'promoter_ratio', 'peak_ratio', 'tsse', "n_fragment","frac_mito","frac_dup",'XIST'),
    'Trophoblast' = c('KRT7', 'GATA3', 'TFAP2A'),
    'CTB' = c('DNMT1', 'CDH1', 'PPARG', 'TEAD4', 'TEAD3', 'MKI67', 'TP53','TP63', 'TP73', 'BCAM'),
    'CTB fusion' = c('ERVFRD-1', 'GCM1', 'OVOL1','PPARD'), 
    'STB nascent' = c('SH3TC2', 'BACE2', 'ESRRG'), 
    'STB general' = c('PSG8', 'CGA', 'PSG2', 'PSG5', 'LEP'), 
    'STB PAPPA' = c('PAPPA', 'ADAMTSL1','ADAMTS6', 'GH2', 'GHR','JAK1', 'JAK2', 'LAMA3', 'AR', 'VDR', 'CSHL1', 'CSH1', 'CSH2', 'STAT5A', 'STAT5B', 'STAT4','FOS', 'FOSB', 'JUNB', 'JUN'), 
    'STB FLT1' = c('FLT1', 'ENG', 'ANGPTL4', 'FSTL3', 'INHBA', 'INHA','MYCN', 'POU2F3', 'LVRN', 'TGFB1', 'FOSL2', 'JUND'), 
    'STB apoptosis' = c('DDX60', 'DDX58', 'MAP4K4', 'SPATA5', 'GDF15', 'CROT', 'CDKN1A', 'ADCY5'),
    'EVT' = c('HLA-G', 'LAIR2', 'PLAC8', 'MMP2'), 
    'STR general' = c('VIM', 'DLK1'), 
    'Vascular Endothelial Cell' = c('PECAM1'), 
    'STR' = c('HIVEP3', 'HLA-A', 'HLA-DPA1', 'HLA-DPB1'), 
    'Mesenchymal STR'= c('THY1'), 
    'Hofbauer Cell'= c('CD68', 'CD14'), 
    'Red blood' = c('HBA1', 'HBZ')
)

marker.genes <- unlist(marker.genes.full)

marker.genes <- list(marker_gene = c("DNMT1","ERVFRD-1",'PSG8','CGA','SH3TC2','LEP','PAPPA','FLT1'))



###########customized FeaturePlot for marker genes ####


###method 0: quick plot
res.marker <- list()
for (gene in marker.genes){
    options(repr.plot.height=5,repr.plot.width=5.5)
    p<- FeaturePlot(placenta, features = gene, reduction = "umap",slot = 'data',cols = c("lightgrey", "darkred") )+ #c("lightgrey", "#ff0000", "#00ff00") c("lightgrey", "darkred")
      #scale_color_gradientn(colours = color_ga)+
      #scale_color_gradientn(colours = color_peak)+
    
      #scale_color_gradientn(colours = color_gradient_my)+
      #scale_color_gradientn(colours = color_pyreds)+
      theme(
            legend.position = 'right',
            axis.text=element_blank(), 
            axis.title = element_text(size = 15, face = "bold"),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(color="black", fill = NA,size=0.8),
    #         panel.background = element_rect(fill = "white", colour = "white", 
    #                 size = rel(1)),
            #panel.border = element_blank(),
            plot.title = element_text(size = 15, face = "bold"),
            #complete = TRUE
            plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
           ) 
    print(p)
    #ggsave(filename = paste('marker_gene.red_blue',gene,'.pdf',sep=''),height=5,width=5.5,useDingbats=FALSE  )
     res.marker[[gene]] <- p

}


##arrange plot by patchwork
options(repr.plot.height=13.5,repr.plot.width=18)
res.marker[['PSG8']] + res.marker[['DNMT1']] + res.marker[['ERVFRD-1']] + res.marker[['PECAM1']] +
res.marker[['CSHL1']] + res.marker[['CDH1']] + res.marker[['LAIR2']] + res.marker[['VIM']] +
res.marker[['FLT1']] + res.marker[['MKI67']] + res.marker[['PLAC8']] + res.marker[['CD14']] +
plot_layout(ncol=4,nrow=3)
#print(res.marker[['FLT1']], vp=viewport(angle=-185))
ggsave(filename = 'marker.gene.umap.pdf',height=13.5,width=18,useDingbats=FALSE) #can save grid multiple plots



#####method 1: use seurat FeaturePlot#####
res.wrap <- list()
for(id in names(marker.genes) ){
    
    cat('do for group ',id,'\n',sep='')
    
#    obj <- placenta.list.obj[[id]]
    

    ###########customized FeaturePlot for marker genes ####
    res.marker <- list()
    marker.genes.sel <- marker.genes[[id]][marker.genes[[id]] %in% rownames(placenta)]
    
    for (gene in marker.genes.sel){
        options(repr.plot.height=7.5,repr.plot.width=8)
        p<- FeaturePlot(placenta, features = gene, reduction = "umap_rotate",slot = 'scale.data',
                        #cols = as.character(color_use), #3+ colors will ignored the >3  colors!!
                        #cols = c("lightgrey", "darkred"),
                        ##cols = c("lightgrey", "darkgreen"), #use this?
                        #cols = c('darkgreen',"lightgrey", "darkred"),
                        #cols = c("lightgrey", "darkred"),
                        ##cols = color_use,
                        #cols = c("lightgrey", "purple"),
                        pt.size = .1,
                        min.cutoff = "q10",
                        max.cutoff = "q99"
                       )+ #c("lightgrey", "#ff0000", "#00ff00") c("lightgrey", "darkred")
        #p<- FeaturePlot(placenta_filter, features = gene, reduction = "umap",slot = 'data',cols = c("lightgrey", "darkred") )+ #c("lightgrey", "#ff0000", "#00ff00") c("lightgrey", "darkred")
          #scale_color_gradientn(colours = color_ga)+
          
          ##scale_color_gradientn(colours = color_set)+

          #scale_color_gradientn(colours = color_gradient_my)+
          scale_color_gradientn(colours = color_use)+
          #scale_color_gradientn(colours = color_pyreds)+
          #ggtitle(paste0('marker gene group ',id)) +
          theme(
                legend.position = 'right',
                axis.line = element_blank(),
                axis.text=element_blank(), 
                axis.title = element_text(size = 15, face = "bold"),
                axis.ticks = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),#element_rect(color="black", fill = NA,size=0.8),
        #         panel.background = element_rect(fill = "white", colour = "white", 
        #                 size = rel(1)),
                #panel.border = element_blank(),
                plot.title = element_text(size = 15, face = "bold"),
                #complete = TRUE
                plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
               ) +
          labs(x = "UMAP1", y = "UMAP2")
        #print(p)
        #ggsave(filename = paste('pdfs/marker_gene.RNA.red_blue',gene,'.pdf',sep=''),height=5,width=5.5,useDingbats=FALSE  )
         res.marker[[gene]] <- p

    }
    
    
    ##arrange plot by patchwork
    
    nplots <- length(res.marker)
    nrow <- ceiling(nplots / 2)
    message(paste0( 'nplots: ', nplots, "  nrow: " ,nrow,'\n')  )
    ##options(repr.plot.height=nrow * 5.5,repr.plot.width=18)
    options(repr.plot.height=nrow * 5.5,repr.plot.width=12)
#     res.marker[['PSG8']] + res.marker[['DNMT1']] + res.marker[['ERVFRD-1']] + res.marker[['PECAM1']] +
#     res.marker[['CSHL1']] + res.marker[['CDH1']] + res.marker[['LAIR2']] + res.marker[['VIM']] +
#     res.marker[['FLT1']] + res.marker[['MKI67']] + res.marker[['PLAC8']] + res.marker[['CD14']] +
#     plot_layout(ncol=4,nrow=3)
    #print(res.marker[['FLT1']], vp=viewport(angle=-185))
    res.wrap[[id]] <- wrap_plots(res.marker,byrow = TRUE,ncol = 4) + plot_annotation(title = paste0('marker gene group ',id),theme = theme(plot.title = element_text(size = 16)))
    ##options(repr.plot.height=nrow * 5.5,repr.plot.width=12)
    ##print(res.wrap)
    #ggsave(filename = paste0('pdfs/marker.gene.plots/',id,'_',marker.genes.name,'.marker.gene.umap.pdf'),height=nrow * 5.5,width=18,useDingbats=FALSE) #can save grid multiple plots
    
    
}


for(i in names(res.wrap)){
    print(res.wrap[[i]])
    
}

##plot all grid
options(repr.plot.height=10,repr.plot.width=20)
print(res.wrap[[1]])



#####method 2: customized ggplot from scale.data mat (expr.mat.z), use this

marker.genes <-  c("DNMT1","ERVFRD-1",'PSG8','CGA','SH3TC2','LEP','PAPPA','FLT1')

exprMat.z <- GetAssayData(placenta,slot = 'scale.data')
exprMat.data <- GetAssayData(placenta,slot = 'data')
exprMat.count <- GetAssayData(placenta,slot = 'counts')


saveRDS(exprMat.data,'exprMat.data.rds')
saveRDS(exprMat.count,'exprMat.count.rds')


# ###i. smooth by sliding window with ArchR function?? (with data.table::frollingmean, k=5, apply to row of all cell barcode, along trajectory)
# exprMat.z.smooth <- as.matrix(t(apply(exprMat.z, 1, function(x) {ArchR:::.centerRollMean(x, 
#     k = 5)}  )   ))
# colnames(exprMat.z.smooth) <- paste0(colnames(mat.all.sel))


##ii. use Magic package to get imputation
library(Rmagic) #v2.0.3.99 from git pip install magic-impute==3.0.0 #v1.4.0 need python magic-imputate, install magic-imputate package in r-reticular conda env:pip install magic-impute==1.4.0
exprMat.z.magic <- Rmagic::magic(t(exprMat.z),
                                 genes = NULL,#marker.genes,
                                 n.jobs=1, #broken with > 1!!
                                ) #need gene as column mat #batch effect?

#start "2023-02-02 5:18 CST"

timetag <- Sys.time()

timetag
"2023-02-02 05:23:17 CST" #5 min

#"2023-02-02 04:57:08 CST" #9 min, quick in fact


saveRDS(exprMat.z.magic,'exprMat.z.magic.rds')



##iii. use figR smoothScoreNN in FigR_shareseq_tutorial.html (use this)
library(FigR)
library(doParallel)

set.seed(123)
cellkNN <- FNN::get.knn(Embeddings(placenta,reduction = 'pca'),  k= 30)$nn.index
#23702 × 30

rownames(cellkNN) <- colnames(exprMat.z)

exprMat.z.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = exprMat.z,nCores = 10) #need gene x cell mat
#29132 x 23702



#saveRDS(exprMat.z.s,'exprMat.z.s.rds') #big! 6.3g
exprMat.z.s <- readRDS('exprMat.z.s.rds')


##plotting with imputated exprMat.z

exprMat.z <- exprMat.z.s
#exprMat.z <- t(exprMat.z.magic$result)


marker.genes <- marker.genes[ marker.genes %in% rownames(exprMat.z)  ]

marker.mat <- t(as.data.frame(as(exprMat.z[marker.genes,],'matrix')))
#rownames(marker.mat) <- rownames()

all.equal(rownames(marker.mat),rownames(cluster.df.add)) #TRUE

marker.mat.df <- cbind(cluster.df.add[,c('cluster','UMAP_1','UMAP_2')],marker.mat )

##cutoff quantiles
cutoff.thresh = list(
    'ERVFRD-1'=list(0.3,1),
    'PAPPA'=list(0.3,.99),
    'LAMA3'=list(0,1),
    'STAT5A'=list(0,.99),
    'ESRRG'=list(0,.99),
    'FLT1'=list(0.6,.99),
    'ENG'=list(0,1),
    'FOSL1'=list(0,.99),
    'JUNB'=list(0,.99),
    'FOS'=list(0,.99),
    'HLA-G'=list(0,.99),
    'PLAC8'=list(0,.99),
    
    'DNMT1'=list(0.3,.99),
    'PSG8'=list(0.1,.99),#'PSG8'=list(0.3,.99),
    'CGA'=list(0.1,1),#list(0.5,.99),
    'SH3TC2'=list(0.5,1),
    'LEP'=list(0.5,.99)
)

marker.mat.df.cutoff <- marker.mat.df[,1:3]
for(gene in marker.genes){
    #cat('gene ',gene)
    feature.value <- marker.mat.df[,gene]
    #cutoff <- quantile(feature.value,probs = c(0,1) )
    low.q <- cutoff.thresh[[gene]][[1]]
    high.q <- cutoff.thresh[[gene]][[2]]
    if(is.null(low.q)){ low.q = 0    }
    if(is.null(high.q)){ high.q = 1    }
    cutoff <- quantile(feature.value,probs = c(low.q,high.q)) 
    cutoff.low <- cutoff[1]
    cutoff.high <- cutoff[2]
    feature.value[feature.value < cutoff.low] <- cutoff.low
    feature.value[feature.value > cutoff.high] <- cutoff.high
    marker.mat.df.cutoff <- cbind(marker.mat.df.cutoff,feature.value)
}
colnames(marker.mat.df.cutoff)[4:ncol(marker.mat.df.cutoff)] <- marker.genes


######

xmin = min(marker.mat.df.cutoff[,'UMAP_1']) 
xmax = max(marker.mat.df.cutoff[,'UMAP_1'])
ymin = min(marker.mat.df.cutoff[,'UMAP_2']) 
ymax = max(marker.mat.df.cutoff[,'UMAP_2']) 

shrink.x = 1.2
shrink.y = 1.05 

#options(repr.plot.height=5.8,repr.plot.width=5)
res.marker <- list()
for(gene in marker.genes){
    
    ##to label title
    low.q <- cutoff.thresh[[gene]][[1]]
    high.q <- cutoff.thresh[[gene]][[2]]
    if(is.null(low.q)){ low.q = 0    }
    if(is.null(high.q)){ high.q = 1    }
    
    p <- ggplot(marker.mat.df.cutoff,aes_string(x="UMAP_1",y="UMAP_2",col= paste0("`",gene,"`") )) +
    #p <- ggplot(marker.mat.df,aes_string(x="UMAP_1",y="UMAP_2",col= paste0("`",gene,"`") )) +
      geom_point(size = .1,show.legend = TRUE,alpha= 1 ) +
      #scale_colour_manual(values = c('red','navy'))  +
      ##scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
      #scale_colour_gradientn(colors = color_tfdev )  +
      scale_colour_gradientn(colors = color_use )  +
      #theme_classic() +
      #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
      theme(
            #legend.position = 'none',
            legend.position = 'right',
            axis.text=element_blank(), 
            axis.title = element_text(size = 15, face = "bold"),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),#element_rect(color="black", fill = NA,size=1),
            plot.title = element_text(size = 15, face = "bold"),
            #complete = TRUE
            plot.margin = unit(c(1,1,1,1), "lines")
           )+
      ggtitle(paste0(gene, ', cutoff: ' ,low.q, ', ' ,high.q )) +
      #ggtitle(gene) +
      #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
      #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
      #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
      ylim(shrink.y*ymin,shrink.y*ymax) +
      xlim(shrink.x*xmin,shrink.x*xmax) +
      #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
      labs(x = "UMAP1", y = "UMAP2")
      #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
    #ggsave(filename = "PLA-8w-RNA-depth.pdf",height=5.5,width=5)
    res.marker[[gene]] <- p
    #print(p)
}



##arrange plot by patchwork
options(repr.plot.height=10.5,repr.plot.width=18)
res.marker[['DNMT1']] + res.marker[['ERVFRD-1']] + res.marker[['PSG8']] + res.marker[['CGA']] +
res.marker[['SH3TC2']] + res.marker[['LEP']] + res.marker[['PAPPA']] + res.marker[['FLT1']]  +
plot_layout(ncol=4,nrow=2)
#print(res.marker[['FLT1']], vp=viewport(angle=-185))
ggsave(filename = 'pdfs/marker_genes/marker.gene.umap.pdf',height=10.5,width=18,useDingbats=FALSE) #can save grid multiple plots

ggsave(filename = 'pdfs/marker_genes/marker.gene.umap.usemagic.pdf',height=10.5,width=18,useDingbats=FALSE)


##plot all grid
options(repr.plot.height=5.5,repr.plot.width=5.5)
for(i in names(res.marker)){
    print(res.marker[[i]])
    ggsave(filename = paste('pdfs/marker_genes/all/marker.gene.',i,'.umap.pdf',sep=''),height=5.5,width=5.5,useDingbats=FALSE)
}







##############do all DEG identification (test on data slot by default)#############

placenta.markers <- FindAllMarkers(placenta, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

placenta.markers.loose <- FindAllMarkers(placenta, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

placenta.markers.veryloose <- FindAllMarkers(placenta, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
##use this, then filter from result

#set min.pct to a small value to detect ERVFRD-1

#saveRDS(placenta.markers,'DEGs/placenta.markers.rds')
#saveRDS(placenta.markers.loose,'DEGs/placenta.markers.loose.rds')
#saveRDS(placenta.markers.veryloose,'DEGs/placenta.markers.veryloose.rds')

placenta.markers <- readRDS('DEGs/placenta.markers.rds') #strict
placenta.markers.loose <- readRDS('DEGs/placenta.markers.loose.rds')
placenta.markers.veryloose <- readRDS('DEGs/placenta.markers.veryloose.rds')


#start at 2023-02-01 5:04
tick <- paste0('Done at ',Sys.time())

'Done at 2023-02-01 06:14:16'
#'Done at 2023-01-18 04:38:19'



table(placenta.markers$cluster)
 1    2    3    4    5    6    7    8    9   10   11 
 117  184  208   66  186 1536 1974  135 1467  223  763

table(placenta.markers.loose$cluster)
  1    2    3    4    5    6    7    8    9   10   11 
 117  188  222   66  192 1683 2117  135 1554  242  814

table(placenta.markers.veryloose$cluster)
  1    2    3    4    5    6    7    8    9   10   11 
 467  596  733  422  665 4225 4785  701 4082  790 1550 



#de.markers.all <- placenta.markers #6859
#de.markers.all <- placenta.markers.loose #7330
de.markers.all <- placenta.markers.veryloose #19016


##save deg table
saveRDS(de.markers.all,'DEGs/de.markers.all.rds') #1,3,4,5


###select top/max n

##use this, sorted by pvalue then weighted(sorted) by avg_logFC
top25 <- de.markers.all %>% group_by(cluster) %>% dplyr::top_n(n = 25, wt = avg_logFC) #sort by pvalue, weight by logFC
top50 <- de.markers.all %>% group_by(cluster) %>% dplyr::top_n(n = 50, wt = avg_logFC)
top100 <- de.markers.all %>% group_by(cluster) %>% dplyr::top_n(n = 100, wt = avg_logFC)
top1000 <- de.markers.all %>% group_by(cluster) %>% dplyr::top_n(n = 1000, wt = avg_logFC)

#not use this, pvalue lost order
max25 <- de.markers.all %>% group_by(cluster) %>% dplyr::slice_max(n = 25, order_by = avg_logFC)
max50 <- de.markers.all %>% group_by(cluster) %>% dplyr::slice_max(n = 50, order_by = avg_logFC)
max100 <- de.markers.all %>% group_by(cluster) %>% dplyr::slice_max(n = 100, order_by = avg_logFC)
max1000 <- de.markers.all %>% group_by(cluster) %>% dplyr::slice_max(n = 1000, order_by = avg_logFC)

##true head n without any sorting
head25 <- de.markers.all %>% group_by(cluster) %>% dplyr::slice_head(n = 25) #only head, no sorting
head50 <- de.markers.all %>% group_by(cluster) %>% dplyr::slice_head(n = 50)
head100 <- de.markers.all %>% group_by(cluster) %>% dplyr::slice_head(n = 100)
head1000 <- de.markers.all %>% group_by(cluster) %>% dplyr::slice_head(n = 1000)



subset(de.markers.all,cluster == '3')[1:25,'gene',drop=TRUE]
'MYCNUT''FSTL3''PFKP''PLCL2''IL1RAP''MIR193BHG''FLT1''MYCN''SH3BP5''GPRC5C''MYCNOS''EGLN3''GDPD4''PCED1B''LIMD1''TGFB1''GALNT2''MTUS1''EZR''MTSS2''NEURL1''SASH1''IGF2BP2''MALAT1''POU2F3'


subset(head25,cluster == '3')[1:25,'gene',drop=TRUE] #use this
'MYCNUT''FSTL3''PFKP''PLCL2''IL1RAP''MIR193BHG''FLT1''MYCN''SH3BP5''GPRC5C''MYCNOS''EGLN3''GDPD4''PCED1B''LIMD1''TGFB1''GALNT2''MTUS1''EZR''MTSS2''NEURL1''SASH1''IGF2BP2''MALAT1''POU2F3'


subset(top25,cluster == '3')[1:25,'gene',drop=TRUE]
'MYCNUT''FSTL3''PFKP''PLCL2''IL1RAP''MIR193BHG''FLT1''MYCN''SH3BP5''EGLN3''PCED1B''LIMD1''TGFB1''GALNT2''SASH1''HDAC2''NEK11''COL27A1''TCL6''LVRN''AJ009632.2''GPC5''MME''LEP''XACT'

subset(max25,cluster == '3')[1:25,'gene',drop=TRUE]
'MYCNUT''COL27A1''EGLN3''FSTL3''PFKP''XACT''AJ009632.2''PLCL2''HDAC2''IL1RAP''TCL6''MIR193BHG''MME''LVRN''TGFB1''NEK11''FLT1''LEP''MYCN''SH3BP5''SASH1''GALNT2''GPC5''LIMD1''PCED1B'



table(top50$cluster)
1  2  3  4  5  6  7  8  9 10 11 
50 50 50 50 50 50 50 50 50 50 50


table(top100$cluster)
1   2   3   4   5   6   7   8   9  10  11 
100 100 100  66 100 100 100 100 100 100 100

#   1   2   3   4   5   7   8   9  10  11  12  13  14  16  17  18 
# 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100

#  0   1   2   3   4   5   6   7   8   9 
# 100 100 100 100 100 100 100 100 100 100


table(top1000$cluster)
  1    2    3    4    5    6    7    8    9   10   11 
 117  188  222   66  192 1000 1000  135 1000  242  814 

#  1    2    3    4    5    7    8    9   10   11   12   13   14   16   17   18 
# 134  213  240  103  214 1000 1000  160 1000  257  810 1000 1000 1000 1000 1000

table(max1000$cluster)

1    2    3    4    5    6    7    8    9   10   11 
 117  188  222   66  192 1000 1000  135 1000  242  814

#  1    2    3    4    5    7    8    9   10   11   12   13   14   16   17   18 
# 134  213  240  103  214 1000 1000  160 1000  257  810 1000 1000 1000 1000 1000 



####test deg for c10 and c3 pairwisely

placenta.markers.loose.c10c3 <- FindMarkers(placenta, ident.1 = '10', ident.2 = '3', only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)

#saveRDS(placenta.markers.loose.c10c3,'DEGs/placenta.markers.loose.c10c3.rds')
placenta.markers.loose.c10c3 <- readRDS('DEGs/placenta.markers.loose.c10c3.rds')


table(placenta.markers.loose.c10c3$avg_logFC > 0) #537
FALSE  TRUE 
  216   321


placenta.markers.loose.c10c3.sort <- placenta.markers.loose.c10c3[order(placenta.markers.loose.c10c3$avg_logFC,decreasing = TRUE),]

deg.c10 <- rownames(placenta.markers.loose.c10c3.sort[placenta.markers.loose.c10c3.sort$avg_logFC >0,] )
'SLC26A7''PAPPA''LINC01483''ADAMTSL1''ZNF117''SCIN''AC022872.1''LINC00578''ABCB1''AC022146.2''ADCY7''AUTS2''LAMA3''ISM1''FBN2''SLC27A6''EPHA1-AS1''KIF6''CATSPERB''ANK3''NAALADL2''SVEP1''AC093817.2''RALGAPA2''IQGAP2''STARD13''TMEM108''KMO''AC004704.1''RNF103-CHMP3''ATP8A2''AC009126.1''APBB2''CDYL2''SLC45A4''NRK''TRPV6''BCKDHB''PTPRG''NID1''NPAS2''RRAS2''CSHL1''ARHGAP17''CLRN1-AS1''ARHGAP26''CSH1''CCSER1''ADGRL3''THSD4''KATNBL1''AC069277.1''POSTN''LINC00882''ZC2HC1A''IL1RAPL2''AC016831.7''AC119674.1''PCDH11X''CERS4''LNX1''ST5''MYO9A''FER1L5''RAB17''GHR''LSS''MORC4''FBLN1''TGFBR3''GABRE''GLCCI1''AC079160.1''NMT2''AL162718.1''KIAA1211''LINC02860''XDH''CPS1''PTPRJ''ITPK1''CLIC5''ARHGAP42''NLRC5''GPC3''BCL2''CSH2''DTNB''DDX5''MIR100HG''KIAA0319''PAG1''PLXDC2''LRCH1''SLC15A1''TCHH''CBR3-AS1''OTOGL''PLCB1''ELMO1''LINC00474''CPXM2''GNG7''GNG12''NFU1''MICU3''PTPRM''STS''LINC02365''PHYHIPL''PLPP3''MKLN1''SGPP2''IMMP2L''AIG1''PCDH11Y''SCMH1''MAPK8''FNDC3A''LINC01091''TWSG1''ERV3-1''DOCK1''MPP5''ICA1''STAT4''MGAT5''SLCO2A1''MB21D2''STON2''TGM2''RIN2''LINC01807''ABCD3''NMNAT2''SLC25A35''PTPRQ''TXNRD1''TTBK2''KLRD1''TMEM164''KLHL5''SIPA1L1''FOXP1''WWC1''CLDN14''AHR''PDE4D''MTSS1''COMMD1''LUZP1''AP003181.1''AR''AC005532.1''SEMA6A''SPRY4-AS1''SUPT3H''PRKCE''ACSS1''RIMKLB''LINC02109''SLC4A4''AP000331.1''FRMD6-AS2''SGSM1''MTMR7''LINC01456''GRIP1''SLC23A2''LINC01949''FHIT''COBL''LINC01119''GH2''RAB3B''GNGT1''AL136962.1''ST3GAL3''PIP5K1B''SLC39A8''SLC7A2''AC073569.2''GREM2''ANO10''AC011287.1''LINC02267''HSD3B1''AC099520.1''WWOX''LRP8''ADAM10''XPO6''ANGPT2''IKZF2''TRIM5''AL138828.1''SSH2''ABCG2''SH3TC2-DT''DDX17''HIBADH''CCDC171''GRM7''GRAMD2B''GSTA4''TC2N''RNGTT''CARD18''AC013287.1''TTC28''PLEKHG1''GULP1''ANKDD1A''AL023574.1''COX10-AS1''COLEC12''PRKCA''CYTH3''THSD7A''PPP2R3A''LYST''CASTOR3''GSK3B''RANBP17''AGAP1''ZNF554''PTPN1''CHST11''CFAP299''AL365295.1''POU6F2''AL589787.1''NEK7''BCR''FAM234B''ZHX2''HSD11B2''GRK5''AC139718.1''MAP3K4''RNF217''TRAPPC9''NFIA''CAST''MAP3K5''COMT''PBX1''PLPP1''OGFRL1''TULP4''SPIDR''TBC1D5''TNS1''SYT12''EXT1''SLC19A3''KIAA1671''CHKA''MED12L''FAM126A''LITAF''AL589740.1''GLUL''ACACA''ARNTL2''GAB2''TGFBR2''MYO1B''CSGALNACT1''ARRB1''ATF7IP2''FAM172A''TWIST1''DLGAP4''PLEKHA7''AMPD3''BPGM''MIR4435-2HG''MAP4K3''NCF1''MAP4K3-DT''CTDSP1''PSG2''TIAM1''ADORA1''TBC1D30''ATXN7L1''HLCS''AC105411.1''TBC1D1''MGST1''NHSL1''LIMS1''LINC00278''LNPEP''CLMN''ARHGAP10''PTH2R''ZDHHC14''AC079760.1''TJP1''AC006378.1''SUCLG2''TMEM120B''GAB1''AC087854.1''ARFGEF3''PPP1R12B''DENND1A''USP46''SRSF6''MAGI1''TSGA10''PDLIM5''SELENOI''HMGB1''NCOA2''FAM171B''AMFR''RERG''AC009264.1'

deg.c3 <- rev(rownames(placenta.markers.loose.c10c3.sort[placenta.markers.loose.c10c3.sort$avg_logFC <0,] ))
'MYCNUT''EGLN3''COL27A1''FSTL3''MME''ARHGAP24''XACT''MIR193BHG''PFKP''BNIP3L''SASH1''FLT1''FAM83D''GPC5''GASK1A''SH3BP5''LIMD1''SLC2A1''SERPINE1''SLC6A8''HDAC2''IL1RAP''ARID3A''MTHFD1L''ANKRD37''PLCL2''ANXA1''MYCN''SIGLEC6''MTUS1''ATP10D''RASSF8''NPEPPS''NECTIN4''SLC2A3''AL513164.1''C4orf47''TBL1XR1''TCL6''PFKFB4''NABP1''MYOF''AC108690.1''ERRFI1''MTSS2''SYNJ2''AC087857.1''FRMD4A''TNFRSF1B''AC020916.1''ZNF331''GATA2''GMCL1''MFSD12''PLIN2''ITPRID2''LPCAT1''ATP2C2''FLNB''GPRC5C''FOCAD''LVRN''EZR''NDRG1''TDRP''NEURL1''MIR29B2CHG''PLEKHM2''MALAT1''GSE1''ADAMTS20''SEMA3B''CMIP''HILPDA''ZNF626''FAM153CP''MGAT3''CORO1C''MAPK8IP3''NEK11''BASP1''HPCAL1''JUP''ARHGAP45''EGLN1''MYCNOS''ARHGAP30''ERO1A''UBC''KRBOX1''INPP5B''FAM43A''AC026167.1''CSNK1E''MID1''DIAPH2''P4HA1''JPT2''SPATA13''TPP2''PNCK''FAM120A''PTPRR''GALNT2''HIST1H2AC''AHNAK''ESYT2''MLLT10''NEAT1''DHX35''ZFP36L1''RABL6''CCDC138''INHBA''KMT2E''CGB3''SHANK2''TGFB1''HIF1A-AS3''CHD2''CALM1''PPP1R12C''HIST1H4H''TTC3''RAB8B''GPR78''ACOT11''MFSD2A''GDPD4''ZMYND8''HMG20B''MXD1''CYP11A1''DDB1''NCOR2''EBP''SLC7A5''SP3''GTF2IRD1''SEMA7A''ITGA5''HTRA4''SH3PXD2A''UCA1''LHB''PCNP''BCL6''NRIP1''FOXJ3''ERGIC1''UBASH3B''KIAA1217''HK2''DLEU2''DGKI''PRDX6''CBLB''DGKZ''TREML2''TPD52''TRIM14''MXI1''AC005786.3''AL645568.1''FHL2''POF1B''PHIP''NGLY1''WDR60''TTLL5''CCNY''TMEM91''MIR503HG''RAB3GAP1''EGFR''DENND4A''FMN1''PHKA2''PLOD2''AC018754.1''PKM''GMDS''PHACTR2''PIM3''FZR1''EFHD1''REEP3''CAPN7''AFAP1''PRKD3''DNAJB6''SCARB1''POMP''SLC16A3''TPBG''AC011453.1''TET3''TNS3''LINC02832''RIPK2''ZNF175''PICALM''MTMR4''SERTAD2''PTDSS1''HES2''INHA''ZCCHC2''TET1''GATA3''PCOLCE2''PORCN''RNMT''DOCK5''SLC1A6''ANGPTL4'


###whether or not important gene in deg list ?


top50[grep('^DNMT1$',top50$gene),] #6 7  9

top50[grep('^ERV',top50$gene),] #7 11 
top50[grep('^OVOL',top50$gene),] #no

top50[grep('^SH3TC2$',top50$gene),] #8
top50[grep('^GCM',top50$gene),] #no

top50[grep('^PAPPA$',top50$gene),] #2, 10
top50[grep('^LAMA3$',top50$gene),] #5,10
top50[grep('^CSHL1',top50$gene),] #2
top50[grep('^CSH1',top50$gene),] #2
top50[grep('^ADAMTS6',top50$gene),] #2

top50[grep('^STAT4$',top50$gene),] #no
top50[grep('^STAT',top50$gene),] #no
top50[grep('^AR$',top50$gene),] #no

top50[grep('^FLT1$',top50$gene),] #1,3,4,5
top50[grep('^ENG$',top50$gene),] #1
top50[grep('^INHBA$',top50$gene),] #3,4,5
top50[grep('^INHA$',top50$gene),] #no

top50[grep('^MYCN$',top50$gene),] #3
top50[grep('^LVRN$',top50$gene),] #3



top100[grep('^DNMT1$',top100$gene),] #6 7  9

top100[grep('^ERVFRD-1',top100$gene),] #11
top100[grep('^OVO',top100$gene),] 

top100[grep('^SH3TC2$',top100$gene),] #8,2
top100[grep('^GCM',top100$gene),] #2

top100[grep('^PAPPA$',top100$gene),] #2,5, 10
top100[grep('^LAMA3$',top100$gene),] #5,11
top100[grep('^CSHL1',top100$gene),] #2
top100[grep('^CSH1',top100$gene),] #2
top100[grep('^ADAMTS6',top100$gene),]  #2

top100[grep('^STAT4$',top100$gene),] #no
top100[grep('^STAT',top100$gene),] #no
top100[grep('^AR$',top100$gene),] #no

top100[grep('^FLT1$',top100$gene),] #1,3,4,5
top100[grep('^ENG$',top100$gene),] #1
top100[grep('^INHBA$',top100$gene),]#1,3,4,5

top100[grep('^MYCN$',top100$gene),] #3
top100[grep('^LVRN$',top100$gene),] #3,5


top1000[grep('^DNMT1$',top1000$gene),] #6 7  9

top1000[grep('^ERVFR',top1000$gene),] #11, many if ERV! 
top1000[grep('^OVOL',top1000$gene),] #no

top1000[grep('^SH3TC2$',top1000$gene),] #8,2
top1000[grep('^GCM',top1000$gene),] #n1,2,8

top1000[grep('^PAPPA$',top1000$gene),] #2,5, 10
top1000[grep('^LAMA3$',top1000$gene),] #5,10
top1000[grep('^CSHL1',top1000$gene),] #2, 10
top1000[grep('^CSH1',top1000$gene),] #2,10
top1000[grep('^ADAMTS6',top1000$gene),] #2

top1000[grep('^STAT4$',top1000$gene),] #10
top1000[grep('^STAT5A',top1000$gene),] #no (need very loose!,0.178 avg_logFC)
top1000[grep('^AR$',top1000$gene),] #10

top1000[grep('^FLT1$',top1000$gene),] #1,3,4,5
top1000[grep('^ENG$',top1000$gene),] #1,3
top1000[grep('^INHBA$',top1000$gene),]#1,3,4,5

top1000[grep('^MYCN$',top1000$gene),] #3 #14 (EVT)
top1000[grep('^LVRN$',top1000$gene),] #3,5,10



##############compare to previous deg table , and compare deg within clusters###########


# shared_gene_pappa <- intersect(subset(top100, cluster == '10')[1:50,'gene',drop=TRUE], 
#                          subset(top50, cluster == '10')[,'gene',drop=TRUE]
#                         ) #34

# length(shared_gene) #33 of 50



# shared_gene_flt1 <- intersect(subset(top100, cluster == '10')[1:50,'gene',drop=TRUE], 
#                          subset(top50, cluster == '10')[,'gene',drop=TRUE]
#                         ) #34

# length(shared_gene) #33 of 50


#read in previous result
top50_previous <- read.table('/home/mjwang/pwdex/placenta_10X_combine/02.seurat_harmony/DEGs/top50.de.gene.txt',header = TRUE,sep = '\t')


#length(    intersect(subset(top50_previous, cluster == '3')[,'gene',drop=TRUE], subset(top50, cluster == '10')[,'gene',drop=TRUE]) ) #38 of 50


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
    intersect(datalist[[1]],datalist[[2]])
    #intersect(datalist[[2]],datalist[[1]])
    
    rowid <- vector()
    ##get combined id (order by ref if set, if not set, use the sortedd allid.combined as row id)
    if( is.null(refid) ){
#       ref <- datalist[[1]]
        rowid <- sort(allid.combined)
    }else{ 
           ref = datalist[[refid]]
           notref <- allid.combined[!allid.combined %in% ref]
           rowid <- c(ref,notref)
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
                               na_col = 'white', 
                               color = color_use,
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


###compare deg with previous deg result
de.list <- list(
    #'de.PAPPA'= subset(top50, cluster == '10')[,'gene',drop=TRUE], 
    #'de.PAPPA.pre' = as.character(subset(top50_previous, cluster == '3')[,'gene',drop=TRUE]) 
    
    #'de.FLT1'= subset(top50, cluster == '3')[,'gene',drop=TRUE], 
    #'de.FLT1.pre' = as.character(subset(top50_previous, cluster == '7')[,'gene',drop=TRUE]) 

    #'de.FLT1'= subset(top100, cluster == '3')[,'gene',drop=TRUE], 
    #'de.FLT1.pre' = as.character(subset(top50_previous, cluster == '7')[,'gene',drop=TRUE]) 
    
    #'de.FLT1'= unique(subset(top100, cluster %in% c('3','1') )[,'gene',drop=TRUE]), 
    #'de.FLT1.pre' = unique(as.character(subset(top50_previous, cluster %in% c('7','4') )[,'gene',drop=TRUE]))
    
    'de.nascent'= subset(top50, cluster == '8')[,'gene',drop=TRUE], 
    'de.nascent.pre' = as.character(subset(top50_previous, cluster == '6')[,'gene',drop=TRUE]) 
    
    
    
)

df.compare.pappa <- list2DF_keepref(datalist = de.list, ref = 'de.PAPPA',plot =TRUE , title = 'PAPPA with \nprevious', height = 10.5, width = 3.5 )
# shared ids length is: 38
# shared ids are: SLC26A7,PAPPA,SLCO2A1,ZNF117,LAMA3,EPHA1-AS1,ISM1,AC022146.2,GPC3,NPAS2,SLC4A4,POSTN,PTPRJ,CATSPERB,KIF6,ADCY7,FER1L5,CLIC5,ARHGAP26,KIAA1211,CPS1,KATNBL1,ADAMTSL1,THSD4,AC009126.1,ST5,IQGAP2,ABCB1,SCIN,HOPX,AC022872.1,LINC01483,AJ009632.2,AC016831.7,AC006378.1,ARHGEF28,AC004704.1,RNF103-CHMP3

# combined id length is 62

df.compare.flt1 <- list2DF_keepref(datalist = de.list, ref = 'de.FLT1',plot =TRUE , title = 'FLT1 with \nprevious', height = 25.5, width = 3.5 )

# shared ids length is: 28
# shared ids are: MYCNUT,FSTL3,PFKP,PLCL2,IL1RAP,MIR193BHG,FLT1,MYCN,SH3BP5,EGLN3,PCED1B,LIMD1,TGFB1,GALNT2,SASH1,IGF2BP2,INHBA,HDAC2,NEK11,COL27A1,PHACTR2,FAM83D,AJ009632.2,ANKRD37,MME,TNFRSF1B,LEP,SERPINE1

# combined id length is 72

df.compare.nascent <- list2DF_keepref(datalist = de.list, ref = 'de.nascent',plot =TRUE , title = 'Nascent with \nprevious', height = 10.5, width = 3.5 )



##compare with in cluster de top50 gene
de.list <- lapply(split(top50, f = top50$cluster), function(x){ x$gene } )

df.compare.c10vsc3 <- list2DF_keepref(datalist = de.list[c('10','3')], ref = '10',plot =TRUE , title = 'PAPPA(c10) with \nFLT1(c3) deg', height = 12.5, width = 3.5 )


df.compare.c10vsc3 <- list2DF_keepref(datalist = de.list[c('10','2','3','1')], ref = '10',plot =TRUE , title = 'PAPPA(c10,c2) with \nFLT1(c3,c1) deg', height = 20.5, width = 3.5 )




##do unexpected genes (sex-related, MT, Ribosomal,cell-cycling  gene, ) in deg?#######


#read chrY genes and chrX genes
#get sex information , read in cellranger 7 gene.gtf chrY genes
genes_chrXY.df <- read.table("data/genes.chrXY.bed") #
genes_chrXY.split.df <- genes_chrXY.df %>% separate(col = V4, sep = "\\|", into = c('geneid','gene_type') )
colnames(genes_chrXY.split.df) <- c('chr','start','end','geneid','gene_type','strand')

sex_genes <- genes_chrXY.split.df$geneid #all chrX chrY genes include protein_coding
sex_genes.X <- subset(genes_chrXY.split.df, chr == 'chrX')$geneid
sex_genes.Y <- subset(genes_chrXY.split.df, chr == 'chrY')$geneid

flag_sex_X <- which(deg.c3 %in% sex_genes.X)
deg.c3[flag_sex_X]
'XACT''SLC6A8''MID1''DIAPH2''PNCK''EBP''POF1B''MIR503HG''PHKA2''PORCN'

flag_sex_Y <- which(deg.c3 %in% sex_genes.Y)
deg.c3[flag_sex_Y]
#no

flag_sex_X <- which(deg.c10 %in% sex_genes.X)
deg.c10[flag_sex_X]
'NRK''IL1RAPL2''PCDH11X''MORC4''GABRE''GPC3''STS''TMEM164''AR''LINC01456''AL023574.1'

flag_sex_Y <- which(deg.c10 %in% sex_genes.Y)
deg.c10[flag_sex_Y]
'PCDH11Y''LINC00278'


flag_sex_X <- which(top1000$gene %in% sex_genes.X)
top1000$gene[flag_sex_X]
#274

##do these X deg gene expressed both in female and male samples?

options(repr.plot.height = 7.5, repr.plot.width = 15)
FeaturePlot(placenta, features = 'PAPPA', reduction = "umap_rotate",slot = 'data',cols = c("lightgrey", "darkred"), split.by = 'sex' )
FeaturePlot(placenta, features = 'FLT1', reduction = "umap_rotate",slot = 'data',cols = c("lightgrey", "darkred"), split.by = 'sex' )

FeaturePlot(placenta, features = 'XACT', reduction = "umap_rotate",slot = 'data',cols = c("lightgrey", "darkred"), split.by = 'sex' ) #both expressed in STB Mature 1 (c3)!

FeaturePlot(placenta, features = 'XIST', reduction = "umap_rotate",slot = 'data',cols = c("lightgrey", "darkred"), split.by = 'sex' )

flag_sex_Y <- which(top1000$gene %in% sex_genes.Y)
unique(top1000$gene[flag_sex_Y])
#'PCDH11Y' 'TTTY14' 'UTY' 'LINC00278' 4 of 111 only ! generally few sex-bias except three chrY protein coding gene.


FeaturePlot(placenta, features = 'PCDH11Y', reduction = "umap_rotate",slot = 'data',cols = c("lightgrey", "darkred"), split.by = 'sex' ) #both expressed in STB Mature 1 (c3)!


##cell cycling gene 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

flag_cellcycle <- which(top1000$gene %in% c(s.genes,g2m.genes))
top1000$gene[flag_cellcycle]
'RANGAP1''CDCA7''POLA1''TOP2A''ATAD2''SMC4''BRIP1''HELLS''CENPE''DTL''POLA1''CENPF''RRM2''NUSAP1''KIF23''GTSE1''ANLN''ECT2''KIF11''TACC3''KIF20B''MKI67''BLM''TMPO''TPX2''CDCA2''NASP''WDR76''KIF2C''NUF2''CLSPN''NCAPD2''CDK1''EXO1''GAS2L3''HMMR''HMGB2''CDC25C''MSH2''NDC80''CKAP2''GMNN''TYMS''RRM1''TTK''BUB1''CDC45''AURKB''CDCA8''CDCA7''UHRF1''CKAP2L''SLBP''DLGAP5''AURKA''POLD3''G2E3''CKAP5''CASP8AP2''CTCF''POLA1''RANGAP1''GAS2L3'

table(top1000[flag_cellcycle,'cluster']) #main in c7 CTB proliferation
 1  2  3  4  5  6  7  8  9 10 11 
 0  0  1  0  0  2 57  0  1  1  1 

##MT and ribosomal



##sort df by cluster of level##


#topx <- max50 #do not use this
topx <- top50

topx$cluster <- factor(topx$cluster, levels = c('7','9','6','11','8','4','2','10','5','1','3' ) )

marker.genes.de <- data.frame()
n = 25
#for(i in levels(topx$cluster) ){
for(i in c('7','9','6','11','8','2','10','5','1','3' ) ){ #remove deg of c4
    topx.sel <- subset(topx,cluster == i)
    marker.genes.de <- rbind.data.frame(marker.genes.de,topx.sel[1:n,] )

}


##save to rds, txt and xlsx

saveRDS(marker.genes.de,'DEGs/marker.genes.de.top25.sorted.rds')

write.table(file='DEGs/marker.genes.de.top25.sorted.txt',x=marker.genes.de,quote = FALSE,sep = '\t',row.names = TRUE,col.names = TRUE)

write.xlsx( marker.genes.de, 
            file = paste('DEGs/marker.genes.de.top25.sorted.xlsx',sep=''),
            rowNames = FALSE,
            overwrite = TRUE
          )



###quick look for c3 and c10
for (i in subset(marker.genes.de,cluster == '10')$gene){
#for (i in subset(marker.genes.de,cluster == '3')$gene){
    #cat(i,'\n',sep='')
    cat("\'",i,"\'",',',sep='')
}
#c10
'SLC26A7','PAPPA','SLCO2A1','ZNF117','LAMA3','EPHA1-AS1','ISM1','AC022146.2','GPC3','NPAS2','SLC4A4','POSTN','ATP8A2','PTPRJ','CATSPERB','KIF6','ADCY7','FER1L5','CLIC5','PAG1','CLRN1-AS1','ARHGAP26','KIAA1211','IGF2BP2','CPS1',

#c3
'MYCNUT','FSTL3','PFKP','PLCL2','IL1RAP','MIR193BHG','FLT1','MYCN','SH3BP5','EGLN3','PCED1B','LIMD1','TGFB1','GALNT2','MTUS1','MTSS2','NEURL1','SASH1','IGF2BP2','INHBA','INPP5B','ARHGAP24','TBL1XR1','HDAC2','HPCAL1'


########plot DAG heatmap (DoHeatmap, use scale.data by default) in a long heatmap####

##set column order

placenta@meta.data$cluster <- factor(placenta@meta.data$cluster, levels = c('7','9','6','11','8','4','2','10','5','1','3' ) )


##modify assay data by overwrite
placenta.smooth <- placenta

all.equal( dimnames(GetAssayData(placenta.smooth,slot = 'scale.data')), dimnames(exprMat.z.s) ) #TRUE

placenta.smooth@assays$RNA@data <- exprMat.z.s #dcgmatix can not assign scale.data (a matrix)


##cellbarcode random ?

cellbarcode.random <- vector()
cluster <- Idents(placenta.smooth)

#for(i in levels(topx$cluster) ){
for(i in c('7','9','6','11','8','2','10','5','1','3' ) ){ #remove deg of c4
    
    cluster.sel <- cluster[cluster == i]
    cluster.sel.id <- names(cluster.sel)
    cluster.sel.id.random <- sample(cluster.sel.id,size = length(cluster.sel.id),replace=FALSE)
    #all.equal(sort(cluster.sel.id),sort(cluster.sel.id.random)) #TRUE
    #all.equal(cluster.sel.id,cluster.sel.id.random)#no
    cellbarcode.random <- c(cellbarcode.random,cluster.sel.id.random)
}

saveRDS(cellbarcode.random,'cellbarcode.random.rds') #no c4 cellbarcode
cellbarcode.random <- readRDS('cellbarcode.random.rds')


##get random obj
placenta.smooth.random <- placenta.smooth[,cellbarcode.random]
placenta.smooth.random[['cluster']]  <- droplevels(placenta.smooth.random[['cluster']])

# ##choose genes to show, and set them with name, others be empty
# #https://github.com/satijalab/seurat/issues/1395
# #not ok, use complexHeatmap instead

# all.genes <- marker.genes.de$gene#[1:10]

# #genes.to.label <- sample(x = all.genes, size = 5) 
# genes.to.label <- c('CENPP','DNMT1','ERVFRD-1','SH3TC2','LAMA3','PAPPA','FLT1','CSHL1')

# #genes.to.label <- c('SCLT1','CENPP','SDK1','NSD2')

# # labels <- rep(x = "transparent", times = length(x = all.genes))
# # labels[match(x = genes.to.label, table = all.genes)] <- "black"

# labels <- rep(x = "", times = length(x = all.genes))
# #idx <- match(x = genes.to.label, table = all.genes) #only return the first hit
# idx <- vector()
# for(i in genes.to.label){
#     idx <- c(idx,grep(pattern = paste0('^',i,'$') , x = all.genes) )
    
# }


# labels[idx] <- all.genes[idx]


pdf(file = 'DEGs/de.marker.top50.heatmap.withlegend.pdf',width = 5.5,height = 7.5,useDingbats = FALSE)
pdf(file = 'DEGs/de.marker.top50.heatmap.nolegend.pdf',width = 5.5,height = 7.5,useDingbats = FALSE)

###plotting deg heatmap at single cell level, quick and ok but use complexheatmap with label annotation will be better
options(repr.plot.height = 7.5, repr.plot.width = 5.5)
DoHeatmap(placenta.smooth.random, 
          #features = top10$gene,
          #features = top25$gene,
          ##features = top25.dedup$gene,
          #features = top50$gene,
          features = marker.genes.de$gene,
          #features = c(deg.c10[1:50],deg.c3[1:50]),
          slot='data',#'scale.data' in fact,
          group.by = 'cluster',
          group.colors = unlist(rna_map_cellcolor),#color_good,#unlist(map_cellcolor)#,
          disp.min = -2,
          disp.max=2,
          group.bar = TRUE,
          label = TRUE
          
         ) +
  scale_fill_gradientn(colors =  color_gradient_my)+
  #scale_fill_gradientn(colors =  c('blue','white','red'))#+
  #scale_fill_gradientn(colors =  viridis(256, option = "D") )#+
  #scale_fill_gradientn(colors =  mapal)#+
  #scale_fill_gradientn(colors =  color_pyreds)#+ #not this
  ##scale_y_discrete(labels = rev(labels) ) + #use this! DoHeatmap use geom_raster/geom_tile, but shift too much,#not ok, use complexHeatmap instead
  #coord_flip() + 
  NoLegend()  
  #theme(axis.text.y = element_text(size = 0)) #remove genename if size = 0
  ##theme(axis.text.y = element_text(size = 6,hjust=1,vjust=1.0) ) 
  ##theme(axis.text.y = element_text(size = 6, color = rev(x = labels) )) #unexpected result, as warned by message

dev.off()





###########use ComplexHeatmap with row annotation at single cell level (similar to DoHeatmap but with row label)##########


#create mat and row labels
all.genes <- marker.genes.de$gene#[1:100]

all.genes <- all.genes[all.genes %in% rownames(exprMat.z.s) ]

mat <- as.matrix(exprMat.z.s[all.genes,cellbarcode.random])

#all.equal(placenta.smooth[['cluster']],placenta[['cluster']])
#all.equal(placenta.smooth[['cluster']],cluster.df.add['cluster'],check.attributes = FALSE) #TRUE
#all.equal(placenta.smooth.random[['cluster']],cluster.df.add['cluster'],check.attributes = FALSE) #no

cluster <- placenta.smooth.random[['cluster']] #a dataframe

all.equal(rownames(cluster),colnames(mat)) #TRUE


##column cluster annotation

colnames(cluster) <- 'Cell type'
ha_top <- columnAnnotation(
     df = cluster, #the dataframe, each column will annotated
     col = list('Cell type'=c(unlist(rna_map_cellcolor),'na'='white') ), #a color list with name
     #simple_anno_size = unit(0.1, "mm")
     show_annotation_name = FALSE
)


###left deg group belonging annotation

DEG_type <- as.data.frame(marker.genes.de[,c('cluster','gene')])
colnames(DEG_type) <- c('DEGs','gene')

##for ID_map add dar type
#DAR_type.bk <- DAR_type
#sum(duplicated(DAR_type.bk$peak )) #0
#rownames(DAR_type.bk) <- DAR_type.bk$peak

##DAR_type <- DAR_type[,2,drop=FALSE]

ha_left <- rowAnnotation(
     df = DEG_type[,1,drop=FALSE],
     col = list('DEGs'=c(unlist(rna_map_cellcolor),'na'='white') )
)





###right gene label annotation
genes.to.label <- c('MKI67','TEAD4','DNMT1','ERVFRD-1','SH3TC2','LAMA3','PAPPA','CSHL1','FLT1','ENG','MYCN')

idx <- vector()
for(i in genes.to.label){
    idx <- c(idx,grep(pattern = paste0('^',i,'$') , x = all.genes) )
    
}

idx.id <- all.genes[idx]

ha_right <- rowAnnotation(foo = anno_mark(at = idx,  #HeatmapAnnotation(... which = 'row')
                                    labels = idx.id,
                                    #labels = as.character(ID_map.use$Gene),
                                    labels_gp = gpar(fontsize = 12),
                                    link_width = unit(2, "mm"),
                                    extend = unit(5, "mm"),
                                    padding = unit(1, "mm") #important
                                    #annotation_width= unit(4, "cm")
                                   )

) 



res.hp = Heatmap(mat, name = "expr", 
         cluster_rows = FALSE, 
         cluster_columns = FALSE, 
         show_row_names = FALSE,
         show_column_names = FALSE,
         use_raster = FALSE,#will use raster if >2000 row or cols, however rstudio do not support raster
         ##col = circlize::colorRamp2(seq(-1.5,1.5,by=3/10), viridis(n = 11,option = "C")),
         #col = circlize::colorRamp2(seq(-1.5,1.5,by=3/(length(color_peak)-1)), color_peak),
         col = circlize::colorRamp2(seq(-2,2,by=4/(length(color_gradient_my)-1)), color_gradient_my),
         #col = circlize::colorRamp2(seq(-1.5,1.5,by=3/(length(color_tfdev1)-1)), color_tfdev1),
         #col = circlize::colorRamp2(seq(-1.5,1.5,by=3/255), color_peak),
         na_col = 'white',
         #heatmap_legend_param = list(color_bar = "continuous"),
         #right_annotation = ha#,heatmap_width=unit(8, "cm"),
         ##clustering_distance_rows  = 'pearson', 
         ##clustering_distance_columns  = 'pearson',
         column_split = cluster[,1],
         column_gap = unit(0.05,'cm'),
         column_labels = levels(cluster[,1]),
         column_names_side = 'top',
         top_annotation = ha_top,#trajectory bin cluster belonging
         #bottom_annotation = , #trajectory arrow and text, no use decorate_heatmap_body
         left_annotation = ha_left, #peak dar cluster belonging
         right_annotation = ha_right #peak annotation gene label
)


pdf(file = 'DEGs/pdfs/de.marker.top50.heatmap.withlegend.withArialMT.pdf',width = 7,height = 7.5,useDingbats = FALSE,fonts = NULL) #family to get linux fonts, fonts to get extral fonts file 
#pdf(file = 'DEGs/de.marker.top50.heatmap.nolegend.pdf',width = 5.5,height = 7.5,useDingbats = FALSE)

#"AvantGarde", "Bookman", "Courier", "Helvetica","Helvetica-Narrow", "NewCenturySchoolbook", "Palatino", "Times", Arial is a microsoft font, not availabel in linux


options(repr.plot.height = 7.5, repr.plot.width = 7)
draw(res.hp, heatmap_legend_side = "right",gap = unit(0.1, "cm")) 

dev.off()



pdf(file = 'DEGs/pdfs/test.withArialMT.pdf',width = 7,height = 7.5,useDingbats = FALSE,family = "ArialMT",fonts = NULL) 

plot(1:100)

dev.off()




##############plot aggre heatmap for deg #####


##get aggregated matrix##
#exprMat.ave <- AverageExpression(placenta_filter, slot = 'data')[[1]]
exprMat.ave <- AverageExpression(placenta, slot = 'data')[[1]]



##zscore normalize
#exprMat.ave.z <-  AverageExpression(placenta_filter,slot = "scale.data")[[1]]#,return.seurat = TRUE)
exprMat.ave.z <-  AverageExpression(placenta,slot = "scale.data")[[1]] #use this

#exprMat.ave.z <- AverageExpression(placenta.smooth, slot = 'data')[[1]] #add from scale.data in fact

# saveRDS(object = exprMat.ave,'exprMat.ave.rds')
# saveRDS(object = exprMat.ave.z,'exprMat.ave.z.rds')
# saveRDS(object = exprMat.ave.z,'exprMat.ave.z.s.rds')

exprMat.ave <- readRDS('exprMat.ave.rds')
exprMat.ave.z <- readRDS('exprMat.ave.z.rds')


#######################get gene pct.exp with Seurat::DotPlot for all rowid of exprMat.ave.z
exprMat.perc <- DotPlot(placenta, features = rownames(exprMat.ave.z) )
exprMat.perc <- exprMat.perc$data


rownames(exprMat.perc) <- NULL

table(exprMat.perc$id)
1     2     3     4     5     6     7     8     9    10    11 
29132 29132 29132 29132 29132 29132 29132 29132 29132 29132 29132 

#  1     2     3     4     5     6     7     8     9    10 
# 24307 24307 24307 24307 24307 24307 24307 24307 24307 24307 


all.equal ( as.character(unique( exprMat.perc$features.plot )), rownames(exprMat.ave.z) ,check.attributes = FALSE) #TRUE


subset(exprMat.perc, features.plot == "AL669831.5")
subset(exprMat.perc, features.plot == "PAPPA")
subset(exprMat.perc, features.plot == "FLT1")
subset(exprMat.perc, features.plot == "PSG8")


saveRDS(exprMat.perc,'exprMat.perc.rds')

#https://bioinformatics.stackexchange.com/questions/7026/how-can-i-obtain-the-percentage-gene-expression-per-identity-class-in-seurat-as







#use seurat DoHeatmap
##DoHeatmap(exprMat.ave.z, features = unlist(TopFeatures(placenta_filter[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE)

# marker.genes = c(
#     "DNMT1", "CDH1", "PAGE4",'ITGA6','MKI67','PCNA',
#     "FLT1", "CSHL1", "PSG8", "PAPPA",'GCM1','CGA','CGB8',
#     "ERVFRD-1", "LAIR2", "PLAC8",'HLA-G','MMP2',
#     'VIM','PECAM1','THY1','DLK1','HIVEP3','CD68','CD14','HLA-A','HLA-DPB1'
#   )


# ##STR markers
# marker.genes = c(
#     "VIM", "PECAM1", "THY1",
#     "NR2F1", "DLK1", "HIVEP3", 
#     "CD68", "CD14", "HLA-A",
#     'HLA-DPA1','HLA-DPB1','MKI67'
#   );


# marker.genes = c("LINC01605", #naive STB denovo spearman correlation detected gene signature
# "KLRD1",
# "BCL6",
# "NABP1",
# "AC108062.1",
# "HLF",
# "TENT5A",
# "CBLB"
# )

# marker.genes = c(#"TNFRSF1", #naive STB gene by liger RNA D1 D2 ATAC D1 D2
# "SH3TC2",
# "EGFR",
# "STS",
# "ELMO1",
# "TNFAIP8",
# #"NAALADL",
# "PREP",
# "DYSF",
# "SMS",
# "DCP2",
# "RYBP")


# marker.genes = c("AJ009632.2",  #LEP, FLT1 stb
# "ARHGEF28",
# "NEK11",
# "PAEP",
# "FAM153C",
# "LEP",
# "TGFB1",
# "PCED1B",
# "SLC7A1",
# "MOCOS",
# "GALNT2",
# "AC100802.1",
# "GRK3",
# "ALOX5",
# "FGD4",
# "PTPRM",
# "OAF",
# "IGF2BP2",
# "FAM167A",
# "FCGR3B",
# "INHBA",
# "SH3PXD2A"
# )

# marker.genes = c( #PAPPA stb
# "LINC01483",
# "AC006378.1",
# "SLC26A7",
# "SLCO2A1",
# "CPS1",
# "KIAA1671",
# "PAPPA",
# "THSD4",
# "KIF6",
# "AC119674.1",
# "ISM1",
# "HOPX",
# "ANGPT2",
# "SLC4A4",
# "LAMA3",
# "POSTN",
# "GNG7",
# "TRIM40",
# "SGSM1",
# "PTPRQ",
# "RNF103-CHMP3",
# "CLIC5"

# )

marker.genes <- marker.genes.de$gene #the same with Doheatmap

marker.genes <- subset(marker.genes.de,cluster %in% c('7','9','6','11','8','2','10','5','1','3' ) )$gene #remove deg of c4

marker.genes = c('FLT1','PAPPA','CGA','PSG1','ATF3','ID3','TCF4','INHBA','ANGPT2')

marker.genes = c('SUMO2','PIAS1','PIAS2','PIAS4','SENP5','SENP6','SENP2') #xiaozy

marker.genes = c("CSH2",'PSG1','PSG2','PSG3','PSG4','PSG5','PSG6','PSG7','PSG8','PSG9')



table(marker.genes %in% rownames(exprMat.ave.z))
TRUE 
 250

TRUE 
 110

#marker.mat <- t(exprMat.ave[marker.genes,])
marker.mat.z <- t(exprMat.ave.z[marker.genes,])


#marker.mat <- exprMat.ave[marker.genes,]
marker.mat.z <- exprMat.ave.z[marker.genes,c('7','9','6','11','8','2','10','5','1','3' )]


round(quantile(as.matrix(marker.mat.z)),digits = 2)
0% -0.79 25% -0.32 50% -0.12 75% 0.46 100% 255.

0% -1.63 25% -0.36 50% -0.17 75% 0.29 100% 3.46

#0% -1.12 25% -0.13 50% -0.05 75% 0.03 100% 7.14

##use complexHeatmap to visualize 
##plot the marker gene matrix heatmap


# bks <- seq(0, 100, by = 100/255) #breaks
options(repr.plot.height=7.5,repr.plot.width=5)
hp = Heatmap(marker.mat.z, name = "marker.mat.z", 
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             show_row_names = FALSE,
             use_raster = FALSE,
             ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             column_title = "DEG of early pregnancy",
             clustering_method_columns = "ward.D", ##ward.D,complete
             col = circlize::colorRamp2(c(-2, 0, 2), c("grey", "white", "red")) #set color and value cutoff at the same time
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
#width = max(grobWidth(textGrob(labels))))
draw(hp, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750



saveRDS(hp,file="hp.Rds",compress = TRUE)
##add a outline box frame
#decorate_heatmap_body("mostVar", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})
##get row order to avoid repeatly hclust
hp.roworder <- row_order(hp)
hp.colorder <- column_order(hp)

##draw hclust for columns directly
cluster.dist <- dist(t(top50.data),method = "euclidean")
plot(as.dendrogram(hclust(cluster.dist,method = "complete") ),main="hclust: complete euclidean" )
####






###save/load snapshot####

#saveRDS(placenta.cstb,'snapshot_rds/placenta.cstb.rds')
#saveRDS(cluster.df.add.cstb,'snapshot_rds/cluster.df.add.cstb.rds')




#saveRDS(placenta,'snapshot_rds/placenta.rds') #placenta.cstb.rds in fact
#saveRDS(cluster.df.add,'snapshot_rds/cluster.df.add.rds') #cluster.df.add in fact

#write.table(file='snapshot_rds/cluster.df.add.txt',x=cluster.df.add,quote = FALSE,sep = '\t',row.names = TRUE,col.names = TRUE)


placenta <- readRDS('snapshot_rds/placenta.rds') #placenta.cstb.rds in fact
cluster.df.add <- readRDS('snapshot_rds/cluster.df.add.rds') #cluster.df.add in fact


#saveRDS(placenta.smooth,'snapshot_rds/placenta.smooth.rds')


##replot
options(repr.plot.width = 7.5, repr.plot.height=7.5)
quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', title= 'early combined', shrink.x = 2, shrink.y = 1.05)


DimPlot(object = placenta, label = TRUE,cols=c(color_good),pt.size = 0.5,label.size = 10,reduction = "umap_rotate") + NoLegend() 




###############scripts end###############################

























# ####plot dot-heatmap for sequence depth
# options(repr.plot.height=7.5,repr.plot.width=7.5)
# ggplot(cluster.df.add ,aes(x=UMAP_1,y=UMAP_2,col=log10(nCount_RNA))  ) +
# #ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=TSS.enrichment)  ) +
#   geom_point(size = 1,show.legend = TRUE,alpha= 1 ) +
#   scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
#   theme_classic() +
#   theme(legend.position = 'bottom',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
#   ggtitle(paste(sample, "nCount_RNA",  sep=" ") ) +
#   #guides(shape = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
#   #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
#   #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
#   labs(x = "UMAP_1", y = "UMAP_2")
#   #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")


# #ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=cluster)  ) +
# ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=sample)  ) +
#   geom_point(size = 0.3,show.legend = TRUE,alpha= 0.3 ) +
#   #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
#   scale_color_manual(values = c('red','blue')) +
#   theme_classic() +
#   theme(legend.position = 'bottom',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
#   ggtitle(paste(sample, "Samples",  sep=" ") ) +
#   #guides(shape = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
#   #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
#   #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
#   labs(x = "UMAP_1", y = "UMAP_2") +
#   guides(col = guide_legend(override.aes = list(size = 8)))
#   #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")










# ####





# res.tsne <- placenta@reductions$tsne@cell.embeddings[,c(1,2)]
# res.umap <- placenta@reductions$umap@cell.embeddings[,c(1,2)]

# cl <- Idents(placenta)

# all.equal( names(cl)  , row.names(res.tsne) )  #true
# res.tsne.cl <- cbind(cluster=as.numeric(as.character(cl)),res.tsne)
# res.tsne.cl <- as.data.frame(res.tsne.cl)
# res.tsne.cl$cluster <- as.factor(res.tsne.cl$cluster)
# #all.equal(as.character(cl), as.character(res.tsne.cl$cluster) ) ##TRUE
# #write.csv(res.tsne,quote = FALSE, file="placenta.FCA7196220.tsne.csv", row.names = TRUE)
# write.table(res.tsne.cl,quote = FALSE, file="placenta.PLA-8w-RNA-new.tsne.cl.txt", row.names = TRUE,col.names = NA, sep='\t')


# res.umap.cl <- cbind(cluster=cl,res.umap)
# res.umap.cl = as.data.frame(res.umap.cl)
# write.table(res.umap.cl,quote = FALSE, file="placenta.PLA-8w-RNA-new.umap.cl.txt", row.names = TRUE,col.names = NA, sep='\t')



# #saveRDS(placenta,'placenta.PLA-8w-RNA-new.rds')
# placenta <- readRDS('placenta.PLA-8w-RNA-new.rds')

