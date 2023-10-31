

##integration and projection with seurat 4 dev
#https://satijalab.org/seurat/articles/integration_mapping.html


##1 BL and CT30 integration
##2 placenta villus CT30, BL integration 
##3 TF gene expression similar in vivo and in vitro (?)
##4 projection

#remotes::install_github(repo = 'satijalab/seurat', ref = 'develop') #to skip the spatstat problem

library('Seurat')
library('SummarizedExperiment')

#remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')

library('ggplot2')
library('pheatmap')
library('ComplexHeatmap')


#install metap for FindConservedMarkers
#install.packages('qqconf') #need R >= 4.0.0
#install.packages('metap')


##CCA can carry on integration and transfer cluster ids (without umap)

# devtools::install_github('satijalab/seurat-data')
# library(SeuratData)
# InstallData("panc8")
# AvailableData()

# data('panc8')




#modified for CTB with dark red colors
color_snap_mod1 = c('1'='#777711','2'='#E31A1C','3'='#68228B','4'='#771122','5'='grey','6'='#1F78B4','7'='#FFD700','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')

color_gradient_my <- c(
    rgb(5,48,97,maxColorValue = 255),
    rgb(42,113,178,maxColorValue = 255),
    rgb(147,198,222,maxColorValue = 255),
    rgb(239,243,245,maxColorValue = 255),
    rgb(253,219,199,maxColorValue = 255),
    rgb(214,96,77,maxColorValue = 255),
    rgb(121,5,34,maxColorValue = 255)

)

color_cellranger <-c('#820610','#C50F1E','#F42428','#F86D30','#FBB33D','#FCFB4E','#C0FB61','#87FB8F','#41FEF9','#2BAED7','#155CB1','#08238D') #for depth


color_signac = c(
'0'='#E6D55E','1'='#792B8A','2'='#DA703D','3'='#9DC8E5','4'='#BA273C','5'='#C2C184','6'='#7F8084','7'='#65AB53','8'='#D082AF','9'='#496EAB','10'='#DE896D','11'='#491F8B','12'='#E1AD49','13'='#8E1B85','14'='#E7EE77','15'='#7D1A1D','16'='#96B355')
names(color_signac) <- NULL

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

color_use <- colorRampPalette(BuenColors::jdb_palette('brewer_red'))(10)#(256)
color_use1 <- c('grey95',colorRampPalette(BuenColors::jdb_palette('brewer_red'))(256) )
color_use2 <- c('grey95',colorRampPalette(BuenColors::jdb_palette('brewer_blue'))(256) )

barplot(1:length(color_use),col = color_use,cex.main=2)


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




color_set <- color_list_archr[['grove']]
barplot(1:length(color_set),col=color_set )


##1
color_set <- c('11'='#1a1334',
               '9'='#01545a',
               '1'='#017351',
               '6'='#03c383',
               '8'='#aad962',
               '2'='#fbbf45',
               '10'='#ef6a32',
               '3'='#ed0345',
               '7'='#a12a5e',
               '5'='#710162',
               '4'='#3B9AB2')

##2
colorset_go <- rev(color_set_yellowbrick.flat[['YlGnBu.9']])
colorset_pathway <- rev(color_set_yellowbrick.flat[['RdPu.9']])
barplot(1:length(colorset_go),col=colorset_go )
barplot(1:length(colorset_pathway),col=colorset_pathway )

color_set <- c(colorset_go,rev(colorset_pathway))
barplot(1:length(color_set),col=color_set )

color_set <- c(
    
    '5' = '#081d58',
    '4' = '#253494',
    '2' = '#225ea8',
    '8' = '#1d91c0',
    '11' = '#FCC140',#fa9fb5
    '10' = '#fa9fb5',#'#dd3497',
    '7' = '#ae017e',
    '1' = '#f768a1',
    '3' = '#ae017e',
    '9' = '#49006a',
    '6' = '#7a0177'
    #'9' = '#41b6c4',
    #'8' = '#7fcdbb',
    #'2' = '#c7e9b4',
    #'7' = '#edf8b1',
    #'1' = '#ffffd9',
    #'3' = '#fff7f3',
    #'11' = '#fde0dd',
    #'12' = '#fcc5c0',
    
    

  


)


# color_set <- c(
#     '6' = '#081d58',
#     '5' = '#253494',
#     '4' = '#225ea8',
#     '9' = '#1d91c0',
#     #'9' = '#41b6c4',
#     #'8' = '#7fcdbb',
#     #'2' = '#c7e9b4',
#     #'7' = '#edf8b1',
#     #'1' = '#ffffd9',
#     #'3' = '#fff7f3',
#     #'11' = '#fde0dd',
#     #'12' = '#fcc5c0',
#     '8' = '#fa9fb5',
#     '2' = '#f768a1',
#     '10' = '#FCC140',#'#dd3497',
#     '7' = '#ae017e',
#     '1' = '#7a0177',
#     '3' = '#49006a'


# )


##3
map_cellcolor_cca <- list(
    '5' = '#710162',
  '7' = '#67001f',
 '3'= '#b2182b',
 '10'= '#FCC140',
 '2'= '#f4a582',
 '8'= '#FB8D3C',
 '1'= '#d6604d',
 '4'='#d6604d', 
 '6'= '#92c5de',
 '0'= '#4393c3',
 '9'= 'darkgreen'
 #'9'= '#053061'
)


library(magrittr)

source("quickDimPlot_labelon.r")
#quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', title= 'in vivo vs in vitro ',color_use = map_cellcolor_cca,shrink.x = .8, shrink.y = .2,shuffle = FALSE,pt.size = .1)

quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', title= 'in vivo vs in vitro ',color_use = color_set,shrink.x = .8, shrink.y = .2,shuffle = FALSE,pt.size = .1)

#quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', title= 'in vivo vs in vitro ',color_use = color_set,shrink.x = .8, shrink.y = .2,shuffle = TRUE,pt.size = .3)



#ggsave('pdfs/PLA-in-vivo-in-vitro.UMAP.coldpurple.pdf',height = 7.5, width = 8.5)

ggsave('pdfs/PLA-in-vivo-in-vitro.UMAP.coldpurple.new.pdf',height = 7.5, width = 8.5)



sample <- 'integration-invivo-vitroSTB'

#sample <- 'placenta(yawei ss2)'

#sample <- 'TS_cluster1_integrate_to_Tanglab'

##################################
#    prepare seurat objects      #
##################################



##STB in vivo

#early 
placenta.villus.early <- readRDS('data/placenta.cstb.rds') ##need to split as samples??
29132 x 23702 

#late
placenta.villus.late <- readRDS('data/PLA-late-combine-RNA.cstb.rds')
28686 x 23981 

##CT30 TS-derived STB
STB.CT30 <- readRDS('data/placenta.CT30.final.rds')
28264 x 5135 

#STB.CT30 <- readRDS('data/placenta.CT30.filter.rds')
#28264 x 5651 

cluster.df.add_ct30 <- readRDS('data/cluster.df.add.CT30.final.rds')
5135 × 16

all.equal( rownames(cluster.df.add_ct30),  colnames(STB.CT30)) #TRUE
all.equal( cluster.df.add_ct30$cluster, Idents(STB.CT30) ,check.attributes = FALSE) #TRUE



##BL TS-derived STB
#STB.BL <- readRDS('data/placenta.GRPTS.filter.filter.reclusterc4.rds')
STB.BL <- readRDS('data/placenta.GRPTS.final.rds')
27138 x 4474


cluster.df.add_grpts <- readRDS('data/cluster.df.add.GRPTS.final.rds')
4474 × 23

all.equal( rownames(cluster.df.add_grpts),  colnames(STB.BL)) #TRUE
all.equal( cluster.df.add_grpts$cluster, Idents(STB.BL) ,check.attributes = FALSE) #TRUE




##TS of CT30
#TS.CT30 <- readRDS('data/placenta.TS-CT30.filter.rds') #not ok
TS.CT30 <- readRDS('data/placenta.hTSC-CT30.rds')

##TS of BL

#TS.BL <- readRDS('data/placenta.TS-BL.filter.rds')
TS.BL <- readRDS('data/placenta.hTSC-BL.rds')

##TS of BT1

#TS.BT1 <- readRDS('data/placenta.TS-BT1.filter.rds')
TS.BT1 <- readRDS('data/placenta.hTSC-BT1.rds')



##dimplot each


options(repr.plot.height=7.5,repr.plot.width=8.5)

DimPlot(object = placenta.villus.early, reduction = "umap_rotate" ,label = TRUE,label.size = 12)
DimPlot(object = placenta.villus.late, reduction = "umap_rotate" ,label = TRUE,label.size = 12)


DimPlot(object = STB.CT30, reduction = "umap" ,label = TRUE,label.size = 8)
DimPlot(object = STB.BL, reduction = "umap_rotate" ,label = TRUE,label.size = 8)


DimPlot(object = TS.CT30, reduction = "umap" ,label = TRUE,label.size = 6)
DimPlot(object = TS.BL, reduction = "umap" ,label = TRUE,label.size = 6)
DimPlot(object = TS.BT1, reduction = "umap" ,label = TRUE,label.size = 6)








##############get average exprmat.z for DEG/TF/hormone heatmap comparison#########


exprMat.ave.z.villus.early <-  AverageExpression(placenta.villus.early,slot = "scale.data")[[1]] #29132 × 11
exprMat.ave.z.villus.late <-  AverageExpression(placenta.villus.late,slot = "scale.data")[[1]]

exprMat.ave.z.CT30 <-  AverageExpression(STB.CT30,slot = "scale.data")[[1]] #28264 × 12
exprMat.ave.z.BL <-  AverageExpression(STB.BL,slot = "scale.data")[[1]] #27429 × 12

deg.c10c3 <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.seurat_harmony/DEGs/placenta.markers.loose.c10c3.rds')



table(deg.c10c3$avg_logFC > 0) #537
FALSE  TRUE 
  216   321

deg.c10c3.sort <- deg.c10c3[order(deg.c10c3$avg_logFC,decreasing = TRUE),]

deg.c10 <- rownames(deg.c10c3.sort[deg.c10c3.sort$avg_logFC >0,] )
'SLC26A7''PAPPA''LINC01483''ADAMTSL1''ZNF117''SCIN''AC022872.1''LINC00578''ABCB1''AC022146.2''ADCY7''AUTS2''LAMA3''ISM1''FBN2''SLC27A6''EPHA1-AS1''KIF6''CATSPERB''ANK3''NAALADL2''SVEP1''AC093817.2''RALGAPA2''IQGAP2''STARD13''TMEM108''KMO''AC004704.1''RNF103-CHMP3''ATP8A2''AC009126.1''APBB2''CDYL2''SLC45A4''NRK''TRPV6''BCKDHB''PTPRG''NID1''NPAS2''RRAS2''CSHL1''ARHGAP17''CLRN1-AS1''ARHGAP26''CSH1''CCSER1''ADGRL3''THSD4''KATNBL1''AC069277.1''POSTN''LINC00882''ZC2HC1A''IL1RAPL2''AC016831.7''AC119674.1''PCDH11X''CERS4''LNX1''ST5''MYO9A''FER1L5''RAB17''GHR''LSS''MORC4''FBLN1''TGFBR3''GABRE''GLCCI1''AC079160.1''NMT2''AL162718.1''KIAA1211''LINC02860''XDH''CPS1''PTPRJ''ITPK1''CLIC5''ARHGAP42''NLRC5''GPC3''BCL2''CSH2''DTNB''DDX5''MIR100HG''KIAA0319''PAG1''PLXDC2''LRCH1''SLC15A1''TCHH''CBR3-AS1''OTOGL''PLCB1''ELMO1''LINC00474''CPXM2''GNG7''GNG12''NFU1''MICU3''PTPRM''STS''LINC02365''PHYHIPL''PLPP3''MKLN1''SGPP2''IMMP2L''AIG1''PCDH11Y''SCMH1''MAPK8''FNDC3A''LINC01091''TWSG1''ERV3-1''DOCK1''MPP5''ICA1''STAT4''MGAT5''SLCO2A1''MB21D2''STON2''TGM2''RIN2''LINC01807''ABCD3''NMNAT2''SLC25A35''PTPRQ''TXNRD1''TTBK2''KLRD1''TMEM164''KLHL5''SIPA1L1''FOXP1''WWC1''CLDN14''AHR''PDE4D''MTSS1''COMMD1''LUZP1''AP003181.1''AR''AC005532.1''SEMA6A''SPRY4-AS1''SUPT3H''PRKCE''ACSS1''RIMKLB''LINC02109''SLC4A4''AP000331.1''FRMD6-AS2''SGSM1''MTMR7''LINC01456''GRIP1''SLC23A2''LINC01949''FHIT''COBL''LINC01119''GH2''RAB3B''GNGT1''AL136962.1''ST3GAL3''PIP5K1B''SLC39A8''SLC7A2''AC073569.2''GREM2''ANO10''AC011287.1''LINC02267''HSD3B1''AC099520.1''WWOX''LRP8''ADAM10''XPO6''ANGPT2''IKZF2''TRIM5''AL138828.1''SSH2''ABCG2''SH3TC2-DT''DDX17''HIBADH''CCDC171''GRM7''GRAMD2B''GSTA4''TC2N''RNGTT''CARD18''AC013287.1''TTC28''PLEKHG1''GULP1''ANKDD1A''AL023574.1''COX10-AS1''COLEC12''PRKCA''CYTH3''THSD7A''PPP2R3A''LYST''CASTOR3''GSK3B''RANBP17''AGAP1''ZNF554''PTPN1''CHST11''CFAP299''AL365295.1''POU6F2''AL589787.1''NEK7''BCR''FAM234B''ZHX2''HSD11B2''GRK5''AC139718.1''MAP3K4''RNF217''TRAPPC9''NFIA''CAST''MAP3K5''COMT''PBX1''PLPP1''OGFRL1''TULP4''SPIDR''TBC1D5''TNS1''SYT12''EXT1''SLC19A3''KIAA1671''CHKA''MED12L''FAM126A''LITAF''AL589740.1''GLUL''ACACA''ARNTL2''GAB2''TGFBR2''MYO1B''CSGALNACT1''ARRB1''ATF7IP2''FAM172A''TWIST1''DLGAP4''PLEKHA7''AMPD3''BPGM''MIR4435-2HG''MAP4K3''NCF1''MAP4K3-DT''CTDSP1''PSG2''TIAM1''ADORA1''TBC1D30''ATXN7L1''HLCS''AC105411.1''TBC1D1''MGST1''NHSL1''LIMS1''LINC00278''LNPEP''CLMN''ARHGAP10''PTH2R''ZDHHC14''AC079760.1''TJP1''AC006378.1''SUCLG2''TMEM120B''GAB1''AC087854.1''ARFGEF3''PPP1R12B''DENND1A''USP46''SRSF6''MAGI1''TSGA10''PDLIM5''SELENOI''HMGB1''NCOA2''FAM171B''AMFR''RERG''AC009264.1'

deg.c3 <- rev(rownames(deg.c10c3.sort[deg.c10c3.sort$avg_logFC <0,] ))
'MYCNUT''EGLN3''COL27A1''FSTL3''MME''ARHGAP24''XACT''MIR193BHG''PFKP''BNIP3L''SASH1''FLT1''FAM83D''GPC5''GASK1A''SH3BP5''LIMD1''SLC2A1''SERPINE1''SLC6A8''HDAC2''IL1RAP''ARID3A''MTHFD1L''ANKRD37''PLCL2''ANXA1''MYCN''SIGLEC6''MTUS1''ATP10D''RASSF8''NPEPPS''NECTIN4''SLC2A3''AL513164.1''C4orf47''TBL1XR1''TCL6''PFKFB4''NABP1''MYOF''AC108690.1''ERRFI1''MTSS2''SYNJ2''AC087857.1''FRMD4A''TNFRSF1B''AC020916.1''ZNF331''GATA2''GMCL1''MFSD12''PLIN2''ITPRID2''LPCAT1''ATP2C2''FLNB''GPRC5C''FOCAD''LVRN''EZR''NDRG1''TDRP''NEURL1''MIR29B2CHG''PLEKHM2''MALAT1''GSE1''ADAMTS20''SEMA3B''CMIP''HILPDA''ZNF626''FAM153CP''MGAT3''CORO1C''MAPK8IP3''NEK11''BASP1''HPCAL1''JUP''ARHGAP45''EGLN1''MYCNOS''ARHGAP30''ERO1A''UBC''KRBOX1''INPP5B''FAM43A''AC026167.1''CSNK1E''MID1''DIAPH2''P4HA1''JPT2''SPATA13''TPP2''PNCK''FAM120A''PTPRR''GALNT2''HIST1H2AC''AHNAK''ESYT2''MLLT10''NEAT1''DHX35''ZFP36L1''RABL6''CCDC138''INHBA''KMT2E''CGB3''SHANK2''TGFB1''HIF1A-AS3''CHD2''CALM1''PPP1R12C''HIST1H4H''TTC3''RAB8B''GPR78''ACOT11''MFSD2A''GDPD4''ZMYND8''HMG20B''MXD1''CYP11A1''DDB1''NCOR2''EBP''SLC7A5''SP3''GTF2IRD1''SEMA7A''ITGA5''HTRA4''SH3PXD2A''UCA1''LHB''PCNP''BCL6''NRIP1''FOXJ3''ERGIC1''UBASH3B''KIAA1217''HK2''DLEU2''DGKI''PRDX6''CBLB''DGKZ''TREML2''TPD52''TRIM14''MXI1''AC005786.3''AL645568.1''FHL2''POF1B''PHIP''NGLY1''WDR60''TTLL5''CCNY''TMEM91''MIR503HG''RAB3GAP1''EGFR''DENND4A''FMN1''PHKA2''PLOD2''AC018754.1''PKM''GMDS''PHACTR2''PIM3''FZR1''EFHD1''REEP3''CAPN7''AFAP1''PRKD3''DNAJB6''SCARB1''POMP''SLC16A3''TPBG''AC011453.1''TET3''TNS3''LINC02832''RIPK2''ZNF175''PICALM''MTMR4''SERTAD2''PTDSS1''HES2''INHA''ZCCHC2''TET1''GATA3''PCOLCE2''PORCN''RNMT''DOCK5''SLC1A6''ANGPTL4'


###prepare dataframe to plot
deg.c10.bk <- deg.c10
deg.c3.bk <- deg.c3


deg.c10 <- deg.c10[1:20]
deg.c3 <- deg.c3[1:20]

table(deg.c10 %in% rownames(exprMat.ave.z.villus)) #321 TRUE
table(deg.c3 %in% rownames(exprMat.ave.z.villus)) #216 TRUE

table(deg.c10 %in% rownames(exprMat.ave.z.CT30)) #319 TRUE 2 FALSE
table(deg.c3 %in% rownames(exprMat.ave.z.CT30)) #216 TRUE


table(deg.c10 %in% rownames(exprMat.ave.z.BL)) #315 TRUE 6 FALSE
table(deg.c3 %in% rownames(exprMat.ave.z.BL)) #216 TRUE


exprMat.ave.z.villus.sel <- exprMat.ave.z.villus[c(deg.c10,deg.c3),c('10','3')]
all.equal(rownames(exprMat.ave.z.villus.sel) ,c(deg.c10,deg.c3)) #TRUE
colnames(exprMat.ave.z.villus.sel) <- c('villus_Mature1','villus_Mature2')


exprMat.ave.z.CT30.sel <- exprMat.ave.z.CT30[c(deg.c10,deg.c3),c('recluster_3','recluster_1')]
all.equal(rownames(exprMat.ave.z.CT30.sel) ,c(deg.c10,deg.c3)) #two string mismatch
idx <- which(rownames(exprMat.ave.z.CT30.sel) != c(deg.c10,deg.c3))
rownames(exprMat.ave.z.CT30.sel)[idx] <- c(deg.c10,deg.c3)[idx]
colnames(exprMat.ave.z.CT30.sel) <- c('CT30_Mature1','CT30_Mature2')


exprMat.ave.z.BL.sel <- exprMat.ave.z.BL[c(deg.c10,deg.c3),c('recluster_2','recluster_5')]
all.equal(rownames(exprMat.ave.z.BL.sel) ,c(deg.c10,deg.c3)) #six string mismatch
idx <- which(rownames(exprMat.ave.z.BL.sel) != c(deg.c10,deg.c3))
#colnames(exprMat.ave.z.BL.sel) <- paste0('BL_',colnames(exprMat.ave.z.BL.sel))
colnames(exprMat.ave.z.BL.sel) <- c('BL_Mature1?','BL_Mature2')

all.equal(rownames(exprMat.ave.z.villus.sel) ,c(deg.c10,deg.c3))
all.equal(rownames(exprMat.ave.z.CT30.sel) ,c(deg.c10,deg.c3))
all.equal(rownames(exprMat.ave.z.BL.sel) ,c(deg.c10,deg.c3))
#TRUE


tflist_rowid <- readRDS('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.seurat_harmony/tflist_row_order.rds')


table(tflist_rowid %in% rownames(exprMat.ave.z.villus) )

TRUE 
 194

table(tflist_rowid %in% rownames(exprMat.ave.z.CT30) )

TRUE 
 194 

table(tflist_rowid %in% rownames(exprMat.ave.z.BL) )

TRUE 
 194



exprMat.ave.z.villus.sel <- exprMat.ave.z.villus[tflist_rowid,c('7','9','6','11','8','1','3','10','2','4','5')]
all.equal(rownames(exprMat.ave.z.villus.sel) ,tflist_rowid) #TRUE
#colnames(exprMat.ave.z.villus.sel) <- c('villus_Mature1','villus_Mature2')

colnames(exprMat.ave.z.villus.sel) <- paste0('villus_',colnames(exprMat.ave.z.villus.sel))

exprMat.ave.z.CT30.sel <- exprMat.ave.z.CT30[tflist_rowid,c('6','2','0','3','recluster_1','recluster_3','recluster_2','recluster_0','recluster_4','7','8')]
all.equal(rownames(exprMat.ave.z.CT30.sel) ,tflist_rowid) #TRUE

colnames(exprMat.ave.z.CT30.sel) <- paste0('CT30_',colnames(exprMat.ave.z.CT30.sel))


#idx <- which(rownames(exprMat.ave.z.CT30.sel) != c(deg.c10,deg.c3))
#rownames(exprMat.ave.z.CT30.sel)[idx] <- c(deg.c10,deg.c3)[idx]
#colnames(exprMat.ave.z.CT30.sel) <- c('CT30_Mature1','CT30_Mature2')


exprMat.ave.z.BL.sel <- exprMat.ave.z.BL[tflist_rowid,c('7','2','1','6','recluster_2','recluster_3','recluster_5','recluster_4','recluster_6','recluster_1','recluster_0','recluster_7')]
all.equal(rownames(exprMat.ave.z.BL.sel) ,tflist_rowid) #TRUE
#idx <- which(rownames(exprMat.ave.z.BL.sel) != c(deg.c10,deg.c3))
colnames(exprMat.ave.z.BL.sel) <- paste0('BL_',colnames(exprMat.ave.z.BL.sel))
#colnames(exprMat.ave.z.BL.sel) <- c('BL_Mature1?','BL_Mature2')

all.equal(rownames(exprMat.ave.z.villus.sel) ,tflist_rowid)
all.equal(rownames(exprMat.ave.z.CT30.sel) ,tflist_rowid)
all.equal(rownames(exprMat.ave.z.BL.sel) ,tflist_rowid)


exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel)
title <- 'villus only'

exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel,exprMat.ave.z.CT30.sel)#,exprMat.ave.z.BL.sel)

title <- 'villus vs CT30'

exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel,exprMat.ave.z.BL.sel)

title <- 'villus vs BL'

exprMat.ave.z.sel.combine <- cbind.data.frame(exprMat.ave.z.villus.sel,exprMat.ave.z.CT30.sel,exprMat.ave.z.BL.sel)#,exprMat.ave.z.BL.sel)

title <- 'villus vs CT30 vs BL'

###split plotting & combine plotting


# options(repr.plot.width = 3.5, repr.plot.height = 10)
# res.p <- pheatmap::pheatmap(
#                    #exprMat.ave.z.villus.sel, 
#                    #exprMat.ave.z.CT30.sel,
#                    #exprMat.ave.z.BL.sel,
#                    exprMat.ave.z.sel.combine,
#                    cluster_cols = FALSE,
#                    cluster_rows = FALSE, 
#                    treeheight_col = 2.5,
#                    treeheight_row = 2.5,
#                    #cellwidth = 20, 
#                    #cellheight = 2.8,
#                    na_col = 'white', 
#                    color = color_use,
#                    border =TRUE,
#                    border_color = 'black',
#                    #labels_col = colid, 
#     #                labels_row = make_bold_names(marker.mat.z,
#     #                                             rownames, 
#     #                                             rowid_hi
#     #                                            ),
#                    show_rownames = TRUE,

#                    angle_col = 270,
#                    #display_numbers = data.df.text,
#                    #number_color = 'white',
#                    #fontsize_number = 10,
#                    fontsize = 13,
#                    fontsize_col = 13,
#                    fontsize_row = 10,
#                    main = paste('top20 DEGs compare\n',title,sep=''),
#                   # main = paste(title,' gene expression zscore',sep=''),
#                    silent = FALSE,
#                    #scale = 'none',
#                    scale = 'row',
#                    #scale = 'column',
#                    #legend_breaks = c(2,4,6,8,10),
#                    #legend_labels = c(2,4,6,8,10),
#                    legend = FALSE
#                    #legend = TRUE
#     #                        height = height,
#     #                        width = width,
#     #                        filename = paste('result_do_GO_Pathway_quickonestep/ck.GO_BP.pvalue_',
#     #                                         pvalue,'.qvalue_', qvalue,'.heatmap.pdf' ,sep = '' 
#     #                                        )

#                   )



# print(res.p)





###zscore row data
scale_rows = function(x){ #pheatmap code
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

#zscore by row
exprMat.ave.z.sel.combine.z <- scale_rows(exprMat.ave.z.sel.combine)


options(repr.plot.height=30,repr.plot.width=18)
#options(repr.plot.height=65,repr.plot.width=6)
hp = ComplexHeatmap::Heatmap(exprMat.ave.z.sel.combine.z, name = "TF gene expression", 
##hp = Heatmap(exprMat.ave.data.z, name = "exprMat.ave.data.z", 
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             show_row_names = TRUE,
             use_raster = FALSE,
             ##col = circlize::colorRamp2(seq(-1,1,by=2/10), viridis(n = 11,option = "C")),
             row_title = "hclust: ward.D",
             clustering_method_row = "ward.D", ##ward.D,complete
             #clustering_distance_rows  = "pearson",
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 18),
             #clustering_distance_columns  = function (m) dist(m,method="manhattan") #euclidean, manhattan,minkowski  #clustering_distance_columns  = "pearson",
             #heatmap_legend_param = list(color_bar = "continuous")
             ##right_annotation = ha#,heatmap_width=unit(8, "cm")
) 
#width = max(grobWidth(textGrob(labels))))

##pdf('TF_gene_expression.heatmap.pdf',height=15,width=4,useDingbats=FALSE)
set.seed(1);draw(hp, heatmap_legend_side = "left",gap = unit(0.1, "cm")) #8.895833 8.843750




# #######################techmann###############
# techmann.tenx <- readRDS("/home/mjwang/pwdex/placenta_10X/03.seurate/FCA7196220/data.TE.rds") #a mat
# techmann.tenx <- Seurat::CreateSeuratObject(techmann.tenx, project='techmann.tenx', 
#                                            min.cells=0, min.features=0,
#                                            meta.data = readRDS("/home/mjwang/pwdex/placenta_10X/03.seurate/FCA7196220/metadata.TE.rds")
#                                           )


# #plotDR(obj = techmann.tenx,sample = 'techmann.tenx',feature = 'seurat_clusters',color = color_good)

# #####################PLA 8w tenx###############
# PLA.8w.tenx <- readRDS("data/scRNA-LSI-PLA-8w-RNA.TE.techmann.rds") #se object of 8w placenta
# PLA.8w.tenx <- Seurat::CreateSeuratObject(assay(PLA.8w.tenx), project='PLA.8w.tenx', 
#                                            min.cells=0, min.features=0,
#                                            meta.data = as.data.frame(colData(PLA.8w.tenx))
#                                           )

# #plotDR(obj = PLA.8w.tenx,sample = 'PLA-8w',feature = 'res.0.6',color = color_good)

# ####################PLA 8w and 24w ss2, from yawei paper################
# PLA.ss2 <- readRDS("data/TR599.rds") 
# # PLA.ss2.obj <- Seurat::CreateSeuratObject(PLA.ss2, project='PLA.ss2', 
# #                                            min.cells=0, min.features=0,
# #                                            meta.data = as.data.frame(cPLA.ss2))
# #                                           )

# PLA.ss2.count <- GetAssayData(PLA.ss2,slot = 'counts')



# ##read in the raw exprMat
# PLA.ss2.mat <- read.table('data/data-scRNASEQ-LiuYW/attachment1-single-cell-RNA-data.txt',header = TRUE,sep='\t') #20866 x 1806
# PLA.ss2.dr <- read.table('data/data-scRNASEQ-LiuYW//attachmetn2-tsne_df_All_0816.txt',header = TRUE,sep = '\t')
# #

# #plot(PLA.ss2.dr[,c('tSNE1','tSNE2')],pch = 19)


# cluster.df.add <- PLA.ss2.dr
# table(cluster.df.add$Type)

# HE24W_EVT  HE8W_CTB  HE8W_EVT  HE8W_STB  HE8W_STR 
#       193       239       430        58       551 

# HE8W_CTB + HE8W_EVT + HE8W_STB = 727

# table(cluster.df.add$Embryo)

#  E1  E2  E3  E4  E5  E6  E7  E8 
# 371  30  87  23 193 182 310 275 

# table(cluster.df.add$APcluster)
#  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17 
# 101  88  47 109 100  91  63 111  84  93  91  69 105  90 125  47  57 

# table(cluster.df.add$APcluster_Merge)
#   1   2   3   4   5   6   7   8  11  12  13  15  16  17 
# 101  88  47 109 100  91 246 195  91  69 105 125  47  57


# cluster.df.add$APcluster <- factor(cluster.df.add$APcluster, levels = 1:17 )
# cluster.df.add$APcluster_Merge <- factor(cluster.df.add$APcluster_Merge, levels = names(table(cluster.df.add$APcluster_Merge)) )



plotDR <- function(obj = NULL, sample = NULL,feature = NULL, color = NULL){

    ##get tne final cluster.df.add
    umap <- Embeddings(obj,reduction = 'umap')
    cluster <- Idents(obj)
    cluster.df = data.frame(cluster=cluster,UMAP_1=umap[,1],UMAP_2=umap[,2]) 
    rownames(cluster.df) = rownames(obj@meta.data)
    
    stopifnot(all.equal (rownames(cluster.df),rownames(obj@meta.data) ) )
    cluster.df.add <- cbind(cluster.df,obj@meta.data)

    colnames(cluster.df.add) <- make.unique( colnames(cluster.df.add) )

    stopifnot(all.equal(Idents(obj),cluster.df.add$cluster,check.attributes = FALSE) )


    ##label on cluster
    options(repr.plot.height=5,repr.plot.width=6.5,repr.plot.res = 150)
    res.p <- ggplot(cluster.df.add,aes_string(x="UMAP_1",y="UMAP_2",col=`feature`  )) +
    #ggplot(cluster.df.add,aes(x=tSNE1,y=tSNE2,col=APcluster_Merge  )) +
      geom_point(size = 0.6,show.legend = TRUE,alpha= 1 ) +

      scale_colour_manual(values = color)  +
      #scale_colour_manual(values = color_good)  +
      #scale_colour_manual(values = unlist(map_cellcolor) )  +
      ##scale_colour_manual(values = color_snap_mod1)  +
      #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
      #theme_classic() +
      #theme_bw() +
      theme(
            #legend.position = 'right',
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
      ggtitle(paste(sample, "total cells:",nrow(cluster.df.add),  sep=" ") ) +
    #   geom_text(data = centers_shift, #the halo
    #             mapping = aes(x=x,y=y,label = cluster), 
    #             colour = "white", 
    #             size = 4.5) +
    #             ##size = 6.5) +
    #   geom_text(data = centers, 
    #             mapping = aes(x=x,y=y,label = cellname), 
    #             ##mapping = aes(x=x,y=y,label = cluster), 
    #             colour = "black", 
    #             size = 4.5) +
    #             ##size = 6.5) +
       guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
    #   xlim(left,right) + ylim(bottom,top) +
      #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
      labs(x = "UMAP1", y = "UMAP2")
      #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

   return(res.p)
}


plotDR(obj = placenta.villus.early,sample = 'placenta.villus.early',feature = 'cluster',color = color_good)
plotDR(obj = placenta.villus.late,sample = 'placenta.villus.late',feature = 'cluster',color = color_good)

plotDR(obj = STB.CT30,sample = 'STB.CT30',feature = 'cluster',color = color_good)
plotDR(obj = STB.BL,sample = 'STB.BL',feature = 'cluster',color = color_good)
plotDR(obj = TS.CT30,sample = 'TS.CT30',feature = 'cluster',color = color_good)
plotDR(obj = TS.BL,sample = 'TS.BL',feature = 'cluster',color = color_good)
#plotDR(obj = TS.BT1,sample = 'TS.BT1',feature = 'cluster',color = color_good)

# plotDR(obj = PLA.ss2,sample = 'placenta_yawei',feature = 'sample',color = color_good)
# plotDR(obj = PLA.ss2,sample = 'placenta_yawei',feature = 'type',color = color_good)

#FeaturePlot(object = PLA.ss2,features = c('MKI67','DNMT1','PSG8','HLA-G') )
#FeaturePlot(object = PLA.ss2,features = c('ERVFRD-1') )

# #############TS line  from Yiming
# TS.ss2 <- readRDS("data/HTS630-OK.rds") #seurat object 

# # TS.ss2.obj <- Seurat::CreateSeuratObject(TS.ss2@assays$RNA@data, project='TS.ss2', 
# #                                            min.cells=0, min.features=0,
# #                                            meta.data = TS.ss2@meta.data
# #                                           )

# plotDR(obj = TS.ss2,sample = 'TS_yiming',feature = 'C5',color = color_good)

# options(repr.plot.height=15,repr.plot.width=15)
# FeaturePlot(object = TS.ss2,features = c('MKI67','DNMT1','PSG8','HLA-G') )
# FeaturePlot(object = TS.ss2,features = c('FLT1','ENG','PAPPA','LEP','CSHL1','STAT4','STAT5A','TP53','EPCAM') )



# ##############TE from TangFC lab in vitro data ss2
# TE.ss2 <- readRDS('data/TE3859.rds')


# plotDR(obj = TE.ss2,sample = 'TE_tanglab',feature = 'type',color = color_good)


##################################
#    integrate by CCA            #
##################################


table(placenta.villus.early[['sample']])
PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
          4959           3610           4418           3436           2799 
PLA-early6-RNA 
          4480

# PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
#           5069           3687           4489           3515           2846 
# PLA-early6-RNA 
#           4583


table(placenta.villus.late[['sample']])

PLA-late1-RNA PLA-late2-RNA PLA-late3-RNA PLA-late5-RNA PLA-late6-RNA 
         2563          3566          3720          3886          4545 
PLA-late9-RNA 
         5701


##split by sample (remove pca, neighbor slot??)


placenta.villus.early.list <- SplitObject(object = placenta.villus.early, split.by = 'sample')
'PLA-early1-RNA''PLA-early2-RNA''PLA-early3-RNA''PLA-early4-RNA''PLA-early5-RNA''PLA-early6-RNA'

placenta.villus.late.list <- SplitObject(object = placenta.villus.late, split.by = 'sample')
'PLA-late1-RNA''PLA-late2-RNA''PLA-late3-RNA''PLA-late5-RNA''PLA-late6-RNA''PLA-late9-RNA'

#placenta.list <- list( 'placenta.villus.early' = placenta.villus.early, 'STB.CT30' = STB.CT30 , 'STB.BL' = STB.BL)#, 'TS.CT30' = TS.CT30, 'TS.BL' = TS.BL)

#placenta.list <- list( 'placenta.villus' = placenta.villus, 'STB.CT30' = STB.CT30 , 'STB.BL' = STB.BL, 'TS.CT30' = TS.CT30, 'TS.BL' = TS.BL, 'TS.BT1' = TS.BT1 )

#placenta.list <- list( 'placenta.villus' = placenta.villus, 'STB.CT30' = STB.CT30 , 'STB.BL' = STB.BL)

#placenta.list <- list(  'STB.CT30' = STB.CT30 , 'STB.BL' = STB.BL )

placenta.list <- c(placenta.villus.early.list,list(  'STB.CT30' = STB.CT30, 'STB.BL' = STB.BL))



names(placenta.list)
'PLA-early1-RNA''PLA-early2-RNA''PLA-early3-RNA''PLA-early4-RNA''PLA-early5-RNA''PLA-early6-RNA''STB.CT30''STB.BL'

#'PLA-early1-RNA''PLA-early2-RNA''PLA-early3-RNA''PLA-early4-RNA''PLA-early5-RNA''PLA-early6-RNA''STB.CT30''STB.BL'

#placenta.list <- c(placenta.villus.list,list(placenta.villus.list,  'STB.CT30' = STB.CT30))
#placenta.list <- c(placenta.villus.list,list(placenta.villus.list[['PLA-early1-RNA']],  'STB.CT30' = STB.CT30))



#placenta.list <- list( 'TS.ss2' = TS.ss2, 'PLA.8w.tenx' = PLA.8w.tenx ,'PLA.ss2' = PLA.ss2 , 'TE.ss2' = TE.ss2)
#placenta.list <- list( 'TS.ss2' = TS.ss2, 'PLA.ss2' = PLA.ss2 , 'TE.ss2' = TE.ss2)

#placenta.list <- list( 'techmann' = techmann.obj, 'tenx' = tenx.obj ,'ss2' = ss2.obj  )


# ##preprocess: log transformation, find variable genes
# for (i in 1:length(placenta.list)) {
#     #placenta.list[[i]] <- NormalizeData(placenta.list[[i]], verbose = FALSE)
#     placenta.list[[i]] <- FindVariableFeatures(placenta.list[[i]], selection.method = "vst", 
#         nfeatures = 2000, verbose = FALSE)
# }



#############part 1, integration for a reference (not so good)########
#reference.list <- placenta.list[c("tenx", "techmann")]
#reference.list <- placenta.list[c("tenx", "ss2")]
#reference.list <- placenta.list[c("tenx", "techmann","ss2")]
#reference.list <- placenta.list[c("techmann", "ss2")]

reference.list <- placenta.list
placenta.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30) #multiple object list, anchor.feature 2000 by defaut


##8 sample new 

# Warning message in CheckDuplicateCellNames(object.list = object.list):
# “Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.”
# Computing 2000 integration features

# Scaling features for provided objects

# Finding all pairwise anchors

# Running CCA

# Merging objects

# Finding neighborhoods

# Finding anchors

# 	Found 12579 anchors

# Filtering anchors

# 	Retained 4531 anchors

# Running CCA

# Merging objects

# Finding neighborhoods

# Finding anchors

# 	Found 13576 anchors

# Filtering anchors



###8 sample old BL
# Scaling features for provided objects

# Finding all pairwise anchors

# Running CCA

# Merging objects

# Finding neighborhoods

# Finding anchors

# 	Found 23776 anchors

# Filtering anchors

#	Retained 3034 anchors

#...many times


saveRDS(placenta.anchors,'placenta.anchors.8.samples.rds')


placenta.integrated <- IntegrateData(anchorset = placenta.anchors, dims = 1:30)
# Warning message:
# “Adding a command log without an assay associated with it”







##########rerun this step (will get integrate data with normalization method LogNormalize or SCT, not used)######
#https://github.com/satijalab/seurat/issues/1447
all.genes <- rownames(placenta.integrated@assays$RNA@data) #31478 #rownames(placenta.integrated)
placenta.integrated.full <- IntegrateData(anchorset = placenta.anchors, dims = 1:30,features.to.integrate = all.genes)

##slow and mem big



DefaultAssay(placenta.integrated.full) <- "integrated"

# Run the standard workflow for visualization and clustering

all.genes <- rownames(placenta.integrated.full)
placenta.integrated.full <- ScaleData(placenta.integrated.full, verbose = FALSE,features = all.genes) #from 2000 gene of data slot
placenta.integrated.full <- RunPCA(placenta.integrated.full, npcs = 30, verbose = FALSE)
placenta.integrated.full <- RunUMAP(placenta.integrated.full, reduction = "pca", dims = 1:30)


placenta.integrated.full <- AddMetaData(placenta.integrated.full,metadata = Idents(placenta.integrated.full), col.name = 'cluster.ori')

placenta.integrated.full <- FindNeighbors(object = placenta.integrated.full, reduction = 'pca', dims = 1:30)
placenta.integrated.full<- FindClusters(object = placenta.integrated.full, verbose = FALSE, algorithm = 1,resolution=0.9)


options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(object = placenta.integrated.full, label = TRUE,cols=color_good,pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() #+ xlim(-10,10) + ylim(-10,6)


##copy placenta.integrated umap and cluster and metatable

DefaultAssay(placenta.integrated.full)
DefaultAssay(placenta.integrated)

all.equal(colnames(placenta.integrated.full),colnames(placenta.integrated) ) #no

table(colnames(placenta.integrated)  %in% colnames(placenta.integrated.full) )
 TRUE 
31369

table(rownames(placenta.integrated)  %in% rownames(placenta.integrated.full) )
TRUE
2000


#saveRDS(placenta.integrated.full,"placenta.integrated.full.rds")


##output slot data
exprMat.integrated_full.data <- GetAssayData(object = placenta.integrated.full,assay = 'integrated',slot = 'data')
25150 x 33827

exprMat.integrated_full.scale.data <- GetAssayData(object = placenta.integrated.full,assay = 'integrated',slot = 'scale.data')

dim(exprMat.integrated_full.scale.data)
25150 x 33827

saveRDS(exprMat.integrated_full.data,'exprMat.integrated_full.data.rds')
saveRDS(exprMat.integrated_full.scale.data, 'exprMat.integrated_full.scale.data.rds')


#########rerun done################







###plot the integrated object
library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(placenta.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
placenta.integrated <- ScaleData(placenta.integrated, verbose = FALSE) #from 2000 gene of data slot
placenta.integrated <- RunPCA(placenta.integrated, npcs = 30, verbose = FALSE)
placenta.integrated <- RunUMAP(placenta.integrated, reduction = "pca", dims = 1:30)



##tuning


#######tuning UMAP embedding#####

seed.use = 123
min.dist = 0.5
spread=1.2


# random_state = 0
# n_comps = 2
# min.dist = 0.
# spread = 1.0

umap.bk <- Embeddings(placenta.integrated,reduction = 'umap')

for(seed.use in c(123,sample( ) ) ){
    for (min.dist in c(0.1,0.2,0.3,0.4,0.5)){
    #for mdist in [0.3]:
        for (spread in c(0.5,1.0,1.5,2.0)){
    #    for spread in [1.0]:
          placenta.integrated = RunUMAP(placenta.integrated,dims=1:30,
                             seed.use = seed.use,
                             min.dist = min.dist,
                             spread=spread,
                             return.model=FALSE,
                             verbose=0
                            ) 
         options(repr.plot.height=7.5,repr.plot.width=7.5)
         res.p <- DimPlot(object = placenta.integrated, label = TRUE,cols=color_good,pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() + ggtitle(paste( c('umap tunning: mdist ',min.dist,' spread ',spread,' seed.use ',seed.use) ,collapse =' '))
         print(res.p)

         res.p <- FeaturePlot(placenta.integrated, features = 'FLT1', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 
         print(res.p)
         res.p <- FeaturePlot(placenta.integrated, features = 'PAPPA', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 
         print(res.p)
         res.p <- FeaturePlot(placenta.integrated, features = 'SH3TC2', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 
         print(res.p)
         res.p <- FeaturePlot(placenta.integrated, features = 'LAMA3', reduction = "umap",pt.size = .1,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 
         print(res.p)


       }
    }

}



##fix umap?


seed.use = 123
min.dist = 0.5
spread= 0.5

placenta.integrated = RunUMAP(placenta.integrated,dims=1:30,
                 seed.use = seed.use,
                 min.dist = min.dist,
                 spread=spread,
                 return.model=TRUE,
                 verbose=0
                ) 


options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(object = placenta.integrated, label = TRUE,cols=color_good,pt.size = 0.5,label.size = 8,reduction = "umap") + NoLegend() + xlim(-10,10) + ylim(-10,6)



#########rotate UMAP###
umap.df <- data.frame(UMAP_1=Embeddings(placenta.integrated,reduction = 'umap')[,1],UMAP_2=Embeddings(placenta.integrated,reduction = 'umap')[,2])
#rownames(umap.df) <- colnames(placenta)

#iterative rotate the umap direction
#for(degree in seq(from = 50,90,5) ) { 
#for(degree in seq(from = 65,90,2) ) { 
for(degree in c(55) ) { 
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

fix degree = 55
#fix degree = 175
#fix degree = 150
##fix degree = 175


#umap.df[,1] <- -1 * umap.df[,1]
#UMAP.rotate <- umap.df


plot(UMAP.rotate)

#rewrite the umap slot
rownames(UMAP.rotate) <- rownames(Embeddings(placenta.integrated,reduction = 'umap'))
#UMAP.bk = Embeddings(placenta,reduction = 'umap')
#Embeddings(placenta,reduction = 'umap_rotate') = as.matrix(UMAP.rotate) #error


placenta.integrated[['umap_rotate']] <- CreateDimReducObject(embeddings = as.matrix(UMAP.rotate), key = 'umaprotate_', assay = 'integrated')

#saveRDS(object=UMAP.rotate,file = 'uwot.umap.seed177.rotate175.rds')

options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(placenta.integrated, reduction = "umap_rotate",label=TRUE,cols=color_good[],label.size = 5,pt.size = .1) + NoLegend() + xlim(-10,10) + ylim(-12,10)




###do clustering

#backup original clusterid

placenta.integrated <- AddMetaData(placenta.integrated,metadata = Idents(placenta.integrated), col.name = 'cluster.ori')


placenta.integrated <- FindNeighbors(object = placenta.integrated, reduction = 'pca', dims = 1:30)
placenta.integrated<- FindClusters(object = placenta.integrated, verbose = FALSE, algorithm = 1,resolution=0.9)

placenta.integrated<- FindClusters(object = placenta.integrated, verbose = FALSE, algorithm = 1,resolution=0.3)



#####tuning clustering###########

res.p.list <- list()

DefaultAssay(placenta.integrated) <- "integrated"

#louvain , leiden (slow)
for (res in seq(0.6,1.2,0.1) ){
#for (res in seq(0.1,0.5,0.1) ){
    placenta.integrated<- FindClusters(object = placenta.integrated, verbose = FALSE, algorithm = 1,resolution=res) #louvain
    #placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 4,resolution=res) #slow, mem big!
    
    options(repr.plot.height=7.5,repr.plot.width=7.5)
    res.p <- DimPlot(object = placenta.integrated, label = TRUE,cols=c(color_good,color_good),pt.size = 0.5,label.size = 8,reduction = "umap_rotate") + NoLegend() + ggtitle(paste( c('resolution ',res) ,collapse =' ')) #+ ylim(-10,10)
     print(res.p)
     res.p.list[[paste( c('res',res) ,collapse ='')]] <- res.p
    
}

for(i in names(res.p.list)){print(res.p.list[[i]]) } #+ ylim(-8,9) + xlim(-5,10))}



#fix cluster with res
res = 0.8#0.6


#placenta.integrated <- FindNeighbors(placenta.integrated, dims = 1:40) 
placenta.integrated <- FindClusters(object = placenta.integrated, verbose = FALSE, algorithm = 1,resolution=res)


DimPlot(placenta.integrated, reduction = "umap_rotate",label=TRUE,cols=color_good,label.size = 8,pt.size = .1,shuffle = TRUE) + NoLegend() + xlim(-10,10) + ylim(-12,10)



#saveRDS(placenta.integrated,"placenta.integrated.tuning.before_filter.rds")
#33478 x 33827

saveRDS(placenta.integrated,"placenta.integrated.reclustering_res_0.8.rds")
33478 x 31369

###get integration cluster.df.add only#####


cluster <- Idents(placenta.integrated)
#umap <- placenta@reductions$umap@cell.embeddings
umap <- placenta.integrated@reductions$umap_rotate@cell.embeddings

#umap <- placenta@reductions$umap_rotate@cell.embeddings

colnames(umap) <- c('UMAP_1','UMAP_2')

all.equal(names(cluster),rownames(umap)) #TRUE


cluster.df <- data.frame(cluster=cluster,umap)

metadata <- placenta.integrated@meta.data
all.equal(rownames(cluster.df),rownames(metadata))#TRUE

intersect(colnames(cluster.df), colnames(metadata) ) #cluster

grep('cluster$',colnames(metadata),value=FALSE )
49

#12 28
#49,50
#13,16

metadata[,c(49)] <- NULL


###
all.equal(metadata,placenta.integrated@meta.data[,-c(49)]) #TRUE

##update cluster integration in saved cluster.df.add$cluster
all.equal(rownames(cluster.df.add),rownames(cluster.df) ) #TRUE

cluster.df.add$cluster <- cluster.df$cluster
###

cluster.df.add <- cbind(cluster.df, metadata)

table(cluster.df.add$cluster)
 0    1    2    3    4    5    6    7    8    9   10   11   12 
5032 3550 3535 2924 2407 2356 2249 2024 1836 1523 1393 1361 1179

#  0    1    2    3    4    5    6    7    8    9   10 
# 6401 5529 4847 3797 2905 2497 2131 2037 1394 1184 1105



idx <- grep("RNA_snn_res",colnames(cluster.df.add) )
cluster.df.add[,colnames(cluster.df.add)[idx]] <- NULL


library(magrittr)

source("quickDimPlot_labelon.r")
#quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', title= 'in vivo vitro ',color_use = map_cellcolor_cca,shrink.x = .8, shrink.y = .2,shuffle = FALSE,pt.size = .3)

quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', title= 'in vivo vitro ',color_use = color_set,shrink.x = .8, shrink.y = .2,shuffle = FALSE,pt.size = .3)




###dotdistri and merge?


par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('4','3','1','11','12','9','0','6','10','2','7','8','5') ){
  dotDistri(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i)
  
}


##merge c7 c8, c2c10?

table(cluster.df.add$cluster,cluster.df.add$cluster.sample.new)
   STB.BL_1 STB.BL_2 STB.BL_3 STB.BL_4 STB.BL_5 STB.BL_6 STB.CT30_1
  0       118      119      104        0       53        0          0
  1         1        7        0      288        2        0        507
  2       352       57        9        0      110        0          2
  3         0       36        1      278        6        0        665
  4         0        0        0       37        0       47         40
  5        14        3       21        0       10        0          0
  6        60       36       51        0       11        0          0
  7       112        3        1        0        7        0          0
  8        69        2        4        0        4        0          0
  9       494      174      418        0      279        0          0
  10       15        8       15        0        2        0          0
  11       25      425       24        0        7        0          1
  12        6        4        0        0        1        0          0
    
     STB.CT30_2 STB.CT30_3 STB.CT30_4 STB.CT30_5 STB.CT30_6 STB.CT30_7
  0           0          2          8         77          0          5
  1         394         39          2          0         10          0
  2           0          1          5         23          0         58
  3         447        645         12          0        159          0
  4         126         27          0          0        506          0
  5           0          0          0          6          0          3
  6           0          1          1        108          0         34
  7           0          0          0         15          0         70
  8           0          0          0         44          0          6
  9           0          3          4         27          0         34
  10          0          0          0         16          0          9
  11          0          7        173          0          0         10
  12          0          0          0          2          0          1
    
     STB.CT30_8 STB.CT30_9 villus_1 villus_10 villus_11 villus_2 villus_3
  0          44         56     1994        87         1      118      278
  1           0          1        0         0         0        0        0
  2          33         74       58       350         1     1154       41
  3           0          0        1         0         2        0        0
  4           0          1        0         0         0        0        0
  5           5          2       42       218         0       12       92
  6           4         11      141        28         2       26     1545
  7           0         28       15      1067         0      522       17
  8           0          0        8       803         0        5       79
  9          47         41        0         0         0        1        0
  10          2          5       91        96         0      299       20
  11        182          8        1         0       480        0        4
  12          5          1       56         0        22       51        5
    
     villus_4 villus_5 villus_6 villus_7 villus_8 villus_9
  0      1560      183        0        0      225        0
  1         0        0     1177       14        0     1108
  2      1031       46        1        0      129        0
  3         2        0      547       33        1       89
  4         3        0        6     1606        0        8
  5      1190      735        0        0        3        0
  6        79      103        0        0        8        0
  7       143       23        0        0        1        0
  8       101      711        0        0        0        0
  9         0        1        0        0        0        0
  10      785        8        0        0       22        0
  11        5        0        2        1        3        3
  12       42        0        0        0      983        0





####merge and rename to cluster before reclustering?


cluster.int <- cluster.df.add$cluster

table(cluster.int)
   0    1    2    3    4    5    6    7    8    9   10   11   12 
5032 3550 3535 2924 2407 2356 2249 2024 1836 1523 1393 1361 1179


cluster.int.mod <- plyr::revalue(x = cluster.int, replace = c(
     '0' = '1',
     '1' = '2',
    '2' = '3',
    '3' = '4',
    '4' = '5',
    '5' = '6',
    '6' = '7',
    '7' = '9',
    '8' = '9',
    '9' = '10',
    '10' = '3',
    '11' = '8',
    '12' = '11'

))


table(cluster.int.mod)
 1    2    3    4    5    6    7    9   10    8   11 
5032 3550 4928 2924 2407 2356 2249 3860 1523 1361 1179


cluster.df.add$cluster_int_mod <- cluster.int.mod

quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster_int_mod', title= 'in vivo vitro ',color_use = color_set,shrink.x = .8, shrink.y = .2,shuffle = FALSE,pt.size = .3)



par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('5','4','2','8','11','10','1','7','3','9','6') ){
  dotDistri(cluster = cluster.df.add[,c('cluster_int_mod','UMAP_1','UMAP_2')], id = i)
  
}

cluster.df.add.bk <- cluster.df.add

cluster.df.add$cluster <- factor(cluster.df.add$cluster_int_mod,levels = c('5','4','2','8','11','10','1','7','3','9','6') )

table(cluster.df.add$cluster)
 5    4    2    8   11   10    1    7    3    9    6 
2407 2924 3550 1361 1179 1523 5032 2249 4928 3860 2356 



saveRDS(cluster.df.add,"cluster.df.add.final.final.reclustering.rds")


quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', title= 'in vivo vitro ',color_use = color_set,shrink.x = .8, shrink.y = .2,shuffle = FALSE,pt.size = .3)




##############dot clean##############


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




par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('5','4','3','8','9','7','1','10','6','0','2') ){
  dotDistri(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i)
  
}


table(cluster.df.add$cluster)
  0    1    2    3    4    5    6    7    8    9   10 
6401 5529 4847 3797 2905 2497 2131 2037 1394 1184 1105


cluster.df.add.filterc10 <- subset(cluster.df.add, cluster != '10')

par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('5','4','3','8','9','7','1','6','0','2') ){
  dotDistri(cluster = cluster.df.add.filterc10[,c('cluster','UMAP_1','UMAP_2')], id = i)
  
}


##look for distance distribution

centers <- cluster.df.add.filterc10 %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))


q.cutoff.list <- list()

par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('5','4','3','8','9','7','1','6','0','2' ) ){
  qx <- dotDist(cluster = cluster.df.add.filterc10[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers, q = 0.98)
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
'5''4''3''8''9''7''1''6''0''2'
#'10''9''11''6''4''1''3''2''7''5''8''15''12''13''14'

#'7''9''6''11''8''4''2''10''5''1''3'

for(i in round(d.filter ,3)){cat(i,',',sep='')}
2.195,1.976,1.806,3.407,1.885,2.86,2.068,2.284,2.475,2.465
#0.877,0.668,1.047,1.806,1.596,2.214,2.219,2.002,2.521,2.126,3.31,1.45,1.162,1.393,2.084

#1.232,1.576,1.148,2.263,1.732,1.616,1.372,1.759,1.291,1.36,2.278 #q98
#1.28,1.686,1.228,2.378,2.433,1.715,1.477,1.873,1.395,1.461,3.01 #q99


# d.filter <- c( '1'=2.5, '2'=0,'3'=2,'4'=2,'5'=3,'6'=0,
#               '7'=2,'8'=2,'9'=0,'10'=0,'11'=0,'12'=0,
#               '13'=0,'14'=0,'15'=0 )


dotDist(cluster = cluster.df.add.filterc10[,c('cluster','UMAP_1','UMAP_2')], id = '7',center = centers, q = 0.8)#60%: 1.275 #80%:2.111
d.filter['7'] <- 1.39

##do cleaning
par(mfrow=c(3,3))
res.flag.dist <- list()
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('5','4','3','8','9','7','1','6','0','2') ){
  res.flag.dist[[i]] <- dotClean(cluster = cluster.df.add.filterc10[,c('cluster','UMAP_1','UMAP_2')], id = i,
           center = centers,d.filter = d.filter)
}


flag.dist.combine <- apply(do.call(cbind,res.flag.dist),1,any)
table(flag.dist.combine)
FALSE  TRUE 
31698  1024

# FALSE  TRUE 
# 32065   657

# FALSE  TRUE 
# 25309   779 

# FALSE  TRUE 
# 23702   487

# FALSE  TRUE 
# 13722   634 

cluster.df.add.filterc10.sel <- cluster.df.add.filterc10[!flag.dist.combine,]
cluster.df.add.filterc10.rm <- cluster.df.add.filterc10[flag.dist.combine,]

###quick plot cluster distribution again
#quickDimPlot(data = cluster.df.add, feature = 'cluster', title= 'late combined')
quickDimPlot_labelon(data = cluster.df.add.filterc10.sel, feature = 'cluster', title= 'in_vivo_vitro (kept)', color_use = color_good,shrink.x = .2, shrink.y = .2,shuffle = TRUE,pt.size = .3)
quickDimPlot_labelon(data = cluster.df.add.filterc10.rm, feature = 'cluster', title= 'in_vivo_vitro (remove)',color_use = color_good,shrink.x = .2, shrink.y = .2,shuffle = FALSE,pt.size = .3)



par(mfrow=c(3,3))
options(repr.plot.height=15,repr.plot.width=15)
for(i in c('5','4','3','8','9','7','1','6','0','2') ){
  dotDistri(cluster = cluster.df.add.filterc10.sel[,c('cluster','UMAP_1','UMAP_2')], id = i)
  
}



cluster.df.add.bk <- cluster.df.add
cluster.df.add <- cluster.df.add.filterc10.sel


placenta.integrated.filter  <- subset(placenta.integrated, cells = rownames(cluster.df.add.filterc10.sel) )
33478 x 31698


placenta.integrated <- placenta.integrated.filter
33478 x 31698

all.equal(rownames(cluster.df.add),colnames(placenta.integrated) ) #TRUE


placenta.integrated[['orig.ident']][,1]  <- factor(placenta.integrated[['orig.ident']][,1], levels = c(   'PLA-early1-RNA' , 'PLA-early2-RNA', 'PLA-early3-RNA', 'PLA-early4-RNA','PLA-early5-RNA','PLA-early6-RNA', 'PLA-BL_new-RNA','PLA-STBline2-RNA')   )


table(placenta.integrated[['sample']])


 PLA-BL_new-RNA   PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA 
            3930             4889             3529             4310 
  PLA-early4-RNA   PLA-early5-RNA   PLA-early6-RNA PLA-STBline2-RNA 
            3357             2705             3764             5214

sampleid <- placenta.integrated[['sample']][,1]
table(sampleid)
 PLA-BL_new-RNA   PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA 
            3930             4889             3529             4310 
  PLA-early4-RNA   PLA-early5-RNA   PLA-early6-RNA PLA-STBline2-RNA 
            3357             2705             3764             5214 


sampleid_mod <- plyr::revalue(x = sampleid, replace = c(
                                          'PLA-BL_new-RNA' = 'hTSC.BL',
                                          'PLA-STBline2-RNA' = 'hTSC.CT30'
    
                                       )
             
             )

sampleid_mod <- factor(sampleid_mod, levels = c('PLA-early1-RNA','PLA-early2-RNA','PLA-early3-RNA','PLA-early4-RNA','PLA-early5-RNA','PLA-early6-RNA', 'hTSC.BL', 'hTSC.CT30'))
table(sampleid_mod)
PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
          4889           3529           4310           3357           2705 
PLA-early6-RNA        hTSC.BL      hTSC.CT30 
          3764           3930           5214

placenta.integrated[['sample']][,1] <- sampleid_mod
table(placenta.integrated[['sample']][,1])

PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
          4889           3529           4310           3357           2705 
PLA-early6-RNA        hTSC.BL      hTSC.CT30 
          3764           3930           5214 



saveRDS(placenta.integrated,"placenta.integrated.early6-BL-CT30.final.rds")

rm(placenta.integrated.filter)

saveRDS(cluster.df.add,'cluster.df.add.integration.rds')


########dot clean ok##



placenta.integrated <- AddMetaData(placenta.integrated,metadata = Idents(placenta.integrated), col.name = 'cluster')




#by cluster
#options(repr.plot.height=7.5,repr.plot.width=12)
#options(repr.plot.width=15.5,repr.plot.height=10.5)
options(repr.plot.width=7.5,repr.plot.height=7.5)
DimPlot(placenta.integrated, reduction = "umap_rotate",label = TRUE,label.size = 6,
        cols = color_good,pt.size = .1,
       )


#by sample
#options(repr.plot.height=7.5,repr.plot.width=12)
#options(repr.plot.width=15.5,repr.plot.height=10.5)
options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(placenta.integrated, reduction = "umap_rotate",label = FALSE,label.size = 5,
        cols = c(color_good[c(1,2,3,4,5,6,8,9)],''),
        #cols = c(color_good,color_good),
        pt.size = .1,
        group.by = 'sample'
        #group.by = 'orig.ident'
        #group.by = 'cluster.ori'
       )


######plot cluster.ori distribution###

#cluster <- Idents(placenta.integrated)
cluster <- placenta.integrated[['cluster.ori']][,1]
names(table(cluster) )
'0''1''10''11''2''3''4''5''6''7''8''9''recluster_0''recluster_1''recluster_2''recluster_3''recluster_4'

'0''1''10''11''2''3''4''5''6''7''8''9''recluster_0''recluster_1''recluster_2''recluster_3''recluster_4''recluster_5''recluster_6''recluster_7'

#'1''2''3''4''5''CTB''early-MTB''early-pCTB''early-pSTB''EVT''fusion competent cell''late-pCTB''mid-pCTB''MTB''primary EVT''proliferating CTB''proliferating EVT''pSTB''STB''TE'

# levels(cluster) <- c(levels(Idents(TS.ss2)),levels(Idents(PLA.ss2)),levels(Idents(TE.ss2)))
# names(table(cluster) )
# '1''2''3''4''5''proliferating CTB''CTB''fusion competent cell''STB''proliferating EVT''primary EVT''EVT''TE''early-pCTB''mid-pCTB''late-pCTB''early-MTB''MTB''early-pSTB''pSTB'
# Idents(placenta.integrated) <- cluster


cluster <- as.character(cluster)
#names(cluster) <- names(Idents(placenta.integrated))

names(cluster) <- rownames(placenta.integrated[['cluster.ori']])

all.equal(names(Idents(placenta.integrated)),rownames(placenta.integrated[['cluster.ori']]) )
#TRUE


#cluster1 <- Idents(STB.CT30)
#cluster2 <- Idents(STB.BL)

# cluster1 <- Idents(placenta.villus)
# cluster2 <- Idents(STB.CT30)
# cluster3 <- Idents(STB.BL)
# cluster4 <- Idents(TS.CT30)
# cluster5 <- Idents(TS.BL)
# cluster6 <- Idents(TS.BT1)

names(placenta.list)
'PLA-early1-RNA''PLA-early2-RNA''PLA-early3-RNA''PLA-early4-RNA''PLA-early5-RNA''PLA-early6-RNA''STB.CT30''STB.BL'

#'PLA-early1-RNA''PLA-early2-RNA''PLA-early3-RNA''PLA-early4-RNA''PLA-early5-RNA''PLA-early6-RNA''STB.CT30''STB.BL'

cluster1 <- Idents(placenta.list[['PLA-early1-RNA']])
cluster2 <- Idents(placenta.list[['PLA-early2-RNA']])
cluster3 <- Idents(placenta.list[['PLA-early3-RNA']])
cluster4 <- Idents(placenta.list[['PLA-early4-RNA']])
cluster5 <- Idents(placenta.list[['PLA-early5-RNA']])
cluster6 <- Idents(placenta.list[['PLA-early6-RNA']])
cluster7 <- Idents(placenta.list[['STB.CT30']]) #PLA-STBline2-RNA
cluster8 <- Idents(placenta.list[['STB.BL']]) #PLA-STBline1-RNA


table(placenta.integrated[['orig.ident']][,1])
 PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA   PLA-early4-RNA 
            4889             3529             4310             3357 
  PLA-early5-RNA   PLA-early6-RNA   PLA-BL_new-RNA PLA-STBline2-RNA 
            2705             3764             3930             5214 

# PLA-BL_new-RNA   PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA 
#             3930             4889             3529             4310 
#   PLA-early4-RNA   PLA-early5-RNA   PLA-early6-RNA PLA-STBline2-RNA 
#             3357             2705             3764             5214

# PLA-BL_new-RNA   PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA 
#             4474             4959             3610             4418 
#   PLA-early4-RNA   PLA-early5-RNA   PLA-early6-RNA PLA-STBline2-RNA 
#             3436             2799             4480             5651 

unique(placenta.integrated[['orig.ident']][,1])
'PLA-early1-RNA''PLA-early2-RNA''PLA-early3-RNA''PLA-early4-RNA''PLA-early5-RNA''PLA-early6-RNA''PLA-STBline2-RNA''PLA-BL_new-RNA'

#'PLA-early1-RNA''PLA-early2-RNA''PLA-early3-RNA''PLA-early4-RNA''PLA-early5-RNA''PLA-early6-RNA''PLA-STBline2-RNA''PLA-BL_new-RNA'

#'PLA-early1-RNA''PLA-early2-RNA''PLA-early3-RNA''PLA-early4-RNA''PLA-early5-RNA''PLA-early6-RNA''PLA-STBline2-RNA''PLA-STBline1-RNA'


table(placenta.integrated[['sample']][,1])
PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
          4889           3529           4310           3357           2705 
PLA-early6-RNA        hTSC.BL      hTSC.CT30 
          3764           3930           5214



idx1 <- grep("_1$",names(cluster))
idx2 <- grep("_2$",names(cluster))
idx3 <- grep("_3$",names(cluster))
idx4 <- grep("_4$",names(cluster))
idx5 <- grep("_5$",names(cluster))
idx6 <- grep("_6$",names(cluster))
idx7 <- grep("_7$",names(cluster)) ##PLA-STBline2-RNA, STB.CT30
idx8 <- grep("_8$",names(cluster)) #PLA-BL_new-RNA, STB.BL #PLA-STBline1-RNA, STB.BL

#length(c(idx1,idx2)) == length(cluster)

#length(c(idx1,idx2,idx3)) == length(cluster) #TRUE

#length(c(idx1,idx2,idx3,idx4,idx5,idx6)) == length(cluster) #TRUE

length(c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8)) == length(cluster) #TRUE

id_idx1to6 <- gsub("_[1-6]$",replacement = "",names(cluster)[c(idx1,idx2,idx3,idx4,idx5,idx6)])

all.equal(id_idx1to6, c(names(cluster1),names(cluster2),names(cluster3),names(cluster4),names(cluster5),names(cluster6)  ) )#TRUE
cluster[c(idx1,idx2,idx3,idx4,idx5,idx6)] <- paste0('villus_',cluster[c(idx1,idx2,idx3,idx4,idx5,idx6)])


id_idx7 <- gsub("_7$",replacement = "",names(cluster)[idx7])
all.equal(id_idx7, names(cluster7))#TRUE
cluster[idx7] <- paste0('STB.CT30_',cluster[idx7])


id_idx8 <- gsub("_8$",replacement = "",names(cluster)[idx8])
all.equal(id_idx8, names(cluster8))#TRUE
cluster[idx8] <- paste0('STB.BL_',cluster[idx8])


# id_idx1 <- gsub("_1$",replacement = "",names(cluster)[idx1])
# all.equal(id_idx1, names(cluster1))#TRUE
# cluster[idx1] <- paste0('STB.CT30_',cluster[idx1])


# id_idx2 <- gsub("_2$",replacement = "",names(cluster)[idx2])
# all.equal(id_idx2, names(cluster2))#TRUE
# cluster[idx2] <- paste0('STB.BL_',cluster[idx2])


# id_idx1 <- gsub("_1$",replacement = "",names(cluster)[idx1])
# all.equal(id_idx1, names(cluster1))#TRUE
# cluster[idx1] <- paste0('placenta.villus_',cluster[idx1])


# id_idx2 <- gsub("_2$",replacement = "",names(cluster)[idx2])
# all.equal(id_idx2, names(cluster2))#TRUE
# cluster[idx2] <- paste0('STB.CT30_',cluster[idx2])



# id_idx3 <- gsub("_3$",replacement = "",names(cluster)[idx3])
# all.equal(id_idx3, names(cluster3))#TRUE
# cluster[idx3] <- paste0('STB.BL_',cluster[idx3])



# id_idx4 <- gsub("_4$",replacement = "",names(cluster)[idx4])
# all.equal(id_idx4, names(cluster4))#TRUE
# cluster[idx4] <- paste0('TS.CT30_',cluster[idx4])


# id_idx5 <- gsub("_5$",replacement = "",names(cluster)[idx5])
# all.equal(id_idx5, names(cluster5))#TRUE
# cluster[idx5] <- paste0('TS.BL_',cluster[idx5])



# id_idx6 <- gsub("_6$",replacement = "",names(cluster)[idx6])
# all.equal(id_idx6, names(cluster6))#TRUE
# cluster[idx6] <- paste0('TS.BT1_',cluster[idx6])


table(cluster) #35

            STB.BL_1             STB.BL_2             STB.BL_3 
                1541                  964                  716 
            STB.BL_4             STB.BL_5             STB.BL_6 
                 629                  577                   47 
          STB.CT30_0          STB.CT30_10           STB.CT30_2 
                1229                   43                 1044 
          STB.CT30_3           STB.CT30_6           STB.CT30_7 
                 760                  680                  251 
          STB.CT30_8 STB.CT30_recluster_0 STB.CT30_recluster_1 
                 222                  346                  326 
STB.CT30_recluster_2 STB.CT30_recluster_3 STB.CT30_recluster_4 
                 282                  240                  228 
            villus_1            villus_10            villus_11 
                2504                 2838                  525 
            villus_2             villus_3             villus_4 
                2238                 2229                 5336 
            villus_5             villus_6             villus_7 
                1857                 1740                 1659 
            villus_8             villus_9 
                1501                 1275


#  STB.BL_1             STB.BL_2             STB.BL_6 
#                 1220                  935                  338 
#             STB.BL_7   STB.BL_recluster_0   STB.BL_recluster_1 
#                  250                  136                   97 
#   STB.BL_recluster_2   STB.BL_recluster_3   STB.BL_recluster_4 
#                   90                   80                   73 
#   STB.BL_recluster_5   STB.BL_recluster_6   STB.BL_recluster_7 
#                   68                   57                   44 
#           STB.CT30_0          STB.CT30_10           STB.CT30_2 
#                 1229                   43                 1044 
#           STB.CT30_3           STB.CT30_6           STB.CT30_7 
#                  760                  680                  251 
#           STB.CT30_8 STB.CT30_recluster_0 STB.CT30_recluster_1 
#                  222                  346                  326 
# STB.CT30_recluster_2 STB.CT30_recluster_3 STB.CT30_recluster_4 
#                  282                  240                  228 
#             villus_1            villus_10            villus_11 
#                 2556                 2896                  536 
#             villus_2             villus_3             villus_4 
#                 2284                 2275                 5445 
#             villus_5             villus_6             villus_7 
#                 1895                 1776                 1693 
#             villus_8             villus_9 
#                 1532                 1301



#  STB.BL_1             STB.BL_2             STB.BL_6 
#                 1220                  935                  338 
#             STB.BL_7   STB.BL_recluster_0   STB.BL_recluster_1 
#                  250                  136                   97 
#   STB.BL_recluster_2   STB.BL_recluster_3   STB.BL_recluster_4 
#                   90                   80                   73 
#   STB.BL_recluster_5   STB.BL_recluster_6   STB.BL_recluster_7 
#                   68                   57                   44 
#           STB.CT30_0          STB.CT30_10           STB.CT30_2 
#                 1229                   43                 1044 
#           STB.CT30_3           STB.CT30_6           STB.CT30_7 
#                  760                  680                  251 
#           STB.CT30_8 STB.CT30_recluster_0 STB.CT30_recluster_1 
#                  222                  346                  326 
# STB.CT30_recluster_2 STB.CT30_recluster_3 STB.CT30_recluster_4 
#                  282                  240                  228




# placenta.villus_1   placenta.villus_10   placenta.villus_11 
#                 2556                 2896                  536 
#    placenta.villus_2    placenta.villus_3    placenta.villus_4 
#                 2284                 2275                 5445 
#    placenta.villus_5    placenta.villus_6    placenta.villus_7 
#                 1895                 1776                 1693 
#    placenta.villus_8    placenta.villus_9             STB.BL_1 
#                 1532                 1301                 1220 
#             STB.BL_2             STB.BL_6             STB.BL_7 
#                  935                  338                  250 
#   STB.BL_recluster_0   STB.BL_recluster_1   STB.BL_recluster_2 
#                  136                   97                   90 
#   STB.BL_recluster_3   STB.BL_recluster_4   STB.BL_recluster_5 
#                   80                   73                   68 
#   STB.BL_recluster_6   STB.BL_recluster_7           STB.CT30_0 
#                   57                   44                 1229 
#          STB.CT30_10           STB.CT30_2           STB.CT30_3 
#                   43                 1044                  760 
#           STB.CT30_6           STB.CT30_7           STB.CT30_8 
#                  680                  251                  222 
# STB.CT30_recluster_0 STB.CT30_recluster_1 STB.CT30_recluster_2 
#                  346                  326                  282 
# STB.CT30_recluster_3 STB.CT30_recluster_4              TS.BL_0 
#                  240                  228                  591 
#              TS.BL_1             TS.BT1_0             TS.BT1_1 
#                  243                  142                  102 
#             TS.BT1_2             TS.BT1_3            TS.CT30_0 
#                  100                    8                 1412 
#            TS.CT30_1            TS.CT30_2 
#                  130                   18


all.equal(names(cluster),colnames(placenta.integrated)) #TRUE


placenta.integrated <- AddMetaData(placenta.integrated,metadata = cluster, col.name = 'cluster.sample')

table(placenta.integrated[['cluster.sample']])

            STB.BL_1             STB.BL_2             STB.BL_3 
                1266                  874                  648 
            STB.BL_4             STB.BL_5             STB.BL_6 
                 603                  492                   47 
          STB.CT30_0          STB.CT30_10           STB.CT30_2 
                1215                   33                  967 
          STB.CT30_3           STB.CT30_6           STB.CT30_7 
                 725                  675                  134 
          STB.CT30_8 STB.CT30_recluster_0 STB.CT30_recluster_1 
                 162                  322                  318 
STB.CT30_recluster_2 STB.CT30_recluster_3 STB.CT30_recluster_4 
                 228                  230                  205 
            villus_1            villus_10            villus_11 
                2407                 2649                  508 
            villus_2             villus_3             villus_4 
                2188                 2081                 4941 
            villus_5             villus_6             villus_7 
                1810                 1733                 1654 
            villus_8             villus_9 
                1375                 1208 



all.equal(rownames(cluster.df.add),colnames(placenta.integrated)) #TRUE


placenta.integrated <- AddMetaData(placenta.integrated,metadata = cluster.df.add$cluster.sample.new, col.name = 'cluster.sample.new')

table(placenta.integrated[['cluster.sample.new']])

  STB.BL_1   STB.BL_2   STB.BL_3   STB.BL_4   STB.BL_5   STB.BL_6 STB.CT30_1 
      1266        874        648        603        492         47       1215 
STB.CT30_2 STB.CT30_3 STB.CT30_4 STB.CT30_5 STB.CT30_6 STB.CT30_7 STB.CT30_8 
       967        725        205        318        675        230        322 
STB.CT30_9   villus_1  villus_10  villus_11   villus_2   villus_3   villus_4 
       228       2407       2649        508       2188       2081       4941 
  villus_5   villus_6   villus_7   villus_8   villus_9 
      1810       1733       1654       1375       1208 




placenta.integrated <- AddMetaData(placenta.integrated,metadata = Idents(placenta.integrated), col.name = 'cluster.ori')

Idents(placenta.integrated) <- cluster.df.add$cluster

placenta.integrated <- AddMetaData(placenta.integrated,metadata = cluster.df.add$cluster, col.name = 'cluster')



##modify sample column

#placenta.integrated[['sample']] <- placenta.integrated[['orig.ident']]

table(placenta.integrated[['sample']])

PLA-early1-RNA PLA-early2-RNA PLA-early3-RNA PLA-early4-RNA PLA-early5-RNA 
          4889           3529           4310           3357           2705 
PLA-early6-RNA        hTSC.BL      hTSC.CT30 
          3764           3930           5214

#  PLA-BL_new-RNA   PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA 
#             3930             4889             3529             4310 
#   PLA-early4-RNA   PLA-early5-RNA   PLA-early6-RNA PLA-STBline2-RNA 
#             3357             2705             3764             5214 

# PLA-BL_new-RNA   PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA 
#             4474             4959             3610             4418 
#   PLA-early4-RNA   PLA-early5-RNA   PLA-early6-RNA PLA-STBline2-RNA 
#             3436             2799             4480             5651

#  PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA   PLA-early4-RNA 
#             5069             3687             4489             3515 
#   PLA-early5-RNA   PLA-early6-RNA PLA-STBline1-RNA PLA-STBline2-RNA 
#             2846             4583             3388             5651


##add sample_meta

all.equal(cluster,placenta.integrated[['cluster.sample']][,1],check.attributes = FALSE) #TRUE

placenta.integrated[['sample_meta']] <- sapply(stringr::str_split(string = cluster, pattern = "_", n = 2),function(x){x[1]} )


table(placenta.integrated[['sample_meta']])
 STB.BL STB.CT30   villus 
    3930     5214    22554 

# STB.BL STB.CT30   villus 
#     4474     5651    23702 

# STB.BL STB.CT30   villus 
#     3388     5651    24189# villus before dot clean






#Idents(placenta.integrated) <- cluster



####plot with repel and label box, for original cluster with sample id
options(repr.plot.width=22.5,repr.plot.height=10.5)


left <- -10
right <- 10
bottom <- -12
top <- 10


options(repr.plot.width=15.5,repr.plot.height=10.5)
res.p <- DimPlot(placenta.integrated, reduction = "umap_rotate",label = FALSE,label.size = 5,
        cols = c(pretty_color,color_good[1:15]),
        #cols = color_good,#unlist(map_cellcolor),
        pt.size = 1,
        repel = FALSE,
        label.box = FALSE,
        #max.overlaps = 10,
        #group.by = 'cluster.sample',
        group.by = 'cluster.sample.new',
        shuffle = FALSE
       ) +
  #scale_colour_manual(values = unlist(map_cellcolor) )  +
  #scale_colour_manual(values = color_good)  +
  #scale_colour_manual(values = unlist(map_cellcolor) )  +
  ##scale_colour_manual(values = color_snap_mod1)  +
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
  #NoLegend()+
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "total cells:",ncol(placenta.integrated),  sep=" ") ) +
#   geom_text(data = centers_shift, #the halo
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "white", 
#             size = 4.5) +
#             ##size = 6.5) +
#   geom_text(data = centers, 
#             mapping = aes(x=x,y=y,label = cellname), 
#             ##mapping = aes(x=x,y=y,label = cluster), 
#             colour = "black", 
#             size = 4.5) +
#             ##size = 6.5) +
   guides(col = guide_legend(override.aes = list(size = 4.5))) +  ##no effect ??
  xlim(left,right) + ylim(bottom,top) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")

##do labeling with repel
LabelClusters(res.p, id = "cluster.sample.new", color = 'black',fill = 'white',size = 8, repel = T,  box.padding = 1.3, box = TRUE,max.overlaps = 30, segment.color = 'black',segment.linetype = 2) #+ 
#   scale_fill_discrete(
#     name = "cluster.sample",
#     c(pretty_color,color_good[1:15]),
#     # The same color scall will apply to both of these aesthetics.
#     aesthetics = c("fill", "segment.color")
#   )


#(res.p, id = "cluster.sample", color = unique(ggplot_build(res.p)$data[[1]]$colour), size = 8, repel = T,  box.padding = 1,box = TRUE)


ggsave('pdfs/invivo-invitro-cluster-distribution.pdf',width=15.5,height=10.5)


#######
table(placenta.integrated$orig.ident)

  PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA   PLA-early4-RNA 
            4889             3529             4310             3357 
  PLA-early5-RNA   PLA-early6-RNA   PLA-BL_new-RNA PLA-STBline2-RNA 
            2705             3764             3930             4885 

 PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA   PLA-early4-RNA 
            4889             3529             4310             3357 
  PLA-early5-RNA   PLA-early6-RNA   PLA-BL_new-RNA PLA-STBline2-RNA 
            2705             3764             3930             5214

# PLA-BL_new-RNA   PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA 
#             3930             4889             3529             4310 
#   PLA-early4-RNA   PLA-early5-RNA   PLA-early6-RNA PLA-STBline2-RNA 
#             3357             2705             3764             5214


#  PLA-BL_new-RNA   PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA 
#             4474             4959             3610             4418 
#   PLA-early4-RNA   PLA-early5-RNA   PLA-early6-RNA PLA-STBline2-RNA 
#             3436             2799             4480             5651

# PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA   PLA-early4-RNA 
#             5069             3687             4489             3515 
#   PLA-early5-RNA   PLA-early6-RNA PLA-STBline1-RNA PLA-STBline2-RNA 
#             2846             4583             3388             5651


# PLA-STBline1-RNA PLA-STBline2-RNA 
#             3388             5651


# PLA-early1-RNA   PLA-early2-RNA   PLA-early3-RNA   PLA-early4-RNA 
#             5069             3687             4489             3515 
#   PLA-early5-RNA   PLA-early6-RNA PLA-STBline1-RNA PLA-STBline2-RNA 
#             2846             4583             3388             5651 
#    PLA-TS_BL-RNA   PLA-TS_BT1-RNA  PLA-TS_CT30-RNA 
#              834              352             1560






##quality check for each clusters
options(repr.plot.height=7.5,repr.plot.width=12.5)
VlnPlot(placenta.integrated,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0.1 ,group.by = 'orig.ident')

options(repr.plot.height=5.5, repr.plot.width = 30)
#options(repr.plot.height=5.5, repr.plot.width = 10)
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'sample',split.by = 'sample') + NoLegend()

options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'sample',shuffle = TRUE) #+ NoLegend()

ggsave('pdfs/PLA-in-vivo-in-vitro-source-of-sample.pdf',height=7.5,width=7.5)


# options(repr.plot.height=5.5, repr.plot.width = 10)
# DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap", group.by = 'sex',split.by = 'sex') + NoLegend()




###look at marker gene expression


DefaultAssay(placenta.integrated) <- "RNA"

placenta.integrated <- ScaleData(placenta.integrated, verbose = FALSE,features = row.names(placenta.integrated) )







####quick look at marker gene with/without splitting


DefaultAssay(placenta.integrated) <- "RNA"


options(repr.plot.height=7.5,repr.plot.width=7.5)

##CTB proliferation/hTSCs and CTB
FeaturePlot(placenta.integrated, features = 'TEAD4', reduction = "umap_rotate",pt.size = .8,slot = 'data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 

#FeaturePlot(placenta.integrated, features = 'BRCA1', reduction = "umap_rotate",pt.size = .8,slot = 'data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 


FeaturePlot(placenta.integrated, features = 'MKI67', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 

FeaturePlot(placenta.integrated, features = 'TOP2A', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 

FeaturePlot(placenta.integrated, features = 'CDH1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 

FeaturePlot(placenta.integrated, features = 'DNMT1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 


FeaturePlot(placenta.integrated, features = 'TP63', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 

FeaturePlot(placenta.integrated, features = 'ITGA6', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 

FeaturePlot(placenta.integrated, features = 'PSG8', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 



#CTB fusion

FeaturePlot(placenta.integrated, features = 'ERVFRD-1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 


##STB nascent

FeaturePlot(placenta.integrated, features = 'SH3TC2', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 


##STB Mature 1
DefaultAssay(placenta.integrated) <- "RNA"
options(repr.plot.height=7.5,repr.plot.width=7.5)
FeaturePlot(placenta.integrated, features = 'PAPPA', reduction = "umap_rotate",pt.size = .8,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 

FeaturePlot(placenta.integrated, features = 'CSHL1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 


DefaultAssay(placenta.integrated) <- "integrated"
options(repr.plot.height=5.5,repr.plot.width=35)
FeaturePlot(placenta.integrated, features = 'CSHL1', reduction = "umap_rotate",pt.size = .8,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'orig.ident') + scale_color_gradientn(colours = c('grey',color_use)) 



DefaultAssay(placenta.integrated) <- "integrated"
#DefaultAssay(placenta.integrated) <- "RNA"
options(repr.plot.height=5.5,repr.plot.width=35)
FeaturePlot(placenta.integrated, features = 'PAPPA', reduction = "umap_rotate",pt.size = .8,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample') #+ scale_color_gradientn(colours = c('grey',color_use)) 

options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'PAPPA', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ NoAxes()
ggsave('pdfs/marker_gene_PAPPA.pdf',height=5.5,width=15)

options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'LAMA3', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ NoAxes()
ggsave('pdfs/marker_gene_LAMA3.pdf',height=5.5,width=15)


options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'SH3TC2', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ NoAxes()
ggsave('pdfs/marker_gene_SH3TC2.pdf',height=5.5,width=15)


FeaturePlot(placenta.integrated, features = 'ERVFRD-1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta')
ggsave('pdfs/marker_gene_ERVFRD-1.pdf',height=5.5,width=15)

FeaturePlot(placenta.integrated, features = 'TOP2A', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta')
ggsave('pdfs/marker_gene_TOP2A.pdf',height=5.5,width=15)


FeaturePlot(placenta.integrated, features = 'CDH1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta')
ggsave('pdfs/marker_gene_CDH1.pdf',height=5.5,width=15)



FeaturePlot(placenta.integrated, features = 'ITGA6', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta')





FeaturePlot(placenta.integrated, features = 'STAT5A', reduction = "umap_rotate",pt.size = .8,slot = 'data',min.cutoff = 'q30',max.cutoff = 'q99',split.by = 'orig.ident') + scale_color_gradientn(colours = c('grey',color_use)) 



##STB Mature 2

DefaultAssay(placenta.integrated) <- "RNA"


options(repr.plot.height=7.5,repr.plot.width=7.5)
FeaturePlot(placenta.integrated, features = 'FLT1', reduction = "umap_rotate",pt.size = .8,slot = 'scale.data',min.cutoff = 'q30',max.cutoff = 'q99') + scale_color_gradientn(colours = c('grey',color_use)) 

DefaultAssay(placenta.integrated) <- "integrated"
options(repr.plot.height=5.5,repr.plot.width=35)
FeaturePlot(placenta.integrated, features = 'FLT1', reduction = "umap_rotate",pt.size = .8,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample') #+ scale_color_gradientn(colours = c('grey',color_use)) 

options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'FLT1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ NoAxes()


FeaturePlot(placenta.integrated, features = 'HDAC2', reduction = "umap_rotate",pt.size = .8,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'orig.ident') + scale_color_gradientn(colours = c('grey',color_use)) 


FeaturePlot(placenta.integrated, features = 'ENG', reduction = "umap_rotate",pt.size = .8,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'orig.ident') + scale_color_gradientn(colours = c('grey',color_use)) 


FeaturePlot(placenta.integrated, features = 'CEBPB', reduction = "umap_rotate",pt.size = .8,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'orig.ident') + scale_color_gradientn(colours = c('grey',color_use)) 






###############


#######look at original cluster.sample distribution  in common umap_rotate



options(repr.plot.height=7.5, repr.plot.width = 15)
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'sample_meta',split.by = 'sample_meta') + NoLegend()

ggsave('pdfs/in-vivo-in-vitro-sample_meta.distribution.pdf',height=7.5, width = 15)


all.equal(placenta.integrated[['cluster']][,1],Idents(placenta.integrated), check.attributes = FALSE ) #TRUE
#cluster is the new integration cluster

placenta.integrated[['cluster']] <- Idents(placenta.integrated)


#look at new cluster of integration with splitting
options(repr.plot.height=5.5, repr.plot.width = 15)
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'cluster',split.by = 'sample_meta') + NoLegend()

#look at cluster.sample with splitting
DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'cluster.sample',split.by = 'sample_meta') #+ NoLegend()

DimPlot(object = placenta.integrated, label = FALSE,cols=c(color_good,color_good),pt.size = 0.1,label.size = 8,reduction = "umap_rotate", group.by = 'cluster.sample.new',split.by = 'sample_meta') #+ NoLegend()


ggsave('pdfs/in-vivo-in-vitro-cluster_sample_new.distribution.pdf',height=7.5, width = 15)






####stat cluster class################

##stat sample_meta and cluster.sample
mat_sample_with_cluster_lib <- table(placenta.integrated[['sample_meta']][,1], placenta.integrated[['cluster.sample']][,1])
mat_sample_with_cluster_lib.df <- as.data.frame.matrix(mat_sample_with_cluster_lib)

##stat sample_meta and cluster (new of integration)

mat_sample_with_cluster <- table(placenta.integrated[['sample_meta']][,1], placenta.integrated[['cluster']][,1])
mat_sample_with_cluster.df <- as.data.frame.matrix(mat_sample_with_cluster)

	      1	     2	   3	 4	     5	  6	     7	  8	     9  10

STB.BL	  597	464	  174	312	    298	  94	136	 1367	479	9
STB.CT30  348	250	  97	1076	1910  727	147	 260	390	9
villus	  5328	4704  4479	2333	638	  1626	1805  2	    497	1142



mat_cluster_lib_with_cluster <- table(placenta.integrated[['cluster.sample']][,1], placenta.integrated[['cluster']][,1])
mat_cluster_lib_with_cluster.df <- as.data.frame.matrix(mat_cluster_lib_with_cluster)

options(repr.plot.width=15, repr.plot.height=15)
res.p <- pheatmap(mat_cluster_lib_with_cluster.df,fontsize_col = 16,fontsize_row = 12,cluster_rows = FALSE, kmeans_k = NA,
         cluster_cols = FALSE,main = "summary of cluster_lib with cluster",scale = "none",
         show_rownames = TRUE,show_colnames = TRUE,
         clustering_method = "complete", #complete,ward.D2
         clustering_distance_cols  = "euclidean",#euclidean,correlation 
         clustering_distance_rows  = "euclidean",#euclidean,correlation
         
         #border_color = "black",
         color = colorRampPalette(c("navy","white","firebrick"))(20),
         #color = colorRampPalette(c("grey","white","firebrick"))(20),
         #color = c(colorRampPalette(colors = c("white","red" ))(20)),
         ##breaks = seq(0,1,by=1/20),
         display_numbers = TRUE,
         number_color = 'red',
         fontsize_number = 20,
         number_format = "%i"
)






####get marker gene plot: integrated and split (use assay of integrated ?)


##integrated umap
DefaultAssay(placenta.integrated) <- "integrated"
#DefaultAssay(placenta.integrated) <- "RNA"

res.p.list <- list()
marker.gene <- c('TEAD4','DNMT1','ERVFRD-1','SH3TC2','CSHL1','PAPPA','EGLN3','FLT1')

for (i in marker.gene){
    options(repr.plot.height=5.5,repr.plot.width=5.5)

    res.p.list[[i]] <- FeaturePlot(placenta.integrated, features = i, reduction = "umap_rotate",pt.size = .3,slot = 'scale.data',min.cutoff = 'q10',max.cutoff = 'q99') + NoAxes() + scale_color_gradientn(colours = c(color_use)) + xlim(-6,6) + ylim(-8,6) 

    print(res.p.list[[i]] )
    ggsave(paste0('pdfs/marker_gene_integration/marker_gene.integration.in_one_umap.',i,'.pdf'),height=5.5,width=5.5)
    
}

res.wrap <- wrap_plots(res.p.list[marker.gene],nrow=2,ncol=4, heights = 8.5,byrow=TRUE) + plot_annotation(title = 'marker gene of integration')


options(repr.plot.height=8.5,repr.plot.width=18.5)
print(res.wrap)


ggsave('pdfs/marker_gene_integration/marker_gene.integration.in_one_umap.grid2_by_4.pdf',height=8.5,width=18.5)



##split
DefaultAssay(placenta.integrated) <- "integrated"
#DefaultAssay(placenta.integrated) <- "RNA"


options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'PAPPA', reduction = "umap_rotate",pt.size = .5,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ NoAxes()
ggsave('pdfs/marker_gene_integration/marker_gene_PAPPA.split_3.pdf',height=5.5,width=15)


FeaturePlot(placenta.integrated, features = 'FLT1', reduction = "umap_rotate",pt.size = .5,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ NoAxes()
ggsave('pdfs/marker_gene_integration/marker_gene_FLT1.split_3.pdf',height=5.5,width=15)


FeaturePlot(placenta.integrated, features = 'SH3TC2', reduction = "umap_rotate",pt.size = .5,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ NoAxes()
ggsave('pdfs/marker_gene_integration/marker_gene_SH3TC2.split_3.pdf',height=5.5,width=15)


options(repr.plot.height=5.5,repr.plot.width=15)
FeaturePlot(placenta.integrated, features = 'LAMA3', reduction = "umap_rotate",pt.size = .5,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta') #+ NoAxes()
ggsave('pdfs/marker_gene_integration/marker_gene_LAMA3.split_3.pdf',height=5.5,width=15)


#
FeaturePlot(placenta.integrated, features = 'TOP2A', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta')
ggsave('pdfs/marker_gene_TOP2A.pdf',height=5.5,width=15)


FeaturePlot(placenta.integrated, features = 'CDH1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta')
ggsave('pdfs/marker_gene_CDH1.pdf',height=5.5,width=15)


FeaturePlot(placenta.integrated, features = 'ITGA6', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta')


FeaturePlot(placenta.integrated, features = 'ERVFRD-1', reduction = "umap_rotate",pt.size = .1,slot = 'scale.data',min.cutoff = 'q0',max.cutoff = 'q99',split.by = 'sample_meta')
ggsave('pdfs/marker_gene_ERVFRD-1.pdf',height=5.5,width=15)








#####get cluster.df.add again with cluster.sample ####

cluster <- Idents(placenta.integrated) #placenta.integrated[['cluster_lib']]$cluster_lib#
#names(cluster) <- rownames(placenta.integrated[['cluster_lib']])

umap <- placenta.integrated@reductions$umap_rotate@cell.embeddings
#umap <- placenta@reductions$umap_rotate@cell.embeddings

colnames(umap) <- c('UMAP_1','UMAP_2')

all.equal(names(cluster),rownames(umap)) #TRUE


cluster.df <- data.frame(cluster=cluster,umap)

metadata <- placenta.integrated@meta.data

all.equal(rownames(cluster.df),rownames(metadata))#TRUE

intersect(colnames(cluster.df), colnames(metadata) ) #cluster

grep('^cluster$',colnames(metadata) )
#49

#49,50
#13,16

all.equal(cluster.df$cluster,metadata$cluster)#no


cluster_mod <- plyr::revalue(x = cluster.df$cluster, replace = c(
                                          '0' = '1',
                                          '1' = '2',
                                          '2' = '3',
                                          '3' = '4',
                                          '4' = '5',
                                          '5' = '6',
                                          '6' = '7',
                                          '7' = '8',
                                          '8' = '9',
                                          '9' = '10'
    
                                       )
             
             )


all.equal(cluster_mod, metadata$cluster,check.attributes = FALSE) #TRUE

cluster.df$cluster <- cluster_mod
all.equal(cluster.df$cluster , metadata$cluster,check.attributes = FALSE) #TRUE


metadata$cluster <- NULL

#metadata[,c(49,50)] <- NULL
#metadata[,c(13,16)] <- NULL

#all.equal(metadata,placenta@meta.data[,-c(49,50)])
#all.equal(metadata,placenta@meta.data[,-c(13,16)]) #TRUE


cluster.df.add <- cbind(cluster.df, metadata)

table(cluster.df.add$cluster)

   1    2    3    4    5    6    7    8    9   10 
6273 5418 4750 3721 2846 2447 2088 1629 1366 1160 

#   0    1    2    3    4    5    6    7    8    9 
# 6273 5418 4750 3721 2846 2447 2088 1629 1366 1160 

# 0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 4180 3948 3804 3748 3730 2449 2440 2075 2060 1645 1104 1049  576  420 

all.equal(cluster.df.add$cluster , placenta.integrated[['cluster']][,1]) #TRUE


#placenta.integrated[['cluster']][,1] <- cluster_mod

###remove unwanted column

idx <- grep("RNA_snn_res",colnames(cluster.df.add) )
cluster.df.add[,colnames(cluster.df.add)[idx]] <- NULL

idx <- grep("integrated_snn",colnames(cluster.df.add) )
cluster.df.add[,colnames(cluster.df.add)[idx]] <- NULL


cluster.df.add$cluster.sample <- factor(cluster.df.add$cluster.sample, names(table(cluster.df.add$cluster.sample) ) )
cluster.df.add$sample_meta <- factor(cluster.df.add$sample_meta, names(table(cluster.df.add$sample_meta) ) )




cluster.df.add$cluster.sample <- factor(cluster.df.add$cluster.sample)


all.equal(Idents(placenta.integrated),cluster.df.add$cluster) #no
all.equal(placenta.integrated[['cluster']][,1], cluster.df.add$cluster) #TRUE

all.equal(rownames(cluster.df.add),colnames(placenta.integrated)) #TRUE


table(cluster.df.add$cluster.sample )
            STB.BL_1             STB.BL_2             STB.BL_3 
                1266                  874                  648 
            STB.BL_4             STB.BL_5             STB.BL_6 
                 603                  492                   47 
          STB.CT30_0          STB.CT30_10           STB.CT30_2 
                1215                   33                  967 
          STB.CT30_3           STB.CT30_6           STB.CT30_7 
                 725                  675                  134 
          STB.CT30_8 STB.CT30_recluster_0 STB.CT30_recluster_1 
                 162                  322                  318 
STB.CT30_recluster_2 STB.CT30_recluster_3 STB.CT30_recluster_4 
                 228                  230                  205 
            villus_1            villus_10            villus_11 
                2407                 2649                  508 
            villus_2             villus_3             villus_4 
                2188                 2081                 4941 
            villus_5             villus_6             villus_7 
                1810                 1733                 1654 
            villus_8             villus_9 
                1375                 1208





saveRDS(cluster.df.add,'cluster.df.add.final.rds') #31698 x 36
saveRDS(placenta.integrated,'placenta.integrated.early6-BL-CT30.final.rds') #33478 x 31698 



#############filter placenta.integrated and cluster.df.add by new cluster.df.add_ct30 and rename cluster_sample

##rename for BL cluster (noname if filtered by qc)
table(cluster.df.add$sample_meta)
   STB.BL STB.CT30   villus 
    3930     5214    22554


table( Idents(placenta.villus.early) )
   1    2    3    4    5    6    7    8    9   10   11 
2504 2238 2229 5336 1857 1740 1659 1501 1275 2838  525

29132 x 23702 


cellid_villus <- colnames(subset(placenta.integrated, sample_meta == "villus" ) )

cellid_villus_mod <- gsub(pattern = "_[123456]$",replacement = "",cellid_villus  )

table(cellid_villus_mod %in% colnames(placenta.villus.early))
TRUE 
22554


cellid_villus_keep <- colnames(subset(placenta.integrated, sample_meta == "villus" ) )
#22554




cluster.df.add_ct30 #5135 x 16

table(cluster.df.add_ct30$cluster)
 1    2    3    4    5    6    7    8    9 
1229 1044  760  228  326  680  240  346  282 

rownames(cluster.df.add_ct30) <- paste0( rownames(cluster.df.add_ct30), "_7" )

cluster.df.add_ct30$cluster <- paste0("STB.CT30_", cluster.df.add_ct30$cluster  )

table(cluster.df.add_ct30$cluster)
STB.CT30_1 STB.CT30_2 STB.CT30_3 STB.CT30_4 STB.CT30_5 STB.CT30_6 STB.CT30_7 
      1229       1044        760        228        326        680        240 
STB.CT30_8 STB.CT30_9 
       346        282



table(rownames(cluster.df.add_ct30) %in% colnames(placenta.integrated) )
FALSE  TRUE 
  250  4885

table(colnames(subset(placenta.integrated, sample_meta == "STB.CT30" )) %in% rownames(cluster.df.add_ct30)   )
FALSE  TRUE 
  329  4885


cellid_ct30_keep <- intersect( rownames(cluster.df.add_ct30),colnames(placenta.integrated)    )
#4885




cluster.df.add_grpts #4474 × 23

table(cluster.df.add_grpts$cluster)
  1    2    3    4    5    6 
1541  964  716  629  577   47

rownames(cluster.df.add_grpts) <- paste0( rownames(cluster.df.add_grpts), "_8" )

cluster.df.add_grpts$cluster <- paste0("STB.BL_", cluster.df.add_grpts$cluster  )


table(rownames(cluster.df.add_grpts) %in% colnames(placenta.integrated) )
FALSE  TRUE 
  544  3930

table(colnames(subset(placenta.integrated, sample_meta == "STB.BL" )) %in% rownames(cluster.df.add_grpts)   )
TRUE 
3930 

cellid_grpts_keep <- intersect( rownames(cluster.df.add_grpts),colnames(placenta.integrated)    )
#3930

placenta.integrated.filter <- subset(placenta.integrated, cells = c(cellid_villus_keep,cellid_ct30_keep,cellid_grpts_keep)  )

33478 x 31369 


table(rownames(cluster.df.add) %in% c(cellid_villus_keep,cellid_ct30_keep,cellid_grpts_keep) )
FALSE  TRUE 
  329 31369

cluster.df.add.filter <- cluster.df.add[rownames(cluster.df.add) %in% c(cellid_villus_keep,cellid_ct30_keep,cellid_grpts_keep),]
31369 × 36

all.equal(rownames(cluster.df.add.filter), colnames(placenta.integrated.filter)  )
#TRUE


cluster.df.add_ct30.filter <- cluster.df.add_ct30[cellid_ct30_keep,]
cluster.df.add_grpts.filter <- cluster.df.add_grpts[cellid_grpts_keep,]


table(rownames(cluster.df.add_ct30.filter) %in% rownames(cluster.df.add.filter) )
TRUE 
4885


table( rownames(cluster.df.add_grpts.filter) %in% rownames(cluster.df.add.filter) )
TRUE 
3930

length(cellid_villus_keep) + nrow(cluster.df.add_ct30.filter) + nrow(cluster.df.add_grpts.filter) == nrow(cluster.df.add.filter) #TRUE


all.equal(sort(rownames(cluster.df.add.filter)),sort(c(cellid_villus_keep,rownames(cluster.df.add_ct30.filter),rownames(cluster.df.add_grpts.filter)))  ) #TRUE


all.equal(rownames(cluster.df.add.filter),c(cellid_villus_keep,rownames(cluster.df.add_ct30.filter),rownames(cluster.df.add_grpts.filter)) ) #TRUE




###substitue new grpts cluster with old one in cluster.df.add.filte with check
cluster.df.add.filter$cluster.sample.new  <- as.character(cluster.df.add.filter$cluster.sample)


table(cluster.df.add.filter[rownames(cluster.df.add_ct30.filter),'cluster.sample.new'] == paste0("STB.CT30_", cluster.df.add_ct30.filter$cluster_ori))
TRUE 
4885

cluster.df.add.filter[rownames(cluster.df.add_ct30.filter),'cluster.sample.new'] <- cluster.df.add_ct30.filter$cluster




table(cluster.df.add.filter[rownames(cluster.df.add_ct30.filter),'cluster.sample.new'] == paste0("STB.CT30_", cluster.df.add_ct30.filter$cluster_ori))




all.equal(as.character(cluster.df.add.filter[rownames(cluster.df.add_grpts.filter),'cluster.sample.new']) ,cluster.df.add_grpts.filter$cluster , check.attributes = FALSE) #TRUE




all.equal( rownames(cluster.df.add.filter), colnames(placenta.integrated.filter)   ) #TRUE



table(cluster.df.add.filter$cluster) #cca integration cluster

   1    2    3    4    5    6    7    8    9   10 
6264 5375 4729 3677 2776 2428 2079 1524 1357 1160 

table(cluster.df.add.filter$cluster.sample) #old ct30 cluster
             STB.BL_1             STB.BL_2             STB.BL_3 
                1266                  874                  648 
            STB.BL_4             STB.BL_5             STB.BL_6 
                 603                  492                   47 
          STB.CT30_0          STB.CT30_10           STB.CT30_2 
                1215                    0                  967 
          STB.CT30_3           STB.CT30_6           STB.CT30_7 
                 725                  675                    0 
          STB.CT30_8 STB.CT30_recluster_0 STB.CT30_recluster_1 
                   0                  322                  318 
STB.CT30_recluster_2 STB.CT30_recluster_3 STB.CT30_recluster_4 
                 228                  230                  205 
            villus_1            villus_10            villus_11 
                2407                 2649                  508 
            villus_2             villus_3             villus_4 
                2188                 2081                 4941 
            villus_5             villus_6             villus_7 
                1810                 1733                 1654 
            villus_8             villus_9 
                1375                 1208 


table(cluster.df.add.filter$cluster.sample.new) #new ct30 cluster
STB.BL_1   STB.BL_2   STB.BL_3   STB.BL_4   STB.BL_5   STB.BL_6 STB.CT30_1 
      1266        874        648        603        492         47       1215 
STB.CT30_2 STB.CT30_3 STB.CT30_4 STB.CT30_5 STB.CT30_6 STB.CT30_7 STB.CT30_8 
       967        725        205        318        675        230        322 
STB.CT30_9   villus_1  villus_10  villus_11   villus_2   villus_3   villus_4 
       228       2407       2649        508       2188       2081       4941 
  villus_5   villus_6   villus_7   villus_8   villus_9 
      1810       1733       1654       1375       1208


table(Idents(placenta.integrated.filter))
 0    1    2    3    4    5    6    7    8    9 
6264 5375 4729 3677 2776 2428 2079 1524 1357 1160


saveRDS(cluster.df.add.filter,'cluster.df.add.final.final.rds')
saveRDS(placenta.integrated.filter, 'placenta.integrated.early6-BL-CT30.final.final.rds')



cluster.df.add <- cluster.df.add.filter
31369 × 37

placenta.integrated <- placenta.integrated.filter
33478 x 31369 



###river plot with cluster.df.add.filter$cluster and cluster.df.add.filter$cluster.sample.new


# edges <- read.csv(text="N1,N2,Value
# From EU,To China,170.4
# From China,To EU,350.6
# From EU,To US,426.0
# From US,To EU,272.7
# From China,To US,482
# From US,To China,116",
# header=T, stringsAsFactors=FALSE)
# print(edges)


# nodes <- data.frame(ID=unique(c(edges$N1, edges$N2)),x=rep(c(1,2), each=3), y=c(1,2,3,2,1,3))

# cols <- c(China="#00990066",EU ="#00009966",US ="#99000066")
# style <- sapply(nodes$ID, function(id){ list(col=cols[ gsub("(From|To) ", "", id) ], simplify=FALSE) } )
                
# r <- riverplot::makeRiver(nodes=nodes, edges=edges, styles=style)
# par(bg="grey98")
# d <- list(srt=0, textcex=1.5) # default style
# plot(r, plot_area=1, nodewidth=10, default_style=d)


source("makeRiverPlot_liger.r")
library(riverplot)


cluster.df.add.villus <- subset(cluster.df.add, sample_meta == 'villus')
cluster.villus <- cluster.df.add.villus$cluster.sample.new
names(cluster.villus) <- rownames(cluster.df.add.villus)
cluster.villus <- factor(cluster.villus)

cluster.df.add.grpts <- subset(cluster.df.add, sample_meta == 'STB.BL')
cluster.grpts <- cluster.df.add.grpts$cluster.sample.new
names(cluster.grpts) <- rownames(cluster.df.add.grpts)
cluster.grpts <- factor(cluster.grpts)

cluster.df.add.ct30 <- subset(cluster.df.add, sample_meta == 'STB.CT30')
cluster.ct30 <- cluster.df.add.ct30$cluster.sample.new
names(cluster.ct30) <- rownames(cluster.df.add.ct30)
cluster.ct30 <- factor(cluster.ct30)

cluster.all  <- cluster.df.add$cluster
names(cluster.all) <- rownames(cluster.df.add)


pdf('pdfs/villus_vs_STB_BL-integration-riverplot.pdf',width=9.5,height=8.5)
options(repr.plot.width=9.5,repr.plot.height=8.5)
makeRiverPlot(cluster1 = cluster.villus, cluster2 = cluster.grpts, cluster_consensus = cluster.all,
              min.frac = 0.1,#0.2,
              river.yscale = 3.5, #node.order = set_node_order,
              river.usr = c(0, 1, -0.6, 1.6),
              river.lty = 0,
              label.cex = 1.2
             )

#ggsave('pdfs/villus_vs_STB_BL-integration-riverplot.pdf',width=9.5,height=8.5)
dev.off()

pdf('pdfs/villus_vs_STB_CT30-integration-riverplot.pdf',width=9.5,height=8.5)
options(repr.plot.width=9.5,repr.plot.height=8.5)
makeRiverPlot(cluster1 = cluster.villus, cluster2 = cluster.ct30, cluster_consensus = cluster.all,
             min.frac = 0.2,river.yscale = 5, #node.order = set_node_order,
              river.usr = c(0, 1, -0.6, 1.6),
              river.lty = 0,
              label.cex = 1.0
             )
#ggsave('pdfs/villus_vs_STB_CT30-integration-riverplot.pdf',width=9.5,height=8.5)

dev.off()

options(repr.plot.width=9.5,repr.plot.height=7.5)
makeRiverPlot(cluster1 = cluster.ct30, cluster2 = cluster.grpts, cluster_consensus = cluster.all,
              min.frac = 0.2,river.yscale = 5, #node.order = set_node_order,
              river.usr = c(0, 1, -0.6, 1.6),
              river.lty = 0,
              label.cex = 1.0
             )

ggsave('pdfs/CTB_CT30_vs_STB_BL-integration-riverplot.pdf',width=9.5,height=7.5)


###confusion map
cM <- as.matrix(ArchR::confusionMatrix(cluster.df.add$cluster, cluster.df.add$cluster.sample.new))



options(repr.plot.width=15.5,repr.plot.height=15)
res.p <- pheatmap::pheatmap(
    as.data.frame(log2(cM+1)),#[

      # c('4','5','6','7','10','9','8','2','1','3'), c('STB.BL_1','STB.BL_2','STB.BL_3','STB.BL_4','STB.BL_5','STB.BL_6','STB.CT30_0','STB.CT30_10','STB.CT30_2','STB.CT30_3','STB.CT30_6','STB.CT30_7','STB.CT30_8','STB.CT30_recluster_0','STB.CT30_recluster_1','STB.CT30_recluster_2','STB.CT30_recluster_3','STB.CT30_recluster_4','villus_1','villus_10','villus_11','villus_2','villus_3','villus_4','villus_5','villus_6','villus_7','villus_8','villus_9')
                      
                    # ],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = color_set_yellowbrick.flat[['Oranges.9']],#color_tfdev,
    border_color = 'white',fontsize = 24,
    display_numbers = TRUE,
    fontsize_number = 12
)

pdf('pdfs/villus_STB_BL-STB_CT30-integration-confusion_map.pdf',width=15.5,height=15)
print(res.p)

dev.off()

#ggsave('pdfs/villus_STB_BL-STB_CT30-integration-confusion_map.pdf',width=15.5,height=15)



#####confusion map or anuven plot?


#cM <- as.matrix(ArchR::confusionMatrix(cluster.df.add$cluster, cluster.df.add$cluster.sample))
cM <- as.matrix(ArchR::confusionMatrix(cluster.df.add$cluster, cluster.df.add$cluster.sample.new))

preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

pheatmap::pheatmap(as.data.frame(cM))


options(repr.plot.width=15.5,repr.plot.height=15)
#pheatmap::pheatmap(as.data.frame(log(cM+1))
res.p <- pheatmap::pheatmap(
    as.data.frame(log2(cM+1))[
                      #c('8','2','1','4','5','3','6','7'),
                      #c('9','7','10','2','4','1','3','5','6','8')
                      #c('8','7','6','5','4','3','2','1'),
                      #c('1','2','3','4','5','6','7','8','9','10')
                      #c('9','5','6','4','1','3','7','2','8'),
                      #c('10','9','11','6','4','1','2','3','7','5','8')
       #c('4','5','6','7','10','9','8','2','1','3'), c('STB.BL_1','STB.BL_2','STB.BL_3','STB.BL_4','STB.BL_5','STB.BL_6','STB.CT30_0','STB.CT30_10','STB.CT30_2','STB.CT30_3','STB.CT30_6','STB.CT30_7','STB.CT30_8','STB.CT30_recluster_0','STB.CT30_recluster_1','STB.CT30_recluster_2','STB.CT30_recluster_3','STB.CT30_recluster_4','villus_1','villus_10','villus_11','villus_2','villus_3','villus_4','villus_5','villus_6','villus_7','villus_8','villus_9')
                      
                     ],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = color_set_yellowbrick.flat[['Oranges.9']],#color_tfdev,
    border_color = 'white',fontsize = 24,
    display_numbers = TRUE,
    fontsize_number = 12
)

pdf('pdfs/rna-invivo-invitro-match.pdf',width=15.5,height=15)
print(res.p)
dev.off()






##display cluster distribution for each sample_meta
par(mfrow=c(3,3))
options(repr.plot.height=25,repr.plot.width=25)

#for(i in grep('^STB.CT30',levels(cluster.df.add$cluster.sample),value=TRUE)   ){
#for(i in grep('^villus',levels(cluster.df.add$cluster.sample),value=TRUE)   ){
for(i in grep('^STB.BL',levels(cluster.df.add$cluster.sample),value=TRUE)   ){
 options(repr.plot.height=7.5,repr.plot.width=7.5)
  dotDistri(cluster = cluster.df.add[,c('cluster.sample','UMAP_1','UMAP_2')], id = i)
}


##the integration cluster distribution
par(mfrow=c(3,3))
for(i in levels(cluster.df.add$cluster)   ){
 options(repr.plot.height=7.5,repr.plot.width=7.5)
  dotDistri(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i)
}







#display cluster from both dataset
options(repr.plot.width=15,repr.plot.height=15)
par(mfrow=c(2,2))
dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_7',id2='STB.CT30_6',title='TS (proliferative CTB)') 


dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_2',id2='STB.CT30_2',title='CTB-1') 


dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_1',id2='STB.CT30_0',title='CTB-2') 


dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_6',id2='STB.CT30_3',title='CTB Fusion') 



options(repr.plot.width=15,repr.plot.height=15)
par(mfrow=c(2,2))
dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_recluster_1',id2='STB.CT30_recluster_0',title='Nascent STB') 


dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_recluster_6',id2='STB.CT30_recluster_0',title='Nascent STB') 


dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_recluster_0',id2='STB.CT30_recluster_2',title='STB Mature1 (PAPPA+)') 


dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_recluster_0',id2='STB.CT30_recluster_3',title='STB Mature1 (PAPPA+)') 


options(repr.plot.width=15,repr.plot.height=15)
par(mfrow=c(2,2))
dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_recluster_0',id2='STB.CT30_recluster_1',title='STB Mature2 (FLT1+)') 


dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_recluster_5',id2='STB.CT30_8',title='STB apoptosis(?)') 


dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_recluster_2',id2='STB.CT30_8',title='STB apoptosis(?)') 

dotDistri_both(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2','orig.ident')], 
               id1 = 'STB.BL_recluster_3',id2='STB.CT30_8',title='STB apoptosis(?)') 



    
#display cluster from three dataset ?
    






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





######
dotDistri_both = function (cluster = NULL, id1 = NULL, id2 = NULL,title = NULL){
    #id1 for atac, id2 for rna
    cluster.sel1 = cluster[(cluster$orig.ident == 'PLA-STBline1-RNA' & cluster$cluster == id1),] #STB.BL
    n_sel1 = nrow(cluster.sel1)
    cluster.sel2 = cluster[(cluster$orig.ident == 'PLA-STBline2-RNA' & cluster$cluster == id2),] #STB.CT30
    n_sel2 = nrow(cluster.sel2)
    
    ##stopifnot(sum(is.na(cluster.sel1)) == 0)
    ##stopifnot(sum(is.na(cluster.sel2)) == 0)
    
    ##write output for kdeplot
    ##fileout <- paste( "atac_c",id1,'_rna_c',id2,'.txt',sep=''   )
    ##write.table(x = rbind(cluster.sel1,cluster.sel2),file = fileout,sep='\t',row.names = TRUE,col.names = TRUE,quote = FALSE)
     
    #color_atac = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel1)),'pink' ,'red'   )
    #color_rna = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel2)),'lightblue' ,'navy'   ) 
   
    #stopifnot(length(color_atac) == nrow(cluster.sel1))
    #stopifnot(length(color_rna) == nrow(cluster.sel2))
    
    plot(cluster$UMAP_1,cluster$UMAP_2,pch = 16, type='p',col='grey',cex=0.5,xlab='UMAP1',ylab='UMAP2',main=paste("STB.BL cells cluster ",id1,", n = ",n_sel1,"\nSTB.CT30 cells cluster ",id2,", n = ",n_sel2,sep=''),cex.main=2,cex.axis=2,cex.lab=2 ) 
    #points(cluster.sel2$UMAP1,cluster.sel2$UMAP2,pch = 16, cex=0.2,col='navy') #rna
    #points(cluster.sel1$UMAP1,cluster.sel1$UMAP2,pch = 16, cex=0.2,col='red') #atac

    points(cluster.sel1$UMAP_1,cluster.sel1$UMAP_2,pch = 16, cex=0.5,col='red') #atac
    points(cluster.sel2$UMAP_1,cluster.sel2$UMAP_2,pch = 16, cex=0.5,col='blue' ) #rna
    
    text(x = 9,y=5.5,labels = title,cex=2,pos = 2,adj = 0.5)
    legend('bottomleft',legend=c('STB.BL','STB.CT30'),fill = c('STB.BL'='red','STB.CT30'='blue'),box.lwd = 0,bg=NULL,cex=1.2,pt.cex=1)
    #legend(x = -6.5,y=6,legend=c('STB.BL','STB.CT30'),fill = c('STB.BL'='red','STB.CT30'='blue'),box.lwd = 0,bg=NULL)
    
    return(paste("cluster ",id1,", ",id2," ok. ",sep='') )
    #return(paste("cluster ",id1,", ",id2," ok. coords output to ",fileout,sep='') )
}


######?
dotDistri_three = function (cluster = NULL, id1 = NULL, id2 = NULL,title = NULL){
    #id1 for atac, id2 for rna
    cluster.sel1 = cluster[(cluster$orig.ident == 'PLA-STBline1-RNA' & cluster$cluster == id1),] #STB.BL
    n_sel1 = nrow(cluster.sel1)
    cluster.sel2 = cluster[(cluster$orig.ident == 'PLA-STBline2-RNA' & cluster$cluster == id2),] #STB.CT30
    n_sel2 = nrow(cluster.sel2)
    
    ##stopifnot(sum(is.na(cluster.sel1)) == 0)
    ##stopifnot(sum(is.na(cluster.sel2)) == 0)
    
    ##write output for kdeplot
    ##fileout <- paste( "atac_c",id1,'_rna_c',id2,'.txt',sep=''   )
    ##write.table(x = rbind(cluster.sel1,cluster.sel2),file = fileout,sep='\t',row.names = TRUE,col.names = TRUE,quote = FALSE)
     
    #color_atac = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel1)),'pink' ,'red'   )
    #color_rna = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel2)),'lightblue' ,'navy'   ) 
   
    #stopifnot(length(color_atac) == nrow(cluster.sel1))
    #stopifnot(length(color_rna) == nrow(cluster.sel2))
    
    plot(cluster$UMAP_1,cluster$UMAP_2,pch = 16, type='p',col='grey',cex=0.5,xlab='UMAP1',ylab='UMAP2',main=paste("STB.BL cells cluster ",id1,", n = ",n_sel1,"\nSTB.CT30 cells cluster ",id2,", n = ",n_sel2,sep=''),cex.main=2,cex.axis=2,cex.lab=2 ) 
    #points(cluster.sel2$UMAP1,cluster.sel2$UMAP2,pch = 16, cex=0.2,col='navy') #rna
    #points(cluster.sel1$UMAP1,cluster.sel1$UMAP2,pch = 16, cex=0.2,col='red') #atac

    points(cluster.sel1$UMAP_1,cluster.sel1$UMAP_2,pch = 16, cex=0.5,col='red') #atac
    points(cluster.sel2$UMAP_1,cluster.sel2$UMAP_2,pch = 16, cex=0.5,col='blue' ) #rna
    
    text(x = 9,y=5.5,labels = title,cex=2,pos = 2,adj = 0.5)
    legend('bottomleft',legend=c('STB.BL','STB.CT30'),fill = c('STB.BL'='red','STB.CT30'='blue'),box.lwd = 0,bg=NULL,cex=1.2,pt.cex=1)
    #legend(x = -6.5,y=6,legend=c('STB.BL','STB.CT30'),fill = c('STB.BL'='red','STB.CT30'='blue'),box.lwd = 0,bg=NULL)
    
    return(paste("cluster ",id1,", ",id2," ok. ",sep='') )
    #return(paste("cluster ",id1,", ",id2," ok. coords output to ",fileout,sep='') )
}



##:::::::::::::on map differential expressed gene test (see work_postintegrate_seuratv4.r in r413 conda):::::::::::###### 


# ######identify conserved marker gene in the integrated space
# #https://satijalab.org/seurat/articles/integration_introduction.html
# # For performing differential expression after integration, we switch back to the original
# # data
# DefaultAssay(placenta.integrated) <- "RNA" #full gene set, instead of 2000 anchor.feature gene

# table(Idents(placenta.integrated))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 4180 3948 3804 3748 3730 2449 2440 2075 2060 1645 1104 1049  576  420



# # ##merge CT30 and BL?
# # sample_meta <- placenta.integrated[['sample_meta']]

# # table(sample_meta)
# #  STB.BL STB.CT30   villus 
# #   3388     5651    24189

# ##or filter and keep only placenta.villus and CT30 ?
# placenta.integrated.villus_CT30 <- subset(placenta.integrated, subset =  sample_meta %in% c('STB.CT30', 'villus')  )

# table(placenta.integrated.villus_CT30[['sample_meta']])

# STB.CT30   villus 
#     5651    24189



# share.markers.villus_vs_CT30 <- FindConservedMarkers(placenta.integrated, ident.1 = '6', grouping.var = "sample_meta", verbose = TRUE)
# head(share.markers)



# ######compare different genes within sample_meta
# ######identify and plot TF heatmap pairwise heatmap with CT30 vs placental.villus, along cluster.sample of placenta.villus?











######save and reload

#saveRDS(object = placenta.integrated,file = 'placenta.integrated.rds')
#placenta.integrated <- readRDS('placenta.integrated.rds')
##saveRDS(object = placenta.integrated,file = 'placenta.integrated.BL-CT30.rds')


    
##saveRDS(object = placenta.integrated,file = 'placenta.integrated.early6-BL-CT30.rds')
#placenta.integrated <- readRDS(file = 'placenta.integrated.early6-BL-CT30.rds')

#saveRDS(object = placenta.integrated,file = 'placenta.integrated.early6-BL-CT30.final.rds')


saveRDS(cluster.df.add,'cluster.df.add.final.final.reclustering.rds')
saveRDS(placenta.integrated, 'placenta.integrated.early6-BL-CT30.final.final.reclustering.rds')




placenta.integrated <- readRDS('placenta.integrated.early6-BL-CT30.final.final.rds')

cluster.df.add <- readRDS('cluster.df.add.final.final.rds')



####integration done########


















##get tne final cluster.df.add
umap <- Embeddings(placenta.integrated,reduction = 'umap')
cluster <- Idents(placenta.integrated)
names(table(cluster) )

cluster.df = data.frame(cluster=cluster,UMAP_1=umap[,1],UMAP_2=umap[,2]) 
rownames(cluster.df) = rownames(placenta.integrated@meta.data)

stopifnot(all.equal (rownames(cluster.df),rownames(placenta.integrated@meta.data) ) )
cluster.df.add <- cbind(cluster.df,placenta.integrated@meta.data)


stopifnot(all.equal(Idents(placenta.integrated),cluster.df.add$cluster,check.attributes = FALSE) )


##annotate source TS.ss2,  PLA.ss2 , TE.ss2
levels(Idents(TS.ss2)) #'1''2''3''4''5'  ##5
levels(Idents(PLA.ss2)) #'proliferating CTB' 'CTB' 'fusion competent cell' 'STB' 'proliferating EVT' 'primary EVT ''EVT' ##7
levels(Idents(TE.ss2)) #'TE' 'early-pCTB' 'mid-pCTB' 'late-pCTB' 'early-MTB' 'MTB' 'early-pSTB' 'pSTB' #8


levels(cluster.df.add$cluster)
levels(cluster.df.add$cluster) <- c(levels(Idents(TS.ss2)),levels(Idents(PLA.ss2)),levels(Idents(TE.ss2)))


reds <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf')
blues <- c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','darkblue')


oranges <- c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6')
purples <- c('#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')

greys <- c('#ffffff',
 '#f0f0f0',
 '#d9d9d9',
 '#bdbdbd',
 '#969696',
 '#737373',
 '#525252',
 '#252525')


greens <- c('#f7fcf5',
 '#e5f5e0',
 '#c7e9c0',
 '#a1d99b',
 '#74c476',
 '#41ab5d',
 '#238b45',
 '#006d2c',
 '#00441b')

map_cellcolor <- list(
    
    '1'=reds[5],
    '2'=reds[4],
    '3'=reds[3],
    '4'=reds[2],
    '5'=reds[1],
    
    'proliferating CTB' =blues[2],
    'CTB' =blues[2],
    "fusion competent cell"=blues[5],
    'STB'=blues[6], 
    'proliferating EVT' =blues[3],
    'primary EVT'=blues[3],
    'EVT'=blues[4],
    
    'TE'=greens[2], 
    'early-pCTB' =greens[3],
    'mid-pCTB'=greens[4], 
    'late-pCTB'=greens[5], 
    'early-MTB'=greens[6], 
    'MTB'=greens[7], 
    'early-pSTB'=greens[8], 
    'pSTB'=greens[9]
 )   
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
    
# #     '1'=purples[4], #STB4
# #     '2'=purples[5],  #STB5
# #     '3'=purples[2],  #STB2
# #     '4'=purples[1], #STB1
# #     '5'=purples[3], #STB3
# #     '6'='#8B0000', #Syncytial knot
# #     '7'='#d8daeb', #naive STB
# #     '8'='#7f3b08', #STB-new
# #     '9'='darkgreen', #CTB
# #     '10'=''
#  )
# map_cellname <- list(
    
#     '1'='TS_c1',
#     '2'='TS_c2',
#     '3'='TS_c3',
#     '4'='TS_c4',
#     '5'='TS_c5',
#     'proliferating CTB' ='proliferating_CTB',
#     'CTB' =blues[2],
#     'fusion competent cell'=blues[3],
#     'STB'=blues[4], 
#     'proliferating EVT' =blues[5],
#     'primary EVT '=blues[6],
#     'EVT'='navy',
    
#     'TE'='black', 
#     'early-pCTB' =purples[1],
#     'mid-pCTB'=purples[2], 
#     'late-pCTB'=purples[3], 
#     'early-MTB'=oranges[3], 
#     'MTB'=oranges[5], 
#     'early-pSTB'='brown', 
#     'pSTB'='darkbrown'
#  ) 

#cellname <- as.character(Idents(placenta.integrated))


# cellcolor <- as.character(cluster.df.add$cluster)


# #for(i in 1:length(cellname)){ cellname[i] = map_cellname[[ cellname[i] ]] }
# for(i in 1:length(cellcolor)){   
#     ##cat('cellcolor i:',i," ",cellcolor[i],'\n',sep='')
#     ##if(map_cellcolor[[ as.character(cellcolor[i]) ]] == ""){cat('cellcolor i:',cellcolor[i],'\n',sep='')}
#     cellcolor[i] = map_cellcolor[[ as.character(cellcolor[i]) ]] 

# }

# #cluster.df.add[,'cellname'] <- factor(cellname)
# cluster.df.add[,'cellcolor'] <- factor(cellcolor)


sample = "TS, TE, placenta(ywLiu) integration"
##label on cluster
options(repr.plot.height=5,repr.plot.width=6.5,repr.plot.res = 150)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster )) +
#ggplot(cluster.df.add,aes(x=tSNE1,y=tSNE2,col=APcluster_Merge  )) +
  geom_point(size = 0.6,show.legend = TRUE,alpha= 1 ) +

  scale_colour_manual(values = unlist(map_cellcolor) )  +
  #scale_colour_manual(values = color_good)  +
  #scale_colour_manual(values = unlist(map_cellcolor) )  +
  ##scale_colour_manual(values = color_snap_mod1)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  theme(
        #legend.position = 'right',
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
  #ggtitle(paste(sample, "total cells:",nrow(cluster.df.add),  sep=" ") ) +
#   geom_text(data = centers_shift, #the halo
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "white", 
#             size = 4.5) +
#             ##size = 6.5) +
#   geom_text(data = centers, 
#             mapping = aes(x=x,y=y,label = cellname), 
#             ##mapping = aes(x=x,y=y,label = cluster), 
#             colour = "black", 
#             size = 4.5) +
#             ##size = 6.5) +
   guides(col = guide_legend(override.aes = list(size = 4.5))) +  ##no effect ??
#   xlim(left,right) + ylim(bottom,top) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")


####









##::::::::::::::::projection::::::::::::::::::::

# ############
# p1 <- DimPlot(placenta.integrated, reduction = "umap", group.by = "tech")
# p2 <- DimPlot(placenta.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
#     repel = TRUE) + NoLegend()

# options(repr.plot.height = 7.5, repr.plot.width = 15)
# p1 + p2



# ######################################################

# placenta.list <- SplitObject(panc8, split.by = "tech") #five techs
# placenta.list <- placenta.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")] #remove indrop


# ##preprocess: log transformation, find variable genes
# for (i in 1:length(placenta.list)) {
#     placenta.list[[i]] <- NormalizeData(placenta.list[[i]], verbose = FALSE)
#     placenta.list[[i]] <- FindVariableFeatures(placenta.list[[i]], selection.method = "vst", 
#         nfeatures = 2000, verbose = FALSE)
# }



# #############part 1, integration for a reference ########
# reference.list <- placenta.list[c("celseq", "celseq2", "smartseq2")]
# placenta.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30) #multiple object list
# placenta.integrated <- IntegrateData(anchorset = placenta.anchors, dims = 1:30)


# #     celseq    celseq2 fluidigmc1     indrop  smartseq2 
# #      1004       2285        638       8569       2394 


# ###analysis the integrated object
# library(ggplot2)
# library(cowplot)
# library(patchwork)
# # switch to integrated assay. The variable features of this assay are automatically
# # set during IntegrateData
# DefaultAssay(placenta.integrated) <- "integrated"

# # Run the standard workflow for visualization and clustering
# placenta.integrated <- ScaleData(placenta.integrated, verbose = FALSE)
# placenta.integrated <- RunPCA(placenta.integrated, npcs = 30, verbose = FALSE)
# placenta.integrated <- RunUMAP(placenta.integrated, reduction = "pca", dims = 1:30)
# p1 <- DimPlot(placenta.integrated, reduction = "umap", group.by = "tech")
# p2 <- DimPlot(placenta.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
#     repel = TRUE) + NoLegend()

# options(repr.plot.height = 7.5, repr.plot.width = 15)
# p1 + p2



#########part 2, projection (label transfer) reference data onto a query object (use this, but only id transfered, no umap) ####

#query = ss2
#reference = tenx
cluster.df <- data.frame( metadata(tenx)$umap, metadata(tenx)$Clusters )
colnames(cluster.df) <- c('UMAP1','UMAP2','clusters')
cluster.df$clusters <- gsub(pattern = '^Cluster',replacement = '',x = cluster.df$clusters)
cluster.df$clusters <- factor(cluster.df$clusters,
                              levels=sort(as.numeric(unique(cluster.df$clusters)))
                             )

#start to transfer
placenta.query <- placenta.list[["ss2"]]
placenta.reference <- placenta.list[["tenx"]] #not a integrated object also ok
placenta.anchors <- FindTransferAnchors(reference = placenta.reference, query = placenta.query, 
    dims = 1:30) #find anchors between an already integrated reference and a query object



predictions <- TransferData(anchorset = placenta.anchors, 
                            #refdata = placenta.reference@meta.data$res.1.5,
                            refdata = cluster.df$clusters,
                            dims = 1:30) #630, only transfer id no coordinates

placenta.query <- AddMetaData(placenta.query, metadata = predictions)

#quick look for query (ss2) marker genes
options(repr.plot.height = 5, repr.plot.width = 15)
VlnPlot(placenta.query,features = c('CGA','KRT7','HLA-G','MKI67'),#'ERVFRD-1'),
        group.by='predicted.id',
        cols = color_good[c(1,11,4,9)],
        ncol=4
       )


saveRDS(predictions,'predictions.ss2_transferto_tenx.rds')


all.equal( rownames(placenta.reference@meta.data), rownames(cluster.df) )#TRUE
placenta.reference <- AddMetaData(placenta.reference,metadata = cluster.df)
VlnPlot(placenta.reference,features = c('CGA','ERVFRD-1','KRT7','HLA-G'),#'MKI67'),
        group.by='clusters',
        cols = color_good,
        ncol=4
       )

#quick look for reference (10x scRNA lsi) marker genes



#placenta.query <- AddMetaData(placenta.query, metadata = predictions)
ref.id <- table( cluster.df$clusters)
# 1    2    3    4    5    6    7    8    9   10   11 
# 222  443 1072  504  761  674  988  813  371  496   93 
proj.id <- table(predictions$predicted.id)

##plot reference tenx 
# 1  11   4   9 
#269  31   6 324 




color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" )


##tenx lsi umap and clusters
options(repr.plot.height = 7.5, repr.plot.width = 10)
plot(cluster.df[,1:2],
     col=color_good[as.numeric(cluster.df$cluster)],
     pch=19,cex=0.5)
legend("topright",
       fill = color_good[as.numeric(levels(cluster.df$cluster))],
       legend=levels(cluster.df$cluster),
       ncol=2,
       x.intersp = 0.25 #the space between legend columns
      )



##draw barplot ##barplot for projection id classification

tmp <- rep(0,11)
names(tmp) <- 1:11
tmp[names(proj.id)] <- proj.id
proj.id <- tmp

options(repr.plot.height = 5, repr.plot.width = 6)
barplot(proj.id,col=color_good,border = color_good,main = 'Seuratv3 CCA projection \nSmartSeq2 TS data to 10x scRNA')



############
idents <- ss2@active.ident
predict.df <- data.frame(annotation=idents,
                         ss2@reductions$umap@cell.embeddings,
                         predictions = predictions[names(idents),'predicted.id']
                        )

table(predict.df$annotation)
#  1   2   3   4   5 
# 91 140 175 130  94 
table(predict.df$predictions)
#  1  11   4   9 
#269  31   6 324 


###plot ss2 seurat result umap and clusters to verify cluster.df is correct
DimPlot(object = ss2, reduction = "umap" ,label = TRUE,label.size = 12)

##ss2 original id
options(repr.plot.height = 7.5, repr.plot.width = 7.5)
plot(predict.df[,2:3],
     col=color_good[as.numeric(predict.df$annotation)],
     pch=19,cex=0.5)
legend("topright",
       fill = color_good[as.numeric(levels(predict.df$annotation))],
       legend=levels(predict.df$annotation),
       ncol=2,
       x.intersp = 0.25 #the space between legend columns
      )

##the predicted id (use this)
plot(predict.df[,2:3],
     col=color_good[as.numeric(as.character(predict.df$predictions))],
     pch=19,cex=1)
legend("topright",
       fill = color_good[sort(as.numeric(levels(predict.df$predictions)))],
       legend=sort(as.numeric(levels(predict.df$predictions)))
       #ncol=2,
       #x.intersp = 0.25 #the space between legend columns
      )




##########the prediction heatmap to visualize prediction power

df <- predictions[,-c(1,ncol(predictions))]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

df.z <- t(apply( df,1,cal_z_score ))
#df.z.1 <- scale(df,scale=TRUE,center=FALSE) #scale is different with cal_z_score!!


anno_row = data.frame(row.names = rownames(df.z) )
anno_row$predicted = predictions$predicted.id

anno_color = list(
       'predicted' = c('1'=color_good[1],'11'=color_good[11],'4'=color_good[4],'9' = color_good[9])

)

library('viridisLite')
options(repr.plot.width = 15,repr.plot.height = 15)
pheatmap::pheatmap(df.z,show_rownames = FALSE,fontsize_col = 25,
                   annotation_colors = anno_color,
                   annotation_row = anno_row,
                   #clustering_distance_rows="euclidean",clustering_method="complete",
                   #viridis(n = 20,option = "C"),kmeans_k = NA,
                   border_color = "NA"#,breaks=seq(-1,3,by=4/20)
                  )











saveRDS(predict.df,'predict.df.ss2_transferto_tenx.rds')



####compare cluster ids by riverplot ###

groups <- predict.df[,c(1,4)]
colnames(groups) <- c('ss2_clusters','tenx_clusters')
#groups$cluster.snap <- droplevels(groups$cluster.snap)
#groups$cluster.signac <- droplevels(groups$cluster.signac)

data <- reshape2::melt(table(groups)) #melt from table result
data.gather <- ggforce::gather_set_data(data,1:2)

library('cowplot')
library(ggforce)

options(repr.plot.height = 25, repr.plot.width = 15)
ggplot(data.gather,aes(x,id = id,split=y,value=value) ) +
  geom_parallel_sets(aes(fill = ss2_clusters),alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(color='red',size = 15,family=2,angle = 0) +
  #scale_color_manual(values = color_good_new) + #no use
  #scale_fill_manual(values = color_good_new) +
  #scale_color_brewer(type='qual',palette='Set3') +
  #scale_fill_brewer(type='qual',palette='Set3') +
  theme_cowplot() +
  theme(
         axis.title = element_text(size =20)

       )


# placenta.query <- placenta.list[["fluidigmc1"]]
# placenta.anchors <- FindTransferAnchors(reference = placenta.integrated, query = placenta.query, 
#     dims = 1:30) #find anchors between an already integrated reference and a query object

# predictions <- TransferData(anchorset = placenta.anchors, refdata = placenta.integrated$celltype, 
#     dims = 1:30) #use TransferData instead of IntegrateData to project
# placenta.query <- AddMetaData(placenta.query, metadata = predictions)


#########part 3, projection (label transfer) tenx reference data onto a query object techmann ####


#start to transfer
placenta.query <- placenta.list[["techmann"]]
placenta.reference <- placenta.list[["tenx"]] #not a integrated object also ok
placenta.anchors <- FindTransferAnchors(reference = placenta.reference, query = placenta.query, 
    dims = 1:30) #find anchors between an already integrated reference and a query object



predictions <- TransferData(anchorset = placenta.anchors, 
                            #refdata = placenta.reference@meta.data$res.1.5,
                            refdata = cluster.df$clusters, #tenx clusters
                            dims = 1:30) #630, only transfer id no coordinates


saveRDS(predictions,'predictions.techmann_transferto_tenx.rds')


proj.id <- table(predictions$predicted.id)

#   1   11    2    8    9  ##no function matured CTB?
# 227   61 1805  178 1126 


##check the prediction score table by heatmap with kmeans or hclust



df <- predictions[,-c(1,ncol(predictions))]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

df.z <- t(apply( df,1,cal_z_score ))
#df.z.1 <- scale(df,scale=TRUE,center=FALSE) #scale is different with cal_z_score!!


anno_row = data.frame(row.names = rownames(df.z) )
anno_row$predicted = predictions$predicted.id

anno_color = list(
       'predicted' = c('1'=color_good[1],'11'=color_good[11],'2'=color_good[2],'8'=color_good[8],
                      '9' = color_good[9])

)

library('viridisLite')
options(repr.plot.width = 15,repr.plot.height = 15)
pheatmap::pheatmap(df.z,show_rownames = FALSE,fontsize_col = 25,
                   annotation_colors = anno_color,
                   annotation_row = anno_row,
                   #clustering_distance_rows="euclidean",clustering_method="complete",
                   #viridis(n = 20,option = "C"),kmeans_k = NA,
                   border_color = "NA"#,breaks=seq(-1,3,by=4/20)
                  )


