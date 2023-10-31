
##plot integration cluster as in liger "work_liger.allgene.te.r"


library('ggplot2')

#library(harmony)

library('magrittr')

#library('RANN')

library('patchwork')

library(hexbin)
library(RColorBrewer)

#library("ggpointdensityplot")
library(viridis)


library(grid)
library(gridExtra)





# ###look up color sets in ArchR###
# library('ArchR')
# color_list <- ArchR::ArchRPalettes

# options(repr.plot.width=12,repr.plot.height=6)
# par(mfrow=c(3,3),oma=c(1,1,1,1))
# for(name in names(color_list)){ 
#   len = length(color_list[[name]])
#   color = color_list[[name]]
#   #barplot(rep(5,len),col = color,main = name,cex.main=2)
#   barplot(1:len,col = color,main = name,cex.main=2)
# }

# #saveRDS(color_list,'ArchR.color_list.rds')

# ##


# ##python yellowbrick color palettes##
# library(rlist)
# color_set_yellowbrick <- readRDS('color_set_yellowbrick.rds')
# color_set_yellowbrick.flat <- list.flatten(color_set_yellowbrick)
# options(repr.plot.width=12,repr.plot.height=6)
# par(mfrow=c(3,3),oma=c(1,1,1,1))
# for(name in names(color_set_yellowbrick.flat)){ 
#   len = length(color_set_yellowbrick.flat[[name]])
#   color = color_set_yellowbrick.flat[[name]]
#   #barplot(rep(5,len),col = color,main = name,cex.main=2)
#   barplot(1:len,col = color,main = name,cex.main=2)
# }

# saveRDS(color_set_yellowbrick.flat,'color_set_yellowbrick.flat.rds')

# ############Buen colors########
# library(BuenColors)
# color_set0 <- jdb_color_maps #17 different colors
# names(color_set0) <- NULL
# #plot(1:17,1:17,pch = 19, cex = 5,col=jdb_color_maps)

# #discrete colors
# color_set1 <- jdb_palette("solar_extra") #9 discrete but gradient colors
# color_set2 <- jdb_palette("brewer_spectra") #9 discrete but gradient colors
# color_set3 <- jdb_palette("flame_light") #9 discrete but gradient colors, good!

# color_set3_ext12 <- colorRampPalette(colors = as.character(color_set3))(12)
# color_set3_ext17 <- colorRampPalette(colors = as.character(color_set3))(17)

# #############ArchR colors############
# #hmcols <- colorRamps::blue2green2red(length(bks) ) #colors
# color_peak <- ArchR::paletteContinuous(set = 'solarExtra',n=256,reverse=FALSE)  
# color_tfdev = ArchR::paletteContinuous(set = 'blueYellow',n=257,reverse=FALSE)                       
# #color_ga <- paletteContinuous(set='solarExtra',n=257,reverse=FALSE) 
# #color_ga <- paletteContinuous(set='horizon',n=257,reverse=FALSE)                      
# #color_ga <- paletteContinuous(set='horizonExtra',n=257,reverse=FALSE)  #good        
# color_rna <- ArchR::paletteContinuous(set='greenBlue',n=256,reverse=FALSE)
# #color_ga <- paletteContinuous(set='blueYellow',n=257,reverse=FALSE)
# color_ga <- ArchR::paletteContinuous(set='greyMagma',n=257,reverse=FALSE)

# color_rna <- colorRampPalette(c('grey','red'))(10) 

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


color_good <- c("#E7D654", "#6F1482" ,"#DC7035", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", "#63AC4E", "#D181B0" ,
                "#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,"#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,
                "#CA362E" ,"#2B3918","#1E1E1E" )


color_set <- c('0'   =    '#2679b4',
'1'   =    '#fd7f28',
'2'   =    '#2f9d6a',
'3'    =   '#d42a2f',
'4'    =   '#a94af8',
'5'    =   '#8b564c',
'6'    =   '#e17ac1',
'7'   =    '#b5bc66',
'8'   =    '#29bece',
'9'    =   '#b1cae9',
'10'  =    '#fdba7d',
'11'   =   '#9ade8d',
'12'   =   '#fd9998',
'13'  =    'grey',
'14'  =    'navy',
'15'  =  'brown',
'16'   =   'black',
'17'   =   'yellow'
              )



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


########
map_cellcolor_liger <- list(
  '7' = 'black', 
 '3'= '#b2182b',
 '8'= 'darkgreen',#'#FCC140',#'#FAC669',#'#EABFE1',#'#FBC04D',#'#d6604d',
 '1'= 'darkred',
 '4'= '#fedcbd',#'#FFA300',#'#fddbc7',
 '2'= '#FB8D3C',#''#AB855A',#'#7A5BA1',#'#f7f7f7',
 '5'= '#d6604d',#'#6F1482',#'navy',#purples[4],#'darkgreen',#'#d6604d', #'#d1e5f0',
 '10'= 'grey',
 '6'= '#67001f',#'#4393c3',
 '11'= '#053061',#'#2166ac',
 '9'= '#92c5de'
)


options(repr.plot.height = 7.5, repr.plot.width = 7.5)
barplot(1:length(map_cellcolor_rna),col = unlist(map_cellcolor_rna),main ='rna_map_cellcolor' ,cex.main=2,names.arg = names(map_cellcolor_rna))



######quick look and save cluster.df.add_glue with color
source('quickDimPlot_labelon.r')

cluster.df.add <- cluster.df.add_glue

colnames(cluster.df.add)[3:4] <- c("UMAP_1",'UMAP_2')

table(cluster.df.add_glue$cluster)


options(repr.plot.width = 8.5, repr.plot.height=7.5)
quickDimPlot_labelon(data = cluster.df.add, feature = 'cluster', color_use = map_cellcolor_liger, title= 'late GLUE integration', shrink.x = 0.2, shrink.y = 0.05,shuffle = FALSE,pt.size = .1,width = 8.5,height = 7.5 )

ggsave(filename = "pdfs/PLA-term-RNA-ATAC-glue.UMAP.labelon.hotcolor.pdf",height=7.5,width=8.5,useDingbats=FALSE)






###select one global color set###
color <- color_good



sample <- 'placenta full term RNA integration with ATAC'

######



##readin GLUE cluster data


cluster.df.add_glue <- read.table('cluster_df_add.glue.txt',sep='\t', header = TRUE, row.names = 1)
#use res = 0.7

##add 1 to louvain and leiden cluster


#cluster.df.add_glue$louvain <- cluster.df.add_glue$louvain + 1
cluster.df.add_glue$leiden <- cluster.df.add_glue$leiden + 1

#table(cluster.df.add_glue$louvain)
# 1    2    3    4    5    6    7    8    9   10 
#3080 2791 2685 2623 1843 1677 1461  956  310  239

table(cluster.df.add_glue$leiden)
# 1     2     3     4     5     6     7     8     9    10    11 
# 11371  8784  8091  6959  4444  4436  2374  1588   376   234    16 

 1    2    3    4    5    6    7    8    9   10   11   12   13 
7223 6957 6821 6516 6100 4438 3308 3095 2038 1549  377  235   16

#  1    2    3    4    5    6    7    8    9   10 
# 3233 2894 2526 2182 1874 1800 1466  745  636  309 


table(cluster.df.add_glue$cluster_lib,cluster.df.add_glue$domain)
#   atac  rna
#   c1  3457 5731
#   c10    0  206
#   c11    0  270
#   c2  3583 4248
#   c3  3186 4116
#   c4  5910 3390
#   c5  2983 2594
#   c6  2651 1822
#   c7  1449 1112
#   c8  1050  393
#   c9   423   99

 atac  rna
  c1  3457 5731
  c10    0  206
  c11    0  270
  c2  3583 4248
  c3  3186 4116
  c4  5910 3390
  c5  2983 2594
  c6  2651 1822
  c7  1449 1112
  c8  1050  393
  c9   423   99

#      atac  rna
#   c1  1801 1858
#   c10    0 1360
#   c2  1749  599
#   c3  1539 1716
#   c4  1044 1638
#   c5  1048  994
#   c6   613  492
#   c7   346  405
#   c8   194  170
#   c9     0   99


cluster.df.add_rna <- readRDS('cluster.df.add.rna.rds')
cluster.df.add_atac <- readRDS('cluster.df.add.atac.rds')

table(cluster.df.add_rna$cluster)
  1    2    3    4    5    6    7    8    9   10   11 
5731 4248 4116 3390 2594 1822 1112  393   99  206  270

#  1    2    3    4    5    6    7    8    9   10 
# 1858  599 1716 1638  994  492  405  170   99 1360

table(cluster.df.add_atac$cluster)
  1    2    3    4    5    6    7    8    9 
3457 3583 3186 5910 2983 2651 1449 1050  423

# 1    2    3    4    5    6    7    8 
# 1801 1749 1539 1044 1048  613  346  194


all.equal(rownames(cluster.df.add_glue) ,c(rownames(cluster.df.add_rna) , rownames(cluster.df.add_atac))  )#TRUE, glue = rna + atac


cluster.df.add_glue$cluster <- cluster.df.add_glue$leiden #use this?
#cluster.df.add_glue$cluster <- cluster.df.add_glue$louvain
cluster.df.add_glue$cluster <- factor(cluster.df.add_glue$cluster)


###plot glue integration clusters on umap######
centers <- cluster.df.add_glue %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP1), 
        y = median(x = UMAP2))

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

###

cluster.df.add_glue.shuffle <- cluster.df.add_glue[sample(rownames(cluster.df.add_glue)),]

##label on 
options(repr.plot.height=5,repr.plot.width=5.5)
#options(repr.plot.height=5,repr.plot.width=5)
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
ggplot(cluster.df.add_glue.shuffle,aes(x=UMAP1,y=UMAP2,col=cluster  )) +
  geom_point(size = 0.2,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
  ##xlim(left,right) + ylim(bottom,top) +
  #scale_colour_manual(values = color_snap_mod1)  +
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
  ggtitle(paste(sample, "\ntotal cells:",nrow(cluster.df.add_glue),  sep=" ") ) +
  geom_text(data = centers_shift, #the halo
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "white", 
            size = 6.5) +
  geom_text(data = centers, 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            size = 6.5) +
  guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/PLA-term-RNA-ATAC-glue.UMAP.labelon.leidencluster.pdf",height=5,width=5,useDingbats=FALSE)

ggsave(filename = "pdfs/PLA-term-RNA-ATAC-glue.UMAP.labelon.leidencluster.withlegend.pdf",height=5,width=5.5,useDingbats=FALSE)

##ggsave(filename = "pdfs/PLA-term-RNA-ATAC-glue.UMAP.labelon.louvaincluster.pdf",height=5,width=5,useDingbats=FALSE)
#ggsave(filename = "pdfs/PLA-8w-RNA-ATAC-liger.UMAP.labelon.final.pdf",height=5,width=5,useDingbats=FALSE)


##by data type
table(cluster.df.add_glue$batch)

atac_D1 atac_D2 atac_D3 atac_D5 atac_D6 atac_D9  rna_D1  rna_D2  rna_D3  rna_D5 
   4999    3550    3897    3248    4625    4373    2563    3566    3720    3886 
 rna_D6  rna_D9 
   4545    5701

# atac_D1 atac_D2  rna_D1  rna_D2 
#    4913    3421    3647    5684

table(cluster.df.add_glue$domain)

 atac   rna 
24692 23981


options(repr.plot.height=5,repr.plot.width=5.5)
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
ggplot(cluster.df.add_glue,aes(x=UMAP1,y=UMAP2,col= domain  )) +
#ggplot(cluster.df.add_glue,aes(x=UMAP1,y=UMAP2,col= batch  )) +
  geom_point(size = 0.1,show.legend = TRUE,alpha= .3 ) +
  #scale_colour_manual(values = c('atac_D1' = 'pink', 'atac_D2'='red','rna_D1'='lightblue', 'rna_D2'='navy')  )  +
  scale_colour_manual(values = c('rna' = 'red', 'atac'='navy')  )  +
  ##xlim(left,right) + ylim(bottom,top) +
  #scale_colour_manual(values = color_snap_mod1)  +
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
  ggtitle(paste(sample, "\ntotal cells:",nrow(cluster.df.add_glue),  sep=" ") ) +
  # geom_text(data = centers_shift, #the halo
  #           mapping = aes(x=x,y=y,label = cluster), 
  #           colour = "white", 
  #           size = 6.5) +
  # geom_text(data = centers, 
  #           mapping = aes(x=x,y=y,label = cluster), 
  #           colour = "black", 
  #           size = 6.5) +
  guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/PLA-term-RNA-ATAC-glue.UMAP.by_data.pdf",height=5,width=5.5,useDingbats=FALSE)


###############



# ##########density 2d with contour plot UMAP#########
# ggplot(cluster.df.add_glue, aes(UMAP1, UMAP2)) + 
#   geom_point(color = "lightgray") +
#   geom_density_2d(color='red') +
#   #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
#   scale_fill_gradientn(colors = c("#FFEDA0", "#FEB24C", "#F03B20")) +
#   #scale_fill_distiller(palette = "Spectral", direction = -1)+
#   theme(
#         legend.position = 'none',
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
#        )
# #########



#########pixel style hexbin density scatterplot (grid)#######
##https://cran.r-project.org/web/packages/hexbin/vignettes/hexagon_binning.pdf
pdf( "pdfs/pixel_style.cluster.density.pdf",height=7.5,width=7.5,useDingbats = FALSE)

options(repr.plot.height=7.5,repr.plot.width=7.5)

bin<-hexbin(cluster.df.add_glue$UMAP1, cluster.df.add_glue$UMAP2, xbins=40,xbnds = range(cluster.df.add_glue$UMAP1)*1.3,ybnds = range(cluster.df.add_glue$UMAP2)*1.3 )

my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))

hb <- plot(bin,main="Nuclei pileup density" , colramp=my_colors , legend=F,xlab='UMAP1',ylab='UMAP2')

dev.off()


#####



########cluster distribution plot ####
dotDistri = function (cluster = NULL, id = NULL, type = NULL){
    if(type == 'all'){
        cluster.sel = cluster[(cluster$cluster == id),]
    }else{
      cluster.sel = cluster[(cluster$domain == type & cluster$cluster == id),]
    }
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    color = ''
    #color = ifelse(type == 'atac','red','navy')
    
    if(type == 'atac'){
      color = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel)),'black' ,'red'   )
    }else if (type == 'rna'){
        color = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel)),'lightblue' ,'navy'   ) 
    }else{
        color = ifelse(test=grepl(pattern='^placenta_atac',x=rownames(cluster.sel)),'navy' ,'red'   )
        
    }
    
    plot(cluster$UMAP1,cluster$UMAP2,pch = 16, type='p',col='grey',cex=0.3,xlab='UMAP1',ylab='UMAP2',main=paste('integration'," cells cluster ",id,"\nn = ",n_sel,sep=''),cex.main=2,cex.axis=2,cex.lab=2 ) 
    points(cluster.sel$UMAP1,cluster.sel$UMAP2,pch = 16, cex=0.3,col=color)
    return(paste("cluster ",id," ok",sep='') )
}


#leiden cluster (res 0.6)

# options(repr.plot.height=15,repr.plot.width=15)
# par(mfrow=c(3,3))
# dotDistri(cluster = cluster.df.add_glue, id = '9', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '11', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '8', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '1', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '2', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '5', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '4', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '3', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '6', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '7', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '10', type = 'all')

pdf("pdfs/PLA-term-RNA-ATAC-glue.integration_cluster_distribution_part1.pdf",height=15,width=15,useDingbats = FALSE)

options(repr.plot.height=15,repr.plot.width=15)
par(mfrow=c(3,3))
dotDistri(cluster = cluster.df.add_glue, id = '9', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '11', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '8', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '4', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '2', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '5', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '1', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '3', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '6', type = 'all')

dev.off()



pdf("pdfs/PLA-term-RNA-ATAC-glue.integration_cluster_distribution_part2.pdf",height=15,width=15,useDingbats = FALSE)

options(repr.plot.height=15,repr.plot.width=15)
par(mfrow=c(3,3))

dotDistri(cluster = cluster.df.add_glue, id = '7', type = 'all')
dotDistri(cluster = cluster.df.add_glue, id = '10', type = 'all')

dev.off()


# options(repr.plot.height=15,repr.plot.width=15)
# par(mfrow=c(3,3))
# dotDistri(cluster = cluster.df.add_glue, id = '11', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '10', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '4', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '8', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '5', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '2', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '7', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '1', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '3', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '6', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '9', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '12', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '13', type = 'all')


# #louvain clusters
# options(repr.plot.height=15,repr.plot.width=15)
# par(mfrow=c(3,3))
# dotDistri(cluster = cluster.df.add_glue, id = '4', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '5', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '2', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '3', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '9', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '2', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '6', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '1', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '8', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '7', type = 'all')
# dotDistri(cluster = cluster.df.add_glue, id = '10', type = 'all')



######










###############use confusion map method to predict integration cluster major cluster matching name in atac and rna ######

cluster.df.add_glue

cluster.df.add_glue$cluster_lib_new <- paste(cluster.df.add_glue$domain,'_',cluster.df.add_glue$cluster_lib)

cM <- table(cluster.df.add_glue$cluster_lib_new, cluster.df.add_glue$cluster)
         
#                1    2    3    4    5    6    7    8    9   10   11
#   atac _ c1  574 1640   27  163  561  452   29    6    3    1    1
#   atac _ c2    9   29 2530  813   24  119    5   14   40    0    0
#   atac _ c3  342  333   48 1140 1233   69    7    4    9    1    0
#   atac _ c4  891 1999  757 1090  239  759  118    4   38    0   15
#   atac _ c5 2241  111   12   10   98  118  237  129   17   10    0
#   atac _ c6 1548  119  409  138   31  239  102   50   14    1    0
#   atac _ c7  140   88   60   43   66  809   13    3    7  220    0
#   atac _ c8   41   15  128   62   16   26  724    6   31    1    0
#   atac _ c9    6    0    3    0    0    0    0  404   10    0    0
#   rna _ c1   866 2159  933  732  465  302  238   13   23    0    0
#   rna _ c10    7    2    8   14    7    3    9    0  156    0    0
#   rna _ c11    4    3    5    0    3    1    3  251    0    0    0
#   rna _ c2    10  701  809 2185  391   74   78    0    0    0    0
#   rna _ c3   893 1370   23  343 1157  212  112    4    2    0    0
#   rna _ c4  2732  194    7   12  137  154  112   35    7    0    0
#   rna _ c5     1   10 2301  194    0   10   76    0    2    0    0
#   rna _ c6  1053    0    4    0    2   21   77  662    3    0    0
#   rna _ c7     5    5   16    7    3 1065    8    1    2    0    0
#   rna _ c8     4    5    1    8    1    0  373    1    0    0    0
#   rna _ c9     4    1   10    5   10    3   53    1   12    0    0

             1    2    3    4    5    6    7    8    9   10   11
  atac _ c1  330 1479   18  134  996  459   30    6    3    1    1
  atac _ c2  933  119 2318    4   27  121    5   16   40    0    0
  atac _ c3 1519  305   43  120 1111   66    7    5    9    1    0
  atac _ c4  648 2622  475  223 1048  776   60    4   39    0   15
  atac _ c5    8   69   10 1591  864  102  190  122   17   10    0
  atac _ c6   79  296  367 1056  522  236   32   48   14    1    0
  atac _ c7   53   95   47   92  113  807   11    3    7  221    0
  atac _ c8   53   56   96   49   13   28  717    6   31    1    0
  atac _ c9    0    1    2    7    0    0    0  403   10    0    0
  rna _ c1   493 2586  665  245 1242  307  158   12   23    0    0
  rna _ c10   15    7    5    1    9    3   10    0  156    0    0
  rna _ c11    1    5    3    2    3    1    3  252    0    0    0
  rna _ c2  2229 1194  566    9   90   79   81    0    0    0    0
  rna _ c3   550 1146    6  465 1639  199  104    5    2    0    0
  rna _ c4    16  205    6 1394 1486  162   81   33    7    0    0
  rna _ c5   274   53 2165    2    0   11   87    0    2    0    0
  rna _ c6     0    2    1 1106   19   17   43  631    3    0    0
  rna _ c7     8    5   19    7    0 1062    8    1    2    0    0
  rna _ c8     5   16    1    8    1    0  361    1    0    0    0
  rna _ c9     9    4    8    1   12    2   50    1   12    0    0


#             1    2    3    4    5    6    7    8    9   10   11   12   13
#   atac _ c1  330 1471   18  134  637  459    8  359   30    6    3    1    1
#   atac _ c2  933   54 2318    4   27  121   65    0    5   16   40    0    0
#   atac _ c3 1519  190   43  120 1061   66  115   50    7    5    9    1    0
#   atac _ c4  648 1469  475  223  633  776 1153  415   60    4   39    0   15
#   atac _ c5    8   45   10 1591  225  102   24  639  190  122   17   10    0
#   atac _ c6   79  101  367 1056  420  236  195  102   32   48   14    1    0
#   atac _ c7   53   62   47   92   83  807   33   30   11    3    7  221    0
#   atac _ c8   53   11   96   49    3   28   45   10  717    6   31    1    0
#   atac _ c9    0    1    2    7    0    0    0    0    0  403   10    0    0
#   rna _ c1   493 1748  665  245  938  307  838  304  158   12   23    0    0
#   rna _ c10   15    1    5    1    7    3    6    2   10    0  156    0    0
#   rna _ c11    1    3    3    2    3    1    2    0    3  252    0    0    0
#   rna _ c2  2229  660  566    9   88   79  534    2   81    0    0    0    0
#   rna _ c3   550  985    6  465 1165  199  161  474  104    5    2    0    0
#   rna _ c4    16  133    6 1394  789  162   72  697   81   33    7    0    0
#   rna _ c5   274   13 2165    2    0   11   40    0   87    0    2    0    0
#   rna _ c6     0    1    1 1106   14   17    1    5   43  631    3    0    0
#   rna _ c7     8    4   19    7    0 1062    1    0    8    1    2    0    0
#   rna _ c8     5    5    1    8    0    0   11    1  361    1    0    0    0
#   rna _ c9     9    0    8    1    7    2    4    5   50    1   12    0    0


#            1   2   3   4   5   6   7   8   9  10
#   atac _ c1 366  94 591 449  52 150  27  19  52   1
#   atac _ c2  12 159  42 370 605 454   6  79  13   9
#   atac _ c3 578 424 154   6  41 193  29  42  70   2
#   atac _ c4 429  53 265  44   8   9 172   6  58   0
#   atac _ c5  43 533  95 130 152   3  25   1  66   0
#   atac _ c6  31  45  38  16  37  31 398   2  14   1
#   atac _ c7  15  10  24  13  10   7  41 219   7   0
#   atac _ c8   3   1   0   1   0   0   0   1   0 188
#   rna _ c1  142 234 707 271 156 113  53  68 106   8
#   rna _ c10   2 107 216 392 527   8  25  26  55   2
#   rna _ c2    0   0   1 406 177   0   2  11   1   1
#   rna _ c3  518 499 212  34  13 201  74  38 124   3
#   rna _ c4  963  33 103   9   2 339  90  61  21  17
#   rna _ c5   25 693  73  26  79   5  36   7  48   2
#   rna _ c6    1   0   1   2   2   2 477   7   0   0
#   rna _ c7  101   2   0   5   0 281   9   7   0   0
#   rna _ c8    1   4   1   2   7   2   0 150   1   2
#   rna _ c9    3   3   3   6   6   2   2   1   0  73


options(repr.plot.height=7.5,repr.plot.width=7.5)
pheatmap::pheatmap(cM,border_color = 'white',scale = 'none',fontsize_row = 15, fontsize_col = 15, display_numbers = TRUE,number_format = '%i')


#############use archR confusionMatrix method to predict integration cluster major cluster matching name in atac and rna ######
confusionMatrix <- function (i = NULL, j = NULL) #archR 
{
    ui <- unique(i)
    uj <- unique(j)
    m <- Matrix::sparseMatrix(i = match(i, ui), j = match(j, 
        uj), x = rep(1, length(i)), dims = c(length(ui), length(uj)))
    rownames(m) <- ui
    colnames(m) <- uj
    m
}

cM_new <- as.matrix(confusionMatrix(cluster.df.add_glue$cluster_lib_new, cluster.df.add_glue$cluster))


	6	1	8	5	4	3	2	9	7	10	11
rna _ c7	1062	8	1	0	7	19	5	2	8	0	0
rna _ c2	79	2229	0	90	9	566	1194	0	81	0	0
rna _ c4	162	16	33	1486	1394	6	205	7	81	0	0
rna _ c1	307	493	12	1242	245	665	2586	23	158	0	0
rna _ c6	17	0	631	19	1106	1	2	3	43	0	0
rna _ c9	2	9	1	12	1	8	4	12	50	0	0
rna _ c10	3	15	0	9	1	5	7	156	10	0	0
rna _ c3	199	550	5	1639	465	6	1146	2	104	0	0
rna _ c8	0	5	1	1	8	1	16	0	361	0	0
rna _ c11	1	1	252	3	2	3	5	0	3	0	0
rna _ c5	11	274	0	0	2	2165	53	2	87	0	0
atac _ c4	776	648	4	1048	223	475	2622	39	60	0	15
atac _ c3	66	1519	5	1111	120	43	305	9	7	1	0
atac _ c5	102	8	122	864	1591	10	69	17	190	10	0
atac _ c2	121	933	16	27	4	2318	119	40	5	0	0
atac _ c9	0	0	403	0	7	2	1	10	0	0	0
atac _ c6	236	79	48	522	1056	367	296	14	32	1	0
atac _ c8	28	53	6	13	49	96	56	31	717	1	0
atac _ c7	807	53	3	113	92	47	95	7	11	221	0
atac _ c1	459	330	6	996	134	18	1479	3	30	1	1


#     A matrix: 20 × 11 of type dbl
#              6	4	8	5	1	3	2	9	7	10	11
# rna _ c7	1065	7	1	3	5	16	5	2	8	0	0
# rna _ c2	74	2185	0	391	10	809	701	0	78	0	0
# rna _ c4	154	12	35	137	2732	7	194	7	112	0	0
# rna _ c1	302	732	13	465	866	933	2159	23	238	0	0
# rna _ c6	21	0	662	2	1053	4	0	3	77	0	0
# rna _ c9	3	5	1	10	4	10	1	12	53	0	0
# rna _ c10	3	14	0	7	7	8	2	156	9	0	0
# rna _ c3	212	343	4	1157	893	23	1370	2	112	0	0
# rna _ c8	0	8	1	1	4	1	5	0	373	0	0
# rna _ c11	1	0	251	3	4	5	3	0	3	0	0
# rna _ c5	10	194	0	0	1	2301	10	2	76	0	0
# atac _ c4	759	1090	4	239	891	757	1999	38	118	0	15
# atac _ c3	69	1140	4	1233	342	48	333	9	7	1	0
# atac _ c5	118	10	129	98	2241	12	111	17	237	10	0
# atac _ c2	119	813	14	24	9	2530	29	40	5	0	0
# atac _ c9	0	0	404	0	6	3	0	10	0	0	0
# atac _ c6	239	138	50	31	1548	409	119	14	102	1	0
# atac _ c8	26	62	6	16	41	128	15	31	724	1	0
# atac _ c7	809	43	3	66	140	60	88	7	13	220	0
# atac _ c1	452	163	6	561	574	27	1640	3	29	1	1


# 	           6	1	10	5	4	3	2	8	7	11	9	12	13
# rna _ c7	1062	8	1	0	7	19	4	0	1	2	8	0	0
# rna _ c2	79	2229	0	88	9	566	660	2	534	0	81	0	0
# rna _ c4	162	16	33	789	1394	6	133	697	72	7	81	0	0
# rna _ c1	307	493	12	938	245	665	1748	304	838	23	158	0	0
# rna _ c6	17	0	631	14	1106	1	1	5	1	3	43	0	0
# rna _ c9	2	9	1	7	1	8	0	5	4	12	50	0	0
# rna _ c10	3	15	0	7	1	5	1	2	6	156	10	0	0
# rna _ c3	199	550	5	1165	465	6	985	474	161	2	104	0	0
# rna _ c8	0	5	1	0	8	1	5	1	11	0	361	0	0
# rna _ c11	1	1	252	3	2	3	3	0	2	0	3	0	0
# rna _ c5	11	274	0	0	2	2165	13	0	40	2	87	0	0
# atac _ c4	776	648	4	633	223	475	1469	415	1153	39	60	0	15
# atac _ c3	66	1519	5	1061	120	43	190	50	115	9	7	1	0
# atac _ c5	102	8	122	225	1591	10	45	639	24	17	190	10	0
# atac _ c2	121	933	16	27	4	2318	54	0	65	40	5	0	0
# atac _ c9	0	0	403	0	7	2	1	0	0	10	0	0	0
# atac _ c6	236	79	48	420	1056	367	101	102	195	14	32	1	0
# atac _ c8	28	53	6	3	49	96	11	10	45	31	717	1	0
# atac _ c7	807	53	3	83	92	47	62	30	33	7	11	221	0
# atac _ c1	459	330	6	637	134	18	1471	359	8	3	30	1	1


# A matrix: 18 × 10 of type dbl
#             3	4	7	5	2	6	9	1	10	8
# rna _ c1	707	271	53	156	234	113	106	142	8	68
# rna _ c6	1	2	477	2	0	2	0	1	0	7
# rna _ c10	216	392	25	527	107	8	55	2	2	26
# rna _ c4	103	9	90	2	33	339	21	963	17	61
# rna _ c3	212	34	74	13	499	201	124	518	3	38
# rna _ c5	73	26	36	79	693	5	48	25	2	7
# rna _ c9	3	6	2	6	3	2	0	3	73	1
# rna _ c7	0	5	9	0	2	281	0	101	0	7
# rna _ c8	1	2	0	7	4	2	1	1	2	150
# rna _ c2	1	406	2	177	0	0	1	0	1	11
# atac _ c3	154	6	29	41	424	193	70	578	2	42
# atac _ c1	591	449	27	52	94	150	52	366	1	19
# atac _ c5	95	130	25	152	533	3	66	43	0	1
# atac _ c2	42	370	6	605	159	454	13	12	9	79
# atac _ c8	0	1	0	0	1	0	0	3	188	1
# atac _ c7	24	13	41	10	10	7	7	15	0	219
# atac _ c4	265	44	172	8	53	9	58	429	0	6
# atac _ c6	38	16	398	37	45	31	14	31	1	2


options(repr.plot.height=7.5,repr.plot.width=7.5)
res.p <- pheatmap::pheatmap(cM_new,border_color = 'white',scale = 'none',fontsize_row = 15, fontsize_col = 15, display_numbers = TRUE,number_format = '%i')
##the same with above method


pdf("pdfs/PLA-term-RNA-ATAC-glue.integration_cluster_match_rna-atac-cluster.pdf",height=7.5,width=7.5,useDingbats = FALSE)
res.p
dev.off()



#ggsave("pdfs/PLA-term-RNA-ATAC-glue.integration_cluster_match_rna-atac-cluster.pdf",height=7.5,width=7.5,useDingbats = FALSE)


preClust <- colnames(cM_new)[apply(cM_new, 1 , which.max)]
cbind(preClust, rownames(cM_new)) #Assignments

#6	rna _ c7
#1	rna _ c2
#5	rna _ c4
#2	rna _ c1
#4	rna _ c6
#7	rna _ c9
9	rna _ c10
#5	rna _ c3
#7	rna _ c8
#8	rna _ c11
#3	rna _ c5

#2	atac _ c4
#1	atac _ c3
#4	atac _ c5
#3	atac _ c2
#8	atac _ c9
#4	atac _ c6
#7	rna _ c8
#6	atac _ c7
#2	atac _ c1

#6	rna _ c7
#4	rna _ c2
#1	rna _ c4
#2	rna _ c1
#1	rna _ c6
#7	rna _ c9
#9	rna _ c10
#2	rna _ c3
#7	rna _ c8
#8	rna _ c11
#3	rna _ c5
#2	atac _ c4
#5	atac _ c3
#1	atac _ c5
#3	atac _ c2
#8	atac _ c9
#1	atac _ c6
#7	atac _ c8
#6	atac _ c7
#2	atac _ c1


#6	rna _ c7
#1	rna _ c2
#4	rna _ c4
#2	rna _ c1
#4	rna _ c6
#9	rna _ c9
#11	rna _ c10
#5	rna _ c3
#9	rna _ c8
#10	rna _ c11
#3	rna _ c5
#2	atac _ c4
#1	atac _ c3
#4	atac _ c5
#3	atac _ c2
#10	atac _ c9
#4	atac _ c6
#9	atac _ c8
#6	atac _ c7
#2	atac _ c1

# 3	rna _ c1
# 7	rna _ c6
# 5	rna _ c10
# 1	rna _ c4
# 1	rna _ c3
# 2	rna _ c5
# 10	rna _ c9
# 6	rna _ c7
# 8	rna _ c8
# 4	rna _ c2
# 1	atac _ c3
# 3	atac _ c1
# 2	atac _ c5
# 5	atac _ c2
# 10	atac _ c8
# 8	atac _ c7
# 1	atac _ c4
# 7	atac _ c6



##merged

#res = 0.7 (merge c7 to c2 and c8 to c5)
8 rna _ c11,atac _ c9
4 rna _ c6,atac _ c5,atac _ c6
2 rna _ c1,atac _ c4,atac _ c1
5 rna _ c4,rna _ c3
1 rna _ c2,atac _ c3
3 rna _ c5,atac _ c2
6 rna _ c7,atac _ c7
7 rna _ c9,rna _ c8,rna _ c8
9 rna _ c10

# 8 rna _ c11,atac _ c9
# 1 rna _ c6,rna _ c4,atac _ c5,atac _ c6
# 2 rna _ c1,rna _ c3,atac _ c4,atac _ c1
# 5 atac _ c3
# 4 rna _ c2,
# 3 rna _ c5,atac _ c2
# 6 rna _ c7,atac _ c7
# 7 rna _ c9,rna _ c8,atac _ c8
# 9 rna _ c10



#res = 0.7
3  rna _ c5,atac _ c2 
10 rna _ c11,atac _ c9
4  rna _ c6,rna _ c4,atac _ c5,atac _ c6
8  unmatched
2  rna _ c1,atac _ c4,atac _ c1 (part)
7   unmatched
5  rna _ c3,
1  rna _ c2,atac _ c3

11  rna _ c10,

6 rna _ c7,atac _ c7

9 rna _ c9,rna _ c8,atac _ c8

13 unmatched
12 unmathced


# 3	rna _ c1,atac _ c1
# 7	rna _ c6, atac _ c6
# 5	rna _ c10,atac _ c2
# 1	rna _ c4,rna _ c3, atac _ c3,atac _ c4
# 2	rna _ c5,atac _ c5
# 10	rna _ c9,atac _ c8
# 6	rna _ c7
# 8	rna _ c8,atac _ c7
# 4	rna _ c2

#9 ?	



###merge c7 to c2 and c8 to c5 and rename



cluster.df.add_glue$cluster_bk <- cluster.df.add_glue$cluster

clusterid <- cluster.df.add_glue$cluster

clusterid <- as.character(clusterid )

table(clusterid) #48673
 1   10   11   12   13    2    3    4    5    6    7    8    9 
7223 1549  377  235   16 6957 6821 6516 6100 4438 3308 3095 2038 

   1    2    3    4    5    6    7    8    9   10   11   12   13 
7223 6957 6821 6516 6100 4438 3308 3095 2038 1549  377  235   16 



map_list <- list(
   '1' = '1',
    '2' = '2',
    '3' = '3',
    '4' = '4',
    '5' = '5',
    '6' = '6',
    '7' = '2',
    '8' = '5',
    '9' = '9',
    '10' = '10',
    '11' = '11',
    '12' = '12',
    '13' = '13'



)

for(i in seq_len(length(clusterid)) ){
    
    clusterid[i] <- map_list[[ clusterid[i] ]]
    
    
}

table(clusterid) 
   1    10    11    12    13     2     3     4     5     6     9 
 7223  1549   377   235    16 10265  6821  6516  9195  4438  2038


map_list <- list(
   '1' = '1',
    '2' = '2',
    '3' = '3',
    '4' = '4',
    '5' = '5',
    '6' = '6',
    #'7' = '2',
    #'8' = '5',
    '9' = '7',
    '10' = '8',
    '11' = '9',
    '12' = '10',
    '13' = '11'



)


for(i in seq_len(length(clusterid)) ){
    
    clusterid[i] <- map_list[[ clusterid[i] ]]
    
    
}

table(clusterid) 
 1    10    11     2     3     4     5     6     7     8     9 
 7223   235    16 10265  6821  6516  9195  4438  2038  1549   377


clusterid <- factor(clusterid,levels = c('1','2','3','4','5','6','7','8','9','10','11'))

table(clusterid)
  1     2     3     4     5     6     7     8     9    10    11 
 7223 10265  6821  6516  9195  4438  2038  1549   377   235    16 

table(cluster.df.add_glue$cluster)
1    2    3    4    5    6    7    8    9   10   11   12   13 
7223 6957 6821 6516 6100 4438 3308 3095 2038 1549  377  235   16 

cluster.df.add_glue$cluster <- clusterid



#####


##################plot atac and rna cluster alignment######################
dotDistri_both = function (cluster = NULL, id1 = NULL, id2 = NULL){
    cluster$cluster_lib <- paste0('c',cluster$cluster_lib)
    #id1 for atac, id2 for rna
    cluster.sel1 = cluster[(cluster$domain == 'atac' & cluster$cluster_lib == id1),]
    n_sel1 = nrow(cluster.sel1)
    cluster.sel2 = cluster[(cluster$domain == 'rna' & cluster$cluster_lib == id2),]
    n_sel2 = nrow(cluster.sel2)
    
    ##stopifnot(sum(is.na(cluster.sel1)) == 0)
    ##stopifnot(sum(is.na(cluster.sel2)) == 0)
    
    ##write output for kdeplot
    ##fileout <- paste( "atac_c",id1,'_rna_c',id2,'.txt',sep=''   )
    ##write.table(x = rbind(cluster.sel1,cluster.sel2),file = fileout,sep='\t',row.names = TRUE,col.names = TRUE,quote = FALSE)
     
    color_atac = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel1)),'pink' ,'red'   )
    color_rna = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel2)),'lightblue' ,'navy'   ) 
   
    stopifnot(length(color_atac) == nrow(cluster.sel1))
    stopifnot(length(color_rna) == nrow(cluster.sel2))
    
    plot(cluster$UMAP1,cluster$UMAP2,pch = 16, type='p',col='grey',cex=0.5,xlab='UMAP1',ylab='UMAP2',main=paste("atac cells cluster ",id1,", n = ",n_sel1,"\nrna cells cluster ",id2,", n = ",n_sel2,sep=''),cex.main=2,cex.axis=2,cex.lab=2 ) 
    #points(cluster.sel2$UMAP1,cluster.sel2$UMAP2,pch = 16, cex=0.2,col='navy') #rna
    #points(cluster.sel1$UMAP1,cluster.sel1$UMAP2,pch = 16, cex=0.2,col='red') #atac

    points(cluster.sel1$UMAP1,cluster.sel1$UMAP2,pch = 16, cex=0.5,col=color_atac) #atac
    points(cluster.sel2$UMAP1,cluster.sel2$UMAP2,pch = 16, cex=0.5,col=color_rna ) #rna
    
    return(paste("cluster ",id1,", ",id2," ok. ",sep='') )
    #return(paste("cluster ",id1,", ",id2," ok. coords output to ",fileout,sep='') )
}



##result of FigR nuclei pairing result, confusion map with rna cluster vs atac cluster


pairing_cluster <- list(
#snRNA-seq      snATAC-atac
    #'10' = '-',
    'c11' = 'c9',
    'c6' = 'c5',    
    'c4' = 'c6',
    'c1' = 'c4',
    'c3' = 'c1',
    'c2' = 'c3',
    'c7' = 'c7',
    'c5' = 'c2',
    'c8' = 'c8'
    
    
# '8' = '7',
# '5' = '5',
# '10' = '2',
# '1' = '1',
# '4' = '4',
# '6' = '6',
# '3' = '3',
# '9' = '8'


)



pdf("pdfs/PLA-term-RNA-ATAC_rna_and_atac_cluster_pairwise_distribution..pdf",height=15,width=15,useDingbats = FALSE)

options(repr.plot.width=15,repr.plot.height=15)
par(mfrow=c(3,3)) #atac vs rna cluster_lib

dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c9',id2='c11') #CTB

dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c5',id2='c6') #STB nascent

dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c6',id2='c4') #STB FLT1 residual (intermediate)
dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c4',id2='c1') #STB FLT1 residual


dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c1',id2='c3') #STB PAPPA+ BMP1+
dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c3',id2='c2') #STB PAPPA+ BMP1+

dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c7',id2='c7') #STB PAPPA+ LAMA3+
dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c2',id2='c5') #STB PAPPA+ LAMA3+

dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c8',id2='c8') #syncytial knot like
#dotDistri_both(cluster = cluster.df.add_glue, id1 = 'c7',id2='c8') #STB FLT1+ residual

dev.off()



#####################evaluate cluster overlap percent by knn#################
library(FNN)
library(igraph)
##get_knnx method use k = 5



doKNN <- function(cluster = NULL,id1=NULL,id2=NULL,k=5,leg.pos = 'bottomright' ){
    ###start the real drive##
    #id1='5' #'8' #'4'
    #id2='3'#'4' #'1'
    #cluster = cluster.all.coord.substitute_snap
    cluster.sel1 = cluster[(cluster$domain == 'atac' & cluster$cluster_lib == id1),]
    #n_sel1 = nrow(cluster.sel1)
    if(grepl(pattern = ',',x = id2)){  
        ids <- unlist(strsplit(x = id2, split = ','))
        cluster.sel2 = cluster[(cluster$domain == 'rna' & cluster$cluster_lib %in% ids),]
    }else{
      cluster.sel2 = cluster[(cluster$domain == 'rna' & cluster$cluster_lib == id2),]
    }
    #n_sel2 = nrow(cluster.sel2)

    stopifnot(sum(is.na(cluster.sel1)) == 0)
    stopifnot(sum(is.na(cluster.sel2)) == 0)

    x1 <- cluster.sel1[,c('UMAP1','UMAP2')] #atac, reference
    #x1 <- cluster[cluster$type == 'atac',c('dim1','dim2')] #all atac, reference, not so good
    len1 = nrow(x1) #504 #424 #765 #605 #7394
    x2 <- cluster.sel2[,c('UMAP1','UMAP2')] #rna, query one cluster per time
    len2=nrow(x2) #1009

    
#     #####do some statistic test ?###
#     x <- x1
#     y <-  x2
#     ks <- peacock2(x, y)
#     cat ('two-dimensional Kolmogorov-Smirnov test,ks=:',ks,'\n')

#     par(mfrow=c(1,1))
#     options(repr.plot.width=2.5,repr.plot.height=2.5)
#     plot(x,pch=19,cex = 0.8, col = 'red', main = paste('two-dimensional Kolmogorov-Smirnov test,ks=:',ks,sep="") )
#     points(y,pch=19,cex = 0.8, col = 'navy')

    ##plot in xy coordinate to view dots nearest directly####
    par(mfrow=c(2,2))
   # cat('Before setting options:width=',options()$repr.plot.width,' height=',options()$repr.plot.height,'\n')
    options(repr.plot.width=10,repr.plot.height=10)
    #cat('After setting options:width=',options()$repr.plot.width,' height=',options()$repr.plot.height,'\n')
    
    plot(x1,cex=0.5,pch  = 19,col='red',main=paste('reference-atac','(c',id1,")"),cex.main=2,cex.lab=2)
    #text(x1,labels = 1:len1,col='grey',cex = 2 )
    plot(x2,cex=0.5,pch  = 19,col='navy',main=paste('query-rna','(c',id2,")"),cex.main=2)
    #text(x2,labels = 1:len2,col='red',cex = 2 )
    xall <- rbind(x1,x2)
    lenall <- nrow(xall)
    plot(xall,cex=0.5,pch  = 19,col=c(rep('red',len1),rep('navy',len2)),main='merge',cex.main=2 )#  c('black','grey'),each=10 ) )
    #text(xall,labels = c(1:len1,1:len2),col=c(rep('grey',len1),rep('red',len2)),cex = 2 )#rep(c('grey','red'),each=10 ) )


    ##query x2 from x1
    #set.seed(100) #graph varies even set seed, but topology is similar
    ###k = 5
    #par(mfrow=c(1,1))
    #options(repr.plot.width=5,repr.plot.height=5)
    nn = get.knnx(x1,x2,k)
    #nn = get.knnx(x1,x2,k) #for one point query
    # nn.df = data.frame(from = rep(1:nrow(nn$nn.index), k), 
    #                     to = as.vector(nn$nn.index), 
    #                     weight =  as.vector(nn$nn.dist)
    #                     #weight = 1/(1 + as.vector(nn$nn.dist))
    #                    )
    #test.df.nw = igraph::graph_from_data_frame(test.df, directed = TRUE)
    #plot(test.df.nw,main='graph of qeury to reference',cex.main=3)


    #all nearest from reference
    NN_index <- unique(as.vector(nn$nn.index) ) #644 #545 #1702
    #344/504 = 68% #396/424 = 93% 644/765 = 0.84
    perc_hit <- 100*(length(NN_index)/len1 )
    perc_hit <- round(perc_hit)
    #see hits dot cluster distribution in atac
    #par(mfrow=c(1,1))
    #options(repr.plot.width=8,repr.plot.height=8)
    plot(x1,cex=0.5,pch  = 19,col='grey',main=paste('query rna in \nreference atac. (',perc_hit,"%)",sep = ''),cex.main=2)
    points(x1[NN_index,],cex=0.5,pch  = 19,col='red')
    points(x2,cex=0.5,pch  = 19,col='navy')
    #legend("topright",legend = c('scATAC','scATAC(outside)','scRNA'),
    legend(leg.pos,legend = c('scATAC','scATAC(outside)','scRNA'),
           fill = c('red','grey','navy'),
           x.intersp = 0.8, #x space
           y.intersp = 0.8, #y space
           cex = 1,
           text.width = 3.2 #overall text width
          )
    

    return(paste('atac:c',id1,' ','rna:c',id2,' ok. ','knn hits: ',perc_hit,"%",sep='') )
}


cluster.df.add_glue.bk <- cluster.df.add_glue

##remove cluster prefix
cluster.df.add_glue$cluster_lib_bk <- cluster.df.add_glue$cluster_lib
cluster.df.add_glue$cluster_lib <- gsub('c','',cluster.df.add_glue$cluster_lib)

# table(cluster.df.add_glue$cluster_lib)
# 1   10    2    3    4    5    6    7    8    9 
# 3659 1360 2348 3255 2682 2042 1105  751  364   99 

# table(cluster.df.add_glue$cluster_lib_bk)
# c1  c10   c2   c3   c4   c5   c6   c7   c8   c9 
# 3659 1360 2348 3255 2682 2042 1105  751  364   99


##for STB subtype, atac vs rna
options(repr.plot.width = 7.5, repr.plot.height = 7.5)
    #id2:rna      id1:atac
    'c11' = 'c9',
    'c6' = 'c5',    
    'c4' = 'c6',
    'c1' = 'c4',
    'c3' = 'c1',
    'c2' = 'c3',
    'c7' = 'c7',
    'c5' = 'c2',
    'c8' = 'c8'



##CTB
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c9-vs-rna_c11.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='9',id2='11',k=5,leg.pos = 'topright' ) #77% 
dev.off()


#STB Nascent
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c5-vs-rna_c6.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='5',id2='6',k=5,leg.pos = 'topright' ) #84% 
dev.off()


#STB PAPPA+ CSH1+
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c6-vs-rna_c4.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='6',id2='4',k=5,leg.pos = 'topright' ) #82% 
dev.off()

#STB PAPPA+ CSH1+
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c4-vs-rna_c1.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='4',id2='1',k=5,leg.pos = 'topright' ) #72% 
dev.off()

##STB PAPPA+ LAMA3+
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c1-vs-rna_c3.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='1',id2='3',k=5,leg.pos = 'topright' ) #83% 
dev.off()

##STB PAPPA+ ADAMTSL1+
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c3-vs-rna_c2.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='3',id2='2',k=5,leg.pos = 'topright' ) #88% 
dev.off()


##syncytial knot like
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c7-vs-rna_c7.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='7',id2='7',k=5,leg.pos = 'topright' ) #71% 
dev.off()

##STB FLT1 residual
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c2-vs-rna_c5.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='2',id2='5',k=5,leg.pos = 'topright' ) #66% 
dev.off()

##STB FLT1 residual apoptosis
pdf('pdfs/knn_quantify/nearest-neighbor-ratio.atac_c8-vs-rna_c8.pdf',width=7.5,height=7.5,useDingbats = FALSE)
  doKNN(cluster = cluster.df.add_glue,id1='8',id2='8',k=5,leg.pos = 'topright' ) #66% 
dev.off()


##summary

#atac    rna    integration   overlap-ratio
atac:c9 rna:c11      c8            65% 
atac:c5 rna:c6       c4            52%
atac:c6 rna:c4       c4            75%
atac:c4 rna:c1       c2            95%
atac:c1 rna:c3       c5            95%
atac:c3 rna:c2       c1            63%
atac:c7 rna:c7       c6            65%
atac:c2 rna:c5       c3            81%
atac:c8 rna:c8       c7            60%
-       rna:c10      c11            -
#-       rna:c9      c7             -


