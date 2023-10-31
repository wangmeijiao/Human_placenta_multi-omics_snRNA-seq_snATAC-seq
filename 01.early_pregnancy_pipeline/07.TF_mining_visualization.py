


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import matplotlib.patches as mpatches

import matplotlib.colors as mcolors

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import math



##https://stackoverflow.com/questions/59695334/custom-color-palette-in-seaborn
##use yellowbrick palette for customized colorset in seaborn
from yellowbrick.style.palettes import PALETTES, SEQUENCES, color_palette

color_use = SEQUENCES['Spectral'][9]
#color_palette(color_use)

#cmap = plt.get_cmap('Reds')
#cmap = plt.get_cmap('RdBu_r')
cmap = sns.color_palette(color_use,9)

# import sys
# sys.path.append("/home/mjwang/lib/python/")
# import color_utils

# color_utils.draw_cmap(cmap)


# data = [24, 24, 24, 16, 16, 2, 2, 2]
# x = list(range(0, len(data)))
# y = list(range(0, 25))

# plt.grid(True)
# plt.scatter(x, data, marker='o', c='g', s=100)
# plt.yticks(y)
# plt.xticks(x)


# plt.show()

##add cell name and cell color (only te clusters)
# map_cellname_rna = { #snRNA-seq
#     '8':'eCTB_Pro',
#     '9':'eCTB1', 
#     '5':'eCTB2',
#     '11':'eCTB_Fusion', #in full cell clusters
#     '6':'eSTB_Naive SH3TC2',
#     '1':'eSTB1',
#     '3':'eSTB2 PAPPA+',
#     '4':'eSTB3',
#     '7':'eSTB4 FLT1+',
#     '2':'eSTB5',
#     '12' : 'EVT',
#     '10' : 'STR1',
#     '13' : 'STR2',
#     '14' : 'STR3',
#     '15' : 'STR4',
#     '16' : 'STR5'
    
# }

map_cellname_rna_te = { #snRNA-seq
  '9' : 'CTB1',
  '6': 'CTB2',
  '11': 'CTB Fusion',
  '8': 'STB Nascent',
  '1': 'STB Premature2',
  '3': 'STB Mature2',
  '4':'STB mixed',
  '2': 'STB Premature1',
  '10': 'STB Mature1'
#   '11'= 'darkgreen',
#   '9'= '#053061'
#)

    
#     '8':'eCTB_Pro',
#     '9':'eCTB1', 
#     '5':'eCTB2',
#     '10':'eCTB_Fusion',
#     '6':'eSTB_Naive SH3TC2',
#     '1':'eSTB1',
#     '3':'eSTB2 PAPPA+',
#     '4':'eSTB3',
#     '7':'eSTB4 FLT1+',
#     '2':'eSTB5'
    
}

# map_cellname_atac = { #snATAC-seq
#     '6':'eCTB1(6)',
#     '3':'eCTB2(3)', 
#     '9':'eCTB_Fusion(9)',
#     '5':'eSTB_Naive SH3TC2(5)',
#     '7':'eSTB3(7)',
#     '1': 'eSTB4 FLT1+(1)',
#     '8': 'eSTB5(8)',
#     '2':'eSTB1 PAPPA+(2)',
#     '4':'eSTB2(4)',
#     '12' : 'EVT(12)',
#     '10' : 'STR1(10)',
#     '13' : 'STR2(13)',
#     '15' : 'STR3(15)',
#     '11' : 'STR4(11)',
#     '14' : 'STR5(14)'
    
    
# }


map_cellname_atac_te = { #snATAC-seq
    

    '1' : 'CTB1',
    '7' : 'CTB2',
    '9' : 'CTB Fusion',
    '2' : 'STB Nascent',
    '8' : 'STB Premature2',
    '3' : 'STB Mature2',
    #'#ECDB68',
    #'#AFD16A',
    #'#60C589',
    '6' : 'STB mixed',#'#2BB9AD',
    '4' : 'STB Premature1',#'#337C99',
    '5' : 'STB Mature1'


    
#     '6':'eCTB1(c6)',
#     '3':'eCTB2(c3)', 
#     '9':'eCTB_Fusion(c9)',
#     '5':'eSTB_Naive SH3TC2(c5)',
#     '7':'eSTB3(c7)',
#     '1': 'eSTB4 FLT1+(c1)',
#     '8': 'eSTB5(8)',
#     '2':'eSTB1 PAPPA+(c2)',
#     '4':'eSTB2(c4)'
    
    
}



#'6','3','9','5','7','1','8','2','4'
## '#74add1','#4575b4','darkgreen','#ffffbf','#fdae61','#f46d43','#fee090','#d73027','#a50026'

color_set = {
    'c1' : '#3D3F69',
    'c7' : '#4776B2',
    'c9' : 'darkgreen',
    'c2' : '#FCC140',
    'c8' : '#F95944',
    'c3' : '#8D1541',
    'c6' : '#FB8D3C',
    'c4' : '#C31240',
    'c5' : '#4F1C47'
#            'c6'='#74add1',
#            'c3'='#4575b4',
#            'c9'='darkgreen',
#            'c5'='#ffffbf',
#            'c7'='#fdae61',
#            'c1'='#f46d43',
#            'c8'='#fee090',
#            'c2'='#d73027',
#            'c4'='#a50026'
}


# color_set = {'c6':'#74add1',
#            'c3':'#4575b4',
#            'c9':'darkgreen',
#            'c5':'#ffffbf',
#            'c7':'#fdae61',
#            'c1':'#f46d43',
#            'c8':'#fee090',
#            'c2':'#d73027',
#            'c4':'#a50026'}

outdir = 'result_pycisTarget_darfull_auc0.005_icisTarget_darp2g_auc0.005'

#tfmotif_data = pd.read_csv(outdir+"/tfmotif_align_long.txt",sep = '\t' )
tfmotif_data = pd.read_csv(outdir+"/tfmotif_align_long.cutoff.q40.txt",sep = '\t' ) 

#color_set_order = [color_set[i]  for i in ['c6','c3','c9','c5','c1','c2','c4'] ]
color_set_order = [color_set[i]  for i in ['c1','c7','c9','c2','c8','c3','c6','c4','c5'] ]

tf_len = len(tfmotif_data['y_id'].value_counts() )

#cell_names =  [map_cellname_atac_te[i] for i in ['6','3','9','5','1','2','4'] ]
cell_names =  [map_cellname_atac_te[i] for i in ['1','7','9','2','8','3','6','4','5'] ]

#########

# d = { 'x': x,  'y' : data, 'x_id' : range(1,len(x)+1) , 'y_id': range(1,len(x)+1)  ,'nes' : range(1,len(x)+1) , 'expr' : range(1,len(x)+1), 'keep' : range(1,len(x)+1)  }
# data_df = pd.DataFrame(d,index = [f"{'rowid' + str(i)}" for i in range(1,len(x)+1)]  )

#tips = sns.load_dataset("tips")
#tips.head()

# def dotHeatmap(data):
    
#     ax = sns.scatterplot(
#         data=data_df, x="x", y="y", hue="size", size="size",
#         #sizes=(20, 200), 
#         legend=False
#     )


#     ax.grid(b=True, which='major',color='black',linewidth=0.2)
#     ax.grid(b=True, which='minor',color='black',linewidth=0.2)

#################

plt.rcParams["figure.figsize"] = [3.5,7.5]
##plt.rcParams["figure.figsize"] = [3.5,20.5]

plt.figure(figsize=(3.5,7.5))
##plt.figure(figsize=(3.5,20.5))
kwargs  =   {'edgecolor':tfmotif_data['edge'], # for edge color
             'linewidth':1, # line width of spot
            }

ax = sns.scatterplot(
    #data=tfmotif_data, x="x_id", y="y_id", size="NES", hue="expression.transform",
    data=tfmotif_data, x="x_id", y="y_id", hue="NES", size="expression.transform",
    sizes=(0, 150), 
    #marker = 'o',
    #linewidth=2,
    #hue_norm=(-0.8,0.8),
    palette = 'RdBu_r',#'Reds',
    #edgecolors = 'blue',
    alpha = 1,
    legend='brief',
    #clip_on = False,
    **kwargs
)


ax.grid(b=True, which='major',color='black',linewidth=0.15,alpha = 0.25)
ax.grid(b=True, which='minor',color='black',linewidth=0.15,alpha = 0.25)

ax.xaxis.set_ticks_position('top')
ax.yaxis.set_ticks_position('right')

ax.set_xticklabels( labels = cell_names,
                    #fontdict = {'horizontalalignment':'center','verticalalignment':'bottom'},
                    rotation = 90, 
                    ha = 'center',
                    va = 'bottom'
                  )

#ax.xaxis.set_tick_params( labelrotation = 45, grid_color = 'grey',length = 0,pad = 15, labelsize='small')
ax.xaxis.set_tick_params(  grid_color = 'grey',length = 0,pad = 18, labelsize='small')

ax.xaxis.label.set_text("")
ax.yaxis.label.set_text("")



legd = ax.legend(bbox_to_anchor=(2,1), loc='upper right',ncol=1)
legd.remove()

##try to add annotation rects below colnames
#https://datavizpyr.com/how-to-draw-a-rectangle-on-a-plot-in-matplotlib/

xlim = ax.get_xlim() #(-0.30000000000000004, 6.3)
ylim = ax.get_ylim() #(32.55, -1.55)

(x_min,x_max) = xlim
(y_min,y_max) = ylim

xticks = [(tick - x_min)/(x_max - x_min) for tick in ax.get_xticks()]
yticks = [(tick - y_min)/(y_max - y_min) for tick in ax.get_yticks()]

# left, bottom, width, height = (-0.3, -2.3, 0.7, 0.8)
# rect=mpatches.Rectangle((left,bottom),width,height, 
#                         fill=True,
#                         color=color_set['c6'],
#                         alpha = 0.8,
#                        linewidth=0,
#                        clip_on=False,
#                        #facecolor="red"
#                        )
# plt.gca().add_patch(rect)


# left, bottom, width, height = (0.5, -2.3, 1, 1)
# rect=mpatches.Rectangle((left,bottom),width,height, 
#                         fill=True,
#                         color="red",
#                         alpha = 0.8,
#                        linewidth=0,
#                        clip_on=False,
#                        facecolor="red"
#                        )
# plt.gca().add_patch(rect)

##column annotation
for i in range(1,10):
    print(i)

    left, bottom, width, height = (i-1.5, -2.5, 1, 1)
    if left < xlim[0]:
        left = xlim[0]
        width = width-(xlim[0]-(i-1.5))
    
    if left + width > xlim[1]:
        width = xlim[1] - left
    
    rect=mpatches.Rectangle((left,bottom),width,height, 
                            fill=True,
                            color=color_set_order[i-1],
                            alpha = 0.8,
                           linewidth=0,
                           clip_on=False,
                           #facecolor="red"
                           )
    plt.gca().add_patch(rect)

# # ##legend marker size 
# plt.gca().set_aspect('equal', adjustable='box')

# cirs = mpatches.Circle((xlim[1]*0.8,ylim[0]*0.8),radius = 1,fill = True, color = 'black',hatch = 'o')
# plt.gca().add_patch(cirs)

###legend: map data to colorbar

data_norm =plt.Normalize( vmin = tfmotif_data['NES'].min(), vmax = math.ceil(tfmotif_data['NES'].max()))

cmap=plt.cm.get_cmap('RdBu_r',6)

scalarmappable = plt.cm.ScalarMappable(norm=data_norm, cmap=cmap)

scalarmappable.set_array(tfmotif_data['NES'])


#cbar = plt.colorbar(scalarmappable,ax=ax, orientation = 'horizontal',shrink = 0.8,pad = 0.02, label = 'Normalized expression')


axins = inset_axes(ax,
                    width="80%",  
                    height="3%",
                    loc='lower center',
                    borderpad=-2
                   )
cbar = plt.colorbar(scalarmappable, cax=axins, orientation="horizontal",label = 'TF motif normalized enrichment score (NES)')



#plt.tight_layout()

#https://stackoverflow.com/questions/4042192/reduce-left-and-right-margins-in-matplotlib-plot
plt.subplots_adjust( #need this if output pdf, turn off this if only need png of juputerlab
    #top=0.85, #scale the edge position related to fraction of height
    top=0.95,
    bottom=0.1,
    left=0.042,
    right=0.8,
    #hspace=0.5,
    #wspace=0.5
)


#plt.savefig(outdir+'/TFmotif-nes-expr.dotheatmap.pdf',format = 'pdf',transparent = False)
plt.savefig(outdir+'/TFmotif-nes-expr.dotheatmap.cutoff.q40.pdf',format = 'pdf',transparent = False)

#plt.savefig('TFmotif-nes-expr.dotheatmap.with_dotsize_legend.pdf',format = 'pdf',transparent = False)




