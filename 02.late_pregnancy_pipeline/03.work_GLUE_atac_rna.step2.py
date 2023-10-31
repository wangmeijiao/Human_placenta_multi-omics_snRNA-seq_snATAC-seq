
import anndata
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams
from matplotlib import pyplot as plt

anndata.__version__ #0.8.0

import pickle

import scanpy.external as sce #for harmony

import re

from umap import UMAP

import plotly.express as px
#import plotly
#plotly.plot()
import plotly.io as pio
pio.renderers.default = "png"

import plotly.graph_objects as go



color_good = ["#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" ]




#Stage 2: Model training


pkg_list = ['anndata','networkx','scanpy','scglue','numpy','pandas','leidenalg','louvain']
def print_header(pkg_list):
    for pkg in pkg_list:
            try:
                imp = __import__(pkg)
                print (pkg + '==' + imp.__version__, end = ", ")
            except (ImportError, AttributeError):
                pass


print_header(pkg_list)
anndata==0.8.0, networkx==2.8.4, scanpy==1.8.2, scglue==0.2.3, numpy==1.22.4, pandas==1.5.3,louvain==0.7.1

import leidenalg
leidenalg.version #0.9.1


scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


#Read preprocessed data¶


rna = anndata.read_h5ad("rna_preprocessed.h5ad")
atac = anndata.read_h5ad("atac_preprocessed.h5ad")
graph = nx.read_graphml("prior.graphml.gz")

#Configure data

scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)

# For the scRNA-seq data, we use the previously backed up raw counts in the “counts” layer, and use the PCA embedding as the first encoder transformation.

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

# For the scATAC-seq data, the raw counts are just atac.X, so it’s unnecessary to specify use_layer. We use the LSI embedding as the first encoder transformation.


##subset graph to retain highly variable features only
graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))


##Build and train GLUE model

glue = scglue.models.fit_SCGLUE( ##slow #from ~  1:00 to 13:00, 12h!
    {"rna": rna, "atac": atac}, graph,
    fit_kws={"directory": "glue"}
)

glue.save("glue.dill")

##
# SCGLUE model with the following network and trainer:
# SCGLUE(
#   (g2v): GraphEncoder(
#     (conv): GraphConv()
#     (loc): Linear(in_features=50, out_features=50, bias=True)
#     (std_lin): Linear(in_features=50, out_features=50, bias=True)
#   )
#   (v2g): GraphDecoder()
#   (x2u): ModuleDict(
#     (rna): NBDataEncoder(
#       (linear_0): Linear(in_features=100, out_features=256, bias=True)
#       (act_0): LeakyReLU(negative_slope=0.2)
#       (bn_0): BatchNorm1d(256, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
#       (dropout_0): Dropout(p=0.2, inplace=False)
#       (linear_1): Linear(in_features=256, out_features=256, bias=True)
#       (act_1): LeakyReLU(negative_slope=0.2)
#       (bn_1): BatchNorm1d(256, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
#       (dropout_1): Dropout(p=0.2, inplace=False)
#       (loc): Linear(in_features=256, out_features=50, bias=True)
#       (std_lin): Linear(in_features=256, out_features=50, bias=True)
#     )
#     (atac): NBDataEncoder(
#       (linear_0): Linear(in_features=100, out_features=256, bias=True)
#       (act_0): LeakyReLU(negative_slope=0.2)
#       (bn_0): BatchNorm1d(256, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
#       (dropout_0): Dropout(p=0.2, inplace=False)
#       (linear_1): Linear(in_features=256, out_features=256, bias=True)
#       (act_1): LeakyReLU(negative_slope=0.2)
#       (bn_1): BatchNorm1d(256, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
#       (dropout_1): Dropout(p=0.2, inplace=False)
#       (loc): Linear(in_features=256, out_features=50, bias=True)
#       (std_lin): Linear(in_features=256, out_features=50, bias=True)
#     )
#   )
#   (u2x): ModuleDict(
#     (rna): NBDataDecoder()
#     (atac): NBDataDecoder()
#   )
#   (du): Discriminator(
#     (linear_0): Linear(in_features=50, out_features=256, bias=True)
#     (act_0): LeakyReLU(negative_slope=0.2)
#     (dropout_0): Dropout(p=0.2, inplace=False)
#     (linear_1): Linear(in_features=256, out_features=256, bias=True)
#     (act_1): LeakyReLU(negative_slope=0.2)
#     (dropout_1): Dropout(p=0.2, inplace=False)
#     (pred): Linear(in_features=256, out_features=2, bias=True)
#   )
#   (prior): Prior()
# )

# SCGLUETrainer(
#   lam_graph: 0.02
#   lam_align: 0.05
#   vae_optim: RMSprop (
#   Parameter Group 0
#     alpha: 0.99
#     centered: False
#     eps: 1e-08
#     lr: 0.0002
#     momentum: 0
#     weight_decay: 0
#   )
#   dsc_optim: RMSprop (
#   Parameter Group 0
#     alpha: 0.99
#     centered: False
#     eps: 1e-08
#     lr: 0.0002
#     momentum: 0
#     weight_decay: 0
#   )
#   freeze_u: False
# )


# import dill

# with open('glue.dill','r') as fh: #error
#     gl = dill.load(fh)


##reload glue model
glue = scglue.models.load_model("glue.dill") 

glue.compile() #register trainer



SCGLUE model with the following network and trainer:

# SCGLUE(
#   (g2v): GraphEncoder(
#     (conv): GraphConv()
#     (loc): Linear(in_features=50, out_features=50, bias=True)
#     (std_lin): Linear(in_features=50, out_features=50, bias=True)
#   )
#   (v2g): GraphDecoder()
#   (x2u): ModuleDict(
#     (rna): NBDataEncoder(
#       (linear_0): Linear(in_features=100, out_features=256, bias=True)
#       (act_0): LeakyReLU(negative_slope=0.2)
#       (bn_0): BatchNorm1d(256, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
#       (dropout_0): Dropout(p=0.2, inplace=False)
#       (linear_1): Linear(in_features=256, out_features=256, bias=True)
#       (act_1): LeakyReLU(negative_slope=0.2)
#       (bn_1): BatchNorm1d(256, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
#       (dropout_1): Dropout(p=0.2, inplace=False)
#       (loc): Linear(in_features=256, out_features=50, bias=True)
#       (std_lin): Linear(in_features=256, out_features=50, bias=True)
#     )
#     (atac): NBDataEncoder(
#       (linear_0): Linear(in_features=100, out_features=256, bias=True)
#       (act_0): LeakyReLU(negative_slope=0.2)
#       (bn_0): BatchNorm1d(256, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
#       (dropout_0): Dropout(p=0.2, inplace=False)
#       (linear_1): Linear(in_features=256, out_features=256, bias=True)
#       (act_1): LeakyReLU(negative_slope=0.2)
#       (bn_1): BatchNorm1d(256, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True)
#       (dropout_1): Dropout(p=0.2, inplace=False)
#       (loc): Linear(in_features=256, out_features=50, bias=True)
#       (std_lin): Linear(in_features=256, out_features=50, bias=True)
#     )
#   )
#   (u2x): ModuleDict(
#     (rna): NBDataDecoder()
#     (atac): NBDataDecoder()
#   )
#   (du): Discriminator(
#     (linear_0): Linear(in_features=50, out_features=256, bias=True)
#     (act_0): LeakyReLU(negative_slope=0.2)
#     (dropout_0): Dropout(p=0.2, inplace=False)
#     (linear_1): Linear(in_features=256, out_features=256, bias=True)
#     (act_1): LeakyReLU(negative_slope=0.2)
#     (dropout_1): Dropout(p=0.2, inplace=False)
#     (pred): Linear(in_features=256, out_features=2, bias=True)
#   )
#   (prior): Prior()
# )

# SCGLUETrainer(
#   lam_graph: 0.02
#   lam_align: 0.05
#   vae_optim: RMSprop (
#   Parameter Group 0
#     alpha: 0.99
#     centered: False
#     eps: 1e-08
#     lr: 0.002
#     momentum: 0
#     weight_decay: 0
#   )
#   dsc_optim: RMSprop (
#   Parameter Group 0
#     alpha: 0.99
#     centered: False
#     eps: 1e-08
#     lr: 0.002
#     momentum: 0
#     weight_decay: 0
#   )
#   freeze_u: False
# )







#Check integration diagnostics
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, graph,
    count_layers={"rna": "counts"}
)
dx
	n_meta	consistency
0	10	0.127772
1	20	0.099819
2	50	0.102143
3	100	0.085907
4	200	0.073060



_ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")
#Empirically, it is safe to assume that the integration is reliable if the curve is above the 0.05 line



#Apply model for cell and feature embedding (encode_data with glue model), data imputation?

rna.obsm["X_glue"] = glue.encode_data("rna", rna) #23981 x 50 #9331 x 50
atac.obsm["X_glue"] = glue.encode_data("atac", atac) #24692 x 50 #8334 x 50

##merge

##add batch of rna
#rna.obs['sample'] =  rna.obs.index
batch =['foo'] * len(rna.obs.index)

for i in range(len(rna.obs.index)):#0-9330
    #print(i,rna.obs.index[i])
    cellid = rna.obs.index[i]
    if bool(re.search('-1$',cellid)):
        batch[i] = 'rna_D1'
    elif bool(re.search('-2$',cellid)):
        batch[i] = 'rna_D2'
    elif bool(re.search('-3$',cellid)):
        batch[i] = 'rna_D3'
    elif bool(re.search('-5$',cellid)):
        batch[i] = 'rna_D5'
    elif bool(re.search('-6$',cellid)):
        batch[i] = 'rna_D6'
    elif bool(re.search('-9$',cellid)):
        batch[i] = 'rna_D9'
    else:
        print('error at ' + str(i))
    
rna.obs['batch'] = batch

#all(rna.obs['batch'] == rna.obs['sample']) #True


rna.obs['batch'].value_counts().sort_index()
rna_D1    2563
rna_D2    3566
rna_D3    3720
rna_D5    3886
rna_D6    4545
rna_D9    5701
#checked



##add batch of atac
batch =['foo'] * len(atac.obs.index)

for i in range(len(atac.obs.index)):#0-9330
    #print(i,rna.obs.index[i])
    cellid = atac.obs.index[i]
    if bool(re.search('-1$',cellid)):
        batch[i] = 'atac_D1'
    elif bool(re.search('-2$',cellid)):
        batch[i] = 'atac_D2'
    elif bool(re.search('-3$',cellid)):
        batch[i] = 'atac_D3'
    elif bool(re.search('-5$',cellid)):
        batch[i] = 'atac_D5'
    elif bool(re.search('-6$',cellid)):
        batch[i] = 'atac_D6'
    elif bool(re.search('-9$',cellid)):
        batch[i] = 'atac_D9'
    else:
        print('error at ' + str(i))
    
atac.obs['batch'] = batch

#all(atac.obs['batch'] == atac.obs['sample']) #True

atac.obs['batch'].value_counts().sort_index()
atac_D1    4999
atac_D2    3550
atac_D3    3897
atac_D5    3248
atac_D6    4625
atac_D9    4373
#checked


combined = anndata.concat([rna, atac]) #23981 + 24692 = 48673 #9331 () + 8334 () = 17665

AnnData object with n_obs × n_vars = 48673 × 0
    obs: 'cell_type', 'domain', 'balancing_weight', 'batch'
    obsm: 'X_umap', 'X_glue'
    
    
combined.obs
combined.obs.batch.value_counts().sort_index()
atac_D1    4999
atac_D2    3550
atac_D3    3897
atac_D5    3248
atac_D6    4625
atac_D9    4373

rna_D1     2563
rna_D2     3566
rna_D3     3720
rna_D5     3886
rna_D6     4545
rna_D9     5701

# rna_D2     5684
# atac_D1    4913
# rna_D1     3647
# atac_D2    3421

combined.to_df() #cell x gene, empty
combined.obs.to_csv('combined.obs.txt',sep='\t')
combined.var_names #empty



##batch effect removal

sce.pp.harmony_integrate(combined, 'batch', basis='X_glue', adjusted_basis='X_glue_harmony') #after performing PCA but before computing the neighbor graph, rewrite the X_pca


#Converged after 9 iterations



#harmonypy - INFO - Converged after 3 iterations
#Converged after 5 iterations
#sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
#sc.tl.louvain(adata,resolution=1.4)
#sc.tl.leiden(adata,resolution=1.25)


##plot
sc.pp.neighbors(combined, use_rep="X_glue_harmony", metric="cosine")
sc.tl.umap(combined)
sc.pl.umap(combined, color=["cell_type", "domain"], wspace=0.65)



#sc.pp.neighbors(combined, n_neighbors=15, n_pcs=30)
#sc.tl.louvain(combined,resolution=1)
sc.tl.leiden(combined,resolution=0.7)
#sc.pl.umap(combined, color=["leiden", "louvain"], wspace=0.65)
sc.pl.umap(combined, color=["leiden"], wspace=0.65)




##########tuning umap######
###umap tunning with parameters: min_dist (default 0.1), spread (default 1.0), random_state (default 0)#######
dim_use = combined.obsm['X_glue_harmony'] #29721 x 18 #11114 x 18
groups = combined.obs['leiden']

random_state = 0 #123
n_comps = 2
#min_dist = 0.1 
#spread = 1.0

config = dict({'staticPlot': True})

from random import sample

###choose:
mdist_sel = 0.3
spread_sel = 1.5
#random_state = 0
#n_comps = 2
for random_state in sample(range(10000),5):
#    for mdist in [0.1,0.2,0.3,0.4,0.5]:
     for mdist in [mdist_sel]:
#        for spread in [0.5,1.0,1.5,2.0]:
         for spread in [spread_sel]:
            umap = UMAP(random_state=random_state, #array
                        n_components=n_comps,
                        min_dist = mdist,
                        spread = spread
                       ).fit_transform(dim_use)

            if umap.shape[0] != len(groups):
                print('error: umap length diff with groups')


            df = pd.DataFrame({
                "UMAP-1": umap[:, 0],
                "UMAP-2": umap[:, 1],
                'leiden': groups,
            })
            fig = px.scatter(
                df, x="UMAP-1", y="UMAP-2", color='leiden',
                color_discrete_sequence=color_good,#px.colors.qualitative.Dark24,
                title = 'umap tunning: mdist ' + str(mdist) + ' spread ' + str(spread) + ' random_state ' + str(random_state),
                width = 550,
                height = 450,
                category_orders = {'leiden':groups.value_counts().index.to_list(),
                                   #'cluster':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
                                  }

            )

            fig.update_traces(
              marker_size=2,
              marker={"opacity": 1},
              #legendgrouptitle_font_color='red', 
              #legendgrouptitle_text= 'cluster',
              #selector=dict(type='scatter')
            )

            fig.update_layout(
              template="simple_white",
              legend= {'itemsizing': 'constant'},
              #title_text='test' #ok!
              margin=dict(l=20, r=20, t=60, b=60),
              xaxis={'visible': False, 'showticklabels': False},
              yaxis={'visible': False, 'showticklabels': False}
            )

            fig.update_annotations(
              #'font': {'color':'black','size':45},
                x=10,
                y=15,
                #'showarrow':False,
                text="A very clear explanation",
                #'textangle':0,
                #'xanchor':'left',
                xref="paper",
                yref="paper"

              # {'font': {'color':'black','size':45},
              #   'x':10,
              #   'y':5,
              #   'showarrow':False,
              #   'text':"A very clear explanation",
              #   'textangle':0,
              #   'xanchor':'left',
              #   'xref':"paper",
              #   'yref':"paper"
              # }

            )

            fig.show(interactive=False,config=config,width = 550,height = 450)#will cause very large mem in chrome!

            
            
            
        #img_bytes = fig.to_image(format="png", width=350, height=350, scale=2)
        #Image(img_bytes)
        
        #fig.write_image("pdfs/umap.pdf")

        
        
#         ##quick annotation within tunning
#         if data.obsm['X_umap'].shape[0] == umap.shape[0]:#True
#             gene_matrix.obsm["X_umap"] = umap
#             gene_matrix.obs['leiden'] = data.obs['leiden']
            
#         else:
#             print('umap length diff')
        
#         thresh = {
#             'ERVFRD-1':['p30','p100'],
#             'PAPPA':['p65','p100'],
#             'LAMA3':['p50','p100'],
#             'STAT5A':['p30','p99'],
#             'ESRRG':['p50','p99'],
#             'FLT1':['p50','p100'],
#             'ENG':['p60','p100'],
#             'FOSL1':['p10','p99'],
#             'JUNB':['p30','p99'],
#             'FOS':['p30','p99'],
#             'HLA-G':['p30','p100'],
#             'PLAC8':['p50','p100'],
#         }
        
        
#         marker_genes = ["ERVFRD-1","PAPPA",'LAMA3','STAT5A','ESRRG','FLT1','ENG','FOSL1','JUNB','FOS','HLA-G','PLAC8']
        
#         marker_genes_sel = []
#         [ marker_genes_sel.append(i) for i in marker_genes if i in thresh.keys()]
        
#         plt.rcParams['figure.figsize'] =  (20,15)
#         plt.rcParams['legend.loc'] = 'upper left'
#         sc.set_figure_params(fontsize = 15, frameon = False, dpi = 100, )
#         res_p = sc.pl.umap(gene_matrix, 
#                    use_raw=False, 
#                    ncols=4,
#                    color= marker_genes_sel,
#                    vmin = [ thresh[i][0]  for i in marker_genes_sel],
#                    vmax= [ thresh[i][1]  for i in marker_genes_sel],
#                    size= 6,#15,
#                    color_map = my_cmap,
#                    #colorbar_loc = None, only in scanpy 1.9
#                    show = False,
#                    return_fig = False, #return as a whole fig for save?
#                    #legend_fontsize = 20,
#                    #legend_fontweight = 10,
#                    legend_loc = 'left margin' #unwork
#                    #legend_loc = ['left margin','right margin','right margin','right margin','right margin','right margin','right margin','right margin','right margin','none','right margin','right margin'] #not work, do not know why
#                   )
        
#         #res_p[0]._remove_legend()
#         #res_p[0].legend(loc='upper left')

#         #plt.colorbar(res_p[0].contourf,shrink=0.5)#fraction=0.046, pad=0.04)
#         #plt.colorbar(res_f,ax=res_p[0], shrink=0.5)#,pad=0.01, fraction=0.08, aspect=30)

#         plt.suptitle('umap tunning: mdist ' + str(mdist) + ' spread ' + str(spread) + ' random_state ' + str(random_state),y = 1.0, x = 0.45)

#         #plt.legend([res_p[0],res_p[1],res_p[2],res_p[3]],['test1','test2','test3','test4'])
        
#         ##save figure
#         #plt.show()
#         ##res_p.savefig()


###rotate umap by any degree###

del range

def rotateXY(x,y,degree):
    rotate_d = degree*(pi/180)
    xx = []
    yy = []
    for i in range(len(x)):
        xx.append(x[i]*cos(rotate_d)-y[i]*sin(rotate_d) )
        yy.append ( y[i]*cos(rotate_d)+x[i]*sin(rotate_d) )
    return xx,yy


(xx,yy) = rotateXY(x=umap[:,0],y=umap[:,1],degree=-15)

len(xx) == len(yy) #True

umap_rotate = []
for i in range(len(xx)):
    umap_rotate.append([xx[i],yy[i]])

umap_rotate = np.array(umap_rotate)



##fix and overwrite obj?
# combine.obsm['X_umap'].shape[0] == umap.shape[0]#True
# combined.obsm['X_umap'] = umap
# combined.obsm['X_umap_rotate'] = umap_rotate


# plt.rcParams['figure.figsize'] =  (20,15)
# sc.pl.umap(combined, color="leiden", interactive=False, width=600,use_rep='X_umap')
# sc.pl.umap(combined, color="cluster", interactive=False, width=600,use_rep='X_umap_rotate')



###tuning umap code ok##




#######tuning leiden cluster
combined.obs['leiden_bk'] = combined.obs['leiden']

import numpy as np
for i in np.arange(0.1,1,0.1):
#for i in np.arange(1,1.5,0.1):
#for i in np.arange(1.5,2.0,0.1):
    print('tuning cluster leiden with res: ' + str(i))
    sc.tl.leiden(combined, resolution = i)
    sc.pl.umap(combined, color=['leiden'], size=26,ncols=2,legend_loc='on data',frameon=True,legend_fontsize=10,title = "leiden " + str(i))


#fix leiden cluster
sc.tl.leiden(combined, resolution = 0.6)
sc.pl.umap(combined, color=['leiden'], size=26,ncols=2,legend_loc='on data',frameon=True,legend_fontsize=10)


combined.obs['leiden'].value_counts()
0     11371
1      8784
2      8091
3      6959
4      4444
5      4436
6      2374
7      1588
8       376
9       234
10       16


##tuning leiden ok


##encode_graph to get feature embedding (the integration embedding)
feature_embeddings = glue.encode_graph(graph)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5] #27768 x 50  #27488 x 50

##GLUE embedding

combined.obsm['X_glue_harmony'].shape #48673 x 50 #17665 x 50

X_glue_harmony = pd.DataFrame(combined.obsm['X_glue_harmony'],index = combined.obs.index, columns = ['GLUE'+str(x+1) for x in range(50) ] )
X_glue_harmony.to_csv('x_glue_harmony.txt',sep = '\t')#,index=True,columns=True)




##save/reload
with open('rna.final.pickle','wb') as fh:
    pickle.dump(rna,fh)
    
with open('atac.final.pickle','wb') as fh:
    pickle.dump(atac,fh)
    
# with open('graph.final.pickle','wb') as fh:
#     pickle.dump(graph,fh)
    
with open('combined.final.pickle','wb') as fh:
    pickle.dump(combined,fh)

with open('feature_embeddings.final.pickle','wb') as fh:
    pickle.dump(feature_embeddings,fh)


#X_map = pd.DataFrame( combined.obsm['X_umap'],columns=['UMAP1','UMAP2'])#,index=combined.var_names )
    



