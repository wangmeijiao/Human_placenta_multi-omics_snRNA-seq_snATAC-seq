
#https://kzhang.org/SnapATAC2/tutorials/pbmc.html


#!pip install snapatac2

# tensorflow 2.5.0 requires h5py~=3.1.0, but you have h5py 2.10.0 which is incompatible.
# tensorflow 2.5.0 requires typing-extensions~=3.7.4, but you have typing-extensions 4.4.0 which is incompatible.

# Successfully installed igraph-0.10.2 kaleido-0.2.1 plotly-5.11.0 polars-0.14.27 pooch-1.6.0 rustworkx-0.12.0 snapatac2-2.1.3 tenacity-8.1.0 tqdm-4.64.1 typing-extensions-4.4.0


#!pip install h5py==3.1.0
#tensorflow 2.5.0 requires typing-extensions~=3.7.4, but you have typing-extensions 4.4.0 which is incompatible.
#cellphonedb 3.0.0 requires h5py<3.0.0, but you have h5py 3.1.0 which is incompatible.
#Successfully installed h5py-3.1.0

#!pip install typing-extensions==3.7.4
#polars 0.14.27 requires typing_extensions>=4.0.0; python_version < "3.10", but you have typing-extensions 3.7.4 which is incompatible.
#Successfully installed typing-extensions-3.7.4

#!pip install tensorflow==2.10.0
#Successfully installed absl-py-1.3.0 flatbuffers-22.10.26 keras-2.10.0 libclang-14.0.6 numpy-1.23.4 tensorboard-2.10.1 tensorflow-2.10.0 tensorflow-estimator-2.10.0 tensorflow-io-gcs-filesystem-0.27.0

#!pip install polars
#Successfully installed typing-extensions-4.4.0

#!pip install scanpy==1.8.2

#!pip install anndata==0.8.0
#!pip install --upgrade magic-impute


import snapatac2 as snap
snap.__version__ #2.2.0  #2.2.0.dev0 #2.1.3

from snapatac2.tools._misc import aggregate_X


import os
os.environ['NUMEXPR_MAX_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'


##solve for "BLAS : Program is Terminated. Because you tried to allocate too many memory regions." error
os.environ['OPENBLAS_NUM_THREADS'] = '1' #important!!
#os.environ['GOTO_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

os.environ['MKL_DYNAMIC'] = 'FALSE'
os.environ['MKL_NUM_THREADS'] = '1'


#import scanpy #1.8.2 #1.6.1
import scanpy as sc #1.8.2

import pickle
import dill

import numpy as np #1.23.4

import polars as pl #0.14.27

import matplotlib
from matplotlib import pyplot as plt

import plotly.express as px
#import plotly
#plotly.plot()
import plotly.io as pio
pio.renderers.default = "png"

import plotly.graph_objects as go

from natsort import index_natsorted
import pandas as pd




from umap import UMAP


import color_utils

from yellowbrick.style.palettes import PALETTES, SEQUENCES, color_palette 


from IPython.display import Image

import scipy
from scipy.sparse import coo_matrix #to change csr to coo then mmwrite will ouput coordinate instead of array
from scipy.io import mmwrite
from scipy.sparse import save_npz
from scipy.stats import zscore

from math import log10
from math import ceil
from math import pi
from math import cos,sin
#from math import sqrt

import seaborn as sns


##show all available backends
plt.rcsetup.all_backends

plt.rcParams['backend'] #current backend
#'module://matplotlib_inline.backend_inline'

plt.get_backend()
#plt.use('Agg')



pl.Config.set_fmt_str_lengths(22)


#sc.logging.print_header()
#scanpy==1.8.2 anndata==0.8.0 umap==0.5.3 numpy==1.22.4 scipy==1.8.1 pandas==1.3.5 scikit-learn==1.1.1 statsmodels==0.13.2 python-igraph==0.10.2 louvain==0.7.1 leidenalg==0.7.0 pynndescent==0.5.7

#######check package version
pkg_list = ['snapatac2','scanpy','scrublet','pandas','numpy','numba','loompy','anndata','umap','polars']
def print_header(pkg_list):
    for pkg in pkg_list:
            try:
                imp = __import__(pkg)
                print (pkg + '==' + imp.__version__, end = ", ")
            except (ImportError, AttributeError):
                pass


print_header(pkg_list)

snapatac2==2.2.0, scanpy==1.8.2, pandas==1.3.5, numpy==1.22.4, numba==0.55.2, loompy==3.0.7, anndata==0.8.0, umap==0.5.3, polars==0.14.12
#snapatac2==2.2.0.dev0, scanpy==1.8.2, pandas==1.3.5, numpy==1.22.4, numba==0.55.2, loompy==3.0.7, anndata==0.8.0





###colors


#fh = open("/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/color_tfdev.txt",'r')

color_tfdev = []
with open("/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/color_tfdev.txt",'r') as fh:
    lines =  fh.readlines()
    for i in lines:
        color_tfdev.append(i.strip("\n"))

#fh.close()


my_cmap_tfdev = matplotlib.colors.LinearSegmentedColormap.from_list('Z score', color_tfdev)
my_cmap_tfdev




palcolor = ['lightgrey',
 '#FFF7F3',
 '#FDE0DD',
 '#FCC5C0',
 '#FA9FB5',
 '#F768A1',
 '#DD3497',
 '#AE017E',
 '#7A0177',
 '#49006A']


my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('expression',palcolor,N=10)
my_cmap


#######################predefined color #######################
color_snap ={'1':'grey','2':'#E31A1C','3':'#FFD700','4':'#771122','5':'#777711','6':'#1F78B4','7':'#68228B','8':'#AAAA44','9':'#60CC52','10':'#771155','11':'#DDDD77','12':'#774411','13':'#AA7744','14':'#AA4455','15':'#117744'}

color_snap_mod1 = {'1':'#777711','2':'#E31A1C','3':'#68228B','4':'#771122','5':'grey','6':'#1F78B4','7':'#FFD700','8':'#AAAA44','9':'#60CC52','10':'#771155','11':'#DDDD77','12':'#774411','13':'#AA7744','14':'#AA4455','15':'#117744'}


color_signac = {
'0':'#E6D55E','1':'#792B8A','2':'#DA703D','3':'#9DC8E5','4':'#BA273C','5':'#C2C184','6':'#7F8084','7':'#65AB53','8':'#D082AF','9':'#496EAB','10':'#DE896D','11':'#491F8B','12':'#E1AD49','13':'#8E1B85','14':'#E7EE77','15':'#7D1A1D','16':'#96B355'}


color_good = ["#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" ]


#####customized build colors#####

#######0. color list for categorical usage
#color_palette = adata.uns['louvain_colors']
#color_palette = adata.uns['leiden_colors']
color_palette_preset = color_good #list
color_palette_preset = [x for x in color_signac.values()] #dict to list
color_palette_preset = [x for x in color_snap_mod1.values()] #dict to list

#color_utils.hex_to_rgb_color_list(color_palette)

##preset pretty categorial colors in yellowbrick package
#https://www.scikit-yb.org/en/latest/api/palettes.html#color-palettes
color_palette_preset = SEQUENCES['Spectral'][5]
color_palette_preset = SEQUENCES['Spectral'][11]
color_palette_preset = SEQUENCES['RdYlBu'][11]
color_palette_preset = SEQUENCES['Reds'][9]
color_palette_preset = SEQUENCES['RdBu'][11]
#color_palette_preset.reverse()

#view 
color_palette(color_palette_preset).plot()



#########1 matplotlib predefined cmap color set ,color_map for continuous data
#https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
color_map = 'PuBu' ##good?
color_map = 'RdBu_r' #reverse of RdBu, use this?
#color_map = 'RdYlGn_r'
#color_map = 'PuOr_r' #good?
#color_map = 'Reds'
#color_map = 'viridis' #scanpy default colors
#color_map = 'RdGy_r' #good?

color_utils.draw_cmap(color_map) #default 255 colors?

#get n hex color code from rgb of cmap name
n = 20
cmap_palette = plt.cm.get_cmap(color_map, n) #extend to n colors

#get color code
cmap_palette = plt.cm.get_cmap(color_map)(np.arange(n)).tolist()
#[print( '#%02x%02x%02x' %  (int(x*255), int(y*255),int(z*255) ) )  for x,y,z,f in color_map_palette]
#the same as followings
for x,y,z,f in cmap_palette:
    print( '#%02x%02x%02x' %  (int(x*255), int(y*255),int(z*255) ) )

#view
color_utils.draw_cmap(cmap_palette) #'RdBu_r', 20 colors predefined matplotlib cmap


######2 seaborn self-present color bar######
#https://seaborn.pydata.org/tutorial/color_palettes.html

color_set = 'vlag'
color_set = 'Spectral_r'
color_set = 'coolwarm'

cmap_sns_palette = sns.color_palette(color_set, as_cmap=True) #need sns 0.11

#get color code
cmap_sns_palette = sns.color_palette(color_set, as_cmap=False)
for x,y,z in color_sns_palette:
    print( '#%02x%02x%02x' %  (int(x*255), int(y*255),int(z*255) ) )

#view
color_utils.draw_cmap(cmap_sns_palette ) #seaborn palettes object

######3 customized color bar (of unbalanced) with color_utils#######

#https://github.com/stefmolin/Custom-Colormaps

#rgbs = [[0, 0, 1], [0,0,1],[1, 1, 1],[1, 0, 0]]

#rgbs = [[45,4,74],[255,255,255],[102,3,32]]
rgbs = [[45,4,74],[255,255,255],[198,70,68],[111,6,34]]
rgbs = [ [i/255  for i in x  ]  for x in rgbs] #255 to ratio
cmap_my_palette = color_utils.blended_cmap(rgbs) #use np.linspace to get equal spaced numbers for r, g, b channels

##blend from preset categorial color

rgbs = color_utils.hex_to_rgb_color_list(color_palette_preset)  #build from yellow brick preset colors
cmap_my_palette = color_utils.blended_cmap(rgbs)

#view
color_utils.draw_cmap(cmap_my_palette ) #a customized Colorbar object of color_utilts

##########




###create anndata for each library

 #for i in {1..6};do echo -ne "(\"placenta_10X_early$i\",\"../../placenta_10X_early${i}/02.snapatac2/PLA-early${i}-atac.h5ad\"),\n"; done

#filtered by tss score and fragments
h5ad_files = [
             ("placenta_10X_early1","../../placenta_10X_early1/02.snapatac2/PLA-early1-atac.h5ad"),
            ("placenta_10X_early2","../../placenta_10X_early2/02.snapatac2/PLA-early2-atac.h5ad"),
            ("placenta_10X_early3","../../placenta_10X_early3/02.snapatac2/PLA-early3-atac.h5ad"),
            ("placenta_10X_early4","../../placenta_10X_early4/02.snapatac2/PLA-early4-atac.h5ad"),
            ("placenta_10X_early5","../../placenta_10X_early5/02.snapatac2/PLA-early5-atac.h5ad"),
            ("placenta_10X_early6","../../placenta_10X_early6/02.snapatac2/PLA-early6-atac.h5ad"),


             ]
#name, fl = files[0]


frag_files = [
             ("placenta_10X_early1","../../placenta_10X_early1/01.data_cellranger_atac/fragments.tsv.gz"),
            ("placenta_10X_early2","../../placenta_10X_early2/01.data_cellranger_atac/fragments.tsv.gz"),
            ("placenta_10X_early3","../../placenta_10X_early3/01.data_cellranger_atac/fragments.tsv.gz"),
            ("placenta_10X_early4","../../placenta_10X_early4/01.data_cellranger_atac/fragments.tsv.gz"),
            ("placenta_10X_early5","../../placenta_10X_early5/01.data_cellranger_atac/fragments.tsv.gz"),
            ("placenta_10X_early6","../../placenta_10X_early6/01.data_cellranger_atac/fragments.tsv.gz"),


             ]


gmat_files = {
            "placenta_10X_early1":"../../placenta_10X_early1/02.snapatac2/gene_matrix.h5ad",
            "placenta_10X_early2":"../../placenta_10X_early2/02.snapatac2/gene_matrix.h5ad",
            "placenta_10X_early3":"../../placenta_10X_early3/02.snapatac2/gene_matrix.h5ad",
            "placenta_10X_early4":"../../placenta_10X_early4/02.snapatac2/gene_matrix.h5ad",
            "placenta_10X_early5":"../../placenta_10X_early5/02.snapatac2/gene_matrix.h5ad",
            "placenta_10X_early6":"../../placenta_10X_early6/02.snapatac2/gene_matrix.h5ad",


           }



sex_list = {
  # 'PLA-early1-RNA' : 'male',
  # 'PLA-early2-RNA' :  'female',
  # 'PLA-early3-RNA' :  'female',
  # 'PLA-early4-RNA' :  'male',
  # 'PLA-early5-RNA' :  'male',
  # 'PLA-early6-RNA' :  'female'
  'PLA-atac-early1' : 'male',
  'PLA-atac-early2' :  'female',
  'PLA-atac-early3' :  'female',
  'PLA-atac-early4' :  'male',
  'PLA-atac-early5' :  'male',
  'PLA-atac-early6' :  'female'
    
}



###################Creating AnnDataSet object##################


datalist = {}

for  id, file in h5ad_files:
    print(id + ": " + file)
    datalist[id] = snap.read(file)
    datalist[id]=datalist[id].to_memory() #in case not to affect the original file!!
    

datalist.keys()

#'placenta_10X_early1', 'placenta_10X_early2', 'placenta_10X_early3', 'placenta_10X_early4', 'placenta_10X_early5', 'placenta_10X_early6'

datalist['placenta_10X_early1'] # 5334 x 6176550
datalist['placenta_10X_early2'] # 7788 x 6176550
datalist['placenta_10X_early3'] # 5876 x 6176550
datalist['placenta_10X_early4'] # 5596 × 6176550 (diff 2 barcodes because of scrublet ) # 5891 x 6176550
datalist['placenta_10X_early5'] # 3237 x 6176550
datalist['placenta_10X_early6'] # 5968 x 6176550


datalist['placenta_10X_early1'].obs_names # 5334
datalist['placenta_10X_early2'].obs_names # 7788
datalist['placenta_10X_early3'].obs_names # 5876
datalist['placenta_10X_early4'].obs_names # 5596 # 5891
datalist['placenta_10X_early5'].obs_names # 3237
datalist['placenta_10X_early6'].obs_names # 5968



# datalist['placenta_10X_early4'].obs["is_doublet"].value_counts()
# False    5583
# True      308

#datalist['placenta_10X_early4'] = datalist['placenta_10X_early4'][~datalist['placenta_10X_early4'].obs["is_doublet"],:]

datalist['placenta_10X_early4'].obs["is_doublet"].value_counts()
False    5596 #after call doublet and rewrite



###quick look each obj for umap and marker gene

for i in datalist.keys():
    print(i)
    
    data = datalist[i]
    gene_matrix = sc.read(gmat_files[i])
    
    #snap.pl.umap(data, color="leiden", interactive=False, width=500)
    
    marker_genes = ['DNMT1','MKI67','ERVFRD-1','PSG8','CGA','PAPPA','CSHL1','FLT1','ENG','DDX60','HLA-G','VIM','PECAM1','CD68','HBZ'] #marker gene
    # = ['DNMT1','ERVFRD-1','PSG8','CGA','SH3TC2','PAPPA','CSHL1','FLT1','ENG','HLA-G','XIST','DDX3Y','RPS4Y1','USP9Y','PCDH11Y'] #sex gene
    #marker_genes = ['DNMT1','ERVFRD-1','PSG8','CGA','SH3TC2','PAPPA','CSHL1','FLT1','ENG','HLA-G','DDX60','DDX58','CDKN1A','MAP4K4','CROT'] #apoptosis gene?

    
    marker_genes_sel = [ i for i in marker_genes if i in gene_matrix.var_names ]
    
    plt.rcParams['figure.figsize'] =  (20,15)
    #plt.rcParams['font.size'] = 45
    sc.set_figure_params(fontsize = 15, frameon = False, dpi = 100, )
    sc.pl.umap(gene_matrix, use_raw=False, 
               color=["leiden"] + marker_genes_sel ,
               vmin='p30',
               vmax='p99',
               size=15,
               color_map = my_cmap,
               show = False# False
            )
    plt.suptitle(i,y = 0.9, x = 0.45)
    
    #sc.pl.violin(gene_matrix, marker_genes, use_raw=False, groupby='leiden')



###simplize the h5ad file and rewrite to new dir (for AnnDataSet need h5ad with identical attrubution)

for i in datalist.keys():
    print(i)
    
    data = datalist[i]
    
    #modify cellid
    data.obs['library'] = pl.Series('library', [i] * data.shape[0] , dtype = 'string')
    cellid = np.array(data.obs['library']) + ":" + np.array(data.obs_names)
    data.obs_names = cellid
    
    
    #data.uns = {}
    #datalist[i].var = []
    #data.obsm = {}
    #datalist[i].obsp = {}
    data.write(i+".simple.h5ad")
    #datalist[i].close()
    
    
    

    
# datalist['placenta_10X_early4'].obs.dtypes
# datalist['placenta_10X_early5'].obs.dtypes

# # datalist['placenta_10X_early4'].obs['leiden']
# # datalist['placenta_10X_early5'].obs['leiden']
# # datalist['placenta_10X_early3'].obs['leiden']

# datalist['placenta_10X_early4'].obs['is_doublet'].value_counts() #here is the problem!
# datalist['placenta_10X_early5'].obs['is_doublet'].value_counts()




h5ad_files_new = [
             ("PLA-atac-early1","placenta_10X_early1.simple.h5ad"),
            ("PLA-atac-early2","placenta_10X_early2.simple.h5ad"),
            ("PLA-atac-early3","placenta_10X_early3.simple.h5ad"),
            ("PLA-atac-early4","placenta_10X_early4.simple.h5ad"),
            ("PLA-atac-early5","placenta_10X_early5.simple.h5ad"),
            ("PLA-atac-early6","placenta_10X_early6.simple.h5ad"),


             ]
    
    
    

##reload new h5ad files
datalist = {}

for  id, file in h5ad_files_new:
    print(id + ": " + file)
    datalist[id] = snap.read(file)
    datalist[id]=datalist[id].to_memory() #in case not to affect the original file!!
    

datalist.keys()


#'placenta_10X_early1', 'placenta_10X_early2', 'placenta_10X_early3', 'placenta_10X_early4', 'placenta_10X_early5', 'placenta_10X_early6'

datalist['PLA-atac-early1'] # 5334 x 6176550
datalist['PLA-atac-early2'] # 7788 x 6176550
datalist['PLA-atac-early3'] # 5876 x 6176550
datalist['PLA-atac-early4'] # 5583 × 6176550 #5891 x 6176550
datalist['PLA-atac-early5'] # 3237 x 6176550
datalist['PLA-atac-early6'] # 5968 x 6176550


datalist['PLA-atac-early1'].obs_names # 5334
datalist['PLA-atac-early2'].obs_names # 7788
datalist['PLA-atac-early3'].obs_names # 5876
datalist['PLA-atac-early4'].obs_names # 5583 # 5891
datalist['PLA-atac-early5'].obs_names # 3237
datalist['PLA-atac-early6'].obs_names # 5968

datalist['PLA-atac-early1'].var_names #6176550, 500bp per bin
datalist['PLA-atac-early2'].var_names #6176550, 500bp per bin
datalist['PLA-atac-early3'].var_names #6176550, 500bp per bin
datalist['PLA-atac-early4'].var_names #6176550, 500bp per bin
datalist['PLA-atac-early5'].var_names #6176550, 500bp per bin
datalist['PLA-atac-early6'].var_names #6176550, 500bp per bin

all(datalist['PLA-atac-early1'].var_names == datalist['PLA-atac-early2'].var_names) #True
all(datalist['PLA-atac-early2'].var_names == datalist['PLA-atac-early3'].var_names) #True
all(datalist['PLA-atac-early3'].var_names == datalist['PLA-atac-early4'].var_names) #True
all(datalist['PLA-atac-early4'].var_names == datalist['PLA-atac-early5'].var_names) #True
all(datalist['PLA-atac-early5'].var_names == datalist['PLA-atac-early6'].var_names) #True


datalist['PLA-atac-early1'].var['selected'].value_counts()
False    5097643
True     1078907

datalist['PLA-atac-early2'].var['selected'].value_counts()
False    5084176
True     1092374

datalist['PLA-atac-early3'].var['selected'].value_counts()
False    5122851
True     1053699

datalist['PLA-atac-early4'].var['selected'].value_counts()
False    5068220
True     1108330

datalist['PLA-atac-early5'].var['selected'].value_counts()
False    5152154
True     1024396

datalist['PLA-atac-early6'].var['selected'].value_counts()
False    5137926
True     1038624


###################################start to merge all h5ad obj################################



dataset.close()


dataset = snap.AnnDataSet(adatas=h5ad_files_new, filename="PLA_atac_early_combined.h5ad")#,add_key='test') #quick in fact
#data = snap.create_dataset(h5ad_files, "PLA_early_combined.h5ads")

dataset

AnnDataSet object with n_obs x n_vars = 33786 x 6176550 backed at 'PLA_atac_early_combined.h5ad'
contains 6 AnnData objects with keys: 'PLA-atac-early1', 'PLA-atac-early2', 'PLA-atac-early3', 'PLA-atac-early4', 'PLA-atac-early5', 'PLA-atac-early6'
    obs: 'sample'
    var: 'selected'
    uns: 'AnnDataSet'
    
    
# AnnDataSet object with n_obs x n_vars = 56877 x 6176550 backed at 'PLA_early_combined.h5ads'
# contains 6 AnnData objects with keys: 'PLA_atac_early1', 'PLA_atac_early2', 'PLA_atac_early3', 'PLA_atac_early4', 'PLA_atac_early5', 'PLA_atac_early6'
#     obs: 'sample'
#     var: 'selected'
#     uns: 'AnnDataSet'


    
##an in-memory anndata (https://kzhang.org/epigenomics-analysis/anndata.html)
dataset = dataset.to_adata() #merge as one obj, already in mem


##kept in mem only

#dataset = dataset.to_memory()


##merge and add obs and inserstions data to merged obj


##check insertion ccr matrix column length aligned
for  id, file in h5ad_files_new:
    print(id + ": " + file)

    print(datalist[id].uns['reference_sequences'].reference_seq_length.sum()) #3088269832
    print("shape is:" + str(datalist[id].obsm['insertion'].shape)) #5968, 3088269832, bp level insertion site coordinations, csr_matrix
    print('reference chr is ' + datalist[id].uns['reference_sequences'].reference_seq_name)
    print('reference len is ' + str(datalist[id].uns['reference_sequences'].reference_seq_length))

#insertion mat column max length is the genome length (snapatac2 predefined genome version)

    
    

##check  obs row name and colname aligned (pandas df)


import pandas as pd

obs_cat = pd.DataFrame()

for  id, file in h5ad_files_new:
    print(id + ": " + file)
    obs_cat = pd.concat([obs_cat,datalist[id].obs], axis = 0)
    

all(obs_cat.index == dataset.obs.index) #True

dataset.obs = pd.concat( [dataset.obs, obs_cat], axis=1  )



#start to merge insertion slot
from scipy.sparse import hstack,vstack


# insertions1 = datalist['PLA-atac-early1'].obsm['insertion']
# insertions2 = datalist['PLA-atac-early2'].obsm['insertion']


insertions_merge = vstack(
                          (datalist['PLA-atac-early1'].obsm['insertion'],
                           datalist['PLA-atac-early2'].obsm['insertion'],
                           datalist['PLA-atac-early3'].obsm['insertion'],
                           datalist['PLA-atac-early4'].obsm['insertion'],
                           datalist['PLA-atac-early5'].obsm['insertion'],
                           datalist['PLA-atac-early6'].obsm['insertion']
                          )


)

dataset.obsm['insertion'] = insertions_merge



dataset.uns['reference_sequences'] = datalist['PLA-atac-early3'].uns['reference_sequences']



snap.pl.tsse(dataset, interactive=False)
#snap.pp.filter_cells(dataset, min_counts=5000, min_tsse=5, max_counts=50000)



####start to standard batch effect removal, dimension reduction clustering 


data = dataset

del dataset

#data = data.to_memory()


###restart from here after reload h5ad file

data.obs_names

#data.obs['sample'].value_counts()

data.obs['sample'].value_counts().sort_index() #series sort by rowname

PLA-atac-early1    5334
PLA-atac-early2    7788
PLA-atac-early3    5876
PLA-atac-early4    5583
PLA-atac-early5    3237
PLA-atac-early6    5968


data.obs['sex'] = data.obs['sample'].map(sex_list)

data.obs['sex'].value_counts().sort_index()




##write bin.bed
data.var['selected'].value_counts()

#before filter_raw_cluster 
false 5097643
true  1078907

##filtered_raw_cluster
false 5164772
true 1011778

data.var_names #6176550, 500bp per bin
data.X.shape #(33786, 6176550)

with open('bin.bed','w') as fh:
    for i in data.var_names:
        (chr,range) = i.split(":")
        (start, end) = range.split("-")
        fh.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + i + "\n")


        
##write bin.selected.bed
index_sel = [i for i, val in enumerate(data.var['selected'].to_list()) if val == True]        

data.var['selected'][index_sel].value_counts() #all 1078907 True



bin_df = pd.DataFrame(data.var_names,columns=['name']) #6176550 x 1

bin_df_sel = bin_df.iloc[index_sel,:] #1078907 x 1 


with open('bin.select.bed','w') as fh:
    for i in bin_df_sel['name']:
        (chr,range) = i.split(":")
        (start, end) = range.split("-")
        fh.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + i + "\n")

del range
        
#######




#####reselect after qc filtering + raw cluster filtering

data
AnnData object with n_obs × n_vars = 27380 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'selected'
    uns: 'AnnDataSet', 'reference_sequences', 'spectral_eigenvalue'
    obsm: 'insertion', 'X_umap', 'X_spectral', 'X_spectral_mnn', 'X_spectral_harmony'
    obsp: 'distances'
    

snap.pl.tsse(data, interactive=False)
#snap.pp.filter_cells(dataset, min_counts=5000, min_tsse=5, max_counts=50000)

snap.pp.select_features(data,blacklist="/sda/mjwang/pwdex/placenta_10X_combine/03.snRNA_snATAC/scenicplus/data/hg38-blacklist.v2.bed") #must recompute after merge (and filtering?)

##selection method: simple zscore method
# mean = count[selected_features].mean()
# std = math.sqrt(count[selected_features].var())
# zscores = np.absolute((count - mean) / std)
# cutoff = np.sort(zscores)[most_variable - 1]
# selected_features &= zscores <= cutoff

data.var['selected'] = True

data.var['selected'].value_counts()

##filter_raw_cluster_c3_c5 (not re-select_feature )
False    5097643
True     1078907


False    5164772
True     1011778


#np.set_printoptions(threshold = 2000)
test_mat = data.X[5000:7000,5000:7000].todense() #no scale and log in ATAC bin-matrit??

test_mat[~np.all(test_mat == 0, axis = 1)] #view non-zero lines in np matrix 


snap.tl.spectral(data)#, sample_size = 10000) #mem big even sample down, will use selected bin with jaccard distance by default, subset most var bin data.X[:,data.var['selected']]
#snap.tl.spectral(data)

snap.pl.spectral_eigenvalues(data, interactive=False)

snap.tl.umap(data, use_dims = 15) #after filtering raw clusters
#snap.tl.umap(data, use_dims = 18)

snap.pl.umap(data, color="sample", interactive=False, width = 800)



###Batch correction###
#apply two different approaches, Harmony and modified MNNCorrect, to remove donor specific differences


snap.pp.mnc_correct(data, "sample", use_dims = 18)#12)
##snap.pp.harmony(data, "sample", use_dims = 18)#, max_iter_harmony = 20) #error, use standalone method
#harmony Stopped before convergence (increase max_iter_harmony )

import harmonypy
harmony_out = harmonypy.run_harmony(data.obsm['X_spectral'][:,range(18)], data.obs, 'sample', max_iter_harmony = 20)

#2023-01-13 22:24:05 - INFO - Converged after 11 iterations

#Converged after 10 iterations


data.obsm["X_spectral_harmony"] = harmony_out.Z_corr.T
33786 x 18



##do umap embedding
snap.tl.umap(data, use_rep="X_spectral_mnn")
snap.pl.umap(data, color="sample", interactive=False, width=800)


snap.tl.umap(data, use_rep="X_spectral_harmony")
snap.pl.umap(data, color="sample", interactive=False, width=800)



######filter data for raw cluster c3 and c5 ######

data.obs['leiden'].value_counts()
0     6541
1     5761
2     5664
3     4065
4     3024
5     2341
6     2230
7     1026
8      964
9      763
10     526
11     510
12     371


cellid = data.obs_names.to_list() #33786
leiden_sel = data.obs['leiden'].iloc[( (data.obs['leiden'] != '3') & (data.obs['leiden'] != '5') ).to_list()]
#27380
#29721

leiden_sel.value_counts()
0     6541
1     5761
2     5664
4     3024
6     2230
7     1026
8      964
9      763
10     526
11     510
12     371
3        0
5        0

cellid_sel = leiden_sel.index.to_list() #27380

data_new = data[cellid_sel,:]

data_new
View of AnnData object with n_obs × n_vars = 27380 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'selected'
    uns: 'AnnDataSet', 'reference_sequences', 'spectral_eigenvalue'
    obsm: 'insertion', 'X_umap', 'X_spectral', 'X_spectral_mnn', 'X_spectral_harmony'
    obsp: 'distances'
    
View of AnnData object with n_obs × n_vars = 29721 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'selected'
    uns: 'AnnDataSet', 'reference_sequences', 'spectral_eigenvalue'
    obsm: 'insertion', 'X_spectral_harmony', 'X_spectral', 'X_spectral_mnn', 'X_umap'
    obsp: 'distances'
    
data_new.obs['leiden'].value_counts()

0     6541
1     5761
2     5664
4     3024
6     2230
7     1026
8      964
9      763
10     526
11     510
12     371

0     6541
1     5761
2     5664
4     3024
5     2341
6     2230
7     1026
8      964
9      763
10     526
11     510
12     371


snap.pl.umap(data_new, color="leiden", interactive=False, width=500)



data = data_new

del data_new

snap.pl.umap(data, color="leiden", interactive=False, width=500)



###umap tunning with parameters: min_dist (default 0.1), spread (default 1.0), random_state (default 0)#######
dim_use = data.obsm['X_spectral_harmony'] #29721 x 18 #11114 x 18
groups = data.obs['leiden']

random_state = 0 #123
n_comps = 2
#min_dist = 0.1 
#spread = 1.0

config = dict({'staticPlot': True})

###choose:
mdist_sel = 0.3
spread_sel = 1.0
#random_state = 0
#n_comps = 2

#for mdist in [0.1,0.2,0.3,0.4,0.5]:
for mdist in [mdist_sel]:
#    for spread in [0.5,1.0,1.5,2.0]:
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

        ##quick annotation within tunning
        if data.obsm['X_umap'].shape[0] == umap.shape[0]:#True
            gene_matrix.obsm["X_umap"] = umap
            gene_matrix.obs['leiden'] = data.obs['leiden']
            
        else:
            print('umap length diff')
        
        thresh = {
            'ERVFRD-1':['p30','p100'],
            'PAPPA':['p65','p100'],
            'LAMA3':['p50','p100'],
            'STAT5A':['p30','p99'],
            'ESRRG':['p50','p99'],
            'FLT1':['p50','p100'],
            'ENG':['p60','p100'],
            'FOSL1':['p10','p99'],
            'JUNB':['p30','p99'],
            'FOS':['p30','p99'],
            'HLA-G':['p30','p100'],
            'PLAC8':['p50','p100'],
        }
        
        
        marker_genes = ["ERVFRD-1","PAPPA",'LAMA3','STAT5A','ESRRG','FLT1','ENG','FOSL1','JUNB','FOS','HLA-G','PLAC8']
        
        marker_genes_sel = []
        [ marker_genes_sel.append(i) for i in marker_genes if i in thresh.keys()]
        
        plt.rcParams['figure.figsize'] =  (20,15)
        plt.rcParams['legend.loc'] = 'upper left'
        sc.set_figure_params(fontsize = 15, frameon = False, dpi = 100, )
        res_p = sc.pl.umap(gene_matrix, 
                   use_raw=False, 
                   ncols=4,
                   color= marker_genes_sel,
                   vmin = [ thresh[i][0]  for i in marker_genes_sel],
                   vmax= [ thresh[i][1]  for i in marker_genes_sel],
                   size= 6,#15,
                   color_map = my_cmap,
                   #colorbar_loc = None, only in scanpy 1.9
                   show = False,
                   return_fig = False, #return as a whole fig for save?
                   #legend_fontsize = 20,
                   #legend_fontweight = 10,
                   legend_loc = 'left margin' #unwork
                   #legend_loc = ['left margin','right margin','right margin','right margin','right margin','right margin','right margin','right margin','right margin','none','right margin','right margin'] #not work, do not know why
                  )
        
        #res_p[0]._remove_legend()
        #res_p[0].legend(loc='upper left')

        #plt.colorbar(res_p[0].contourf,shrink=0.5)#fraction=0.046, pad=0.04)
        #plt.colorbar(res_f,ax=res_p[0], shrink=0.5)#,pad=0.01, fraction=0.08, aspect=30)

        plt.suptitle('umap tunning: mdist ' + str(mdist) + ' spread ' + str(spread) + ' random_state ' + str(random_state),y = 1.0, x = 0.45)

        #plt.legend([res_p[0],res_p[1],res_p[2],res_p[3]],['test1','test2','test3','test4'])
        
        ##save figure
        #plt.show()
        ##res_p.savefig()

#         sc.pl.umap(gene_matrix, use_raw=False, color=["ERVFRD-1"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
        
#         sc.pl.umap(gene_matrix, use_raw=False, color=["PAPPA"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
#         sc.pl.umap(gene_matrix, use_raw=False, color=["LAMA3"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
#         sc.pl.umap(gene_matrix, use_raw=False, color=["STAT5A"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
#         sc.pl.umap(gene_matrix, use_raw=False, color=["ESRRG"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
        
#         sc.pl.umap(gene_matrix, use_raw=False, color=["FLT1"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
#         sc.pl.umap(gene_matrix, use_raw=False, color=["ENG"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
#         sc.pl.umap(gene_matrix, use_raw=False, color=["FOSL1"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
#         sc.pl.umap(gene_matrix, use_raw=False, color=["JUNB"],vmin='p30',vmax='p99',size=10,color_map = my_cmap)
#         sc.pl.umap(gene_matrix, use_raw=False, color=["FOS"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
        
#         sc.pl.umap(gene_matrix, use_raw=False, color=["HLA-G"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
#         sc.pl.umap(gene_matrix, use_raw=False, color=["PLAC8"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
        

####tuning umap code stop####



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
data.obsm['X_umap'].shape[0] == umap.shape[0]#True

data.obsm['X_umap'] = umap

data.obsm['X_umap_rotate'] = umap_rotate


plt.rcParams['figure.figsize'] =  (20,15)
snap.pl.umap(data, color="leiden", interactive=False, width=600,use_rep='X_umap')
snap.pl.umap(data, color="cluster", interactive=False, width=600,use_rep='X_umap_rotate')



#########Clustering (use harmony corrected dimension matrix)##############

snap.pp.knn(data, use_rep="X_spectral_harmony", use_dims = 18)

data.obs['leiden_bk'] = data.obs['leiden']

snap.tl.leiden(data, resolution = 1.2)
snap.tl.leiden(data, resolution = 0.9) #need knn run first to generate obsp['distances']
#snap.tl.leiden(data, resolution = 1)
#snap.tl.leiden(data, resolution = 1.2)
#snap.tl.leiden(data, resolution = 1.6)

snap.pl.umap(data, color="leiden", interactive=False, width=500)
#sc.pl.umap(data, use_raw=False, color=["leiden"],size=10, legend_loc = 'on data')

res_px = snap.pl.umap(data, color="leiden", interactive=False, marker_size=2,width=500,show=False,out_file=None)
#res_px.layout(yaxis = list( autorange="reversed") )
res_px.update_yaxes(autorange="reversed")
res_px.show(interactive=False)


#####iteratively tuning cluster resolution
#config = dict({'scrollZoom': False,'responsive': False,'staticPlot': True})
config = dict({'staticPlot': False})


##use resolution x for leiden

for i in np.arange(0.1,2,0.1):
#for i in [1.4,1.5,1.6,1.7]:
#for i in [0.8]: #use this
    snap.tl.leiden(data, resolution = i)
    #res_px = snap.pl.umap(data, color="leiden", interactive=False, marker_size=2,width=500,show=False,out_file=None)
    #res_px.show(interactive=False)
    plot_umap(adata=data, res= i)
    #Image(plot_umap(adata=data, res= i))

    
#fix leiden cluster
snap.tl.leiden(data, resolution = 1.4)
#snap.pl.umap(data, color="leiden", interactive=False, width=500)
plot_umap(adata=data, res= 1.4)

data.obs['leiden'].value_counts()
0     5698
1     5442
2     4545
3     4212
4     2663
5     1013
6      997
7      728
8      719
9      534
10     470
11     359

# 0     5725
# 1     5656
# 2     5070
# 3     3930
# 4     2609
# 5     2564
# 6     1023
# 7      976
# 8      760
# 9      526
# 10     510
# 11     372


    

def plot_umap(adata,
              marker_size: float = 2,
              color_use: list = color_good,
              res: float|int = 1.0
             ):
    #code modified from snap.pl.umap
    embedding = adata.obsm['X_umap'] #array
    groups = adata.obs['leiden']
    df = pd.DataFrame({
        "UMAP-1": embedding[:, 0],
        "UMAP-2": embedding[:, 1],
        'leiden': groups,
    })
    fig = px.scatter(
        df, x="UMAP-1", y="UMAP-2", color='leiden',
        color_discrete_sequence=color_use,#px.colors.qualitative.Dark24,
        title = 'resolution ' + str(res),
        width = 550,
        height = 450,
        category_orders = {'leiden':groups.value_counts().index.to_list(),
                           #'cluster':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
                          }

    )
    
    fig.update_traces(
      marker_size=marker_size,
      marker={"opacity": 1},
      #legendgrouptitle_font_color='red', 
      #legendgrouptitle_text= 'cluster',
      #selector=dict(type='scatter')
    )

    fig.update_layout(
      template="simple_white",
      legend= {'itemsizing': 'constant'},
      #title_text='test' #ok!
      margin=dict(l=20, r=20, t=60, b=20),
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

    fig.show(interactive=False,config=config,width = 550,height = 450) #will cause very large mem in chrome!

    #img_bytes = fig.to_image(format="png", width=300, height=350, scale=2)
    #Image(img_bytes)
    
    #fig.write_image("pdfs/leiden.pdf")
    #return img_bytes





##annotation with gene activity score

gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.hg38) #will copy leiden

gene_matrix = snap.pp.make_gene_matrix(data_cstb, snap.genome.hg38) 


#var_names varies of each calling!! but values eqaul
idx = gene_matrix.var_names.get_indexer(  gene_matrix_1.var_names )
(gene_matrix_1.X != gene_matrix.X[:,idx]).sum() #0



##let promoters = Promoters::new(transcripts, 2000, 0, true); use gene as count level by default (the overlap range of all transcripts of this gene, plus 2000 upstream, but use insertion site <start and end position of one fragment bed record>)
##snap.genome.hg38 == snap.genome.GRCh38
##better than snapATAC v1 (use bmat with gene body only with scale by RPM and smooth by magic)

#/home/mjwang/.cache/snapatac2/gencode_v41_GRCh38.gff3.gz



sc.pp.filter_genes(gene_matrix, min_cells= 5)
sc.pp.normalize_total(gene_matrix)
sc.pp.log1p(gene_matrix)



#Imputation (MAGIC to perform imputation and data smoothing)
#pip install magic-impute (use this)
#tmtoolkit 0.11.2 requires pandas>=1.4.0, but you have pandas 1.3.5 which is incompatible.
#Successfully installed Deprecated-1.2.13 graphtools-1.5.2 magic-impute-3.0.0 pandas-1.3.5 pygsp-0.5.1 scprep-1.2.1 tasklogger-1.2.0


sc.external.pp.magic(gene_matrix, solver="approximate")

# Copy umap embedding

gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]
gene_matrix.obs["leiden"] = data.obs["leiden"]

gene_matrix.obsm["X_umap"] = data_cstb.obsm["X_umap_rotate"]
gene_matrix.obs["leiden"] = data_cstb.obs["leiden"]
gene_matrix.obs["cluster"] = data_cstb.obs["cluster"]


###save and reload gene_matrix
gene_matrix


AnnData object with n_obs × n_vars = 22785 × 59264 
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
    var: 'n_cells'
    uns: 'log1p'
    obsm: 'X_umap'

AnnData object with n_obs × n_vars = 22815 × 59265 (cstb)
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden'
    var: 'n_cells'
    uns: 'log1p'
    obsm: 'X_umap'

AnnData object with n_obs × n_vars = 27380 × 59365 (filtered)
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'n_cells'
    uns: 'log1p'
    obsm: 'X_umap'
    
AnnData object with n_obs × n_vars = 29721 × 59393 (final)
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'n_cells'
    uns: 'log1p'
    obsm: 'X_umap'
    
AnnData object with n_obs × n_vars = 33786 × 59459
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'n_cells'
    uns: 'log1p', 'leiden_colors'
    obsm: 'X_umap'

    
#gene_matrix.write( filename="gene_matrix_cstb.h5ad")
gene_matrix = sc.read("gene_matrix_cstb.h5ad")

#gene_matrix.write( filename="gene_matrix.h5ad")
gene_matrix = sc.read("gene_matrix.h5ad")

#gene_matrix.write( filename="gene_matrix.filter_raw_cluster.h5ad")
gene_matrix = sc.read("gene_matrix.filter_raw_cluster.h5ad")

# AnnData object with n_obs × n_vars = 29721 × 59393
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
#     var: 'n_cells'
#     uns: 'log1p'
#     obsm: 'X_umap'
 
####save gmat raw or processed to mm format for r##

gene_matrix.X.shape #(29721, 59393) ,numpy array
gene_matrix.raw.X.shape #(29721, 59393) sparse matrix of numpy.float32,Compressed Sparse Row format


##Tn5 gmat raw count (for test reproducibility only, see atac_activity_Tn5/work_get_gmat_Tn5.py)
mmwrite('toRDS/gmat.raw.mtx.1',gene_matrix.X,field = 'integer') 

##save gmat colname and rowname
with open('toRDS/gmat.raw.rowname.txt.1','w') as fh:
    for i in gene_matrix.obs_names.to_list():
        fh.write(i + "\n")

with open('toRDS/gmat.raw.colname.txt.1','w') as fh:
    for i in gene_matrix.var_names.to_list():
        fh.write(i + "\n")
    
    
##Tn5 insertion count with filter, scale, magic smoothing    
mmwrite('toRDS/gmat.mtx',coo_matrix(gene_matrix.X),precision=3) #37G file size! read in R with Matrix::readMM() and save as sparse matrix rds will reduce file size

#save_npz(file = 'gmat.npz',matrix=gene_matrix.X)

##save gmat colname and rowname
with open('toRDS/gmat.rowname.txt','w') as fh:
    for i in gene_matrix.obs_names.to_list():
        fh.write(i + "\n")

with open('toRDS/gmat.colname.txt','w') as fh:
    for i in gene_matrix.var_names.to_list():
        fh.write(i + "\n")
    

# ##save rds by pyreadr? no, too big and too slow

# #https://github.com/ofajardo/pyreadr#basic-usage--writing-files
# !pip install pyreadr   
# import pyreadr #0.4.7

# ##prepare pd DataFrame to write (will big!)

# gmat_df = pd.DataFrame(gene_matrix.X,index=gene_matrix.obs_names, columns = gene_matrix.var_names)
# #22815 X 59265
# pyreadr.write_rds('gmat_cstb.rds',gmat_df, compress = 'gzip')



######plot marker gene #####

#marker_genes = ['MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ', 'PPBP']

marker_genes = ['PSG8','PAPPA','FLT1','STAT4','ERVFRD-1','DNMT1','HLA-G','VIM','PECAM1','CD68']
marker_genes = ['PSG8','PAPPA','CSHL1','FLT1','ENG','ERVFRD-1','DNMT1','HLA-G','VIM','PECAM1','CD68']


marker_genes = [
    "DNMT1", "CDH1", "MKI67",
    "FLT1", "CSHL1", "PSG8", 
    "ERVFRD-1", "LAIR2", "PLAC8",
    'VIM','PECAM1','CD14'
]

marker_genes = [ #try to distinguish STB terminals 
    "FLT1", 'ENG','INSIG2', #FLT1-enriched
    "CSHL1",'CSH1','PAPPA',  #PAPPA-enriched
    'PSG8','SH3TC2','GCM1' #general

]



marker_genes = [ #try to distinguish STB terminals 
    "DDX60", 'MAP4K4','SH3TC2', #FLT1-enriched
    "STAT4",'STAT5A','AR',  #PAPPA-enriched
    'LVRN','INHBA','MYCN' #general

]




marker_genes_full = { #from seurat code
    'Quality control' : ["tsse","n_fragment","frac_mito","frac_dup",'XIST'],
    'Trophoblast' : ['KRT7', 'GATA3', 'TFAP2A'],
    'CTB' : ['DNMT1', 'CDH1', 'PPARG', 'TEAD4', 'TEAD3', 'MKI67', 'TP53','TP63', 'TP73', 'BCAM'],
    'CTB fusion' : ['ERVFRD-1', 'GCM1', 'OVOL1','PPARD'], 
    'STB nascent' : ['SH3TC2', 'BACE2', 'ESRRG'], 
    'STB general' : ['PSG8', 'CGA', 'PSG2', 'PSG5', 'LEP'], 
    'STB PAPPA' : ['PAPPA', 'ADAMTSL1','ADAMTS6', 'GH2', 'GHR','JAK1', 'JAK2', 'LAMA3', 'AR', 'VDR', 'CSHL1', 'CSH1', 'CSH2', 'STAT5A', 'STAT5B', 'STAT4','FOS', 'FOSB', 'JUNB', 'JUN'], 
    'STB FLT1' : ['FLT1', 'ENG', 'ANGPTL4', 'FSTL3', 'INHBA', 'INHA','MYCN', 'POU2F3', 'LVRN', 'TGFB1', 'FOSL2', 'JUND'], 
    'STB apoptosis' : ['DDX60', 'DDX58', 'MAP4K4', 'SPATA5', 'GDF15', 'CROT', 'CDKN1A', 'ADCY5'],
    'EVT' : ['HLA-G', 'LAIR2', 'PLAC8', 'MMP2'], 
    'STR general' : ['VIM', 'DLK1'], 
    'Vascular Endothelial Cell' : ['PECAM1'], 
    'STR' : ['HIVEP3', 'HLA-A', 'HLA-DPA1', 'HLA-DPB1'], 
    'Mesenchymal STR': ['THY1'], 
    'Hofbauer Cell': ['CD68', 'CD14'], 
    'Red blood' : ['HBA1', 'HBZ']
}




#sc.settings.set_figure_params(dpi=100,figsize= (7.5,7.5) )
#sc.pl.umap(gene_matrix, use_raw=False, color=["leiden"] + marker_genes,vmin='p10',vmax='p99',size=10)
#sc.pl.violin(gene_matrix, marker_genes, use_raw=False, groupby='leiden')

####plot all marker gene list with for loop#######
for id in marker_genes_full.keys():
    genes = marker_genes_full[id]
    sc.settings.set_figure_params(dpi=100,figsize= (5.5,5.5) )
    for gene in genes:
        #sc.pl.umap(gene_matrix, use_raw=False, color=gene,vmin='p60',vmax='p100',size=8,color_map = my_cmap)
        #sc.pl.umap(gene_matrix, use_raw=False, color=gene,vmin='p50',vmax='p99',size=10,color_map = my_cmap)
        sc.pl.umap(gene_matrix, use_raw=False, color=gene,vmin='p10',vmax='p100',size=10,color_map = my_cmap)

        

####plot with marker gene grid with threshold cutoff (code from umap tuning funcgion)
thresh = {
    'ERVFRD-1':['p30','p100'],
    'PAPPA':['p65','p99'],
    'LAMA3':['p50','p100'],
    'STAT5A':['p30','p99'],
    'ESRRG':['p50','p99'],
    'FLT1':['p60','p100'],
    'ENG':['p60','p100'],
    'FOSL1':['p10','p99'],
    'JUNB':['p30','p99'],
    'FOS':['p30','p99'],
    'HLA-G':['p30','p100'],
    'PLAC8':['p50','p100'],
    
    'DNMT1':['p10','p100'],
    'PSG8':['p10','p100'],
    'CGA':['p10','p100'],
    'SH3TC2':['p60','p100'],
    'LEP':['p10','p100'],
    
}


#marker_genes = ["ERVFRD-1","PAPPA",'LAMA3','STAT5A','ESRRG','FLT1','ENG','FOSL1','JUNB','FOS','HLA-G','PLAC8']
marker_genes = ["DNMT1","ERVFRD-1",'PSG8','CGA','SH3TC2','LEP','PAPPA','FLT1']

marker_genes_sel = []
[ marker_genes_sel.append(i) for i in marker_genes if i in thresh.keys()]

##for plotting
plt.rcParams['figure.figsize'] =  (4,6)
plt.rcParams['legend.loc'] = 'upper left'
plt.rcParams['figure.dpi'] = 100

##for save
#plt.rcParams['savefig.directory'] = "pdfs/marker_genes"
plt.rcParams['savefig.dpi'] = 100

#sc.set_figure_params(fontsize = 15, frameon = False, dpi = 100, )

##get and set umap xlim ylim to make plotting slim

xmin = gene_matrix.obsm['X_umap'][:,0].min()
xmax = gene_matrix.obsm['X_umap'][:,0].max()
ymin = gene_matrix.obsm['X_umap'][:,1].min()
ymax = gene_matrix.obsm['X_umap'][:,1].max()

shrink = 2

res_p = sc.pl.umap(gene_matrix, 
           use_raw=False, 
           ncols=4,
           color= marker_genes_sel,
           vmin = [ thresh[i][0]  for i in marker_genes_sel],
           vmax= [ thresh[i][1]  for i in marker_genes_sel],
           size= 10,#15,
           color_map = my_cmap,#my_cmap_tfdev,#my_cmap,
           #colorbar_loc = None, only in scanpy 1.9
           show = False,
           title = [ i+ " cutoff: " + thresh[i][0] + " - " + thresh[i][1] for i in marker_genes_sel],
           return_fig = False, #return as a whole fig for save?
           #legend_fontsize = 20,
           #legend_fontweight = 10,
           legend_loc = 'left margin' #unwork
           #legend_loc = ['left margin','right margin','right margin','right margin','right margin','right margin','right margin','right margin','right margin','none','right margin','right margin'] #not work, do not know why
          )

#res_p[0]._remove_legend()
#res_p[0].legend(loc='upper left')

#plt.colorbar(res_p[0].contourf,shrink=0.5)#fraction=0.046, pad=0.04)
#plt.colorbar(res_f,ax=res_p[0], shrink=0.5)#,pad=0.01, fraction=0.08, aspect=30)

#######make slim#######
#plt.xlim([shrink*xmin,shrink*xmax]) #only the last one
#plt.ylim([shrink*ymin,shrink*ymax])

# plt.all_set_xlim(shrink*xmin,shrink*xmax) #error
# plt.all_set_ylim(shrink*ymin,shrink*ymax)

for i in range(len(res_p)):
    res_p[i].set_xlim(shrink*xmin,shrink*xmax)
    res_p[i].set_ylim(shrink*ymin,shrink*ymax)
    
#plt.margins(x = [0.7]*8, y = [0.8]*8) #error
#plt.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.3) #not work

plt.suptitle('Marker gene activity score',y = 1.0, x = 0.45)
#plt.tight_layout()

###save pdf ####
plt.savefig("pdfs/marker_genes/marker_genes.pdf",format='pdf',dpi = 100)

      
#########plotting by splitting of libraries########

gene_matrix.obs['library'].value_counts().sort_index()

placenta_10X_early1    3725
placenta_10X_early2    5262
placenta_10X_early3    4024
placenta_10X_early4    3619
placenta_10X_early5    1977
placenta_10X_early6    4208

# placenta_10X_early1    5334
# placenta_10X_early2    7788
# placenta_10X_early3    5876
# placenta_10X_early4    5583
# placenta_10X_early5    3237
# placenta_10X_early6    5968


for id in ['placenta_10X_early1', 'placenta_10X_early2', 'placenta_10X_early3', 'placenta_10X_early4', 'placenta_10X_early5', 'placenta_10X_early6']:
    sc.pl.umap(gene_matrix[gene_matrix.obs['library'] == id,], use_raw=False, color=["leiden"],size=10,title=id)
    sc.pl.umap(gene_matrix, use_raw=False, color=["PAPPA"],vmin='p10',vmax='p99',size=10)
    sc.pl.umap(gene_matrix, use_raw=False, color=["CSHL1"],vmin='p10',vmax='p99',size=10)
    sc.pl.umap(gene_matrix, use_raw=False, color=["FLT1"],vmin='p30',vmax='p99',size=10)
    sc.pl.umap(gene_matrix, use_raw=False, color=["ENG"],vmin='p50',vmax='p99',size=10)
    sc.pl.umap(gene_matrix, use_raw=False, color=["DDX60"],vmin='p50',vmax='p100',size=10)
    

######plotting by one gene #####   

sc.pl.umap(gene_matrix, use_raw=False, color=["leiden"],size=10, legend_loc = 'on data', )
#snap.pl.umap(gene_matrix, color="leiden", interactive=False, width = 500)


sc.pl.umap(gene_matrix, use_raw=False, color=["DNMT1"],vmin='p50',vmax='p99',size=10)
sc.pl.umap(gene_matrix, use_raw=False, color=["TP63"],vmin='p50',vmax='p99',size=10)


sc.pl.umap(gene_matrix, use_raw=False, color=["ERVFRD-1"],vmin='p30',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["ANXA1"],vmin='p10',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["OVOL1"],vmin='p30',vmax='p99',size=10,color_map = my_cmap)

sc.pl.umap(gene_matrix, use_raw=False, color=["ESRRG"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["PPARD"],vmin='p10',vmax='p99',size=10,color_map = my_cmap)


sc.pl.umap(gene_matrix, use_raw=False, color=["SH3TC2"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)


sc.pl.umap(gene_matrix, use_raw=False, color=["PSG8"],vmin='p10',vmax='p100',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["LEP"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["CGA"],vmin='p30',vmax='p100',size=10,color_map = my_cmap)


sc.pl.umap(gene_matrix, use_raw=False, color=["PAPPA"],vmin='p65',vmax='p100',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["CSHL1"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["CSH1"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)

 
sc.pl.umap(gene_matrix, use_raw=False, color=["FLT1"],vmin='p65',vmax='p100',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["ENG"],vmin='p60',vmax='p100',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["INHA"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["INHBA"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["TGFB1"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)


sc.pl.umap(gene_matrix, use_raw=False, color=["FOSL1"],vmin='p50',vmax='p100',size=20,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["FOSL2"],vmin='p20',vmax='p100',size=20,color_map = my_cmap)


sc.pl.umap(gene_matrix, use_raw=False, color=["CD68"],vmin='p60',vmax='p100',size=20,color_map = my_cmap)


sc.pl.umap(gene_matrix[~gene_matrix.obs['leiden'].isin(['20','19','17','18','13','14','1']) ,], use_raw=False, color=["FOSL1"],vmin='p20',vmax='p100',size=20,color_map = my_cmap)

sc.pl.umap(gene_matrix[~gene_matrix.obs['leiden'].isin(['20','19','17','18','13','14','1']) ,], use_raw=False, color=["AR"],vmin='p20',vmax='p100',size=20,color_map = my_cmap)


sc.pl.umap(gene_matrix, use_raw=False, color=["DDX60"],vmin='p50',vmax='p100',size=10,color_map = my_cmap)
#sc.pl.umap(gene_matrix, use_raw=False, color=["SPATA5"],vmin='p10',vmax='p99',size=10)

#sc.pl.umap(gene_matrix, use_raw=False, color=["DDX58"],vmin='p10',vmax='p99',size=10)
#sc.pl.umap(gene_matrix, use_raw=False, color=["MAP4K4"],vmin='p10',vmax='p99',size=10)

#sc.pl.umap(gene_matrix, use_raw=False, color=["PDCD7"],vmin='p50',vmax='p100',size=10)
#sc.pl.umap(gene_matrix, use_raw=False, color=["PDE4D"],vmin='p50',vmax='p100',size=10)
#sc.pl.umap(gene_matrix, use_raw=False, color=["CASP3"],vmin='p50',vmax='p100',size=10)
sc.pl.umap(gene_matrix, use_raw=False, color=["CASP9"],vmin='p50',vmax='p100',size=10)
#sc.pl.umap(gene_matrix, use_raw=False, color=["BCL2"],vmin='p50',vmax='p100',size=10)
#sc.pl.umap(gene_matrix, use_raw=False, color=["BIRC2"],vmin='p50',vmax='p100',size=10)
sc.pl.umap(gene_matrix, use_raw=False, color=["NAT2"],vmin='p50',vmax='p100',size=10)


sc.pl.umap(gene_matrix, use_raw=False, color=["HLA-G"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["PLAC8"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)


sc.pl.umap(gene_matrix, use_raw=False, color=["HBZ"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["HBA1"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color=["HBB"],vmin='p50',vmax='p99',size=10,color_map = my_cmap)


sc.pl.umap(gene_matrix, use_raw=False, color=["tsse"],vmin='p10',vmax='p99',size=10)
sc.pl.umap(gene_matrix, use_raw=False, color=["n_fragment"],vmin='p10',vmax='p99',size=10)
sc.pl.umap(gene_matrix, use_raw=True, color=["frac_mito"],vmin='p10',vmax='p99',size=10)
sc.pl.umap(gene_matrix, use_raw=True, color=["frac_dup"],vmin='p10',vmax='p99',size=10)

sc.pl.umap(gene_matrix, use_raw=False, color=["doublet_score"])

###plot any gene of interest

plt.rcParams['figure.figsize'] =  (4.5,4.5)
#plt.rcParams['legend.loc'] = 'upper left'
plt.rcParams['figure.dpi'] = 150


sc.pl.umap(gene_matrix, use_raw=False, color='CLIC5',vmin='p10',vmax='p100',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color='MAFG',vmin='p10',vmax='p100',size=10,color_map = my_cmap)

sc.pl.umap(gene_matrix, use_raw=False, color='SIRT5',vmin='p30',vmax='p100',size=10,color_map = my_cmap)
sc.pl.umap(gene_matrix, use_raw=False, color='GFOD1',vmin='p30',vmax='p100',size=10,color_map = my_cmap)


###save/reload gmat
gene_matrix.write('gene_matrix.h5ad')

gene_matrix = sc.read("gene_matrix.h5ad")



#############################output cluster.add.df#####################

#import plotly.express as px
#from natsort import index_natsorted
#import pandas as pd

embedding = data.obsm['X_umap'] 
embedding = data.obsm['X_umap_rotate'] 

groups = data.obs['leiden'].to_numpy()

# idx = index_natsorted(groups)
# embedding = embedding[idx, :]
# groups = [groups[i] for i in idx]


cluster_df_add = pd.DataFrame({
    "UMAP-1": embedding[:, 0],
    "UMAP-2": embedding[:, 1],
    "cluster":groups,
    #'color': groups,
   },
   index = data.obs_names

)
#27380 x 14
#33786 x 4

obs_df = data.obs#[:,:].to_pandas() #data.obs[:,:].to_pandas()
#obs_df.index = data.obs_names #33786 x 9

all(cluster_df_add.index == obs_df.index) #True

[True if i in obs_df.columns else False for i in cluster_df_add.columns ]
#Not in obs_df

cluster_df_add = pd.concat([cluster_df_add, obs_df], axis=1)
#27380 x 14
#33786 x 13


sum(cluster_df_add.columns.value_counts() >1) #0, no duplicated

(cluster_df_add.cluster == cluster_df_add.leiden ).all() #True


cluster_df_add.cluster = cluster_df_add.cluster.astype('category')

cluster_df_add.cluster.cat.categories
['0', '1', '10', '11', '12', '13', '14', '15', '2', '3', '4', '5', '6',
       '7', '8', '9']

cluster_df_add.cluster = cluster_df_add.cluster.cat.set_categories( ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']   )



# with open('cluster_df_add.pkl','wb') as fh:
#     dill.dump(cluster_df_add,file=fh)

with open('cluster_df_add.txt','w') as fh:
    cluster_df_add.to_csv(fh,index=True,sep="\t",header=True)


with open('cluster_df_add.filter_rawcluster.txt','w') as fh:
    cluster_df_add.to_csv(fh,index=True,sep="\t",header=True)



# sc.settings.set_figure_params(dpi=300,figsize= (5,5) ) #no use for plotly (js based)
# plt.rcParams['figure.figsize'] =  (5,5)

fig = px.scatter(
    cluster_df_add, x="UMAP-1", y="UMAP-2", color='cluster',
    color_discrete_sequence=color_good,#px.colors.qualitative.Dark24,
    title = 'SnapATAC2 cluster',
    #width = 550,
    #height = 450,
    category_orders = {#'leiden':groups.value_counts().index.to_list(),
                       'cluster':cluster_df_add.cluster.value_counts().index.to_list()
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

fig.show(width=450,height=450)#, interactive=False)

#fig.write_image("umap.pdf")


    
####################doDistri of cluster and merge cluster 5 and cluster 7?##################################  

# ##collect cluster_df_add #####

# cluster_df = data.obs.loc[:,['louvain','cell_type_louvain','leiden','cell_type_leiden']]

# umap = pd.DataFrame(adata.obsm['X_umap'],columns=['UMAP1','UMAP2'] )
# umap.index = cluster_df.index

# #cluster_df = pd.concat([cluster_df,umap],axis=1)

# idx = np.r_[0:8] #choose column idx
# metadata = adata.obs.iloc[:,idx]

# all(cluster_df.index==metadata.index) #True

# cluster_df_add = pd.concat([cluster_df,umap,metadata],axis=1)
# #save
# #cluster_df_add.to_csv('cluster_df_add.txt',index=True,sep='\t')

# ########



###plot cluster distribution and filter####


def dotDistri(data = cluster_df_add,cl_type='leiden'):
    #use seaborn scatterplot to visualize cluster distribution on UMAP
    
    cluster = data[cl_type]
    levels = cluster.cat.categories
    n_cl = cluster.cat.categories.size
    n_row = ceil (n_cl / 3)
    
    fig, (ax_list) = plt.subplots(n_row, 3, figsize=(20, 42), dpi=150, sharey=False)

    kwargs = {'edgecolor':'none', #for edge color
              #'linewidth':0, #line width of spot
              #'linestyle':'--',#line style of spot 
              's':10, #size? a matlibplot parameter

    }
    
    for i in range(0,n_cl):
        cl = levels[i]
        
        rid = int(i / 3)
        cid = i % 3
        
        print(cl,rid,cid)
        
        #the grey background
        g1 =sns.scatterplot( data['UMAP-1'], data['UMAP-2'],ax=ax_list[rid][cid], color='grey',**kwargs)
        
        #the selected cluster
        data_sel = data.loc[data[cl_type] == cl,]
        g2 =sns.scatterplot( data_sel['UMAP-1'], data_sel['UMAP-2'],ax=ax_list[rid][cid], color='red',**kwargs)
    
        #ax_list[rid][cid].title.set_text('cluster' + str(i) )
        ax_list[rid][cid].set_title('cluster' + cl,fontsize = 25 )
        ax_list[rid][cid].set_xlabel('UMAP-1',fontsize = 25 )
        ax_list[rid][cid].set_ylabel('UMAP-2',fontsize = 25 )
        
        ax_list[rid][cid].grid(False)
 
        #g1.legend(loc=0,ncol=1,bbox_to_anchor=(1,1)) #loc='upper right', will legend inside plot

        #fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')

    fig.tight_layout()
        
    return(1)




def clusterPlot(df ,feature1 = 'cluster',feature2=  'leiden',color_use = color_good, title = ['cluster','leiden'],size = 10, fontsize = 25, width = 13, height = 8, dpi = 100,shrink = 1.2, save = None):
    #use seaborn scatterplot to visualize cluster on UMAP, label on cluster
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(width, height), dpi= dpi, sharey=False)

    kwargs = {'edgecolor':'none', #for dot edge color
              #'linewidth':0, #line width of spot
              #'linestyle':'--',#line style of spot 
              's':size, #dot size, a matlibplot parameter
              #'legend': False

    }

    ##plt.rcParams['figure.figsize'] =  (20,15)
    #plt.rcParams['font.size'] = 45
    ##sc.set_figure_params(fontsize = 15, frameon = False, dpi = 100, )
    
    
    ##for plotting
    ##plt.rcParams['figure.figsize'] =  (4,6)
    #plt.rcParams['legend.loc'] = 'on data'
    ##plt.rcParams['figure.dpi'] = dpi

    ##for save
    #plt.rcParams['savefig.directory'] = "pdfs/marker_genes"
    plt.rcParams['savefig.dpi'] = dpi

    #sc.set_figure_params(fontsize = 15, frameon = False, dpi = 100, )

    if df[feature1].dtype != 'category':
        print('feature1 not category')
        return('Error')
    if df[feature2].dtype != 'category':
        print('feature2 not category')
        return('Error')
    
    n_color1 = df[feature1].cat.categories.size
    n_color2 = df[feature2].cat.categories.size
    
    
    if n_color1 > len(color_use) or n_color2 > len(color_use):
        print('not enough colors')
        return('fail')
    
    ##get centroid of each cluster
    center_x = df.groupby(feature1)['UMAP-1'].mean()
    center_y = df.groupby(feature1)['UMAP-2'].mean()
    
    if (center_x.size == center_y.size) & (center_x.index == center_y.index).all():
        centers = pd.DataFrame( { 'x': center_x, 'y': center_y, 'cluster': center_x.index  }  )
    
    ########add legend outline (text halo) (code from R plotting )
    centers_shift = pd.DataFrame({},columns = ['cluster','x','y'])
    
    #add little shift for text x and y, to plot text halo, borrow from snapATAC
    theta= np.linspace(0, 2*pi, num=50)
    r=0.1
    strwidth = 0.5 
    strheight = 0.5
    xo = r*strwidth # r*strwidth('A')
    yo = r*strheight #r*strheight('A')
    for i in range(len(centers)):
        for j in theta :
            centers_shift = pd.concat([centers_shift,
                                      pd.DataFrame(
                                          {'cluster':[ centers['cluster'][i] ],
                                          'x':[ centers['x'][i] + cos(j)*xo ], 
                                          'y':[ centers['y'][i] + sin(j)*yo ]
                                          },

                                         )],
                                      #axis = 0
                           )
    centers_shift.index = [  "row"+str(i+1) for i in range(centers_shift.shape[0]) ]
    
    
    g1 =sns.scatterplot( df['UMAP-1'], df['UMAP-2'],ax=ax1, hue=df[feature1],palette=color_use[:n_color1],legend = False,**kwargs)
    g2 = sns.scatterplot(df['UMAP-1'], df['UMAP-2'],ax=ax2, hue=df[feature2],palette=color_use[:n_color2],legend = 'full',**kwargs)
    #g3 = sns.scatterplot( center_x, center_y,ax=ax2, color = 'black',size= 25)
    
 
    for i in range(centers_shift.shape[0]):
            g1.text(centers_shift['x'].to_list()[i],
                    centers_shift['y'].to_list()[i],
                    centers_shift['cluster'].to_list()[i],
                    horizontalalignment='center',
                    size=fontsize + 0.5,
                    color='white',
                    weight = 'bold',
                    alpha = 1
                   )
            
    if center_x.size == center_y.size:
        for i in range(center_x.size):
            g1.text(center_x[i],
                    center_y[i],
                    center_x.index[i],
                    horizontalalignment='center',
                    size=fontsize,
                    color='black',
                    weight = 'bold'
                   )

    #g3 = sns.scatterplot( data['UMAP_1_scanpy'], data['UMAP_2_scanpy'],ax=ax3, hue=data["cluster"],palette=color_set,hue_order=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12','13','14','15','16','17'],**kwargs)

    

    
    
    #ax1.title.set_text('louvain clusters')
    #ax2.title.set_text('leiden clusters')
    #ax3.title.set_text('Scanpy with Seurat clusters')

    ax1.set_title(title[0],fontsize = 25 )
    ax1.set_axis_off()
    ax1.set_xlabel('UMAP-1',fontsize = 25 )
    ax1.set_ylabel('UMAP-2',fontsize = 25 )
    
    ax2.set_title(title[1],fontsize = 25 ,fontdict = {'horizontalalignment':'left'} )
    ax2.set_axis_off()
    ax2.set_xlabel('UMAP-1',fontsize = 25 )
    ax2.set_ylabel('UMAP-2',fontsize = 25 )
    ax2.legend(loc = 'upper right',fontsize = 15,markerscale = 2, title = 'Cell Type',title_fontsize = 'x-large')
    #move legend to right
    #sns.move_legend(ax2, 'upper right', bbox_to_anchor = (1,1))
    
    ax1.grid(False)
    ax2.grid(False)
    #ax3.grid(False)

    ##get and set umap xlim ylim to make plotting slim
    xmin = df['UMAP-1'].min()
    xmax = df['UMAP-1'].max()
    ymin = df['UMAP-2'].min()
    ymax = df['UMAP-2'].max()

    shrink = 2
    
    ##shrink/shift with xlim/ylim
    ax1.set_xlim(shrink*xmin,shrink*xmax)
    ax1.set_ylim(shrink*ymin,shrink*ymax)
    
    ax2.set_xlim(shrink*xmin,shrink*xmax)
    ax2.set_ylim(shrink*ymin,shrink*ymax)
    
    #plt.legend(fontsize = 55)
    
    #g1.legend() #loc='upper right', will legend inside plot
    #g2.legend(loc='upper right',ncol=1,bbox_to_anchor=(1,1))
    ##g3.legend(loc=0,ncol=1,bbox_to_anchor=(1,1))
    #fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='x-large')

    #fig.tight_layout()

    #fig.show()
    if save != None:
        fig.savefig(save,dpi = dpi)
    
    
    return('done.')

##



###
dotDistri(data = cluster_df_add,cl_type = 'cluster')

clusterPlot(df = cluster_df_add, feature1 = 'cluster', feature2 = 'leiden',color_use = color_good, title = 'placenta early pregnancy snATAC-seq'+"\n"+"total nuclei:" + str(cluster_df_add.shape[0]) )


######filter louvain cluster ? in cluster_df_add


cluster_df_add.cluster.value_counts()
0     4018
1     3668
2     2981
3     2545
4     2276
5     2140
6     1670
7     1649
8     1477
9     1017
10     994
11     862
12     720
13     534
14     470
15     359


##merge c7 to c5
map_cl_dict = {
#     for i in range(16):
#     print ("\'%i\': \'%i\'," % (i, i))
    
    '0': '0',
    '1': '1',
    '2': '2',
    '3': '3',
    '4': '4',
    '5': '5',
    '6': '6',
    '7': '5', #merge to c5
    '8': '8',
    '9': '9',
    '10': '10',
    '11': '11',
    '12': '12',
    '13': '13',
    '14': '14',
    '15': '15'
   
}


sum(cluster_df_add.cluster == '7') #1649

cluster_df_add.cluster = cluster_df_add.cluster.map(map_cl_dict)

cluster_df_add.cluster.value_counts().index.astype('int').sort_values().astype('string')
'0', '1', '2', '3', '4', '5', '6', '8', '9', '10', '11', '12', '13','14', '15'

map_cl_dict_mod = {
#    for i in cluster_df_add.cluster.value_counts().index.astype('int').sort_values().astype('string'):
#        print ("\'%s\': \'%s\'," % (i, i))
'0': '1',
'1': '2',
'2': '3',
'3': '4',
'4': '5',
'5': '6',
'6': '7',
'8': '8',
'9': '9',
'10': '10',
'11': '11',
'12': '12',
'13': '13',
'14': '14',
'15': '15',
   
}


cluster_df_add.cluster = cluster_df_add.cluster.map(map_cl_dict_mod)

cluster_df_add.cluster = cluster_df_add.cluster.astype('category').cat.set_categories( ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']   )

cluster_df_add.cluster.cat.categories
'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
       '14', '15'
    
cluster_df_add.cluster.value_counts().sort_index()
1     4018
2     3668
3     2981
4     2545
5     2276
6     3789
7     1670
8     1477
9     1017
10     994
11     862
12     720
13     534
14     470
15     359


cluster_df_add.cluster.cat.remove_unused_categories(inplace=True)

dotDistri(data = cluster_df_add,cl_type = 'cluster')



# cluster_df_add_filter = cluster_df_add.loc[ ~((cluster_df_add.louvain == '3') | (cluster_df_add.louvain == '14') | ( cluster_df_add.louvain == '20') ),] #12419

# ###cluster_df_add_filter.louvain.droplevel(['3','14','20']), #not this
# cluster_df_add_filter.louvain.cat.remove_unused_categories(inplace=True)
# cluster_df_add_filter.louvain.value_counts() 
    

# cluster_df_add_filter.leiden.cat.remove_unused_categories(inplace=True)
# cluster_df_add_filter.leiden.value_counts()



#######de-noise by removing dots far from 98% quantile from the centroid (copy from SnapATAC R code)####


def dotDist(cluster = cluster_df_add.loc[:,['cluster','UMAP-1','UMAP-2']], id = '5', center = centers, q = 0.98):#, ax = None):
    colids = cluster.columns #must 'cluster','dim1','dim2' and rowname
    cluster_sel = cluster.loc[ cluster.loc[:,'cluster'] == id,:]
    n_sel = cluster_sel.shape[0]
    dx = cluster_sel.loc[:,'UMAP-1'] - center.loc[center['cluster'] == id,'x'].item()
    dy = cluster_sel.loc[:,'UMAP-2'] - center.loc[center['cluster'] == id,'y'].item()
    d = np.sqrt(dx**2 + dy**2)
    #d.name = 'dist'
    
    plt.rcParams['figure.figsize'] = (5.5,4.5)
    res_h = sns.histplot(data=d,x=None,bins=100)#, ax = ax)
    d_quantile = d.quantile(np.arange(0,1,0.05)).to_list()
    d_quantile_q = d.quantile(q)
    
    for i in d_quantile:
        res_h.axvline(x=i,ymin=0,ymax=1,color='red',linewidth=0.5,linestyle='dashed')
    res_h.grid(False)
    (ymin,ymax) = res_h.get_ylim()
    res_h.text(x=d_quantile_q,y=ymax*q,s="%%%.1f percentile is %.3f in cl %s" % (q*100,d_quantile_q,id)  )
    
    print("cluster %s ok\n" % id )
    plt.show()
    
    
    return(d_quantile_q)


# ##plot distance distribution hist
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




##get centroid of each cluster
center_x = cluster_df_add.groupby('cluster')['UMAP-1'].mean()
center_y = cluster_df_add.groupby('cluster')['UMAP-2'].mean()

if (center_x.size == center_y.size) & (center_x.index == center_y.index).all():
    centers = pd.DataFrame( { 'x': center_x, 'y': center_y, 'cluster': center_x.index  }  )

#centers



####plot dot distance distribution##
#fig, axes = plt.subplots(5, 3, figsize=(15, 20), dpi=150, sharey=False)

kwargs = {'edgecolor':'none', #for dot edge color
          #'linewidth':0, #line width of spot
          #'linestyle':'--',#line style of spot 
          's':10, #size? a matlibplot parameter
          #'legend': False

}

levels = cluster_df_add.cluster.cat.categories
n_cl = cluster_df_add.cluster.cat.categories.size
#n_row = ceil (n_cl / 3)

dist_q_list = {}

for i in range(n_cl):
    cl = levels[i]
    print("do dotDist for cluster %s  " % cl)
    
    
    #get row and col id
    #rid = int(i / 3)
    #cid = i % 3

    print( "cluster is %s" % cl)
    print( "cluster is %s, row is %i, column is %i" % (cl,rid,cid))

    
    dist_q = dotDist(cluster = cluster_df_add.loc[:,['cluster','UMAP-1','UMAP-2']], id = cl, center = centers, q = 0.98)#, ax = axes[rid][cid])
    dist_q_list[cl] = round(dist_q,3)

#fig.show()
#fig.tight_layout()   



#dist_q_list
# q95.5
{'1': 2.148,
 '2': 2.38,
 '3': 2.481,
 '4': 2.482,
 '5': 4.1,
 '6': 2.137,
 '7': 2.0,
 '8': 1.912,
 '9': 1.725,
 '10': 7.561,
 '11': 4.582,
 '12': 1.429,
 '13': 1.442,
 '14': 1.802,
 '15': 1.248}

## q96
{'1': 2.174,
 '2': 2.436,
 '3': 2.517,
 '4': 2.509,
 '5': 4.292,
 '6': 2.173,
 '7': 2.023,
 '8': 1.964,
 '9': 1.746,
 '10': 7.583,
 '11': 4.616,
 '12': 1.437,
 '13': 1.466,
 '14': 1.822,
 '15': 1.264}


## q98
{'1': 2.313,
 '2': 2.739,
 '3': 2.788,
 '4': 2.628,
 '5': 4.794,
 '6': 2.387,
 '7': 2.333,
 '8': 3.95,
 '9': 1.855,
 '10': 7.773,
 '11': 4.717,
 '12': 1.504,
 '13': 1.611,
 '14': 1.924,
 '15': 1.344}

## q99
{'1': 2.492,
 '2': 3.143,
 '3': 3.105,
 '4': 2.818,
 '5': 5.104,
 '6': 2.541,
 '7': 8.58,
 '8': 4.718,
 '9': 1.945,
 '10': 7.947,
 '11': 4.766,
 '12': 1.536,
 '13': 1.711,
 '14': 1.976,
 '15': 14.212}


## q95
{'1': 2.127,
 '2': 2.331,
 '3': 2.428,
 '4': 2.447,
 '5': 3.917,
 '6': 2.108,
 '7': 1.975,
 '8': 1.85,
 '9': 1.692,
 '10': 7.53,
 '11': 4.534,
 '12': 1.398,
 '13': 1.437,
 '14': 1.774,
 '15': 1.23
}

# dist_q99_list_manual = {'1': 2.127,
#                          '2': 2.331,
#                          '3': 2.428,
#                          '4': 2.447,
#                          '5': 3.917,
#                          '6': 2.108,
#                          '7': 1.975,
#                          '8': 1.85,
#                          '9': 1.692,
#                          '10': 7.53,
#                          '11': 4.534,
#                          '12': 1.398,
#                          '13': 1.437,
#                          '14': 1.774,
#                          '15': 1.23
#                         }



def dotClean(cluster = None, id = None, center = None, d_filter = dist_q_list):#, ax = None):
    colids = cluster.columns #must 'cluster','dim1','dim2' and rowname
    cluster_sel = cluster.loc[ cluster.loc[:,'cluster'] == id,:]
    d_filter_sel = d_filter[id]
    n_sel = cluster_sel.shape[0]
    dx = cluster_sel.loc[:,'UMAP-1'] - center.loc[center['cluster'] == id,'x'].item()
    dy = cluster_sel.loc[:,'UMAP-2'] - center.loc[center['cluster'] == id,'y'].item()
    d = np.sqrt(dx**2 + dy**2)
    #d.name = 'dist'
    
    kwargs = {'edgecolor':'none', #for dot edge color
          #'linewidth':0, #line width of spot
          #'linestyle':'--',#line style of spot 
          's':2, #size? a matlibplot parameter
          #'legend': False

    }
    
    flag_dist = d > d_filter_sel
    flag_dist_all = cluster.index.isin(cluster_sel.index[flag_dist])
    
    plt.rcParams['figure.figsize'] = (5.5,4.5)
    #the grey background (all clusters)
    g0 =sns.scatterplot( cluster['UMAP-1'], cluster['UMAP-2'], color='grey',**kwargs)

    #the selected cluster
    g1 =sns.scatterplot( cluster_sel['UMAP-1'], cluster_sel['UMAP-2'], color='red',**kwargs)
    
    #filtered dots in the selected cluster
    g2 =sns.scatterplot( cluster_sel.loc[flag_dist,'UMAP-1'], cluster_sel.loc[flag_dist,'UMAP-2'],color='blue',**kwargs)

    #kept dots in the selected cluster
    #g3 =sns.scatterplot( cluster_sel[~flag_dist,'UMAP-1'], cluster_sel[~flag_dist,'UMAP-2'], color='blue',**kwargs)
    
    #the centroid
    g3 =sns.scatterplot( center.loc[center['cluster'] == id,'x'], center.loc[center['cluster'] == id,'y'],color='black',s = 8)
    
    #the filtered dots
    
    
    #ax_list[rid][cid].title.set_text('cluster' + str(i) )
    g1.set_title('cluster' + id,fontsize = 15 )
    g1.set_xlabel('UMAP-1',fontsize = 15 )
    g1.set_ylabel('UMAP-2',fontsize = 15 )

    g1.grid(False)
    
    
    print("cluster %s ok\n" % id )
    plt.show()
    
    
    return(flag_dist_all)



res_flag_dist = {}

for i in levels:
    res_flag_dist[i] = dotClean(cluster = cluster_df_add.loc[:,['cluster','UMAP-1','UMAP-2']], 
                                id = i,
                                center = centers,
                                d_filter = dist_q_list
                               )



##merge dist_flag to filter cluster_df_add
flag_dist_combine = pd.DataFrame(res_flag_dist).any(axis = 1)#  apply(do.call(cbind,res.flag.dist),1,any)

flag_dist_combine.value_counts()
False    26824
True       556


False    26140
True      1240


len(flag_dist_combine) == cluster_df_add.shape[0] #True

cluster_df_add_filterdist = cluster_df_add.loc[(~flag_dist_combine).to_list(),:]
#26824 x 14 #q98
#26140 x 14 #q95.5


clusterPlot(df = cluster_df_add_filterdist, feature1 = 'cluster', feature2 = 'leiden',color_use = color_good, title = 'placenta early pregnancy snATAC-seq'+"\n"+"total nuclei:" + str(cluster_df_add_filterdist.shape[0]), size = 5 )


# table(flag.dist.combine)
# #flag.dist.combine
# FALSE  TRUE 
# 13722   634 

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
# for(i in c('8','3','4','5','7','1','4') ){
#   dotDist(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers)
# }

# # par(mfrow=c(2,3))
# # options(repr.plot.height=10,repr.plot.width=15)
# # for(i in c('12','9','13','15','11','14') ){
# #   dotDist(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i,center = centers)
# # }

# d.filter <- c( '1'=2.5, '2'=0,'3'=2,'4'=2,'5'=3,'6'=0,
#               '7'=2,'8'=2,'9'=0,'10'=0,'11'=0,'12'=0,
#               '13'=0,'14'=0,'15'=0 )


# par(mfrow=c(3,3))
# res.flag.dist <- list()
# options(repr.plot.height=15,repr.plot.width=15)
# for(i in c('8','3','4','5','7','1','4') ){
#   res.flag.dist[[i]] <- dotClean(cluster = cluster.df.add[,c('cluster','UMAP_1','UMAP_2')], id = i,
#            center = centers,d.filter = d.filter)
# }


# flag.dist.combine <- apply(do.call(cbind,res.flag.dist),1,any)
# table(flag.dist.combine)
# #flag.dist.combine
# FALSE  TRUE 
# 13722   634 






########## create a dictionary to map cluster to annotation label
cluster2annotation = {
     '1': '1: CTB-1',
     '2': '2: STB Nascent',
     '3': '3: STB Mature 2 (FLT1)',
     '4': '4: STB Premature 1 (PAPPA)',
     '5': '5: STB Mature 1 (PAPPA)',
     '6': '6: STB Mixed',
     '7': '7: CTB-2',
     '8': '8: STB Premature 1 (FLT1)',
     '9': '9: STR-5',
     '10': '10: EVT',
     '11': '11: CTB Fusion',
     '12': '12: STR-2',
     '13': '13: STR-4',
     '14': '14: STR-3',
     '15': '15: STR-1',

}

# add a new `.obs` column called `cell type` by mapping clusters to annotation using pandas `map` function

cluster_df_add['cell_type_leiden'] = cluster_df_add.cluster.map(cluster2annotation)

##plot the annotated umap
clusterPlot(df = cluster_df_add, feature1 = 'cluster' ,feature2 = 'cell_type_leiden', color_use = color_good, title = 'placenta early pregnancy snATAC-seq'+"\n"+"total nuclei:" + str(cluster_df_add.shape[0]))



#####write the cell annotated custer_df_add ####
with open('cluster_df_add.final.final.txt','w') as fh:
    cluster_df_add.to_csv(fh,index=True,sep="\t",header=True)



##filter and write a new meta-column to snapatac2 obj
data.obs['cell_type_leiden'] = cluster_df_add.cell_type_leiden

sc.pl.umap(data, color='cell_type_leiden', legend_loc='on data', title='cell type of earyly pregenancy six donors leiden cluster', frameon=False,palette=color_palette, size=26,legend_fontsize=5)#, save='.pdf')




########Peak calling at the cluster-level#######

snap.tl.call_peaks(data, groupby="leiden")


####get peak matrix for all cell type#####

peak_mat = snap.pp.make_peak_matrix(data)#, file="peak_matrix.h5ad")

AnnData object with n_obs × n_vars = 27380 × 671430
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
data


AnnData object with n_obs × n_vars = 33786 × 684610
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library'



peak_mat.write('peak_matrix.final.h5ad')

peak_mat = sc.read("peak_matrix.final.h5ad")



peak_mat.write('peak_matrix.h5ad')

peak_mat = sc.read("peak_matrix.h5ad")


##output combined peaks 

#with open('peak.combine.bed','w') as fh:

with open('peak.combine.final.bed','w') as fh:
    for i in peak_mat.var_names:
        (chr,range) = i.split(":")
        (start, end) = range.split("-")
        fh.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + i + "\n")


        
        
##output bw tracks (from Tn5 insertion tracks, in obsm), without any filtering of bin range
data.obs['leiden'].value_counts()
0     4018
1     3668
2     2981
3     2545
4     2276
5     2140
6     1670
7     1649
8     1477
9     1017
10     994
11     862
12     720
13     534
14     470
15     359

cluster_df_add.leiden.value_counts()
0     4018
1     3668
2     2981
3     2545
4     2276
5     2140
6     1670
7     1649
8     1477
9     1017
10     994
11     862
12     720
13     534
14     470
15     359

cluster_df_add.cluster.value_counts().sort_index()
1     4018
2     3668
3     2981
4     2545
5     2276
6     3789
7     1670
8     1477
9     1017
10     994
11     862
12     720
13     534
14     470
15     359


(data.obs_names == cluster_df_add.index).all() #True, 27380


data.obs['cluster'] = cluster_df_add.cluster

#snap.ex.export_bigwig(data,out_dir = "./bw", prefix = "cluster",groupby = 'leiden' )
#snap.ex.export_bed(data,out_dir = "./bed", prefix = "cluster",groupby = 'leiden' )
 
snap.ex.export_bigwig(data,out_dir = "./bw", prefix = "cluster",groupby = 'cluster' )
#snap.ex.export_bed(data,out_dir = "./MACS", prefix = "cluster",groupby = 'cluster' )

##will normalized! in snapatac2-core/src/export.rs::export_insertions_as_bigwig
#norm_factor = total_count * resolution



#################get cstb obj#####################3

clusterPlot(df = cluster_df_add_filterdist, feature1 = 'cluster', feature2 = 'leiden',color_use = color_good, title = 'placenta early pregnancy snATAC-seq'+"\n"+"total nuclei:" + str(cluster_df_add_filterdist.shape[0]), size = 5 )
cluster_df_add_cstb = cluster_df_add_filterdist.loc[cluster_df_add_filterdist.cluster.isin(['7','1','11','2','6','8','3','4','5']),:]
#22815

cluster_df_add_cstb.cluster.value_counts().sort_index()
1     3937
6     3713
2     3594
3     2921
4     2493
5     2230
7     1636
8     1447
11     844
9        0
10       0
12       0
13       0
14       0
15       0

cluster_df_add_cstb.cluster = cluster_df_add_cstb.cluster.cat.remove_unused_categories()

cluster_df_add_cstb.cluster.value_counts().sort_index()
1     3937
2     3594
3     2921
4     2493
5     2230
6     3713
7     1636
8     1447
11     844

map_dict_cstb = {
    '1': '1',#     3937
    '2': '2',#     3594
    '3': '3',#     2921
    '4' : '4',#    2493
    '5': '5',#     2230
    '6' : '6',#    3713
    '7': '7',#     1636
    '8' : '8',#    1447
    '11' : '9',#    844


    
    
}

cluster_df_add_cstb.cluster = cluster_df_add_cstb.cluster.map(map_dict_cstb)

cluster_df_add_cstb.cluster.value_counts().sort_index()
1    3937
2    3594
3    2921
4    2493
5    2230
6    3713
7    1636
8    1447
9     844



cluster2annotation_cstb = {
     '1': '1: CTB-1',
     '2': '2: STB Nascent',
     '3': '3: STB Mature 2 (FLT1)',
     '4': '4: STB Premature 1-2 (PAPPA)',
     '5': '5: STB Mature 1 (PAPPA)',
     '6': '6: STB Premature 1-1 (PAPPA)',
     '7': '7: CTB-2',
     '8': '8: STB Premature 2 (FLT1)',
     '9': '9: CTB Fusion'

}


cluster_df_add_cstb['cell_type_leiden'] = cluster_df_add_cstb.cluster.map(cluster2annotation_cstb)


clusterPlot(df = cluster_df_add_cstb, feature1 = 'cluster', feature2 = 'cell_type_leiden',color_use = color_good, title = 'placenta early pregnancy snATAC-seq (cstb)'+"\n"+"total nuclei:" + str(cluster_df_add_cstb.shape[0]), size =6, width = 10, height = 8, dpi =150 )


all(cluster_df_add_cstb.index.isin(data.obs_names))#True

data_cstb = data[cluster_df_add_cstb.index,:]

View of AnnData object with n_obs × n_vars = 22815 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster'
    var: 'selected'
    uns: 'reference_sequences', 'AnnDataSet', 'spectral_eigenvalue', 'peaks'
    obsm: 'X_spectral', 'insertion', 'X_umap', 'X_spectral_harmony', 'X_spectral_mnn', 'X_umap_rotate'
    obsp: 'distances'
    
data_cstb.obs.cluster.value_counts().sort_index()


data_cstb.obs.cluster = cluster_df_add_cstb.cluster
data_cstb.obs['cell_type_leiden'] = cluster_df_add_cstb.cell_type_leiden


#####write the cell annotated custer_df_add_cstb ####
with open('cluster_df_add_cstb.txt','w') as fh:
    cluster_df_add_cstb.to_csv(fh,index=True,sep="\t",header=True)

    
(cluster_df_add_cstb.index == data_cstb.obs_names).all() #True

(cluster_df_add_cstb.cluster == data_cstb.obs['cluster']).all() #True
(cluster_df_add_cstb['UMAP-1'] == data_cstb.obsm['X_umap_rotate'][:,0]).all()#True
(cluster_df_add_cstb['UMAP-2'] == data_cstb.obsm['X_umap_rotate'][:,1]).all()#True


####save/reload data_cstb
data_cstb.write('PLA_atac_early_combined_cstb.h5ad')



##save X_spectral_harmony###

X_spectral_harmony_df = pd.DataFrame(data_cstb.obsm['X_spectral_harmony'],index=data_cstb.obs_names, columns=[ 'DC' + str(i) for i in np.arange(1,19,1)] )
#22815 x 18

with open('X_spectral_harmony_cstb.txt','w') as fh:
    X_spectral_harmony_df.to_csv(fh,index=True,sep="\t",header=True)


    



##################get cluster_df_add and data_cstb again with filtering c9 ervfrd-1 low cell###################


##filter c11 fusion for ERVFRD-1 low accessible cell###

gene_matrix

#AnnData object with n_obs × n_vars = 22815 × 59265
# obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
# var: 'n_cells'
# uns: 'log1p', 'leiden_colors', 'cluster_colors', 'sample_colors', 'sex_colors'
# obsm: 'X_umap'
    
gene_mat_ervfrd1 = gene_matrix[gene_matrix.obs['cluster'] == '9','ERVFRD-1']
#844 x 1

sc.pl.umap(gene_mat_ervfrd1, use_raw=False, color=["cluster"],size=10, legend_loc = 'on data', legend_fontoutline = 2,add_outline = True, )

#sc.pl.umap(gene_matrix, use_raw=False, color=["ERVFRD-1"],size=10 )
sc.pl.umap(gene_mat_ervfrd1, use_raw=False, color=["ERVFRD-1"],size=20 )


# #############filterby ervftd-1 accessibilitye#######
# # value = [i[0] for i in gene_mat_ervfrd1.X.tolist()]
# # #844, 1
# # value = np.array(value)

# value = gene_mat_ervfrd1.X.flatten()

# value_q = np.quantile(value,np.arange(0,1,0.05))

# # ArrayView([-0.80848403, -0.59739709, -0.51320473, -0.42741213,
# #            -0.3331489 , -0.25505869, -0.15756946, -0.06712278,
# #            -0.01581164,  0.04769488,  0.12311837,  0.16288468,
# #             0.21341889,  0.26237091,  0.32247899,  0.36806871,
# #             0.44418396,  0.54453369,  0.64108562,  0.97443484])

# # array([0.  , 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 ,
# #        0.55, 0.6 , 0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95])

# # ArrayView([0.15134259, 0.24904549, 0.29763578, 0.36896615, 0.43872666,
# #        0.51344039, 0.56673064, 0.61303063, 0.66800017, 0.6973906 ,
# #        0.73032653, 0.76164656, 0.78583355, 0.80562185, 0.82963589,
# #        0.85225077, 0.86984333, 0.91109843, 0.96803587, 1.07474329])

# flag = value > value_q[1] #5%
# flag = value > value_q[2] #10%
# flag = value > value_q[3] #15%

# gene_mat_ervfrd1_sel = gene_mat_ervfrd1[flag,:]

# #717 x 1 
# #759 x 1
# #801 x 1

# sc.pl.umap(gene_mat_ervfrd1_sel, use_raw=False, color=["cluster"],size=10, legend_loc = 'on data', legend_fontoutline = 2,add_outline = True, )
# sc.pl.umap(gene_mat_ervfrd1_sel, use_raw=False, color=["ERVFRD-1"],size=20 )



#################filter by umap xmax####
value = gene_mat_ervfrd1.obsm['X_umap'][:,0]


value.max() #4.13

flag = value < 3.8


gene_mat_ervfrd1_sel = gene_mat_ervfrd1[flag,:]
#30 of 814 removed

sc.pl.umap(gene_mat_ervfrd1_sel, use_raw=False, color=["cluster"],size=10, legend_loc = 'on data', legend_fontoutline = 2,add_outline = True, )
sc.pl.umap(gene_mat_ervfrd1_sel, use_raw=False, color=["ERVFRD-1"],size=20 )


flag = value >= 3.8

barcode_c9_rm = gene_mat_ervfrd1.obs_names[flag] #30

for i in barcode_c9_rm:
    print (i)

# placenta_10X_early1:TTACGGAAGAGCGAAA-1
# placenta_10X_early2:GATTGACAGCTCGTTA-1
# placenta_10X_early2:GCGGGTTGTGTCCCAG-1
# placenta_10X_early2:GCGTAGCCAAAGGTCG-1
# placenta_10X_early2:GCTCGAGAGCTACGCC-1
# placenta_10X_early2:TCGATTTTCATTGCCC-1
# placenta_10X_early4:AGCCTCTTCGTTGTAG-1
# placenta_10X_early4:CACAACATCGACTCGG-1
# placenta_10X_early4:TTGCACCTCGATGCAT-1
# placenta_10X_early5:AACAAAGAGGAGAACA-1
# placenta_10X_early6:ACAGAAAGTTCTACGA-1
# placenta_10X_early6:ACAGACTCACTACACA-1
# placenta_10X_early6:ACTGTCCGTCAGGTGA-1
# placenta_10X_early6:AGATAGACAACTACTG-1
# placenta_10X_early6:AGTGCCGGTTTCGTTT-1
# placenta_10X_early6:CAAGAAATCAAATGGA-1
# placenta_10X_early6:CCGCATTTCCTTTGCG-1
# placenta_10X_early6:CGCACAGAGTGATATG-1
# placenta_10X_early6:CGCGCAATCTTGTCGC-1
# placenta_10X_early6:CGTAAACAGCGTTGCC-1
# placenta_10X_early6:GCAGCCACAAACTACC-1
# placenta_10X_early6:GGTAGGAGTGCAAGAC-1
# placenta_10X_early6:GTGCACGGTCTCAAAC-1
# placenta_10X_early6:TAAACCGCAAGCAATA-1
# placenta_10X_early6:TCAAGACGTACGTATC-1
# placenta_10X_early6:TCAGCTCCACTAAACC-1
# placenta_10X_early6:TCTAGTTCATCATGTG-1
# placenta_10X_early6:TGATTTCGTATACGCT-1
# placenta_10X_early6:TGCCTGTGTTTCTCTA-1
# placenta_10X_early6:TTGCTTATCAGGAATA-1    
    
    
###get cluster_df_add and data_cstb

all(cluster_df_add.index == data_cstb.obs_names) #True
all(cluster_df_add.index == gene_matrix.obs_names) #TRUE
all(data_cstb.obs_names == gene_matrix.obs_names) #TRUE


[True if i in cluster_df_add.index else False for i in barcode_c9_rm] #all True


flag = [False if i in barcode_c9_rm else True for i in cluster_df_add.index]
cluster_df_add_filter_c9 = cluster_df_add.iloc[flag,:] #22785 of 22815


data_cstb_filter_c9 = data_cstb[cluster_df_add_filter_c9.index,:]
gene_matrix_filter_c9 = gene_matrix[cluster_df_add_filter_c9.index,:]



all(cluster_df_add_filter_c9.index == data_cstb_filter_c9.obs_names) #True
all(cluster_df_add_filter_c9.index == gene_matrix_filter_c9.obs_names) #TRUE
all(data_cstb_filter_c9.obs_names == gene_matrix_filter_c9.obs_names) #TRUE



sc.pl.umap(gene_matrix_filter_c9, use_raw=False, color=["cluster"],size=10, legend_loc = 'on data', legend_fontoutline = 2,add_outline = True, )

#sc.pl.umap(gene_matrix, use_raw=False, color=["ERVFRD-1"],size=10 )
sc.pl.umap(gene_matrix_filter_c9, use_raw=False, color=["ERVFRD-1"],size=20 )


##
cluster_df_add = cluster_df_add_filter_c9

data_cstb = data_cstb_filter_c9

gene_matrix = gene_matrix_filter_c9

del cluster_df_add_filter_c9
del  data_cstb_filter_c9
del gene_matrix_filter_c9


data_cstb.write('snapshot_h5ad/data_cstb.h5ad')
gene_matrix.write('snapshot_h5ad/gene_matrix_cstb.h5ad')


cluster_df_add.to_csv('snapshot_h5ad/cluster_df_add_cstb.txt',index=True,sep="\t",header=True)



    
#############customize plot cluster_df_add for umap and quality control plotting###############



###plot pretty umap with clusterPlot with cluster_df_add_cstb 

sc.pl.umap(gene_matrix, use_raw=False, color=["leiden"],size=10, legend_loc = 'on data', )

#with shrink/shift and text halo
clusterPlot(cluster_df_add ,
            feature1 = 'cluster',
            feature2=  'cell_type_leiden',
            color_use = color_good, 
            title = ['cluster','cell type'],
            size = 1,
            fontsize = 10,
            width = 13, height = 8, 
            dpi = 100,
            shrink = 2, 
            save = 'pdfs/UMAP/PLA-early-combined-atac.UMAP.pdf'
           )


################plot marker gene/quality control/sample_sex distribution  with sc.pl.umap###############

thresh = {
    'ERVFRD-1':['p30','p100'],
    'PAPPA':['p65','p99'],
    'LAMA3':['p50','p100'],
    'STAT5A':['p30','p99'],
    'ESRRG':['p50','p99'],
    'FLT1':['p60','p100'],
    'ENG':['p60','p100'],
    'FOSL1':['p10','p99'],
    'JUNB':['p30','p99'],
    'FOS':['p30','p99'],
    'HLA-G':['p30','p100'],
    'PLAC8':['p50','p100'],
    
    'DNMT1':['p10','p100'],
    'PSG8':['p10','p100'],
    'CGA':['p10','p100'],
    'SH3TC2':['p60','p100'],
    'LEP':['p10','p100'],
    
    
    'tsse': ['p0','p100'],
    'n_fragment' : ['p0','p100']
}


#marker_genes = ["ERVFRD-1","PAPPA",'LAMA3','STAT5A','ESRRG','FLT1','ENG','FOSL1','JUNB','FOS','HLA-G','PLAC8']
marker_genes = ["DNMT1","ERVFRD-1",'PSG8','CGA','SH3TC2','LEP','PAPPA','FLT1']

marker_genes = ['tsse', 'n_fragment']

marker_genes_sel = []
[ marker_genes_sel.append(i) for i in marker_genes if i in thresh.keys()]


##for plotting
plt.rcParams['figure.figsize'] =  (4,6)
plt.rcParams['legend.loc'] = 'upper left'
plt.rcParams['figure.dpi'] = 100

##for save
#plt.rcParams['savefig.directory'] = "pdfs/marker_genes"
plt.rcParams['savefig.dpi'] = 100
#sc.set_figure_params(fontsize = 15, frameon = False, dpi = 100, )

##get and set umap xlim ylim to make plotting slim
xmin = gene_matrix.obsm['X_umap'][:,0].min()
xmax = gene_matrix.obsm['X_umap'][:,0].max()
ymin = gene_matrix.obsm['X_umap'][:,1].min()
ymax = gene_matrix.obsm['X_umap'][:,1].max()

shrink = 2

##for marker gene/quality control metric
res_p = sc.pl.umap(gene_matrix, 
           use_raw=False, 
           ncols=4,
           color= marker_genes_sel,
           vmin = [ thresh[i][0]  for i in marker_genes_sel],
           vmax= [ thresh[i][1]  for i in marker_genes_sel],
           size= 10,#15,
           color_map = my_cmap,#my_cmap_tfdev,#my_cmap,
           #colorbar_loc = None, only in scanpy 1.9
           show = False,
           title = [ i+ " cutoff: " + thresh[i][0] + " - " + thresh[i][1] for i in marker_genes_sel],
           return_fig = False, #return as a whole fig for save?
           legend_loc = 'left margin' #unwork

          )

##for sample, sex distribution
res_p = sc.pl.umap(gene_matrix, 
           #use_raw=False, 
           #ncols=4,
           color= ['sample','sex'],
           #vmin = [ thresh[i][0]  for i in marker_genes_sel],
           #vmax= [ thresh[i][1]  for i in marker_genes_sel],
           size= 10,#15,
           #color_map = my_cmap_tfdev,#my_cmap,
           palette = color_good,
           #colorbar_loc = None, only in scanpy 1.9
           show = False,
           title = ['sample','sex'],
           return_fig = False, #return as a whole fig for save?
           legend_loc = 'right margin',
           #legend_fontsize = 20,
           #legend_fontweight = 10,
           #legend_loc = 'left margin' #unwork
           #legend_loc = ['left margin','right margin','right margin','right margin','right margin','right margin','right margin','right margin','right margin','none','right margin','right margin'] #not work, do not know why
          )

for i in range(len(res_p)):
    res_p[i].set_xlim(shrink*xmin,shrink*xmax)
    res_p[i].set_ylim(shrink*ymin,shrink*ymax)

#plt.suptitle('Marker gene activity score',y = 1.0, x = 0.45)
#plt.suptitle('Quality control metric',y = 1.0, x = 0.45)
#plt.tight_layout()

###save pdf ####

plt.savefig("pdfs/marker_genes/marker_genes.pdf",format='pdf',dpi = 100)
#plt.savefig("pdfs/QC/quality_control.pdf",format='pdf',dpi = 100)
#plt.savefig("pdfs/QC/sample_sex_distribution.pdf",format='pdf',dpi = 100)


######split plotting with sample and sex (by subset loop or by seaborn faceGrid )####
color_sex = {'female': 'magent','male': 'darkblue'}

# for i in cluster_df_add['sex'].cat.categories:

#     cluster_df_add_sel = cluster_df_add[cluster_df_add['sex'] == i]
#     cluster_df_add_sel['sex'] = cluster_df_add_sel['sex'].cat.remove_unused_categories()
#     clusterPlot(cluster_df_add_sel ,
#                 feature1 = 'cluster',
#                 feature2=  'sex',
#                 color_use = color_sex[i], 
#                 title = ['cluster','sex'],
#                 size = 5,
#                 fontsize = 15,
#                 width = 13, height = 8, 
#                 dpi = 100,
#                 shrink = 1.2, 
#                 save = None#'pdfs/UMAP/PLA-early-combined-atac.UMAP.pdf'
#                )

##########seaborn FacetGrid scatter for cluster_df_add column######
#######https://cduvallet.github.io/posts/2018/11/facetgrid-ylabel-access##


df = cluster_df_add

columnid = 'sex'
#columnid = 'sample'

##get and set umap xlim ylim to make plotting slim
xmin = df['UMAP-1'].min()
xmax = df['UMAP-1'].max()
ymin = df['UMAP-2'].min()
ymax = df['UMAP-2'].max()

shrink = 1.5

# plt.rcParams['figure.figsize'] =  (4,6)
# plt.rcParams['legend.loc'] = 'upper left'
# plt.rcParams['figure.dpi'] = 100

sns.set_style(style = 'white')

g = sns.FacetGrid(df, 
                  col = columnid,#'sex',
                  height = 3.5, 
                  aspect = 0.65, 
                  palette = color_sex, 
                  #hue_kws = {'size':1},
                  #sharex = True, 
                  #sharey = True,
                  #xlim = (shrink*xmin,shrink*xmax),
                  #ylim = (shrink*ymin,shrink*ymax),
                  #subplot_kws = { 'grid':False   }
                  
                 )

kws = {'s':.1, 'linewidth' : .5, 'edgecolor' :  None}
g.map_dataframe(sns.scatterplot,x="UMAP-1",y="UMAP-2",hue = columnid, **kws)

g.despine(left=True, bottom=True, right = True,top = True)

g.set_xlabels('UMAP-1',fontsize = 10 )
g.set_ylabels('UMAP-2',fontsize = 10 )

g.set( xlim=(shrink*xmin,shrink*xmax),ylim=(shrink*ymin,shrink*ymax),
       xticklabels = [], yticklabels = []
     )

g.add_legend() 
#g.tight_layout()

#
g.savefig("pdfs/QC/sex_distribution.split.pdf",format='pdf',dpi = 100)
#g.savefig("pdfs/QC/sample_distribution.split.pdf",format='pdf',dpi = 100)



######call peaks with data_cstb with strict pval == 0.01 by snapatac2 method (only for summit)######

#snap.tl.call_peaks(data_cstb, groupby="cluster",q_value = 0.01, out_dir = 'peaks_cstb')
#snap.tl.call_peaks(data_cstb, groupby="cluster",q_value = 0.001, out_dir = 'peaks_cstb_q0.001')

snap.tl.call_peaks(data_cstb, groupby="cluster",q_value = 0.01, out_dir = 'MACS_snapatac2_use_no_lambda')


peak_mat_cstb = snap.pp.make_peak_matrix(data_cstb)
AnnData object with n_obs × n_vars = 22815 × 404198
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster'

peak_mat_cstb_strict = snap.pp.make_peak_matrix(data_cstb)
AnnData object with n_obs × n_vars = 22815 × 368809
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden'


peak_mat_cstb.write('peak_matrix_cstb.h5ad')
peak_mat_cstb_strict.write('peak_matrix_cstb_q0.001.h5ad')


##output combined peaks 

with open('peak.combine.cstb.bed','w') as fh:
    for i in peak_mat_cstb.var_names:
        (chr,range) = i.split(":")
        (start, end) = range.split("-")
        fh.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + i + "\n")
        

with open('peak.combine.cstb.q0.001.bed','w') as fh:
    for i in peak_mat_cstb_strict.var_names:
        (chr,range) = i.split(":")
        (start, end) = range.split("-")
        fh.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + i + "\n")


#snap.ex.export_bigwig(data_cstb,out_dir = "./bw_cstb_q0.001", prefix = "cluster",groupby = 'cluster' )     



######call peak by standalone macs2 with SnapATAC v1 parameters and make_peak_matrix manually#######

##snapatac2 call peak parameters (use --call-summits then link afterwards)
#/home/mjwang/anaconda3/envs/r413/bin/macs2 #v2.2.7.1
#-shift -100 -extsize 200 -nomodel -callsummits -nolambda -keep-dup all

#1  code in snapatac2-core/src/export.rs
"callpeak",
        "-f", "BED",
        "-t", bed_file.as_ref().to_str().unwrap(),
        "--keep-dup", "all",
        "--outdir", format!("{}", dir.path().display()).as_str(),
        "--qvalue", format!("{}", q_value).as_str(),
        "-g", format!("{}", (genome_size as f64 * 0.9).round()).as_str(), genome_size = self.read_chrom_sizes()?.into_iter().map(|(_, v)| v).sum()
        "--call-summits", #??
        "--nomodel", "--shift", "-100", "--extsize", "200",
        "--nolambda", ##?
        "--tempdir", format!("{}", dir.path().display()).as_str(),

##2 snapatac call peak parameters
--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR
# -B --bdg save extended fragment pileup
# --SPMR SAVE signal per million reads for fragment pileup profiles

"callpeak", 
			   "-t", combined.bed, 
			   "-f", "BED",
			   "-g", gsize,
			   macs.options,
			   "-n", output.prefix


#use Tn5 insertions in pl.call_peak

##or output Tn5 insertion bed 
#snap.ex.export_bed(data_cstb,out_dir = "./insertions_cstb", prefix = "insertion_cluster",groupby = 'cluster', ids = False,suffix='insertion' )  #the same with call_peak output insertions, 1-8 cluster checked

#for i in {1..8};do echo $i; diff -s <(zcat MACS_cstb/${i}_insertion.bed.gz) <(zcat insertions_cstb/cluster${i}.bed.gz |cut -f 1-3 ); done #all identical


#run macs2 in MACS_cstb dir

peaks_bed = pd.read_csv('MACS_cstb/peaks.combined.bed',sep='\t',header=None) #merged peak
#274189, checked in igv

peaks_bed_join = peaks_bed.iloc[:,0] + ":" + peaks_bed.iloc[:,1].astype('string') + "-" + peaks_bed.iloc[:,2].astype('string')

peaks_bed_join_list = peaks_bed_join.to_list()

peak_mat_cstb_manual = snap.pp.make_peak_matrix(adata=data_cstb,  #will get obsm insertion to count, see fn raw_count_iter
                                                use_rep= peaks_bed_join_list, 
                                                #peak_file='peaks_cstb_q0.001/peaks.combined.bed'
                                               ) 

AnnData object with n_obs × n_vars = 22815 × 274189
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden'


(peaks_bed_join_list == peak_mat_cstb_manual.var_names).all() #True

#peak_mat_cstb_manual.write('peak_mat_cstb_manual.h5ad')
#peak_mat_cstb_manual = sc.read("peak_matrix_cstb.h5ad")


# with open('peak.combine.cstb.manual.bed','w') as fh:
#     for i in peak_mat_cstb_manual.var_names:
#         (chr,range) = i.split(":")
#         (start, end) = range.split("-")
#         fh.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + i + "\n")




###output bw track file
snap.ex.export_bigwig(data_cstb,out_dir = "./bw_cstb", prefix = "cluster",groupby = 'cluster' )     



####save pmat raw or processed to mm format for r##

peak_mat = peak_mat_cstb_manual


peak_mat.X.shape #22815 x 274189 sparse matrix of numpy.float64, compressed sparse row format
peak_mat.raw.X.shape #22815 x 274189 sparse matrix of numpy.float32,Compressed Sparse Row format

mem()
#Total memory usage: 193.1791 GB

from datetime import datetime

timetag1 = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
print(timetag1)

mmwrite('toRDS/pmat.mtx',peak_mat.X,precision=3) #5.1G, sparse matrix, float, no need to use coo_matrix to transform as coo format
mmwrite('toRDS/pmat.raw.mtx',peak_mat.raw.X,field = 'integer') #sparse matrix integer, or will assign unsigned integer


timetag2 = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
print(timetag2) #quick, 10 min

mem()

##save gmat colname and rowname
with open('toRDS/pmat.rowname.txt','w') as fh:
    for i in peak_mat.obs_names.to_list():
        fh.write(i + "\n")

with open('toRDS/pmat.colname.txt','w') as fh:
    for i in peak_mat.var_names.to_list():
        fh.write(i + "\n")

        
###add peak logic table to uns####

peak_logic_df = pd.read_csv('MACS_cstb/peak.combine.logic.txt',sep='\t')

peak_logic_df['Peaks'] = peak_logic_df.index
peak_logic_df = peak_logic_df.loc[:,['Peaks','1','2','3','4','5','6','7','8','9']]


##substitute uns['peaks']

data_cstb.uns['peaks'] = peak_logic_df
        

    
##filter for removal of c9 30 cell



all(data_cstb.obs_names ==  cluster_df_add.index) #True

all(cluster_df_add.index.isin(peak_mat.obs_names))

peak_mat_filter_c9 = peak_mat[cluster_df_add.index,:]

all(peak_mat_filter_c9.obs_names == data_cstb.obs_names) #True


peak_mat = peak_mat_filter_c9

del peak_mat_filter_c9




###calculate summit_mat

summits_bed = pd.read_csv('MACS_cstb/summits.combined.merge100.extend200.bed',sep='\t',header=None) 
#387932, first merge within 100bp summit then extend to 200bp if length < 200bp
#263723, first merge within 500bp summit then extend to 200bp if length < 200bp
#1333735 , 200bp extend from highest peak point (-100, +100)

summits_bed_join = summits_bed.iloc[:,0] + ":" + summits_bed.iloc[:,1].astype('string') + "-" + summits_bed.iloc[:,2].astype('string')

summits_bed_join_list = summits_bed_join.to_list()

summit_mat_cstb_manual = snap.pp.make_peak_matrix(adata=data_cstb, 
                                                use_rep= summits_bed_join_list, 
                                                #peak_file='peaks_cstb_q0.001/peaks.combined.bed'
                                               ) 

AnnData object with n_obs × n_vars = 22785 × 387932
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
    
# AnnData object with n_obs × n_vars = 22815 × 274189
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden'


(summits_bed_join_list == summit_mat_cstb_manual.var_names).all() #True


summit_mat = summit_mat_cstb_manual

del summit_mat_cstb_manual


summit_mat.write('snapshot_h5ad/summit_mat_cstb.h5ad')



####save smat raw or processed to mm format for r##

summit_mat.X.shape #22785, 387932  sparse matrix of numpy.float64, compressed sparse row format
#summit_mat.raw.X.shape 

mem()
#Total memory usage: 108.8138 GB

from datetime import datetime

timetag1 = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
print(timetag1)

mmwrite('toRDS/smat.mtx',summit_mat.X,field = 'integer') #5.1G, sparse matrix, float, no need to use coo_matrix to transform as coo format
#mmwrite('toRDS/smat.raw.mtx',summit_mat.raw.X,field = 'integer') #sparse matrix integer, or will assign unsigned integer


timetag2 = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
print(timetag2) #quick, 10 min

mem()

##save gmat colname and rowname
with open('toRDS/smat.rowname.txt','w') as fh:
    for i in summit_mat.obs_names.to_list():
        fh.write(i + "\n")

with open('toRDS/smat.colname.txt','w') as fh:
    for i in summit_mat.var_names.to_list():
        fh.write(i + "\n")



#############Identify differentially accessible regions (DAR) by a quick and dirty method (aggregate signal across cells and utilizes z-scores to identify specifically enriched peaks) #############


#data = snap.read(str(snap.datasets.pbmc5k(type="annotated_h5ad")))

#data = snap.read('data.snapatac2.h5ad')

    
#data
#5977 x 6176550
#4362 x 6176550 


##marker_peaks = snap.tl.marker_regions(peak_mat, groupby='leiden', pvalue=0.01) #will cause problem
##marker_peaks = snap.tl.marker_regions(peak_mat, groupby='cluster', pvalue=0.01)
##marker_peaks = snap.tl.marker_regions(peak_mat_cstb_manual, groupby='cluster', pvalue=0.01)





#########do marker peaks step by step with snapatac v2.2.0 and scanpy rank_gene method##########

#count = pl.DataFrame(aggregate_X(peak_mat, groupby='leiden', normalize="RPKM"))

(peak_mat.obs_names == data_cstb.obs_names ).all()#True

peak_mat.obsm['X_umap'] = data_cstb.obsm['X_umap_rotate']



peak_mat.obs['cluster'].value_counts().sort_index()
1    3937
2    3594
3    2921
4    2493
5    2230
6    3713
7    1636
8    1447
9     844

peak_mat_aggre = aggregate_X(peak_mat, groupby='cluster', normalize="RPKM")
9 × 274189
10 × 483312

peak_mat_aggre.obs_names
#'1', '2', '3', '4', '5', '6', '7', '8', '9'

peak_mat_aggre.obs_names = 'c' + peak_mat_aggre.obs_names


peak_mat_aggre_df = pd.DataFrame(peak_mat_aggre.X,index=peak_mat_aggre.obs_names,columns=peak_mat_aggre.var_names)#.T

peak_mat_aggre_df = peak_mat_aggre_df.loc[['c1','c7','c9','c2','c8','c3','c6','c4','c5'],:]


count = pl.DataFrame(peak_mat_aggre_df)
#aggre_group x peak  , 10 x 483312


count = count[:,['c1','c7','c9','c2','c8','c3','c6','c4','c5']]


#######snapatac v2.2.0 method####

names = np.array(peak_mat.var_names)

#names = np.array(peak_mat_aggre_df.index)

#len(names) == count.shape[1]# True

#first zscore data
z = scipy.stats.zscore(
    np.log2(1 + count.to_numpy()),
    axis = 1, #by column
)
#10 x 483312, np.array

#then test by column
peaks = {}
pvalue=0.01
del range
for i in range(z.shape[1]): #test by each column (peaks)
    pvals = scipy.stats.norm.sf(z[:, i])
    select = pvals < pvalue
    if np.where(select)[0].size >= 1:
        peaks[count.columns[i]] = names[select]

with open('marker_peaks.pkl','wb') as fh:
    dill.dump(marker_peaks,file=fh)


marker_peaks = peaks

for i in marker_peaks.keys():
    print("length of cluster %s is %i" % (i,len(marker_peaks[i])) )

length of cluster 1 is 4789
length of cluster 2 is 239
length of cluster 3 is 25
length of cluster 4 is 137
length of cluster 5 is 3180
length of cluster 6 is 13
length of cluster 7 is 1409
length of cluster 8 is 47
length of cluster 9 is 3729


marker_peaks_merge = np.concatenate(list(marker_peaks.values())).tolist()
#13568

marker_peaks_sel = marker_peaks['5'] #3180
marker_peaks_sel = marker_peaks['3'] #25
marker_peaks_sel = marker_peaks['2']

marker_peaks_sel = np.concatenate([marker_peaks['2'],marker_peaks['5']])

marker_peaks_sel = marker_peaks_merge



    
#######use scanpy difftest method at single cell level (use this)#########

##do diff test for raw peak count (no logarithemize)

peak_mat.X.todense()[1:1000,1:1000]
peak_mat.raw.X.todense()[1:1000,1:1000]
peak_mat.raw = peak_mat

sc.tl.rank_genes_groups(peak_mat, 'cluster', method='t-test', key_added = "t-test")
sc.pl.rank_genes_groups(peak_mat, n_genes=25, sharey=False, key = "t-test")


#peak_mat.uns['t-test'].keys()
#peak_mat.uns['t-test']['names']


##do diff test after logarithemize and normalize

#sc.pp.normalize_per_cell(peak_mat, counts_per_cell_after=1e1) #will be error
sc.pp.log1p(peak_mat)
sc.tl.rank_genes_groups(peak_mat, 'cluster', method='t-test', key_added = "t-test-log")
sc.pl.rank_genes_groups(peak_mat, n_genes=25, sharey=False, key = "t-test-log")



##extract to df
marker_peaks_df = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 't-test',pval_cutoff = 0.001, log2fc_min = 1)
marker_peaks_df = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 't-test-log',pval_cutoff = 0.001, log2fc_min = 1)
marker_peaks_df.group.value_counts().sort_index()
1    62380
2     6341
3     4342
4     2401
5    40169
6       50
7    25678
8      158
9     5541

(marker_peaks_df.pvals <= 0.001).all() #True



##save dar full table
marker_peaks_df.to_csv('DARs_ttest/marker_peaks.scanpy_rank.ttest.txt',sep='\t',index=False)


##save as dict
marker_peaks_df.group
#Categories (9, object): ['1', '2', '3', '4', ..., '6', '7', '8', '9']

marker_peaks_dict = marker_peaks_df.groupby('group').apply(lambda x: x.names.to_list()).to_dict() #apply to dataframe, select series to dict
type(marker_peaks_dict)
marker_peaks_dict.keys()
['1', '2', '3', '4', '5', '6', '7', '8', '9']

marker_peaks_dict['3'][0:10]
marker_peaks_dict['3'][-10:]

for i in marker_peaks_dict.keys():
    print(i + ': ' + str(len(marker_peaks_dict[i])))
1: 62380
2: 6341
3: 4342
4: 2401
5: 40169
6: 50
7: 25678
8: 158
9: 5541


##select top 1000 marker peak to plot

#marker_peaks_sel_df = marker_peaks_df.groupby('group').head( 1000 )
marker_peaks_sel_df = marker_peaks_df.groupby('group').apply(lambda x: x.sort_values(['logfoldchanges'], ascending = False).head(1000)).reset_index(drop = True) #sort by log2foldchange top1000 use this?

marker_peaks_sel_df.group.value_counts().sort_index()
1    1000
2    1000
3    1000
4    1000
5    1000
6      50
7    1000
8     158
9    1000

marker_peaks_sel_df.to_csv('DARs_ttest/marker_peaks.scanpy_rank.ttest.top1000.sortby_logfc.txt',sep='\t',index=False)

marker_peaks_sel = marker_peaks_sel_df.names



marker_peaks_sel_dict = marker_peaks_sel_df.groupby('group').apply(lambda x: x.names.to_list()).to_dict()

for i in marker_peaks_sel_dict.keys():
    print(i + ': ' + str(len(marker_peaks_sel_dict[i])))
1: 1000
2: 1000
3: 1000
4: 1000
5: 1000
6: 50
7: 1000
8: 158
9: 1000

#marker_peaks_df.groupby(['group'])['logfoldchanges'].nlargest(100)


# marker_peaks_merge = marker_peaks_merge_df.names
# marker_peaks_sel = marker_peaks_merge





#########visuzalize marker peak ##############


# def marker_regions(
#     data: AnnData | AnnDataSet,
#     groupby: str | list[str],
#     pvalue: float = 0.01,
# ) -> dict[str, list[str]]:
#     """
#     A quick-and-dirty way to get marker regions.
#     """
#     import scipy.stats
#     import polars as pl

#     count = pl.DataFrame(aggregate_X(data, groupby, normalize="RPKM"))
#     names = np.array(data.var_names)
#     z = scipy.stats.zscore(
#         np.log2(1 + count.to_numpy()),
#         axis = 1,
#     )
#     peaks = {}
#     for i in range(z.shape[1]):
#         pvals = scipy.stats.norm.sf(z[:, i])
#         select = pvals < pvalue
#         if np.where(select)[0].size >= 1:
#             peaks[count.columns[i]] = names[select]
#     return peaks



####plot dar peaks

##snap.pl.regions(peak_mat, groupby='cluster', peaks=marker_peaks, interactive=False) #cause problem


#count = pl.DataFrame(aggregate_X(data, groupby=groupby, normalize="RPKM"))


peakid = peak_mat.var_names.to_list()

# idx = peakid.index('chr1:8396853-8397361')# marker_peaks_merge)
# peakid[idx]

# idx = peakid.index(['chr1:8396853-8397361','chr1:8392538-8392940']) #wrong for list index matching, only one element for one time!
# peakid[idx]


# idx = []
# for i,v in enumerate(marker_peaks_sel):#slow but ok
# #    if v in marker_peaks_c5:
#     if v in peakid:
#         idx.append(i)

# idx = []
# for i in marker_peaks_sel:
#     for j,v in enumerate(peakid):
#         if v == i:
#             idx.append(j)


idx = []#quick! use this!
for i in marker_peaks_sel:
    if i in peakid:
        idx.append(peakid.index(i))#list index will only return the first one match idx
            

            

mat = np.log2(1 + count.to_numpy()[idx, :])


mat.shape
13568 x 9
239 x 9 #c2
3180 x  9 #c5 marker peak
13568 x 9 #all marker peak


##calculate zscore for dataframe by column, borrow from scenic code

# mat_Z = pd.DataFrame( index=peak_mat_aggre_df.index ) #create an empty DataFrame

# for col in list(peak_mat_aggre_df.columns):#very slow
#     peak_mat_aggre_df_Z[ col ] = ( peak_mat_aggre_df[col] - peak_mat_aggre_df[col].mean()) / peak_mat_aggre_df[col].std(ddof=0)


# z = scipy.stats.zscore(
#     np.log2(1 + count.to_numpy()),
#     axis = 1, #by column
# )

mat_Z = zscore(mat,axis = 1,ddof = 0 )




# for col in list(peak_mat_aggre_df.columns):#very slow
#     peak_mat_aggre_df_Z[ col ] = zscore(peak_mat_aggre_df[ col ])




trace = go.Heatmap(
    x=count.columns,
    y=[peakid[i] for i in idx], #np.concatenate(list(peaks.values()))[::-1],
    z=mat_Z,
    zmin = -1,
    zmax = 1,
    zmid=0,
    type='heatmap',
    colorscale='Blues',#'Reds',#'RdBu_r',#'YlGnBu',#'Viridis',
    colorbar={ "title": "zscore log2(1 + RPKM)" },
)
#data = [trace]
layout = {
    "yaxis": { "visible": False, "autorange": "reversed" },
    "xaxis": { "title":'cluster' },
}

fig = go.Figure(data=[trace], layout=layout)


fig.update_layout(
    template="simple_white",
    legend= {'itemsizing': 'constant'},
)


fig.show(width=500,height=450)

#render_plot(fig, width, height, interactive, show, out_file)

#return render_plot(fig, width, height, interactive, show, out_file)






    
# def regions(
#     data: AnnData | AnnDataSet,
#     groupby: str | list[str],
#     peaks: dict[str, list[str]],
#     width: float = 600,
#     height: float = 400,
#     show: bool = True,
#     interactive: bool = True,
#     out_file: str | None = None,
# ) -> 'plotly.graph_objects.Figure' | None:
#     """
#     Parameters
#     ----------
#     data
#         Annotated data matrix.
#     groupby
#         Group the cells into different groups. If a `str`, groups are obtained from
#         `.obs[groupby]`.
#     peaks
#         Peaks of each group.
#     width
#         The width of the plot
#     height
#         The height of the plot
#     show
#         Show the figure
#     interactive
#         Whether to make interactive plot
#     out_file
#         Path of the output file for saving the output image, end with
#         '.svg' or '.pdf' or '.png' or '.html'.

#     Returns
#     -------
#     'plotly.graph_objects.Figure' | None
#         If `show=False` and `out_file=None`, an `plotly.graph_objects.Figure` will be 
#         returned, which can then be further customized using the plotly API.
#     """
#     import polars as pl
#     import plotly.graph_objects as go

#     count = pl.DataFrame(aggregate_X(data, groupby=groupby, normalize="RPKM"))
#     idx = data.var_ix(np.concatenate(list(peaks.values())).tolist())
#     mat = np.log2(1 + count.to_numpy()[idx, :])

#     trace = go.Heatmap(
#         x=count.columns,
#         y=np.concatenate(list(peaks.values()))[::-1],
#         z=mat,
#         type='heatmap',
#         colorscale='Viridis',
#         colorbar={ "title": "log2(1 + RPKM)" },
#     )
#     data = [trace]
#     layout = {
#         "yaxis": { "visible": False, "autorange": "reversed" },
#         "xaxis": { "title": groupby },
#     }
#     fig = go.Figure(data=data, layout=layout)
#     return render_plot(fig, width, height, interactive, show, out_file)




####visulize maker peaks in one umap####

#https://github.com/scverse/scanpy/issues/532
#...That wouldn't give you the sum of the gene expression values, but the average expression minus an average expression of a random gene set.

#for i in marker_peaks_sel_dict.keys():
##for i in marker_peaks_dict.keys():
for i in ['1','7','9','2','8','3','6','4','5']:
    print('do dar visualization for c%s ' % i)
    
    #sc.tl.score_genes(peak_mat, marker_peaks_sel_dict[i], score_name='c'+i+'_dar_score') #method1 use score_genes
    sc.tl.score_genes(peak_mat, marker_peaks_dict[i], score_name='c'+i+'_dar_score') #method1 use score_genes
    #peak_mat.obs['c3_dar_score_new'] = peak_mat.X[:,marker_peaks_sel_dict['3']].sum(1) #method2 simple sum

    plt.rcParams['figure.figsize'] =  (5,5)
    #plt.rcParams['font.size'] = 45
    sc.set_figure_params(fontsize = 10, frameon = False, dpi = 150, )

    sc.pl.umap(peak_mat, use_raw=False, 
               color=['c'+i+'_dar_score'],
               #color=["cluster",'c'+i+'_dar_score'],
               vmin='p30',
               vmax='p99',
               size=10,
               color_map = my_cmap,
               title = ['c'+i+'_dar_score'+'(number of DARs: ' + str(len(marker_peaks_dict[i])) + '  )'],
               #title = ['cluster','c'+i+'_dar_score'+'(number of DARs:' + str(len(marker_peaks_sel_dict[i])) + '  )'],
               show = True# False
            )


#peak_mat.write('peak_mat_cstb_manual.h5ad')



#######see tf motif enrichment in these dars(Homer enrichment)####
motifs = snap.tl.motif_enrichment( 
    motifs=snap.datasets.cis_bp(unique=True),
    regions=marker_peaks_dict,#116001 also quick #marker_peaks_sel_dict quick for 7150 peaks
    genome_fasta=snap.genome.hg38,
)

snap.pl.motif_enrichment(motifs, min_log_fc= 2, max_fdr=0.0001, height=2200, interactive=False)

snap.pl.motif_enrichment(motifs, min_log_fc= 5, max_fdr=0.0001, height=5200, interactive=True)





##Identify differentially accessible regions by regression-based differential test, taking variations across cells into consideration (single cell level) ##use this



#group1 = "Naive/Intermediate B"
#group2 = "Memory B"

group1 = '3'
group2 = '5'



##select peak
peaks_selected = np.logical_or( #274189
    data_cstb.uns["peaks"][group1].to_numpy(),
    data_cstb.uns["peaks"][group2].to_numpy(),
)

peaks_selected.sum() #192182 of 274189, only group1 and group2 specific peak

##select cell
flag1 = data_cstb.obs['cluster'] == group1
flag2 = data_cstb.obs['cluster'] == group2

cell_1 = data_cstb.obs_names[flag1].to_list() #2921
cell_2 = data_cstb.obs_names[flag2].to_list() #2230


idx = [ i  for i in range(0,len(data_cstb.obs['cluster']) ) if data_cstb.obs['cluster'].to_list()[i] == group1  ]
#cell1 = data_cstb.obs['Cell'][idx].to_list()
cell1 = data_cstb.obs_names[idx].to_list()

idx = [ i  for i in range(0,len(data_cstb.obs['cluster']) ) if data_cstb.obs['cluster'].to_list()[i] == group2  ]
#cell2 = data_cstb.obs['Cell'][idx].to_list()
cell2 = data_cstb.obs_names[idx].to_list()

(peak_mat.obs_names == data_cstb.obs_names).all() #True
cell_1 == cell1 #True
cell_2 == cell2 #True

diff_peaks = snap.tl.diff_test(
    peak_mat,
    cell_group1= cell1, #naive_B,
    cell_group2= cell2 ,#memory_B,
    features=peaks_selected,
)

##min_log_fc = 0.25, min_pct = 0.05

diff_peaks.head() #polars
45379 x 4, 45379 of 192182

diff_peaks.write_csv("DARs_pair_regression/diff_peaks_c5c3.txt",sep='\t')


# snap.pl.regions(
#     peak_mat,
#     groupby = 'leiden',#'cell_type',
#     peaks = {
#         group1: diff_peaks.filter(pl.col("log2(fold_change)") > 0)['feature name'].to_numpy(),
#         group2: diff_peaks.filter(pl.col("log2(fold_change)") < 0)['feature name'].to_numpy(),
#     },
#     interactive = False,
# )


##########rank gene by scanpy  peak_mat_c5c3##########


peak_mat_c5c3 = peak_mat[]

sc.tl.rank_genes_groups(peak_mat_c5c3, 'cluster', method='t-test', key_added = "t-test-log")
sc.pl.rank_genes_groups(peak_mat_c5c3, n_genes=25, sharey=False, key = "t-test-log")



##extract to df
marker_peaks_df_c5c3 = sc.get.rank_genes_groups_df(peak_mat_c5c3,group = None, key = 't-test',pval_cutoff = 0.001, log2fc_min = 1)
marker_peaks_df_c5c3 = sc.get.rank_genes_groups_df(peak_mat_c5c3,group = None, key = 't-test-log',pval_cutoff = 0.001, log2fc_min = 1)
marker_peaks_df_c5c3.group.value_counts().sort_index()
1    62380
2     6341
3     4342
4     2401
5    40169
6       50
7    25678
8      158
9     5541

(marker_peaks_df_c5c3.pvals <= 0.001).all() #True



##save dar full table
marker_peaks_df_c5c3.to_csv('DARs_ttest/marker_peaks_c5c3.scanpy_rank.ttest.txt',sep='\t',index=False)









# ##diff test one-with-rest


# barcodes = np.array(data.obs_names)
# background = []
# for i in np.unique(data.obs['leiden']):
#     if i != group2:
        
#         idx = [ j  for j in range(0,len(data.obs['leiden']) ) if data.obs['leiden'].to_list()[j] == i  ]
#         cells = barcodes[idx]
#         if np.unique(data.obs['leiden'][idx]) != i:
#             #print('wrong with leiden subsetting')
#             raise SystemExit('wrong with leiden subsetting')
#         cells = np.random.choice(
#             #barcodes[data.obs['leiden'] == i],
#             cells,
#             size = 50, #gather randomly 50 cells as background
#             replace = False,
#         )
#         background.append(cells)
        
# background = np.concatenate(background)


# diff_peaks = snap.tl.diff_test(
#     peak_mat,
#     cell_group1 = cell1,#memory_B,
#     cell_group2 = background,
#     features = data.uns["peaks"]['0'].to_numpy(),
#     direction = "positive",
# )


# snap.pl.regions(
#     peak_mat,
#     groupby = 'leiden',
#     peaks = {
#         group2: diff_peaks['feature name'].to_numpy(),
#     },
#     interactive = False,
# )

# #27410 dar fo cluster 0


##see tf motif enrichment in these dars
motifs = snap.tl.motif_enrichment(
    motifs=snap.datasets.cis_bp(unique=True),
    regions=marker_peaks,
    genome_fasta=snap.genome.hg38,
)

snap.pl.motif_enrichment(motifs, max_fdr=0.0001, height=1200, interactive=False)






#############save and reload obj#####

data
AnnData object with n_obs × n_vars = 33786 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'selected'
    uns: 'reference_sequences', 'spectral_eigenvalue', 'AnnDataSet'
    obsm: 'X_spectral_mnn', 'X_spectral', 'X_umap', 'X_spectral_harmony', 'insertion'
    obsp: 'distances'
    
# AnnData object with n_obs × n_vars = 33786 × 6176550
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library'
#     var: 'selected'
#     uns: 'AnnDataSet', 'spectral_eigenvalue', 'reference_sequences'
#     obsm: 'insertion', 'X_spectral', 'X_umap', 'X_spectral_mnn', 'X_spectral_harmony'
#     obsp: 'distances'
    
# AnnData object with n_obs × n_vars = 7788 × 6176550
#     obs: 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden'
#     var: 'selected'
#     uns: 'reference_sequences', 'scrublet_sim_doublet_score', 'scrublet_threshold', 'spectral_eigenvalue'
#     obsm: 'insertion', 'X_spectral', 'X_umap'
#     obsp: 'distances'
    
    
# AnnData object with n_obs × n_vars = 5334 × 6176550
#     obs: 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden'
#     var: 'selected'
#     uns: 'reference_sequences', 'scrublet_sim_doublet_score', 'scrublet_threshold', 'spectral_eigenvalue'
#     obsm: 'insertion', 'X_spectral', 'X_umap'
#     obsp: 'distances'
    

    
#data.write('PLA_atac_early_combined.h5ad') #very quick


##data.write('PLA_atac_early_combined.new.h5ad')
##!mv PLA_atac_early_combined.new.h5ad PLA_atac_early_combined.h5ad


#data.write('PLA_atac_early_combined.filtered_raw_cluster.h5ad')
data = snap.read('PLA_atac_early_combined.filtered_raw_cluster.h5ad')

#output barcode only
# with open('barcode.filtered_rawcluster.txt','w') as fh:
# for i in data.obs_names:
#     fh.write(i + "\n")
#data.close()

data = snap.read('PLA_atac_early_combined.h5ad') #33786 x 6176550


####
#data.write('PLA_atac_early_combined.final.h5ad')
data = snap.read('PLA_atac_early_combined.final.h5ad') #29721 x 6176550


data.write('PLA_atac_early_combined.final.final.h5ad')
data = snap.read('PLA_atac_early_combined.final.final.h5ad')

##to memory
data = data.to_memory() #up to 17.9G (h5ad file 17G)

data

AnnData object with n_obs x n_vars = 27380 x 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'selected'
    uns: 'AnnDataSet', 'reference_sequences', 'spectral_eigenvalue'
    obsm: 'X_spectral', 'insertion', 'X_spectral_harmony', 'X_spectral_mnn', 'X_umap'
    obsp: 'distances'
    
    
AnnData object with n_obs × n_vars = 29721 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'selected'
    uns: 'reference_sequences', 'AnnDataSet', 'spectral_eigenvalue'
    obsm: 'X_spectral_harmony', 'X_spectral', 'X_spectral_mnn', 'X_umap', 'insertion'
    obsp: 'distances'
    

AnnData object with n_obs × n_vars = 33786 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk'
    var: 'selected'
    uns: 'reference_sequences', 'AnnDataSet', 'spectral_eigenvalue'
    obsm: 'X_spectral', 'insertion', 'X_umap', 'X_spectral_mnn', 'X_spectral_harmony'
    obsp: 'distances'


    
##################save/reload a snapshot#######################

######save as h5ad
data.write('snapshot_h5ad/data.h5ad')

data_cstb.write('snapshot_h5ad/data_cstb.h5ad')
gene_matrix.write('snapshot_h5ad/gene_matrix_cstb.h5ad')
peak_mat.write('snapshot_h5ad/peak_mat_cstb.h5ad')
    
cluster_df_add.to_csv('snapshot_h5ad/cluster_df_add.txt',index=True,sep="\t",header=True)




######save as dill pkl (will avoid use snap.read or sc.read_h5ad), in work_doDARs.py


# with open('snapshot_pkl/data_cstb.pkl','wb') as fh:#to avoid use snapatac2 read
#     dill.dump(data_cstb,file=fh)



# with open('snapshot_pkl/cluster_df_add_cstb.pkl','rb') as fh:
#     cluster_df_add = dill.load(file=fh)


# with open('snapshot_pkl/peak_mat.pkl','wb') as fh:
#     dill.dump(peak_mat,file=fh)
    
    
# with open('snapshot_pkl/gene_matrix.pkl','wb') as fh:
#     dill.dump(gene_matrix,file=fh)




# with open('snapshot_pkl/summit_mat_cstb.pkl','wb') as fh:
#     dill.dump(summit_mat,file=fh)



####reload

data = snap.read('snapshot_h5ad/data.h5ad').to_memory()

# AnnData object with n_obs × n_vars = 27380 × 6176550
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster'
#     var: 'selected'
#     uns: 'AnnDataSet', 'spectral_eigenvalue', 'peaks', 'reference_sequences'
#     obsm: 'X_spectral_harmony', 'X_umap', 'X_spectral_mnn', 'X_umap_rotate', 'insertion', 'X_spectral'
#     obsp: 'distances'
    
    
data_cstb = snap.read('snapshot_h5ad/data_cstb.h5ad').to_memory()

AnnData object with n_obs × n_vars = 22785 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
    var: 'selected'
    uns: 'peaks', 'reference_sequences', 'AnnDataSet', 'spectral_eigenvalue'
    obsm: 'X_spectral_harmony', 'insertion', 'X_spectral_mnn', 'X_umap', 'X_spectral', 'X_umap_rotate'
    obsp: 'distances'
    
    
# AnnData object with n_obs × n_vars = 22815 × 6176550
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden'
#     var: 'selected'
#     uns: 'peaks', 'reference_sequences', 'AnnDataSet', 'spectral_eigenvalue'
#     obsm: 'X_umap', 'X_umap_rotate', 'X_spectral', 'X_spectral_harmony', 'insertion', 'X_spectral_mnn'
#     obsp: 'distances'


cluster_df_add = pd.read_csv('snapshot_h5ad/cluster_df_add_cstb.txt',sep='\t',index_col = 0)

(cluster_df_add.index == data_cstb.obs_names).all() #TRUE

cluster_df_add.dtypes
# UMAP-1              float64
# UMAP-2              float64
# cluster               int64
# sample               object
# tsse                float64
# n_fragment            int64
# frac_dup            float64
# frac_mito           float64
# doublet_score       float64
# is_doublet             bool
# leiden                int64
# library              object
# leiden_bk             int64
# cell_type_leiden     object
# sex                  object


cluster_df_add['cluster'] = cluster_df_add['cluster'].astype('str').astype('category')
cluster_df_add['leiden'] = cluster_df_add['leiden'].astype('str').astype('category')
cluster_df_add['cell_type_leiden'] = cluster_df_add['cell_type_leiden'].astype('category')

cluster_df_add['sample'] = cluster_df_add['sample'].astype('category')
cluster_df_add['library'] = cluster_df_add['library'].astype('category')

cluster_df_add['sex'] = cluster_df_add['sex'].astype('category')

cluster_df_add.dtypes
# UMAP-1               float64
# UMAP-2               float64
# cluster             category
# sample              category
# tsse                 float64
# n_fragment             int64
# frac_dup             float64
# frac_mito            float64
# doublet_score        float64
# is_doublet              bool
# leiden              category
# library             category
# leiden_bk              int64
# cell_type_leiden    category
# sex                 category




gene_matrix = snap.read('snapshot_h5ad/gene_matrix_cstb.h5ad').to_memory()

AnnData object with n_obs x n_vars = 22785 x 59265 backed at 'snapshot_h5ad/gene_matrix_cstb.h5ad'
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
    var: 'n_cells'
    uns: 'leiden_colors', 'sex_colors', 'sample_colors', 'log1p', 'cluster_colors'
    obsm: 'X_umap'

gene_matrix.to_memory()

# AnnData object with n_obs × n_vars = 22815 × 59265
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden'
#     var: 'n_cells'
#     uns: 'log1p'
#     obsm: 'X_umap'


peak_mat = snap.read('snapshot_h5ad/peak_mat_cstb.h5ad').to_memory()

##dar in uns load in 'snapshot_pkl/peak_mat_cstb.pkl'

AnnData object with n_obs × n_vars = 22785 × 274189
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'n_counts', 'c3_dar_score', 'c1_dar_score', 'c2_dar_score', 'c4_dar_score', 'c5_dar_score', 'c6_dar_score', 'c7_dar_score', 'c8_dar_score', 'c9_dar_score', 'sex'
    obsm: 'X_umap'
    

# AnnData object with n_obs x n_vars = 22815 x 274189 backed at 'snapshot_h5ad/peak_mat.h5ad'
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'n_counts', 'c3_dar_score', 'c1_dar_score', 'c2_dar_score', 'c4_dar_score', 'c5_dar_score', 'c6_dar_score', 'c7_dar_score', 'c8_dar_score', 'c9_dar_score'
#     uns: 't-test', 'log1p', 'cluster_colors', 't-test-log'
#     obsm: 'X_umap'

##error to process uns, del 
    
#peak_mat.uns = {}

#peak_mat = peak_mat.to_memory()

# AnnData object with n_obs × n_vars = 22815 × 274189
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'n_counts', 'c3_dar_score', 'c1_dar_score', 'c2_dar_score', 'c4_dar_score', 'c5_dar_score', 'c6_dar_score', 'c7_dar_score', 'c8_dar_score', 'c9_dar_score'
#     obsm: 'X_umap'
    



summit_mat = snap.read('snapshot_h5ad/summit_mat_cstb.h5ad').to_memory()  

AnnData object with n_obs × n_vars = 22785 × 387932
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
    
