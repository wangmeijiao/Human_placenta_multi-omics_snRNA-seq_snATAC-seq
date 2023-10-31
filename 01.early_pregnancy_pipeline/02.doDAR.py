
###do dar test and visualization for snapatac2 peak_mat obj (an extend analysis of the standard snapatac2 pipeline)#########



import scanpy as sc

#import snapatac2 as snap
from snapatac2.tools._misc import aggregate_X

# import os
# os.environ['NUMEXPR_MAX_THREADS'] = '1'
# os.environ['NUMEXPR_NUM_THREADS'] = '1'


# ##solve for "BLAS : Program is Terminated. Because you tried to allocate too many memory regions." error
# os.environ['OPENBLAS_NUM_THREADS'] = '1' #important!!
# #os.environ['GOTO_NUM_THREADS'] = '1'
# os.environ['OMP_NUM_THREADS'] = '1'

# os.environ['MKL_DYNAMIC'] = 'FALSE'
# os.environ['MKL_NUM_THREADS'] = '1'




import pandas as pd
import numpy as np

from scipy.stats import zscore
from scipy.sparse import coo_matrix #to change csr to coo then mmwrite will ouput coordinate instead of array

import polars as pl

import dill

from getMem import mem

from scipy.io import mmwrite

import plotly.graph_objects as go

import plotly.io as pio
pio.renderers.default = "png"


from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



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


my_cmap = LinearSegmentedColormap.from_list('expression',palcolor,N=10)
my_cmap




##read in obj and re-calculate peak mat

# data_cstb = snap.read('snapshot_h5ad/data_cstb.h5ad').to_memory()

# AnnData object with n_obs × n_vars = 22785 × 6176550
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
#     var: 'selected'
#     uns: 'peaks', 'reference_sequences', 'spectral_eigenvalue', 'AnnDataSet'
#     obsm: 'X_umap', 'X_spectral_mnn', 'insertion', 'X_umap_rotate', 'X_spectral', 'X_spectral_harmony'
#     obsp: 'distances'


# with open('snapshot_pkl/data_cstb.pkl','wb') as fh:#to avoid use snapatac2 read
#     dill.dump(data_cstb,file=fh)

with open('snapshot_pkl/data_cstb.pkl','rb') as fh:#to avoid use snapatac2 read
    data_cstb = dill.load(file=fh)
    #data_cstb_new = dill.load(file=fh)

data_cstb

AnnData object with n_obs × n_vars = 22785 × 6176550
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
    var: 'selected'
    uns: 'peaks', 'reference_sequences', 'spectral_eigenvalue', 'AnnDataSet'
    obsm: 'X_umap', 'X_spectral_mnn', 'insertion', 'X_umap_rotate', 'X_spectral', 'X_spectral_harmony'
    obsp: 'distances'


#all(data_cstb_new.obs == data_cstb.obs) #True
#all(data_cstb_new.var == data_cstb.var) #True
#all(data_cstb_new.uns['peaks'] == data_cstb.uns['peaks']) #True
#(data_cstb_new.obsm['X_umap_rotate'] == data_cstb.obsm['X_umap_rotate']).all()
#(data_cstb_new.obsp['distances'][1:10000,1:10000].todense() == data_cstb.obsp['distances'][1:10000,1:10000].todense()).all() #True
##all True


# cluster_df_add = pd.read_csv('snapshot_h5ad/cluster_df_add_cstb.txt',sep='\t',index_col = 0)

# (cluster_df_add.index == data_cstb.obs_names).all() #TRUE

# cluster_df_add.dtypes
# # UMAP-1              float64
# # UMAP-2              float64
# # cluster               int64
# # sample               object
# # tsse                float64
# # n_fragment            int64
# # frac_dup            float64
# # frac_mito           float64
# # doublet_score       float64
# # is_doublet             bool
# # leiden                int64
# # library              object
# # leiden_bk             int64
# # cell_type_leiden     object
# # sex                  object


# cluster_df_add['cluster'] = cluster_df_add['cluster'].astype('str').astype('category')
# cluster_df_add['leiden'] = cluster_df_add['leiden'].astype('str').astype('category')
# cluster_df_add['cell_type_leiden'] = cluster_df_add['cell_type_leiden'].astype('category')

# cluster_df_add['sample'] = cluster_df_add['sample'].astype('category')
# cluster_df_add['library'] = cluster_df_add['library'].astype('category')

# cluster_df_add['sex'] = cluster_df_add['sex'].astype('category')

# cluster_df_add.dtypes


# with open('snapshot_pkl/cluster_df_add_cstb.pkl','wb') as fh:
#     dill.dump(cluster_df_add,file=fh)


with open('snapshot_pkl/cluster_df_add_cstb.pkl','rb') as fh:
    cluster_df_add = dill.load(file=fh)

cluster_df_add
#22785 x 15


# with open('snapshot_pkl/gene_matrix.pkl','rb') as fh:#to avoid use snapatac2 read
#     gene_matrix = dill.load(file=fh)

# gene_matrix
# AnnData object with n_obs × n_vars = 22785 × 59265
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
#     var: 'n_cells'
#     uns: 'log1p', 'leiden_colors', 'cluster_colors', 'sample_colors', 'sex_colors'
#     obsm: 'X_umap'


# gene_matrix = snap.read('snapshot_h5ad/gene_matrix_cstb.h5ad')#.to_memory()

# gene_matrix.obs
# gene_matrix.var.dtypes

# gene_matrix = gene_matrix.to_memory()

# with open('snapshot_pkl/gene_matrix.pkl','wb') as fh:
#     dill.dump(gene_matrix,file=fh)
    

########recalculate peak mat (keep raw count to raw.X)############
peaks_bed = pd.read_csv('MACS_cstb/peaks.combined.bed',sep='\t',header=None) #merged peak
#274189, checked in igv

peaks_bed_join = peaks_bed.iloc[:,0] + ":" + peaks_bed.iloc[:,1].astype('string') + "-" + peaks_bed.iloc[:,2].astype('string')

peaks_bed_join_list = peaks_bed_join.to_list()

peak_mat_cstb_manual = snap.pp.make_peak_matrix(adata=data_cstb, 
                                                use_rep= peaks_bed_join_list, 
                                                #peak_file='peaks_cstb_q0.001/peaks.combined.bed'
                                               ) 


AnnData object with n_obs × n_vars = 22785 × 274189
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'

# AnnData object with n_obs × n_vars = 22815 × 274189
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden'


(peaks_bed_join_list == peak_mat_cstb_manual.var_names).all() #True


peak_mat = peak_mat_cstb_manual
del peak_mat_cstb_manual


with open('snapshot_pkl/peak_mat.pkl','wb') as fh:
    dill.dump(peak_mat,file=fh)

    
summit_mat = snap.read('snapshot_h5ad/gene_matrix_cstb.h5ad').to_memory()
    

####done   
    

with open('snapshot_pkl/peak_mat.pkl','rb') as fh:
    peak_mat = dill.load(file=fh)

peak_mat

AnnData object with n_obs × n_vars = 22785 × 274189
    obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'
    uns: 't-test', 'wilcoxon', 'logreg'
    obsm: 'X_umap'
    
# AnnData object with n_obs × n_vars = 22785 × 274189
#     obs: 'sample', 'tsse', 'n_fragment', 'frac_dup', 'frac_mito', 'doublet_score', 'is_doublet', 'leiden', 'library', 'leiden_bk', 'cluster', 'cell_type_leiden', 'sex'


peak_mat.X[1:1000,1:1000].todense()    

#peak_mat.raw = peak_mat   #make a backup 
    
peak_mat.raw.X[1:1000,1:1000].todense() 




###add some data_cstb columns

(peak_mat.obs_names == data_cstb.obs_names ).all()#True

peak_mat.obsm['X_umap'] = data_cstb.obsm['X_umap_rotate']

peak_mat.obs['cluster'].dtype

peak_mat.obs['cluster'] = peak_mat.obs['cluster'].cat.reorder_categories(['1','7','9','2','8','3','6','4','5'])


peak_mat.obs['cluster'].value_counts().sort_index()
1    3937
7    1636
9     814 #rm 30
2    3594
8    1447
3    2921
6    3713
4    2493
5    2230

1    3937
2    3594
3    2921
4    2493
5    2230
6    3713
7    1636
8    1447
9     844


          
# ###add peak logic table to uns####

# peak_logic_df = pd.read_csv('MACS_cstb/peaks.combined.logic.txt',sep='\t')

# peak_logic_df['Peaks'] = peak_logic_df.index
# peak_logic_df = peak_logic_df.loc[:,['Peaks','1','2','3','4','5','6','7','8','9']]

# ##substitute uns['peaks']

# all(data_cstb.uns['peaks'] == peak_logic_df) #True
        

    

###use scanpy rank deg test for dar identification####


##use raw peak_mat
peak_mat.X[1:1000,1:1000].todense()
peak_mat.raw.X[1:1000,1:1000].todense()


mmwrite('pmat_DAR/pmat.raw.mtx',peak_mat.raw.X,field = 'integer') 


##save gmat colname and rowname
with open('pmat_DAR//pmat.rowname.txt','w') as fh:
    for i in peak_mat.obs_names.to_list():
        fh.write(i + "\n")

with open('pmat_DAR/pmat.colname.txt','w') as fh:
    for i in peak_mat.var_names.to_list():
        fh.write(i + "\n")
    

##make binary
pmat = peak_mat.raw.X #22785x274189, Compressed Sparse Row format, numpy.uint32

pmat[1:1000,1:1000].todense()#integer indeed

pmat.min() #0
pmat.max() #112

pmat.dtype #uint32

pmat.count_nonzero() #260129883

pmat_bin = np.where( pmat.todense() > 0, 1, 0)
#pmat_bin = pmat_bin.flatten() #mem big and slow!!

pmat_bin.dtype #int64
#pmat_bin.count_nonzero() #int64

pmat_bin[1:1000,1:1000] #correct

pmat_bin.shape #22785, 274189 #no need to reshape
#pmat_bin_reshape = pmat_bin.reshape(-1)


mmwrite('pmat_DAR/pmat.bin.mtx',coo_matrix(pmat_bin),field = 'integer') 


######restart from here#########


#make logarithm
#sc.pp.normalize_per_cell(peak_mat, counts_per_cell_after=1e1) #will be error
sc.pp.log1p(peak_mat) #Natural logarithm by default,X = log(X + 1)

peak_mat.uns['log1p']['base'] = 2.71828


mmwrite('pmat_DAR/pmat.log1p.mtx',peak_mat.X,precision=3) 

peak_mat.X[1:1000,1:1000].todense()
peak_mat.raw.X[1:1000,1:1000].todense()






####do dar test with rank_genes_group (raw peak count by default  )
sc.tl.rank_genes_groups(peak_mat, 'cluster', method='t-test', key_added = "t-test")
sc.pl.rank_genes_groups(peak_mat, n_genes=25, sharey=False, key = "t-test")
##WARNING: It seems you use rank_genes_groups on the raw count data. Please logarithmize your data before calling rank_genes_groups.

sc.tl.rank_genes_groups(peak_mat, 'cluster', method='wilcoxon', key_added = "wilcoxon") #slow but ok
sc.pl.rank_genes_groups(peak_mat, n_genes=25, sharey=False, key = "wilcoxon")

#WARNING: It seems you use rank_genes_groups on the raw count data. Please logarithmize your data before calling rank_genes_groups.

sc.tl.rank_genes_groups(peak_mat, 'cluster', method='logreg', key_added = "logreg") #slow, will use multiple cpus, but ok
sc.pl.rank_genes_groups(peak_mat, n_genes=25, sharey=False, key = "logreg")

#https://github.com/scverse/scanpy/issues/95
#https://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters

#clf = sklearn.linear_model.LogisticRegressionCV()
#clf.fit(peak_mat.X, peak_mat.obs['cluster'])


####do dar test with rank_genes_group (with logarithmized peak count)

sc.tl.rank_genes_groups(peak_mat, 'cluster', method='t-test', key_added = "t-test-log", use_raw = False)
sc.pl.rank_genes_groups(peak_mat, n_genes=25, sharey=False, key = "t-test-log")


sc.tl.rank_genes_groups(peak_mat, 'cluster', method='wilcoxon', key_added = "wilcoxon-log", use_raw = False) #slow but ok
sc.pl.rank_genes_groups(peak_mat, n_genes=25, sharey=False, key = "wilcoxon-log")


sc.tl.rank_genes_groups(peak_mat, 'cluster', method='logreg', key_added = "logreg-log", use_raw = False) #slow, will use multiple cpus, but ok
sc.pl.rank_genes_groups(peak_mat, n_genes=25, sharey=False, key = "logreg-log")



#save with the full dar test table?
with open('snapshot_pkl/peak_mat.pkl','wb') as fh:
    dill.dump(peak_mat,file=fh)

with open('snapshot_pkl/peak_mat.log1p.pkl','wb') as fh:
    dill.dump(peak_mat,file=fh)

    

##extract to df (rawcount-strict)
marker_peaks_df_ttest_rawstrict = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 't-test',pval_cutoff = 0.001, log2fc_min = 1)
marker_peaks_df_wilcoxon_rawstrict = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 'wilcoxon',pval_cutoff = 0.001, log2fc_min = 1)

#marker_peaks_df_logreg = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 'logreg') #error!


###get logreg dar df manually
#peak_mat.uns['t-test'].keys()
#'params', 'names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges'

peak_mat.uns['logreg'].keys()
#'params', 'names', 'scores'

group = list(peak_mat.uns['logreg']['names'].dtype.names)
colnames = ['names', 'scores']
d = [pd.DataFrame(peak_mat.uns['logreg'][c])[group] for c in colnames]
d = pd.concat(d, axis=1, names=[None, 'group'], keys=colnames)
d = d.stack(level=1).reset_index()
d['group'] = pd.Categorical(d['group'], categories=group)
d = d.sort_values(['group', 'level_0']).drop(columns='level_0')


#marker_peaks_df_logreg = d
marker_peaks_df_logreg_rawstrict = d.groupby('group').head(100)



# ##code from scanpy/get/get.py
# d = [pd.DataFrame(adata.uns[key][c])[group] for c in colnames]
# d = pd.concat(d, axis=1, names=[None, 'group'], keys=colnames)
# d = d.stack(level=1).reset_index()
# d['group'] = pd.Categorical(d['group'], categories=group)
# d = d.sort_values(['group', 'level_0']).drop(columns='level_0')

# if pval_cutoff is not None:
#     d = d[d["pvals_adj"] < pval_cutoff]
# if log2fc_min is not None:
#     d = d[d["logfoldchanges"] > log2fc_min]
# if log2fc_max is not None:
#     d = d[d["logfoldchanges"] < log2fc_max]
# if gene_symbols is not None:
#     d = d.join(adata.var[gene_symbols], on="names")

# for pts, name in {'pts': 'pct_nz_group', 'pts_rest': 'pct_nz_reference'}.items():
#     if pts in adata.uns[key]:
#         pts_df = (
#             adata.uns[key][pts][group]
#             .rename_axis(index='names')
#             .reset_index()
#             .melt(id_vars='names', var_name='group', value_name=name)
#         )
#         d = d.merge(pts_df)

# # remove group column for backward compat if len(group) == 1
# if len(group) == 1:
#     d.drop(columns='group', inplace=True)


##save to txt
marker_peaks_df_ttest_rawstrict.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.ttest_rawstrict.txt',sep='\t',index=False)
marker_peaks_df_wilcoxon_rawstrict.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.wilcoxon_rawstrict.txt',sep='\t',index=False)
marker_peaks_df_logreg_rawstrict.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.logreg_rawstrict.txt',sep='\t',index=False)


##count dar by cluster
marker_peaks_df_ttest_rawstrict.group.value_counts().sort_index()
1    62393
7    25666
9     5847
2     6337
8      155
3     4301
6       48
4     2379
5    40118


marker_peaks_df_wilcoxon_rawstrict.group.value_counts().sort_index()
1    15768
7     6420
9      454
2      302
8        3
3     2427
6        6
4       85
5     5324

marker_peaks_df_logreg_rawstrict.group.value_counts().sort_index()
1    100
7    100
9    100
2    100
8    100
3    100
6    100
4    100
5    100


#########rawcount-loose table
marker_peaks_df_ttest_rawloose = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 't-test',pval_cutoff = 0.01, log2fc_min = .2)
marker_peaks_df_wilcoxon_rawloose = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 'wilcoxon',pval_cutoff = 0.01, log2fc_min = .2)
marker_peaks_df_logreg_rawloose = d.groupby('group').head(1000)


##save to txt
marker_peaks_df_ttest_rawloose.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.ttest_rawloose.txt',sep='\t',index=False)
marker_peaks_df_wilcoxon_rawloose.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.wilcoxon_rawloose.txt',sep='\t',index=False)
marker_peaks_df_logreg_rawloose.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.logreg_rawloose.txt',sep='\t',index=False)


marker_peaks_df_ttest_rawloose.group.value_counts().sort_index()
1    94077
7    50117
9    13115
2    37479
8     9962
3    49402
6    20307
4    23082
5    89466

marker_peaks_df_wilcoxon_rawloose.group.value_counts().sort_index()
1    27565
7    13259
9     1003
2     3587
8      766
3    20079
6     4652
4     1753
5    15485

marker_peaks_df_logreg_rawloose.group.value_counts().sort_index()
1    1000
7    1000
9    1000
2    1000
8    1000
3    1000
6    1000
4    1000
5    1000




##extract to df (logcount-strict)
marker_peaks_df_ttest_logstrict = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 't-test-log',pval_cutoff = 0.001, log2fc_min = 1)
marker_peaks_df_wilcoxon_logstrict = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 'wilcoxon-log',pval_cutoff = 0.001, log2fc_min = 1)

#marker_peaks_df_logreg = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 'logreg') #error!


###get logreg dar df manually
#peak_mat.uns['t-test'].keys()
#'params', 'names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges'

peak_mat.uns['logreg-log'].keys()
#'params', 'names', 'scores'

group = list(peak_mat.uns['logreg-log']['names'].dtype.names)
colnames = ['names', 'scores']
d = [pd.DataFrame(peak_mat.uns['logreg-log'][c])[group] for c in colnames]
d = pd.concat(d, axis=1, names=[None, 'group'], keys=colnames)
d = d.stack(level=1).reset_index()
d['group'] = pd.Categorical(d['group'], categories=group)
d = d.sort_values(['group', 'level_0']).drop(columns='level_0')


#marker_peaks_df_logreg = d
marker_peaks_df_logreg_logstrict = d.groupby('group').head(100)



##save to txt
marker_peaks_df_ttest_logstrict.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.ttest_logstrict.txt',sep='\t',index=False)
marker_peaks_df_wilcoxon_logstrict.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.wilcoxon_logstrict.txt',sep='\t',index=False)
marker_peaks_df_logreg_logstrict.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.logreg_logstrict.txt',sep='\t',index=False)



marker_peaks_df_ttest_logstrict.group.value_counts().sort_index()
1    62064
7    25370
9     5773
2     6050
8      113
3     2540
6       34
4     2322
5    39883


marker_peaks_df_wilcoxon_logstrict.group.value_counts().sort_index()
1    15057
7     5744
9      384
2      229
8        0
3      997
6        1
4       55
5     4663


marker_peaks_df_logreg_logstrict.group.value_counts().sort_index()

1    100
7    100
9    100
2    100
8    100
3    100
6    100
4    100
5    100



#########logcount-loose table
marker_peaks_df_ttest_logloose = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 't-test-log',pval_cutoff = 0.01, log2fc_min = .2)
marker_peaks_df_wilcoxon_logloose = sc.get.rank_genes_groups_df(peak_mat,group = None, key = 'wilcoxon-log',pval_cutoff = 0.01, log2fc_min = .2)
marker_peaks_df_logreg_logloose = d.groupby('group').head(1000)


##save to txt
marker_peaks_df_ttest_logloose.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.ttest_logloose.txt',sep='\t',index=False)
marker_peaks_df_wilcoxon_logloose.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.wilcoxon_logloose.txt',sep='\t',index=False)
marker_peaks_df_logreg_logloose.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.logreg_logloose.txt',sep='\t',index=False)


marker_peaks_df_ttest_logloose.group.value_counts().sort_index()
1    94245
7    51584
9    13821
2    38672
8    10558
3    49620
6    20489
4    25252
5    92325




marker_peaks_df_wilcoxon_logloose.group.value_counts().sort_index()
1    27316
7    13252
9     1003
2     3576
8      766
3    20027
6     4633
4     1753
5    15483

marker_peaks_df_logreg_logloose.group.value_counts().sort_index()
1    1000
7    1000
9    1000
2    1000
8    1000
3    1000
6    1000
4    1000
5    1000



##########visualize dars by scanpy plot functions########

marker_peaks_list = {'ttest':{},'wilcoxon':{},'logreg':{}}

marker_peaks_list['ttest']['ttest_rawstrict'] = marker_peaks_df_ttest_rawstrict
marker_peaks_list['ttest']['ttest_rawloose'] = marker_peaks_df_ttest_rawloose
marker_peaks_list['ttest']['ttest_logstrict'] = marker_peaks_df_ttest_logstrict
marker_peaks_list['ttest']['ttest_logloose'] = marker_peaks_df_ttest_logloose


marker_peaks_list['wilcoxon']['wilcoxon_rawstrict'] = marker_peaks_df_wilcoxon_rawstrict
marker_peaks_list['wilcoxon']['wilcoxon_rawloose'] = marker_peaks_df_wilcoxon_rawloose
marker_peaks_list['wilcoxon']['wilcoxon_logstrict'] = marker_peaks_df_wilcoxon_logstrict
marker_peaks_list['wilcoxon']['wilcoxon_logloose'] = marker_peaks_df_wilcoxon_logloose



marker_peaks_list['logreg']['logreg_rawstrict'] = marker_peaks_df_logreg_rawstrict
marker_peaks_list['logreg']['logreg_rawloose'] = marker_peaks_df_logreg_rawloose
marker_peaks_list['logreg']['logreg_logstrict'] = marker_peaks_df_logreg_logstrict
marker_peaks_list['logreg']['logreg_logloose'] = marker_peaks_df_logreg_logloose



##save dar full dar list
with open('DARs_doDAR/marker_peaks_list.pkl','wb') as fh:
    dill.dump(marker_peaks_list,file=fh)


    
    
##compare dar similarity


#for i in marker_peaks_list.keys():
for i in ['rawstrict','rawloose','logstrict','logloose']:
    print('compare condition ' + i)
    
    dar_df_ttest = marker_peaks_list['ttest']['ttest_'+i]
    dar_df_wilcoxon = marker_peaks_list['wilcoxon']['wilcoxon_'+i]
    
    dar_ttest = set(dar_df_ttest['names'].to_list())
    dar_wilcoxon = set(dar_df_wilcoxon['names'].to_list())
    
    dar_share = dar_ttest & dar_wilcoxon
    
    len_ttest = len(dar_ttest)
    len_wilcoxon = len(dar_wilcoxon)
    len_share = len(dar_share)
    
    print('share len: ' + str(len_share) + ' of ttest ' + str(100*len_share/len_ttest) + "%" + "; of wilcoxon " + str(100*len_share/len_wilcoxon) + " %" )
    

####all wilcox result included by ttest


    
    
    
    
##choose one dar set
marker_peaks_df = marker_peaks_list['ttest']['ttest_rawstrict']

#marker_peaks_df = marker_peaks_list['wilcoxon']['wilcoxon_logstrict']


##save dar full table
marker_peaks_df.to_csv('DARs_doDAR/marker_peaks_df_use.txt',sep='\t',index=False)



##compute RPKM from raw count
peak_mat.obs['cluster'].value_counts().sort_index()
1    3937
7    1636
9     814
2    3594
8    1447
3    2921
6    3713
4    2493
5    2230

# 1    3937
# 2    3594
# 3    2921
# 4    2493
# 5    2230
# 6    3713
# 7    1636
# 8    1447
# 9     844





peak_mat_aggre = aggregate_X(peak_mat, groupby='cluster', normalize="RPKM") #Aggregate values in adata.X in a row-wise fashion
9 × 274189
#10 × 483312

peak_mat_aggre.obs_names
#'1', '2', '3', '4', '5', '6', '7', '8', '9'

peak_mat_aggre.obs_names = 'c' + peak_mat_aggre.obs_names
#'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9'


peak_mat_aggre_df = pd.DataFrame(peak_mat_aggre.X,index=peak_mat_aggre.obs_names,columns=peak_mat_aggre.var_names)#.T

peak_mat_aggre_df = peak_mat_aggre_df.loc[['c1','c7','c9','c2','c8','c3','c6','c4','c5'],:]


#count = pl.DataFrame(peak_mat_aggre_df)
#aggre_group x peak  , 10 x 483312

#count = count[:,['c1','c7','c9','c2','c8','c3','c6','c4','c5']]


count = peak_mat_aggre_df
#9 x 274189


# count_z = zscore(
#     np.log2(1 + count.to_numpy()),
#     axis = 1, #by column (vertically)
# )
# ##9 x 274189


#dat = count.iloc[:,0] 
#(dat-dat.mean())/dat.std(ddof=0)


#count_z = pd.DataFrame(count_z,index=count.index,columns=count.columns)


with open('DARs_doDAR/count.pkl','wb') as fh:
    dill.dump(count,file=fh)

with open('DARs_doDAR/count.log1p.pkl','wb') as fh:
    dill.dump(count,file=fh)



##save as dict
marker_peaks_df.group
#Categories (9, object): ['1', '2', '3', '4', ..., '6', '7', '8', '9']

marker_peaks_dict = marker_peaks_df.groupby('group').apply(lambda x: x.names.to_list()).to_dict() #apply to dataframe, select series to dict
type(marker_peaks_dict)
marker_peaks_dict.keys()
['1', '7', '9', '2', '8', '3', '6', '4', '5']
#['1', '2', '3', '4', '5', '6', '7', '8', '9']

marker_peaks_dict['3'][0:10]
marker_peaks_dict['3'][-10:]

for i in marker_peaks_dict.keys():
    print(i + ': ' + str(len(marker_peaks_dict[i])))
1: 62393
7: 25666
9: 5847
2: 6337
8: 155
3: 4301
6: 48
4: 2379
5: 40118

# 1: 62380
# 2: 6341
# 3: 4342
# 4: 2401
# 5: 40169
# 6: 50
# 7: 25678
# 8: 158
# 9: 5541


marker_peaks_c8c3 = set([j for i in ['8','3'] for j in marker_peaks_dict[i]]) 
marker_peaks_c4c5 = set([j for i in ['4','5'] for j in marker_peaks_dict[i]]) 

len(marker_peaks_c8c3) #4450
len(marker_peaks_c4c5) #40830
len(marker_peaks_c8c3 & marker_peaks_c4c5) #1511 shared!


##select top n marker peak to plot
n = 1000
n = 3000
n = 5000

n=10000 #for pairwise logregression method of snapatac2

#marker_peaks_sel_df = marker_peaks_df.groupby('group').head(n).reset_index(drop = True) #keep pvalue sorted but not log2foldchange
marker_peaks_sel_df = marker_peaks_df.groupby('group').apply(lambda x: x.sort_values(['logfoldchanges'], ascending = False).head(n)).reset_index(drop = True) #sort by log2foldchange top1000 use this?

marker_peaks_sel_df = marker_peaks_df.groupby('group').apply(lambda x: x.sort_values(['log2fc'], ascending = False).head(n)).reset_index(drop = True)  #for pairwise logregresson method of snapatac2
marker_peaks_sel_df = marker_peaks_df.groupby('group').head(n).reset_index(drop = True)

marker_peaks_sel_df.group.value_counts().sort_index()

c3    10000
c5    10000

1    5000
7    5000
9    5000
2    5000
8     155
3    4301
6      48
4    2379
5    5000

1    3000
7    3000
9    3000
2    3000
8     155
3    3000
6      48
4    2379
5    3000

1    1000
7    1000
9    1000
2    1000
8     155
3    1000
6      48
4    1000
5    1000

# 1    1000
# 2    1000
# 3    1000
# 4    1000
# 5    1000
# 6      50
# 7    1000
# 8     158
# 9    1000

#marker_peaks_sel_df.to_csv('DARs_doDAR/marker_peaks.scanpy_rank.ttest.top1000.sortby_logfc.txt',sep='\t',index=False)

#marker_peaks_sel_df.to_csv('DARs_doDAR/marker_peaks_sel.scanpy_rank.ttest_rawstrict.top3000.sortby_logfc.txt',sep='\t',index=False)

marker_peaks_sel_df.to_csv('DARs_doDAR/marker_peaks_sel.scanpy_rank.ttest_rawstrict.top5000.sortby_logfc.txt',sep='\t',index=False)

marker_peaks_sel_df.to_csv('DARs_pair_regression/marker_peaks_sel.snapatac2_logiregression.top10000.sortby_logfc.txt',sep='\t',index=False)

#marker_peaks_sel_df = pd.read_csv('DARs_ttest/marker_peaks.scanpy_rank.ttest.top1000.sortby_logfc.txt',sep='\t')


marker_peaks_sel_dict = marker_peaks_sel_df.groupby('group').apply(lambda x: x.names.to_list()).to_dict()

for i in marker_peaks_sel_dict.keys():
    print(i + ': ' + str(len(marker_peaks_sel_dict[i])))

c3: 10000
c5: 10000

1: 5000
7: 5000
9: 5000
2: 5000
8: 155
3: 4301
6: 48
4: 2379
5: 5000


1: 3000
7: 3000
9: 3000
2: 3000
8: 155
3: 3000
6: 48
4: 2379
5: 3000

1: 1000
7: 1000
9: 1000
2: 1000
8: 155
3: 1000
6: 48
4: 1000
5: 1000

# 1: 1000
# 2: 1000
# 3: 1000
# 4: 1000
# 5: 1000
# 6: 50
# 7: 1000
# 8: 158
# 9: 1000


marker_peaks_sel = marker_peaks_sel_df.names 


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



########plot dar peaks######

# ##snap.pl.regions(peak_mat, groupby='cluster', peaks=marker_peaks, interactive=False) #cause problem
# #count = pl.DataFrame(aggregate_X(data, groupby=groupby, normalize="RPKM"))

# peakid = peak_mat.var_names.to_list()

# # idx = peakid.index('chr1:8396853-8397361')# marker_peaks_merge)
# # peakid[idx]

# # idx = peakid.index(['chr1:8396853-8397361','chr1:8392538-8392940']) #wrong for list index matching, only one element for one time!
# # peakid[idx]


# # idx = []
# # for i,v in enumerate(marker_peaks_sel):#slow but ok
# # #    if v in marker_peaks_c5:
# #     if v in peakid:
# #         idx.append(i)

# # idx = []
# # for i in marker_peaks_sel:
# #     for j,v in enumerate(peakid):
# #         if v == i:
# #             idx.append(j)


# idx = []#quick! use this!
# for i in marker_peaks_sel:
#     if i in peakid:
#         idx.append(peakid.index(i))#list index will only return the first one match idx
            

# peakid_sel =   [peakid[i] for i in idx]

# all(count.T.index == peakid) #True


# mat = np.log2(1 + count.T.to_numpy()[idx, :])


# mat.shape
# 7203 x 9

# 13568 x 9
# 239 x 9 #c2
# 3180 x  9 #c5 marker peak
# 13568 x 9 #all marker peak


# ##calculate zscore for dataframe by column, borrow from scenic code

# # mat_Z = pd.DataFrame( index=peak_mat_aggre_df.index ) #create an empty DataFrame

# # for col in list(peak_mat_aggre_df.columns):#very slow
# #     peak_mat_aggre_df_Z[ col ] = ( peak_mat_aggre_df[col] - peak_mat_aggre_df[col].mean()) / peak_mat_aggre_df[col].std(ddof=0)


# # z = scipy.stats.zscore(
# #     np.log2(1 + count.to_numpy()),
# #     axis = 1, #by column
# # )

# mat_Z = zscore(mat,axis = 1,ddof = 0 ) #subset then zscore by column (within all peak sel of one cluster)
# #7203 x 9

# # for col in list(peak_mat_aggre_df.columns):#very slow
# #     peak_mat_aggre_df_Z[ col ] = zscore(peak_mat_aggre_df[ col ])


# #marker_peaks_sel = marker_peaks_df_logreg['names']

# #marker_peaks_sel = peakid_sel


#########directly get mat_Z
mat = count.loc[:,marker_peaks_sel] ##9 x 7203
#mat = count.loc[:,marker_peaks_sel.tolist()] #9 x 7203
mat_Z = zscore(
    np.log2(1 + mat.to_numpy()),
    #axis = 1, #by column (vertically)
    axis = 0, #by row (horizontal, within cluster, use this)
)
##9 x 7203

mat_Z = pd.DataFrame(mat_Z,index=mat.index,columns=mat.columns).T
#7203 x 9

# plt.rcParams['figure.figsize'] =  (4,6)
# plt.rcParams['legend.loc'] = 'upper left'
# plt.rcParams['figure.dpi'] = 100


trace = go.Heatmap(
    #x=count.T.columns,
    #y=[peakid[i] for i in idx], #np.concatenate(list(peaks.values()))[::-1],
    x = mat_Z.columns,
    y = mat_Z.index,
    z=mat_Z,
    zmin = -1.2,
    zmax = 1.2,
    zmid=0,
    type='heatmap',
    colorscale='Blues',#'Blues',#'Reds',#'RdBu_r',#'YlGnBu',#'Viridis',
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
    autosize = False,
    width = 500,
    height = 450,
    # margin=dict(
    #     l=10,
    #     r=0,
    #     b=0,
    #     t=10,
    #     pad=0
    # ),
    # paper_bgcolor="white",
)

fig.show()



###########plot automaticlly###############

##choose one dar set

marker_peaks_df = marker_peaks_list['ttest']['ttest_rawstrict'] #= marker_peaks_df_ttest_rawstrict
#marker_peaks_df = marker_peaks_list['ttest']['ttest_rawloose'] #= marker_peaks_df_ttest_rawloose
#marker_peaks_df = marker_peaks_list['ttest']['ttest_logstrict'] #= marker_peaks_df_ttest_logstrict
#marker_peaks_df = marker_peaks_list['ttest']['ttest_logloose'] #= marker_peaks_df_ttest_logloose


title = 'ttest_rawstrict top5000 sortby logFC'
#title = 'ttest_rawloose top10000 sortby logFC'
#title = 'ttest_logstrict top5000 sortby logFC'
#title = 'ttest_logloose top10000 sortby logFC'


######pairwise c5 vs c3, snapatac2 regression method
#marker_peaks_df = pd.read_csv("DARs_pair_regression/diff_peaks_c5c3.txt",sep='\t')
marker_peaks_df_c3 = pd.read_csv("DARs_pair_regression/diff_peaks_c3.sortlogFC.txt",sep='\t',header= None) #22472
marker_peaks_df_c3.columns = ['names','log2fc','pvalue','adjpvalue']
marker_peaks_df_c3['group'] = 'c3'


marker_peaks_df_c5 = pd.read_csv("DARs_pair_regression/diff_peaks_c5.sortlogFC.txt",sep='\t', header = None)
#22906 
marker_peaks_df_c5.columns = ['names','log2fc','pvalue','adjpvalue']
marker_peaks_df_c5['group'] = 'c5'

marker_peaks_df = pd.concat([marker_peaks_df_c3,marker_peaks_df_c5],axis=0)
#45379

title = 'regresion_pairwise_c3vsc5'

######



####COSG detected dar??#####

marker_peaks_df_cosg = pd.read_csv("/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/DARs_COSG/ranks_stats.mu.100.out3000.txt",sep='\t',header= 0) #top3000
marker_peaks_df_cosg.columns = 'c' + marker_peaks_df_cosg.columns

marker_peaks_df_cosg['idx'] = [str(i) for i in range(0,3000,1)]

marker_peaks_df = pd.wide_to_long(marker_peaks_df_cosg,stubnames='c',i = ['idx'], j = 'group')


marker_peaks_df.columns = ['names']
[  str(i)+' '+str(j) for i,j in marker_peaks_df.index] #multiple index

marker_peaks_df['group'] = [  str(j) for i,j in marker_peaks_df.index] 

marker_peaks_df.index = [  str(i) for i,j in marker_peaks_df.index] 


title = 'cosine_similarity_dar'





########
marker_peaks_df.group.value_counts().sort_index()

c3    22473 #snapatac2 pairwise logistic regression method
c5    22906

1    62393 #ttest_rawstrict
7    25666
9     5847
2     6337
8      155
3     4301
6       48
4     2379
5    40118

7    50117 #ttest_rawloose
9    13115
2    37479
8     9962
3    49402
6    20307
4    23082
5    89466


1    62064 #ttest_logstrict
7    25370
9     5773
2     6050
8      113
3     2540
6       34
4     2322
5    39883

1    94245 #ttest_logloose
7    51584
9    13821
2    38672
8    10558
3    49620
6    20489
4    25252
5    92325


marker_peaks_dict = marker_peaks_df.groupby('group').apply(lambda x: x.names.to_list()).to_dict() 

for i in marker_peaks_dict.keys():
    print(i + ': ' + str(len(marker_peaks_dict[i])))



##select topn
#n = 1000
#n = 3000

#n = 5000
n = 10000

#marker_peaks_sel_df = marker_peaks_df.groupby('group').head(n).reset_index(drop = True) #keep pvalue sorted but not log2foldchange
marker_peaks_sel_df = marker_peaks_df.groupby('group').apply(lambda x: x.sort_values(['logfoldchanges'], ascending = False).head(n)).reset_index(drop = True) #sort by log2foldchange top1000 use this?

marker_peaks_sel_df.group.value_counts().sort_index()

1    5000 #ttest_rawstrict
7    5000
9    5000
2    5000
8     155
3    4301
6      48
4    2379
5    5000

1    10000 #ttest_rawloose
7    10000
9    10000
2    10000
8     9962
3    10000
6    10000
4    10000
5    10000


1    5000 #ttest_logstrict
7    5000
9    5000
2    5000
8     113
3    2540
6      34
4    2322
5    5000

1    10000 #ttest_logloose
7    10000
9    10000
2    10000
8    10000
3    10000
6    10000
4    10000
5    10000


marker_peaks_sel_dict = marker_peaks_sel_df.groupby('group').apply(lambda x: x.names.to_list()).to_dict()

for i in marker_peaks_sel_dict.keys():
    print(i + ': ' + str(len(marker_peaks_sel_dict[i])))
1: 10000
7: 10000
9: 5773
2: 6050
8: 113
3: 2540
6: 34
4: 2322
5: 10000

    
marker_peaks_sel = marker_peaks_sel_df.names


shared_peakid_c3c5 = pd.read_csv('DARs_doDAR/DARs_bed_ttest_rawstrict/shared_peakid_c3c5.txt',sep='\t',header=None)

marker_peaks_sel = shared_peakid_c3c5.loc[:,0]
title = 'shared_peakid_c3c5'

plotDAR_heatmap(marker_peaks_sel = marker_peaks_sel, 
                count = count, #must have rowname and colname, cluster  x peak
                color_use = 'Blues', 
                width= 500, 
                height = 450,
                title = title ,
                save = "DARs_doDAR/pdfs/"+ title.replace(" ","_") + ".pdf"#None
               )



###
def plotDAR_heatmap(marker_peaks_sel = None, count = None, color_use = None, width= None, height = None, title = None, save = None):
    
    #########directly get mat_Z
    
    if not all([True if i in count.columns else False for i in marker_peaks_sel ]):
        print('error: peakid not in count mat: ' + i)
    
    mat = count.loc[:,marker_peaks_sel] ##9 x 7203
    #mat = count.loc[:,marker_peaks_sel.tolist()] #9 x 7203
    mat_Z = zscore(
        np.log2(1 + mat.to_numpy()),
        #axis = 1, #by column (vertically)
        axis = 0, #by row (horizontal, within cluster, use this)
    )
    ##9 x 7203

    mat_Z = pd.DataFrame(mat_Z,index=mat.index,columns=mat.columns).T
    #7203 x 9

    trace = go.Heatmap(
        #x=count.T.columns,
        #y=[peakid[i] for i in idx], #np.concatenate(list(peaks.values()))[::-1],
        x = mat_Z.columns,
        y = mat_Z.index,
        z=mat_Z,
        zmin = -1.2,
        zmax = 1.2,
        zmid=0,
        type='heatmap',
        colorscale= color_use,#'Blues',#'Blues',#'Reds',#'RdBu_r',#'YlGnBu',#'Viridis',
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
        autosize = False,
        width = width,#500,
        height = height,#450,
        # margin=dict(
        #     l=10,
        #     r=0,
        #     b=0,
        #     t=10,
        #     pad=0
        # ),
        # paper_bgcolor="white",
        title_text = title
    )

    fig.show()

    
    if save != None:
        fig.write_image(save)
    
    return('Done')



# #####visualize dar with scanpy heatmap at single cell level (very slow!)

# sc.pl.rank_genes_groups_heatmap(peak_mat,  #slow and need dendrogram from pca
#                                 key = 'wilcoxon-log',
#                                 groupby = 'cluster',
#                                 n_genes=10, 
#                                 use_raw=False, 
#                                 swap_axes=True, 
#                                 show_gene_labels=False,
#                                 vmin=-3, vmax=3,
#                                 cmap='bwr'
#                                )





######quick peakid umap##


plt.rcParams['figure.figsize'] = [4.5,4.5]
plt.rcParams['figure.dpi'] = 100
sc.pl.umap(peak_mat, color=["cluster"],size=10, legend_loc = 'on data', legend_fontoutline = 2,add_outline = True, )


sc.pl.umap(peak_mat,
           use_raw=False, 
           color=["chrX:66638346-66639510",'chr15:64004500-64005140','chr1:109453896-109454696'],
           size=10, 
           color_map = my_cmap,
           vmin = .25,
           #vmax = 1.75,
           show = False,
           title = ['chrX:66638346-66639510','chr15:64004500-64005140','chr1:109453896-109454696'],
          )


##quick&better look the umap

# thresh = {
    
#     'ERVFRD-1':['p30','p100'],
#     'PAPPA':['p65','p99'],
#     'LAMA3':['p50','p100'],
#     'STAT5A':['p10','p99'],
#     'ESRRG':['p50','p99'],
#     'FLT1':['p60','p100'],
#     'ENG':['p60','p100'],
#     'FOSL1':['p50','p99'],
#     'JUNB':['p50','p99'],
#     'FOS':['p30','p99'],
#     'HLA-G':['p30','p100'],
#     'PLAC8':['p50','p100'],
    
#     'DNMT1':['p10','p100'],
#     'PSG8':['p10','p100'],
#     'CGA':['p10','p100'],
#     'SH3TC2':['p60','p100'],
#     'LEP':['p10','p100'],

#     'MYCN':['p60','p100'],    
#     'MYCNUT':['p60','p100'],

#     'FOSL2':['p50','p99'],
#     'JUND':['p30','p99'],
#     'JUN':['p10','p99'],
#     'VDR':['p50','p99'],
#     'AR':['p50','p99'],
#     'MITF':['p50','p99'],
#     'VDR':['p50','p99'],
    
#     'STAT5B':['p50','p99'],
#     'STAT4':['p50','p99']
    
# }


#marker_genes = ["ERVFRD-1","PAPPA",'LAMA3','STAT5A','ESRRG','FLT1','ENG','FOSL1','JUNB','FOS','HLA-G','PLAC8']
marker_genes = ["chrX:66638346-66639510",'chr15:64004500-64005140','chr1:109453896-109454696']




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

xmin = peak_mat.obsm['X_umap'][:,0].min()
xmax = peak_mat.obsm['X_umap'][:,0].max()
ymin = peak_mat.obsm['X_umap'][:,1].min()
ymax = peak_mat.obsm['X_umap'][:,1].max()

shrink = 2

res_p = sc.pl.umap(peak_mat, 
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


for i in range(len(res_p)):
    res_p[i].set_xlim(shrink*xmin,shrink*xmax)
    res_p[i].set_ylim(shrink*ymin,shrink*ymax)
    
#plt.margins(x = [0.7]*8, y = [0.8]*8) #error
#plt.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.3) #not work

plt.suptitle('Marker gene activity score',y = 1.0, x = 0.45)
#plt.tight_layout()




##score dar peakid and plot in one umap

#https://github.com/scverse/scanpy/issues/532
#...That wouldn't give you the sum of the gene expression values, but the average expression minus an average expression of a random gene set.

#for i in marker_peaks_sel_dict.keys():
##for i in marker_peaks_dict.keys():
#for i in ['1','7','9','2','8','3','6','4','5']:
for i in ['c3','c5']:#pairwise snapatac2
    print('do dar visualization for c%s ' % i)
    
    #sc.tl.score_genes(peak_mat, marker_peaks_sel_dict[i], score_name='select c'+i+'_dar_score') #method1 use score_genes
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




###compare within dar dict overlap

marker_peaks_set = {}
for i in marker_peaks_dict.keys():
    marker_peaks_set[i] = set(marker_peaks_dict[i])


#from matplotlib_venn import venn2
#venn2(subsets = marker_peaks_set, set_labels = marker_peaks_set.keys())
#plt.show()


shareMat = pd.DataFrame(index=marker_peaks_set.keys(), columns=marker_peaks_set.keys())

for i in marker_peaks_set.keys():
    for j in marker_peaks_set.keys():
        shareid = marker_peaks_set[i] & marker_peaks_set[j]
        print('shareid of ' + i + ":len " + str(len(marker_peaks_set[i])) + ' vs ' + j + ":len " + str(len(marker_peaks_set[j])) + ' sharelen: ' + str(len(shareid)))
        shareMat.loc[i,j] = len(shareid)

shareMat.dtypes #all object type
#shareMat = shareMat.astype('int') #sns.heatmap need float type




# ##pairwise c3 vs c5, no peakid in common        
# shareid of c3:len 22473 vs c3:len 22473 sharelen: 22473
# shareid of c3:len 22473 vs c5:len 22906 sharelen: 0
# shareid of c5:len 22906 vs c3:len 22473 sharelen: 0
# shareid of c5:len 22906 vs c5:len 22906 sharelen: 22906

shareMat.to_csv('DARs_doDAR/DARs_shareMat//shareMat.ttest.rawstrict.txt',sep='\t')


##quick plot heatmap

from seaborn import heatmap
import seaborn as sns

shareMat = shareMat.astype('float32') #sns.heatmap need float type

sns.heatmap(shareMat, linewidth = 0.5 ,annot=True, )
plt.title( "shared peakid count within marker_peaks_dict" )
plt.show()

#glue = sns.load_dataset("glue").pivot("Model", "Task", "Score")
#np.array(shareMat)

#data_set = glue
#data_set = np.random.rand( 10 , 10 )
#data_set =np.array(shareMat)

# ax = sns.heatmap( data_set , linewidth = 0.5 ,annot=True )
# plt.title( "2-D Heat Map" )
# plt.show()


###save peak_mat with dar table again


with open('snapshot_pkl/peak_mat.pkl','wb') as fh:
    dill.dump(peak_mat,file=fh)

    






######################use diffxpy package for more method??######################








