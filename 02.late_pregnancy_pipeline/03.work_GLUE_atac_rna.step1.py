

##a quick look for GLUE integration of term RNA and ATAC data (use peak - gene interaction, not gene activity score )

##GLUE is developed for unpaired and paired single-celll multi-omics data

import anndata
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams

from ipywidgets import FloatProgress
#! pip install ipywidgets
#!jupyter nbextension enable --py widgetsnbextension


import pickle
import numpy as np

import pandas as pd

#!pip install -U leidenalg #0.7.0 leidenalg will be get problem for scanpy 1.8.2? update to newest 0.9.1

pkg_list = ['anndata','networkx','scanpy','scglue','numpy','leidenalg']
def print_header(pkg_list):
    for pkg in pkg_list:
            try:
                imp = __import__(pkg)
                print (pkg + '==' + imp.__version__, end = ", ")
            except (ImportError, AttributeError):
                pass


print_header(pkg_list)
anndata==0.8.0, networkx==2.8.4, scanpy==1.8.2, scglue==0.2.3, numpy==1.22.4,leidenalg==0.7.0

import leidenalg
leidenalg.version #0.9.1

#Stage 1: Data preprocessing (rna and atac come from one nucleus or two split nuclei batches)

#Read data
rna = anndata.read_h5ad("data/rna.h5ad") #23981 × 28686 #9331 x 22331
rna
rna.obs.domain
rna.obs.cell_type

rna.to_df() #cell x gene
#rna.obs.cell_type.to_csv('rna.cellid.txt',sep='\t')
rna.var_names

atac = anndata.read_h5ad("data/atac.snapatac.h5ad") #24692 × 191816 #8334 x 173870
atac
atac.obs.domain
atac.obs.cell_type
atac.var_names

atac.to_df() #cell x peak 24692 x 191816 #8334 x 173870

#atac.obs.cell_type.to_csv('atac.cellid.txt',sep='\t')


#Preprocess scRNA-seq data
rna.X
rna.X.data


rna.layers["counts"] = rna.X.copy()

sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")

##plot rna
sc.pp.neighbors(rna, metric="cosine")
sc.tl.umap(rna)
sc.pl.umap(rna, color="cell_type")

#sc.tl.leiden(rna, resolution = 0.9)# error

#Preprocess scATAC-seq data
atac.X
atac.X.data
scglue.data.lsi(atac, n_components=100, n_iter=15)

##plot atac
sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.umap(atac)
sc.pl.umap(atac, color="cell_type")

#sc.tl.leiden(atac, resolution = 0.9)

###Construct prior regulatory graph .  the most commonly used prior information linking ATAC peaks with genes is genomic proximity
#Obtain genomic coordinates

rna.var.head()


scglue.data.get_gene_annotation(
    #rna, gtf='data/genes.rna.rename.filter.gtf', #cellranger rna 3.0 gtf with removal of MT and add 'chr'
    #rna, gtf="data/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz",
    rna, gtf='data/genes.gtf',
    gtf_by="gene_name"
)
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].tail()
rna.var.dtypes

flag = np.isfinite(rna.var.loc[:, [ "chromStart", "chromEnd"]]).chromStart
#all(flag) #False
flag.value_counts()
True     27767 #use cellranger 7.0 2020 ref gtf
False      919

True     25615 #use cellranger 3.0 ref gtf
False     3071

#True    22331

idx = np.where(~flag)[0]
rna.var.iloc[idx,].loc[:,["chrom", "chromStart", "chromEnd"]]

pd.Series(rna.var.iloc[idx,].index)


##filter gene without coordinate
rna_filter = rna[:,flag]
23981 × 28686

flag = np.isfinite(rna_filter.var.loc[:, [ "chromStart", "chromEnd"]]).chromStart
#all(flag) #False
flag.value_counts()
True    27767

rna_filter.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
rna_filter.var.loc[:, ["chrom", "chromStart", "chromEnd"]].tail()

rna = rna_filter

del rna_filter

all(rna.var.index == rna.var.genes) #True


rna
23981 × 27767


rna.var.highly_variable.value_counts()
False    25815
True      1952




atac.var_names[:5]

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1])
atac.var["chromEnd"] = split.map(lambda x: x[2])
atac.var.head()

atac.var.dtypes

atac

atac.var.highly_variable.value_counts()
False    160842
True      30974



with open('rna.pickle','wb') as fh:
    pickle.dump(rna,fh)
    
with open('atac.pickle','wb') as fh:
    pickle.dump(atac,fh)
    

# with open('rna.pickle','rb') as fh:
#     rna = pickle.load(fh)
    
# with open('atac.pickle','rb') as fh:
#     atac = pickle.load(fh)
    
    
#stage1: Graph construction

##vertic: peak and gene, edge peak gene proximity
#ATAC peak is connected to a gene if they overlap in either the gene body or promoter region

graph = scglue.genomics.rna_anchored_prior_graph(rna, atac) #very quick in fact
#
   
graph

##check graph
graph.number_of_nodes(), graph.number_of_edges()
(219583, 510819)
#(196201, 439991)

#(270687, 590267)

# Graph node covers all omic features
all(graph.has_node(gene) for gene in rna.var_names), \
all(graph.has_node(peak) for peak in atac.var_names)
#(True, True)

# Edge attributes contain weights and signs
for _, e in zip(range(5), graph.edges):
    print(f"{e}: {graph.edges[e]}")

# Each node has a self-loop
all(graph.has_edge(gene, gene) for gene in rna.var_names), \
all(graph.has_edge(peak, peak) for peak in atac.var_names)
#(True,True)

# Graph is symmetric
all(graph.has_edge(j, i) for i, j, _ in graph.edges)
#True

atac.var.head()

##Save preprocessed data files

rna.write("rna_preprocessed.h5ad", compression="gzip")
atac.write("atac_preprocessed.h5ad", compression="gzip")
nx.write_graphml(graph, "prior.graphml.gz")






