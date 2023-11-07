

####1, use v10nr_clust motif ranking db for dar_full and dar_regression_pairwise_c5c3 (run_pycistarget)
####2, use v10nr_clust (homer only) motif ranking db for dar_regression_pairwise_c5c3 (run_cistarget, not used)

####3, use homer in-build motif for dar_regression_pairwise_c5c3 (run_homer wrapper <-size given> but with cistrome ) (problem occured)

##parameter: auc, ratio_overlap_region_to_db, nes

#%matplotlib inline
#import scenicplus #0.1.dev447+gd4fd733
from scenicplus.wrappers.run_pycistarget import run_pycistarget


import pandas as pd #1.5.3

import pyranges as pr

#pycistarget 1.0.2.dev9+g0ed8289'

#pycistarget 1.0.2.dev9+g0ed8289
from pycistarget.motif_enrichment_cistarget import * #cisTargetDatabase, run_cistarget, run_homer
from pycistarget.motif_enrichment_dem import *

from pycistarget.motif_enrichment_homer import *


from pycistarget.utils import *
from pycistarget.utils import region_names_to_coordinates

#import pycistarget

import os
import sys
import dill

##Motif enrichment analysis using pycistarget


##read in dar

#dar_full
#region_df = pd.read_csv('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/DARs_doDAR/DARs_bed_ttest_rawstrict/DARs.clusterALL.ttest_rawstrict.filtersharedc3c5.bed',sep='\t',header=None,)



##1 read in cut&tag / ChIP-seq peaks
#region_df = pd.read_csv('CS1_vs_IgG_peaks.narrowPeak',sep='\t',header=None) #50377

region_df = pd.read_csv('peaks.combined.CEBPB.bed',sep='\t',header=None) #86242
region_df = pd.read_csv('peaks.combined.FOSL2.bed',sep='\t',header=None)  #178367


##2 read in eGRN
CEBPB_cistrom_target = pd.read_csv("/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/cisTarget/combine_compare_eGRN_cistarget/tf_footprint/CEBPB.cistrome.target.bed",sep='\t',header=None)
#428 x 6

FOSL2_cistrom_target = pd.read_csv("/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/cisTarget/combine_compare_eGRN_cistarget/tf_footprint/FOSL2.cistrome.target.bed",sep='\t',header=None)
#461 x 6



#pairwise c5c3
#region_df = pd.read_csv('/sda/mjwang/pwdex/placenta_10X_revision/placenta_10X_early_combine/02.snapATAC_harmony/DARs_pair_regression/diff_peaks_pairwiseALL.sortlogFC.txt.addid.bed',sep='\t',header=None,)

#http://genome.ucsc.edu/FAQ/FAQformat.html#format12
#narrowpeak format
##region_df.columns = ['chr','start','end','peakname','score','peakid','signalValue','pValue','qValue','peak']
region_df.columns = ['chr','start','end']
#region_df.columns = ['chr','start','end','peakid_cluster']
#temp = region_df.loc[:,'peakid_cluster'].str.split('_',expand = True)


region_df['peakid'] =  [ "chr"+region_df.iloc[i,0] + ":" + str(region_df.iloc[i,1]) +"-" + str(region_df.iloc[i,2]) for i in range(len(region_df)) ]  #'chr'+region_df['chr']+":"+str(region_df['start'].to_list() )+ str(region_df['end'])

CEBPB_cistrom_target.columns = ['chr','start','end','peak-target','score','strand']

CEBPB_cistrom_target['peakid'] =  [ CEBPB_cistrom_target.iloc[i,0] + ":" + str(CEBPB_cistrom_target.iloc[i,1]) +"-" + str(CEBPB_cistrom_target.iloc[i,2]) for i in range(len(CEBPB_cistrom_target)) ] 


FOSL2_cistrom_target.columns = ['chr','start','end','peak-target','score','strand']

FOSL2_cistrom_target['peakid'] =  [ FOSL2_cistrom_target.iloc[i,0] + ":" + str(FOSL2_cistrom_target.iloc[i,1]) +"-" + str(FOSL2_cistrom_target.iloc[i,2]) for i in range(len(FOSL2_cistrom_target)) ] 





#region_df.loc[:,'peakid'] = temp.iloc[:,0]
#region_df.loc[:,'cluster'] = temp.iloc[:,1]
# region_df.loc[:,'cid'].value_counts()
# cluster1    24066
# cluster7    16831
# cluster5     4444
# cluster9     2568
# cluster3     1428
# cluster4      654
# cluster2      328
# cluster8       41
# cluster6       17

# cluster5    22906
# cluster3    22473




#region_set = region_df.groupby(region_df.cid)

region_set = {}

#for i in region_df.cluster.value_counts().sort_index().index:
# for i in region_df.cid.value_counts().sort_index().index:
#     #df_sel = region_df.loc[region_df.cluster == i,:]
#     df_sel = region_df.loc[region_df.cid == i,:]
#     #region_set[i] = df_sel.drop(columns=['peakid_cluster','cluster'])
#     region_set[i] = df_sel.drop(columns='cid')
#     print(i+ ' = ' + str(df_sel.shape[0]) )

# cluster1 = 24066 #dar overlap p2g
# cluster2 = 328
# cluster3 = 1428
# cluster4 = 654
# cluster5 = 4444
# cluster6 = 17
# cluster7 = 16831
# cluster8 = 41
# cluster9 = 2568

    
# cluster3 = 22473
# cluster5 = 22906


# cluster1 = 62393 #full dar
# cluster2 = 6337
# cluster3 = 2806
# cluster4 = 2379
# cluster5 = 36863
# cluster6 = 48
# cluster7 = 25666
# cluster8 = 155
# cluster9 = 5847


##to see whether or not peakid of interest are in ctx db metaregion

##peak near CSH1, CSH2 CSHL1, GH loci
# region_test = ['chr17:63882671-63883416','chr17:63873729-63873957','chr17:63912030-63912428','chr17:63909698-63910045','chr17:63911267-63911561','chr17:63896549-63896943','chr17:63873729-63873957']

# ##promotor peak and upstream downstream peaks of ERVFRD-1
#                                              #promoter
# region_test = ['chr6:11110923-11112506','chr6:11102469-11103550','chr6:11116034-11116522','chr6:11129830-11130560']


# for i in region_set.keys():
#     print(i)
#     print(str(region_set[i].peakid.isin(region_test).sum()))
# #ERVFRD-1 cluster 9 found 3

# pd.Series(region_test).isin(region_set['cluster9'].peakid)
# 0     True
# 1     True
# 2    False
# 3     True


    
###3 overlap summary for cut&tag / ChIPseq peaks with eGRN peaks

region_df_pr = pr.PyRanges(region_names_to_coordinates(region_df['peakid']))
#178367 FOSL2
#86242 CEBPB merge 3reps
#51429 CEBPB rep1

CEBPB_cistrom_gr= pr.PyRanges(region_names_to_coordinates(CEBPB_cistrom_target['peakid']))
#428

CEBPB_cistrom_gr.intersect(region_df_pr)
302 (302/428, 70.6%)
#248 (57.9%)

#CEBPB_cistrom_gr.overlap(region_df_pr)
##238



FOSL2_cistrom_gr= pr.PyRanges(region_names_to_coordinates(FOSL2_cistrom_target['peakid']))
#461 

FOSL2_cistrom_gr.intersect(region_df_pr)
344 (344/461 74.6%)



region_set['CEBPB_cuttag'] = region_df
region_set['CEBPB_cistrom'] = CEBPB_cistrom_target



region_set_pr = {}

for k in region_set.keys():
    region_set_pr[k] = pr.PyRanges(region_names_to_coordinates(region_set[k].peakid))
51,429
428

# region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
# region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
# markers_dict = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'rb'))


# region_sets = {}
# region_sets['topics_otsu'] = {}
# region_sets['topics_top_3'] = {}
# region_sets['DARs'] = {}

# #wiil do for three region sets for CTX and DEM method with/without promoters !!
# ##total 12 dirs


# for topic in region_bin_topics_otsu.keys():
#     regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
#     region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

# for topic in region_bin_topics_top3k.keys():
#     regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
#     region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

# for DAR in markers_dict.keys():
#     regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
#     region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))


##4 verify cut&tag/ChIPseq peak motif entrichment

for key in region_set_pr.keys():
    print(f'{key}: {region_set_pr[key].keys()}')


# db_fpath = "/staging/leuven/stg_00002/lcb/icistarget/data/make_rankings/v10_clust/CTX_hg38"
# motif_annot_fpath = "/staging/leuven/stg_00002/lcb/cbravo/cluster_motif_collection_V10_no_desso_no_factorbook/snapshots"


# rankings_db = os.path.join(db_fpath, 'cluster_SCREEN.regions_vs_motifs.rankings.v2.feather')
# scores_db =  os.path.join(db_fpath, 'cluster_SCREEN.regions_vs_motifs.scores.v2.feather')
# motif_annotation = os.path.join(motif_annot_fpath, 'motifs-v10-nr.hgnc-m0.00001-o0.0.tbl')


db_fpath = "/home/mjwang/dataex/cisTargetDb_v10/"
motif_annot_fpath = "/home/mjwang/dataex/cisTargetDb_v10/motif2tf/"

rankings_db = os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather' )
motif_annotation = os.path.join(motif_annot_fpath, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")


#do pycisTarget with wrapper run_pycistarget 

#work_dir = 'result_pycisTarget'
#work_dir = 'result_pycisTarget_p2g_auc0.001_metaoverlap0.2'
work_dir = 'result_pycisTarget_CEBPB_CS1-vs-IgG_auc0.005'
#work_dir = 'result_pycisTarget_p2g_auc0.001'
#work_dir = 'result_pycisTarget_pairwise_c5c3'
#work_dir = 'result_pycisTarget_pairwise_c5c3_homer'


# tmp_dir = './'

# if not os.path.exists(os.path.join(work_dir, 'motifs')):
#     os.makedirs(os.path.join(work_dir, 'motifs'))

##read in cisTargetDatabase

##ctx_db = cisTargetDatabase(fname = rankings_db) #Loading complete cistarget database slow! and mem big!
ctx_db = cisTargetDatabase(fname = rankings_db, region_sets = region_set_pr)#,fraction_overlap = 0.2)
#dem_db = DEMDatabase(scores_db, region_set_pr) 

ctx_db.db_rankings.shape
5876, 54488

#5876, 42430, p2g, 0.2 overlap ratio
#5876, 43595, p2g, 0.1 overlap ratio
#5876, 40480, p2g
#5876, 84150, pairwise dar
#5876, 1837304, motif x metaregion, full database
#5876, 111964, motif x region, only overlap with region_set_gr

ctx_db.db_rankings.index
ctx_db.db_rankings.columns


#####check metaregion overlap peak of interest####

##ERVFRD-1 promoter peak (chr6:11110923-11112506) overlap meta-region (chr6:11111459-11111710)
#ctx_db.db_rankings.columns.isin(['chr6:11111460-11111710']).sum() #not in
#ctx_db_rankings_pr = pr.PyRanges(region_names_to_coordinates(ctx_db.db_rankings.columns))
#ctx_db_rankings_pr.overlap( pr.PyRanges(region_names_to_coordinates(['chr6:11111460-11111710'])) )

#chr6:11111459-11111710
#ctx_db.db_rankings.columns.isin(['chr6:11111459-11111710']).sum() #1



#ERVFRD-1  peak  overlap meta-region
#ctx_db_rankings_pr.overlap( pr.PyRanges(region_names_to_coordinates(region_test)) )
#0	chr6	11111459	11111710
#1	chr6	11129769	11130111
#2	chr6	11130259	11130492
    

###check done



# metaregion_df = pd.Series(ctx_db.db_rankings.columns).str.split(":|-",expand=True)
# metaregion_df.columns = ['chr','start','end']

# metaregion_df.to_csv('metaregion.bed',sep='\t',header=False,index=False)





# homer_motif = pd.read_csv('/home/mjwang/dataex/cisTargetDb_v10/motif2tf/motif.homer.addH.tbl',sep='\t',header = 0)
# 3109 x 13

# homer_motif.motif_id.isin( ctx_db.db_rankings.index).value_counts()
# True     3013
# False      96

# len(set(homer_motif.motif_id)) #308


# pd.Series(ctx_db.db_rankings.index).isin(homer_motif.motif_id).value_counts()
# False    5630
# True      246

# flag = pd.Series(ctx_db.db_rankings.index).isin(homer_motif.motif_id).to_list()

# sum(flag) #246 of 5876

# ctx_db_homer = ctx_db
# ctx_db_homer.db_rankings = ctx_db.db_rankings.iloc[flag,:]

# #ranking db:
# 246 x 84150
# #regions_to_db:
# 52750 x 2
# 31402  x 2


# ##to see whether or not peakid of interest are in ctx db metaregion

# ##CSH1 CSH2 CSHL1 GH2 near peak
# region_test = ['chr17:63882671-63883416','chr17:63873729-63873957','chr17:63912030-63912428','chr17:63909698-63910045','chr17:63911267-63911561','chr17:63896549-63896943','chr17:63873729-63873957']

# region_test = ['chr6:11110923-11112506','chr6:11102469-11103550','chr6:11116034-11116522','chr6:11129830-11130560'] #ERVFRD-1 near peak

# for i in ctx_db.regions_to_db.keys():
#     print(i)
#     print(str(ctx_db.regions_to_db[i].Target.isin(region_test).sum()))
#     print(str(ctx_db.regions_to_db[i].Query.isin(region_test).sum() ))

# ##CSH1 CSH2 CSHL1 GH2  near peak

# #cluster5 found 1, cluster 3 found 1 (p2g overlap0.2)
# #cluster3 found 1, cluster5 found 1 (p2g)
# #cluster5 found 1 (pairwise)
# #cluster5 found 3, cluster 3 found 1 (full dar)

# #ERVFRD-1
# #cluster9 found 3 (p2g overlap0.2)

# #pd.concat([ctx_db.regions_to_db[x] for x in ctx_db.regions_to_db.keys()])['Query']

# # region_set_pr_test = { 'dar_ttest_rawstrict.filtersharedc3c5': {  'c1':  region_set_pr['cluster1'], 
# #                                                                   'c5':  region_set_pr['cluster5'] 
                       
# #                                                                 }
#                      }

region_set_pr_dict = {'peak_CEBPB_target_cuttag': region_set_pr}
#region_set_pr_dict = {'dar_regression_pairwise_c5c3': region_set_pr} #two levels
#region_set_pr_dict = {'dar_ttest_rawstrict.filtersharedc3c5': region_set_pr} #two levels


##will run cistarget and DEM analysis (wrapper for run_cistarget and DEM)
run_pycistarget( #will take more than 1h!
    region_sets = region_set_pr_dict,
    species = 'homo_sapiens',
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,#run ctx by default
    #dem_db_path = scores_db, #run dem if not empty
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True, #will both run with/without promoter
    n_cpu = 10, #> or eq to dar cluster number ?
    #_temp_dir = tmp_dir,#os.path.join(tmp_dir, 'ray_spill'), #better set None
    annotation_version = 'v10nr_clust',
    save_partial = True, #will keep temp pkl #trun this switch on, or dill will cause fail, do not know why
    ctx_auc_threshold = 0.001 #0.005
)

##additional useful parameters
#ctx_auc_threshold: float = 0.005
#promoter_space: int = 500
#ctx_nes_threshold: float = 3.0
#annotation : List[str] = ['Direct_annot', 'Orthology_annot']
#exclude_motifs: path of csv file of motifs to exclude from the analysis 
#exclude_collection: list of strings identifying which motif collections to exclude from analysis


# ###directly run cistarget with subset of rank database (only homer motifs with region_sets)
# cistarget_dict = run_cistarget(ctx_db = ctx_db_homer,
#                                region_sets = region_set_pr_dict['dar_regression_pairwise_c5c3'],
#                                specie = 'homo_sapiens',
#                                annotation_version = 'v10nr_clust', #Use v10nr_clust for the latest motif collection
#                                auc_threshold = 0.005,
#                                nes_threshold = 3.0,
#                                rank_threshold = 0.05,
#                                path_to_motif_annotations = motif_annotation,
#                                annotation = ['Direct_annot', 'Orthology_annot'],
#                                n_cpu = 10,
#                                #motifs_to_use = motif_list,
#                                #_temp_dir='/scratch/leuven/313/vsc31305/ray_spill'
#                               )


# ##additional useful parameters
# #ctx_db: Path to the cisTarget database to use, or a preloaded cisTargetDatabase object (using the same region sets to be analyzed



# menr_ctx = cistarget_dict


#####





menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))
#'CTX_dar_ttest_rawstrict.filtersharedc3c5_All', 'CTX_dar_ttest_rawstrict.filtersharedc3c5_No_promoters', 'DEM_dar_ttest_rawstrict.filtersharedc3c5_All', 'DEM_dar_ttest_rawstrict.filtersharedc3c5_No_promoters'

'CTX_peak_CEBPB_target_cuttag_All', 'CTX_peak_CEBPB_target_cuttag_No_promoters'



#'CTX_dar_ttest_rawstrict.filtersharedc3c5.p2g_All', 'CTX_dar_ttest_rawstrict.filtersharedc3c5.p2g_No_promoters'
#'CTX_dar_regression_pairwise_c5c3_All', 'CTX_dar_regression_pairwise_c5c3_No_promoters'
#'CTX_dar_ttest_rawstrict.filtersharedc3c5_All', 'CTX_dar_ttest_rawstrict.filtersharedc3c5_No_promoters'

menr_ctx = menr['CTX_peak_CEBPB_target_cuttag_No_promoters']
#menr_ctx = menr['CTX_dar_regression_pairwise_c5c3_No_promoters']
#menr_ctx = menr['CTX_dar_ttest_rawstrict.filtersharedc3c5_No_promoters']
#menr_ctx = menr['CTX_dar_ttest_rawstrict.filtersharedc3c5_All']

CEBPB_cuttag , CEBPB_cistrom





#menr = dill.load(open(os.path.join(work_dir, #'motifs/CTX_dar_ttest_rawstrict.filtersharedc3c5_No_promoters.pkl'), 'rb'))




##########show as html result table###########3


cistarget_results(menr_ctx,name = 'CEBPB_cuttag')
cistarget_results(menr_ctx,name = 'CEBPB_cistrom')


##verified!!




#####pipeline stop#######















# # menr_ctx_c5 = cistarget_results(menr_ctx,name = 'cluster5')
# # menr_ctx_c3 = cistarget_results(menr_ctx,name = 'cluster3')

# # menr_ctx_c5
# # menr_ctx_c3




# # menr_ctx_c5 = cistarget_results(menr_ctx,name = 'cluster5')
# # menr_ctx_c4 = cistarget_results(menr_ctx,name = 'cluster4')

# # menr_ctx_c5
# # menr_ctx_c4

# # menr_ctx_c3 = cistarget_results(menr_ctx,name = 'cluster3')
# # menr_ctx_c8 = cistarget_results(menr_ctx,name = 'cluster8')

# # menr_ctx_c3
# # menr_ctx_c8



# # menr_ctx_c9 = cistarget_results(menr_ctx,name = 'cluster9')
# # menr_ctx_c2 = cistarget_results(menr_ctx,name = 'cluster2')

# # menr_ctx_c9
# # menr_ctx_c2


# # menr_ctx_c1 = cistarget_results(menr_ctx,name = 'cluster1')
# # menr_ctx_c7 = cistarget_results(menr_ctx,name = 'cluster7')

# # menr_ctx_c1
# # menr_ctx_c7


# #########:::::::::::::::::::access result, output enrichment table ::::::::::::::######

# cid = 'cluster3'
# cid = 'cluster8'

# cid = 'cluster5'
# cid = 'cluster4'

# dir(menr_ctx[cid])
# #  'add_motif_annotation_cistarget',
# #  'annotation',
# #  'annotation_version',
# #  'auc_threshold',
# #  'cistromes',
# #  'motif_enrichment',
# #  'motif_hits',
# #  'motif_similarity_fdr',
# #  'motifs_to_use',
# #  'name',
# #  'nes_threshold',
# #  'orthologous_identity_threshold',
# #  'path_to_motif_annotations',
# #  'rank_threshold',
# #  'region_set',
# #  'regions_to_db',
# #  'run_ctx',
# #  'specie'
    

# ##add motif peak hits (and near target gene candidate?) for each motif row in enrichment table ?

# ##0 pycistarget parameters
# menr_ctx[cid].annotation_version #v10nr_clust
# menr_ctx[cid].auc_threshold #0.001 #0.005 
# menr_ctx[cid].nes_threshold #3.0


# ####1 map region to region_db################
# menr_ctx[cid].region_set #2674 of total 4301 found db neta region
# menr_ctx[cid].regions_to_db #3071

# menr_ctx[cid].region_set.peakid = menr_ctx[cid].region_set.Chromosome.astype('string') + ":" +	menr_ctx[cid].region_set.Start.map(str) + "-" +	menr_ctx[cid].region_set.End.map(str) 


# #Target store the dar region, query is the db region
# menr_ctx[cid].regions_to_db.Target.isin(menr_ctx[cid].region_set.peakid).value_counts()
# True    3071
# menr_ctx[cid].regions_to_db.Query.isin(menr_ctx[cid].region_set.peakid).value_counts()
# False    3071

# #check dar of this cluster
# dar_full = pd.read_csv('../../cicero_Granja/peakAnno_merged/dar_all_pressel_2018.txt',sep='\t',header = 0)
# dar_full = region_df
# all(menr_ctx[cid].region_set.peakid.isin(dar_full.loc[dar_full.cid == cid,:].peakid)) #True
# #all(menr_ctx[cid].region_set.peakid.isin(dar_full.loc[dar_full.cluster == cid,:].peakid))
# #True



# ##2 motif hits peak region
# menr_ctx[cid].motif_hits #to bed, and assign near gene with ChIPseeker and GREAT result

# menr_ctx[cid].motif_hits['Region_set'].keys() == menr_ctx[cid].motif_hits['Database'].keys() #TRUE

# menr_ctx[cid].motif_hits['Region_set']['metacluster_46.4']
# menr_ctx[cid].motif_hits['Database']['metacluster_46.4']



# ##3 motif level enrichment table (similar to i-cistarget but no target gene assignment)
# menr_ctx[cid].motif_enrichment

# menr_ctx[cid].motif_enrichment.columns
# 'Logo', 'Region_set', 'Direct_annot', 'Orthology_annot', 'NES', 'AUC','Rank_at_max', 'Motif_hits'



# ##quick look up weather or not TF of interest exist
# tf_anno_direct = menr_ctx[cid].motif_enrichment.Direct_annot
# tf_anno_direct = tf_anno_direct[tf_anno_direct.notna().to_list()]

# sum(tf_anno_direct.str.contains('FOSL2').dropna()) #6
# sum(tf_anno_direct.str.contains('CEBPB').dropna()) #5
# sum(tf_anno_direct.str.contains('HDAC').dropna()) #0

# sum(tf_anno_direct.str.contains('STAT5A').dropna()) #4 in c5, 4 in c4
# sum(tf_anno_direct.str.contains('MITF').dropna()) #1 in c4, 0 in c5
# sum(tf_anno_direct.str.contains(' AR,').dropna()) #0 in c5, 1 in c8?


# idx = tf_anno_direct.str.contains('STAT5A')
# idx = tf_anno_direct.str.contains('MITF')
# idx = tf_anno_direct.str.contains(' AR,')

# tf_anno_direct[idx.to_list()]




# ##quick look done




# ###read and choose peakAnno to add ChIPSeeker and GREAT precalculated peak annotation to  motif hits

# ##chipseeker all peak annotation
# peakAnno_chipseeker = pd.read_csv('../../cicero_Granja/ChIPSeeker/peaks.gr.anno.filter.txt',sep='\t',header = 0)
# 274147 x 19 #all peaks

# ##great dar annotation
# peakAnno_great_dar = pd.read_csv('../../cicero_Granja/peakAnno_merged/peakAnno_great_dar.all.pressel_2018.txt',sep='\t',header = 0)

# 147244 × 7 # equal to dar full


# ##merge two method chipseeker and great annotation overlap highconf annotation

# peakAnno_merge_chipseekergreat = pd.read_csv('../../cicero_Granja/peakAnno_merged/peakAnno.all.merge.txt',sep='\t',header = 0, dtype={'annotation':'string','gene ':'string','gene1': 'string','gene2':'string','dist1': 'string', 'dist2': 'string'},na_filter = False)#keep_default_na = True)#,na_filter = False,na_values = ['NA'] )
# 274147 x 23 #all peaks, use this

# peakAnno_merge_chipseekergreat.dtypes


# ##check if NA is correctly read
# peakAnno_merge_chipseekergreat.loc['chrY:56881388-56881989',].isna()
# #NA treat as string!! good





# peakAnno_merge_chipseekergreat_highconf = pd.read_csv('../../cicero_Granja/peakAnno_merged/peakAnno.all.merge.highconf.txt',sep='\t',header = 0 )

# 196511 × 23 #peak annotation agree between two method


# peakAnno_merge_chipseekergreat_dar_highconf = pd.read_csv('../../cicero_Granja/peakAnno_merged/peakAnno.merge_chipseekergreat.dar.highconf.txt',sep='\t',header = 0)
# 79546 x 23 #dar annotation with two method overlap


# ##choose one of them
# #peakAnno =  peakAnno_merge_chipseekergreat_dar_highconf
# peakAnno =  peakAnno_merge_chipseekergreat #use this ?
# len(set(peakAnno.index)) == len(peakAnno.index) #True, no dups

# ###read and choose peakAnno done



# ########output motif level result table########


# menr_ctx.keys()
# #['cluster1', 'cluster2', 'cluster3', 'cluster4', 'cluster5', 'cluster6', 'cluster7', 'cluster8', 'cluster9']

# cid = 'cluster5'#'cluster1'#'cluster3'

# list(menr_ctx[cid].motif_hits['Region_set'].keys() ) == menr_ctx[cid].motif_enrichment.index.to_list() #True

# motifid = 'metacluster_29.2'#'metacluster_171.1'#'metacluster_137.2'


# menr_ctx[cid].motif_enrichment.iloc[0:20,0:4]





# #for cid in ['cluster3']:
# for cid in menr_ctx.keys():
#     print('processing cluster ' + cid)
    
#     #make bed dir for this cluster
#     if not os.path.isdir(work_dir + '/motif_beds/' + cid):
#         os.makedirs(work_dir + '/motif_beds/' + cid, )#will mkdir recursively
#         #os.mkdir(work_dir + '/motif_beds/' + cid, )
    
#     outdir = work_dir + '/motif_beds/' + cid + '/'
    
#     tf_motif_map = {} 
    
#     for motifid in menr_ctx[cid].motif_enrichment.index:
#         print('processing motif ' + motifid)
#         #len(menr_ctx[cid].motif_hits['Region_set'][motifid]) #321
#         #len(menr_ctx[cid].motif_hits['Database'][motifid]) #336

#         if not len(menr_ctx[cid].motif_hits['Database'][motifid]) == menr_ctx[cid].motif_enrichment.loc[motifid,'Motif_hits']:
#             print('motif_hits number not recorded rightly in motif_enrichment table')
#             sys.exit() 
#         #motif_enrichment record db_region hits in Motif_hits column

#         region_motif_hits =  menr_ctx[cid].motif_hits['Region_set'][motifid]
#         region_motif_hits_db =  menr_ctx[cid].motif_hits['Database'][motifid]

#         if not len(set(region_motif_hits)) == len(region_motif_hits): #no duplication
#             print('duplication in region_motif_hits')
#             sys.exit()

#         if not len(set(region_motif_hits_db)) == len(region_motif_hits_db): #no duplication
#             print('duplication in region_motif_hits_db')
#             sys.exit()
            
            
#         if not pd.Series(region_motif_hits).isin(dar_full.peakid).all(): #True,is dar region
#             print('region_motif_hits must in dar_full region')
#             sys.exit()

#         ##output bed
#         region_motif_hits_df = region_names_to_coordinates(region_motif_hits)
#         region_motif_hits_df.loc[:,'peakid'] = region_motif_hits_df.loc[:,'Chromosome'] + ":" + region_motif_hits_df.loc[:,'Start'].map(str) + '-' + region_motif_hits_df.loc[:,'End'].map(str)
#         region_motif_hits_df.to_csv(outdir + motifid + '.peaks.bed',sep='\t',header=False,index = False)
        
#         region_motif_hits_db_df = region_names_to_coordinates(region_motif_hits_db)
#         region_motif_hits_db_df.loc[:,'metaid'] = region_motif_hits_db_df.loc[:,'Chromosome'] + ":" + region_motif_hits_db_df.loc[:,'Start'].map(str) + '-' + region_motif_hits_db_df.loc[:,'End'].map(str)
#         region_motif_hits_db_df.to_csv(outdir + motifid + '.metaregion.bed',sep='\t',header=False,index = False)
    
            
#         if not len(set(region_motif_hits).intersection( set(peakAnno.index) ) ) == len(region_motif_hits): #all has annotation in peakAnno table
#             print('not all peak in region_motif_hits has annotation,will filter regions without annotation')
#             print ( str(len(region_motif_hits) - len(set(region_motif_hits).intersection( set(peakAnno.index) ) )) + ' has no annotaton, filtering...' )
            
#             region_motif_hits = [v for(v,i) in zip(region_motif_hits ,pd.Series(region_motif_hits).isin(peakAnno.index).to_list()) if i]
            
#             #sys.exit()

#         region_motif_hits_anno_df = peakAnno.loc[region_motif_hits,:]

#         if not all(region_motif_hits_anno_df.index == region_motif_hits): #True
#             print('not all peak in region_motif_hits has annotation')
#             sys.exit()

#         ##region_motif_hits_anno_df.columns
#         # 'seqnames', 'start', 'end', 'width', 'strand', 'annotation', 'geneChr',
#         #        'geneStart', 'geneEnd', 'geneLength', 'geneStrand', 'geneId',
#         #        'transcriptId', 'distanceToTSS', 'ENSEMBL', 'SYMBOL', 'GENENAME',
#         #        'gene', 'site_name', 'gene1', 'dist1', 'gene2', 'dist2'

#         if not (region_motif_hits_anno_df.gene == region_motif_hits_anno_df.SYMBOL).all(): #True
#             print('in region_motif_hit_anno_df, gene column != SYNBOL')
#             sys.exit()
#         #(region_motif_hits_anno_df.index == region_motif_hits_anno_df.site_name).all()

#         ##gene = [i for i in region_motif_hits_anno_df.gene.value_counts().sort_index().index] #chipseeker assigned gene
#         ##gene1 = [i for i in region_motif_hits_anno_df.gene1.value_counts().sort_index().index] #great assigned gene1
#         ##gene2 = [i for i in region_motif_hits_anno_df.gene2.value_counts().sort_index().index] #great assigned gene2


#         ##add a formated annotation string

#         region_motif_hits_anno_df.loc[:,'region_set'] = region_motif_hits_anno_df.index.to_list()

#         region_motif_hits_anno_df.loc[:,'annotype'] = region_motif_hits_anno_df.annotation.str.split(' ',expand = True).iloc[:,0]

#         #region,chipseeker gene|dist|type,great gene1|dist1|gene2|dist2
#         region_motif_hits_anno_df.loc[:,'anno_string'] = region_motif_hits_anno_df.loc[:,'region_set'] + "," + region_motif_hits_anno_df.loc[:,'gene'] + '|' + region_motif_hits_anno_df.loc[:,'distanceToTSS'].map(str) + '|' + region_motif_hits_anno_df.loc[:,'annotype'] + ',' + region_motif_hits_anno_df.loc[:,'gene1'] + '|' + region_motif_hits_anno_df.loc[:,'dist1'].map(str) + '|' + region_motif_hits_anno_df.loc[:,'gene2'] + '|' + region_motif_hits_anno_df.loc[:,'dist2'].map(str)

#         ##add region_hits, region_anno, region_db hits save as a i-cistarget compatible table######

#         #motif, tf (, joined), NES, topTargetMeta(, joined), peaks(, joined), nearGene(, joined)
#         # enrichTable <- rbind.data.frame(enrichTable, data.frame(  
#         #             'motif' = motifid, 
#         #             'tf' = paste(tfid,collapse = ','),
#         #             'NES' = nes,
#         #             'topTargetMeta' = paste(topTargetMeta,collapse = ','),
#         #             'peaks' = paste(peaks,collapse = ','),
#         #             'nearGene' = paste(nearGene,collapse = ',')


#         #          ) )

        
#         print ('adding region_hits, region_anno, region_db to motif_enrichment table')
#         menr_ctx[cid].motif_enrichment.loc[motifid,'region_hits'] = ','.join(region_motif_hits_anno_df.index.to_list()) #join annodf one column to a string

#         menr_ctx[cid].motif_enrichment.loc[motifid,'region_anno'] = ' '.join(region_motif_hits_anno_df.anno_string)

#         menr_ctx[cid].motif_enrichment.loc[motifid,'region_db'] =  ','.join(menr_ctx[cid].motif_hits['Database'][motifid])
#         print('processing motif ' + motifid + ' done')
        
#         #output bed with annotation
#         temp = region_motif_hits_anno_df.loc[:,'anno_string'].str.split(',',expand = True)
#         temp.columns = ['peakid','chipseeker_anno','great_anno']
#         temp.loc[:,'annotation'] = temp.loc[:,'peakid'] + ',' + temp.loc[:,'chipseeker_anno'] + ',' + temp.loc[:,'great_anno']
        
#         temp1 = temp.loc[:,'peakid'].str.split(":|-",expand = True)
#         temp1.columns = ['chr','start','end']
#         temp.loc[:,['chr','start','end']] = temp1 # index and column must aggree
#         temp.loc[:,['chr','start','end','annotation']].to_csv(outdir + motifid + '.peaks.annotation.bed',sep='\t',header=False,index = False)

#         #get motif - tf relationship
#         if isinstance(menr_ctx[cid].motif_enrichment.loc[motifid,'Direct_annot'],str):
#             tf_anno_direct = menr_ctx[cid].motif_enrichment.loc[motifid,'Direct_annot'].split(', ')
#         else:
#             tf_anno_direct = []
        
#         if len(tf_anno_direct) != 0:
#             for j in tf_anno_direct:
#                 if j not in tf_motif_map.keys():
#                     tf_motif_map[j] = [motifid]
#                 else:
#                     tf_motif_map[j].append(motifid)
        
    
#     ###write tf-motif to file for this cluster
    
#     tf_motif_map_df = pd.DataFrame(columns = ['tf','motif'])
#     for k in tf_motif_map.keys():
#         v = ','.join(tf_motif_map[k])
#         tf_motif_map_df.loc[k,'tf'] = k
#         tf_motif_map_df.loc[k,'motif'] = v
    
#     print('#write tf-motif table file for cluster ' + cid )
#     tf_motif_map_df.to_csv(outdir  + 'tf-motif-table.txt',sep='\t',index=False)
    
#     ###write to file
#     print('#write result file for cluster ' + cid )
#     menr_ctx[cid].motif_enrichment.to_csv(work_dir + '/motif_enrichment.' + cid + '.txt',sep='\t')



# #save new menr_ctx


# with open(work_dir + '/menr_ctx.modify_motf_enrichment.dill','wb') as fh:
#     dill.dump(menr_ctx,fh)



    
    
    
# ###4 TF level cistrome: TF cistrome peaks (combine similar annotated motif to the same TF)

# #verify motif to TF??
# motif2tf = pd.read_csv('/home/mjwang/dataex/cisTargetDb_v10/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl',sep='\t',header=0,quoting= False)
# 253096 x 13 ##motif2tf annotation table



# # menr_ctx[cid].cistromes #combine of motif hits ?https://pycistarget.readthedocs.io/en/latest/tools.html
# # menr_ctx[cid].cistromes['Region_set'].keys() #TF gene name dict
# # menr_ctx[cid].cistromes['Region_set'].keys() == menr_ctx[cid].cistromes['Database'].keys() #TF gene name identical but length diff



# # menr_ctx[cid].cistromes['Region_set']['CEBPB_extended_(619r)'] #combine all CEBPB motifs to TF cistrome 
# # menr_ctx[cid].cistromes['Region_set']['GCM1_(3993r)']


# # pr.PyRanges(region_names_to_coordinates(menr_ctx[cid].cistromes['Region_set']['CEBPB_extended_(619r)'])).to_bed('test.bed')




# # ###add ChIPSeeker and GREAT precalculated peak annotation to tf-cistrome

# # menr_ctx[cid].cistromes['Region_set'].keys()

# # region_cistrome = menr_ctx[cid].cistromes['Region_set']['CEBPB_extended_(619r)']
# # #619
# # region_cistrome = menr_ctx[cid].cistromes['Region_set']['GCM1_(3993r)']
# # len(set(region_cistrome)) == len(region_cistrome) #619, no duplication


# # len(set(region_cistrome).intersection( set(peakAnno.index) ) ) == len(region_cistrome)
# # #619 #use peakAnno_merge_chipseekergreat
# # #432 #use peakAnno_merge_chipseekergreat_dar_highconf



# # region_cistrome_anno_df = peakAnno.loc[region_cistrome,:]

# # all(region_cistrome_anno_df.index == region_cistrome) #True


# # region_cistrome_anno_df.columns
# # # 'seqnames', 'start', 'end', 'width', 'strand', 'annotation', 'geneChr',
# # #        'geneStart', 'geneEnd', 'geneLength', 'geneStrand', 'geneId',
# # #        'transcriptId', 'distanceToTSS', 'ENSEMBL', 'SYMBOL', 'GENENAME',
# # #        'gene', 'site_name', 'gene1', 'dist1', 'gene2', 'dist2'
    
# # gene = [i for i in region_cistrome_anno_df.gene.value_counts().sort_index().index] #chipseeker assigned gene
# # gene1 = [i for i in region_cistrome_anno_df.gene1.value_counts().sort_index().index] #great assigned gene1
# # gene2 = [i for i in region_cistrome_anno_df.gene2.value_counts().sort_index().index] #great assigned gene2




# #for cid in ['cluster5']:
# for cid in menr_ctx.keys():
#     print('processing cluster ' + cid)
    
#     tf_cistrome_table = pd.DataFrame(columns=['TF','cistrome_id','motif_group','region','region_annotation','region_db'])
    
#     #make bed dir for this cluster
#     if not os.path.isdir(work_dir+'/cistrome_beds/' + cid):
#         #os.mkdir('result_pycisTarget/cistrome_beds/' + cid )
#         os.makedirs(work_dir + '/cistrome_beds/' + cid )
        
#     outdir = work_dir + '/cistrome_beds/' + cid + '/'
    
#     #menr_ctx[cid].cistromes.keys() #'Database', 'Region_set'
#     cistromeids = menr_ctx[cid].cistromes['Region_set'].keys() #277
#     cistromeids_db = menr_ctx[cid].cistromes['Database'].keys() #277
#     if not [i.split('_')[0] for i in cistromeids] == [i.split('_')[0] for i in cistromeids_db]: #True
#         print('cistromeids ne cistromeids_db')
#         sys.exit()
    

#     for cistromeid in cistromeids:
#     #for cistromeid in ['GCM1_(3993r)']:
#         print('processing tf cistrome ' + cistromeid)
 
#         #cistromeid = 'GCM1_(3993r)'
#         tfid = cistromeid.split('_')[0]
#         cistromeid_db = list(cistromeids_db)[list(cistromeids).index(cistromeid)]
#         tf_cistrome_table.loc[cistromeid,'TF'] = tfid
#         tf_cistrome_table.loc[cistromeid,'cistrome_id'] = cistromeid
#         tf_cistrome_table.loc[cistromeid,'motif_group'] = 'unknown'
        
#         print('tf ' + tfid + ' region_set ' + cistromeid + ' match to database metaregion ' + cistromeid_db)
        
#         region_cistrome = menr_ctx[cid].cistromes['Region_set'][cistromeid]
#         region_cistrome_db = menr_ctx[cid].cistromes['Database'][cistromeid_db]

#         tf_cistrome_table.loc[cistromeid,'region'] = ','.join(region_cistrome)
#         tf_cistrome_table.loc[cistromeid,'region_db'] = ','.join(region_cistrome_db)
        
        
#         if not len(set(region_cistrome)) == len(region_cistrome): #no duplication
#             print('duplication in region_cistrome')
#             sys.exit()

#         if not len(set(region_cistrome_db)) == len(region_cistrome_db): #no duplication
#             print('duplication in region_cistrome_db')
#             sys.exit()

#         if not pd.Series(region_cistrome).isin(dar_full.peakid).all(): #True,is dar region
#             print('region_cistrome must in dar_full region')
#             sys.exit()


#         ##output bed
#         region_cistrome_df = region_names_to_coordinates(region_cistrome)
#         region_cistrome_df.loc[:,'peakid'] = region_cistrome_df.loc[:,'Chromosome'] + ":" + region_cistrome_df.loc[:,'Start'].map(str) + '-' + region_cistrome_df.loc[:,'End'].map(str)
#         region_cistrome_df.to_csv(outdir + cistromeid + '.peaks.bed',sep='\t',header=False,index = False)
        
#         region_cistrome_db_df = region_names_to_coordinates(region_cistrome_db)
#         region_cistrome_db_df.loc[:,'peakid'] = region_cistrome_db_df.loc[:,'Chromosome'] + ":" + region_cistrome_db_df.loc[:,'Start'].map(str) + '-' + region_cistrome_db_df.loc[:,'End'].map(str)
#         region_cistrome_db_df.to_csv(outdir + cistromeid+ '_'+cistromeid_db + '.metaregion.bed',sep='\t',header=False,index = False)

#          ##add peak annotation for region_cistrome   
#         if not len(set(region_cistrome).intersection( set(peakAnno.index) ) ) == len(region_cistrome): #all has annotation in peakAnno table
#             print('not all peak in region_motif_hits has annotation,will filter regions without annotation')
#             print ( str(len(region_cistrome) - len(set(region_cistrome).intersection( set(peakAnno.index) ) )) + ' has no annotaton, filtering...' )
            
#             region_cistrome = [v for(v,i) in zip(region_cistrome ,pd.Series(region_cistrome).isin(peakAnno.index).to_list()) if i]
            
#             #sys.exit()

#         region_cistrome_anno_df = peakAnno.loc[region_cistrome,:]

#         if not all(region_cistrome_anno_df.index == region_cistrome): #True
#             print('not all peak in region_cistrome has annotation')
#             sys.exit()

#         ##region_motif_hits_anno_df.columns
#         # 'seqnames', 'start', 'end', 'width', 'strand', 'annotation', 'geneChr',
#         #        'geneStart', 'geneEnd', 'geneLength', 'geneStrand', 'geneId',
#         #        'transcriptId', 'distanceToTSS', 'ENSEMBL', 'SYMBOL', 'GENENAME',
#         #        'gene', 'site_name', 'gene1', 'dist1', 'gene2', 'dist2'

#         if not (region_cistrome_anno_df.gene == region_cistrome_anno_df.SYMBOL).all(): #True
#             print('in region_cistrome_anno_df, gene column != SYNBOL')
#             sys.exit()

#         ##add a formated annotation string
#         region_cistrome_anno_df.loc[:,'region_set'] = region_cistrome_anno_df.index.to_list()

#         region_cistrome_anno_df.loc[:,'annotype'] = region_cistrome_anno_df.annotation.str.split(' ',expand = True).iloc[:,0]

#         #region,chipseeker gene|dist|type,great gene1|dist1|gene2|dist2
#         region_cistrome_anno_df.loc[:,'anno_string'] = region_cistrome_anno_df.loc[:,'region_set'] + "," + region_cistrome_anno_df.loc[:,'gene'] + '|' + region_cistrome_anno_df.loc[:,'distanceToTSS'].map(str) + '|' + region_cistrome_anno_df.loc[:,'annotype'] + ',' + region_cistrome_anno_df.loc[:,'gene1'] + '|' + region_cistrome_anno_df.loc[:,'dist1'].map(str) + '|' + region_cistrome_anno_df.loc[:,'gene2'] + '|' + region_cistrome_anno_df.loc[:,'dist2'].map(str)

        
#         tf_cistrome_table.loc[cistromeid,'region_annotation'] = ' '.join(region_cistrome_anno_df.loc[:,'anno_string'] )
        
        
#         #output bed with annotation
#         temp = region_cistrome_anno_df.loc[:,'anno_string'].str.split(',',expand = True)
#         temp.columns = ['peakid','chipseeker_anno','great_anno']
#         temp.loc[:,'annotation'] = temp.loc[:,'peakid'] + ',' + temp.loc[:,'chipseeker_anno'] + ',' + temp.loc[:,'great_anno']
        
#         temp1 = temp.loc[:,'peakid'].str.split(":|-",expand = True)
#         temp1.columns = ['chr','start','end']
#         temp.loc[:,['chr','start','end']] = temp1 # index and column must aggree
#         temp.loc[:,['chr','start','end','annotation']].to_csv(outdir + cistromeid + '.peaks.annotation.bed',sep='\t',header=False,index = False)

#     ###write to file
#     print('#write result file for cluster ' + cid )
#     tf_cistrome_table.to_csv(work_dir + '/tfcistrome_enrichment.' + cid + '.txt',sep='\t',index = False)

    


# ##1 output tf-cistrome-peaks-annotation
# ##2 output a summary table file by cluster
    
# ##output tf cistrome result done
    
    
    
    
    
    
    
    
    
    
# ####run homer tf motif enrichment (problem with run_homer)

# homer_path = '/home/mjwang/anaconda3/envs/r413/share/homer/bin/'#'/home/mjwang/anaconda3/envs/myenv_ori/bin/'
# #need /home/mjwang/anaconda3/envs/r413/share/homer/bin/../motifs/extras/motifTable.txt


# #meme_collection_path = '/sda/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/chromVAR_TF_specific/homer_motif_clustering/known.motif.meme'
# #1006 motifs from all species, with redundancy
# ##reformat from /home/mjwang/anaconda3/envs/myenv_ori/share/homer/data/knownTFs/known.motifs

# meme_collection_path = '/home/mjwang/dataex/cisTargetDb_v10/motif2tf/v10nr_clust_public/homer_meme_from_cb/tf_motif_homer.meme'
# #homer_meme_from_cb from singleton dir
# meme_path = '/home/mjwang/anaconda3/envs/pacbio/bin/' #for denoe motif ?


# #path_to_motif_annotations = '/home/mjwang/dataex/cisTargetDb_v10/motif2tf/motif.homer.addH.tbl'
# path_to_motif_annotations = '/home/mjwang/dataex/cisTargetDb_v10/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl' #homer id stored as source_name in motif annotation table

# work_dir = 'result_pycisTarget_pairwise_c5c3_runhomer'
# #work_dir = 'result_pycisTarget_pairwise_c5c3_runhomer_new'
# outdir = work_dir + '/homer'



# #can not found SCENIC+ tf meme files
# ## use homer meme file instead ?

# homer_dict=run_homer(homer_path,
#                      region_set_pr_dict['dar_regression_pairwise_c5c3'],
#                      outdir,
#                      genome = 'hg38',
#                      size='given', #or will find motif center +/- nbp
#                      mask=False, #homer use smsk seq or not
#                      denovo=False, #off
#                      length='8,10,12',
#                      n_cpu=1,
#                      meme_path = None, #skip annotation of de novo motifs
#                      meme_collection_path = None, #skip annotation of de novo motifs
#                      annotation_version = 'v10_nrclust',
#                      path_to_motif_annotations = path_to_motif_annotations,
#                      #path_to_motif_annotations = '/staging/leuven/stg_00002/lcb/icistarget/data/motif2tf_project/motif_to_tf_db_data/snapshots/motifs-v10-nr.mgi-m0.00001-o0.0.tbl',
#                      cistrome_annotation = ['Direct_annot'],
#                      #cistrome_annotation = ['Direct_annot', 'Orthology_annot'],
#                      #_temp_dir='/scratch/leuven/313/vsc31305/ray_spill'
#                     )



# ####run homer commands info


# (homer_ray pid=165030) 2023-04-21 20:13:10,227 Homer        INFO     Running cluster3
# (homer_ray pid=165030) 2023-04-21 20:13:10,227 Homer        INFO     Running Homer for cluster3 with /home/mjwang/anaconda3/envs/r413/share/homer/bin/findMotifsGenome.pl result_pycisTarget_pairwise_c5c3_runhomer/homer/regions_bed/cluster3.bed hg38 result_pycisTarget_pairwise_c5c3_runhomer/homercluster3 -preparsedDir result_pycisTarget_pairwise_c5c3_runhomer/homercluster3 -size given -len 8,10,12 -nomotif -keepFiles
# (homer_ray pid=165031) 2023-04-21 20:13:10,250 Homer        INFO     Running cluster5
# (homer_ray pid=165031) 2023-04-21 20:13:10,251 Homer        INFO     Running Homer for cluster5 with /home/mjwang/anaconda3/envs/r413/share/homer/bin/findMotifsGenome.pl result_pycisTarget_pairwise_c5c3_runhomer/homer/regions_bed/cluster5.bed hg38 result_pycisTarget_pairwise_c5c3_runhomer/homercluster5 -preparsedDir result_pycisTarget_pairwise_c5c3_runhomer/homercluster5 -size given -len 8,10,12 -nomotif -keepFiles


# # ####running in command line
# # ! /home/mjwang/anaconda3/envs/r413/bin/findMotifsGenome.pl result_pycisTarget_pairwise_c5c3_runhomer/homer/regions_bed/cluster3.bed hg38 result_pycisTarget_pairwise_c5c3_runhomer/homercluster3 -preparsedDir result_pycisTarget_pairwise_c5c3_runhomer/homercluster3 -size given -len 8,10,12 -mask -nomotif -keepFiles

# # ! /home/mjwang/anaconda3/envs/r413/bin/findMotifsGenome.pl result_pycisTarget_pairwise_c5c3_runhomer/homer/regions_bed/cluster5.bed hg38 result_pycisTarget_pairwise_c5c3_runhomer/homercluster5 -preparsedDir result_pycisTarget_pairwise_c5c3_runhomer/homercluster5 -size given -len 8,10,12 -mask -nomotif -keepFiles

# # #success


# ##SnapATAC runHomer parameters:
# -len 10 -size 300 -S 2 -p 5 -cache 100 -fdr 5
# #-S <#> (Number of motifs to optimize, default: 25)
# #-p threads 
# #-len motif length (default 8,10,12)




# known_motifs = pd.read_csv('result_pycisTarget_pairwise_c5c3_runhomer/homercluster5/knownResults.txt', sep='\t')
# ctx_motif_annotation = load_motif_annotations('hg38',
#                      version = 'v10',
#                      fname = path_to_motif_annotations,
#                      motif_similarity_fdr= 0.1,
#                      orthologous_identity_threshold=0.1)

# add_motif_annotation_homer()


# find_motif_hits(n_cpu=1)
# ptint("Getting cistromes for " + self.name) 
# self.get_cistromes(self.cistrome_annotation)

# ##problem here 
        

        
    

# ##save run_homer result table

# with open(work_dir + '/enr_runhomer.dill','wb') as fh:
#     dill.dump(homer_dict,fh)


    

    
# ##Exploring Homer results


# #see html table
# homer_results(homer_dict, 'cluster5', results='known')
# homer_results(homer_dict, 'cluster3', results='known')

# #get cistrome region




