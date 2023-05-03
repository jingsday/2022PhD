#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import csv
import gzip
import anndata as ad
from pathlib import Path
import glob

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)


# In[2]:


def sce(sample):

    #set corresponding folders which the results of each specimen are written to
  

    f_db_glob = "/home/jing/scratch/resource/*feather"
    f_db_names = ' '.join( glob.glob(f_db_glob) )

    # motif databases
    f_motif_path = "/home/jing/scratch/resource/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
    f_tfs = "/home/jing/scratch/resource/allTFs_hg38.txt"
    
    output_folder = "/home/jing/projects/samples/"
    output_path = os.path.join(output_folder, sample)
    
    f_loom_path_scenic = f"{output_path}/filtered_scenic.loom"
    os.chdir(output_path)

    # path to  loom file

    f_anndata_path = f"{output_path}/anndata.h5ad"
    f_pyscenic_output = f"{output_path}/output.loom"
    f_final_loom = f"{output_path}/scenic_integrated-output.loom"
    
    #I
    get_ipython().system('pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 8')
    adjacencies = pd.read_csv("adj.csv", index_col=False, sep='\t')

    adjacencies.head()
    
    #II
   
    get_ipython().system('pyscenic ctx adj.csv         {f_db_names}         --annotations_fname {f_motif_path}         --expression_mtx_fname {f_loom_path_scenic}         --output reg.csv         --mask_dropouts         --num_workers 8')
    
    
   #III     
    get_ipython().system('pyscenic aucell         {f_loom_path_scenic}         reg.csv         --output {f_pyscenic_output}         --num_workers 8')
    
    lf = lp.connect(f_pyscenic_output, mode='r+', validate=False)

    auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
    lf.close()
    auc_mtx.to_pickle('auc_mtx.pickle')
    
    
    


# In[3]:


sce('SF3448')


# In[ ]:




