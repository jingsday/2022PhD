Two approaches attemped. Please also find pyscenic usage at the end of this document.
1. Pyscenic command lines:
nohup pyscenic ctx /home/jyang/projects/GBM/specimens/individuals/SF2777/adj.csv --annotations_fname /home/jyang/projects/GBM/resource/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
--expression_mtx_fname /home/jyang/projects/GBM/specimens/individuals/SF2777/SF2777tumor.loom --output /home/jyang/projects/GBM/specimens/tumor_merged/SF2777/reg.csv --mask_dropouts
--num_workers 4 /home/jyang/projects/GBM/resource/hg19-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather
/home/jyang/projects/GBM/resource/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather
/home/jyang/projects/GBM/resource/hg19-tss-centered-5kb-10species.mc9nr.genes_vs_motifs.rankings.feather
  

3. Command generated using ipynb(please refer to the notebook scenic_step2.ipynb under hpcrun.




usage: pyscenic ctx [-h] [-o OUTPUT] [-n] [--chunk_size CHUNK_SIZE]
                    [--mode {custom_multiprocessing,dask_multiprocessing,dask_cluster}] [-a] [-t]
                    [--rank_threshold RANK_THRESHOLD] [--auc_threshold AUC_THRESHOLD]
                    [--nes_threshold NES_THRESHOLD] [--min_orthologous_identity MIN_ORTHOLOGOUS_IDENTITY]
                    [--max_similarity_fdr MAX_SIMILARITY_FDR] --annotations_fname ANNOTATIONS_FNAME
                    [--num_workers NUM_WORKERS] [--client_or_address CLIENT_OR_ADDRESS]
                    [--thresholds THRESHOLDS [THRESHOLDS ...]]
                    [--top_n_targets TOP_N_TARGETS [TOP_N_TARGETS ...]]
                    [--top_n_regulators TOP_N_REGULATORS [TOP_N_REGULATORS ...]] [--min_genes MIN_GENES]
                    [--expression_mtx_fname EXPRESSION_MTX_FNAME] [--mask_dropouts]
                    [--cell_id_attribute CELL_ID_ATTRIBUTE] [--gene_attribute GENE_ATTRIBUTE] [--sparse]
                    module_fname database_fname [database_fname ...]
