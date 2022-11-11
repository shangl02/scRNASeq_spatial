''' This file has modules for analyzing visium data
'''

import scanpy as sc
import squidpy as sq
import pandas as pd
import anndata as ad
import os
import glob



def get_adata_list(file_pres):
    ''' given file prefix of spaceranger result folder, readin the data and do qc analysis
    '''
    adata_list = []
    samples = []
    for pre in file_pres:
        adata_tmp = sc.read_visium(pre)
        adata_tmp.var_names_make_unique()
        adata_tmp.var.index = adata_tmp.var.index.str.upper()
        sample = f'{os.path.basename(pre)}'
        adata_tmp.obs['sample'] = sample
        # store raw data
        adata_tmp.layers['counts'] = adata_tmp.X.copy()
        # add cell ranger cluster
        cr_clst_fn = glob.glob(f'{pre}/analysis/clustering/*graph*/clusters.csv')[0]
        cr_clst_df = pd.read_csv(cr_clst_fn, header=0, index_col=0, names=['cr_cluster']).astype(str)
        adata_tmp.obs = pd.merge(adata_tmp.obs, cr_clst_df, left_index=True, right_index=True,how='left')
        samples.append(sample)
        # get metrics
        sc.pp.calculate_qc_metrics(adata_tmp,inplace=True)
        adata_tmp = adata_tmp[adata_tmp.obs['total_counts'] > 0]
        adata_list.append(adata_tmp)
    return adata_list, samples


def norm_hv_gene(adata_list):
    '''normalize data and find highly variable genes
    '''
    for ad in adata_list:
        sc.pp.normalize_total(ad, inplace=True)
        sc.pp.log1p(ad)
        # store raw data
        ad.raw = ad
        sc.pp.highly_variable_genes(ad, flavor="seurat", n_top_genes=2000, inplace=True)
        

def filter_spots_by_gene(adata, genes, log_threshold=0, rings=3):
    '''filter spots based on genes expressoin
    * log_threshold: usually raw data of anndata stores log normalized data.
    * rings: number of rings surrond the target spots, will not be considered control spots
    
    return: 
        * target spot indexes.
        * ring spot indexes
        * the rest control indexes
    '''
    sub = adata.copy()
    for g in genes:
        sub = sub[sub.raw[:,g].X > 0]
    all_cells = adata.obs.index.tolist()
    target_cells = sub.obs.index.tolist()
    target_index = [all_cells.index(c) for c in target_cells]
    # get ring spots
    sq.gr.spatial_neighbors(adata, n_rings=rings, coord_type="grid", n_neighs=6)
    _, idx = adata.obsp["spatial_connectivities"][target_index, :].nonzero()
    close_ring_idx = list(set(idx)) + target_index
    ring_cells = [all_cells[i] for i in close_ring_idx if all_cells[i] not in target_cells]
    norm_cells = list(set(all_cells) - set(ring_cells) - set(target_cells))
    
    return target_cells, ring_cells, norm_cells 


def use_geneid_as_index(adata_list):
    '''
    This function change ensembl gene id as index
    '''
    for ad in adata_list:
        ad.var['SYMBOL'] = ad.var_names
        ad.var.set_index('gene_ids', drop=True, inplace=True)
        
        
def rm_mt_genes(adata_list):
    '''
    This function remove mitochondria genes
    '''
    for ad in adata_list:
        # find mitochondria-encoded (MT) genes
        ad.var['MT_gene'] = [gene.startswith('MT-') for gene in ad.var.index]

        # remove MT genes for spatial mapping (keeping their counts in the object)
        ad.obsm['MT'] = ad[:, ad.var['MT_gene'].values].X.toarray()
        ad = ad[:, ~ad.var['MT_gene'].values]

    