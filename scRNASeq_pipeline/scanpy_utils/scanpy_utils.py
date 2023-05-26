import scanpy as sc
import pandas as pd
import anndata
import os
import numpy as np
import seaborn as sb
from functools import reduce
from bioinfokit import analys, visuz

# process/readin raw count files
def check_gene_file_format(file_pres, version_10x):
    '''
    This function makes sure feature/gene file has 3 columns.
    * file_pres: prefix of 10X data files.
    * version: 2 (file ends with genes.tsv) or 3 (file ends with features.tsv.gz).
    '''
    for pre in file_pres:
        if version_10x == 2:
            feature = 'genes.tsv'
            pref_fn = f'{pre}{feature}'
            pref_df = pd.read_csv(pref_fn, sep='\t', header=None)
        elif version_10x == 3:
            feature = 'features.tsv.gz'
            pref_fn = f'{pre}{feature}'
            pref_df = pd.read_csv(pref_fn, sep='\t', header=None, compression='gzip')
        # make column unique
        pref_df[0] = anndata.utils.make_index_unique(pd.Index(pref_df[0]))
        # add column if missing
        if pref_df.shape[1] == 1:
            pref_df[1] = pref_df[0].str.upper()
        if pref_df.shape[1] == 2:
            pref_df[2] = 'Gene Expression'
        if version_10x == 2:
            pref_df.to_csv(pref_fn,sep='\t',index=False,header=None)
        elif version_10x == 3:
            pref_df.to_csv(pref_fn,sep='\t',index=False,header=None,compression='gzip')


def read_10x_raw_files(file_pres, file_tree, **kwargs):
    '''read 10x raw data in format of 3 files
    '''
    # read data using scanpy
    adata_list = []
    for pre in file_pres:
        if file_tree == 'in_1_folder':
            adata_tmp = sc.read_10x_mtx(os.path.dirname(pre),prefix=os.path.basename(pre), cache=True)
        elif file_tree == 'in_n_folder':
            adata_tmp = sc.read_10x_mtx(pre)
        adata_tmp.var_names_make_unique()
        adata_tmp.obs.index = adata_tmp.obs.index.map(lambda x: x.split('-')[0])
        adata_tmp.obs['sample'] = f'{os.path.basename(pre[:-1])}'
        adata_list.append(adata_tmp)
    adata = sc.concat(adata_list)
    adata.obs_names_make_unique()
    adata.var = adata_tmp.var
    del adata_list
    return adata
            

def read_10x_h5_files(h5_fns, **kwargs):
    ''' read 10x h5 data and merge into one anndata object
    '''
    adata_list = []
    samples = []
    for h5 in h5_fns:
        adata_tmp = sc.read_10x_h5(h5, gex_only=False)
        adata_tmp.var_names_make_unique()
        adata_tmp.obs.index = adata_tmp.obs.index.map(lambda x: x.split('-')[0])
        sample = os.path.splitext(os.path.basename(h5))[0]
        samples.append(sample)
        adata_tmp.obs['sample'] = sample
        adata_list.append(adata_tmp)
    adata = sc.concat(adata_list,keys=samples,index_unique='-')
    adata.obs_names_make_unique()
    adata.var = adata_tmp.var
    del adata_list
    return adata


# preprocessing qc
def qc_plots(adata):
    sc.pp.calculate_qc_metrics(adata,inplace=True)
    adata.obs['n_genes'] = (adata.X > 0).sum(1)
    mt_gene_mask = [gene.startswith('MT-') for gene in adata.var_names]
    ribo_gene_mask = [gene[:3] in ['RPL','RPS'] for gene in adata.var_names]
    adata.obs['mt_frac'] = np.array(adata.X[:, mt_gene_mask].sum(1)).flatten() / adata.obs['total_counts']
    adata.obs['ribo_frac'] = np.array(adata.X[:, ribo_gene_mask].sum(1)).flatten() / adata.obs['total_counts']
    # get largest gene
    adata_pp = adata[:, adata.var.index != 'MALAT1'].copy()
    adata.obs['largest_gene'] = adata_pp.var.index[np.argmax(adata_pp.X, axis=1).reshape(-1,1)]
    adata.obs['largest_count_frac'] = np.max(adata_pp.X.A, axis=1) / adata.obs['total_counts']
    del adata_pp
    # Sample quality plots
    t1 = sc.pl.violin(adata, 'total_counts', groupby='sample', size=1, log=True,cut=0,rotation=90, figsize=(2,2))
    t2 = sc.pl.violin(adata, ['mt_frac','n_genes','ribo_frac', 'largest_count_frac'], groupby='sample', rotation=90)
    # ncount and ngenes summary, you can change the axis region to zoom in
    p1 = sc.pl.scatter(adata,'total_counts','n_genes', color='mt_frac')
    p2 = sc.pl.scatter(adata[adata.obs['total_counts']<5000], 'total_counts', 'n_genes', color='mt_frac')
    #Thresholding decision: counts, use height and aspect to change figure size
    p3 = sb.displot(adata.obs['total_counts'], kde=False,height=5, aspect=1.2)
    p4 = sb.displot(adata.obs['total_counts'][adata.obs['total_counts']<10000], kde=False, bins=60,height=5,aspect=1.2)
    p5 = sb.displot(adata.obs['total_counts'][adata.obs['total_counts']>10000], kde=False, bins=60,height=5,aspect=1.2)
    #Thresholding decision: genes
    p5 = sb.displot(adata.obs['n_genes'], kde=False,height=5, aspect=1.2)
    p6 = sb.displot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=60,height=5,aspect=1.2)
    


# cell type annotation
def scsa_cell_anno(adata, rank_key, clust_col, species='Human',tissue='Blood', tmp_dir='/tmp'):
    '''
    * adata: scanpy anndata object
    * rank_key: key name in uns for ranking test
    * clust_col: column name in obs dataframe
    You can run the following command to check which species and tissue is avaialble in the database
            !python /opt/SCSA/SCSA.py -d /opt/SCSA/whole.db -i test -l
    '''
    # 1. format files
    result = adata.uns[rank_key]
    groups = result['names'].dtype.names
    dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
    cluster_fn = f'{tmp_dir}/{rank_key}.tsv'
    dat.to_csv(cluster_fn)
    # 2. run scsa
    cell_anno_db_fn = '/opt/SCSA/whole.db'
    scsa_anno_fn = f'{tmp_dir}/scsa_anno_{rank_key}.txt'
    cmd = f'python3 /opt/SCSA/SCSA.py -d {cell_anno_db_fn} -i {cluster_fn} -s scanpy -E -f1.5 \
       -p 0.01 -o {scsa_anno_fn} -m txt -b -g {species} -k {tissue}'
    log = os.popen(cmd).read()
    print(log)
    # 3. build inital cell type annotation
    n_clust = adata.obs[clust_col].unique().shape[0]
    init_ct = log.split('\n')[-(n_clust+1):-1]
    init_celltype = {}
    for c in init_ct:
        c = c.replace('\'','')
        c = c[1:-1].split(',')
        init_celltype[c[0]] = c[2][1:]
    return init_celltype


def get_cluster_markers(adata, key):
    '''
    this function gets biomarkers for each cluster, with group name append before the columns
    '''
    result = adata.uns[key]
    groups = result['names'].dtype.names
    de_df = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals','pvals_adj']})
    if 'pts' in result:
        pts_df = result['pts']
        pts_df.columns = ['pts_' + c for c in pts_df.columns]
        de_df = pd.merge(de_df, pts_df, how='left', left_on=groups[0] + '_names', right_index=True)
    return de_df



# DE analysis
def clust_de_wilcoxon(adata, clust_col, condition, ref, pre_key):
    '''
    * clust_col: column name to indicate the cluster of cells
    * condition: column name in obs to indicate all conditions.
    * ref: reference group as control in DE analysis
    * pre_key: key added to the uns to indicate what comparisons were did
    '''
    clusts = adata.obs[clust_col].cat.categories
    de_dfs = []
    for clust in clusts:
        cl = adata[adata.obs.query(f'`{clust_col}` == @clust').index,:]
        sc.pp.filter_genes(cl, min_cells = cl.shape[0] * 0.1)
        key = f'{pre_key}_clust{clust}' # f'IgA+_Bcell_igan_vs_ctrl'
        sc.tl.rank_genes_groups(cl, condition, reference=ref, method='wilcoxon', key_added=key, pts=True, use_raw=False)
        adata.uns[key] = cl.uns[key].copy()
        # extract results
        result = cl.uns[key]
        groups = result['names'].dtype.names
        de_df = pd.DataFrame({f'clust{clust}' + '_' + k: result[k][group] for group in groups for k in ['names', 'logfoldchanges','scores','pvals','pvals_adj']})
        pts_df = result['pts'].copy()
        pts_df.columns = [f'clust{clust}_{c}_pts' for c in pts_df.columns]
        de_df = pd.merge(de_df, pts_df, how='left', left_on=f'clust{clust}_names', right_index=True)
        de_dfs.append(de_df)
    # merge data using index, keep ranking of each cluster
    concat_df = pd.concat(de_dfs, axis=1)
    for de in de_dfs:
        de.index = de.iloc[:,0]
        del de[de.columns[0]]
    # merge data using the same index
    # merge data using the cluster 0 gene ranking
    merge_df = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how='left'), de_dfs)
    return concat_df, merge_df


def volcano(df, gene_col, lfc_col, qval_col, lfc, qval, out_path, title, topN=0):
    '''
    * gene_col: column name indicating gene name
    * lfc_col: column name indicating log2foldchange
    * qval_col: column name indicating qvalue
    * lfc: threshold of log2foldchange
    * qval: threshold of qvalue
    * out_path: path to save the output
    '''
    os.chdir(out_path)
    df = df[~df[qval_col].isna()].reset_index(drop=True)
    df = df[~df[lfc_col].isna()].reset_index(drop=True)
    
    up = df.query(f'{qval_col} < {qval} and {lfc_col} > {lfc}')[gene_col].tolist()
    dn = df.query(f'{qval_col} < {qval} and {lfc_col} < -{lfc}')[gene_col].tolist()
    if topN == 0:
        plot_genes = up + dn
    else:
        plot_genes = up[:topN] + dn[:topN]
    

    visuz.GeneExpression.volcano(df=df, lfc=lfc_col, pv=qval_col,
            color=("#E10600FF","grey","#00239CFF"),figname=title,
            sign_line=True, lfc_thr=(lfc, lfc),geneid=gene_col,
                                genenames=tuple(plot_genes))
    
    
# plot
def plot_rank_gene_de_num(merge_df, title, ncol_clust, lfc=1, padj=0.05):
    '''
    plot DE gene number for each cluster using results extracted from the results
    * merge_df: each cluster has ncol_clust columns.
    * title: title of the plot
    * ncol_clust: number of columns for each subcluster in the dataframe. 
            eg: cluster 0 has the following columns: [clust0_logfoldchanges	clust0_scores	clust0_pvals	clust0_pvals_adj	clust0_IgAN_pts	clust0_control_pts],
            there are 6 columns, then the paramter is set as 6.
    
    return: 
        * ax: object of plot
        * other_dict: dictionary {cluster: list of non-DE genes}
    '''
    up_n = []; dn_n = []
    up_dict = {}; dn_dict = {}; other_dict = {}
    for i in [str(i) for i in range(merge_df.shape[1]//ncol_clust)]:
        up = merge_df.query(f'clust{i}_logfoldchanges > {lfc} and clust{i}_pvals_adj < {padj}')
        dn = merge_df.query(f'clust{i}_logfoldchanges < -{lfc} and clust{i}_pvals_adj < {padj}')
        up_n.append(up.shape[0])
        dn_n.append(dn.shape[0])
        up_dict[i] = up.index
        dn_dict[i] = dn.index
        other_dict[i] = list(set(merge_df.index) - set(up.index) - set(dn.index))
    df = pd.DataFrame({'up':up_n, 'down':dn_n})
    ax = df.plot(kind='bar')
    ax.set_title(title)
    ax.set_xlabel('cluster')
    ax.set_ylabel('number of DE genes')
    return ax, (up_dict, dn_dict, other_dict)