import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")
warnings.filterwarnings("ignore", message=".*DataFrame is highly fragmented*")

import os
import sys
import re
from scipy.cluster.hierarchy import fcluster
os.environ['OPENBLAS_NUM_THREADS'] = '4'

if __name__ == '__main__':
    try:
        emb_data_path = sys.argv[1]
        genes_to_macrogenes_path = sys.argv[2]
        alignedterm_path = sys.argv[3]
        annot_map_path = sys.argv[4]
    except:
        print("Usage: python3 " + os.path.basename(__file__) + " emb_data_path genes_to_macrogenes.pkl alignedterm_path annot_map_path")
        sys.exit(1)

emb_basename = os.path.basename(emb_data_path)

import pickle
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

if emb_data_path.endswith('loom'):
    adata_emb = sc.read_loom(emb_data_path,validate=False)
elif emb_data_path.endswith('h5ad'):
    adata_emb = sc.read_h5ad(emb_data_path)

genes_to_macrogenes = pickle.load(open(genes_to_macrogenes_path, "rb"))
genes_to_macrogenes_df = pd.DataFrame(genes_to_macrogenes)

adata_emb.obs['OrigID'] = adata_emb.obs.index
adata_emb.obs_names_make_unique()
adata_emb.var_names_make_unique()

adata_emb.layers["counts"] = adata_emb.X.copy()

if annot_map_path != 'NA':
    annot_map = pd.read_table(annot_map_path, header=None)
    annot_map_dict = annot_map.set_index(0)[1].to_dict()
    adata_emb.obs['superterm'] = adata_emb.obs['labels'].apply(lambda x: x+"_"+annot_map_dict[x] if x in annot_map_dict else x)
else:
    adata_emb.obs['superterm'] = adata_emb.obs['labels']

for sp in adata_emb.obs['species'].unique():
    adata_emb.obs[sp] = adata_emb.obs['species'].apply(lambda x: sp if x == sp else np.nan)
    adata_emb.obs[sp+'_labels'] = np.where(adata_emb.obs['species'] == sp, adata_emb.obs['labels'].astype('str').replace(sp+"_"+sp, '', regex=True), np.nan)

sc.pp.pca(adata_emb, n_comps=50)
sc.pp.neighbors(adata_emb, n_neighbors=15)

sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(24, 6), fontsize=6, facecolor='white')
for cor in ['pearson', 'kendall', 'spearman']:
    #for linkage in ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']:
    for linkage in ['ward']:
        sc.tl.dendrogram(adata_emb, 'superterm', linkage_method=linkage, cor_method=cor, key_added=f"den.{cor}.{linkage}")
        sc.pl.dendrogram(adata_emb, 'superterm', dendrogram_key=f"den.{cor}.{linkage}", save="." + emb_basename + '.' + cor + '.' + linkage + '.pdf')

adata_emb.uns['dendrogram_superterm'] = adata_emb.uns['den.pearson.ward']

if alignedterm_path != 'NA':
    aligned_df = pd.read_table(alignedterm_path, header=None)
    aligned_dict = aligned_df.set_index(0)[1].to_dict()
    adata_emb.obs['alignedterm'] = adata_emb.obs['labels'].apply(lambda x: aligned_dict[x])
else:
    aligned = fcluster(adata_emb.uns['den.pearson.ward']['linkage'], t=0.5, criterion='distance')
    aligned_dict = {value: f"aln_{aligned[i]}" for i, value in enumerate(sorted(adata_emb.uns['den.pearson.ward']['categories_ordered']))}
    adata_emb.obs['alignedterm'] = adata_emb.obs['superterm'].apply(lambda x: aligned_dict[x])

adata_emb.obs['alignedterm'] = adata_emb.obs['alignedterm'].astype('category')

aln_superterm_dict={}
for key, value in aligned_dict.items():
    aln_superterm_dict.setdefault(value, []).append(key)

macrogene_adata = sc.AnnData(adata_emb.obsm["macrogenes"])
macrogene_adata = macrogene_adata[:,genes_to_macrogenes_df.max(axis=1)>0.5]
macrogene_adata.X = macrogene_adata.X.astype(int)
macrogene_adata.layers["counts"] = macrogene_adata.X.copy()
macrogene_adata.obs = adata_emb.obs
macrogene_adata.uns = adata_emb.uns
macrogene_adata.groups_order, macrogene_adata.groups_masks = sc._utils.select_groups(macrogene_adata, 'all', 'labels')
macrogene_adata.pts = np.zeros((macrogene_adata.groups_masks.shape[0], macrogene_adata.shape[1]))
for imask, mask in enumerate(macrogene_adata.groups_masks):
    X_mask = macrogene_adata.X[mask]
    macrogene_adata.pts[imask] = (np.count_nonzero(X_mask, axis=0) / X_mask.shape[0]).round(3)
    macrogene_adata.pts[imask] = macrogene_adata.pts[imask].round(3)

sc.pp.normalize_total(macrogene_adata, target_sum=1e4)
macrogene_adata.layers["norm"] = macrogene_adata.X.copy()
sc.pp.log1p(macrogene_adata)
macrogene_adata.layers["log1p"] = macrogene_adata.X.copy()
#sc.pp.scale(macrogene_adata, max_value=10)
sc.tl.rank_genes_groups(macrogene_adata, 'alignedterm', method='wilcoxon', key_added = "wilcoxon_alignedterm")

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

sorted_alignedterm = sorted(list(set(macrogene_adata.obs['alignedterm'])), key=natural_sort_key)

os.system(f"rm -rf {emb_basename}.wilcoxon_alignedterm.p0.01.lfc0.25.csv")
macrodegs={}
for i in sorted_alignedterm:
    deg = sc.get.rank_genes_groups_df(macrogene_adata, group=i, key='wilcoxon_alignedterm', pval_cutoff=0.01, log2fc_min=0.25)
    deg['group'] = i
    deg['subgroup'] = deg.apply(lambda x: ':'.join(aln_superterm_dict[x.group]), axis=1)
    deg['pct'] = deg.apply(lambda x: ':'.join(str(macrogene_adata.pts[np.where(macrogene_adata.groups_order == str(i))[0][0]][macrogene_adata.var_names.tolist().index(x.names)]) for i in x.subgroup.split(':')), axis=1)
    deg['pct_min'] = deg.apply(lambda x: np.min(list(macrogene_adata.pts[np.where(macrogene_adata.groups_order == str(i))[0][0]][macrogene_adata.var_names.tolist().index(x.names)] for i in x.subgroup.split(':'))), axis=1)
    deg['pct_mean'] = deg.apply(lambda x: np.mean(list(macrogene_adata.pts[np.where(macrogene_adata.groups_order == str(i))[0][0]][macrogene_adata.var_names.tolist().index(x.names)] for i in x.subgroup.split(':'))), axis=1)
    deg_sorted = (
            deg.reset_index()
            .groupby("group")
            .apply(lambda x: x.sort_values(by=["scores", "pct_min"], ascending=False))
            .reset_index(drop=True)
    )
    marker_genes = list(deg_sorted.head(5)["names"])
    macrodegs[i] = marker_genes
    deg_sorted.to_csv(f"{emb_basename}.wilcoxon_alignedterm.p0.01.lfc0.25.csv", mode='a', index=False)

macrogene_adata.X = macrogene_adata.layers["log1p"].copy()
#macrogene_adata_quantile_99 = np.percentile(macrogene_adata.X, 99)
#sc.pl.rank_genes_groups_dotplot(macrogene_adata, var_names=macrodegs['aln_1'], vmax = macrogene_adata_quantile_99, groupby='superterm', key='wilcoxon_alignedterm', save=f".{emb_basename}.wilcoxon_alignedterm.p0.01.lfc0.25.aln_1.pdf")
macrogene_adata_quantile_995 = np.percentile(macrogene_adata.X, 99.5)
sc.pl.rank_genes_groups_dotplot(macrogene_adata, var_names=np.array(list(macrodegs.values())).flatten(), vmax = macrogene_adata_quantile_995, groupby='superterm',key='wilcoxon_alignedterm', save=f".{emb_basename}.wilcoxon_alignedterm.p0.01.lfc0.25.pdf")

sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))

sc.tl.umap(adata_emb)

sc.pl.umap(
        adata_emb,
        color=['species'] + adata_emb.obs['species'].unique().tolist() + [f"{sp}_labels" for sp in adata_emb.obs['species'].unique()],
        na_color="grey",
        ncols=2,
        size=2,
        wspace=1,
        legend_loc="on data",
        legend_fontsize="xx-small",
        legend_fontweight="normal",
        save="." + emb_basename + '.species.pdf',
        )

for feature in ['batch_labels', 'labels']:
    if feature in adata_emb.obs:
        if len(adata_emb.obs[feature].unique()) > 80:
            sc.pl.umap(adata_emb,
                    color=[feature],
                    frameon=False,
                    wspace=0.6,
                    palette=list(colors.CSS4_COLORS.values()),
                    save="." + emb_basename + '.' + feature + '.pdf',
                    )
        else:
            sc.pl.umap(adata_emb,
                    color=[feature],
                    frameon=False,
                    wspace=0.6,
                    save="." + emb_basename + '.' + feature + '.pdf',
                    )

os.system('mkdir -p ' + emb_basename + ".figures")

os.system('\\cp figures/*.' + emb_basename + '.*.pdf ' + emb_basename + ".figures" + '/')

