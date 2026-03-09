import pandas as pd
import numpy as np
from samap.analysis import get_mapping_scores
from samap.analysis import (sankey_plot)
import pickle
import holoviews as hv
from holoviews import dim


#df= pd.read_csv('test.txt', sep='\t', header=0, index_col=0)

samap_pkl = f'/HOMA/samap.pkl'
with open(samap_pkl, 'rb') as f:
    sm1 = pickle.load(f)
keys = {i:'seurat_clusters' for i in sm1.ids}
D1, M1 = get_mapping_scores(sm1, keys, n_top=1000)

samap_pkl = f'/MAMU/samap.pkl'
with open(samap_pkl, 'rb') as f:
    sm2 = pickle.load(f)
keys = {i:'seurat_clusters' for i in sm2.ids}
D2, M2 = get_mapping_scores(sm2, keys, n_top=1000)

samap_pkl = f'/MUTU/samap.pkl'
with open(samap_pkl, 'rb') as f:
    sm7 = pickle.load(f)
keys = {i:'seurat_clusters' for i in sm7.ids}
D7, M7 = get_mapping_scores(sm7, keys, n_top=1000)

samap_pkl = f'/TUDO/samap.pkl'
with open(samap_pkl, 'rb') as f:
    sm3 = pickle.load(f)
keys = {i:'seurat_clusters' for i in sm3.ids}
D3, M3 = get_mapping_scores(sm3, keys, n_top=1000)

samap_pkl = f'/DOOV/samap.pkl'
with open(samap_pkl, 'rb') as f:
    sm4 = pickle.load(f)
keys = {i:'seurat_clusters' for i in sm4.ids}
D4, M4 = get_mapping_scores(sm4, keys, n_top=1000)

samap_pkl = f'/OVGO/samap.pkl'
with open(samap_pkl, 'rb') as f:
    sm5 = pickle.load(f)
keys = {i:'seurat_clusters' for i in sm5.ids}
D5, M5 = get_mapping_scores(sm5, keys, n_top=1000)

samap_pkl = f'/GOSU/samap.pkl'
with open(samap_pkl, 'rb') as f:
    sm6 = pickle.load(f)
keys = {i:'seurat_clusters' for i in sm6.ids}
D6, M6 = get_mapping_scores(sm6, keys, n_top=1000)


#nodes =list(set( M1.index.tolist() + M2.index.tolist() + M7.index.tolist() + M3.index.tolist() + M4.index.tolist() + M5.index.tolist() + M6.index.tolist()))
#na_df = pd.DataFrame(np.nan, index=nodes , columns=nodes )
#df = na_df.fillna(M1).fillna(M1.T).fillna(M2).fillna(M2.T).fillna(M7).fillna(M7.T).fillna(M3).fillna(M3.T).fillna(M4).fillna(M4.T).fillna(M5).fillna(M5.T).fillna(M6).fillna(M6.T).fillna(0)

M1.to_csv('HOMA.results.txt', sep='\t')
M2.to_csv('MAMU.results.txt', sep='\t')
M7.to_csv('MUTU.results.txt', sep='\t')
M3.to_csv('TUDO.results.txt', sep='\t')
M4.to_csv('DOOV.results.txt', sep='\t')
M5.to_csv('OVGO.results.txt', sep='\t')
M6.to_csv('GOSU.results.txt', sep='\t')


fig=sankey_plot(M1, align_thr=0.3)
import holoviews as hv
hv.extension('matplotlib')
hv.save(obj=fig, filename='HOMA', fmt='svg')

fig=sankey_plot(M2, align_thr=0.3)
import holoviews as hv
hv.extension('matplotlib')
hv.save(obj=fig, filename='MAMU', fmt='svg')

fig=sankey_plot(M7, align_thr=0.4)
import holoviews as hv
hv.extension('matplotlib')
hv.save(obj=fig, filename='MUTU', fmt='svg')

fig=sankey_plot(M3, align_thr=0.4)
import holoviews as hv
hv.extension('matplotlib')
hv.save(obj=fig, filename='TUDO', fmt='svg')

fig=sankey_plot(M4, align_thr=0.4)
import holoviews as hv
hv.extension('matplotlib')
hv.save(obj=fig, filename='DOOV', fmt='svg')

fig=sankey_plot(M5, align_thr=0.4)
import holoviews as hv
hv.extension('matplotlib')
hv.save(obj=fig, filename='OVGO', fmt='svg')

fig=sankey_plot(M6, align_thr=0.3)
import holoviews as hv
hv.extension('matplotlib')
hv.save(obj=fig, filename='GOSU', fmt='svg')


#import holoviews as hv
#hv.extension('matplotlib')
#hv.save(obj=fig, filename='7species1', fmt='svg')
