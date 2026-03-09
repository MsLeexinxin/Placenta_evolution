import pandas as pd
import numpy as np
from samap.analysis import get_mapping_scores
from samap.analysis import (sankey_plot)
import pickle

M1= pd.read_csv('HOMA.results.txt', sep='\t', header=0, index_col=0)
M2= pd.read_csv('MAMU.results.txt', sep='\t', header=0, index_col=0)
M3= pd.read_csv('MUTU.results.txt', sep='\t', header=0, index_col=0)
M4= pd.read_csv('TUDO.results.txt', sep='\t', header=0, index_col=0)
M5= pd.read_csv('DOOV.results.txt', sep='\t', header=0, index_col=0)
M6= pd.read_csv('OVGO.results.txt', sep='\t', header=0, index_col=0)
M7= pd.read_csv('GOSU.results.txt', sep='\t', header=0, index_col=0)

nodes =list(set(M1.index.tolist() + M2.index.tolist() + M3.index.tolist() + M4.index.tolist() + M5.index.tolist() + M6.index.tolist() + M7.index.tolist()))


na_df = pd.DataFrame(np.nan, index=nodes , columns=nodes )

df = na_df.fillna(M1).fillna(M1.T).fillna(M2).fillna(M2.T).fillna(M3).fillna(M3.T).fillna(M4).fillna(M4.T).fillna(M5).fillna(M5.T).fillna(M6).fillna(M6.T).fillna(M7).fillna(M7.T).fillna(0)

df.to_csv('HO_MA_MU_TU_DO_OV_GO_SU.merge.results.txt', sep='\t')
