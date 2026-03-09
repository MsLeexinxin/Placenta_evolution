import pandas as pd
import numpy as np
from samap.analysis import get_mapping_scores
from samap.analysis import (sankey_plot)
import pickle
import holoviews as hv
from holoviews import dim

df= pd.read_csv('HO_MA_MU_TU_DO_OV_GO_SU.merge.results.txt', sep='\t', header=0, index_col=0)


def sankey_plot(M,species_order=None,align_thr=0.1,**params):
    """Generate a sankey plot
    
    Parameters
    ----------
    M: pandas.DataFrame
        Mapping table output from `get_mapping_scores` (second output).

    align_thr: float, optional, default 0.1
        The alignment score threshold below which to remove cell type mappings.
    
    species_order: list, optional, default None
        Specify the order of species (left-to-right) in the sankey plot.
        For example, `species_order=['hu','le','ms']`.

    Keyword arguments
    -----------------
    Keyword arguments will be passed to `sankey.opts`.
    """

    def q(arr):
        return np.array(arr)
    
    if species_order is not None:
        ids = np.array(species_order)
    else:
        ids = np.unique([x.split('_')[0] for x in M.index])

    if len(ids)>2:
        d = M.values.copy()
        d[d<align_thr]=0
        x,y = d.nonzero()
        x,y = np.unique(np.sort(np.vstack((x,y)).T,axis=1),axis=0).T
        values = d[x,y]
        nodes = q(M.index)

        node_pairs = nodes[np.vstack((x,y)).T]
        sn1 = q([xi.split('_')[0] for xi in node_pairs[:,0]])
        sn2 = q([xi.split('_')[0] for xi in node_pairs[:,1]])
        filt = np.zeros(len(x), dtype=bool)
#        filt = np.logical_or(
#            np.logical_or(np.logical_and(sn1==ids[0],sn2==ids[1]),np.logical_and(sn1==ids[1],sn2==ids[0])),
#            np.logical_or(np.logical_and(sn1==ids[1],sn2==ids[2]),np.logical_and(sn1==ids[2],sn2==ids[1]))
#        )
        for i in range(len(ids) - 1):
            for j in range(i + 1, len(ids)):
                filt |= np.logical_or(np.logical_and(sn1 == ids[i], sn2 == ids[j]), np.logical_and(sn1 == ids[j], sn2 == ids[i]))

        x,y,values=x[filt],y[filt],values[filt]
        
        d=dict(zip(ids,list(np.arange(len(ids)))))        
        depth_map = dict(zip(nodes,[d[xi.split('_')[0]] for xi in nodes]))
        data =  nodes[np.vstack((x,y))].T
        for i in range(data.shape[0]):
            if d[data[i,0].split('_')[0]] > d[data[i,1].split('_')[0]]:
                data[i,:]=data[i,::-1]
        R = pd.DataFrame(data = data,columns=['source','target'])
        
        R['Value'] = values       
    else:
        d = M.values.copy()
        d[d<align_thr]=0
        x,y = d.nonzero()
        x,y = np.unique(np.sort(np.vstack((x,y)).T,axis=1),axis=0).T
        values = d[x,y]
        nodes = q(M.index)
        R = pd.DataFrame(data = nodes[np.vstack((x,y))].T,columns=['source','target'])
        R['Value'] = values
        depth_map=None
    
    try:
        from holoviews import dim
        #from bokeh.models import Label
        import holoviews as hv
        hv.extension('bokeh',logo=False)
        hv.output(size=100)        
    except:
        raise ImportError('Please install holoviews-samap with `!pip install holoviews-samap`.')

    def f(plot,element):
        plot.handles['plot'].sizing_mode='scale_width'    
        plot.handles['plot'].x_range.start = -600    
        #plot.handles['plot'].add_layout(Label(x=plot.handles['plot'].x_range.end*0.78, y=plot.handles['plot'].y_range.end*0.96, text=id2))
        plot.handles['plot'].x_range.end = 1500    
        #plot.handles['plot'].add_layout(Label(x=0, y=plot.handles['plot'].y_range.end*0.96, text=id1))


    sankey1 = hv.Sankey(R, kdims=["source", "target"])#, vdims=["Value"])

    cmap = params.get('cmap','Colorblind')
    label_position = params.get('label_position','outer')
    edge_line_width = params.get('edge_line_width',0)
    show_values = params.get('show_values',False)
    node_padding = params.get('node_padding',4)
    node_alpha = params.get('node_alpha',1.0)
    node_width = params.get('node_width',40)
    node_sort = params.get('node_sort',True)
    frame_height = params.get('frame_height',2000)
    frame_width = params.get('frame_width',1600)
    bgcolor = params.get('bgcolor','snow')
    apply_ranges = params.get('apply_ranges',True)


    sankey1.opts(cmap=cmap,label_position=label_position, edge_line_width=edge_line_width, show_values=show_values,
                 node_padding=node_padding,node_cmap=depth_map, node_alpha=node_alpha, node_width=node_width,
                 node_sort=node_sort,frame_height=frame_height,frame_width=frame_width,bgcolor=bgcolor,
                 apply_ranges=apply_ranges,hooks=[f])

    return sankey1

fig=sankey_plot(df, align_thr=0.3, species_order = ['HO','MA','MU','TU','DO','OV','GO','SU'])

import holoviews as hv
hv.extension('matplotlib')
hv.save(obj=fig, filename='HO_MA_MU_TU_DO_OV_GO_SU', fmt='svg')


