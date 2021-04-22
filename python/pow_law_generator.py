import hypernetx as hnx
import networkx as nx
import pandas as pd
import numpy as np
from collections import Counter, OrderedDict
import pickle
import matplotlib.pyplot as plt
import datetime as dtm


class Demo():
    """ 
    Simulation Graph generator

    """

    def __init__(self, use_nwhy=False, seed=0, collapsed=True, saved_state='demo_hypergraph.p'):
        
        default_data = {
            "expected_deg": [],
            "expected_size_hist": {}
        }
        
        try:
            self.hypergraph = hnx.Hypergraph.recover_from_state(saved_state)
            self.data = self.hypergraph.edges.data
            self.labels = self.hypergraph.edges.labels
        except:
            np.random.seed(seed)
            n, estubs = hist_to_stubs(default_data["expected_size_hist"])  # edges
            m, vstubs = list_to_stubs(default_data["expected_deg"])  # nodes
            self.data = create_data(estubs, vstubs)  # edges x nodes
            self.labels = create_labels(n, m)  # edges x nodes
            entity = hnx.StaticEntitySet(data=self.data, labels=self.labels)
            self.hypergraph = hnx.Hypergraph(entity, use_nwhy=use_nwhy, filepath=saved_state)
        if collapsed:
            try:
                self.collapsed = hnx.Hypergraph.recover_from_state(f'{saved_state[:-2]}_collapsed.p' )
            except:
                self.collapsed,self.edge_classes = self.hypergraph.collapse_edges(return_equivalence_classes=True)
                self.collapsed.set_state(equivalence_classes=self.edge_classes)
                self.collapsed.filepath = f'{saved_state[:-2]}_collapsed.p'           
                self.collapsed.save_state(fpath=self.collapsed.filepath)
        
# expected_deg, expected_size_hist


def list_to_stubs(degree_sequence):
    m = len(degree_sequence)
    stubs = list()
    for vdx in range(m):
        stubs += [vdx] * degree_sequence[vdx]
    return m, stubs


def hist_to_stubs(degree_histogram):
    stubs = list()
    n = 0
    for k, v in degree_histogram.items():
        for vdx in range(n, n + v):
            stubs += [vdx] * k
        n += v
    return n, stubs


def create_data(estubs, vstubs):
    assert len(estubs) == len(vstubs)
    np.random.shuffle(estubs)
    np.random.shuffle(vstubs)
    return np.array(list(zip(estubs, vstubs)))


def create_labels(num_edges, num_nodes, edgeprefix='e', nodeprefix='v', edgelabel='Edges', nodelabel='Nodes'):
    enames = np.array([f'{edgeprefix}{idx}' for idx in range(num_edges)])
    nnames = np.array([f'{nodeprefix}{jdx}' for jdx in range(num_nodes)])
    return OrderedDict([(edgelabel, enames), (nodelabel, nnames)])

def add_nwhy(h, fpath=None):
    """
    Add nwhy functionality to a hypergraph.

    Parameters
    ----------
    h : hnx.Hypergraph

    Returns
    -------
    hnx.Hypergraph

    """
    if h.nwhy:
        h.filepath = fpath
        return h
    elif not h.isstatic:
        return hnx.Hypergraph(hnx.StaticEntitySet(h.edges), use_nwhy=True, filepath=fpath)
    else:
        data = h.edges.data
        labels = h.edges.labels
        sd = h.state_dict
        sd['sedgelg'] = dict()
        sd['snodelg'] = dict()
        H = hnx.Hypergraph(hnx.StaticEntitySet(data=data, labels=labels), use_nwhy=True, filepath=fpath)
        H.state_dict.update(sd)
        return H
    
    
# Hypergraph statistics


