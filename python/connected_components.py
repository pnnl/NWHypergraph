import numpy as np
import sys
import scipy
import scipy.io
import sys

MI = sys.maxsize   #  MAXINT aka -1 aka unvisited


# Q: how to do without scanning whole graph
def replace(arr, old, new):            # "compress" in bgl17 cc
    temp = (arr == old)
    return new * temp + arr * (~temp)


# edx -- hyperedge numbering
# ndx -- hypernode numbering
def comps(edx, ndx):   #  "struct of arrays"  edge list bipartite graph  
 
    n = max(edx) + 1     # number of hyperedges
    m = max(ndx) + 1     # number of hypernodes
    E = MI * np.ones(n)  # label for each hyperedge -- all initialized to maxint (-1)
    N = MI * np.ones(m)  # label for each hypernode -- all initialized to maxint (-1)

    # for each bi-partite edge

    for idx in range(len(edx)):
        node = ndx[idx]    # hypernode
        edge = edx[idx]    # edge hyperedge

        if E[edge] == MI:        # if hyperedge is unlabeled
            E[edge] = edge       # label it with itself
        if E[edge] == N[node]:   # if hyperedge label is equal to hypernode label
            continue             #   continue
        if N[node] == MI:        # if hypernode is unlabeled
            N[node] = E[edge]    #   give it the hyperedge label
        elif N[node] > E[edge]:             # if hypernode label is larger than hyperedge label
            temp = N[node]                  #   temp gets larger value
            N = replace(N, temp, E[edge])   #   replace all occurrences of higher value with lower value in E
            E = replace(E, temp, E[edge])   #   replace all occurrences of higher value with lower value in N
        elif N[node] < E[edge]:             # if hyperedge label is larger than hypernode label
            temp = E[edge]                  #   temp gets larger value
            E = replace(E, temp, N[node])   #   replace all occurrences of higher value with lower value in E
            N = replace(N, temp, N[node])   #   replace all occurrences of higher value with lower value in E


    ER = np.unique(E[np.argwhere(E != MI)])   # array of unique labels of connected components \ unlabeled hyperedges

    edgecomps = list()

    for udx in np.sort(ER):                                     # for every component, in increasing order
        edgecomps.append(np.argwhere(E == udx).transpose()[0])  #   gather the associated elements in E with that label

    components = list()

    for comp in edgecomps:                                              # for each array of hyperedges in given component
        f = np.frompyfunc(lambda x: x in set(comp), 1, 1)               #   gather the associated bipartite edges
        temp = np.argwhere(f(edx)).transpose()
        components.append([edx[temp].flatten(), ndx[temp].flatten()])

    return components
