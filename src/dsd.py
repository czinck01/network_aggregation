import sys
import numpy as np
import networkx as nx
import os

PATH = '../data/networks/anonymized'

def build_transition_matrix(adjacency_graph):
    n = np.size(adjacency_graph[0])

    transition = np.zeros((n,n))

    for i in range(n):
        s = np.sum(adjacency_graph[i])
        if s:
            transition[i] = np.divide(adjacency_graph[i], s)
        else:
            transition[i][i] = 1
    return transition

def calc_hescotts(transition, iters, v=True, n=None):
    if n is None:
        n = np.size(transition[0])
    m = np.size(transition[0])
    hescotts = np.zeros((n,m))
    for i in range(n): hescotts[i,i] = 1
    for i in range(iters):
        if v: print "Calculating hescotts for step "+str(i+1)+"..."
        hescotts = np.dot(hescotts, transition)
        for i in range(n): hescotts[i,i]+=1
    return hescotts

def calc_dsd(hescotts):
    n = np.size(hescotts[:,0])
    dsd = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            d = np.linalg.norm((hescotts[i,:]-hescotts[j,:]), ord=1)
            dsd[i][j] = d
            dsd[j][i] = d
    return dsd

# All code from this file outside this main function comes from Tufts DSD source
# code
def main(argv):
    for filename in sorted(os.listdir(PATH)):
        print (filename)
        G = nx.read_weighted_edgelist(os.path.join(PATH, filename))
        nlist = map(str, sorted(map(int, list(G.nodes()))))
        adj = nx.to_numpy_matrix(G, nodelist=nlist)
        trans = build_transition_matrix(adj)
        hescotts = calc_hescotts(trans, 7)
        dsd = calc_dsd(hescotts)
        np.save(filename, dsd)
        with open(filename + ".nlist", "w+") as f:
            for i in range(0, len(nlist)):
                f.write("{}\n".format(nlist[i]))

if __name__ == '__main__':
    main(sys.argv)
