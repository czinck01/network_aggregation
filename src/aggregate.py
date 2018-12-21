import sys, os, operator, time, math, datetime
import networkx as nx
import numpy as np

from funcassociate.client import _fc as fc
from sklearn.cluster import KMeans

DSDs = '../data/networks/DSDs'
NODELISTS = '../data/nodelists'
GENEFILE = '../data/ids/gene_ids.txt'
PARAMETERS = 'params'

NUM_GENES = 21115
NUM_GOS = 20270
K = 250
P = .05

# Takes a path and loads all npy files in that directory into np matrices
def load_networks(path):
    start = time.time()
    adjs = []
    for filename in sorted(os.listdir(path)):
        adj = np.load(os.path.join(path, filename))
        adjs.append(adj)
    end = time.time()
    print("Loaded networks in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return adjs

# Takes a path and reads all nodelists files in that directory into lists
def read_nodelists(path):
    nodelists = []
    for filename in sorted(os.listdir(path)):
        with open(os.path.join(path, filename), "r") as file:
            nodelist = []
            for line in file.readlines():
                nodelist.append(int(line.strip()))
            nodelists.append(nodelist)
    return nodelists

# Takes a filename that corresponds to a set of parameters and returns them
def read_parameters(filename):
    params = []
    with open(filename, "r") as file:
        for line in file.readlines():
            params.append(float(line))
    return params

# Takes a filename that corresponds to a list of key value pairs and creates a
# dictionary of gene ids 
def gene_id_dict(filename):
    gene_ids = {}
    with open(filename, "r") as file:
        for line in file.readlines():
            d = line.strip().split("\t")
            gene_ids[int(d[1])] = d[0]
    return gene_ids

# Takes an unnormalized weight vector and creates the normalized weight vector
def normalize_vector(raw):
    sum = np.sum(raw)
    normal = [x/float(sum) for x in raw]
    return normal

# Takes a list of np matrices and a list of nodes of the same length and resizes
# the matrices to the same size, where each gene is the same index in all
# matrices
def resize_networks(networks, nodelists):
    start = time.time()
    len_networks = len(networks)
    adjs = []
    for i in range(0, len_networks):
        len_nodelist = len(nodelists[i])
        adj = np.zeros((NUM_GENES, NUM_GENES))
        for j in range(0, len_nodelist):
            for k in range(0, len_nodelist):
                adj[nodelists[i][j]][nodelists[i][k]] = networks[i][j][k]
        adjs.append(adj)
    end = time.time()
    print("Resized networks in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return adjs

# Takes a list of np matrices and a weight vector and sums all matrices together
# scaled by the weights in the vector
def aggregate(adjs, wvec):
    start = time.time()
    agg = np.zeros_like(adjs[0])
    for i in range(0, len(adjs)):
        agg = np.add(agg, np.multiply(adjs[i], (wvec[i])))
    end = time.time()
    print("Aggregated networks in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return agg

# Takes a k value, an aggregated matrix, and a dictionary of gene ids. Uses
# k-means clustering on the aggregated matrix to produce clusters of genes
def cluster(k, agg, gene_ids):
    start = time.time()
    clustering = KMeans(n_clusters=k)
    labels = clustering.fit(agg).labels_
    clusters = d = [[] for x in xrange(k)]
    for i in range(0, labels.size):
           clusters[labels[i]].append(gene_ids[i])
    end = time.time()
    print("Clustered aggregated network in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return clusters

# Takes a p value and a set of clusters and scores each cluster using Func
# Associate. Returns the percentage of enriched clusters, the percentage of
# genes in enriched clusters, and the best GO term for each enriched cluster
def score(pval_cutoff, clusters):
    start = time.time()
    cluster_score = 0
    gene_score = 0
    cluster_terms = []
    func = fc()
    for i in range(0, len(clusters)):
        try:
            ret = func.functionate(query=clusters[i], species="Homo sapiens",
            namespace="hgnc_symbol", cutoff=(pval_cutoff))["over"]
            p = ret[0][4] * NUM_GOS
	    size = ret[0][1]
            term = ret[0][7]
        except:
            p = 1

        if (p < pval_cutoff):
            cluster_score += 1
	    gene_score += size
            cluster_terms.append(term)
	else:
	    cluster_terms.append("")
    end = time.time()
    print("Scored clusters in " +
        str(datetime.timedelta(seconds=int(end-start))))
    return cluster_score/float(K), gene_score/float(NUM_GENES), cluster_terms

# Takes a set of clusters, their GO terms, 2 scores, and the weight vector used
# to produce these results and writes this information to filename
def to_file(filename, clusters, cluster_terms, cs, gs, wvec):
    with open(filename, "w+") as file:
        file.write("{}% Cluster Enrichment\t{}% Gene Enrichment\n"
		.format(100 * cs, 100 * gs))
	file.write("Weights: {}\n".format(wvec))
        for i in range(0, len(clusters)):
            file.write("Cluster {}: {}\n".format(i+1, cluster_terms[i]))
            for j in range(0, len(clusters[i])):
                file.write("{}\t".format(clusters[i][j]))
            file.write("\n")

if __name__ == '__main__':
    # Load data
    networks = load_networks(DSDs)
    nodelists = read_nodelists(NODELISTS)
    gene_ids = gene_id_dict(GENEFILE)
    adjs = resize_networks(networks, nodelists)

    for filename in sorted(os.listdir(PARAMETERS)):

        # Parameters
        print("Using {} parameters".format(filename))
        params = read_parameters(os.path.join(PARAMETERS, filename))
        wvec = normalize_vector(params)
        
	# Aggregate, cluster, score, and output
        agg = aggregate(adjs, wvec)
        clusters = cluster(K, agg, gene_ids)
	cs, gs, cluster_terms = score(P, clusters)
        to_file(filename + ".out", clusters, cluster_terms, cs, gs, wvec)
