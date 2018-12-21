import sys
import numpy as np
import networkx as nx
import os

PATH = '../data/networks/anonymized'

if __name__ == '__main__':
    for filename in sorted(os.listdir(PATH)):
        print(filename)

	G = nx.read_weighted_edgelist(os.path.join(PATH,filename))

	print("# of Edges:\t{}".format(G.number_of_edges()))	
	print("# of Nodes:\t{}".format(G.number_of_nodes()))	

	print("Assort:\t{}".format(nx.degree_assortativity_coefficient(G)))
	print("Density:\t{}".format(nx.density(G)))
	
	degs = [x[1] for x in G.degree()]
        
        print("Max degree:\t{}".format(np.max(degs)))
	print("Mean degree:\t{}".format(np.mean(degs)))
        print("Median degree:\t{}".format(np.median(degs)))

        cents = nx.eigenvector_centrality(G).values()

	print("Max centrality:\t{}".format(np.max(cents)))
	print("Mean centrality:\t{}".format(np.mean(cents)))
	print("Median centrality:\t{}".format(np.median(cents)))
