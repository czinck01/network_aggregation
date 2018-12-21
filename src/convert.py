import sys
import networkx as nx

def main(argv):
    G = nx.read_weighted_edgelist(argv[1], create_using=nx.MultiGraph())
    with open("output.txt", 'w') as f:
	for e in G.edges():
	    weights = G[e[0]][e[1]]
	    num = 0.0
	    sum = 0.0
	    for k, v in weights.items():
		num += 1
		sum += float(v['weight'])
	    f.write("{}\t{}\t{}\n".format(e[0], e[1], sum/num))
if __name__ == '__main__':
        main(sys.argv)
