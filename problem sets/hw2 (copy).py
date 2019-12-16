import networkx as nx
import copy
import random
import numpy as np
import matplotlib.pyplot as plt
'''
# 5a
G = nx.read_adjlist("adj.adjlist")
print(G.nodes())
print(len(G.edges()))

centrality = nx.eigenvector_centrality(G)
print(sorted(((v, '{:0.2f}'.format(c)) for v, c in centrality.items()), key=lambda x: x[1],reverse=True ))

betweenness = nx.betweenness_centrality(G)
print(sorted(((v, '{:0.2f}'.format(c)) for v, c in betweenness.items()), key=lambda x: x[1],reverse=True ))

harmonic = nx.harmonic_centrality(G)
print(sorted(((v, '{:0.2f}'.format(c)) for v, c in harmonic.items()), key=lambda x: x[1],reverse=True ))

# 5b
degreeSeq = []
for i in G.nodes():
    for j in range(G.degree[i]):
        degreeSeq.append(i)
harmonicDist = [[] for i in range(len(G.nodes))]

for i in range(1000):
    deg = copy.deepcopy(degreeSeq)
    shuffle = True
    while shuffle:
        random.shuffle(deg)
        shuffle = False  
        tuples = zip(*[deg[i::2] for i in range(2)])
        for t in tuples:
            if tuples.count(t) + tuples.count((t[1],t[0])) > 1:
                shuffle = True
                break
            if t[0] == t[1]:
                shuffle = True
                break
    tuples = zip(*[deg[i::2] for i in range(2)])
    Gp=nx.Graph()
    Gp.add_edges_from(tuples)
    Gp.add_nodes_from(G.nodes)
    harmon = nx.harmonic_centrality(Gp)
    for j in Gp.nodes():
        harmonicDist[int(j)].append(harmon[j])
        
quart1 = []
quart2 = []
quart3 = []
x = range(len(G.nodes))
for i in x:
    h = harmonic[str(i)]
    quart1.append(h - np.quantile(harmonicDist[i], .25))
    quart2.append(h - np.quantile(harmonicDist[i], .50))
    quart3.append(h - np.quantile(harmonicDist[i], .75))

plt.figure()
plt.plot(x, [0 for i in range(len(x))], '--k')
plt.plot(x, quart2, 'ro-')
plt.fill_between(x, quart1, quart3, alpha=0.5)
plt.xlabel("Node Label")
plt.ylabel("Difference")
plt.title("Difference of Harmonic Centrality")
plt.savefig("harmonic_dist.png")
'''

# 6a
edgeList = [(1,0), (1,2), (1,3), (1,4), (2,3), (3,4), (5,6), (4,5), (0,2)] #, (0,6)]
G = nx.Graph()
G.add_edges_from(edgeList)
print(G.nodes())
print(len(G.edges()))

degree = nx.degree_centrality(G)
print(sorted(((v, '{:0.3f}'.format(c)) for v, c in degree.items()), key=lambda x: x[1],reverse=True ))

centrality = nx.eigenvector_centrality(G)
#print(sorted(((v, '{:0.3f}'.format(c)) for v, c in centrality.items()), key=lambda x: x[1],reverse=True ))

betweenness = nx.betweenness_centrality(G)
print(sorted(((v, '{:0.3f}'.format(c)) for v, c in betweenness.items()), key=lambda x: x[1],reverse=True ))

harmonic = nx.harmonic_centrality(G)
#print(sorted(((v, '{:0.3f}'.format(c)) for v, c in harmonic.items()), key=lambda x: x[1],reverse=True ))

plt.figure()
#plt.axis('off')
nx.draw_networkx(G)
plt.draw()
plt.annotate("Max Degree", (0.15, -0.2))
plt.annotate("Max Betweenness", (0,0.1))
#plt.axis('off')
plt.savefig("net.png")   
    
