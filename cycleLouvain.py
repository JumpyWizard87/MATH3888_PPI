import networkx as nx
import community as community_louvain

def cycleLouvain(cycleG,n=1,label=False,resolution=1.0):
    """
    Example use: outputGraph=cycleLouvain(G,1,False)
    n iterations of the Louvian algorithm.
    cycleG=initial graph to feed in for Louvain iteration
    n=number of iterations, node after n=20 returns same community of ~7 proteins.
    label=output graph with node labels (True otherwise set to False e.g. if there are lots of nodes)
    note that this outputs labels even if False, might remove if can't fix.
    resolution=favours formation of smaller communities if set to >1, won't need many iterations if so.
    output=plot of network, returns graph object
    uses: get smaller network with sgs1 and use other algorithms to investigate graph.
    to add: we have a list of target nodes, could stop iteration if split into different communities. 
    
    """
    # Base case:
    if n == 0:
        nx.draw(cycleG,node_size=10,with_labels=label)
        return cycleG
    # Computing Louvain algorithm:
    partition=community_louvain.best_partition(cycleG,resolution=resolution)
    # Getting all proteins in SGS1's community:
    targetCommunity={key for key, value in partition.items() if value == partition[protein_target]}
    subG=cycleG.subgraph(targetCommunity) 
    return cycleLouvain(subG,(n-1),resolution)


