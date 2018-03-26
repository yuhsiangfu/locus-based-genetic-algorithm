"""
Measure: modularity (set)
@auth:  Yu-Hsiang Fu
@date   2015/10/09
@update 2016/02/13
"""


# 模塊性: Newman's modularity
def modularity(G, community_list):
    """
    The estimated time complexity of this version (2016/02/13) is approximating
    O(V) + O(E)
    """
    import copy as c
    NODE_DEGREE = 'node_degree'

    # ls, ds variables
    intra_degree = {i: 0 for i in range(0, len(community_list))}  # ds
    intra_edges  = {i: 0 for i in range(0, len(community_list))}  # ls

    # calculate ds, time complexity: O(V)
    community_index = 0
    community_id = {}

    for com in community_list:
        tmp_index = c.copy(community_index)
        for i in com:
            intra_degree[tmp_index] += G.node[i][NODE_DEGREE]
            community_id[i] = tmp_index
        community_index += 1

    # calculate ls, time complexity: O(E)
    for (ei, ej) in G.edges():
        if community_id[ei] == community_id[ej]:
            intra_edges[community_id[ei]] += 1
        else:
            pass

    # calculate modularity Q, time complexity: O(C)
    modularity = 0
    num_edges = G.number_of_edges()
    for i in range(0, len(community_list)):
        ls = intra_edges[i] / num_edges
        ds = pow((intra_degree[i] / (2 * num_edges)), 2)
        modularity += (ls - ds)
    return modularity
