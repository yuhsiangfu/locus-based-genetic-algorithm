"""
Locus-based Genetic Algorithm (LGA)

@auth: Yu-Hsiang Fu
@date: 2016/04/28
@update: 2018/03/26
"""
# --------------------------------------------------------------------------------
# 1.Import modular
# --------------------------------------------------------------------------------
# import modular
import copy as c
import networkx as nx
import numpy as np
import pickle
import random as r

# import custom-modular
import util.data_structure.disjoint_set as djs
import util.handler.edgelist_handler as eh
import util.handler.pickle_handler as ph
import util.handler.pairvalue_handler as pvh

# import folder-constant
from util.constant.constant_folder import FOLDER_EDGELIST
from util.constant.constant_folder import FOLDER_FILE

# import graph-constant
from util.constant.constant_graph import NODE_COMMUNITY
from util.constant.constant_graph import NODE_DEGREE
from util.constant.constant_graph import NODE_LAYOUT_XY


# --------------------------------------------------------------------------------
# 2.Define variable
# --------------------------------------------------------------------------------
# program variable
GENE_TO_NODE = {}
NODE_TO_GENE = {}
NEIGHBOR_LIST = {}

# GA variable
IDV_LENGTH = 0
IDV_FITNESS = 'fitness'
IDV_GENOTYPE = 'genotype'
IDV_PHENOTYPE = 'phenotype'


# --------------------------------------------------------------------------------
# 3.Define function
# --------------------------------------------------------------------------------
def create_network(edge_list):
    g = nx.parse_edgelist(edge_list, nodetype=int)
    g = g.to_undirected()
    g.remove_edges_from(g.selfloop_edges())

    return g


def add_new_attribute(g, new_attr, pair_dict):
    for key in pair_dict.keys():
        g.node[key][new_attr] = pair_dict[key]

    return g


def p_deepcopy(data, dumps=pickle.dumps, loads=pickle.loads):
    return loads(dumps(data, -1))


def modularity(g, community_list):
    # ls, ds variables
    intra_degree = {i: 0 for i in range(0, len(community_list))}  # ds
    intra_edges = {i: 0 for i in range(0, len(community_list))}   # ls

    # calculate ds, time complexity: O(V)
    community_index = 0
    community_id = {}

    for community in community_list:
        tmp_index = c.copy(community_index)
        for i in community:
            intra_degree[tmp_index] += g.node[i][NODE_DEGREE]
            community_id[i] = tmp_index
        community_index += 1

    # calculate ls, time complexity: O(E)
    for (ei, ej) in g.edges():
        if community_id[ei] == community_id[ej]:
            intra_edges[community_id[ei]] += 1
        else:
            pass

    # calculate modularity Q, time complexity: O(C)
    q = 0
    num_edges = g.number_of_edges()

    for i in range(0, len(community_list)):
        ls = intra_edges[i] / num_edges
        ds = pow((intra_degree[i] / (2 * num_edges)), 2)
        q += (ls - ds)

    return q


# --------------------------------------------------
# generate function
# --------------------------------------------------
def generate_genotype():
    # genotype = []
    #
    # for i in range(0, IDV_LENGTH):
    #     node_id = GENE_TO_NODE[i]
    #     neighbor_id = r.choice(NEIGHBOR_LIST[node_id])
    #     neighbor_gene = NODE_TO_GENE[neighbor_id]
    #     genotype.append(c.copy(neighbor_gene))
    #
    genotype = [NODE_TO_GENE[r.choice(NEIGHBOR_LIST[GENE_TO_NODE[i]])] for i in range(0, IDV_LENGTH)]

    return genotype


def generate_phenotype(idv):
    # create disjoint-set
    ds = djs.disjoint_set(IDV_LENGTH)

    for x in range(0, IDV_LENGTH):
        y = idv[IDV_GENOTYPE][x]
        ds.union(x, y)

    # create community list, map node_index to node_id
    community_list = {}

    for i in range(IDV_LENGTH):
        node_id = GENE_TO_NODE[i]
        ri = ds.find_set(i)

        if ri in community_list:
            community_list[ri].append(node_id)
        else:
            community_list[ri] = [node_id]

    return list(community_list.values())


def generate_an_individual():
    individual = dict()
    individual[IDV_GENOTYPE] = generate_genotype()

    return individual


def generate_two_individuals(idv_pool, rate_crossover=0.5):
    # random select parent
    # parent_x = dict()
    # parent_y = dict()
    pool_size = len(idv_pool)

    if pool_size == 1:
        parent_x = p_deepcopy(idv_pool[0][IDV_GENOTYPE])
        parent_y = generate_an_individual()
    else:
        parent_index = r.sample(range(0, pool_size), 2)
        parent_x = p_deepcopy(idv_pool[parent_index[0]][IDV_GENOTYPE])
        parent_y = p_deepcopy(idv_pool[parent_index[1]][IDV_GENOTYPE])

    # uniform crossover
    new_idv1 = dict()
    new_idv2 = dict()
    new_idv1[IDV_GENOTYPE] = []
    new_idv2[IDV_GENOTYPE] = []

    # assign gene value
    if r.random() <= rate_crossover:
        for i in range(0, IDV_LENGTH):
            # mask == 0, px[i] -> idv1[i], py[i] -> idv2[i]
            if r.randint(0, 1) == 0:
                new_idv1[IDV_GENOTYPE].append(parent_x[i])
                new_idv2[IDV_GENOTYPE].append(parent_y[i])
            # mask == 1, px[i] -> idv2[i], py[i] -> idv1[i]
            else:
                new_idv1[IDV_GENOTYPE].append(parent_y[i])
                new_idv2[IDV_GENOTYPE].append(parent_x[i])
    else:
        new_idv1[IDV_GENOTYPE] = parent_x
        new_idv2[IDV_GENOTYPE] = parent_y

    return new_idv1, new_idv2


def generate_population(size_population=100):
    idv_pool = []
    for i in range(0, size_population):
        idv_pool.append(generate_an_individual())

    return idv_pool


# --------------------------------------------------
# GA operator
# --------------------------------------------------
def selection(idv_pool, size_population=100, rate_selection=0.1):
    # truncate selection
    cut_index = int(rate_selection * size_population)
    idv_pool = idv_pool[0: cut_index]

    return idv_pool


def crossover(idv_pool, size_population=100, rate_crossover=0.5):
    pool_size = len(idv_pool)
    empty_size = size_population - len(idv_pool)
    new_pool = p_deepcopy(idv_pool)

    if (empty_size % 2) == 0:
        for i in range(pool_size, size_population, 2):
            new_idv1, new_idv2 = generate_two_individuals(idv_pool, rate_crossover)
            new_pool.append(new_idv1)
            new_pool.append(new_idv2)
    else:
        for i in range(pool_size, size_population - 1, 2):
            new_idv1, new_idv2 = generate_two_individuals(idv_pool, rate_crossover)
            new_pool.append(new_idv1)
            new_pool.append(new_idv2)
        new_pool.append(generate_two_individuals(idv_pool, rate_crossover)[r.randint(0, 1)])

    return new_pool


def mutation(idv_pool, size_population=100, rate_mutation=0.05):
    # bit-by-bit mutation
    for i in range(0, size_population):
        idv = idv_pool[i][IDV_GENOTYPE]

        for j in range(0, IDV_LENGTH):
            if r.random() <= rate_mutation:
                node_id = GENE_TO_NODE[j]
                neighbor_id = r.choice(NEIGHBOR_LIST[node_id])
                idv[j] = c.copy(NODE_TO_GENE[neighbor_id])

    return idv_pool


# --------------------------------------------------
def initialization(g):
    global IDV_LENGTH, GENE_TO_NODE, NODE_TO_GENE, NEIGHBOR_LIST

    # length of individual
    IDV_LENGTH = c.copy(g.number_of_nodes())

    # mapping gene-to-node and node-to-gene
    gene_index = 0

    for i in g:
        GENE_TO_NODE[gene_index] = c.copy(i)
        NODE_TO_GENE[i] = c.copy(gene_index)
        gene_index += 1

    # create neighbor list
    for i in g:
        NEIGHBOR_LIST[i] = list(g.neighbors(i))
        g.node[i][NODE_DEGREE] = len(NEIGHBOR_LIST[i])

    return g


def evaluation(g, idv_pool, size_population=100):
    for i in range(0, size_population):
        idv = idv_pool[i]
        idv[IDV_PHENOTYPE] = generate_phenotype(idv)
        idv[IDV_FITNESS] = modularity(g, idv[IDV_PHENOTYPE])

    idv_pool = sorted(idv_pool, key=lambda x: x[IDV_FITNESS], reverse=True)

    return idv_pool


def locus_based_genetic_algorithm(g,
                                  num_evolution=1,
                                  num_generation=100,
                                  size_population=100,
                                  rate_selection=0.1,
                                  rate_crossover=0.5,
                                  rate_mutation=0.05):
    # initialization
    evo_best = []
    g = initialization(g)

    # evolution
    for i in range(0, num_evolution):
        print(" --- Evolution {0}".format(i + 1))
        fitness_avg = []
        fitness_best = []

        # create individual pool
        idv_pool = generate_population(size_population)

        # generation
        idv_best = dict()

        for j in range(0, num_generation):
            # show progress
            # if (j + 1) % 50 != 0:
            #     print(".", end="")
            # elif (j + 1) == num_generation:
            #     print(".")
            # else:
            #     print(".")

            # evaluate individuals
            idv_pool = evaluation(g, idv_pool, size_population)

            # maintain idv_best
            if not idv_best:
                idv_best = p_deepcopy(idv_pool[0])
            elif idv_pool[0][IDV_FITNESS] > idv_best[IDV_FITNESS]:
                idv_best = p_deepcopy(idv_pool[0])
            else:
                pass

            # record fitness_avg and fitness_best
            fitness_avg.append(np.mean([idv[IDV_FITNESS] for idv in idv_pool]))
            fitness_best.append(idv_best[IDV_FITNESS])

            # stop condition
            if (j + 1) == num_generation:
                break
            else:
                pass

            # genetic operation: selection, crossover and mutation
            idv_pool = selection(idv_pool, size_population, rate_selection)
            idv_pool = crossover(idv_pool, size_population, rate_crossover)
            idv_pool = mutation(idv_pool, size_population, rate_mutation)

        # maintain evo_best
        if not evo_best:
            evo_best = [idv_best, fitness_avg, fitness_best]
        elif idv_best[IDV_FITNESS] > evo_best[0][IDV_FITNESS]:
            evo_best = [idv_best, fitness_avg, fitness_best]
        else:
            pass

    return evo_best


# --------------------------------------------------------------------------------
# 4.Main function
# --------------------------------------------------------------------------------
def main_function():
    # test network
    # filename_list = ["14p_gcc"]
    #
    filename_list = ["14p_gcc",
                     "karate_gcc",
                     "dolphins_gcc",
                     "LFR_benchmark_n=300_u=0.05"]

    # --------------------------------------------------
    # GA variables
    num_evolution = 10
    num_generation = 50
    size_population = 100
    rate_selection = 0.1
    rate_crossover = 0.8
    rate_mutation = 0.05

    # --------------------------------------------------
    # read edge-list file of networks
    print(" Locus-based genetic algorithm (LGA)")
    for net_name in filename_list:
        print(" - [Net] {0}:".format(net_name))
        print(" -- Read edge-list file")
        file_path = "{0}{1}.txt".format(FOLDER_EDGELIST, net_name)
        g = create_network(eh.read_edgelist(file_path))

        print(" -- Add new attribute")
        g = add_new_attribute(g,
                              NODE_LAYOUT_XY,
                              pvh.read_pairvalue("{0}{1}_pos.txt".format(FOLDER_EDGELIST, net_name)))
        g = add_new_attribute(g,
                              NODE_COMMUNITY,
                              pvh.read_pairvalue("{0}{1}_community.txt".format(FOLDER_EDGELIST, net_name), is_int=True))

        file_path = "{0}{1}, analysis.pickle".format(FOLDER_FILE, net_name)
        ph.write_pickle_file(g, file_path)

        print(" -- LGA evolution")
        evo_result = locus_based_genetic_algorithm(g,
                                                   num_evolution,
                                                   num_generation,
                                                   size_population,
                                                   rate_selection,
                                                   rate_crossover,
                                                   rate_mutation)

        print(" -- Save evolution result")
        file_path = "{0}{1}-lga.pickle".format(FOLDER_FILE, net_name)
        ph.write_pickle_file(evo_result, file_path)
        print(" - [/Net]\n")


if __name__ == "__main__":
    main_function()
