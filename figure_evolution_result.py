"""
Figure: convergence of LGA and communities of network

@auth: Yu-Hsiang Fu
@date: 2016/05/08
@ update 2018/03/26
"""
# --------------------------------------------------------------------------------
# 1.Import modular
# --------------------------------------------------------------------------------
# import modular
import matplotlib.cbook
import matplotlib.pyplot as plt
import networkx as nx
import warnings

# import custom-modular
import util.handler.pickle_handler as ph

# import folder-constant
from util.constant.constant_folder import FOLDER_FILE
from util.constant.constant_folder import FOLDER_IMAGE

# import graph-constant
from util.constant.constant_graph import NODE_LAYOUT_XY


# --------------------------------------------------------------------------------
# 2.Define variable
# --------------------------------------------------------------------------------
# GA variable
IDV_LENGTH = 0
IDV_FITNESS = 'fitness'
IDV_GENOTYPE = 'genotype'
IDV_PHENOTYPE = 'phenotype'

# figure variables
# COLOR_LIST = ['gray', 'orange', 'y', 'b', 'c', 'm', 'r', 'k']
# MARKER_LIST = ['^', 'v', '8', 'H', 's', 'D']
PLOT_NET_X_SIZE = 4
PLOT_NET_Y_SIZE = 4
PLOT_X_SIZE = 3
PLOT_Y_SIZE = 3
PLOT_DPI = 300
PLOT_FORMAT = 'png'


# --------------------------------------------------------------------------------
# 3.Define function
# --------------------------------------------------------------------------------
def draw_convergence_figure(net_name, evo_result):
    best_idv, fitness_avg, fitness_best = evo_result

    # create figure
    fig, ax = plt.subplots(figsize=(PLOT_X_SIZE, PLOT_Y_SIZE), facecolor='w')

    # draw plot
    ax.plot(fitness_best, color="r", linewidth=1, marker="o", markersize=4, markevery=1, fillstyle='none')
    ax.plot(fitness_avg, color="b", linewidth=1, marker="o", markersize=4, markevery=1, fillstyle='none')

    # plot setting
    ax.grid(color="gray", linestyle="dotted", linewidth=0.5)
    ax.set_xlim(-1, len(fitness_avg) + 1)
    ax.set_ylim(-0.01, max(fitness_best) + 0.02)
    ax.set_xlabel("Generation", fontdict={'fontsize': 8})
    ax.set_ylabel("Fitness: modularity Q", fontdict={'fontsize': 8})
    ax.tick_params(axis="both", direction="in", which="major", labelsize=6)

    # legend-text
    legend_text = ["fitness-best", "fitness-avg"]
    ax.legend(legend_text, loc=4, fontsize='medium', prop={'size': 6}, ncol=1, framealpha=0.5)

    # save image
    image_path = "{0}{1}, lga-convergence.png".format(FOLDER_IMAGE, net_name)
    plt.tight_layout()
    plt.savefig(image_path, dpi=PLOT_DPI, format=PLOT_FORMAT)
    plt.close()


def draw_network_figure(net_name, g, evo_result):
    best_idv = evo_result[0]

    # community info.
    community_index = 1
    community_list = dict()
    node_community = dict()

    for community in best_idv[IDV_PHENOTYPE]:
        # nodes' community
        for i in community:
            node_community[i] = community_index

        # community list
        community_list[community_index] = community
        community_index += 1

    # --------------------------------------------------
    # edge info.
    inside_edge = dict()
    between_edge = list()

    for (ei, ej) in g.edges():
        ci = node_community[ei]
        cj = node_community[ej]

        if ci == cj:
            # same community
            if ci in inside_edge:
                inside_edge[ci].append((ei, ej))
            else:
                inside_edge[ci] = [(ei, ej)]
        else:
            # different community
            between_edge.append((ei, ej))

    # --------------------------------------------------
    # address warning
    warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

    # create figure
    fig, ax = plt.subplots(figsize=(PLOT_NET_X_SIZE, PLOT_NET_Y_SIZE), facecolor='w')

    # draw plot
    color_max = max(community_list)
    pos = {i: g.node[i][NODE_LAYOUT_XY] for i in g}

    # between-community-edge
    nx.draw_networkx_edges(g, pos=pos, edgelist=between_edge, edge_color='gray', width=1.0, alpha=0.5)

    # inside-community-edge
    for i in inside_edge:
        edge_color_list = [i] * len(inside_edge[i])

        nx.draw_networkx_edges(g,
                               pos=pos,
                               edgelist=inside_edge[i],
                               edge_color=edge_color_list,
                               edge_cmap=plt.cm.jet,
                               width=1.5,
                               edge_vmin=0.9,
                               edge_vmax=color_max+0.1,
                               alpha=0.8)

    # draw community-nodes
    for i in community_list:
        node_color_list = [i] * len(community_list[i])

        nx.draw_networkx_nodes(g,
                               pos=pos,
                               nodelist=community_list[i],
                               node_size=25,
                               node_color=node_color_list,
                               cmap=plt.cm.jet,
                               vmin=0.9,
                               vmax=color_max+0.1)

    # plot setting
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.axis('off')

    # save image
    image_path = "{0}{1}, lga-identified-community.png".format(FOLDER_IMAGE, net_name)
    plt.tight_layout()
    plt.savefig(image_path, dpi=PLOT_DPI, format=PLOT_FORMAT)
    plt.close()


# --------------------------------------------------------------------------------
# 4.Main function
# --------------------------------------------------------------------------------
def main_function():
    # test network
    # filename_list = ["LFR_benchmark_n=300_u=0.05"]
    #
    filename_list = ["14p_gcc",
                     "karate_gcc",
                     "dolphins_gcc",
                     "LFR_benchmark_n=300_u=0.05"]

    # read evolutionary pickle files
    for net_name in filename_list:
        print(" - [Net] {0}:".format(net_name))
        print(" -- Read network pickle file")
        file_path = "{0}{1}, analysis.pickle".format(FOLDER_FILE, net_name)
        g = ph.read_pickle_file(file_path)

        print(" -- Read LGA pickle files")
        file_path = "{0}{1}-lga.pickle".format(FOLDER_FILE, net_name)
        evo_result = ph.read_pickle_file(file_path)

        print(" -- Draw LGA convergence figure")
        draw_convergence_figure(net_name, evo_result)

        print(" -- Draw LGA identified community figure")
        draw_network_figure(net_name, g, evo_result)

        print(" - [/Net]\n")


if __name__ == "__main__":
    main_function()
