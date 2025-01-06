import sys
from matplotlib import pyplot as plt
from pathlib import Path
import networkx as nx
import numpy
import random
import scipy as sp
import csv


def get_two_node_hash_table():
    """
    Get all 2-node graphlet binary edge vectors
    see diagram for reference
    """
    # 2-node graphlet binary edge vectors
    # see diagram for reference
    a_1 = (1, 0, 0)
    a_2 = (1, 0, 0)
    a_hash = hash(a_1 + a_2)

    b_1 = (1, 1, 0)
    b_2 = (1, 0, 1)
    b_hash = hash(b_1 + b_2)

    c_1 = (1, 1, 1)
    c_2 = (1, 1, 1)
    c_hash = hash(c_1 + c_2)

    d_1 = (0, 1, 0)
    d_2 = (0, 0, 1)
    d_hash = hash(d_1 + d_2)

    e_1 = (0, 1, 1)
    e_2 = (0, 1, 1)
    e_hash = hash(e_1 + e_2)

    two_node_hash_table = {a_hash: 0, b_hash: 0, c_hash: 0, d_hash: 0, e_hash: 0}
    return two_node_hash_table


def generate_random_multi_graph(n, edge_probability=0.3, edge_label_probability=0.3):
    """
    Given the number of nodes ,n, and some edge and label probabilites,
    generate a random MultiDiGraph
    """

    print("Generating random graph")
    G = nx.MultiDiGraph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(random.randrange(0, 50)):
            if i != j:
                if random.random() < edge_probability:
                    if random.random() < edge_label_probability:
                        G.add_edge(i, j, label="ppi")
                    else:
                        G.add_edge(i, j, label="reg")
                        if random.random() < 0.7:
                            G.add_edge(j, i, label="reg")
                # want to generate a graph that is connected
                elif G.degree(j) == 0:
                    if random.random() < edge_label_probability:
                        G.add_edge(j, i, label="ppi")
                    else:
                        G.add_edge(j, i, label="reg")
    return G


def update_protein_id_dict(file_path, res_dict, start_index):
    """Helper function to update the protein ID dictionary from a file."""
    with open(file_path, "r") as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            for protein_id in row[:2]:
                res_dict.setdefault(protein_id, start_index)
                if res_dict[protein_id] == start_index:
                    start_index += 1
    return start_index


def get_protein_id_dict(ppi_path, reg_path):
    """Creates a dictionary mapping protein IDs to unique integers."""
    res_dict = {}
    start_index = 0
    start_index = update_protein_id_dict(ppi_path, res_dict, start_index)
    update_protein_id_dict(reg_path, res_dict, start_index)
    return res_dict


def process_edges(
    file_path, G, protein_id_dict, visited_nodes, label, node_limit, edge_limit
):
    """Helper function to process edges and add them to the graph."""
    print(f"currently processing {label} edges")
    with open(file_path, "r") as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        node_count = len(visited_nodes)
        edge_count = G.number_of_edges()

        for row in csv_reader:
            if node_count > node_limit or edge_count > edge_limit:
                break

            id1 = protein_id_dict[row[0]]
            id2 = protein_id_dict[row[1]]

            if id1 not in visited_nodes:
                visited_nodes.add(id1)
                node_count += 1

            if id2 not in visited_nodes:
                visited_nodes.add(id2)
                node_count += 1

            G.add_edge(id1, id2, label=label)
            edge_count += 1


def read_csv(
    ppi_path,
    reg_path,
    protein_id_dict,
    node_size_limit=float("inf"),
    edge_size_limit=float("inf"),
):
    """Reads CSV files and constructs a graph with edges labeled as 'ppi' or 'reg'."""
    G = nx.MultiDiGraph()
    visited_nodes = set()

    process_edges(
        ppi_path,
        G,
        protein_id_dict,
        visited_nodes,
        "ppi",
        node_size_limit,
        edge_size_limit,
    )
    process_edges(
        reg_path,
        G,
        protein_id_dict,
        visited_nodes,
        "reg",
        node_size_limit,
        edge_size_limit,
    )

    return G


def draw_labeled_multigraph(G, attr_name, ax=None):
    """
    Length of connectionstyle must be at least that of a maximum number of edges
    between pair of nodes. This number is maximum one-sided connections
    for directed graph and maximum total connections for undirected graph.
    """

    pos = nx.shell_layout(G)
    nx.draw_networkx_nodes(G, pos, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=20, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color="grey", ax=ax)

    labels = {
        tuple(edge): f"{attr_name}={attrs[attr_name]}"
        for *edge, attrs in G.edges(keys=True, data=True)
    }
    nx.draw_networkx_edge_labels(
        G,
        pos,
        labels,
        label_pos=0.3,
        font_color="blue",
        bbox={"alpha": 0},
        ax=ax,
    )


def get_two_node_graphlet_stats(G, two_node_hash_table):
    G_adj_list = get_two_node_adjacency_list(G)
    orbits_dict = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
    for neighbors in G_adj_list:
        i = neighbors[0]
        # print((i / len(G.nodes())) * 100, end="\r")
        for j in neighbors[1]:
            vectors = G_adj_list[i][1][j] + G_adj_list[j][1][i]
            orbits_dict = get_two_node_orbit_position(
                i, j, G_adj_list[i][1][j], G_adj_list[j][1][i], orbits_dict
            )
            if hash(tuple(vectors)) in two_node_hash_table:
                two_node_hash_table[hash(tuple(vectors))] += 1
    return two_node_hash_table, orbits_dict


def get_two_node_adjacency_list(G):
    """Get the adjacency list for a MultiDiGraph"""
    print("getting adjacency list")

    adj_list_vector = [{} for _ in range(len(G.nodes()))]

    for i, j, data in G.edges(data=True):
        label = data.get("label")
        if label == "ppi":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [1, 0, 0]
            # else:
            #     adj_list_vector[i][j][0] += 1
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [1, 0, 0]
            # else:
            #     adj_list_vector[j][i][0] += 1
        elif label == "reg":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [0, 1, 0]
            else:
                adj_list_vector[i][j][1] += 1
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [0, 0, 1]
            else:
                adj_list_vector[j][i][2] += 1

    final_adj_list_vector = [
        (i, neighbors) for i, neighbors in enumerate(adj_list_vector)
    ]
    return final_adj_list_vector


def get_two_node_orbit_position(a, b, a_vec, b_vec, orbit_dict):
    # check orbit positioning via vectors
    if a_vec == [1, 0, 0] and b_vec == [1, 0, 0]:
        if a not in orbit_dict[1]:
            orbit_dict[1] += [a]
        if b not in orbit_dict[1]:
            orbit_dict[1] += [b]
    elif (
        a_vec == [0, 1, 0]
        and b_vec == [0, 0, 1]
        or a_vec == [0, 0, 1]
        and b_vec == [0, 1, 0]
    ):
        if a_vec == [0, 1, 0] and a not in orbit_dict[2]:
            orbit_dict[2] += [a]
            orbit_dict[3] += [b]
        elif b not in orbit_dict[2]:
            orbit_dict[2] += [b]
            orbit_dict[3] += [a]

    elif (
        a_vec == [1, 1, 0]
        and b_vec == [1, 0, 1]
        or a_vec == [1, 0, 1]
        and b_vec == [1, 1, 0]
    ):
        if a_vec == [1, 1, 0] and a not in orbit_dict[4]:
            orbit_dict[4] += [a]
            orbit_dict[5] += [b]
        elif b not in orbit_dict[4]:
            orbit_dict[4] += [b]
            orbit_dict[5] += [a]

    elif a_vec == [1, 1, 1] and b_vec == [1, 1, 1]:
        if a not in orbit_dict[6]:
            orbit_dict[6] += [a]
        if b not in orbit_dict[6]:
            orbit_dict[6] += [b]

    elif a_vec == [0, 1, 1] and b_vec == [0, 1, 1]:
        if a not in orbit_dict[7]:
            orbit_dict[7] += [a]
        if b not in orbit_dict[6]:
            orbit_dict[7] += [b]

    return orbit_dict


def get_two_node_graphlet_dist_adj_matrix(G, two_node_hash_table):
    G_adj_matrix = get_adjacency_matrix(G)
    for i in range(len(G_adj_matrix)):
        for j in range(len(G_adj_matrix[0])):
            vector = G_adj_matrix[i][j] + G_adj_matrix[j][i]
            if hash(tuple(vector)) in two_node_hash_table:
                two_node_hash_table[hash(tuple(vector))] += 1
    return two_node_hash_table


def get_adjacency_matrix(G):
    """Get the adjacency matrix for a MultiDiGraph"""

    G_adj_matrix = nx.adjacency_matrix(G)
    adj_matrix = [
        [[0, 0, 0] for _ in range(len(G.nodes()))] for _ in range(len(G.nodes()))
    ]

    for edge in G.edges(data=True):
        i = edge[0]
        j = edge[1]
        edge_type = edge[2]["label"]

        if edge_type == "ppi":
            adj_matrix[i][j][0] = 1
            adj_matrix[j][i][0] = 1
        elif edge_type == "reg":
            adj_matrix[i][j][1] = 1
            adj_matrix[j][i][2] = 1

    return G_adj_matrix


def get_three_node_graphlet_dist_adj_list(G: nx.MultiDiGraph):
    three_node_hash = {}

    # create all the binary edge vectors
    adj_list_vector = [{} for _ in range(len(G.nodes()))]

    for i, j, data in G.edges(data=True):
        label = data.get("label")
        if label == "ppi":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [1, 0, 0]
            # else:
            #     adj_list_vector[i][j][0] += 1
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [1, 0, 0]
            # else:
            #     adj_list_vector[j][i][0] += 1
        elif label == "reg":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [0, 1, 0]
            else:
                adj_list_vector[i][j][1] += 1
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [0, 0, 1]
            else:
                adj_list_vector[j][i][2] += 1

    # find all combinations of potential 3 node graphlets
    # pick an edge between A and B
    # for each edge pair, find the union of neighbors between A and B
    three_node_combination = []
    for (
        i,
        j,
    ) in G.edges():
        neighbors = G.edges([i, j])
        for k in neighbors:
            if k[1] != i and k[1] != j:
                three_node_combination.append([i, j, k[1]])

                # for each triplet combination, we have to get all their binary edge vectors
                # for nodes A, B, C, there are six binary edge vectors:
                # A-B
                # A-C
                # B-A
                # B-C
                # C-A
                # C-B
                a1 = a2 = b1 = b2 = c1 = c2 = 0
                if j in adj_list_vector[i]:
                    a1 = adj_list_vector[i][j]
                else:
                    a1 = [0, 0, 0]
                if k[1] in adj_list_vector[i]:
                    a2 = adj_list_vector[i][k[1]]
                else:
                    a2 = [0, 0, 0]
                if i in adj_list_vector[j]:
                    b1 = adj_list_vector[j][i]
                else:
                    b1 = [0, 0, 0]
                if k[1] in adj_list_vector[j]:
                    b2 = adj_list_vector[j][k[1]]
                else:
                    b2 = [0, 0, 0]
                if i in adj_list_vector[k[1]]:
                    c1 = adj_list_vector[k[1]][i]
                else:
                    c1 = [0, 0, 0]
                if j in adj_list_vector[k[1]]:
                    c2 = adj_list_vector[k[1]][j]
                else:
                    c2 = [0, 0, 0]

                vector = a1 + a2 + b1 + b2 + c1 + c2

                if hash(tuple(vector)) not in three_node_hash:
                    three_node_hash[hash(tuple(vector))] = 0
                else:
                    three_node_hash[hash(tuple(vector))] += 1
        print((i / len(G.nodes())) * 100, end="\r")

    return three_node_hash


def main():
    two_node_hash_table = get_two_node_hash_table()
    fly_ppi_path = Path("data/fly_ppi.csv")
    fly_reg_path = Path("data/fly_reg.csv")
    bsub_ppi_path = Path("data/bsub_ppi.csv")
    bsub_reg_path = Path("data/bsub_reg.csv")
    ceravisiae_ppi_path = Path("data/ceravisiae_ppi.csv")
    ceravisiae_reg_path = Path("data/ceravisiae_reg.csv")
    drerio_ppi_path = Path("data/drerio_ppi.csv")
    drerio_reg_path = Path("data/drerio_reg.csv")
    elegans_ppi_path = Path("data/elegans_ppi.csv")
    elegans_reg_path = Path("data/elegans_reg.csv")
    protein_id_dict = get_protein_id_dict(fly_ppi_path, fly_reg_path)
    G = read_csv(
        fly_ppi_path,
        fly_reg_path,
        protein_id_dict,
        node_size_limit=999999999,
        edge_size_limit=999999999,
    )

    print(f"Number of nodes: {len(G.nodes())}")
    print(f"Number of edges: {len(G.edges())}")

    # two_node_hash_table = get_two_node_graphlet_dist_adj_matrix(G, two_node_hash_table)
    two_node_hash_table, two_node_orbit_dict = get_two_node_graphlet_stats(
        G, two_node_hash_table
    )
    print("\ntwo node graphlet counts")
    for key in two_node_hash_table:
        print(f"{key} = {two_node_hash_table[key]}")
    print("\ntwo node graphlet orbit counts")
    for key in two_node_orbit_dict:
        print(f"{key} = {len(two_node_orbit_dict[key])}")
    # three_node_hash_table = get_three_node_graphlet_dist_adj_list(G)
    # print(three_node_hash_table)

    # draw_labeled_multigraph(G, "label")
    plt.show()


if __name__ == "__main__":
    main()
