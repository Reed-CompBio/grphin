import sys
from matplotlib import pyplot as plt
from pathlib import Path
import networkx as nx
import numpy as np
import random
import scipy as sp
import csv
import curses
from itertools import combinations
import time
import csv
import ast
import re
from scipy.stats import norm
import seaborn as sns
from collections import Counter, defaultdict

csv.field_size_limit(sys.maxsize)


def get_two_node_dict():
    """
    Get all 2-node graphlet binary edge vectors
    see diagram for reference
    """
    # 2-node graphlet binary edge vectors
    # see diagram for reference
    a_1 = (1, 0, 0)
    a_2 = (1, 0, 0)
    a_hash = hash(a_1 + a_2)

    b_1 = (0, 1, 0)
    b_2 = (0, 0, 1)
    b_hash = hash(b_1 + b_2)

    c_1 = (1, 1, 0)
    c_2 = (1, 0, 1)
    c_hash = hash(c_1 + c_2)

    d_1 = (1, 1, 1)
    d_2 = (1, 1, 1)
    d_hash = hash(d_1 + d_2)

    e_1 = (0, 1, 1)
    e_2 = (0, 1, 1)
    e_hash = hash(e_1 + e_2)

    two_node_graphlet_dict = {a_hash: 0, b_hash: 0, c_hash: 0, d_hash: 0, e_hash: 0}
    two_node_graphlet_labels = {a_hash: 1, b_hash: 2, c_hash: 3, d_hash: 4, e_hash: 5}
    return two_node_graphlet_dict, two_node_graphlet_labels


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


import networkx as nx


def simplify_graph_to_undirected(G):
    """
    Simplify the given graph to a new undirected graph with only one edge per connected node pair.

    Parameters:
    G (networkx.Graph): The original graph (can be directed or undirected, multigraph or simple).

    Returns:
    G_prime (networkx.Graph): The simplified undirected graph.
    """
    # Create an empty undirected graph
    G_prime = nx.Graph()

    # Add all nodes to G_prime
    G_prime.add_nodes_from(G.nodes())

    # Add a single undirected edge for every connected node pair in G
    for u, v in G.edges():
        G_prime.add_edge(u, v)

    return G_prime


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
    print()

    return G


def load_graphlet_config(file_path):
    """Load graphlet lookup table from a CSV file."""
    graphlet_config = []
    with open(file_path, mode="r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            graphlet_config.append(
                {
                    "key": ast.literal_eval(row["key"]),
                    "a_expected": ast.literal_eval(row["a_expected"]),
                    "b_expected": ast.literal_eval(row["b_expected"]),
                    "c_expected": ast.literal_eval(row["c_expected"]),
                    "orbits": (
                        int(row["orbit1"]),
                        int(row["orbit2"]),
                        int(row.get("orbit3", -1)),
                    ),  # Handle missing orbit3
                }
            )
    return graphlet_config


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


def get_two_node_graphlet_stats(G, two_node_graphlet_dict):
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
            if hash(tuple(vectors)) in two_node_graphlet_dict:
                two_node_graphlet_dict[hash(tuple(vectors))] += 1
    return two_node_graphlet_dict, orbits_dict


def get_two_node_adjacency_list(G):
    """Get the adjacency list for a MultiDiGraph"""
    print("getting adjacency list")

    adj_list_vector = [{} for _ in range(len(G.nodes()))]

    for i, j, data in G.edges(data=True):
        label = data.get("label")
        if label == "ppi":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [1, 0, 0]
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [1, 0, 0]
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


def hash_tuple(xy):
    return hash(tuple(xy))


def get_three_node_graphlet_dict(ab, ac, ba, bc, ca, cb):
    # Process into a hashed tuple
    ab = hash_tuple(ab)
    ac = hash_tuple(ac)
    ba = hash_tuple(ba)
    bc = hash_tuple(bc)
    ca = hash_tuple(ca)
    cb = hash_tuple(cb)

    # Define the dictionary
    my_dict = {
        hash((0, 0, 0)): 0,
        hash((1, 0, 0)): 1,
        hash((0, 1, 0)): 2,
        hash((0, 0, 1)): 3,
        hash((1, 1, 0)): 4,
        hash((1, 0, 1)): 5,
        hash((0, 1, 1)): 6,
        hash((1, 1, 1)): 7,
    }

    return my_dict[ab], my_dict[ac], my_dict[ba], my_dict[bc], my_dict[ca], my_dict[cb]


def get_three_node_graphlet_dist_adj_list_v3(G: nx.MultiDiGraph, G_prime: nx.Graph):
    three_node_graphlet_dict = {}
    graphlet_mapper = {}
    orbit_dict = {}
    orbit_mapper = {}
    graphlet_config = load_graphlet_config("graphlet_config.csv")
    start_time = time.time()

    # create all the binary edge vectors
    adj_list_vector = [{} for _ in range(len(G.nodes()))]
    for i, j, data in G.edges(data=True):
        label = data.get("label")
        if label == "ppi":
            if j not in adj_list_vector[i]:
                adj_list_vector[i][j] = [1, 0, 0]
            if i not in adj_list_vector[j]:
                adj_list_vector[j][i] = [1, 0, 0]
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
    three_node_combination = set()  # Use a set for fast triplet lookups

    # Preprocess neighbors for each node once
    neighbors_dict = {i: set(G_prime.neighbors(i)) for i in G_prime.nodes()}
    completed_i = set()

    for i in G_prime.nodes():
        print(f"Node: {i}/{len(G_prime.nodes)}", end="\r")
        for j in neighbors_dict[i]:
            for k in neighbors_dict[j].difference(completed_i):
                if (
                    (i < k) and (i != j) and (j != k)
                ):  # Ensure no duplicates by enforcing i < k and i != j
                    triplet = tuple(sorted([i, j, k]))
                    if triplet not in three_node_combination:
                        three_node_combination.add(triplet)

                        a = i
                        b = j
                        c = k

                        ab = ac = ba = bc = ca = cb = 0
                        if b in adj_list_vector[a]:
                            ab = adj_list_vector[a][b]
                        else:
                            ab = [0, 0, 0]
                        if c in adj_list_vector[a]:
                            ac = adj_list_vector[a][c]
                        else:
                            ac = [0, 0, 0]
                        if a in adj_list_vector[b]:
                            ba = adj_list_vector[b][a]
                        else:
                            ba = [0, 0, 0]
                        if c in adj_list_vector[b]:
                            bc = adj_list_vector[b][c]
                        else:
                            bc = [0, 0, 0]
                        if a in adj_list_vector[c]:
                            ca = adj_list_vector[c][a]
                        else:
                            ca = [0, 0, 0]
                        if b in adj_list_vector[c]:
                            cb = adj_list_vector[c][b]
                        else:
                            cb = [0, 0, 0]
                        a_b, a_c, b_a, b_c, c_a, c_b = get_three_node_graphlet_dict(
                            ab, ac, ba, bc, ca, cb
                        )

                        # order A, B, C edge values internally
                        a_edges = tuple(sorted([a_b, a_c]))
                        b_edges = tuple(sorted([b_a, b_c]))
                        c_edges = tuple(sorted([c_a, c_b]))

                        # Create a list of tuples in order [A, B, C]
                        tuples_list = [a_edges, b_edges, c_edges]

                        # Sort the tuples first by the first index, then by the second index
                        sorted_tuples = tuple(
                            sorted(tuples_list, key=lambda x: (x[0], x[1]))
                        )

                        # Add the graphlet if it has not been seen yet and update the count for the graphlet
                        if hash(sorted_tuples) not in three_node_graphlet_dict:
                            three_node_graphlet_dict[hash(sorted_tuples)] = 0
                            graphlet_mapper[hash(sorted_tuples)] = sorted_tuples
                        three_node_graphlet_dict[hash(sorted_tuples)] += 1
                        # print(f'Graphlet: {sorted_tuples}')
                        orbit_dict = get_orbit_per_graphlet(
                            orbit_dict,
                            orbit_mapper,
                            sorted_tuples,
                            a_edges,
                            b_edges,
                            c_edges,
                            i,
                            j,
                            k,
                            graphlet_config,
                        )

        # Once we're done processing i, mark it as completed
        completed_i.add(i)
        # print(f"Completed node: {i}")

    run_time = time.time() - start_time
    print("run time : %.3f seconds" % run_time)
    return three_node_graphlet_dict, graphlet_mapper, orbit_dict, orbit_mapper, run_time


def get_orbit_per_graphlet(
    orbit_dict,
    orbit_mapper,
    sorted_tuples,
    a_edges,
    b_edges,
    c_edges,
    i,
    j,
    k,
    graphlet_config,
):
    """Assign orbit information based on the graphlet configuration."""
    for graphlet in graphlet_config:
        if sorted_tuples == graphlet["key"]:
            orbit_change = get_orbit_position_change(
                a_edges,
                b_edges,
                c_edges,
                graphlet["a_expected"],
                graphlet["b_expected"],
                graphlet["c_expected"],
                i,
                j,
                k,
            )
            for idx, orbit in enumerate(graphlet["orbits"]):
                if orbit == -1:  # Skip missing orbits
                    continue
                graphlet_key = (sorted_tuples, orbit)
                if hash(graphlet_key) not in orbit_dict:
                    orbit_dict[hash(graphlet_key)] = []
                    orbit_mapper[hash(graphlet_key)] = f"{graphlet_key}, Orbit {orbit}"
                orbit_dict[hash(graphlet_key)] += [orbit_change[idx]]

    return orbit_dict


def get_orbit_position_change(
    a_edges, b_edges, c_edges, a_expected, b_expected, c_expected, i, j, k
):
    # a, b, c
    if a_edges == a_expected and b_edges == b_expected and c_edges == c_expected:
        return [i, j, k]
    # a, c, b
    if a_edges == a_expected and b_edges == c_expected and c_edges == b_expected:
        return [i, k, j]
    # b, a, c
    if a_edges == b_expected and b_edges == a_expected and c_edges == c_expected:
        return [j, i, k]
    # b, c, a
    if a_edges == b_expected and b_edges == c_expected and c_edges == a_expected:
        return [j, k, i]
    # c, a, b
    if a_edges == c_expected and b_edges == a_expected and c_edges == b_expected:
        return [k, i, j]
    # c, b, a
    if a_edges == c_expected and b_edges == b_expected and c_edges == a_expected:
        return [k, j, i]
    raise ValueError(
        f"Unexpected edges configuration: a_edges={a_edges}, b_edges={b_edges}, c_edges={c_edges}"
    )


def plot_three_node_graphlet_distribution(
    graphlet_dict, graphlet_mapper, indexed_graphlet_dict, selected_network, output_dir
):
    print("plotting graphlet distribution")
    hist_data = []
    x_label = [*range(0, len(indexed_graphlet_dict), 1)]
    for graphlet in indexed_graphlet_dict:
        if graphlet in graphlet_dict:
            print(
                indexed_graphlet_dict[graphlet],
                graphlet_mapper[graphlet],
                graphlet_dict[graphlet],
            )
            hist_data.append(graphlet_dict[graphlet])
        else:
            hist_data.append(0)
    hist_data = [value if value > 0 else 0.1 for value in hist_data]

    fig = plt.figure(figsize=(14, 6))
    plt.bar(x_label, hist_data, color="skyblue", edgecolor="black")
    plt.yscale("log")
    plt.title(f"{selected_network} Graphlet Count Distribution", fontsize=16)
    plt.xlabel("Graphlet Index", fontsize=14)
    plt.ylabel("Count (log scale)", fontsize=14)
    plt.xticks(x_label[::2], fontsize=8)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/three_node_graphlet_dist.pdf")
    plt.show()

    sorted_graphlet_dict = {
        key: value
        for key, value in sorted(
            graphlet_dict.items(), key=lambda item: item[1], reverse=True
        )
    }

    with open(f"{output_dir}/top_graphlet_counts.csv", "w") as f:
        f.write(f"graphlet\tid\tcount\n")
        for graphlet in sorted_graphlet_dict:
            f.write(
                f"{graphlet_mapper[graphlet]}\t{indexed_graphlet_dict[graphlet]}\t{sorted_graphlet_dict[graphlet]}\n"
            )

    return None

def plot_three_node_orbit_dist(orbit_dict, orbit_mapper, indexed_orbit_dict, selected_network, output_dir):
    orbit_counts = {}
    for orbit in orbit_dict:
        count = 0
        for protein_set in orbit_dict[orbit]:
            count += protein_set[1]
        orbit_counts[orbit] = count

    sorted_orbit_dict = {
        key: value
        for key, value in sorted(
            orbit_counts.items(), key=lambda item: item[1], reverse=True
        )
    }

    with open(f"{output_dir}/top_orbit_counts.csv", "w") as f:
        f.write(f"orbit\tid\tcount\n")
        for orbit in sorted_orbit_dict:
            f.write(
                f"{orbit_mapper[orbit]}\t{indexed_orbit_dict[orbit]}\t{orbit_counts[orbit]}\n"
            )
    return None


def get_stress_proteins(protein_id_dict, stress_path, delimiter):
    stress_proteins = []

    with open(stress_path, "r") as file:
        next(file)
        reader = csv.reader(file, delimiter=delimiter)
        for line in reader:
            if line[0] in protein_id_dict:
                stress_proteins.append(protein_id_dict[line[0]])
    return stress_proteins


def plot_stress_orbit_distribution(
    orbit_dict,
    orbit_mapper,
    indexed_orbit_dict,
    stress_proteins_list,
    protein_id_dict,
    selected_network,
    output_dir,
    node_orbit_arr
):

    observed_median_count = {}
    # print("Stress proteins : ", stress_proteins_list)
    stress_protein_orbit_dict = {}
    print("Counting stress proteins\n")
    for orbit in orbit_dict:
        # print(f"orbit {orbit}", end="\r")
        stress_protein_orbit_dict[orbit] = []
        for protein in stress_proteins_list:
            protein_orbit_count = node_orbit_arr[protein][int(indexed_orbit_dict[orbit])]
            stress_protein_orbit_dict[orbit] += [protein_orbit_count]

    print("Counting stress medians\n")
    for orbit in stress_protein_orbit_dict:
        # print(f"orbit {orbit}", end="\r")
        orbit_list = []
        if len(stress_protein_orbit_dict[orbit]) == 0:
            orbit_list = [0]
        else:
            orbit_list = stress_protein_orbit_dict[orbit]
        sorted_list = np.sort(orbit_list)
        median = np.median(sorted_list)
        observed_median_count[orbit] = median

    sample = 1000
    sample_size = len(stress_proteins_list)
    print("getting non stress proteins \n")
    non_stress_proteins = []
    for protein in protein_id_dict:
        # print(f"protein {protein}", end="\r")
        if protein_id_dict[protein] not in stress_proteins_list:
            non_stress_proteins.append(protein_id_dict[protein])

    sample_results = {}

    for i in range(0,sample):
        print("sample :", i, end="\r")
        non_stress_sample = random.sample(non_stress_proteins, sample_size)
        array_stack = []
        for protein in non_stress_sample:
            array_stack.append(node_orbit_arr[protein])
        stacked = np.vstack(array_stack)
        orbit_medians_list = np.median(stacked, axis = 0)

        for orbit in orbit_dict:
            orbit_index = indexed_orbit_dict[orbit]
            if orbit not in sample_results:
                sample_results[orbit] = []
            sample_results[orbit] += [orbit_medians_list[int(orbit_index)]]




    # sample_results = {}
    # for i in range(0, sample):
    #     print("sample : ", i)
    #     # randomly select n_number of non-stress proteins
    #     non_stress_sample = random.sample(non_stress_proteins, sample_size)
    #     sample_non_stress_orbit_dict = {}
    #     # if orbit count is zero, lets skip this,

    #     for orbit in orbit_dict:
    #         print("orbit", orbit)
    #         sample_non_stress_orbit_dict[orbit] = []
    #         for protein_id in non_stress_sample:
    #             sample_non_stress_orbit_dict[orbit] += [
    #                 orbit_dict[orbit].count(protein_id)
    #             ]

    #     for orbit in orbit_dict:
    #         orbit_list = []
    #         if len(sample_non_stress_orbit_dict[orbit]) == 0:
    #             orbit_list = [0]
    #         else:
    #             orbit_list = sample_non_stress_orbit_dict[orbit]
    #         sorted_list = np.sort(orbit_list)
    #         median = np.median(sorted_list)
    #         if orbit not in sample_results:
    #             sample_results[orbit] = []
    #         sample_results[orbit] += [median]
    #         # orbit 1 = [1, 4, 5, 6, 7, 8 .. samples_samples]

    significance = {}

    # Calculate if an orbit is significant
    for orbit in orbit_dict:
        count = 0
        # from the vector of all the random medians at a given orbit
        for random_median in sample_results[orbit]:
            if observed_median_count[orbit] > random_median:
                count += 1

        if count >= float(sample) * 0.99:
            significance[orbit] = 1 
            print(indexed_orbit_dict[orbit], ": significant")
        else:
            significance[orbit] = 0
            print(indexed_orbit_dict[orbit], ": NOT significant")

    with open(f"{output_dir}/stress_orbit_significance.csv", "w") as f:
            f.write(f"orbit_id\tsignificant?\tobserved_median\tvector_random_medians\n")
            for orbit in significance:
                f.write(f"{indexed_orbit_dict[orbit]}\t{significance[orbit]}\t{observed_median_count[orbit]}\t{sample_results[orbit]}\n")


    i = 0
    for orbit in orbit_dict:
        if significance[orbit] == 1:
            hist_data = sample_results[orbit]

            fig = plt.figure(figsize=(14, 6))
            plt.hist(hist_data)
            plt.axvline(x=int(observed_median_count[orbit]), color='r', linestyle='--')
            plt.title(
                f"{selected_network} Random Median Samples at Orbit {int(indexed_orbit_dict[orbit])} Distribution",
                fontsize=16,
            )
            plt.savefig(f"{output_dir}/sig_orbit/orbit{i}.pdf")

            # x_label = [*range(0, sample_size, 1)]
            # fig = plt.figure(figsize=(14, 6))
            # plt.bar(x_label, hist_data, color="skyblue", edgecolor="black")
            # # plt.yscale("log")
            # plt.title(
            #     f"{selected_network} Stress Protein Significance per Orbit Distribution",
            #     fontsize=16,
            # )
            # plt.xlabel("Orbit Index", fontsize=14)
            # plt.ylabel("Count (log scale)", fontsize=14)
            # plt.xticks(x_label[::2], fontsize=8)
            # plt.grid(axis="y", linestyle="--", alpha=0.7)
            # plt.tight_layout()
            # plt.savefig(f"{output_dir}/stres_orbit_significance_dist_v2.pdf")
            # plt.show()
            plt.close()
        i+=1


    # hist_data = []
    # x_label = [*range(0, len(indexed_orbit_dict), 1)]
    # for orbit in indexed_orbit_dict:
    #     hist_data.append(significance[orbit])

    # fig = plt.figure(figsize=(14, 6))
    # plt.bar(x_label, hist_data, color="skyblue", edgecolor="black")
    # plt.yscale("log")
    # plt.title(
    #     f"{selected_network} Stress Protein Significance per Orbit Distribution",
    #     fontsize=16,
    # )
    # plt.xlabel("Orbit Index", fontsize=14)
    # plt.ylabel("Count (log scale)", fontsize=14)
    # plt.xticks(x_label[::2], fontsize=8)
    # plt.grid(axis="y", linestyle="--", alpha=0.7)
    # plt.tight_layout()
    # plt.savefig(f"{output_dir}/stres_orbit_significance_dist_v2.pdf")
    # plt.show()

    # with open(f"{output_dir}/stress_orbit_significance.csv", "w") as f:
    #         for orbit in significance:
    #             f.write(f"{indexed_orbit_dict[orbit]}\t{significance[orbit]}\n")

    return None


def species_wide_two_node_plots():
    species_list = ["bsub", "fly", "cerevisiae", "drerio", "elegans"]
    graphlet_ids = ["G_1", "G_2", "G_3", "G_4", "G_5"]
    species_graphlet_data = {}

    for species in species_list:
        print(f"Processing species: {species}")
        species_graphlet_data[species] = {}
        with open(f"final_output/{species}/two_node_graphlet_counts.csv", "r") as file:
            csv_reader = csv.reader(file, delimiter=",")
            for row in csv_reader:
                species_graphlet_data[species][row[0]] = int(row[1].strip())

    graphlet_counts = {graphlet: [] for graphlet in graphlet_ids}

    for species_name in species_list:
        species_counts = species_graphlet_data[species_name]
        total_count = sum(species_counts[graphlet] for graphlet in graphlet_ids)
        for graphlet in graphlet_ids:
            percentage = (species_counts[graphlet] / total_count) * 100 if total_count > 0 else 0
            graphlet_counts[graphlet].append(percentage)

    index = np.arange(len(species_list))
    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.5
    bottom_values = np.zeros(len(species_list))
    colors = ["blue", "green", "red", "yellow", "pink"]

    for i, graphlet in enumerate(graphlet_ids):
        bar = ax.bar(
            index,
            graphlet_counts[graphlet],
            bar_width,
            bottom=bottom_values,
            label=graphlet,
            color=colors[i]
        )
        bottom_values += np.array(graphlet_counts[graphlet])

    ax.set_xlabel("Species")
    ax.set_ylabel("Percentage of Graphlet Counts")
    ax.set_title("Stacked Bar Chart of Graphlet Percentages per Species")
    ax.set_xticks(index)
    ax.set_xticklabels(species_list)
    ax.legend()

    plt.tight_layout()
    plt.savefig("final_output/species_two_node_graphlet_dist_percentage.pdf")
    plt.show()
    plt.close()


    orbit_ids = ["1", "2", "3", "4", "5", "6", "7"]
    species_orbit_data = {}

    for species in species_list:
        print(f"Processing species: {species}")
        species_orbit_data[species] = {}
        with open(f"final_output/{species}/two_node_orbit_counts.csv", "r") as file:
            csv_reader = csv.reader(file, delimiter=",")
            for row in csv_reader:
                print(row[0], row[1].strip())
                species_orbit_data[species][row[0]] = int(row[1].strip())

    orbit1_count = []
    orbit2_count = []
    orbit3_count = []
    orbit4_count = []
    orbit5_count = []
    orbit6_count = []
    orbit7_count = []

    for species_name in species_list:
        species_counts = species_orbit_data[species_name]
        total_count = sum(species_counts.values())
        orbit1_count.append(species_counts["1"] / total_count * 100)
        orbit2_count.append(species_counts["2"] / total_count * 100)
        orbit3_count.append(species_counts["3"] / total_count * 100)
        orbit4_count.append(species_counts["4"] / total_count * 100)
        orbit5_count.append(species_counts["5"] / total_count * 100)
        orbit6_count.append(species_counts["6"] / total_count * 100)
        orbit7_count.append(species_counts["7"] / total_count * 100)

    orbit1_count = np.array(orbit1_count)
    orbit2_count = np.array(orbit2_count)
    orbit3_count = np.array(orbit3_count)
    orbit4_count = np.array(orbit4_count)
    orbit5_count = np.array(orbit5_count)
    orbit6_count = np.array(orbit6_count)
    orbit7_count = np.array(orbit7_count)

    index = np.arange(len(species_list))

    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.5
    bar1 = ax.bar(index, orbit1_count, bar_width, label="1", color="blue")
    bar2 = ax.bar(
        index, orbit2_count, bar_width, bottom=orbit1_count, label="2", color="green"
    )
    bar3 = ax.bar(
        index,
        orbit3_count,
        bar_width,
        bottom=orbit1_count + orbit2_count,
        label="3",
        color="red",
    )
    bar4 = ax.bar(
        index,
        orbit4_count,
        bar_width,
        bottom=orbit1_count + orbit2_count + orbit3_count,
        label="4",
        color="yellow",
    )
    bar5 = ax.bar(
        index,
        orbit5_count,
        bar_width,
        bottom=orbit1_count + orbit2_count + orbit3_count + orbit4_count,
        label="5",
        color="pink",
    )
    bar6 = ax.bar(
        index,
        orbit6_count,
        bar_width,
        bottom=orbit1_count + orbit2_count + orbit3_count + orbit4_count + orbit5_count,
        label="6",
        color="orange",
    )
    bar7 = ax.bar(
        index,
        orbit7_count,
        bar_width,
        bottom=orbit1_count
        + orbit2_count
        + orbit3_count
        + orbit4_count
        + orbit5_count
        + orbit6_count,
        label="7",
        color="purple",
    )

    ax.set_xlabel("Species")
    ax.set_ylabel("Percentage of orbit ids")
    ax.set_title("Stacked Bar Chart of Orbit Percentages per Species")
    ax.set_xticks(index)
    ax.set_xticklabels(species_list)
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"final_output/species_two_node_orbit_dist.pdf")
    plt.show()
    plt.close()

    return None


def species_wide_3_node_plots():
    species_list = ["bsub", "fly", "cerevisiae", "drerio", "elegans"]
    graphlet_counts = defaultdict(int)
    species_graphlet_counts = {}
    for species in species_list:
        print(f"Processing species: {species}")
        with open(f"final_output/{species}/top_graphlet_counts.csv", "r") as file:
            csv_reader = csv.reader(file, delimiter="\t")
            next(csv_reader)
            for row in csv_reader:
                graphlet_id = row[1]
                count = int(row[2])
                graphlet_counts[graphlet_id] += count
                if species not in species_graphlet_counts:
                    species_graphlet_counts[species] = {}
                if graphlet_id not in species_graphlet_counts[species]:
                    species_graphlet_counts[species][graphlet_id] = count

    top_graphlets = sorted(graphlet_counts.items(), key=lambda x: x[1], reverse=True)[
        :5
    ]

    print("Top 10 graphlets with the highest counts across all species:")
    for graphlet_id, count in top_graphlets:
        print(f"Graphlet ID: {graphlet_id}, Total Count: {count}")
        for species in species_list:
            print(f"{species} Count: {species_graphlet_counts[species][graphlet_id]}")

    with open(f"final_output/species_top_three_node_graphlet_counts.csv", "w") as f:
        f.write(f"graphlet_id\ttotal_count\tbsub\tcerevisiae\tdrerio\telegans\tfly\n")
        for graphlet in top_graphlets:
            f.write(
                f"{graphlet[0]}\t{species_graphlet_counts["bsub"][graphlet_id]}\t{species_graphlet_counts["cerevisiae"][graphlet_id]}\t{species_graphlet_counts["drerio"][graphlet_id]}\t{species_graphlet_counts["elegans"][graphlet_id]}\t{species_graphlet_counts["bsub"][graphlet_id]}\n"
            )
    return None


def dropdown_menu(stdscr, options):
    curses.start_color()
    curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)
    curses.curs_set(0)
    current_index = 0

    while True:
        stdscr.clear()

        for i, option in enumerate(options):
            x = 2
            y = 1 + i
            if i == current_index:
                stdscr.attron(curses.color_pair(1))
                stdscr.addstr(y, x, option.ljust(20))
                stdscr.attroff(curses.color_pair(1))
            else:
                stdscr.addstr(y, x, option.ljust(20))

        stdscr.refresh()

        key = stdscr.getch()

        if key == curses.KEY_UP and current_index > 0:
            current_index -= 1
        elif key == curses.KEY_DOWN and current_index < len(options) - 1:
            current_index += 1
        elif key in [curses.KEY_ENTER, 10, 13]:
            return options[current_index]
        elif key == 27:
            return None


def type_answer(stdscr, step):
    curses.start_color()
    curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)
    curses.curs_set(0)
    x = 0
    y = 1
    typed_number = ""

    while True:
        stdscr.clear()
        if step == 0:
            stdscr.addstr(0, 0, "Node limit? (Press ENTER for no limit)")
        elif step == 1:
            stdscr.addstr(0, 0, "Edge limit? (Press ENTER for no limit)")

        stdscr.addstr(y, x, typed_number)
        stdscr.refresh()

        key = stdscr.getkey()

        # Check if key is a digit (0-9)
        if key.isdigit():
            typed_number += key
        elif key in ("\n", "\r") and typed_number == "" or int(typed_number) == 0:
            return 999999999
        elif key in ("\n", "\r") and int(typed_number) != 0:
            return int(typed_number)


def get_user_inputs(selected_network, selected_graphlet, selected_processing, selected_network_type):
    ppi_path = None
    reg_path = None
    graphlet_option = None
    output_dir = None
    stress_data_path = None
    process_type = None

    txids = ["txid6239", "txid7227", "txid7955", "txid224308", "txid559292"]

    match selected_network:
        case "D. melanogaster":
            if selected_network_type == "Normal":
                ppi_path = Path("data/fly_ppi.csv")
                reg_path = Path("data/fly_reg.csv")
                output_dir = Path("final_output/fly")
            else:
                ppi_path = Path("data/oxidative_stress/txid7227/stress_ppi.csv")
                reg_path = Path("data/oxidative_stress/txid7227/stress_reg.csv")
                output_dir = Path("final_output/fly_stress")
            stress_data_path = Path(
                "data/oxidative_stress/txid7227/txid7227-stress-proteins.csv"
            )
            process_type = get_process_type(selected_processing)
        case "B. subtilis":
            if selected_network_type == "Normal":
                ppi_path = Path("data/bsub_ppi.csv")
                reg_path = Path("data/bsub_reg.csv")
                output_dir = Path("final_output/bsub")
            else:
                ppi_path = Path("data/oxidative_stress/txid224308/stress_ppi.csv")
                reg_path = Path("data/oxidative_stress/txid224308/stress_reg.csv")
                output_dir = Path("final_output/bsub_stress")
            stress_data_path = Path(
                "data/oxidative_stress/txid224308/txid224308-stress-proteins.csv"
            )
            process_type = get_process_type(selected_processing)
        case "S. cerevisiae":
            if selected_network_type == "Normal":
                ppi_path = Path("data/cerevisiae_ppi.csv")
                reg_path = Path("data/cerevisiae_reg.csv")
                output_dir = Path("final_output/cerevisiae")
            else: 
                ppi_path = Path("data/oxidative_stress/txid559292/stress_ppi.csv")
                reg_path = Path("data/oxidative_stress/txid559292/stress_reg.csv")
                output_dir = Path("final_output/cerevisiae_stress")
            stress_data_path = Path(
                "data/oxidative_stress/txid559292/txid559292-stress-proteins.csv"
            )
            process_type = get_process_type(selected_processing)
        case "D. rerio":
            if selected_network_type == "Normal":
                ppi_path = Path("data/drerio_ppi.csv")
                reg_path = Path("data/drerio_reg.csv")
                output_dir = Path("final_output/drerio")
            else:
                ppi_path = Path("data/oxidative_stress/txid7955/stress_ppi.csv")
                reg_path = Path("data/oxidative_stress/txid7955/stress_reg.csv")
                output_dir = Path("final_output/drerio_stress")
            stress_data_path = Path(
                "data/oxidative_stress/txid7955/txid7955-stress-proteins.csv"
            )
            process_type = get_process_type(selected_processing)
        case "C. elegans":
            if selected_network_type == "Normal":
                ppi_path = Path("data/elegans_ppi.csv")
                reg_path = Path("data/elegans_reg.csv")
                output_dir = Path("final_output/elegans")
            else:
                ppi_path = Path("data/oxidative_stress/txid6239/stress_ppi.csv")
                reg_path = Path("data/oxidative_stress/txid6239/stress_reg.csv")
                output_dir = Path("final_output/elegans_stress")
            stress_data_path = Path(
                "data/oxidative_stress/txid6239/txid6239-stress-proteins.csv"
            )
            process_type = get_process_type(selected_processing)
        case "Test network":
            ppi_path = Path("data/test_ppi.csv")
            reg_path = Path("data/test_reg.csv")
            output_dir = Path("final_output/test")
            process_type = get_process_type(selected_processing)
        case "Toy network":
            ppi_path = Path("data/toy_ppi.csv")
            reg_path = Path("data/toy_reg.csv")
            output_dir = Path("final_output/toy")
            process_type = get_process_type(selected_processing)
        case "Shuffled toy network":
            ppi_path = Path("data/toy2_ppi.csv")
            reg_path = Path("data/toy2_reg.csv")
            output_dir = Path("final_output/toy2")
            process_type = get_process_type(selected_processing)
        case "2-node species wide":
            output_dir = Path("final_output")
            process_type = 2
        case "3-node species wide":
            output_dir = Path("final_output")
            process_type = 3

    if selected_graphlet == "2-node":
        graphlet_option = 2
    elif selected_graphlet == "3-node":
        graphlet_option = 3

    return (
        ppi_path,
        reg_path,
        graphlet_option,
        output_dir,
        stress_data_path,
        process_type,
    )

def get_process_type(selected_processing):
    if selected_processing == "all process":
        return 0
    elif selected_processing == "post-processing only":
        return 1

def main(stdscr):
    try:
        network = [
            "D. melanogaster",
            "B. subtilis",
            "S. cerevisiae",
            "D. rerio",
            "C. elegans",
            "2-node species wide",
            "3-node species wide",
            "Test network",
            "Toy network",
            "Shuffled toy network",
            "Exit",
        ]
        selected_network = dropdown_menu(stdscr, network)
        stdscr.clear()

        graphlet = ["2-node", "3-node", "Exit"]
        selected_graphlet = dropdown_menu(stdscr, graphlet)
        stdscr.clear()

        processing = ["all process", "post-processing only"]
        selected_processing = dropdown_menu(stdscr, processing)
        stdscr.clear()

        network_type = ["Normal", "Stress"]
        selected_network_type = dropdown_menu(stdscr, network_type)
        stdscr.clear()

        node_limit = type_answer(stdscr, 0)
        stdscr.clear()

        edge_limit = type_answer(stdscr, 1)
        stdscr.clear()
    except Exception as e:
        stdscr.addstr(0, 0, f"An error occurred: {str(e)}")
        stdscr.getch()
    finally:
        curses.endwin()
    ppi_path, reg_path, graphlet_mode, output_dir, stress_data_path, process_type = (
        get_user_inputs(selected_network, selected_graphlet, selected_processing, selected_network_type)
    )

    if graphlet_mode == 2:
        if process_type == 2:
            print("species wide processing")
            species_wide_two_node_plots()
        else:
            two_node_graphlet_dict, two_node_graphlet_labels = get_two_node_dict()
            protein_id_dict = get_protein_id_dict(ppi_path, reg_path)

            G = read_csv(
                ppi_path,
                reg_path,
                protein_id_dict,
                node_size_limit=node_limit,
                edge_size_limit=edge_limit,
            )
            print(selected_network)
            print(f"Number of nodes: {len(G.nodes())}")
            print(f"Number of edges: {len(G.edges())}")

            # Simplify the graph G to G_prime
            G_prime = simplify_graph_to_undirected(G)

            # Print the simplified graph's details
            print(f"\nSimplified Graph:")
            print(f"Number of nodes: {len(G_prime.nodes())}")
            print(f"Number of edges: {len(G_prime.edges())}")
            two_node_graphlet_dict, two_node_orbit_dict = get_two_node_graphlet_stats(
                G, two_node_graphlet_dict
            )
            print("\ntwo node graphlet counts")
            for key in two_node_graphlet_labels:
                print(
                    f"G_{two_node_graphlet_labels[key]} = {two_node_graphlet_dict[key]}"
                )
            print("\ntwo node graphlet orbit counts")
            for key in two_node_orbit_dict:
                print(f"{key} = {len(two_node_orbit_dict[key])}")

            with open(f"{output_dir}/two_node_graphlet_counts.csv", "w") as f:
                for key in two_node_graphlet_labels:
                    f.write(
                        f"G_{two_node_graphlet_labels[key]}, {two_node_graphlet_dict[key]}\n"
                    )
            f.close()

            with open(f"{output_dir}/two_node_orbit_counts.csv", "w") as f:
                for key in two_node_orbit_dict:
                    f.write(f"{key}, {len(two_node_orbit_dict[key])}\n")
            f.close()

    elif graphlet_mode == 3:
        if process_type == 3:
            print("species wide processing")
            species_wide_3_node_plots()

        else:
            two_node_graphlet_dict, two_node_graphlet_labels = get_two_node_dict()
            protein_id_dict = get_protein_id_dict(ppi_path, reg_path)

            G = read_csv(
                ppi_path,
                reg_path,
                protein_id_dict,
                node_size_limit=node_limit,
                edge_size_limit=edge_limit,
            )
            print(selected_network)
            print(f"Number of nodes: {len(G.nodes())}")
            print(f"Number of edges: {len(G.edges())}")

            # Simplify the graph G to G_prime
            G_prime = simplify_graph_to_undirected(G)

            # Print the simplified graph's details
            print(f"\nSimplified Graph:")
            print(f"Number of nodes: {len(G_prime.nodes())}")
            print(f"Number of edges: {len(G_prime.edges())}")
            graphlet_config = load_graphlet_config("graphlet_config.csv")

            three_node_graphlet_dict = {}
            graphlet_mapper = {}
            orbit_dict = {}
            orbit_mapper = {}
            run_time = 0
            indexed_graphlet_dict = {}
            indexed_orbit_dict = {}
            orbit_id_dict = {}
            node_orbit_arr = np.zeros(
                (len(G.nodes), len(indexed_orbit_dict)), dtype=int
            )
            orbit_count = 0
            graphlet_count = 0
            for graphlet in graphlet_config:
                indexed_graphlet_dict[hash(graphlet["key"])] = graphlet_count
                for i in set(graphlet["orbits"]):
                    indexed_orbit_dict[hash((graphlet["key"], i))] = orbit_count
                    orbit_count += 1
                graphlet_count += 1
            # if we want to compute the whole counts
            if process_type == 0:
                (
                    three_node_graphlet_dict,
                    graphlet_mapper,
                    orbit_dict,
                    orbit_mapper,
                    run_time,
                ) = get_three_node_graphlet_dist_adj_list_v3(G, G_prime)
                node_orbit_arr = np.zeros(
                    (len(G.nodes), len(indexed_orbit_dict)), dtype=int
                )
                for orbit in orbit_dict:
                    match = re.match(
                        r"(\(\(\(.*?\)\), \d+\)), Orbit \d+", orbit_mapper[orbit]
                    )
                    if match:
                        result = match.group(1)
                        try:
                            result_tuple = ast.literal_eval(result)
                        except (SyntaxError, ValueError) as e:
                            print(f"Error converting result to tuple: {e}")
                            result_tuple = None
                        if result_tuple is not None:
                            for node in orbit_dict[orbit]:
                                if hash(result_tuple) in indexed_orbit_dict:
                                    orbit_index = indexed_orbit_dict[hash(result_tuple)]
                                    node_orbit_arr[node][int(orbit_index)] += 1

                # save the output files
                print(node_orbit_arr)
                np.savetxt(
                    f"{output_dir}/node_orbit.csv",
                    node_orbit_arr,
                    delimiter=",",
                    fmt="%d",
                )
                print("INSIDEEE")
                with open(f"{output_dir}/protein_id_mapper.csv", "w") as f:
                    for protein in protein_id_dict:
                        f.write(f"{protein}, {protein_id_dict[protein]}\n")
                f.close()

                with open(f"{output_dir}/orbit_id_mapper.csv", "w") as f:
                    for orbit in indexed_orbit_dict:
                        if orbit in orbit_mapper:
                            f.write(
                                f"{orbit_mapper[orbit]} : {indexed_orbit_dict[orbit]}\n"
                            )

                with open(f"{output_dir}/orbit_hash_mapper.csv", "w") as f:
                    for orbit in indexed_orbit_dict:
                        if orbit in orbit_mapper:
                            f.write(f"{orbit} : {orbit_mapper[orbit]}\n")

                with open(f"{output_dir}/graphlet_id_mapper.csv", "w") as f:
                    for graphlet in indexed_graphlet_dict:
                        if graphlet in graphlet_mapper:
                            f.write(
                                f"{graphlet_mapper[graphlet]} : {indexed_graphlet_dict[graphlet]}\n"
                            )

                with open(f"{output_dir}/graphlet_hash_mapper.csv", "w") as f:
                    for graphlet in indexed_graphlet_dict:
                        if graphlet in graphlet_mapper:
                            f.write(f"{graphlet} : {graphlet_mapper[graphlet]}\n")

                with open(f"{output_dir}/graphlet_counts.csv", "w") as f:
                    for graphlet in indexed_graphlet_dict:
                        if graphlet in graphlet_mapper:
                            f.write(
                                f"{graphlet_mapper[graphlet]} : {three_node_graphlet_dict[graphlet]}\n"
                            )

            # if we already have the output files,
            else:
                print("post-processing only")
                with open(f"{output_dir}/graphlet_hash_mapper.csv", "r") as file:
                    csv_reader = csv.reader(file, delimiter=":")
                    for row in csv_reader:
                        graphlet_hash = int(row[0])
                        graphlet = ast.literal_eval(row[1])
                        graphlet_mapper[graphlet_hash] = graphlet
                print("post processing step 1/6")

                with open(f"{output_dir}/graphlet_counts.csv", "r") as file:
                    csv_reader = csv.reader(file, delimiter=":")
                    for row in csv_reader:
                        graphlet = ast.literal_eval(row[0])
                        count = int(row[1].strip())
                        three_node_graphlet_dict[hash(graphlet)] = count
                print("post processing step 2/6")

                with open(f"{output_dir}/orbit_hash_mapper.csv", "r") as file:
                    csv_reader = csv.reader(file, delimiter=":")
                    for row in csv_reader:
                        orbit_hash = int(row[0])
                        orbit_mapper[orbit_hash] = row[1]
                print("post processing step 3/6")

                with open(f"{output_dir}/graphlet_id_mapper.csv", "r") as file:
                    csv_reader = csv.reader(file, delimiter=":")
                    for row in csv_reader:
                        graphlet = ast.literal_eval(row[0])
                        id = row[1].strip()
                        indexed_graphlet_dict[hash(graphlet)] = id
                print("post processing step 4/6")

                with open(f"{output_dir}/orbit_id_mapper.csv", "r") as file:
                    csv_reader = csv.reader(file, delimiter=":")
                    for row in csv_reader:
                        match = re.search(r"\(\(\(.*?\)\), \d+\)", row[0])
                        if match:
                            result = ast.literal_eval(match.group())
                            id = row[1].strip()
                            # print(result, type(result), row[1].strip())
                            indexed_orbit_dict[hash(result)] = id
                            orbit_id_dict[id] = result
                print("post processing step 5/6")

                node_orbit_arr = np.loadtxt(f"{output_dir}/node_orbit.csv", delimiter=",", dtype=int)
                rows, cols = node_orbit_arr.shape

                node_orbit_arr = np.loadtxt(
                    f"{output_dir}/node_orbit.csv", delimiter=",", dtype=int
                )
                rows, cols = node_orbit_arr.shape
                # print(rows, cols)
                count = 0
                for orbit in orbit_id_dict:
                    print(orbit, orbit_id_dict[orbit])
                    for protein in range(0, rows, 1):
                        # print(orbit, type(orbit), protein, type(protein))
                        orbit_hash = hash(orbit_id_dict[orbit])
                        if orbit_hash not in orbit_dict:
                            orbit_dict[orbit_hash] = []
                        protein_count = node_orbit_arr[protein][int(orbit)]
                        # print(orbit, protein, protein_count)
                        orbit_dict[orbit_hash] += [(protein, protein_count)]
                        print(count, "/", rows * cols, end="\r")
                        count += 1

                # # Precompute orbit hashes for efficiency
                # orbit_hashes = {orbit: hash(orbit_id_dict[orbit]) for orbit in orbit_id_dict}
                # orbit_dict = {orbit_hash: [] for orbit_hash in orbit_hashes.values()}
                # count  = 0
                # # Efficiently populate orbit_dict
                # for orbit, orbit_hash in orbit_hashes.items():
                #     print(count, "/", rows * cols, end="\r")
                #     protein_counts = node_orbit_arr[:, int(orbit)]
                    
                #     # Avoid repeated appending, use list comprehension
                #     expanded_proteins = [
                #         protein - 1 for protein in range(rows) for _ in range(protein_counts[protein])
                #     ]
                    
                #     orbit_dict[orbit_hash].extend(expanded_proteins)
                #     count +=1


                # node_orbit_arr = np.loadtxt(
                #     f"{output_dir}/node_orbit.csv", delimiter=",", dtype=int
                # )
                # rows, cols = node_orbit_arr.shape
                # # print(rows, cols)
                # count = 0
                # for orbit in orbit_id_dict:
                #     print(orbit, orbit_id_dict[orbit])
                #     for protein in range(0, rows, 1):
                #         # print(orbit, type(orbit), protein, type(protein))
                #         orbit_hash = hash(orbit_id_dict[orbit])
                #         if orbit_hash not in orbit_dict:
                #             orbit_dict[orbit_hash] = []
                #         protein_count = node_orbit_arr[protein][int(orbit)]
                #         orbit_dict[orbit_hash] += [
                #             protein - 1 for _ in range(protein_count)
                #         ]
                #         print(count, "/", rows * cols, end="\r")
                #         count += 1
                print("post processing step 6/6")
                print("Completed post processing")

            print("getting stats")
            with open(f"{output_dir}/stats.csv", "w") as f:
                f.write(f"{selected_network}\n")
                f.write(f"Number of nodes: {len(G.nodes())}\n")
                f.write(f"Number of edges: {len(G.edges())}\n")
                f.write(f"Simplified Graph:\n")
                f.write(f"Number of nodes: {len(G_prime.nodes())}\n")
                f.write(f"Number of edges: {len(G_prime.edges())}\n")
                f.write(f"run time : %.3f seconds\n" % run_time)
                f.write("three node graphlet counts\n")
                count = 0
                for key in three_node_graphlet_dict:
                    f.write(
                        f"{graphlet_mapper[key]} = {three_node_graphlet_dict[key]}\n"
                    )
                    count += three_node_graphlet_dict[key]
                f.write(f"Total graphlets found: {count}\n")
                f.write(f"unique graphlet counts : {len(three_node_graphlet_dict)}\n")
                f.write(f"three node orbit counts\n")
                total_orbit_count = 0
                for orbit in orbit_dict:
                    orbit_count = 0
                    for protein_set in orbit_dict[orbit]:
                        orbit_count+=protein_set[1]
                    total_orbit_count += orbit_count
                    f.write(f"{orbit_mapper[orbit]} : {orbit_count}\n")
                f.write(f"Total orbits found: {total_orbit_count}\n")
                f.write(f"unique orbits counts : {len(orbit_dict)}\n")
            f.close()
            print()

            # plot_three_node_graphlet_distribution(
            #     three_node_graphlet_dict,
            #     graphlet_mapper,
            #     indexed_graphlet_dict,
            #     selected_network,
            #     output_dir,
            # )

            # plot_three_node_orbit_dist(
            #     orbit_dict,
            #     orbit_mapper,
            #     indexed_orbit_dict,
            #     selected_network,
            #     output_dir,)

            # print("stress protein stats")

            stress_proteins_list = get_stress_proteins(protein_id_dict, stress_data_path, '\t')
            print(stress_proteins_list)

            significance = plot_stress_orbit_distribution(orbit_dict, orbit_mapper, indexed_orbit_dict, stress_proteins_list, protein_id_dict, selected_network, output_dir, node_orbit_arr)

    # draw_labeled_multigraph(G, "label")
    # plt.show()

    # nx.draw_networkx(G_prime, with_labels=True, font_size=10)
    # plt.show()


if __name__ == "__main__":
    curses.wrapper(main)
