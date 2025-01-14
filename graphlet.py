import sys
from matplotlib import pyplot as plt
from pathlib import Path
import networkx as nx
import numpy
import random
import scipy as sp
import csv
import curses
from itertools import combinations
import time


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


def get_two_node_graphlet_dist_adj_matrix(G, two_node_graphlet_dict):
    G_adj_matrix = get_adjacency_matrix(G)
    for i in range(len(G_adj_matrix)):
        for j in range(len(G_adj_matrix[0])):
            vector = G_adj_matrix[i][j] + G_adj_matrix[j][i]
            if hash(tuple(vector)) in two_node_graphlet_dict:
                two_node_graphlet_dict[hash(tuple(vector))] += 1
    return two_node_graphlet_dict


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

    # Example of a simple dictionary
    my_dict = {
        hash((0, 0, 0)): "0",
        hash((1, 0, 0)): "1",
        hash((0, 1, 0)): "2",
        hash((0, 0, 1)): "3",
        hash((1, 1, 0)): "4",
        hash((1, 0, 1)): "5",
        hash((0, 1, 1)): "6",
        hash((1, 1, 1)): "7",
    }

    return my_dict[ab], my_dict[ac], my_dict[ba], my_dict[bc], my_dict[ca], my_dict[cb]

def get_three_node_graphlet_dist_adj_list_v3(G: nx.MultiDiGraph, G_prime: nx.Graph):
    three_node_graphlet_dict = {}
    graphlet_mapper = {}
    orbit_dict = {}
    orbit_mapper = {}
    start_time = time.time()

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
    three_node_combination = set()  # Use a set for fast triplet lookups

    # Preprocess neighbors for each node once
    neighbors_dict = {i: set(G_prime.neighbors(i)) for i in G_prime.nodes()}
    completed_i = set()

    for i in G_prime.nodes():
        print(f"Node: {i}/{len(G_prime.nodes)}", end="\r")
        for j in neighbors_dict[i]:
            for k in neighbors_dict[j].difference(completed_i):
                if (i < k) and (i != j) and (j != k):  # Ensure no duplicates by enforcing i < k and i != j
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
                        a_b, a_c, b_a, b_c, c_a, c_b = get_three_node_graphlet_dict(ab, ac, ba, bc, ca, cb)

                        # get unsorted edge values for A, B, C
                        a_edges = [a_b, a_c]
                        b_edges = [b_a, b_c]
                        c_edges = [c_a, c_b]
                        
                        # order A, B, C edge values internally
                        a_sort = tuple(sorted(a_edges))
                        b_sort = tuple(sorted(b_edges))
                        c_sort = tuple(sorted(c_edges))

                        # Create a list of tuples in order [A, B, C]
                        tuples_list = [a_sort, b_sort, c_sort]

                        # Sort the tuples first by the first index, then by the second index
                        sorted_tuples = tuple(sorted(tuples_list, key=lambda x: (x[0], x[1])))

                        # Add the graphlet if it has not been seen yet and update the count for the graphlet
                        if hash(sorted_tuples) not in three_node_graphlet_dict:
                            three_node_graphlet_dict[hash(sorted_tuples)] = 0
                            graphlet_mapper[hash(sorted_tuples)] = sorted_tuples
                        three_node_graphlet_dict[hash(sorted_tuples)]+=1
                        print(f'Graphlet: {sorted_tuples}')
                        orbit_dict = get_orbit_per_graphlet(orbit_dict, orbit_mapper, sorted_tuples, a_edges, b_edges, c_edges, i, j, k)

        # Once we're done processing i, mark it as completed
        completed_i.add(i)
        # print(f"Completed node: {i}")

    run_time =  time.time() - start_time
    print("run time : %.3f seconds" % run_time)
    return three_node_graphlet_dict, graphlet_mapper, orbit_dict, orbit_mapper

def get_orbit_per_graphlet(orbit_dict, orbit_mapper, sorted_tuples, a_edges, b_edges, c_edges, i, j, k):

    # # 3-node line basic only ppi example 1
    # we can check which graphlet we are looking at via the sorted tuple
    if sorted_tuples == (('0', '1'), ('0', '1'), ('1', '1')):
        a_expected = ('0', '1')
        b_expected = ('0', '1')
        c_expected = ('1', '1')
        #since we arbitrarily chose what i, j, and k is, 
        # we need to check which of them match the a, b, c in the sorted tuple
        orbit_change = get_orbit_position_change(a_edges, b_edges, c_edges,
                                                 a_expected, b_expected, c_expected,
                                                 i, j, k)

        #store orbits with the key:
        # graphlet_sorted_tuple + orbit_position
        # i.e 
        # (('0', '1'), ('0', '1'), ('1', '1') ('0','0')) this specific graphlet at its 0th orbit
        # (('0', '1'), ('0', '1'), ('1', '1') ('1','1')) this specific graphlet at its 1th orbit
        zero_orbit = (0, 0)
        first_orbit = (1, 1)

        a_orbit = sorted_tuples + zero_orbit
        b_orbit = sorted_tuples + first_orbit
        # c_orbit = graphlet_sorted_tuple + second_orbit

        if hash(a_orbit) not in orbit_dict:
            orbit_dict[hash(a_orbit)] = []
            orbit_mapper[hash(a_orbit)] = "3-node line only ppi orbit 0"
        orbit_dict[hash(a_orbit)] += [orbit_change[0]]
        if hash(b_orbit) not in orbit_dict:
            orbit_dict[hash(b_orbit)] = []
            orbit_mapper[hash(b_orbit)] = "3-node line only ppi orbit 1"
        orbit_dict[hash(b_orbit)] += [orbit_change[1], orbit_change[2]]

    return orbit_dict

def get_orbit_position_change(a_edges, b_edges, c_edges, a_expected, b_expected, c_expected, i, j, k):
    # a, b, c
    if a_edges == a_expected and b_edges == b_expected and c_edges == c_expected:
        return [i,j,k]
    # a, c, b
    if a_edges == a_expected and b_edges == c_expected and c_edges == b_expected:
        return [i,k,j]
    # b, a, c
    if a_edges == b_expected and b_edges == a_expected and c_edges == c_expected:
        return [j,i,k]
    # b, c, a
    if a_edges == b_expected and b_edges == c_expected and c_edges == a_expected:
        return [j,k,i]
    # c, a, b
    if a_edges == c_expected and b_edges == a_expected and c_edges == b_expected:
        return [k,i,j]
    # c, b, a
    if a_edges == c_expected and b_edges == b_expected and c_edges == a_expected:
        return [k,j,i]
    raise ValueError(
        f"Unexpected edges configuration: a_edges={a_edges}, b_edges={b_edges}, c_edges={c_edges}"
    )    

def get_three_node_graphlet_dist_adj_list_v4(G: nx.MultiDiGraph, G_prime: nx.Graph):
    three_node_graphlet_dict = {}
    graphlet_mapper = {}
    start_time = time.time()

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
    three_node_combination = set()  # Use a set for fast triplet lookups

    # Preprocess neighbors for each node once
    neighbors_dict = {i: set(G_prime.neighbors(i)) for i in G_prime.nodes()}
    completed_i = set()

    for i in G_prime.nodes():
        print(f"Node: {i}", end="\r")
        for j in neighbors_dict[i]:
            for k in neighbors_dict[j].difference(completed_i):
                if (i < k) and (i != j) and (j != k):  # Ensure no duplicates by enforcing i < k and i != j
                    triplet = tuple(sorted([i, j, k]))
                    if triplet not in three_node_combination:
                        three_node_combination.add(triplet)
                        # print(f"Triplet: {i}, {j}, {k}")

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
                        a_b, a_c, b_a, b_c, c_a, c_b = get_three_node_graphlet_dict(ab, ac, ba, bc, ca, cb)

                        # get unsorted edge values for A, B, C
                        a_unsort = [a_b, a_c]
                        b_unsort = [b_a, b_c]
                        c_unsort = [c_a, c_b]
                        
                        # order A, B, C edge values internally
                        a_sort = tuple(sorted(a_unsort))
                        b_sort = tuple(sorted(b_unsort))
                        c_sort = tuple(sorted(c_unsort))

                        # Create a list of tuples in order [A, B, C]
                        tuples_list = [a_sort, b_sort, c_sort]

                        unsorted_tuples = tuple(tuples_list)

                        # Sort the tuples first by the first index, then by the second index
                        sorted_tuples = tuple(sorted(tuples_list, key=lambda x: (x[0], x[1])))

                        # Add the graphlet if it has not been seen yet and update the count for the graphlet
                        if hash(sorted_tuples) not in three_node_graphlet_dict:
                            three_node_graphlet_dict[hash(sorted_tuples)] = 0
                            graphlet_mapper[hash(sorted_tuples)] = sorted_tuples
                        three_node_graphlet_dict[hash(sorted_tuples)]+=1

        # Once we're done processing i, mark it as completed
        completed_i.add(i)
        # print(f"Completed node: {i}")

    run_time =  time.time() - start_time
    print("run time : %.3f seconds" % run_time)
    return three_node_graphlet_dict, graphlet_mapper

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


def get_user_inputs(selected_network, selected_graphlet):
    ppi_path = None
    reg_path = None
    graphlet_option = None

    match selected_network:
        case "D. melanogaster":
            ppi_path = Path("data/fly_ppi.csv")
            reg_path = Path("data/fly_reg.csv")
        case "B. subtilis":
            ppi_path = Path("data/bsub_ppi.csv")
            reg_path = Path("data/bsub_reg.csv")
        case "S. cerevisiae":
            ppi_path = Path("data/ceravisiae_ppi.csv")
            reg_path = Path("data/ceravisiae_reg.csv")
        case "D. rerio":
            ppi_path = Path("data/drerio_ppi.csv")
            reg_path = Path("data/drerio_reg.csv")
        case "C. elegans":
            ppi_path = Path("data/elegans_ppi.csv")
            reg_path = Path("data/elegans_reg.csv")
        case "Test network":
            ppi_path = Path("data/test_ppi.csv")
            reg_path = Path("data/test_reg.csv")
        case "Toy network":
            ppi_path = Path("data/toy_ppi.csv")
            reg_path = Path("data/toy_reg.csv")
        case "Shuffled toy network":
            ppi_path = Path("data/toy2_ppi.csv")
            reg_path = Path("data/toy2_reg.csv")

    if selected_graphlet == "2-node":
        graphlet_option = 2
    elif selected_graphlet == "3-node":
        graphlet_option = 3

    return ppi_path, reg_path, graphlet_option


def main(stdscr):
    try:
        network = [
            "D. melanogaster",
            "B. subtilis",
            "S. cerevisiae",
            "D. rerio",
            "C. elegans",
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

        node_limit = type_answer(stdscr, 0)
        stdscr.clear()

        edge_limit = type_answer(stdscr, 1)
        stdscr.clear()
    except Exception as e:
        stdscr.addstr(0, 0, f"An error occurred: {str(e)}")
        stdscr.getch()
    finally:
        curses.endwin()

    ppi_path, reg_path, graphlet_mode = get_user_inputs(
        selected_network, selected_graphlet
    )

    two_node_graphlet_dict, two_node_graphlet_labels = get_two_node_dict()
    protein_id_dict = get_protein_id_dict(ppi_path, reg_path)

    for protein in protein_id_dict:
        print(f"Protein ID dictionary: {protein} {protein_id_dict[protein]}")
    
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

    if graphlet_mode == 2:
        two_node_graphlet_dict, two_node_orbit_dict = get_two_node_graphlet_stats(
            G, two_node_graphlet_dict
        )
        print("\ntwo node graphlet counts")
        for key in two_node_graphlet_labels:
            print(f"G_{two_node_graphlet_labels[key]} = {two_node_graphlet_dict[key]}")
        print("\ntwo node graphlet orbit counts")
        for key in two_node_orbit_dict:
            print(f"{key} = {len(two_node_orbit_dict[key])}")
    elif graphlet_mode == 3:
        # three_node_graphlet_dict = get_three_node_graphlet_dist_adj_list(G)
        # three_node_graphlet_dict, graphlet_mapper = get_three_node_graphlet_dist_adj_list_v4(G, G_prime)
        three_node_graphlet_dict, graphlet_mapper, orbit_dict, orbit_mapper = get_three_node_graphlet_dist_adj_list_v3(G, G_prime)        # three_node_graphlet_dict, graphlet_mapper = get_three_node_graphlet_dist_adj_list_v4(G, G_prime)
        print("\nthree node graphlet counts")
        count = 0 
        for key in three_node_graphlet_dict:
            # print(f"{key} = {three_node_graphlet_dict[key]}")
            print(f"{graphlet_mapper[key]} = {three_node_graphlet_dict[key]}")
            count += three_node_graphlet_dict[key]
        print(f"Total graphlets found: {count}")
        print(f"\n unique graphlet counts : {len(three_node_graphlet_dict)}")
        print(f"\n three node orbits")
        for orbit in orbit_dict:
            print(f"{orbit_mapper[orbit]} : {orbit_dict[orbit]}")

    # draw_labeled_multigraph(G, "label")
    # plt.show()

    # nx.draw_networkx(G_prime, with_labels=True, font_size=10)
    # plt.show()


if __name__ == "__main__":
    curses.wrapper(main)