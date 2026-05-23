import random
import numpy as np


def create_population_overlapping(protein_pairs, population_size):
    num_edges = protein_pairs.shape[0]

    adjacent_edges = [[] for _ in range(num_edges)]

    for i in range(num_edges):
        node_a = protein_pairs[i, 0]
        node_b = protein_pairs[i, 1]

        for j in range(num_edges):
            if i == j:
                continue

            if protein_pairs[j, 0] in [node_a, node_b] or protein_pairs[j, 1] in [node_a, node_b]:
                adjacent_edges[i].append(j + 1)  

    population = []
    for _ in range(population_size):
        chromosome = np.zeros(num_edges, dtype=int)
        for edge_counter in range(num_edges):
            if len(adjacent_edges[edge_counter]) > 0:
                chromosome[edge_counter] = random.choice(adjacent_edges[edge_counter])

        population.append({
            'Chromosome': chromosome,
            'CommsNodes': {'CommsNodes': []},
            'Obj': []
        })

    return population
