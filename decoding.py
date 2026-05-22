import numpy as np


def compute_cmplx_decoding_overlapping(individual, protein_label_pairs):
    chromosome = individual['Chromosome']
    num_edges = len(chromosome)

    adj_list = [[] for _ in range(num_edges)]

    for i in range(num_edges):
        j = int(chromosome[i])
        if j > 0 and j <= num_edges:
            adj_list[i].append(j - 1)
            adj_list[j - 1].append(i)

    visited = [False] * num_edges
    comm_edges = []

    for i in range(num_edges):
        if not visited[i]:
            comm = []
            stack = [i]
            visited[i] = True

            while len(stack) > 0:
                current = stack.pop()
                comm.append(current + 1)

                for neighbor in adj_list[current]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        stack.append(neighbor)

            comm_edges.append(sorted(comm))

    comms_nodes = []
    for comm_e in comm_edges:
        comm_n = []
        for edge_id in comm_e:
            node_source = int(protein_label_pairs[edge_id - 1, 0])
            node_target = int(protein_label_pairs[edge_id - 1, 1])
            comm_n.extend([node_source, node_target])

        comms_nodes.append(sorted(list(set(comm_n))))

    individual['CommsNodes'] = {'CommsNodes': comms_nodes}
    return individual


def individual2cmplx_decoding_overlapping(protein_label_pairs, population, population_size):
    for i in range(population_size):
        population[i] = compute_cmplx_decoding_overlapping(population[i], protein_label_pairs)
    return population
