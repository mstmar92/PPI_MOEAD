import random
import copy


def crossover(parent1, parent2, N, child):
    pc = 0.8
    new_child = copy.deepcopy(child)

    limit = min(N, len(new_child['Chromosome']))

    for i in range(limit):
        if random.random() <= pc:
            new_child['Chromosome'][i] = parent1['Chromosome'][i]
        else:
            new_child['Chromosome'][i] = parent2['Chromosome'][i]

    return new_child
