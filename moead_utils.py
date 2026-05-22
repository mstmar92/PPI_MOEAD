import numpy as np


def init_weights(pop_size, neighbour, objective_dimension):
    sub_problem = []

    for i in range(1, pop_size + 1):
        if objective_dimension == 2:
            weight = np.array([i / pop_size, (pop_size - i) / pop_size], dtype=float)
            sub_problem.append({
                'Weight': weight,
                'Neighbour': [],
                'Optimal': float('inf'),
                'OptPoint': [],
                'Curpoint': []
            })
        elif objective_dimension == 3:
            for j in range(1, pop_size + 1):
                if i + j <= pop_size:
                    k = pop_size - i - j
                    weight = np.array([i / pop_size, j / pop_size, k / pop_size], dtype=float)
                    sub_problem.append({
                        'Weight': weight,
                        'Neighbour': [],
                        'Optimal': float('inf'),
                        'OptPoint': [],
                        'Curpoint': []
                    })

    length = len(sub_problem)
    distance_matrix = np.zeros((length, length))

    for i in range(length):
        for j in range(i + 1, length):
            a = sub_problem[i]['Weight']
            b = sub_problem[j]['Weight']
            distance_matrix[i, j] = np.sum((a - b) * (a - b))
            distance_matrix[j, i] = distance_matrix[i, j]

        sindex = np.argsort(distance_matrix[i, :])
        sub_problem[i]['Neighbour'] = sindex[:neighbour].tolist()

    return sub_problem


def update_reference_overlapping(ideal_point, indiv_point, population, population_size):
    for i in range(population_size):
        if population[i]['Obj'][0] < ideal_point[0]:
            ideal_point[0] = population[i]['Obj'][0]
            indiv_point[0] = {
                'Chromosome': population[i]['Chromosome'].copy(),
                'CmplxID': population[i]['CommsNodes']['CommsNodes'],
                'Ideal': ideal_point[0],
                'Obj': population[i]['Obj'][:]
            }

        if population[i]['Obj'][1] < ideal_point[1]:
            ideal_point[1] = population[i]['Obj'][1]
            indiv_point[1] = {
                'Chromosome': population[i]['Chromosome'].copy(),
                'CmplxID': population[i]['CommsNodes']['CommsNodes'],
                'Ideal': ideal_point[1],
                'Obj': population[i]['Obj'][:]
            }

    return ideal_point, indiv_point


def scalar_function(individual, weight, indiv_point, objective_dimension):
    max_function = -1.0e30

    for d in range(objective_dimension):
        z_d = indiv_point[d]['Ideal']
        f_d = individual['Obj'][d]

        diff = abs(f_d - z_d)
        val = weight[d] * diff

        if val > max_function:
            max_function = val

    return max_function


def update_problem_overlapping(individual, child, sub_problems, params, indiv_point, objective_dimension):
    for problem_counter in range(params['PopSize']):
        for neighbour_counter in range(params['Neighbour']):
            neighbour_index = sub_problems[problem_counter]['Neighbour'][neighbour_counter]
            weight_of_neighbour = sub_problems[neighbour_index]['Weight']

            f1 = scalar_function(individual[neighbour_index], weight_of_neighbour, indiv_point, objective_dimension)
            f2 = scalar_function(child[problem_counter], weight_of_neighbour, indiv_point, objective_dimension)

            if f2 < f1:
                individual[neighbour_index]['Chromosome'] = child[problem_counter]['Chromosome'].copy()
                individual[neighbour_index]['CommsNodes'] = {
                    'CommsNodes': [list(c) for c in child[problem_counter]['CommsNodes']['CommsNodes']]
                }
                individual[neighbour_index]['Obj'] = child[problem_counter]['Obj'][:]
                individual[neighbour_index]['Obj1_Conductance'] = child[problem_counter]['Obj'][0]
                individual[neighbour_index]['Obj2_TriLack'] = child[problem_counter]['Obj'][1]

    return individual
