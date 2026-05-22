import numpy as np


def same_complexes(a, b):
    if len(a) != len(b):
        return False
    for x, y in zip(a, b):
        if list(x) != list(y):
            return False
    return True


def dominate(child, near_pareto_optimal_set, nd, objective_dimension):
    if objective_dimension == 2:
        old = near_pareto_optimal_set[nd]
        if (child['Obj'][0] <= old['Obj'][0] and
            child['Obj'][1] <= old['Obj'][1] and
            (child['Obj'][0] < old['Obj'][0] or child['Obj'][1] < old['Obj'][1])):
            near_pareto_optimal_set[nd]['Dominated'] = 'T'
        else:
            near_pareto_optimal_set[nd]['Dominated'] = 'F'

    return near_pareto_optimal_set


def erase_from_near_pareto_optimal_set_overlapping(near_pareto_optimal_set):
    temp = []

    for item in near_pareto_optimal_set:
        if item.get('Dominated', '') != 'T':
            temp.append({
                'Chromosome': item['Chromosome'].copy(),
                'Dominated': item.get('Dominated', ''),
                'CmplxID': [list(c) for c in item['CmplxID']],
                'Obj': item['Obj'][:],
                'Obj1_Conductance': item['Obj1_Conductance'],
                'Obj2_TriLack': item['Obj2_TriLack']
            })

    return temp


def update_near_pareto_optimal_set_overlapping(near_pareto_optimal_set,
                                               population, pop_size,
                                               objective_dimension):
    for i in range(pop_size):
        flag = 1

        for nd in range(len(near_pareto_optimal_set)):
            if same_complexes(near_pareto_optimal_set[nd]['CmplxID'], population[i]['CommsNodes']['CommsNodes']):
                near_pareto_optimal_set[nd]['Dominated'] = ''
                flag = 0
                break

            if objective_dimension == 2:
                if (near_pareto_optimal_set[nd]['Obj'][0] <= population[i]['Obj'][0] and
                    near_pareto_optimal_set[nd]['Obj'][1] <= population[i]['Obj'][1] and
                    (near_pareto_optimal_set[nd]['Obj'][0] < population[i]['Obj'][0] or
                     near_pareto_optimal_set[nd]['Obj'][1] < population[i]['Obj'][1])):
                    near_pareto_optimal_set[nd]['Dominated'] = ''
                    flag = 0
                    break

            near_pareto_optimal_set = dominate(population[i], near_pareto_optimal_set, nd, objective_dimension)

        if flag:
            near_pareto_optimal_set = erase_from_near_pareto_optimal_set_overlapping(near_pareto_optimal_set)

            near_pareto_optimal_set.append({
                'Chromosome': population[i]['Chromosome'].copy(),
                'Dominated': '',
                'CmplxID': [list(c) for c in population[i]['CommsNodes']['CommsNodes']],
                'Obj': population[i]['Obj'][:],
                'Obj1_Conductance': population[i]['Obj'][0],
                'Obj2_TriLack': population[i]['Obj'][1]
            })

    return near_pareto_optimal_set
