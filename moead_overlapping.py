import random
import copy

from moead_utils import init_weights, update_reference_overlapping, update_problem_overlapping
from pareto import update_near_pareto_optimal_set_overlapping
from operators import crossover
from mutation import mutation_sor_overlapping
from decoding import individual2cmplx_decoding_overlapping
from fitness import compute_fitness_overlapping


def moead_overlapping(A, N, m,
                      indices_interaction_protein, num_interaction_protein, max_num_interaction_protein,
                      number_of_proteins_in_complexes, number_of_known_proteins_in_complexes,
                      population, population_size,
                      heuristic, pm,
                      known_proteins, protein_label_pairs):

    objective_dimension = 2

    ideal_point = [10000000, 10000000]
    indiv_point = [{}, {}]

    params = {}
    params['PopSize'] = population_size
    params['Neighbour'] = 5
    params['Generations'] = 100
    params['DMethod'] = 'ts'

    sub_problem = init_weights(params['PopSize'], params['Neighbour'], objective_dimension)
    params['PopSize'] = len(sub_problem)

    params['Mut'] = {}
    params['Mut']['MinSize'] = 3
    params['Mut']['pAddTriangle'] = 0.6
    params['Mut']['pPrune'] = 0.4
    params['Mut']['MaxComplexAddsPerNode'] = 2
    params['Mut']['MaxMovesPerNode'] = 1

    ideal_point, indiv_point = update_reference_overlapping(ideal_point, indiv_point,
                                                            population, params['PopSize'])

    results = []
    near_pareto_optimal_set = []

    near_pareto_optimal_set = update_near_pareto_optimal_set_overlapping(
        near_pareto_optimal_set, population, params['PopSize'], objective_dimension
    )

    results.append({'NearParetoOptimalSet': near_pareto_optimal_set})

    for generation_counter in range(2, params['Generations'] + 1):
        print('Generation:', generation_counter)
        print('IdealPoint:', ideal_point)

        child = copy.deepcopy(population)
        parent1 = [None] * params['PopSize']
        parent2 = [None] * params['PopSize']

        for problem_counter in range(params['PopSize']):
            rand1 = random.randint(0, params['Neighbour'] - 2)
            rand2 = random.randint(0, params['Neighbour'] - 2)

            p1 = sub_problem[problem_counter]['Neighbour'][rand1]
            p2 = sub_problem[problem_counter]['Neighbour'][rand2]

            parent1[problem_counter] = population[p1]
            parent2[problem_counter] = population[p2]

        for problem_counter in range(params['PopSize']):
            child[problem_counter] = crossover(parent1[problem_counter], parent2[problem_counter], N,
                                               child[problem_counter])

        for problem_counter in range(params['PopSize']):
            child[problem_counter] = mutation_sor_overlapping(child[problem_counter], A,
                                                              indices_interaction_protein,
                                                              num_interaction_protein,
                                                              pm, params['Mut'])

        child = individual2cmplx_decoding_overlapping(protein_label_pairs, child, population_size)
        child = compute_fitness_overlapping(A, N, m,
                                            indices_interaction_protein, num_interaction_protein,
                                            max_num_interaction_protein,
                                            child, population_size)

        ideal_point, indiv_point = update_reference_overlapping(ideal_point, indiv_point,
                                                                child, params['PopSize'])

        population = update_problem_overlapping(population, child,
                                                sub_problem, params,
                                                indiv_point, objective_dimension)

        near_pareto_optimal_set = update_near_pareto_optimal_set_overlapping(
            near_pareto_optimal_set, child, params['PopSize'], objective_dimension
        )

        results.append({'NearParetoOptimalSet': near_pareto_optimal_set})

    return results
