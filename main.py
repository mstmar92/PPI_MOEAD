import os
import pickle
import numpy as np
from scipy.io import loadmat

from population import create_population_overlapping
from decoding import individual2cmplx_decoding_overlapping
from fitness import compute_fitness_overlapping
from moead_overlapping import moead_overlapping


def get_value(data, name, default=None):
    if name in data:
        return data[name]
    return default


def load_dataset(network_number):
    dataset_name = ''


    if network_number == 1:
        ppi_file = 'DataSets/PPI/1-Protein-Yeast-D1-Files.mat'
        complex_file = 'DataSets/Complex/Complex-D1-Files.mat'
    elif network_number == 2:
        ppi_file = 'DataSets/PPI/2-Protein-Yeast-D2-Files.mat'
        complex_file = 'DataSets/Complex/Complex-D2-Files.mat'
    elif network_number == 4:
        dataset_name = 'CollinsCYC'
        ppi_file = 'DataSets/PPI/4-Collins-CYC2008.mat'
        complex_file = 'DataSets/Complex/4-CYC2008.mat'
    elif network_number == 5:
        ppi_file = 'DataSets/PPI/5-Collins-MIPS.mat'
        complex_file = 'DataSets/Complex/5-MIPS.mat'
    else:
        raise ValueError('Wrong network number')

    ppi_data = loadmat(ppi_file)

    data = {}
    data.update(ppi_data)

    if os.path.exists(complex_file):
        complex_data = loadmat(complex_file)
        data.update(complex_data)
    else:
        data['NumberOfProteinsInComplexes'] = None
        data['NumberOfKnownProteinsInComplexes'] = None

    data['DatasetName'] = dataset_name

    return data

def main():
    network_number = int(input('Enter PPINs Network: '))

    data = load_dataset(network_number)

    A = get_value(data, 'A')
    N = int(np.squeeze(get_value(data, 'N')))
    protein_label_pairs = get_value(data, 'ProteinLabelPairs')

    indices_interaction_protein = get_value(data, 'IndicesInteractionProtein')
    num_interaction_protein = np.squeeze(get_value(data, 'NumInteractionProtein')).astype(int)
    max_num_interaction_protein = get_value(data, 'MaxNumInteractionProtein')

    number_of_proteins_in_complexes = get_value(data, 'NumberOfProteinsInComplexes')
    number_of_known_proteins_in_complexes = get_value(data, 'NumberOfKnownProteinsInComplexes')
    known_proteins = get_value(data, 'KnownProteins')

    m = np.sum(A) / 2

    max_run = 10
    pm = 0.2
    heuristic = 1
    population_size = 100

    results_group = []
    population_group = []

    for run_number in range(max_run):
        population = create_population_overlapping(protein_label_pairs, population_size)
        population_group.append({'Population': population})

    for run_number in range(max_run):
        print('RunNumber:', run_number + 1)

        population = population_group[run_number]['Population']

        population = individual2cmplx_decoding_overlapping(protein_label_pairs,
                                                           population, population_size)

        population = compute_fitness_overlapping(A, N, m,
                                                 indices_interaction_protein,
                                                 num_interaction_protein,
                                                 max_num_interaction_protein,
                                                 population, population_size)

        results = moead_overlapping(A, N, m,
                                    indices_interaction_protein,
                                    num_interaction_protein,
                                    max_num_interaction_protein,
                                    number_of_proteins_in_complexes,
                                    number_of_known_proteins_in_complexes,
                                    population,
                                    population_size,
                                    heuristic,
                                    pm,
                                    known_proteins,
                                    protein_label_pairs)

        results_group.append({'Results': results})
        population_group[run_number]['Population'] = population

    os.makedirs('RepOverlap/MOEA', exist_ok=True)
    out_name = f'RepOverlap/MOEA/MOEAD_PPI_{network_number}_Heuristic_{heuristic}_Pm_{pm}NumOfRuns{max_run}.pkl'

    with open(out_name, 'wb') as f:
        pickle.dump({'ResultsGroup': results_group}, f)

    print('Saved:', out_name)


if __name__ == '__main__':
    main()
