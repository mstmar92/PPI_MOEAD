import numpy as np
from math import comb


def compute_fitness_overlapping(A, N, m,
                                indices_interaction_protein, num_interaction_protein, max_num_interaction_protein,
                                population, population_size):
    for i in range(population_size):
        population[i]['Obj'] = []
        population[i]['Obj1_Conductance'] = 0
        population[i]['Obj2_TriLack'] = 0
        population[i]['Kv'] = 0
        population[i]['avgTri'] = 0

    for i in range(population_size):
        population[i] = compute_our_fitness_collection_overlapping(
            A, N, indices_interaction_protein, num_interaction_protein, population[i]
        )

    return population


def compute_our_fitness_collection_overlapping(A, N,
                                               indices_interaction_protein, num_interaction_protein,
                                               individual):
    min_size = 3
    eps_val = 1e-12

    A = np.array(A)
    A = A.copy()
    np.fill_diagonal(A, 0)

    deg = np.sum(A, axis=1)
    total_vol = np.sum(deg)

    cmplx_id = individual['CommsNodes']['CommsNodes']

    cond_sum = 0
    tri_sum = 0
    kv = 0

    for members in cmplx_id:
        members = sorted(list(set([int(x) for x in members if int(x) >= 1 and int(x) <= N])))
        s = len(members)

        if s < min_size:
            continue

        kv += 1
        idx = [x - 1 for x in members]
        A_sub = A[np.ix_(idx, idx)]

        m_in = np.sum(np.triu(A_sub, 1) != 0)
        vol_c = np.sum(deg[idx])
        cut_c = vol_c - 2 * m_in

        denom = min(vol_c, total_vol - vol_c) + eps_val
        phi = cut_c / denom
        cond_sum += phi

        triangles = np.trace(A_sub @ A_sub @ A_sub) / 6
        tri_den = triangles / (comb(s, 3) + eps_val)
        tri_sum += tri_den

    if kv == 0:
        obj1 = 1
        obj2 = 1
        avg_tri = 0
    else:
        obj1 = cond_sum / kv
        avg_tri = tri_sum / kv
        obj2 = 1 - avg_tri

    individual['Obj'] = [float(obj1), float(obj2)]
    individual['Obj1_Conductance'] = float(obj1)
    individual['Obj2_TriLack'] = float(obj2)
    individual['Kv'] = kv
    individual['avgTri'] = float(avg_tri)

    return individual
