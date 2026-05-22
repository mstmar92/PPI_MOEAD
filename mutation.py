import random
import numpy as np


def mutation_sor_overlapping(child, A,
                             indices_interaction_protein, num_interaction_protein,
                             pm, params=None):
    if params is None:
        params = {}

    min_size = params.get('MinSize', 3)
    p_add_triangle = params.get('pAddTriangle', 0.6)
    p_prune = params.get('pPrune', 0.4)
    max_complex_adds_per_node = params.get('MaxComplexAddsPerNode', 2)
    max_moves_per_node = params.get('MaxMovesPerNode', 1)

    A = np.array(A)
    N = A.shape[0]

    comms = [list(set(c)) for c in child['CommsNodes']['CommsNodes']]
    K = len(comms)

    protein_complex_map = [[] for _ in range(N + 1)]
    for c in range(K):
        members = sorted([int(x) for x in comms[c] if int(x) >= 1 and int(x) <= N])
        comms[c] = members
        for u in members:
            protein_complex_map[u].append(c)

    cmask = np.zeros((N + 1, K), dtype=bool)
    for c in range(K):
        for u in comms[c]:
            cmask[u, c] = True

    for u in range(1, N + 1):
        if num_interaction_protein[u - 1] <= 0 or random.random() > pm:
            continue

        curr_cs = protein_complex_map[u][:]
        if len(curr_cs) == 0:
            continue

        move_count = 0
        for c in curr_cs:
            if move_count >= max_moves_per_node:
                break

            kin = 0
            kout = 0

            for t in range(int(num_interaction_protein[u - 1])):
                v = int(indices_interaction_protein[u - 1, t])
                if v < 1 or v > N:
                    continue

                if cmask[v, c]:
                    kin += A[u - 1, v - 1]
                else:
                    kout += A[u - 1, v - 1]

            if kin <= kout:
                best_c = None
                best_score = -float('inf')

                cand = []
                for t in range(int(num_interaction_protein[u - 1])):
                    v = int(indices_interaction_protein[u - 1, t])
                    if v >= 1 and v <= N:
                        cand.extend(protein_complex_map[v])

                cand = sorted(list(set(cand)))
                cand = [x for x in cand if x != c]

                for cc in cand:
                    if cmask[u, cc]:
                        continue

                    kin2 = 0
                    kout2 = 0
                    for t in range(int(num_interaction_protein[u - 1])):
                        v = int(indices_interaction_protein[u - 1, t])
                        if v < 1 or v > N:
                            continue

                        if cmask[v, cc]:
                            kin2 += A[u - 1, v - 1]
                        else:
                            kout2 += A[u - 1, v - 1]

                    tri_gain = 0
                    for t in range(int(num_interaction_protein[u - 1])):
                        v = int(indices_interaction_protein[u - 1, t])
                        if v < 1 or v > N or not cmask[v, cc] or A[u - 1, v - 1] == 0:
                            continue

                        vn = indices_interaction_protein[v - 1, :int(num_interaction_protein[v - 1])]
                        vn = [int(x) for x in vn if int(x) >= 1 and int(x) <= N]
                        tri_gain += sum(cmask[x, cc] for x in vn)

                    score = (kin2 - kout2) + 0.1 * tri_gain

                    if score > best_score:
                        best_score = score
                        best_c = cc

                if best_c is not None:
                    comms[best_c].append(u)
                    cmask[u, best_c] = True
                    protein_complex_map[u].append(best_c)
                    move_count += 1

        if random.random() <= p_add_triangle:
            adds = 0
            cand = []
            for t in range(int(num_interaction_protein[u - 1])):
                v = int(indices_interaction_protein[u - 1, t])
                if v >= 1 and v <= N:
                    cand.extend(protein_complex_map[v])

            cand = sorted(list(set(cand)))

            for cc in cand:
                if adds >= max_complex_adds_per_node:
                    break

                if cmask[u, cc]:
                    continue

                nin = 0
                for t in range(int(num_interaction_protein[u - 1])):
                    v = int(indices_interaction_protein[u - 1, t])
                    if v >= 1 and v <= N and cmask[v, cc] and A[u - 1, v - 1] == 1:
                        nin += 1
                        if nin >= 2:
                            break

                if nin >= 2:
                    comms[cc].append(u)
                    cmask[u, cc] = True
                    protein_complex_map[u].append(cc)
                    adds += 1

        if random.random() <= p_prune:
            curr_cs2 = sorted(list(set(protein_complex_map[u])))

            for c in curr_cs2:
                members = sorted(list(set(comms[c])))
                if len(members) <= min_size:
                    continue

                kin = 0
                kout = 0
                tri_support = 0

                for t in range(int(num_interaction_protein[u - 1])):
                    v = int(indices_interaction_protein[u - 1, t])
                    if v < 1 or v > N:
                        continue

                    if cmask[v, c]:
                        kin += A[u - 1, v - 1]

                        vn = indices_interaction_protein[v - 1, :int(num_interaction_protein[v - 1])]
                        vn = [int(x) for x in vn if int(x) >= 1 and int(x) <= N]
                        tri_support += sum(cmask[x, c] for x in vn)
                    else:
                        kout += A[u - 1, v - 1]

                if kin <= kout and tri_support <= 1:
                    comms[c] = [x for x in members if x != u]
                    cmask[u, c] = False
                    protein_complex_map[u] = [x for x in protein_complex_map[u] if x != c]

    for c in range(K):
        comms[c] = sorted(list(set(comms[c])))

    child['CommsNodes']['CommsNodes'] = comms
    return child
