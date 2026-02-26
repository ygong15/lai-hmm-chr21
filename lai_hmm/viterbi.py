import math
from .model import STATES

def trans_logprob(prev_state, next_state, r):
    # r is per-step recombination probability per haplotype
    dist = sum(a != b for a, b in zip(prev_state, next_state))
    if dist == 0:
        p = (1 - r) ** 2
    elif dist == 1:
        p = r * (1 - r)
    else:
        p = r ** 2
    return math.log10(max(p, 1e-300))

def viterbi(emissions, positions_bp, rho):
    V = []
    back = []

    init_prob = math.log10(1 / len(STATES))
    V.append({s: init_prob + emissions[0][s] for s in STATES})
    back.append({})

    for i in range(1, len(emissions)):
        V.append({})
        back.append({})

        delta = positions_bp[i] - positions_bp[i-1]
        r = 1 - math.exp(-rho * delta)   # distance-aware

        for s in STATES:
            best_prob = None
            best_prev = None
            for sp in STATES:
                prob = V[i-1][sp] + trans_logprob(sp, s, r) + emissions[i][s]
                if best_prob is None or prob > best_prob:
                    best_prob = prob
                    best_prev = sp
            V[i][s] = best_prob
            back[i][s] = best_prev

    last_state = max(V[-1], key=V[-1].get)
    path = [last_state]
    for i in reversed(range(1, len(emissions))):
        last_state = back[i][last_state]
        path.append(last_state)
    path.reverse()
    return path