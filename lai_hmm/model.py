import math

STATES = ["AA", "AB", "BA", "BB"]

def _clamp(p, eps=1e-6):
    # keep away from 0 and 1 to prevent log(0)
    if p < eps:
        return eps
    if p > 1 - eps:
        return 1 - eps
    return p

def emission_logprob(state, hap1, hap2, pA, pB, eps=1e-6):
    pA = _clamp(pA, eps)
    pB = _clamp(pB, eps)

    probs = []
    for hap, pop in zip([hap1, hap2], state):
        p_alt = pA if pop == "A" else pB
        prob = p_alt if hap == "1" else (1 - p_alt)
        probs.append(prob)

    return math.log10(probs[0] * probs[1])