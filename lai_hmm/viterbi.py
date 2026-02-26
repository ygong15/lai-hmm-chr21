from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Mapping, Sequence, Tuple

from .model import STATES

LogProb = float
State = str


def _safe_log10(p: float, floor: float = 1e-300) -> float:
    """
    Just log10 but we clamp really small values so we never hit log(0).
    """
    return math.log10(max(p, floor))


def _hamming_dist(a: str, b: str) -> int:
    """
    Count how many haplotypes switch ancestry between two states.
    Since states are length-2 strings (like "AB"), this is at most 2.
    """
    return sum(x != y for x, y in zip(a, b))


def _step_recomb_prob(delta_bp: int, rho: float) -> float:
    """
    Convert physical distance into recombination probability.

    We use:
        r = 1 - exp(-rho * delta)

    So larger distance means higher chance of switching.
    """
    return 1.0 - math.exp(-rho * delta_bp)


def trans_logprob(prev_state: State, next_state: State, r: float) -> LogProb:
    """
    Transition probability between two diploid ancestry states.

    We treat the two haplotypes independently.
    So:
        0 switches → (1-r)^2
        1 switch   → r(1-r)
        2 switches → r^2
    """
    dist = _hamming_dist(prev_state, next_state)

    if dist == 0:
        p = (1.0 - r) ** 2
    elif dist == 1:
        p = r * (1.0 - r)
    else:
        p = r ** 2

    return _safe_log10(p)


@dataclass
class _DPCell:
    """
    So we can return (score, prev_state) together.
    """
    score: LogProb
    prev: State


def _init_scores(first_emission: Mapping[State, LogProb]) -> Dict[State, LogProb]:
    """
    Initialize DP with uniform prior over states.
    """
    init = _safe_log10(1.0 / len(STATES))
    return {s: init + first_emission[s] for s in STATES}


def _best_prev_state(
    prev_scores: Mapping[State, LogProb],
    curr_state: State,
    curr_emission: LogProb,
    r: float,
) -> _DPCell:
    """
    For a fixed current state, try all previous states
    and pick the one that gives the best score.
    """
    best_prob = None
    best_prev = None

    for sp in STATES:
        prob = prev_scores[sp] + trans_logprob(sp, curr_state, r) + curr_emission
        if best_prob is None or prob > best_prob:
            best_prob = prob
            best_prev = sp

    return _DPCell(score=best_prob, prev=best_prev)  # safe since STATES non-empty


def viterbi(
    emissions: Sequence[Mapping[State, LogProb]],
    positions_bp: Sequence[int],
    rho: float,
) -> List[State]:
    """
    Standard Viterbi decoding.

    positions_bp lets us compute distance dependent recombination.
    """
    n = len(emissions)
    if n == 0:
        return []

    if len(positions_bp) != n:
        raise ValueError("positions and emissions must have same length")

    # back[i][s] = best previous state leading to s at position i
    back: List[Dict[State, State]] = [{} for _ in range(n)]

    # Initialize with uniform prior
    prev_scores = _init_scores(emissions[0])

    # Forward DP
    for i in range(1, n):
        delta = positions_bp[i] - positions_bp[i - 1]
        r = _step_recomb_prob(delta, rho)

        curr_scores: Dict[State, LogProb] = {}

        for s in STATES:
            cell = _best_prev_state(prev_scores, s, emissions[i][s], r)
            curr_scores[s] = cell.score
            back[i][s] = cell.prev

        prev_scores = curr_scores

    # Pick best ending state
    last_state = max(prev_scores, key=prev_scores.get)

    # Backtrack
    path: List[State] = [last_state]
    for i in range(n - 1, 0, -1):
        last_state = back[i][last_state]
        path.append(last_state)

    path.reverse()
    return path


# Extra helper
def viterbi_score_only(
    emissions: Sequence[Mapping[State, LogProb]],
    positions_bp: Sequence[int],
    rho: float,
) -> Tuple[List[State], LogProb]:
    """
    Return both path and final log-score.
    """
    path = viterbi(emissions, positions_bp, rho)
    if not path:
        return path, float("-inf")

    score = _safe_log10(1.0 / len(STATES)) + emissions[0][path[0]]

    for i in range(1, len(path)):
        r = _step_recomb_prob(positions_bp[i] - positions_bp[i - 1], rho)
        score += trans_logprob(path[i - 1], path[i], r) + emissions[i][path[i]]

    return path, score