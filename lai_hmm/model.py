from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, Literal, Mapping, Sequence, Tuple

# Hidden ancestry states: (hap1 ancestry, hap2 ancestry)
# For example
# "AB" means hap1 from pop A, hap2 from pop B
STATES: Sequence[str] = ("AA", "AB", "BA", "BB")

Hap = Literal["0", "1"]
State = Literal["AA", "AB", "BA", "BB"]


@dataclass(frozen=True)
class EmissionParams:
    """
    Holds allele frequencies for a single variant site.

    pA_alt / pB_alt are ALT allele frequencies in the two populations.
    eps is just used to avoid log(0) when frequencies are extreme.
    """
    pA_alt: float
    pB_alt: float
    eps: float = 1e-6


def is_valid_state(s: str) -> bool:
    """Small helper so we don't accidentally use invalid states."""
    return s in STATES


def is_valid_hap(h: str) -> bool:
    """Haplotype alleles should just be '0' or '1'."""
    return h == "0" or h == "1"


def _clamp_prob(p: float, eps: float) -> float:
    """
    Clamp probability away from 0 and 1.
    This keeps things nice when we take logs.
    """
    if p < eps:
        return eps
    if p > 1.0 - eps:
        return 1.0 - eps
    return p


def _alt_freq_for_pop(pop: str, pA: float, pB: float) -> float:
    """Pick ALT allele frequency depending on ancestry label."""
    return pA if pop == "A" else pB


def _allele_prob(hap: Hap, p_alt: float) -> float:
    """
    Probability of seeing hap allele (0/1)
    given ALT allele frequency.
    """
    # hap == "1" means ALT allele
    # hap == "0" means REF allele
    return p_alt if hap == "1" else (1.0 - p_alt)


def emission_prob(
    state: str,
    hap1: str,
    hap2: str,
    pA: float,
    pB: float,
    eps: float = 1e-6,
) -> float:
    """
    Non-log emission probability for a phased genotype.

    We assume the two haplotypes are independent
    given the hidden ancestry state.
    """
    if not is_valid_state(state):
        raise ValueError(f"Invalid state: {state!r}")
    if not (is_valid_hap(hap1) and is_valid_hap(hap2)):
        raise ValueError(f"Invalid haplotypes: {hap1!r}, {hap2!r}")
    if not (0.0 < eps < 0.5):
        raise ValueError(f"eps must be between 0 and 0.5")

    # Clamp allele frequencies to avoid log(0)
    pA_c = _clamp_prob(pA, eps)
    pB_c = _clamp_prob(pB, eps)

    # Pick the right ALT frequency for each haplotype
    p_alt_1 = _alt_freq_for_pop(state[0], pA_c, pB_c)
    p_alt_2 = _alt_freq_for_pop(state[1], pA_c, pB_c)

    # Multiply per-haplotype probabilities
    prob1 = _allele_prob(hap1, p_alt_1)
    prob2 = _allele_prob(hap2, p_alt_2)

    return prob1 * prob2


def emission_logprob(
    state: str,
    hap1: str,
    hap2: str,
    pA: float,
    pB: float,
    eps: float = 1e-6,
) -> float:
    p = emission_prob(state, hap1, hap2, pA, pB, eps)
    return math.log10(p)


def site_emissions(
    hap1: Hap,
    hap2: Hap,
    params: EmissionParams,
) -> Dict[str, float]:
    """
    Compute emission logprobs for all states at a site.
    """
    return {
        s: emission_logprob(s, hap1, hap2, params.pA_alt, params.pB_alt, params.eps)
        for s in STATES
    }


def batch_site_emissions(
    haps: Iterable[Tuple[Hap, Hap]],
    params_list: Iterable[EmissionParams],
) -> Sequence[Dict[str, float]]:
    """
    Helper for building a full emissions list.
    """
    out: list[Dict[str, float]] = []
    for (h1, h2), params in zip(haps, params_list):
        out.append(site_emissions(h1, h2, params))
    return out


def assert_emissions_shape(emissions: Sequence[Mapping[str, float]]) -> None:
    """
    each site should have all states.
    """
    for i, e in enumerate(emissions):
        for s in STATES:
            if s not in e:
                raise ValueError(f"Missing state {s} at site {i}")