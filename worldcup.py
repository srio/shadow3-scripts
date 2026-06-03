"""
FIFA World Cup 2026 — Full Tournament Simulation + Monte Carlo
Imports strength scores from worldcup_features.py

Usage:
    python worldcup.py              # single run, verbose
    python worldcup.py --mc 1000    # Monte Carlo, 1000 runs
"""

import sys
import random
import numpy as np
from collections import defaultdict

# ── Import strength lookup from features module ───────────────
from worldcup_features import (
    nations, GROUPS, DEFAULT_WEIGHTS, RANDOM_WEIGHTS, compute_strength
)

# ── Match simulators ──────────────────────────────────────────

PAIR_IDX = [(i, j) for i in range(4) for j in range(i+1, 4)]

def simulate_match(n1, n2, strength):
    """
    Group phase: draw allowed.
    Weights biased by strength difference.
    Returns (pts1, pts2).
    """
    s1, s2 = strength[n1], strength[n2]
    total  = s1 + s2
    p1     = s1 / total   # prob n1 wins
    p2     = s2 / total   # prob n2 wins
    pd     = 0.50         # draw factor — both sides scaled down

    w_win1  = p1 * (1 - pd)
    w_draw  = pd
    w_win2  = p2 * (1 - pd)

    return random.choices(
        [(3, 0), (1, 1), (0, 3)],
        weights=[w_win1, w_draw, w_win2]
    )[0]

def simulate_match_no_draws(n1, n2, strength):
    """
    Knockout phase: no draws, one winner always.
    Probability proportional to strength scores.
    """
    s1, s2 = strength[n1], strength[n2]
    return random.choices([n1, n2], weights=[s1, s2])[0]

# ── Slot allowed for 3rd-place teams ─────────────────────────

SLOT_ALLOWED = {
    "M74": set("ABCDF"),  "M77": set("CDFGH"),
    "M79": set("CEFHI"),  "M80": set("EHIJK"),
    "M81": set("BEFIJ"),  "M82": set("AEHIJ"),
    "M85": set("EFGIJ"),  "M87": set("DEIJL"),
}

def assign_thirds(advancing):
    """Backtracking to guarantee complete assignment."""
    slots = list(SLOT_ALLOWED.keys())
    remaining = set(advancing)
    assignment = {}
    def backtrack(idx):
        if idx == len(slots): return True
        slot = slots[idx]
        for chosen in sorted(SLOT_ALLOWED[slot] & remaining):
            assignment[slot] = chosen; remaining.discard(chosen)
            if backtrack(idx + 1): return True
            del assignment[slot]; remaining.add(chosen)
        return False
    backtrack(0)
    return assignment
def _assign_thirds_old(advancing):
    remaining  = set(advancing)
    assignment = {}
    slots = sorted(SLOT_ALLOWED.keys(),
                   key=lambda s: len(SLOT_ALLOWED[s] & remaining))
    for slot in slots:
        candidates = SLOT_ALLOWED[slot] & remaining
        if candidates:
            chosen = sorted(candidates)[0]
            assignment[slot] = chosen
            remaining.discard(chosen)
    return assignment

# ── Full tournament ───────────────────────────────────────────

def run_tournament(weights=DEFAULT_WEIGHTS, verbose=False):
    """
    Run one full tournament simulation.
    weights: passed to compute_strength() from worldcup_features.
    Returns (champion, runner_up, third, fourth).
    """

    # Compute strength scores for this run
    strength, _ = compute_strength(weights)

    # ── Group stage ───────────────────────────────────────────
    points = np.zeros((12, 4), dtype=int)
    for g in range(12):
        for i, j in PAIR_IDX:
            pts1, pts2 = simulate_match(
                nations[g, i], nations[g, j], strength)
            points[g, i] += pts1
            points[g, j] += pts2

    standings = np.array([
        nations[g, np.argsort(-points[g])] for g in range(12)
    ])
    pts_sorted = np.array([
        np.sort(points[g])[::-1] for g in range(12)
    ])

    def first(g):  return standings[GROUPS.index(g), 0]
    def second(g): return standings[GROUPS.index(g), 1]
    def third_p(g):  return standings[GROUPS.index(g), 2]

    if verbose:
        print(f"\n{'='*55}")
        print(f"  GROUP STAGE STANDINGS")
        print(f"{'='*55}")
        print(f"  {'Grp':<5} {'1st':<22} {'2nd':<22} {'3rd'}")
        print(f"  {'-'*65}")
        for g in GROUPS:
            print(f"  {g:<5} {first(g):<22} {second(g):<22} {third_p(g)}")

    # ── Best 8 third-placed teams ─────────────────────────────
    third_pts  = pts_sorted[:, 2]
    best8_idx  = set(np.argsort(-third_pts)[:8])
    adv_thirds = sorted([GROUPS[i] for i in best8_idx])
    assignment = assign_thirds(set(adv_thirds))

    def t(slot):
        grp = assignment.get(slot)
        return third_p(grp) if grp else "TBD"

    # ── Round of 32 fixtures ──────────────────────────────────
    r32 = [
        (second("A"), second("B")), (first("E"),   t("M74")),
        (first("F"),  second("C")), (first("C"),   second("F")),
        (first("I"),  t("M77")),    (second("E"),  second("I")),
        (first("A"),  t("M79")),    (first("L"),   t("M80")),
        (first("D"),  t("M81")),    (first("G"),   t("M82")),
        (second("K"), second("L")), (first("H"),   second("J")),
        (first("B"),  t("M85")),    (first("J"),   second("H")),
        (first("K"),  t("M87")),    (second("D"),  second("G")),
    ]

    # ── Knockout rounds ───────────────────────────────────────
    def knockout_round(pairs, label=""):
        winners = []
        losers  = []
        if verbose and label:
            print(f"\n{'='*55}")
            print(f"  {label}  ({len(pairs)} matches)")
            print(f"{'='*55}")
            print(f"  {'Team 1':<24} {'Team 2':<24} {'Winner'}")
            print(f"  {'-'*58}")
        for n1, n2 in pairs:
            w = simulate_match_no_draws(n1, n2, strength)
            l = n2 if w == n1 else n1
            winners.append(w)
            losers.append(l)
            if verbose and label:
                print(f"  {n1:<24} {n2:<24} {w}  (out: {l})")
        return winners, losers

    r32_w, _   = knockout_round(r32,  "Round of 32")
    r16_w, _   = knockout_round(
        [(r32_w[i], r32_w[i+1]) for i in range(0, 16, 2)], "Round of 16")
    qf_w,  _   = knockout_round(
        [(r16_w[i], r16_w[i+1]) for i in range(0,  8, 2)], "Quarter-finals")
    sf_w,  sf_l = knockout_round(
        [(qf_w[i],  qf_w[i+1])  for i in range(0,  4, 2)], "Semi-finals")

    third_w, _ = knockout_round(
        [(sf_l[0], sf_l[1])], "3rd Place Match")
    fourth     = sf_l[1] if third_w[0] == sf_l[0] else sf_l[0]

    final_w, _ = knockout_round(
        [(sf_w[0], sf_w[1])], "FINAL")
    champion   = final_w[0]
    runner_up  = sf_w[1] if champion == sf_w[0] else sf_w[0]

    if verbose:
        print(f"\n{'='*55}")
        print(f"  FINAL STANDINGS")
        print(f"{'='*55}")
        print(f"  🥇  {champion}")
        print(f"  🥈  {runner_up}")
        print(f"  🥉  {third_w[0]}")
        print(f"  4th {fourth}")
        print(f"{'='*55}")

    return champion, runner_up, third_w[0], fourth

# ── Monte Carlo ───────────────────────────────────────────────

def run_monte_carlo(n=1000, weights=DEFAULT_WEIGHTS):
    """
    Run n tournament simulations and print champion frequency.
    If any weight is None/negative, a NEW random value is drawn
    each simulation — full uncertainty mode.
    """
    counts = {rank: defaultdict(int) for rank in [1, 2, 3, 4]}

    for _ in range(n):
        results = run_tournament(weights)
        for rank, team in enumerate(results, 1):
            counts[rank][team] += 1

    all_teams = sorted(nations.flatten())

    print(f"\n{'='*68}")
    print(f"  MONTE CARLO — {n} simulations")
    print(f"{'='*68}")
    print(f"  {'Team':<22} {'1st':>6} {'2nd':>6} {'3rd':>6} {'4th':>6}  {'1st%':>6}  Bar")
    print(f"  {'-'*65}")

    ranked = sorted(all_teams,
                    key=lambda t: counts[1][t], reverse=True)
    for team in ranked:
        c1, c2, c3, c4 = (counts[r][team] for r in [1, 2, 3, 4])
        if c1+c2+c3+c4 == 0:
            continue
        pct = c1 / n * 100
        bar = "█" * int(pct / 1.5)
        print(f"  {team:<22} {c1:>6} {c2:>6} {c3:>6} {c4:>6}  {pct:>5.1f}%  {bar}")

    print(f"\n  Random baseline: 1/48 = {100/48:.1f}%")

# ── Entry point ───────────────────────────────────────────────

if __name__ == "__main__":
    # if "--mc" in sys.argv:
    #     idx = sys.argv.index("--mc")
    #     n   = int(sys.argv[idx + 1]) if idx + 1 < len(sys.argv) else 1000
    #     run_monte_carlo(n=n, weights=DEFAULT_WEIGHTS)
    # else:
    #     run_tournament(weights=DEFAULT_WEIGHTS, verbose=True)

    n = 100000
    weights = RANDOM_WEIGHTS
    if n > 1:
        run_monte_carlo(n=n, weights=weights)
    else:
        run_tournament(weights=weights, verbose=True)