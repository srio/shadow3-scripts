"""
FIFA World Cup 2026 — Monte Carlo Simulation
Runs the full tournament N times and ranks teams by
how often they finish 1st, 2nd, 3rd, 4th.
"""

import numpy as np
import random
from collections import defaultdict

# ── Data ──────────────────────────────────────────────────────

nations = np.array([
    ["Mexico",       "South Korea",  "South Africa", "Czechia"      ],
    ["Canada",       "Bosnia-Herz.", "Qatar",        "Switzerland"  ],
    ["Brazil",       "Morocco",      "Haiti",        "Scotland"     ],
    ["USA",          "Paraguay",     "Australia",    "Curaçao"      ],
    ["Germany",      "Colombia",     "Côte d'Ivoire","Tunisia"      ],
    ["Spain",        "Cape Verde",   "Saudi Arabia", "Uruguay"      ],
    ["Belgium",      "Egypt",        "Iran",         "New Zealand"  ],
    ["Portugal",     "DR Congo",     "Uzbekistan",   "Ecuador"      ],
    ["France",       "Senegal",      "Iraq",         "Norway"       ],
    ["Argentina",    "Algeria",      "Austria",      "Jordan"       ],
    ["Netherlands",  "Japan",        "Türkiye",      "Sweden"       ],
    ["England",      "Croatia",      "Ghana",        "Panama"       ],
])

GROUPS   = list("ABCDEFGHIJKL")
PAIR_IDX = [(i, j) for i in range(4) for j in range(i+1, 4)]

SLOT_ALLOWED = {
    "M74": set("ABCDF"),  "M77": set("CDFGH"),
    "M79": set("CEFHI"),  "M80": set("EHIJK"),
    "M81": set("BEFIJ"),  "M82": set("AEHIJ"),
    "M85": set("EFGIJ"),  "M87": set("DEIJL"),
}

# ── Match simulators ──────────────────────────────────────────

def simulate_match(n1, n2):
    """Group phase: draw allowed. Returns (pts1, pts2)."""
    # if n1 == "Spain": return (3, 0)
    # if n2 == "Spain": return (0, 3)
    out = random.choices([(3,0),(1,1),(0,3)], weights=[0.25,0.50,0.25])[0]
    print(">>", n1, n2, "Points: ", out)
    return out

def simulate_match_no_draws(n1, n2):
    # if n1 == "Spain": return n1
    # if n2 == "Spain": return n2
    """Knockout: no draws, 50/50. One winner always."""
    out = random.choice([n1, n2])
    print(">>>", n1, n2, "Winner is: ", out)
    return out

# ── Full tournament → returns (champion, runner_up, third, fourth) ──

def run_tournament():

    # Group stage
    points = np.zeros((12, 4), dtype=int)
    for g in range(12):
        for i, j in PAIR_IDX:
            pts1, pts2 = simulate_match(nations[g,i], nations[g,j])
            points[g, i] += pts1
            points[g, j] += pts2

    standings = np.array([nations[g, np.argsort(-points[g])] for g in range(12)])
    pts_sorted = np.array([np.sort(points[g])[::-1] for g in range(12)])

    def first(g):  return standings[GROUPS.index(g), 0]
    def second(g): return standings[GROUPS.index(g), 1]
    def third(g):  return standings[GROUPS.index(g), 2]

    # Best 8 third-placed teams
    third_pts  = pts_sorted[:, 2]
    best8_idx  = set(np.argsort(-third_pts)[:8])
    adv_thirds = sorted([GROUPS[i] for i in best8_idx])

    # Assign thirds to slots
    remaining  = set(adv_thirds)
    assignment = {}
    slots = sorted(SLOT_ALLOWED.keys(),
                   key=lambda s: len(SLOT_ALLOWED[s] & remaining))
    for slot in slots:
        candidates = SLOT_ALLOWED[slot] & remaining
        if candidates:
            chosen = sorted(candidates)[0]
            assignment[slot] = chosen
            remaining.discard(chosen)

    def t(slot):
        grp = assignment.get(slot)
        return third(grp) if grp else "TBD"

    # Round of 32 fixtures
    r32 = [
        (second("A"), second("B")), (first("E"),  t("M74")),
        (first("F"),  second("C")), (first("C"),  second("F")),
        (first("I"),  t("M77")),    (second("E"), second("I")),
        (first("A"),  t("M79")),    (first("L"),  t("M80")),
        (first("D"),  t("M81")),    (first("G"),  t("M82")),
        (second("K"), second("L")), (first("H"),  second("J")),
        (first("B"),  t("M85")),    (first("J"),  second("H")),
        (first("K"),  t("M87")),    (second("D"), second("G")),
    ]

    # Knockout rounds — adjacent winners pair up each round
    def knockout_round(pairs):
        winners = [simulate_match_no_draws(n1, n2) for n1, n2 in pairs]
        return winners

    r32_w = knockout_round(r32)
    r16_w = knockout_round([(r32_w[i], r32_w[i+1]) for i in range(0,16,2)])
    qf_w  = knockout_round([(r16_w[i], r16_w[i+1]) for i in range(0, 8,2)])
    sf_w  = []
    sf_l  = []
    for n1, n2 in [(qf_w[i], qf_w[i+1]) for i in range(0, 4, 2)]:
        w = simulate_match_no_draws(n1, n2)
        sf_w.append(w)
        sf_l.append(n2 if w == n1 else n1)

    champion   = simulate_match_no_draws(sf_w[0], sf_w[1])
    runner_up  = sf_w[1] if champion == sf_w[0] else sf_w[0]
    third      = simulate_match_no_draws(sf_l[0], sf_l[1])
    fourth     = sf_l[1] if third == sf_l[0] else sf_l[0]

    return champion, runner_up, third, fourth

# ── Monte Carlo ───────────────────────────────────────────────

N = 1
counts = {rank: defaultdict(int) for rank in [1, 2, 3, 4]}

for _ in range(N):
    results = run_tournament()
    for rank, team in enumerate(results, 1):
        counts[rank][team] += 1

# ── Results ───────────────────────────────────────────────────

all_teams = sorted(nations.flatten())

print(f"\n{'='*65}")
print(f"  MONTE CARLO — {N} simulations")
print(f"{'='*65}")
print(f"  {'Team':<22} {'1st':>7} {'2nd':>7} {'3rd':>7} {'4th':>7}  {'1st %':>7}")
print(f"  {'-'*63}")

# Sort by number of 1st place finishes
ranked = sorted(all_teams,
                key=lambda t: counts[1][t],
                reverse=True)

for team in ranked:
    c1 = counts[1][team]
    c2 = counts[2][team]
    c3 = counts[3][team]
    c4 = counts[4][team]
    if c1 + c2 + c3 + c4 == 0:
        continue
    pct = c1 / N * 100
    bar = "█" * int(pct / 2)
    print(f"  {team:<22} {c1:>7} {c2:>7} {c3:>7} {c4:>7}  {pct:>6.1f}%  {bar}")

print(f"\n  Total simulations: {N}")
print(f"  Each team had equal 50/50 probability in every match")
print(f"  → Expected champion rate ≈ 1/48 = {100/48:.1f}%  (pure random)")