"""
Nation Strength Score
Combines 6 features into a single strength value per team.

DEFAULT_WEIGHTS: set any value you like — they are auto-normalized to sum=1.
                 Set a weight to None or negative → replaced by random value.
"""

import numpy as np
import random

# ── Nations matrix (12x4) ─────────────────────────────────────

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

GROUPS = list("ABCDEFGHIJKL")

# ── Raw feature matrices (12x4) ───────────────────────────────

world_cups_won = np.array([
    [0, 0, 0, 0], [0, 0, 0, 0], [5, 0, 0, 0], [0, 0, 0, 0],
    [4, 0, 0, 0], [1, 0, 0, 2], [0, 0, 0, 0], [0, 0, 0, 0],
    [2, 0, 0, 0], [3, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0],
], dtype=float)

fifa_ranking = np.array([
    [15,  22,  61,  44], [21,  65,  55,  19], [ 6,   8,  83,  43],
    [16,  40,  26,  82], [10,  13,  51,  44], [ 2,  69,  61,  17],
    [ 9,  35,  33,  85], [ 5,  46,  50,  61], [ 1,  14,  57,  30],
    [ 3,  32,  29,  63], [ 7,  18,  28,  24], [ 4,  11,  74,  38],
], dtype=float)

distance_km = np.array([
    [ 1502, 10973, 14685,  8595], [ 2288,  9258, 12739,  8369],
    [ 7490,  8044,  2906,  7253], [ 1902,  7673, 14051,  3651],
    [ 8378,  3912,  9946,  9190], [ 7956,  7572, 12547,  8635],
    [ 7951, 11143, 11712, 12217], [ 7624, 12330, 11655,  4140],
    [ 7934,  8145, 11616,  7788], [ 8499,  8669,  8834, 11275],
    [ 7901, 10397, 10416,  8164], [ 7641,  8975, 10286,  3186],
], dtype=float)

squad_value_mEUR = np.array([
    [164.6, 184.3,  12.0,  85.0], [185.0,  60.0,  23.2, 281.7],
    [1135.0,318.1,   8.0, 115.0], [270.2,  55.0,  41.3,  30.0],
    [775.0, 211.8,  90.0,  54.1], [861.0,  25.0,  14.5, 424.0],
    [549.4,  75.0,  51.4,  15.0], [1000.0, 18.0,  20.0, 236.4],
    [1195.0,211.8,  24.5, 180.0], [820.7,  40.0, 150.0,  23.1],
    [671.7, 284.8, 190.0, 120.0], [1345.0,325.5,  25.0,  45.0],
], dtype=float)

goals_scored_pg = np.array([
    [1.8, 2.1, 1.2, 1.9], [2.0, 1.4, 1.0, 2.2], [2.3, 2.0, 0.8, 1.8],
    [2.1, 1.6, 1.5, 1.4], [2.4, 2.1, 1.9, 1.5], [2.5, 1.3, 1.2, 1.8],
    [2.0, 1.4, 1.3, 1.0], [2.2, 1.0, 1.1, 1.7], [2.8, 1.9, 1.1, 4.6],
    [1.7, 1.5, 2.0, 1.2], [2.3, 3.2, 1.8, 2.1], [2.4, 1.7, 1.5, 1.6],
], dtype=float)

goals_conceded_pg = np.array([
    [1.2, 1.4, 1.8, 1.3], [1.1, 1.5, 1.9, 0.9], [0.8, 0.9, 2.5, 1.2],
    [0.9, 1.7, 1.4, 1.8], [1.0, 1.5, 1.6, 1.7], [0.7, 1.6, 1.8, 1.1],
    [1.1, 1.5, 1.4, 1.5], [0.6, 2.1, 1.7, 1.2], [0.5, 1.0, 2.0, 0.6],
    [0.6, 1.4, 1.0, 1.8], [0.8, 0.7, 1.3, 0.9], [0.7, 1.1, 1.8, 1.4],
], dtype=float)

population_M = np.array([
    [130.0,  51.7,  60.6,  10.9], [ 38.2,   3.3,   2.9,   8.9],
    [215.3,  37.5,  11.7,   5.5], [335.0,   6.8,  26.5,   0.16],
    [ 84.4,  52.2,  28.0,  12.0], [ 47.4,   0.6,  36.9,   3.5],
    [ 11.6, 106.0,  89.8,   5.1], [ 10.3,  99.0,  35.3,  18.1],
    [ 68.0,  17.6,  42.3,   5.5], [ 45.8,  45.6,   9.1,  10.8],
    [ 17.9, 124.5,  85.7,  10.5], [ 56.5,   4.0,  33.5,   4.4],
], dtype=float)

# ── Normalized feature matrices (computed once) ───────────────

def _norm_minmax(m, invert=False):
    mn, mx = m.min(), m.max()
    n = (m - mn) / (mx - mn + 1e-9)
    return 1.0 - n if invert else n

def _norm_log(m, invert=False):
    lg = np.log1p(m)
    mn, mx = lg.min(), lg.max()
    n = (lg - mn) / (mx - mn + 1e-9)
    return 1.0 - n if invert else n

goal_ratio = goals_scored_pg / (goals_scored_pg + goals_conceded_pg + 1e-9)

FEATURES = {
    "rank":     _norm_minmax(fifa_ranking,   invert=True),
    "value":    _norm_log(squad_value_mEUR,  invert=False),
    "form":     _norm_minmax(goal_ratio,     invert=False),
    "pop":      _norm_log(population_M,      invert=False),
    "distance": _norm_minmax(distance_km,    invert=True),
    "wc":       _norm_minmax(world_cups_won, invert=False),
}

# ── DEFAULT WEIGHTS ───────────────────────────────────────────
# Rules:
#   - Any positive number → used as-is (auto-normalized to sum=1)
#   - None or negative    → replaced by a random positive value
#   - No need to sum to 1 yourself — normalization is automatic

DEFAULT_WEIGHTS = {
    "rank":     0.35,
    "value":    0.25,
    "form":     0.20,
    "pop":      0.10,
    "distance": 0.05,
    "wc":       0.05,
}

# ── compute_strength ──────────────────────────────────────────

def compute_strength(weights=DEFAULT_WEIGHTS):
    """
    Compute strength lookup dict {team: score} from weights.

    weights: dict with keys matching FEATURES.
             - Positive value → used as weight.
             - None or negative → replaced by random().
             Auto-normalized to sum=1.
    Returns: (lookup dict, resolved & normalized weights dict)
    """
    # Step 1: resolve None / negative → random
    resolved = {
        k: (random.random() if (v is None or v < 0) else v)
        for k, v in weights.items()
    }

    # Step 2: normalize to sum=1
    total = sum(resolved.values())
    normalized = {k: v / total for k, v in resolved.items()}

    # Step 3: weighted sum of feature matrices
    strength = sum(normalized[k] * FEATURES[k] for k in normalized)

    # Step 4: build lookup dict
    lookup = {nations[g, c]: strength[g, c]
              for g in range(12) for c in range(4)}

    return lookup, normalized


# ── print_strength ────────────────────────────────────────────

def print_strength(lookup, weights):
    print(f"\n{'='*65}")
    print(f"  Strength Ranking")
    print(f"  Weights (normalized):")
    for k, v in weights.items():
        print(f"    {k:<10} {v:.3f}")
    print(f"{'='*65}")
    print(f"  {'#':<4} {'Team':<22} {'Score':>7}  Bar")
    print(f"  {'-'*55}")
    for i, (team, score) in enumerate(
            sorted(lookup.items(), key=lambda x: x[1], reverse=True), 1):
        bar = "█" * int(score * 30)
        print(f"  {i:<4} {team:<22} {score:>7.4f}  {bar}")


# ── Demo ──────────────────────────────────────────────────────

if __name__ == "__main__":

    # 1. Default weights
    print("\n── Default weights ──────────────────────────────────────")
    lookup, w = compute_strength(DEFAULT_WEIGHTS)
    print_strength(lookup, w)

    # 2. Random weights (all None)
    print("\n\n── All-random weights ───────────────────────────────────")
    lookup_r, w_r = compute_strength({k: None for k in DEFAULT_WEIGHTS})
    print_strength(lookup_r, w_r)

    # 3. Mix: some fixed, some random
    print("\n\n── Mixed: rank fixed, rest random ───────────────────────")
    lookup_m, w_m = compute_strength({
        "rank":     0.5,   # fixed
        "value":    None,  # random
        "form":     None,  # random
        "pop":      -1,    # random (negative)
        "distance": None,  # random
        "wc":       None,  # random
    })
    print_strength(lookup_m, w_m)