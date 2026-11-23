"""Regression baselines for phylogenetic diversity measures.

Stores expected outcomes for experiment networks (exp1) and the shared
``small_tree`` fixture so tests can assert against fixed values.

Notes
-----
``EXP1_FIXED_SUBSET_DIVERSITY_BASELINES`` uses fixed taxon sets. They were
originally the APD-optimal sets, but since edges of length 0 are no longer
counted as length 1 (version 2.1.2), most of them differ from the current
``APD_BASELINES`` solutions; they are kept as they are. For networks ``00101``
and ``00401`` with budget 1, ``min_tree.compute_diversity`` cannot induce the
subnetwork for some singletons (PhyloZoo parallel-edge validation); those entries
use the singleton taxa ``t1`` and ``t10`` respectively.

All values were checked against brute-force computations (all-paths diversity
over all taxon subsets, and MaxTreePD/MinTreePD over all switchings).
"""

from __future__ import annotations

from typing import Any

APD_BASELINES: dict[str, dict[int, dict[str, Any]]] = {
    "00001": {
        1: {"solution": ["t4"], "value": 6.61566628},
        2: {"solution": ["t15", "t4"], "value": 13.049230782599999},
        3: {"solution": ["t15", "t3", "t4"], "value": 18.4358022246},
        5: {"solution": ["t15", "t18", "t3", "t4", "t8"], "value": 22.782559121299997},
    },
    "00101": {
        1: {"solution": ["t20"], "value": 3.776914961761},
        2: {"solution": ["t14", "t20"], "value": 5.726578479760999},
        3: {"solution": ["t11", "t14", "t20"], "value": 7.514500559360999},
        5: {"solution": ["t1", "t11", "t14", "t16", "t20"], "value": 10.410450645661},
    },
    "00201": {
        1: {"solution": ["t2"], "value": 6.758689968600001},
        2: {"solution": ["t2", "t9"], "value": 9.33579965345},
        3: {"solution": ["t2", "t8", "t9"], "value": 11.678709709769997},
        5: {"solution": ["t12", "t2", "t6", "t8", "t9"], "value": 14.47370332461},
    },
    "00301": {
        1: {"solution": ["t18"], "value": 9.350581974344},
        2: {"solution": ["t11", "t18"], "value": 14.174457809244},
        3: {"solution": ["t11", "t18", "t21"], "value": 17.825392599044},
        5: {"solution": ["t11", "t18", "t21", "t7", "t8"], "value": 21.633099618044003},
    },
    "00401": {
        1: {"solution": ["t1"], "value": 6.79687361743},
        2: {"solution": ["t1", "t7"], "value": 10.78769663536},
        3: {"solution": ["t1", "t21", "t7"], "value": 13.02042449136},
        5: {"solution": ["t1", "t13", "t18", "t21", "t7"], "value": 14.76729930646},
    },
}

APD_BUDGETS: tuple[int, ...] = (1, 2, 3, 5)

# MaxTreePD under unit costs (may differ from all-paths APD on reticulate networks).
# Only objective values are pinned: optimal taxon sets can tie under the DP tie-breaking.
MAX_TREE_BUDGET_BASELINES: dict[str, dict[int, float]] = {
    "00001": {1: 6.61566628, 2: 13.049230782599999, 3: 18.4358022246, 5: 22.782559121299997},
    "00101": {1: 2.261911132161, 2: 4.493075351760999, 3: 6.442738869660999, 5: 10.016015764760999},
    "00201": {1: 5.70226862625, 2: 9.10160002537, 3: 11.678709709769997, 5: 14.47370332461},
    "00301": {1: 5.568313233674, 2: 9.398525289674, 3: 13.227885574074001, 5: 19.071657187704},
    "00401": {1: 3.6377712514300002, 2: 6.86585837909, 3: 9.1854516432, 5: 12.66698782813},
}

# Fixed taxon sets: all-paths diversity vs MinTreePD (displayed-tree minimum) on that set.
EXP1_FIXED_SUBSET_DIVERSITY_BASELINES: dict[str, dict[int, dict[str, Any]]] = {
    "00001": {
        1: {"taxa": ["t3"], "all_paths": 6.61566628, "min_tree": 6.61566628},
        2: {"taxa": ["t15", "t3"], "all_paths": 13.049230782599999, "min_tree": 13.049230782599999},
        3: {"taxa": ["t15", "t3", "t5"], "all_paths": 18.4358022246, "min_tree": 18.4358022246},
        5: {
            "taxa": ["t13", "t15", "t18", "t3", "t5"],
            "all_paths": 22.7825591213,
            "min_tree": 22.7825591213,
        },
    },
    "00101": {
        1: {"taxa": ["t1"], "all_paths": 2.2619111317609994, "min_tree": 2.261911131761},
        2: {"taxa": ["t14", "t20"], "all_paths": 5.726578479760999, "min_tree": 4.211574649399999},
        3: {
            "taxa": ["t14", "t20", "t23"],
            "all_paths": 6.577806633520998,
            "min_tree": 4.9277350392449994,
        },
        5: {
            "taxa": ["t11", "t14", "t15", "t20", "t23"],
            "all_paths": 10.151083528621001,
            "min_tree": 8.946821338706,
        },
    },
    "00201": {
        1: {"taxa": ["t2"], "all_paths": 6.758689968600001, "min_tree": 5.7022686258},
        2: {"taxa": ["t2", "t9"], "all_paths": 9.335799653450001, "min_tree": 8.27937831065},
        3: {
            "taxa": ["t15", "t2", "t9"],
            "all_paths": 11.678709709769999,
            "min_tree": 10.937423458769999,
        },
        5: {
            "taxa": ["t15", "t2", "t24", "t6", "t9"],
            "all_paths": 14.473703324609998,
            "min_tree": 13.732417073609998,
        },
    },
    "00301": {
        1: {"taxa": ["t18"], "all_paths": 9.350581974344, "min_tree": 5.5683132327700005},
        2: {"taxa": ["t18", "t23"], "all_paths": 14.174457809244002, "min_tree": 9.017098883144001},
        3: {
            "taxa": ["t18", "t20", "t23"],
            "all_paths": 17.825392599043997,
            "min_tree": 12.847310939144,
        },
        5: {
            "taxa": ["t18", "t20", "t23", "t7", "t8"],
            "all_paths": 21.633099618044,
            "min_tree": 16.655017958144,
        },
    },
    "00401": {
        1: {"taxa": ["t10"], "all_paths": 3.6377712501899997, "min_tree": 3.63777125019},
        2: {"taxa": ["t1", "t22"], "all_paths": 9.026861714919999, "min_tree": 3.9354776114199996},
        3: {
            "taxa": ["t1", "t21", "t22"],
            "all_paths": 12.08548716412,
            "min_tree": 5.050392373419999,
        },
        5: {
            "taxa": ["t1", "t13", "t21", "t22", "t7"],
            "all_paths": 14.644732605950004,
            "min_tree": 10.23585942413,
        },
    },
}

EXP1_MIN_TREE_FULL_LEAF_DIVERSITY: dict[str, float] = {
    "00001": 33.18148433767,
    "00101": 15.631596940151,
    "00201": 19.76914009861,
    "00301": 26.972085813173997,
    "00401": 17.134867692660002,
}

SMALL_TREE_BUDGET_COSTS: dict[str, int] = {"a": 2, "b": 1, "c": 2, "d": 1}

SMALL_TREE_BUDGET_MAXIMIZATION: dict[str, Any] = {
    "budget": 3,
    "costs": SMALL_TREE_BUDGET_COSTS,
    "value": 7.5,
    "taxa": frozenset({"a", "d"}),
}

SMALL_TREE_FIXED_DIVERSITY: dict[frozenset[str], float] = {
    frozenset({"a", "c", "d"}): 10.5,
    frozenset({"a", "b", "c", "d"}): 12.5,
}
