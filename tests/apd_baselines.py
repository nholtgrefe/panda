"""Baseline APD regression values for selected experiment networks."""

from __future__ import annotations

APD_BASELINES: dict[str, dict[int, dict[str, object]]] = {'00001': {1: {'solution': ['t3'], 'value': 6.61566628},
           2: {'solution': ['t15', 't3'], 'value': 13.049230782599999},
           3: {'solution': ['t15', 't3', 't5'], 'value': 18.4358022246},
           5: {'solution': ['t13', 't15', 't18', 't3', 't5'], 'value': 22.782559121299997}},
 '00101': {1: {'solution': ['t20'], 'value': 5.7769149617609985},
           2: {'solution': ['t14', 't20'], 'value': 7.726578479760999},
           3: {'solution': ['t14', 't20', 't23'], 'value': 9.577806633521},
           5: {'solution': ['t11', 't14', 't15', 't20', 't23'], 'value': 13.151083528621}},
 '00201': {1: {'solution': ['t2'], 'value': 8.7586899686},
           2: {'solution': ['t2', 't9'], 'value': 11.33579965345},
           3: {'solution': ['t15', 't2', 't9'], 'value': 13.67870970977},
           5: {'solution': ['t15', 't2', 't24', 't6', 't9'], 'value': 16.47370332461}},
 '00301': {1: {'solution': ['t18'], 'value': 13.350581974344001},
           2: {'solution': ['t18', 't23'], 'value': 18.174457809244},
           3: {'solution': ['t18', 't20', 't23'], 'value': 21.825392599044},
           5: {'solution': ['t18', 't20', 't23', 't7', 't8'], 'value': 25.633099618044}},
 '00401': {1: {'solution': ['t1'], 'value': 8.79687361743},
           2: {'solution': ['t1', 't22'], 'value': 13.02686171492},
           3: {'solution': ['t1', 't21', 't22'], 'value': 17.08548716412},
           5: {'solution': ['t1', 't13', 't21', 't22', 't7'], 'value': 19.64473260595}}}

APD_BUDGETS: tuple[int, ...] = (1, 2, 3, 5)
