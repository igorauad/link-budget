# Link Budget

![Tests](https://github.com/igorauad/link_budget/workflows/Tests/badge.svg?branch=master)
[![Docs](https://readthedocs.org/projects/link-budget/badge/?version=latest)](https://link-budget.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/igorauad/link-budget/branch/master/graph/badge.svg?token=72U3BI51OT)](https://codecov.io/gh/igorauad/link_budget)

A Python-based link budget calculator for satellite communications and radar
systems. To install it, run:

```
pip install link-budget
```

Please refer to the [documentation page](https://link-budget.readthedocs.io/)
for further information.

## Running

Example:

```
link-budget \
  --eirp 52 \
  --freq 12.45e9 \
  --if-bw 24e6 \
  --rx-dish-size 0.46 \
  --lnb-noise-fig 0.6 \
  --lnb-gain 40 \
  --coax-length 110 \
  --rx-noise-fig 10 \
  --sat-long -101 \
  --rx-long -82.43 \
  --rx-lat 29.71
```
