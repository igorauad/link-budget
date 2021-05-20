# Link Budget

![Tests](https://github.com/igorauad/link_budget/workflows/Tests/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/igorauad/link_budget/branch/master/graph/badge.svg?token=72U3BI51OT)](https://codecov.io/gh/igorauad/link_budget)

A Python-based link budget calculator for basic satellite communications and
radar systems.

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-generate-toc again -->
**Table of Contents**

- [Link Budget](#link-budget)
    - [Installation](#installation)
    - [Running](#running)

<!-- markdown-toc end -->


## Installation

Package `link-budget` provides the link budget calculator as a command-line
tool. To build the package and install it, run:

```
make && make install
```

## Running

Example:

```
link-budget \
  --eirp 52 \
  --freq 12.45e9 \
  --if-bw 24e6 \
  --rx-dish-size 0.46 \
  --antenna-noise-temp 20 \
  --lnb-noise-fig 0.6 \
  --lnb-gain 40 \
  --coax-length 110 \
  --rx-noise-fig 10 \
  --sat-long -101 \
  --rx-long -82.43 \
  --rx-lat 29.71
```
