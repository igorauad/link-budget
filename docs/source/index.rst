.. Link Budget documentation master file, created by
   sphinx-quickstart on Thu May 20 21:48:04 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python-based Link Budget Calculator
=======================================


A Python-based link budget calculator for satellite communications and radar
systems.

Installation
==================

Package `link-budget` provides the link budget calculator as a command-line
tool. To install it, run:

.. code-block:: bash

   pip install link-budget

Basic Usage
======================

The `link-budget` utility accepts a range of parameters. The following example
computes the link budget for a 52 dBW EIRP at 12.45 GHz, with an
intermediate-frequency (IF) bandwidth of 24 MHz. The Rx dish has a diameter of
0.46 m, while the LNB has a conversion gain of 40 dB and a noise figure of 0.6
dB. Furthermore, the LNB connects to a receiver with a noise figure of 10 dB
over a coaxial cable with 110 ft. Lastly, this example receiver station is at
longitude -82.43 degrees and latitude 29.71 degrees.

.. code-block:: bash

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

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   cli
   calc
   antenna
   pointing
   propagation
   util

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
