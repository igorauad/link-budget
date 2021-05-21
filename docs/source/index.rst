.. Link Budget documentation master file, created by
   sphinx-quickstart on Thu May 20 21:48:04 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python-based Link Budget Calclator
=======================================


A Python-based link budget calculator for basic satellite communications and
radar systems.

Installation
==================

Package `link-budget` provides the link budget calculator as a command-line
tool. To build the package and install it, run:

.. code-block:: bash

   make && make install

Command-line Interface
======================

Example:


.. code-block:: bash

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

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   main
   calc
   antenna
   pointing
   util

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
