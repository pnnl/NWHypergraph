Installing NWhy
===============

The NWhy library provides Pybind11_ APIs for analysis of complex data set intepret as hypergraphs.

.. _Pybind11: https://github.com/pybind/pybind11

To install in an Anaconda environment
-------------------------------------

	>>> conda create -n <env name> python=3.8

Then activate the environment
-----------------------------

	>>> conda activate <env name> 

Install Intel Threading Building Blocks(TBB)
--------------------------------------------

To install TBB_:

.. _TBB: https://github.com/oneapi-src/oneTBB

	>>> conda install tbb

Install using Pip
-----------------

For installation:

	>>> pip install nwhy

For upgrade:

	>>> pip install nwhy --upgrade

or 

	>>> pip install nwhy -U


Quick test with import
----------------------

For quick test:

	>>> python -c "import nwhy"

If there is no import error, then installation is done.