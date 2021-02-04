Installing NWhy
===============

The NWhy library provides Pybind11 APIs for analysis of complex data set intepret as hypergraphs.


To install in an Anaconda environment
-------------------------------------

	>>> conda create -n <env name> python=3.8
	>>> conda activate <env name> 

Install Intel Threading Building Blocks(TBB)
--------------------------------------------

For installation:

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
