# NWHy: Northwest Hypergraph Processing Framework
NWHy is a C++ hypergraph processing framework for shared memory. 
It focuses on constructing s-line graphs, a lower-order approximation of a given hypergraph, and compute different metrics such as s-connected components, s-betweenness centrality, s-closeness centrality, etc. It also provides Python APIs for s-overlap computation. The Python APIs are provided using [Pybind11](https://pybind11.readthedocs.io/en/stable/). Pybind11 is included as a git submodule.

## Organization

The organization of our library is shown as follow:
```
$NWHy_HOME/
├── README.md
├── CMakeLists.txt
├── bench/
|   ├── CMakeLists.txt
|   ├── adjoinbfs.cpp
|   ├── adjoincc.cpp
|   ├── common.hpp
|   ├── hyperbfs.cpp
|   ├── hypercc.cpp
|   ├── imdb.cpp
|   ├── Log.hpp
|   ├── soverlapbc.cpp
|   ├── soverlapbfs.cpp
|   ├── soverlapcc.cpp
|   ├── soverlapsssp.cpp
|   └── toplexes.cpp
├── docker/
│   └── Dockerfile.gcc11
├── include/
│   ├── algorithms/
|   |   └── experimental/
|   ├── containers/
|   ├── generators/
|   ├── io/
|   ├── util/
|   ├── CMakeLists.txt
|   └── s_overlap.hpp
├── python/
│   ├── pybind11/
|   ├── test/
|   ├── CMakeLists.txt
|   ├── nwhy.cpp
|   └── ...
├── test/
├── setup.cfg
├── setup.py
└── LICENSE
```

The algorithms for hypergraph and s-line graphs are under `NWHy_HOME/include/algorithms/` diretory. The applications for hypergraph analytics are under `NWHy_HOME/bench/` diretory. The Python API definitions are under `NWHy_HOME/python/` directory. The Pytest cases are under `NWHy_HOME/python/test/` diretory. The pybind11 is a git module fetched from its Github. The C++ test cases are under `NWHy_HOME/test/` diretory.

## How to compile

NWHy is built upon NWGraph (NWGr) library. It uses NWGr graph abstractions and concepts as the building blocks for hypergraph models, containers and implementations in NWHy. NWHy also uses many graph algorithms in NWGr. Sepecifically, NWHy relies on a the master branch of NWGr. NWHy (as well as NWGr) uses [Intel OneTBB](https://github.com/oneapi-src/oneTBB) as the parallel backend.   



### Requirements

* g++ &gt;= 11 with support for OneTBB as parallel backend
* oneTBB &gt;= 2021
* NWGr master branch

Note that older versions of g++ are supported along with oneTBB 2020. 

### Compilation

```
$ mkdir build; cd build
$ cmake ..
```

### Useful things to know 
To specify compiler:
```
$ cmake .. -DCMAKE_CXX_COMPILER=g++-11
```
To specify build type as Release or Debug, default is Release:
```
$ cmake .. -DCMAKE_BUILD_TYPE=Release (or Debug)
```
To enable test cases and examples under build/test directory:
```
$ cmake .. -DNW_HYPERGRAPH_BUILD_TEST=ON (or OFF)
```
To generate applications under build/bench/ directory:
```
$ cmake .. -DNW_HYPERGRAPH_BUILD_BENCH=ON (or OFF)
```
To enable Python binding modules compilation, and generate Python wheel file:
```
$ cmake .. -DNW_HYPERGRAPH_BUILD_PYBIND=OFF (or ON)
```
To generate tools under build/tools/ directory:
```
$ cmake .. -DNW_HYPERGRAPH_BUILD_TOOLS=OFF (or ON)
```
To see verbose information during compilation:
```
$ make VERBOSE=1
```
To run C++ test cases after compilation, assume the test cases are enabled during cmake:
```
$ make test
```

## Running code in NWHy

NWHy uses command-line interface description language [DOCOPT](http://docopt.org/) to define the interface of our command-line applications and tools.

A typical interface of the binary looks like this:
```
$ hystats.exe: hypergraph stats driver.
  Usage:
      hystats.exe (-h | --help)
      hystats.exe [-f FILE...] [--deg-dist FILE...] [--output FILE...] [-D NUM] [--degree NUM] [--relabel NUM] [--direction DIR] [-d] [--log FILE] [--log-header]

  Options:
      -h, --help            show this screen
      -f FILE               edge list or matrix market input file paths
      --deg-dist FILE       output degree distribution of edges/nodes to FILE (must have same num of input files)
      --output FILE         output matrix market file (with/without relabeling by degree)
      -D NUM                specify column either [0](edge) or [1](node) [default: 0]
      --relabel NUM         relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
      --degree NUM          check the percentile above and below degree [NUM] [default: -1]
      --direction DIR       graph relabeling direction - ascending/descending [default: descending]
      -d, --debug           run in debug mode
```

The applications takes options followed by the arguments of the options as inputs. A minimal example takes a hypergraph as input is as follow:
```
$ hystats.exe -f karate.mtx
```


## Supported file format in NWHy

NWHy recogonizes three types of file format:
* [Matrix Market Exchange Formats](https://math.nist.gov/MatrixMarket/formats.html)
* Comma-separated values (CSV) format
* Hypergraph Adjacency format


### Hypergraph Adjacency Format

The hypergraph adjacency format starts with a sequence of offsets
 one for each vertex, followed by a sequence of incident hyperedges
 (the vertex is an incoming member of the hyperedge) ordered by
 vertex, followed by a sequence of offsets one for each hyperedge, and
 finally a sequence of incident vertices (the vertex is an outgoing
 member of the hyperedge) ordered by hyperedge. All vertices,
 hyperedges, and offsets are 0 based and represented in decimal. For a
 graph with *nv* vertices, *mv* incident hyperedges for the vertices,
 *nh* hyperedges, and *mh* incident vertices for the hyperedges, the
 specific format is as follows:

AdjacencyHypergraph  
&lt;nv>		     
&lt;mv>		     
&lt;nh>		     
&lt;mh>		     
&lt;ov0>  
&lt;ov1>  
...  
&lt;ov(nv-1)>  
&lt;ev0>  
&lt;ev1>  
...  
&lt;ev(mv-1)>  
&lt;oh0>  
&lt;oh1>  
...  
&lt;oh(nh-1)>  
&lt;eh0>  
&lt;eh1>  
...  
&lt;eh(mh-1)>  

This file is represented as plain text.

Weighted hypergraphs: the weights are listed as
another sequence following the sequence of neighbors for vertices or
hyperedges file (i.e., after &lt;ev(mv-1)> and &lt;eh(mh-1)>), and the
first line of the file should store the string
"WeightedAdjacencyHypergraph". Note the weighted hypergraphs are not fully tests.

## Use NWHy as a Python Module

NWHy has been used by [HyperNetX(HNX)](https://pnnl.github.io/HyperNetX/build/index.html) as the C++ backend. HNX is a hypergraph library provides classes and methods for modeling the entities and relationships found in complex networks as hypergraphs, the natural models for multi-dimensional network data. See 
[NWHy](https://pypi.org/project/nwhy/) PyPI release page and 
[User Guide](https://pnnl.github.io/HyperNetX/build/nwhy.html) for more information.

### Requirements

* Python &gt;= 3.9
* oneTBB &gt;= 2021



### Install NWHy

#### To install in an Anaconda environment
First create a conda environment.  
```
$ conda create -n <env name> python=3.9
```
Then activate the environment.

```
$ conda activate <env name> 
```
Next install Intel one Threading Building Blocks(TBB) in current conda environment.
```
$ conda install tbb
```
If oneTBB has been installed locally, we can specify TBBROOT
```
$ export TBBROOT=/opt/tbb/
```
Finally install NWHy using pip.

For installation:
```
$ pip install nwhy
```
For upgrade:
```
$ pip install nwhy --upgrade
```
or 
```
$ pip install nwhy -U
```
For quick test whether the environment is ready:
```
$ python -c "import nwhy"
```
If there is no import error, then installation is done.

### Pytest

NWHy uses Pytest framework for the unit test for our Python APIs. The unit test cases are under python/test/ directory.

## Tools in NWHy

The tools are under build/tools/ directory. The following tools are provided in NWHy:

### Converters

**adj2mmio** converts the hypergraph adjacency format to the Matrix Market format.

**mmio2adj** converts the Matrix Market format to the hypergraph adjacency format.

**mmio2mmio** transposes the Matrix Market format file to its transpose hypergraph.

**imdb2mtx** extracts some information from the raw IDMB tsv format and save it to the Matrix Market format.

### Random generator

**hypergen** takes the hyperedge degree distribution and the hypernodes degree distribution of the input hypergraph to generate a new random hypergraph using configuration model. The hypergraph has the same degree distribution of the old hypergraph. Note the generator is not fully tested.

### Other

**hystat** collects some of the basic statistics of a hypergraph.

Please cite:
-----------
Xu T. Liu, Jesun Firoz, Andrew Lumsdaine, Cliff Joslyn, Sinan G. Aksoy, Brenda Praggastis, Assefaw H. Gebremedhin.Parallel Algorithms for Efficient Computation of High-Order Line Graphs of Hypergraphs. 28th IEEE Int’l Conference on High Performance Computing, Data, and Analytics (HiPC 2021), 2021.

Xu T. Liu, Jesun Firoz, Sinan G. Aksoy, Ilya Amburg, Andrew Lumsdaine, Cliff Joslyn, Brenda Praggastis, Assefaw H. Gebremedhin. High-order Line Graphs of Non-uniform Hypergraphs: Algorithms, Applications, and Experimental Analysis. 36th IEEE Int’l Parallel & Distributed Processing Symposium (IPDPS), 2022.
