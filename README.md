# BnB_CETSP

A branch-and-bound algorithm for the CETSP based the [original work by Walton Pereira Coutinho](https://github.com/waltonpcoutinho/BnB_CETSP) and [Wenda Zhang](https://github.com/UranusR/BnB_CETSP_CBFS)

This code accompanies the paper "Reduce, Reuse, Recycle: Efficient Second-Order Cone Programming for the CLose Enough Traveling Salesman Problem," submitted to WAFR 2024.

## Dependencies

### Clarabel.rs and Clarabel.cpp wrapper
This code relies on a fork of the Clarabel cone programming solver. This fork is available at https://github.com/SURI-Shared/Clarabel.cpp

Clone the Clarabel.cpp repository (including the submodule, which contains the actual modified Clarabel solver!):
```
git clone --recurse-submodules git@github.com:SURI-Shared/Clarabel.cpp.git
```

### Eigen, QSOPT, Concorde, and GMP BigNum
tar files and archives of these dependencies can be obtained by running download_dependencies.sh

### CPLEX
CPLEX is used only by the baseline approach. You will need to obtain a CPLEX install to build the rest of the code, however.

## Building an executable

### Locally
1. Build concorde following the instructions at https://www.math.uwaterloo.ca/tsp/concorde/DOC/README.html
2. Build gmp_bignum (this is SLOW)
3. Build Clarabel.cpp using cmake (We recommend doing a Release build)
4. Install CPLEX
5. Update BnB_CETSP_CBFS/BnB_CETSP/makefile (NOT BnB_CETSP_CBFS/docker/makefile) to reflect the installation paths
For help compiling things, look at the various docker files for hints.

The executable targets are:
* exeCVXHULL: As in (https://pubsonline.informs.org/doi/abs/10.1287/ijoc.2016.0711), with some additional logging
* clarabel_redundant: Replaces CPLEX with Clarabel, but no other changes from exeCVXHULL
* clarabel_reduce: Uses Clarabel to solve RSOCP not MSOCP (see the paper)
* clarabel_reuse: as clarabel_reduce, but reuses the Clarabel solvers between nodes in the B&B tree
* clarabel_recycle: as clarabel_reuse, but warm starts Clarabel by solving a small SOCP.


### Docker
The docker folder contains Dockerfiles for building images that can then be derived from to build the actual codebase.
1. build concorde.Dockerfile
2. build gmp_bignum.Dockerfile (this is SLOW)
3. build clarabelcpp's Docker image using the Dockerfile from https://github.com/SURI-Shared/Clarabel.cpp
4. Update the image tags in BnB_CETSP_CBFS/Dockerfile and build

To run, attach interactively to the running container produced by step 4. Result files can be extracted by using a [bind-mount ](https://docs.docker.com/storage/bind-mounts/)

To build an executable, run the following command from the project root directory "BnB_CETSP":
```
make
```
Please, note that:
 * You will need to edit the makefile in order to fix the correct path to your CPLEX installation.
 * You will need to install the "GMP MP Bignum Library" (https://gmplib.org/) in order to run this code.
 * You will need to install Concorde (http://www.math.uwaterloo.ca/tsp/concorde.html).

## Calling command:

Python scripts to run multiple instances in parallel are provided in BnB_CETSP/exec

A single executable can be run using the following command structure
```
./exeCVXHULL [Path to the instance] [OPTIONS] [OVERLAP FACTOR] [TIME LIMIT] [BRANCHING RULE] [BRANCHING STRATEGY] [ROOT SELECTION] [S.B SIZE] [CONTOUR SELECTION] [LAH SELECTION] [LAH TIME LIMIT] [FSI SELECTION]
```
Where:

 * OPTIONS = "2D" or "3D"
 * OVERLAP FACTOR = One of the following values "{0.1, 0.25, 0.5, 1.0, 1.5}"
 * TIME LIMIT = Any integer value
 * BRANCHING RULE = One of the following values "{V1,SB}"
 * BRANCHING Strategy = One of the following values "{DFS,BFS,BeFS,CBFS}"
 * ROOT SELECTION = One of the following values "{1,2,3}"
 * S.B SIZE = Any integer value
 * CONTOUR SELECTION = One of the following values "{1, 2}"
 * LAH SELECTION = One of the following values "{0, 1}"
 * LAH TIME LIMIT = Any positive real number, but it is suggested to use a number around 1.0.
 * FSI SELECTION = One of the following values "{0, 1}"

Not all combinations of parameters will work for all instances. Check the paper for details (https://pubsonline.informs.org/doi/abs/10.1287/ijoc.2016.0711).
Note that, if not running exeCVXHULL, only V1 branching, BeFS, root selection 1, S.B Size of 1 should be used.

Examples:

Running a 2D instance:
```
./exeCVXHULL 2D/rotatingDiamonds1.txt 2D 1.0 3600 V1 BeFS 1 1
```

Running a 3D instance:
```
./exeCVXHULL 3D/d493.txt 3D 0.5 3600 V1 BeFS 1 1
```

## License

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details


