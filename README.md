# BnB_CETSP

A branch-and-bound algorithm for the CETSP based the [original work by Walton Pereira Coutinho](https://github.com/waltonpcoutinho/BnB_CETSP) and [Wenda Zhang](https://github.com/UranusR/BnB_CETSP_CBFS)

This code accompanies the paper "Efficient Second-Order Cone Programming for the Close Enough Traveling Salesman Problem," submitted to ICRA 2025.

## Dependencies

### Clarabel.rs and Clarabel.cpp wrapper
This code relies on a fork of the Clarabel cone programming solver. This fork is available at https://github.com/SURI-Shared/Clarabel.cpp

Clone the Clarabel.cpp repository (including the submodule, which contains the actual modified Clarabel solver!):
```
git clone --recurse-submodules git@github.com:SURI-Shared/Clarabel.cpp.git
```

### SCS
This code relies on a fork of the SCS cone programming solver. This fork is available at [https://github.com/SURI-Shared/scs/tree/tridiagonal](https://github.com/SURI-Shared/scs/tree/tridiagonal)

Clone the SCS repository (make sure to switch to the tridiagonal branch):
```
git clone git@github.com:SURI-Shared/scs
cd scs
git checkout tridiagonal
```

### Eigen, QSOPT, Concorde, and GMP BigNum
tar files and archives of these dependencies can be obtained by running download_dependencies.sh

### CPLEX
CPLEX is used only by the baseline approach. You will need to obtain a CPLEX install to build the rest of the code, however.

### LAPACK and LAPACKE
SCS depends on both LAPACK for linear algebra and LAPACKE (the C interface to LAPACK).
```
sudo apt-get install liblapack-dev liblapacke-dev
```

## Building an executable

### Locally
1. Build concorde following the instructions at https://www.math.uwaterloo.ca/tsp/concorde/DOC/README.html
2. Build gmp_bignum (this is SLOW)
3. Build Clarabel.cpp using cmake (We recommend doing a Release build)
4. Build SCS using cmake (multiple shared object libraries are produced, one for each choice of linear system solver. We use both scsdir and scstridir).
```
cd scs
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
5. Install CPLEX
6. Update BnB_CETSP_CBFS/BnB_CETSP/makefile (NOT BnB_CETSP_CBFS/docker/makefile) to reflect the paths (everything in the section ###directories)

For help compiling things, look at the various docker files for hints.

The executable targets are:
* exeCVXHULL: As in (https://pubsonline.informs.org/doi/abs/10.1287/ijoc.2016.0711), with some additional logging
* clarabel_redundant: Replaces CPLEX with Clarabel, but no other changes from exeCVXHULL
* clarabel_reduce: Uses Clarabel to solve RSOCP not MSOCP (see the paper)
* clarabel_reuse: as clarabel_reduce, but reuses the Clarabel solvers between nodes in the B&B tree
* clarabel_recycle: as clarabel_reuse, but warm starts Clarabel by solving a small SOCP.
* scs_reduce: Uses SCS to solve RSOCP not MSOCP (see the paper)
* scs_reuse: as scs_reduce, but reuses the SCS solvers between nodes in the B&B tree (and also the matrix factorization)
* scs_recycle: as scs_reuse, but warm starts SCS by solving a small SOCP.
* scs_tridiag_reduce: Uses SCS with a custom tridiagonal solver to solve RSOCP not MSOCP (see the paper)
* scs_tridiag_reuse: as scs_tridiag_reduce, but reuses the SCS solvers between nodes in the B&B tree (and also the matrix factorization)
* scs_tridiag_recycle: as scs_tridiag_reuse, but warm starts SCS by solving a small SOCP.


### Docker
The docker folder contains Dockerfiles for building images that can then be derived from to build the actual codebase.
1. build concorde.Dockerfile
```
cd BnB_CETSP_CBFS/docker
docker build -t concorde -f concorde.Dockerfile .
```
2. build gmp_bignum.Dockerfile (this is SLOW)
```
docker build -t gmp_bignum -f gmp_bignum.Dockerfile .
```
3. build clarabelcpp's Docker image using the Dockerfile from https://github.com/SURI-Shared/Clarabel.cpp
```
cd Clarabel.cpp
docker build . -t clarabel.cpp
```
4. build scs's Docker image using the Dockerfile from [https://github.com/SURI-Shared/scs/tree/tridiagonal](https://github.com/SURI-Shared/scs/tree/tridiagonal)
```
cd scs
docker build . -t scs
```
5. build cplex's Docker image using cplex.Dockerfile. You will need to place cplex_install_response_file.properties into the same folder as cplex_studio2211.linux_x86_64.bin
```
cd BnB_CETSP_CBFS/docker
cp cplex_install_response_file.properties /folder/with/cplex_studio2211.linux_x86_64.bin
docker build -t cplex2211 -f cplex.Dockerfile /folder/with/cplex_studio2211.linux_x86_64.bin
```
6. Update the image tags in BnB_CETSP_CBFS/Dockerfile and build
```
cd BnB_CETSP_CBFS
docker build -t r3cetsp .
```

To run, attach interactively to the container produced by step 6. 
```
docker run -it r3cetsp
```

Result files can be extracted by using a [bind-mount ](https://docs.docker.com/storage/bind-mounts/)

To build an executable, run the following command from the project root directory "BnB_CETSP":
```
make
```

## Calling command:

Scripts to run multiple instances in parallel are provided in BnB_CETSP/exec

To process all instances in a folder in parallel using one of the executables:
```
python exec/run_instances.py [n_processes] [executable] [folder of instances] [folder to save log files] [OPTIONS] [OVERLAP FACTOR] [TIME LIMIT] [BRANCHING RULE] [BRANCHING Strategy] [ROOT SELECTION] [S.B SIZE]
```

To process the 152 problem instances from https://pubsonline.informs.org/doi/abs/10.1287/ijoc.2016.0711:
```
python exec/run_all.py [n_processes] [executable] [folder to save log files]
```

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


