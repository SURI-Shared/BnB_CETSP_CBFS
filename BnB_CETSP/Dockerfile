FROM rust 
LABEL author=ggutow@andrew.cmu.edu 

#system dependencies
RUN apt-get install eigen

#CPLEX
WORKDIR /src
COPY cplex_installer.bin .
RUN cplex_installer.bin

#QSOPT
WORKDIR /src/QSOPT
RUN wget https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/ubuntu/qsopt.h
RUN wget https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/ubuntu/qsopt.a

#Concorde
WORKDIR /src
RUN wget https://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz
RUN gunzip co031219.tgz
RUN tar xvf co031219.tar
WORKDIR /src/concorde/concorde_build
RUN ../configure --with-qsopt=/src/QSOPT
RUN make

#Clarabel
WORKDIR /src
COPY ../Clarabel.cpp .

WORKDIR /src/Clarabel.cpp/build
RUN cmake ..
RUN make

#GMP BigNum
WORKDIR /src
RUN wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
RUN tar xz gmp-6.3.0.tar.xz
WORKDIR /src/gmp-6.3.0
RUN configure
RUN make
RUN make install

#branch and bound
WORKDIR /src
COPY ../BnB_CETSP_CBFS .

WORKDIR /src/BnB_CETSP_CBFS
RUN make exeCVXHULL
RUN make clarabel_redundant
RUN make clarabel_dropin
RUN make clarabel_recycling
RUN make clarabel_warmstart
