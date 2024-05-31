FROM rust 
LABEL author=ggutow@andrew.cmu.edu 

#Eigen
WORKDIR /src
COPY eigen-3.4.0.tar.gz .
RUN tar xvf eigen-3.4.0.tar.gz

#CPLEX
WORKDIR /src
COPY ../cplex_studio2211.linux_x86_64.bin .
RUN cplex_studio2211.linux_x86_64.bin

#QSOPT
WORKDIR /src
COPY QSOPT /src/

#Concorde
WORKDIR /src
COPY co031219.tgz .
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
COPY gmp-6.3.0.tar.xz .
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
