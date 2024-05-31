FROM rust 
LABEL author=ggutow@andrew.cmu.edu 

#branch and bound
WORKDIR /src
COPY BnB_CETSP_CBFS .

#Eigen
WORKDIR /src/BnB_CETSP_CBFS
RUN tar xf eigen-3.4.0.tar.gz

#CPLEX
WORKDIR /src
COPY cplex_studio2211.linux_x86_64.bin .
RUN cplex_studio2211.linux_x86_64.bin

#Concorde
WORKDIR /src/BnB_CETSP_CBFS
RUN gunzip co031219.tgz
RUN tar xf co031219.tar
WORKDIR /src/BnB_CETSP_CBFS/concorde/concorde_build
RUN ../configure --with-qsopt=/src/BnB_CETSP_CBFS/QSOPT
RUN make

#Clarabel
WORKDIR /src
COPY Clarabel.cpp .

WORKDIR /src/Clarabel.cpp/build
RUN cmake ..
RUN make

#GMP BigNum
WORKDIR /src/BnB_CETSP_CBFS
RUN tar xz gmp-6.3.0.tar.xz
WORKDIR /src/gmp-6.3.0
RUN configure
RUN make
RUN make install

WORKDIR /src/BnB_CETSP_CBFS
RUN make exeCVXHULL
RUN make clarabel_redundant
RUN make clarabel_dropin
RUN make clarabel_recycling
RUN make clarabel_warmstart
