FROM rust 
LABEL author=ggutow@andrew.cmu.edu 

#branch and bound
WORKDIR /src
COPY . .

#Eigen
WORKDIR /src/BnB_CETSP_CBFS
RUN tar xf /src/BnB_CETSP_CBFS/eigen-3.4.0.tar.gz

#CPLEX
WORKDIR /src
RUN chmod u+x cplex_studio2211.linux_x86_64.bin
RUN ./cplex_studio2211.linux_x86_64.bin -f cplex_install_response_file.properties

#Concorde
WORKDIR /src/BnB_CETSP_CBFS
RUN gunzip co031219.tgz
RUN tar xf co031219.tar
WORKDIR /src/BnB_CETSP_CBFS/concorde/concorde_build
RUN /src/BnB_CETSP_CBFS/concorde/configure --with-qsopt=/src/BnB_CETSP_CBFS/QSOPT
RUN make

#Clarabel
WORKDIR /src/Clarabel.cpp/build
RUN cmake /src/Clarabel.cpp
RUN make

#GMP BigNum
WORKDIR /src/BnB_CETSP_CBFS
RUN tar xz /src/BnB_CETSP_CBFS/gmp-6.3.0.tar.xz
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
