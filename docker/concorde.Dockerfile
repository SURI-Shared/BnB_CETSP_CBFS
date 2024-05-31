FROM rust

WORKDIR /src/QSOPT
RUN wget https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/ubuntu/qsopt.h
RUN wget https://www.math.uwaterloo.ca/~bico/qsopt/downloads/codes/ubuntu/qsopt.a
WORKDIR /src
RUN wget https://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz
RUN gunzip co031219.tgz
RUN tar xf co031219.tar
WORKDIR /src/concorde/concorde_build
RUN /src/concorde/configure --with-qsopt=/src/QSOPT --host=x86_64-unknown-linux-gnu
RUN make