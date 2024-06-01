FROM rust

LABEL author=ggutow@andrew.cmu.edu

#GMP BigNum
WORKDIR /src
RUN wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
RUN tar xf /src/gmp-6.3.0.tar.xz
WORKDIR /src/gmp-6.3.0
RUN /src/gmp-6.3.0/configure
RUN make
RUN make install