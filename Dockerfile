FROM ggutow/cplex2211:latest as cplex
FROM ggutow/concorde:latest as concorde
FROM ggutow/gmp_bignum:latest as gmp_bignum
FROM ggutow/clarabel.cpp:latest as clarabel
LABEL author=ggutow@andrew.cmu.edu

#copy from dependencies
WORKDIR /src
COPY --from=cplex /opt/ibm/ILOG/CPLEX_Studio2211/ /opt/ibm/ILOG/CPLEX_Studio2211/
COPY --from=gmp_bignum /src/gmp-6.3.0 /src/gmp-6.3.0
WORKDIR /src/gmp-6.3.0
RUN make install
COPY --from=concorde /src/QSOPT /src/QSOPT
COPY --from=concorde /src/concorde/ /src/concorde

#branch and bound
WORKDIR /src/BnB_CETSP_CBFS
COPY . .
WORKDIR /src/BnB_CETSP_CBFS/BnB_CETSP/obj/exec
WORKDIR /src/BnB_CETSP_CBFS/BnB_CETSP/obj/test
WORKDIR /src/BnB_CETSP_CBFS/BnB_CETSP/obj/src
WORKDIR /src/BnB_CETSP_CBFS/BnB_CETSP/run
WORKDIR /src/BnB_CETSP_CBFS/docker
RUN make exeCVXHULL
RUN make clarabel_redundant
RUN make clarabel_dropin
RUN make clarabel_recycling
RUN make clarabel_warmstart

RUN apt-get install -qq pip
RUN pip install --break-system-packages --no-cache-dir git-python