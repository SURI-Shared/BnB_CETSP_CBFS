FROM ggutow/cplex2211:latest AS cplex
FROM ggutow/concorde:latest AS concorde
FROM ggutow/gmp_bignum:latest AS gmp_bignum
FROM ggutow/scs:latest AS scs
FROM ggutow/clarabel.cpp:latest AS clarabel
LABEL author=ggutow@andrew.cmu.edu

#copy from dependencies
WORKDIR /src
COPY --from=cplex /opt/ibm/ILOG/CPLEX_Studio2211/ /opt/ibm/ILOG/CPLEX_Studio2211/
COPY --from=gmp_bignum /src/gmp-6.3.0 /src/gmp-6.3.0
WORKDIR /src/gmp-6.3.0
RUN make install
COPY --from=concorde /src/QSOPT /src/QSOPT
COPY --from=concorde /src/concorde/ /src/concorde
COPY --from=scs /src/scs /src/scs
RUN apt-get install -qq liblapack-dev liblapacke-dev

#branch and bound
WORKDIR /src/BnB_CETSP_CBFS
COPY docker docker
COPY BnB_CETSP/2D BnB_CETSP/2D 
COPY BnB_CETSP/3D BnB_CETSP/3D
COPY BnB_CETSP/Behdani BnB_CETSP/Behdani
COPY BnB_CETSP/medium_2D_Behdani_CETSPs BnB_CETSP/medium_2D_Behdani_CETSPs
COPY BnB_CETSP/exec BnB_CETSP/exec
COPY BnB_CETSP/RND BnB_CETSP/RND
COPY BnB_CETSP/src BnB_CETSP/src
COPY BnB_CETSP/test BnB_CETSP/test
COPY BnB_CETSP/Upper_Bounds BnB_CETSP/Upper_Bounds
WORKDIR /src/BnB_CETSP_CBFS/BnB_CETSP/obj/exec
WORKDIR /src/BnB_CETSP_CBFS/BnB_CETSP/obj/test
WORKDIR /src/BnB_CETSP_CBFS/BnB_CETSP/obj/src
WORKDIR /src/BnB_CETSP_CBFS/BnB_CETSP/run
WORKDIR /src/BnB_CETSP_CBFS/docker
RUN make exeCVXHULL
RUN make clarabel_redundant
RUN make clarabel_reduce
RUN make clarabel_reuse
RUN make clarabel_recycle

RUN apt-get install -qq pip
RUN pip install --break-system-packages --no-cache-dir git-python