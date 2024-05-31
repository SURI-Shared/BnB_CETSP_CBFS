FROM rust

RUN apt-get install cmake

WORKDIR /src
COPY Clarabel.cpp .
WORKDIR /src/Clarabel.cpp/build
RUN cmake /src/Clarabel.cpp
RUN make