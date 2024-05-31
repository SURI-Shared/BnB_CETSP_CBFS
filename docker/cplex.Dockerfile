FROM openjdk:8-jre

WORKDIR /src
COPY cplex_studio2211.linux_x86_64.bin .
COPY cplex_install_response_file.properties .
RUN chmod u+x cplex_studio2211.linux_x86_64.bin
RUN ./cplex_studio2211.linux_x86_64.bin -f cplex_install_response_file.properties