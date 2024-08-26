#makefile

SYSTEM = x86-64_linux
LIBFORMAT = static_pic 

####diretorios com as libs do cplex
INCLUDE 	  = /home/ggutow/eclipse-workspace/BnB_CETSP_CBFS/BnB_CETSP/src
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio2211/cplex
QSOPTDIR	  = /home/ggutow/eclipse-workspace/concorde/QSOPT
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio2211/concert
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCORDEINCDIR = /home/ggutow/eclipse-workspace/concorde/INCLUDE 
CONCORDELIBDIR = /home/ggutow/eclipse-workspace/concorde/install
CLARABEL_INCLUDE = /home/ggutow/eclipse-workspace/Clarabel.cpp/include
CLARABEL_LIBDIR = /home/ggutow/eclipse-workspace/Clarabel.cpp/rust_wrapper/target/release/
EIGEN_INCLUDE = /home/ggutow/miniconda3/envs/gobnilp/include/eigen3/
SCS_INCLUDE = /home/ggutow/eclipse-workspace/scs-python/scs_source/include/
SCS_LIBDIR = /home/ggutow/eclipse-workspace/scs-python/scs_source/build/lib/

#### define o compilador
CPPC = g++
#############################

#### compilation options and includes
CCOPT = $(BITS_OPTION) -O3 -g -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++14
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(INCLUDE) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(CONCORDEINCDIR) -I${CLARABEL_INCLUDE} -I${SCS_INCLUDE} -I${EIGEN_INCLUDE}
#############################

#### flags do linker
CCLNFLAGS = -L$(CONCORDELIBDIR) -lconcorde -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lgmpxx -lgmp -ldl -lm -lpthread -L${QSOPTDIR} -lqsopt
CLARABELLINFLAGS = -L${CLARABEL_LIBDIR} -lclarabel_c
SCSDIRECTLINFLAGS = -L${SCS_LIBDIR} -lscsdir#link to generic direct SCS solver
#############################

#### diretorios com os source files e com os objs files
SRCDIR = src
OBJDIR = obj
TESTDIR = test
EXECDIR = exec
#############################

#### list of all srcs and all objs
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/src/%.o, $(SRCS))
CLARABEL_SRCS = $(wildcard $(SRCDIR)/clarabel_interface/*.cpp)
CLARABEL_OBJS = $(patsubst $(SRCDIR)/clarabel_interface/%.cpp, $(OBJDIR)/src/%.o, $(CLARABEL_SRCS))
SCS_SRCS = $(wildcard $(SRCDIR)/scs_interface/*.cpp)
SCS_OBJS = $(patsubst $(SRCDIR)/scs_interface/%.cpp, $(OBJDIR)/src/%.o, $(SCS_SRCS))
TEST_SRCS = $(wildcard $(TESTDIR)/*.cpp)
TEST_OBJS = $(patsubst $(TESTDIR)/%.cpp, $(OBJDIR)/test/%.o, $(TEST_SRCS))
EXEC_SRCS =  $(wildcard $(EXECDIR)/*.cpp)
TEST_OBJS = $(patsubst $(EXECDIR)/%.cpp, $(OBJDIR)/exec/%.o, $(EXEC_SRCS))
#############################

#### regra principal, gera o executavel
exeCVXHULL: $(OBJDIR)/exec/main.o $(OBJS) $(CLARABEL_OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(CLARABEL_OBJS) $(OBJDIR)/exec/main.o -Wl,-rpath,${CLARABEL_LIBDIR} -o run/$@ $(CCLNFLAGS) $(CLARABELLINFLAGS)
############################

clarabel_dropin: $(OBJDIR)/exec/clarabel_dropin.o $(OBJS) $(CLARABEL_OBJS) $(SCS_OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(SCS_OBJS) $(CLARABEL_OBJS) $(OBJDIR)/exec/clarabel_dropin.o -Wl,-rpath,${SCS_LIBDIR} -Wl,-rpath,${CLARABEL_LIBDIR} -o run/$@ $(CCLNFLAGS) $(CLARABELLINFLAGS) $(SCSDIRECTLINFLAGS)

clarabel_redundant: $(OBJDIR)/exec/clarabel_redundant.o $(OBJS) $(CLARABEL_OBJS) $(SCS_OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(SCS_OBJS) $(CLARABEL_OBJS) $(OBJDIR)/exec/clarabel_redundant.o -Wl,-rpath,${SCS_LIBDIR} -Wl,-rpath,${CLARABEL_LIBDIR} -o run/$@ $(CCLNFLAGS) $(CLARABELLINFLAGS) $(SCSDIRECTLINFLAGS)

clarabel_recycling: $(OBJDIR)/exec/clarabel_recycling.o $(OBJS) $(CLARABEL_OBJS) $(SCS_OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(SCS_OBJS) $(CLARABEL_OBJS) $(OBJDIR)/exec/clarabel_recycling.o -Wl,-rpath,${SCS_LIBDIR} -Wl,-rpath,${CLARABEL_LIBDIR} -o run/$@ $(CCLNFLAGS) $(CLARABELLINFLAGS) $(SCSDIRECTLINFLAGS)

clarabel_warmstart: $(OBJDIR)/exec/clarabel_warmstart.o $(OBJS) $(CLARABEL_OBJS) $(SCS_OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(SCS_OBJS) $(CLARABEL_OBJS) $(OBJDIR)/exec/clarabel_warmstart.o -Wl,-rpath,${SCS_LIBDIR} -Wl,-rpath,${CLARABEL_LIBDIR} -o run/$@ $(CCLNFLAGS) $(CLARABELLINFLAGS) $(SCSDIRECTLINFLAGS)

test/%: $(OBJDIR)/test/%.o $(OBJS) $(CLARABEL_OBJS) $(SCS_OBJS)
	@echo "\033[31m \nLinking test: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(SCS_OBJS) $(CLARABEL_OBJS) -Wl,-rpath,${SCS_LIBDIR} -Wl,-rpath,${CLARABEL_LIBDIR} $< -o run/$@ $(CCLNFLAGS) $(CLARABELLINFLAGS) $(SCSDIRECTLINFLAGS)

scs_dropin: $(OBJDIR)/exec/scs_dropin.o $(OBJS) $(SCS_OBJS) $(CLARABEL_OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(SCS_OBJS) $(CLARABEL_OBJS) $(OBJDIR)/exec/scs_dropin.o -Wl,-rpath,${SCS_LIBDIR} -Wl,-rpath,${CLARABEL_LIBDIR} -o run/$@ $(CCLNFLAGS) $(SCSDIRECTLINFLAGS) $(CLARABELLINFLAGS)

scs_recycling: $(OBJDIR)/exec/scs_recycling.o $(OBJS) $(CLARABEL_OBJS) $(SCS_OBJS)
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) $(SCS_OBJS) $(CLARABEL_OBJS) $(OBJDIR)/exec/scs_recycling.o -Wl,-rpath,${SCS_LIBDIR} -Wl,-rpath,${CLARABEL_LIBDIR} -o run/$@ $(CCLNFLAGS) $(CLARABELLINFLAGS) $(SCSDIRECTLINFLAGS)

#includes dependency files
-include $(OBJS:.o=.d)

#rule for each object file: compiles and generates object file dependencies file
#each object file depends on the .c and the headers (header information is in the dependency file generated by the compiler)
$(OBJDIR)/src/%.o: $(SRCDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC)  -MM -I$(INCLUDE) $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

$(OBJDIR)/src/%.o: $(SRCDIR)/clarabel_interface/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC)  -MM -I$(INCLUDE) $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

$(OBJDIR)/src/%.o: $(SRCDIR)/scs_interface/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC)  -MM -I$(INCLUDE) $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

$(OBJDIR)/test/%.o: $(TESTDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC)  -MM -I$(INCLUDE) $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

$(OBJDIR)/exec/%.o: $(EXECDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC)  -MM -I$(INCLUDE) $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp
#delete objetos e arquivos de dependencia
clean:
	@echo "\033[31mcleaning obj directory \033[0m"
	@rm -f $(OBJDIR)/src/**.o $(OBJDIR)/src/**.d $(OBJDIR)/test/**.o $(OBJDIR)/test/**.d $(OBJDIR)/exec/**.o $(OBJDIR)/exec/**.d

rebuild: clean exeCVXHULL

