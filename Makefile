# Split Operator Code, v 2.0
# First version in C++
# Josh Szekely, May 2013
# Makefile

CC = icpc
CFLAGS = -g -O2 -debug -parallel -openmp -DMKL_ILP64 -I$(MKLROOT)/include

FFTW3_DIR = /projects/p20191/jes060/Codes/fftw-3.3.3_c++/
FFTW_LIB = -L$(FFTW3_DIR)

FFTW_INC = -I$(FFTW3_DIR)include
SPLITOP_INC = -I/projects/p20191/jes060/Codes/ImprovedSplitOp/

INI_LIB = -L/projects/p20191/jes060/Codes/iniparser_c++/
INI_INC = -I/projects/p20191/jes060/Codes/iniparser_c++/src/

BOOST_INC = -I/software/boost/1.41.0/

OBJ  =   obj/main.o  obj/array_structs.o  obj/numerics.o obj/junction_hamiltonian.o obj/imag_time_prop.o obj/quantum_dynamics.o

LIBS =   -L$(FFTW3_DIR)lib -lfftw3_threads -lfftw3  -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -liniparser

HEADS =  include/array_structs.h include/numerics.h include/junction_hamiltonian.h include/imag_time_prop.h include/quantum_dynamics.h
BIN  =   SplitOpPropagator

RM = rm  -f

.PHONY: all all-before all-after clean clean-custom
all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS)  -o $(BIN) -I./include $(FFTW_INC) $(INI_INC) $(SPLITOP_INC) $(OBJ) $(FFTW_LIB) $(INI_LIB) $(LIBS)

obj/main.o: src/main.cpp
	$(CC) $(CFLAGS) -c  src/main.cpp  $(INI_LIB) -o obj/main.o $(BOOST_INC) $(INI_INC) -I./include $(FFTW_INC) $(SPLITOP_INC)

obj/array_structs.o: src/array_structs.cpp
	$(CC) $(CFLAGS)  -c src/array_structs.cpp  -o obj/array_structs.o  -I./include $(FFTW_INC) $(SPLITOP_INC)

obj/numerics.o: src/numerics.cpp
	$(CC) $(CFLAGS) -c src/numerics.cpp  -o obj/numerics.o  -I./include $(FFTW_INC) $(SPLITOP_INC)

obj/junction_hamiltonian.o: src/junction_hamiltonian.cpp
	$(CC) $(CFLAGS) -c src/junction_hamiltonian.cpp  -o obj/junction_hamiltonian.o  -I./include $(FFTW_INC) $(SPLITOP_INC)

obj/imag_time_prop.o: src/imag_time_prop.cpp
	$(CC) $(CFLAGS) -c src/imag_time_prop.cpp -o obj/imag_time_prop.o -I./include $(BOOST_INC) $(FFTW_INC) $(SPLITOP_INC)

obj/quantum_dynamics.o: src/quantum_dynamics.cpp
	$(CC) $(CFLAGS) -c src/quantum_dynamics.cpp -o obj/quantum_dynamics.o -I./include $(FFTW_INC) $(FFTW_INC)  $(SPLITOP_INC)
