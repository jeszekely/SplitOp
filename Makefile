# Split Operator Code
# Author J. Szekely
# Makefile

CC = g++

#MKLROOT = /opt/intel/composer_xe_2013_sp1/mkl
CC = /opt/local/bin/g++-mp-4.8
#CC = clang++ -I/opt/local/include/

CFLAGS = -O3 -I$(MKLROOT)/include -I$(MKLROOT)/include/fftw -std=c++11 -openmp -Wall -Wno-sign-compare -Wno-unused-function -Werror

BOOST_INC = -I/opt/local/include/boost/

OBJ  =   obj/main.o  obj/matrix.o obj/wvfxn.o obj/input_parser.o obj/splitop.o obj/junction.o obj/chebyshev.o obj/surrogate.o

LIBS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_sequential -lpthread -lm

HEADS =  src/arrays.hpp src/wvfxn.hpp src/input_parser.hpp src/splitop.hpp src/junction.hpp src/chebyshev.hpp src/surrogate.hpp src/matrix.hpp
BIN  =   SplitOp

RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) -I./src $(OBJ) $(LIBS)

obj/main.o: src/main.cpp
	$(CC) $(CFLAGS) -c  src/main.cpp -o obj/main.o $(BOOST_INC) -I./src

obj/matrix.o: src/matrix.cpp
	$(CC) $(CFLAGS) -c  src/matrix.cpp -o obj/matrix.o -I./src

obj/wvfxn.o: src/wvfxn.cpp
	$(CC) $(CFLAGS) -c  src/wvfxn.cpp -o obj/wvfxn.o -I./src

obj/input_parser.o: src/input_parser.cpp
	$(CC) $(CFLAGS) -c  src/input_parser.cpp -o obj/input_parser.o $(BOOST_INC) -I./src

obj/splitop.o: src/splitop.cpp
	$(CC) $(CFLAGS) -c  src/splitop.cpp -o obj/splitop.o $(BOOST_INC) -I./src

obj/junction.o: src/junction.cpp
	$(CC) $(CFLAGS) -c  src/junction.cpp -o obj/junction.o $(BOOST_INC) -I./src

obj/chebyshev.o: src/chebyshev.cpp
	$(CC) $(CFLAGS) -c  src/chebyshev.cpp -o obj/chebyshev.o $(BOOST_INC) -I./src

obj/surrogate.o: src/surrogate.cpp
	$(CC) $(CFLAGS) -c  src/surrogate.cpp -o obj/surrogate.o $(BOOST_INC) -I./src
