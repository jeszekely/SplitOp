# Split Operator Code
# Author J. Szekely
# Makefile

CC = g++

#MKLROOT = /opt/intel/composer_xe_2013_sp1/mkl
CC = /opt/local/bin/g++-mp-4.8

CFLAGS = -O3 -I$(MKLROOT)/include -Wall -Wno-sign-compare -Wno-unused-function -Werror -std=c++11 -openmp

BOOST_INC = -I/opt/local/include/

OBJ  =   obj/main.o obj/wvfxn.o

LIBS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_sequential -lpthread -lm

HEADS =  src/arrays.hpp src/wvfxn.hpp
BIN  =   SplitOp

RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) -I./src $(OBJ) $(BOOST_INC) $(LIBS)

obj/main.o: src/main.cpp
	$(CC) $(CFLAGS) -c  src/main.cpp -o obj/main.o $(BOOST_INC) -I./src

obj/wvfxn.o: src/wvfxn.cpp
	$(CC) $(CFLAGS) -c  src/wvfxn.cpp -o obj/wvfxn.o $(BOOST_INC) -I./src