# Droplet Makefile

##################
# Makefile links #
##################

# http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
# http://mrbook.org/blog/tutorials/make/
# http://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules

#########
# Flags #
#########

# Project directories
BASE=.
ROOT=$(BASE)/root
TEST=$(BASE)/test
BIN=$(BASE)/bin
LIB=$(BASE)/lib

YAML_PREFIX=$(BASE)
YAML_INC=$(YAML_PREFIX)/include
YAML_INC=$(YAML_PREFIX)/lib

# C++ Compiler
CC=g++
# Root flags
RF=`root-config --glibs --cflags` #-lMathMore -lMinuit

CFLAGS=$(RF) -lyaml

CLARGS=some/path/inFile.txt some/path/outFile.txt $(ROOT)/droplet_analysis_options.txt

###############
# EXECUTABLES #
###############

$(BIN)/droplet_analysis.out: $(LIB)/droplet_analysis.o $(LIB)/Utils.o $(LIB)/MDBase.o $(LIB)/Droplet.o $(LIB)/Parameters.o $(LIB)/Atoms.o $(LIB)/Readers.o
	$(CC) -o $@ $^ $(CFLAGS)


#############
# UTILITIES #
#############

all: $(LIB)/Utils.o $(LIB)/MDBase.o $(LIB)/Parameters.o $(LIB)/Atoms.o $(LIB)/Readers.o $(LIB)/Droplet.o

clean:
	rm $(LIB)/*.o -f
	rm $(TEST)/*.out -f

run: $(BIN)/droplet_analysis.out
	$(BIN)/droplet_analysis.out

#########
# TESTS #
#########

test: $(TEST)/test_root.out $(TEST)/test_parameters.out
	$(TEST)/test_root.out
	$(TEST)/test_parameters.out $(CLARGS)

$(TEST)/test_inputstream.out: $(ROOT)/test_inputstream.cpp $(ROOT)/Utils.cpp
	$(CC) -o $@ $^ $(CFLAGS)

$(TEST)/test_root.out: $(ROOT)/test_root.cpp
	$(CC) -o $@ $^ $(CFLAGS)

$(TEST)/test_parameters.out: $(LIB)/test_parameters.o $(LIB)/Parameters.o $(LIB)/Utils.o
	$(CC) -o $@ $^ $(CFLAGS)

$(LIB)/test_parameters.o: $(ROOT)/test_parameters.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


#############
# LIBRARIES #
#############

$(LIB)/Utils.o: $(ROOT)/Utils.cpp $(ROOT)/Utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

$(LIB)/MDBase.o: $(ROOT)/MDBase.cpp $(ROOT)/MDBase.h
	$(CC) -o $@ -c $< $(CFLAGS)

$(LIB)/Parameters.o: $(ROOT)/Parameters.cpp $(ROOT)/Parameters.h
	$(CC) -o $@ -c $< $(CFLAGS)

$(LIB)/Atoms.o: $(ROOT)/Atoms.cpp $(ROOT)/Atoms.h
	$(CC) -o $@ -c $< $(CFLAGS)

$(LIB)/Readers.o: $(ROOT)/Readers.cpp $(ROOT)/Readers.h
	$(CC) -o $@ -c $< $(CFLAGS)

$(LIB)/Time.o: $(ROOT)/Time.cpp $(ROOT)/Time.h
	$(CC) -o $@ -c $< $(CFLAGS)

$(LIB)/Droplet.o: $(ROOT)/Droplet.cpp $(ROOT)/Droplet.h
	$(CC) -o $@ -c $< $(CFLAGS)

$(LIB)/droplet_analysis.o: $(ROOT)/droplet_analysis.cpp
	$(CC) -o $@ -c $< $(CFLAGS)
