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


# C++ Compiler
CXX=g++
# Root flags
ROOTFLAGS=`root-config --glibs --cflags` #-lMathMore -lMinuit
CFLAGS=$(ROOTFLAGS) -lyaml

# Project directories
BASE=..
INC=.
SRC=$(BASE)/src/cpp
BIN=$(BASE)/bin
LIB=$(BASE)/lib
YAML=$(BASE)/libyaml
TEST=$(BASE)/test

# All C++ source, header, and object files
cpp_src=$(wildcard $(SRC)/*.cpp)
cpp_hdr=$(cpp_src:$(SRC)/%.cpp=$(INC)/%.h)
cpp_obj=$(cpp_src:$(SRC)/%.cpp=$(INC)/%.o)

#############
# UTILITIES #
#############

.PHONY: all clean run echo

all: $(BIN)/droplet_analysis.out

clean:
	rm $(LIB)/*.o -f
	rm $(TEST)/*.out -f

run: $(BIN)/droplet_analysis.out
	$(BIN)/droplet_analysis.out

echo:
	echo $(cpp_obj)

#############
# LIBRARIES #
#############

$(cpp_obj): %.o: $(SRC)/%.cpp $(INC)/%.h
	$(CXX) -o $@ -c $< $(CFLAGS)

Parameters.o: $(LIB)/libyaml.so

###################
# Compile libYAML #
###################

$(LIB)/libyaml.so:
	cd $(YAML) && ./bootstrap
	cd $(YAML) && ./configure --prefix=$(PWD)/..
	cd $(YAML) && make
	cd $(YAML) && make install

###############
# EXECUTABLES #
###############

$(BIN)/droplet_analysis.out: $(cpp_obj)
	$(CXX) -o $@ $^ $(CFLAGS)

#########
# TESTS #
#########

CLARGS=some/path/inFile.txt some/path/outFile.txt $(ROOT)/droplet_analysis_options.txt

test: $(TEST)/test_root.out $(TEST)/test_parameters.out
	$(TEST)/test_root.out
	$(TEST)/test_parameters.out $(CLARGS)

$(TEST)/test_inputstream.out: $(ROOT)/test_inputstream.cpp $(ROOT)/Utils.cpp
	$(CXX) -o $@ $^ $(CFLAGS)

$(TEST)/test_root.out: $(ROOT)/test_root.cpp
	$(CXX) -o $@ $^ $(CFLAGS)

$(TEST)/test_parameters.out: $(LIB)/test_parameters.o $(LIB)/Parameters.o $(LIB)/Utils.o
	$(CXX) -o $@ $^ $(CFLAGS)

$(LIB)/test_parameters.o: $(ROOT)/test_parameters.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)


