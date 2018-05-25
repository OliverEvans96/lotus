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
BASE=..
INC=.
SRC=$(BASE)/src/cpp
BIN=$(BASE)/bin
LIB=$(PWD)/../lib
YAML=$(BASE)/libyaml
TEST=$(BASE)/test/src

# C++ Compiler
CXX=g++
# Root flags
ROOTFLAGS=`root-config --glibs --cflags` -lMathMore -lMinuit
# Specify rpath for custom shared library path
CFLAGS=$(ROOTFLAGS) -I$(INC) -L$(LIB) -lyaml -Wl,-rpath=$(LIB)

# Main executable
main_src=$(SRC)/droplet_analysis.cpp
main_bin=$(main_src:$(SRC)/%.cpp=%)

# All other C++ source, header, and object files
cpp_src=$(filter-out $(main_src), $(wildcard $(SRC)/*.cpp))
cpp_hdr=$(cpp_src:$(SRC)/%.cpp=$(INC)/%.h)
cpp_obj=$(cpp_src:$(SRC)/%.cpp=$(INC)/%.o)

# Tests
test_src=$(wildcard ../src/test/*)
test_bin=$(test_src:$(TEST)/%.cpp=%)

#############
# UTILITIES #
#############

.PHONY: all clean cleanall run echo $(test_bin) $(main_bin)

all: $(main_bin)

clean:
	rm $(LIB)/*.o -f
	rm $(INC)/*.o -f

cleanall: clean
	rm $(LIB)/* -rf

run: $(BIN)/droplet_analysis.out
	$(BIN)/droplet_analysis.out

echo:
	echo $(main_src)
	echo $(cpp_src)

#############
# LIBRARIES #
#############

$(cpp_obj): %.o: $(SRC)/%.cpp $(INC)/%.h
	$(CXX) -o $@ -c $< $(CFLAGS)

MDBase.o: Utils.o Parameters.o Time.o
Visualization.o: Utils.o Quiver.o FieldViz.o
Readers.o: Utils.o MDBase.o Atoms.o Time.o
Atoms.o: Utils.o
Writers.o: Utils.o
Droplet.o: Utils.o CircleFit.o
Parameters.o: $(LIB)/libyaml.so Utils.o

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

$(main_bin): $(main_src) $(cpp_obj)
	$(CXX) -o $(BIN)/$@ $^ $(CFLAGS)

#########
# TESTS #
#########

test_parameters: $(TEST)/test_parameters.cpp Parameters.o
	$(CXX) -o $(BIN)/$@ $^ $(CFLAGS)
	$(BIN)/$@ $(TEST)/../run/test_config.yaml

$(TEST)/test_inputstream.out: $(ROOT)/test_inputstream.cpp $(ROOT)/Utils.cpp
	$(CXX) -o $@ $^ $(CFLAGS)

$(TEST)/test_root.out: $(ROOT)/test_root.cpp
	$(CXX) -o $@ $^ $(CFLAGS)

$(TEST)/test_parameters.out: $(LIB)/test_parameters.o $(LIB)/Parameters.o $(LIB)/Utils.o
	$(CXX) -o $@ $^ $(CFLAGS)

$(LIB)/test_parameters.o: $(ROOT)/test_parameters.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)

# $(test_bin): %: $(TEST)/%.cpp


