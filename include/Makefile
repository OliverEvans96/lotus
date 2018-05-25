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
TEST=$(BASE)/test
TEST_SRC=$(TEST)/src
TEST_INC=$(TEST)/include

# C++ Compiler
CXX=g++
# Root flags
ROOTFLAGS=`root-config --glibs --cflags` -lMathMore -lMinuit
# Specify rpath for custom shared library path
CFLAGS=$(ROOTFLAGS) -I$(INC) -L$(LIB) -lyaml -Wl,-rpath=$(LIB)
TESTFLAGS=-I$(TEST_INC)

# Test data file
testdata_zipped=$(wildcard $(DATA)/*.tar.gz)
testdata_unzipped=$(testdata_zipped:%.tar.gz=%)

# Main executable
main_src=$(SRC)/droplet_analysis.cpp
main_bin=$(main_src:$(SRC)/%.cpp=%)

# All other C++ source, header, and object files
cpp_src=$(filter-out $(main_src), $(wildcard $(SRC)/*.cpp))
cpp_hdr=$(cpp_src:$(SRC)/%.cpp=$(INC)/%.h)
cpp_obj=$(cpp_src:$(SRC)/%.cpp=$(INC)/%.o)

# Main test
main_test_src=$(TEST_SRC)/test_main.cpp
main_test_obj=$(main_test_src:$(TEST_SRC)/%.cpp=$(TEST_INC)/%.o)
main_test_bin=$(BIN)/test

# Tests
test_src=$(filter-out $(main_test_src), $(wildcard $(TEST_SRC)/*.cpp))

# Eliminate this
test_bin=$(test_src:$(TEST_SRC)/%.cpp=%)

$(main_test_obj): $(main_test_src)
	$(CXX) -o $@ -c $^ $(TESTFLAGS)

$(main_test_bin): $(main_test_obj) $(test_src) $(cpp_obj)
	$(CXX) -o $@ $^ $(CFLAGS) $(TESTFLAGS)

#############
# UTILITIES #
#############

.PHONY: all clean cleanall run echo test $(test_bin) $(main_bin)

all: $(main_bin)

clean:
	rm $(LIB)/*.o -f
	rm $(INC)/*.o -f
	rm $(main_test_obj) -f

cleandata:
	rm $(datafile_unzipped) -f

cleanall: clean cleandata
	rm $(LIB)/* -rf

run: $(BIN)/droplet_analysis.out
	$(BIN)/droplet_analysis.out

echo:
	echo $(test_src)

test: $(main_test_bin)
	$<

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

########
# DATA #
########

$(testdata_unzipped): %: %.tar.gz
	tar -xvzf $@

########################
# TESTS (Remove these) #
########################

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


