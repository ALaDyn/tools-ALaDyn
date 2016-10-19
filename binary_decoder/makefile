EXE      = binary_decoder.exe
OPT      = -std=c++11

OBJECTS  = binning.o filter.o binary_decoder.o grid_file_decoder.o phasespace_file_decoder.o parameters.o swap_tools.o

all: $(EXE)

perf: 	OPT += -O3
perf: 	$(EXE)

debug:	OPT += -Og -g 
debug:	$(EXE)

CINECA:	OPT = -std=c++0x -O3 -DCINECA
CINECA:	$(EXE)


$(EXE): $(OBJECTS)
	$(CXX) $(OPT) -o $(EXE) $(OBJECTS)


parameters.o: parameters.cpp parameters.h binary_decoder.h
	$(CXX) $(OPT) -I. -c parameters.cpp

binning.o: binning.cpp binning.h parameters.o binary_decoder.h
	$(CXX) $(OPT) -I. -c binning.cpp

filter.o: filter.cpp filter.h parameters.o binary_decoder.h
	$(CXX) $(OPT) -I. -c filter.cpp

swap_tools.o: swap_tools.cpp swap_tools.h binary_decoder.h
	$(CXX) $(OPT) -I. -c swap_tools.cpp

grid_file_decoder.o: grid_file_decoder.cpp grid_file_decoder.h binary_decoder.h parameters.h swap_tools.o parameters.o
	$(CXX) $(OPT) -I. -c grid_file_decoder.cpp

phasespace_file_decoder.o: phasespace_file_decoder.cpp phasespace_file_decoder.h binary_decoder.h parameters.h swap_tools.o parameters.o binning.o filter.o
	$(CXX) $(OPT) -I. -c phasespace_file_decoder.cpp

binary_decoder.o: binary_decoder.cpp binary_decoder.h parameters.h phasespace_file_decoder.o grid_file_decoder.o
	$(CXX) $(OPT) -I. -c binary_decoder.cpp

clean:
	rm -f $(OBJECTS)

cleanall: clean
	rm -f $(EXE) *~ .*~ 


