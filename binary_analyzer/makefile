
FILES = binary_analyzer.cpp

EXES = $(addsuffix .exe, $(basename $(FILES)))

OPT = -std=c++11


all: $(EXES)


%.exe: %.cpp
	$(CXX) $(OPT) -o $@ $<


clean:
	rm -f $(EXES) *~ .*~ 

cleanall: clean


