
FILES = analizza_binario.cpp

EXES = $(addsuffix .exe, $(basename $(FILES)))

OPT = -std=c++11


all: $(EXES)


%.exe: %.cpp
	$(CXX) $(OPT) -c -o $@ $<


clean:
	rm -f $(EXES) *~ .*~ 

cleanall: clean


