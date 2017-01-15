
FILES = InGenUO.cpp

EXES = $(addsuffix .exe, $(basename $(FILES)))

OPT = -std=c++11 -I../ext/jsoncons/src


all: $(EXES)


%.exe: %.cpp
	$(CXX) $(OPT) -o $@ $<


clean:
	rm -f $(EXES) *~ .*~ 

cleanall: clean


