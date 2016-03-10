
FILES = InGenUO.cpp

EXES = $(addsuffix .exe, $(basename $(FILES)))

OPT = -std=c++11 -I.


all: $(EXES)


%.exe: %.cpp
	$(CXX) $(OPT) -o $@ $<


clean:
	rm -f $(EXES) *~ .*~ 

cleanall: clean


