COMPILER = g++
EXE      = leggi_bin.exe
OPT      = -std=c++11

OBJECTS  = binnaggio.o filtro.o leggi_binario_ALaDyn_fortran.o leggi_griglia.o leggi_particelle.o parametri.o scrittura.o swap_tools.o

all: $(EXE)

perf: 	OPT += -O3
perf: 	$(EXE)

debug:	OPT += -Og -g 
debug:	$(EXE)

CINECA:	OPT = -std=c++0x -O3 -DCINECA
CINECA:	$(EXE)


$(EXE): $(OBJECTS)
	$(COMPILER) $(OPT) -o $(EXE) $(OBJECTS)


parametri.o: parametri.cpp leggi_binario_ALaDyn_fortran.h
	$(COMPILER) $(OPT) -c parametri.cpp

binnaggio.o: binnaggio.cpp binnaggio.h parametri.o leggi_binario_ALaDyn_fortran.h
	$(COMPILER) $(OPT) -c binnaggio.cpp

filtro.o: filtro.cpp filtro.h parametri.o leggi_binario_ALaDyn_fortran.h
	$(COMPILER) $(OPT) -c filtro.cpp

scrittura.o: scrittura.cpp scrittura.h parametri.o leggi_binario_ALaDyn_fortran.h
	$(COMPILER) $(OPT) -c scrittura.cpp

swap_tools.o: swap_tools.cpp swap_tools.h leggi_binario_ALaDyn_fortran.h
	$(COMPILER) $(OPT) -c swap_tools.cpp

leggi_griglia.o: leggi_griglia.cpp leggi_griglia.h leggi_binario_ALaDyn_fortran.h swap_tools.o parametri.o
	$(COMPILER) $(OPT) -c leggi_griglia.cpp

leggi_particelle.o: leggi_particelle.cpp leggi_particelle.h leggi_binario_ALaDyn_fortran.h swap_tools.o parametri.o binnaggio.o filtro.o scrittura.o
	$(COMPILER) $(OPT) -c leggi_particelle.cpp

leggi_binario_ALaDyn_fortran.o: leggi_binario_ALaDyn_fortran.cpp leggi_binario_ALaDyn_fortran.h leggi_particelle.o leggi_griglia.o parametri.o
	$(COMPILER) $(OPT) -c leggi_binario_ALaDyn_fortran.cpp

clean:
	rm -f $(OBJECTS)

cleanall:
	rm -f $(OBJECTS) $(EXE) *~ .*~ 


