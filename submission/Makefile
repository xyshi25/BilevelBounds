SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /usr/local/cplex/cplex_studio128/cplex/
CONCERTDIR    = /usr/local/cplex/cplex_studio128/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0 -std=c++14 


# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD


# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl



all:
	make all_cpp
	
	
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 



#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
CPP_EX = KnapInterdiction BVertexCover BMST



all_cpp: $(CPP_EX)
	/bin/rm -rf *.o
	/bin/mkdir -p bin
	/bin/mv $(CPP_EX) bin

	 


clean :
	/bin/rm -rf bin
	/bin/rm -rf Log
	/bin/rm -rf Dataset
	/bin/rm -rf *.log *.txt



# ------------------------------------------------------------
#
# The examples
#
SRCDIR = src
KnapInterdiction: $(SRCDIR)/KnapInterdiction.cpp 
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/KnapInterdiction.cpp -o KnapInterdiction.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o KnapInterdiction KnapInterdiction.o  $(CCLNFLAGS)

BVertexCover: BVertexCover.o 
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o BVertexCover BVertexCover.o  $(CCLNFLAGS)
BVertexCover.o: $(SRCDIR)/BVertexCover.cpp 
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/BVertexCover.cpp -o BVertexCover.o

BMST: BMST.o 
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o BMST BMST.o  $(CCLNFLAGS)
BMST.o: $(SRCDIR)/BMST.cpp 
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/BMST.cpp -o BMST.o


# ------------------------------------------------------------





# Local Variables:
# mode: makefile
# End:
