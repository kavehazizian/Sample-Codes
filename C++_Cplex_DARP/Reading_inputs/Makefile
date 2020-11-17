SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /media/kaveh/DATA/Installed_cplex_linux/cplex
CONCERTDIR    = /media/kaveh/DATA/Installed_cplex_linux/concert
MYDIR         = /media/kaveh/DATA/C++_Cplex_DARP/Reading_inputs

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0
CC  = gcc -O0
JAVAC = javac 

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
COPT  = -m64 -fPIC -fno-strict-aliasing
JOPT  = -classpath $(CPLEXDIR)/lib/cplex.jar -O

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread
CLNFLAGS  = -lcplex -lm -lpthread
JAVA      = java  -d64 -Djava.library.path=$(CPLEXDIR)/bin/x86-64_linux -classpath $(CPLEXJARDIR):


all:
	
	make all_cpp
	

execute: all
	
	make execute_cpp
	
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = $(MYDIR)



CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
JCFLAGS = $(JOPT)


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------




CPP_EX = main \
         


all_cpp: $(CPP_EX)





execute_cpp: $(CPP_EX)
	 
	 ./main	 
	 
# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf  $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------
#
# The examples
#




main: main.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o main main.o $(CCLNFLAGS)
main.o: $(MYDIR)/main.cpp
	$(CCC) -c $(CCFLAGS) $(MYDIR)/main.cpp -o main.o




# Local Variables:
# mode: makefile
# End:
