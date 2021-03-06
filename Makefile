## Simple make file, It does NOT check for interdependencies of modules.
## It probably is best to do make clean; make if you change a module.

## Make files are extremely picky about syntax. 
## Indentations have to be done with tabulators!
## Absolutely nothing can follow a continuation backslash.

## List main program module last
## 
SOURCES = precise.f90 \
	io.f90 \
	vtkxmlmod.f90 \
	panelgeometry.f90 \
	constants.f90 \
	influence.f90 \
	slae.f90 \
	flowsolve.f90

## Define name of main program
PROGRAM = flowsolve

# Compiler
FF = gfortran

# Delete program
# Linux
RM = rm -f
# DOSe
#RM = del

## Compiler options
# for the Intel Fortran 90 compiler
# CFLAGS = -c -fast -heap-arrays

# for the g95 compiler
CFLAGS = -c -O3

#Linker Options
# for the Intel Fortran 90 compiler
# LDFLAGS = -fast -heap-arrays

# for the g95 compiler
# LDFLAGS = -02 # -c 03
LDFLAGS = -O3 -march=native -ffast-math -funroll-loops # -g -02


## Probably no changes necessary below this line
OBJECTS = $(SOURCES:.f90=.o)

all: $(SOURCES) $(PROGRAM)


$(PROGRAM): $(OBJECTS)
	$(FF) $(LDFLAGS) $(OBJECTS) -o $@

$(OBJECTS): %.o: %.f90
	$(FF) $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJECTS) *.mod 

realclean:
	$(RM) $(OBJECTS) *.mod $(PROGRAM)