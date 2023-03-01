########################################################################
#
# Makefile for PMI MEX files.
#
########################################################################

########################################################################
# Matlab paths, change to match local system

# Root matlab directory
MATLAB_ROOT = /usr/local/matlab

########################################################################
# Gnu Scientific Library path and library name (without the leading "lib")

GSLDIR = /usr/local/lib
GSLLIB = gsl

########################################################################
# Compiler flags

CC   = gcc 
COPT = -O2 -funroll-all-loops
CWRN = -Wall 
CINC = -I$(MATLAB_ROOT)/extern/include
CDCL =

########################################################################
# Linker flags

LD   = ld

########################################################################
# END OF USER CONFIGURATION

# MEX options, as extracted from the stock matlab scripts

MLIBDIR = $(MATLAB_ROOT)/extern/lib/glnx86
MBINDIR = $(MATLAB_ROOT)/bin/glnx86

MLIB =  --version-script=$(MLIBDIR)/mexFunction.map\
	--rpath-link=$(MBINDIR) --rpath-link=$(MLIBDIR)\
	-L$(MLIBDIR) -L$(MBINDIR) -lmx -lmex -lm

# Repackage above flags into compler flags I can pass to other Makefiles

MEXCC      = $(CC)
MEXLD      = $(LD)

MEXCFLAGS  = $(CFLAGS) $(CWRN) $(COPT) $(CDCL) $(CINC)
MEXLDFLAGS = $(LDFLAGS)

MEXLIBS    = $(MLIB) $(LIBS)
GSL_LIBS   = -L$(GSLDIR) -l$(GSLLIB)

########################################################################

# Make configuration variables availabe to sub-makefiles
export MATLAB_ROOT MEXCC MEXLD MEXCFLAGS MEXLDFLAGS MEXLIBS GSL_LIBS

# Where to find the files I need to build
SUBDIRS =\
	Forward/Born/MEX \
	Forward/Born/FrequencyDomain/MEX \
	Forward/Born/TimeDomain/MEX

all:
	for mdir in $(SUBDIRS) ; do make -C $$mdir   all; done

clean:
	for mdir in $(SUBDIRS) ; do make -C $$mdir clean; done
