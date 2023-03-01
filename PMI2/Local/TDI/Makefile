########################################################################
# Pack everything into individual object files after doing a partial
# link.  This VASTLY simplifies the Matlab mex command, since I can
# link against the single object file and NO external libraries to get
# out a working Mex function.

# Matlab version to link against.  Needs to match executable
TMW_ROOT = /usr/local/matlab

########################################################################

CC   = gcc
COPT = -O2 
CWRN = -Wall 
CINC = -I$(TMW_ROOT)/extern/include
CDCL =

LD   = ld
CLIB = 

# MEX options, as best as I can tell from the matlab scripts

MLIBDIR = $(TMW_ROOT)/extern/lib/glnx86
MBINDIR = $(TMW_ROOT)/bin/glnx86

MLIB = --version-script=$(MLIBDIR)/mexFunction.map\
	--rpath-link=$(MBINDIR) --rpath-link=$(MLIBDIR)\
	-L$(MLIBDIR) -L$(MBINDIR) -lmx -lmex -lm

########################################################################

CFLAGS  = $(CWRN) $(COPT) $(CDCL) $(CINC)

########################################################################

all: defiber.mexglx

defiber.mexglx: defiber.o
	$(LD) -shared $(LOPT) -o $@ $^ $(MLIB)

defiber.o: defiber.c

clean:
	-/bin/rm -f core *.o *.mexglx

###########################################################################

