########################################################################
# Pack everything into individual object files after doing a partial
# link.  This VASTLY simplifies the Matlab mex command, since I can
# link against the single object file and NO external libraries to get
# out a working Mex function.

# Matlab version to link against.  Needs to match executable.  Not
# used if called by top-level makefile.  Change to match local system.

# If not called by top-level makefile, supply usable defaults

ifndef MATLAB_ROOT
  MATLAB_ROOT= /usr/local/matlab
endif

ifndef MEXCC
  MEXCC      = $(CC)
endif

ifndef MEXCFLAGS
  MEXCFLAGS  = $(CFLAGS) -I$(MATLAB_ROOT)/extern/include
endif

ifndef MEXLD
  MEXLD      = $(LD)
endif

ifndef MEXLDFLAGS
  MEXLDFLAGS = $(LDFLAGS)
endif

ifndef MEXLIBS
  MLIBDIR = $(MATLAB_ROOT)/extern/lib/glnx86
  MBINDIR = $(MATLAB_ROOT)/bin/glnx86

  MLIB = --version-script=$(MLIBDIR)/mexFunction.map\
	--rpath-link=$(MBINDIR) --rpath-link=$(MLIBDIR)\
	-L$(MLIBDIR) -L$(MBINDIR) -lmx -lmex -lm

  MEXLIBS    = $(MLIB)
endif

########################################################################

all: FD_GF.mexglx mexversion.o

FD_GF.mexglx: FD_GF.o mexversion.o
	$(MEXLD) -shared  $(MEXLDFLAGS) -o $@ $^ $(MEXLIBS)

FD_GF.o: FD_GF.c
	$(MEXCC) $(MEXCFLAGS) -c FD_GF.c

mexversion.o: $(MATLAB_ROOT)/extern/src/mexversion.c
	$(MEXCC) $(MEXCFLAGS) -c $(MATLAB_ROOT)/extern/src/mexversion.c

.c.o:
	$(MEXCC) $(MEXCFLAGS) -c $*.c

###########################################################################

clean:
	-/bin/rm -f core *.o *.mexglx

###########################################################################

