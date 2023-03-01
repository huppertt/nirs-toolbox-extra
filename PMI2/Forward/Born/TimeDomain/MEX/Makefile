########################################################################
# Pack everything into individual object files after doing a partial
# link.  This VASTLY simplifies the Matlab mex command, since I can
# link against the single object file and NO external libraries to get
# out a working Mex function.

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

ifndef GSL_LIBS
  GSLDIR = /usr/local/lib
  GSLLIB = gsl

  GSL_LIBS   = -L$(GSLDIR) -l$(GSLLIB)
endif

########################################################################

SHARED_ALL = mexversion.o mexStruct.o cierr.o jshack.o tdgf_zbnd.o tdgf.o
SHARED_2PT = mextd2pt.o td2pt.o $(SHARED_ALL) 
SHARED_3PT = mextd3pt.o td3pt.o $(SHARED_ALL) tdggf.o analytic.o samplevol.o

########################################################################

all: mextd2ptNB.mexglx mextd2ptZB.mexglx mextd2ptSB.mexglx\
	mextd3ptNB.mexglx mextd3ptZB.mexglx mextd3ptSB.mexglx\
	mexzbnd.mexglx

########################################################################

mextd2ptNB.mexglx: td2ptinft.o
	$(MEXLD) -shared $(MEXLDFLAGS) -o $@ $^

mextd2ptZB.mexglx: td2ptsemi.o
	$(MEXLD) -shared $(MEXLDFLAGS) -o $@ $^

mextd2ptSB.mexglx: td2ptslab.o
	$(MEXLD) -shared $(MEXLDFLAGS) -o $@ $^

mextd3ptNB.mexglx: td3ptinft.o
	$(MEXLD) -shared $(MEXLDFLAGS) -o $@ $^

mextd3ptZB.mexglx: td3ptsemi.o
	$(MEXLD) -shared $(MEXLDFLAGS) -o $@ $^

mextd3ptSB.mexglx: td3ptslab.o
	$(MEXLD) -shared $(MEXLDFLAGS) -o $@ $^

mexzbnd.mexglx: zbnd.o
	$(MEXLD) -shared $(MEXLDFLAGS) -o $@ $^

########################################################################
# Static linkage against GSL, since the shared library won't be
#  visible from the nodes.  Shared linkage above to libc and libm.

zbnd.o:  mexzbnd.o cierr.o tdgf_zbnd.o mexversion.o
	$(MEXLD) $(MEXLDFLAGS) -r -o $@ $^ $(GSL_LIBS)

td2ptinft.o: $(SHARED_2PT) inft2pt.o
	$(MEXLD) $(MEXLDFLAGS) -r -o $@ $^ $(GSL_LIBS)

td2ptsemi.o: $(SHARED_2PT) semi2pt.o
	$(MEXLD) $(MEXLDFLAGS) -r -o $@ $^ $(GSL_LIBS)

td2ptslab.o: $(SHARED_2PT) slab2pt.o
	$(MEXLD) $(MEXLDFLAGS) -r -o $@ $^ $(GSL_LIBS)

td3ptinft.o: $(SHARED_3PT) inft3pt.o
	$(MEXLD) $(MEXLDFLAGS) -r -o $@ $^ $(GSL_LIBS)

td3ptsemi.o: $(SHARED_3PT) semi3pt.o
	$(MEXLD) $(MEXLDFLAGS) -r -o $@ $^ $(GSL_LIBS)

td3ptslab.o: $(SHARED_3PT) slab3pt.o
	$(MEXLD) $(MEXLDFLAGS) -r -o $@ $^ $(GSL_LIBS)

mexversion.o: $(MATLAB_ROOT)/extern/src/mexversion.c
	$(MEXCC) $(MEXCFLAGS) -c $(MATLAB_ROOT)/extern/src/mexversion.c

# Explicitly declare the non-trivial dependancies

cierr.o:	td2pt.h
inft2pt.o:	td2pt.h
inft3pt.o:	td2pt.h td3pt.h
mextd2pt.o:	td2pt.h
mextd3pt.o:	td2pt.h td3pt.h
mexzbnd.o:	td2pt.h
semi2pt.o:	td2pt.h
semi3pt.o:	td2pt.h td3pt.h
slab2pt.o:	td2pt.h
slab3pt.o:	td2pt.h td3pt.h
td2pt.o:	td2pt.h
td3pt.o:	td2pt.h td3pt.h
tdgf.o:		td2pt.h
tdgf_zbnd.o:	td2pt.h
tdggf.o:	td2pt.h

.c.o:
	$(MEXCC) $(MEXCFLAGS) -c $*.c

###########################################################################

clean:
	-/bin/rm -f core *.o *.mexglx

###########################################################################

