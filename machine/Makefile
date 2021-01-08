## < ====================      DEVELOPER'S NOTE      ==================== > ##
##                                                                          ##
##  Makefile                                                                ##
##                                                                          ##
##  When you compile, "HDF5" and "MPICH2" (or "MPICH3") must be installed.  ##
##                                                                          ##
##  SYSTYPE : Select the computer to compile.                               ##
##  CFILE   : Select the code which you want to compile.                    ##
##                                                                          ##
##  Skip lines starting with "#" or line breaks.                            ##
##                                                                          ##
##  > make clean                                                            ##
##  > make                                                                  ##
##                                                                          ##
## > ==================================================================== < ##


# SYSTYPE="macbook"
SYSTYPE="macbook2"

# CFILE="FX_demo"
CFILE="FFT_FX"
# CFILE="test"


ifeq ($(SYSTYPE),"xc")
CC       = cc
DFLAGS   = -DH5_USE_16_API #-DGRACKLE=3
CFLAGS   = -O3 -Wall $(DFLAGS)
MPICHINC =
MPIGHLIB =
LIBS     = -lm -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"pluto")
CC 	     = mpicc
# DFLAGS   = -DGRACKLE=3
DFLAGS   = -DH5_USE_16_API #-DGRACKLE=3
CFLAGS 	 = -O3 -Wall $(DFLAGS)
INCDIRS  = -I/home/ishimizu/local/include
LIBDIRS  = -L/home/ishimizu/local/lib
MPICHINC =
MPIGHLIB =
LIBS 	   = -lm -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"macbook")
CC 	     = mpicc
DFLAGS   = -DH5_USE_16_API #-DGRACKLE=3
CFLAGS 	 = -O3 -Wall $(DFLAGS)
INCDIRS  = -I/opt/local/include
LIBDIRS  = -L/opt/local/lib
MPICHINC = -I/Users/yuza/local/mpich2/include
MPIGHLIB = -L/Users/yuza/local/mpich2/lib
LIBS     = -lm -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"macbook2")
CC 	     = mpicc
# CC 	     = gcc
# DFLAGS   = -DGRACKLE=3
DFLAGS   = -DH5_USE_16_API #-DGRACKLE=3
CFLAGS 	 = -O3 -Wall $(DFLAGS)
INCDIRS  = -I/opt/local/include
LIBDIRS  = -L/opt/local/lib
MPICHINC = -I/Users/nakatayuuki/local/mpich2/include
MPIGHLIB = -L/Users/nakatayuuki/local/mpich2/lib
# MPICHINC = -I/Users/nakatayuuki/local/mpich3/include
# MPIGHLIB = -L/Users/nakatayuuki/local/mpich3/lib
LIBS 	= -lm -lhdf5 -lgsl -lgslcblas -lfftw3 -lz
endif



ifeq ($(CFILE),"FX_demo")
INCL = common.h fx_demo.h
MAIN = common.o fx_demo6.o main_fx.o
ALL  = main_fx

main: $(MAIN)
	$(CC) $(CFLAGS) -o main_fx $(LIBDIRS) $(LIBS) $(MAIN)
endif

ifeq ($(CFILE),"FFT_FX")
INCL = common.h fft_fx.h
MAIN = common.o fft_fx16.o main_fft.o
ALL  = main_fft

main: $(MAIN)
	$(CC) $(CFLAGS) -o main_fft $(LIBDIRS) $(LIBS) $(MAIN)
endif

ifeq ($(CFILE),"test")
INCL = #common.h fft_fx.h
MAIN = test_fftw3.o
ALL  = test_fftw3

main: $(MAIN)
	$(CC) $(CFLAGS) -o test_fftw3 $(LIBDIRS) $(LIBS) $(MAIN)
endif


## Compile
.c.o: $*.c
	$(CC) -c $*.c $(CFLAGS) $(INCDIRS) $(MPICHINC)


## Clean

clean:
	rm -f $(ALL) *.o *~
