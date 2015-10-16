#OUTNAME = offtrac.cfcs
#CC=gcc

ifeq ($(strip $(OUTNAME)),)
OUTNAME=offtrac
endif

#CC=gcc
#CFLAGS = -O3 -g -fopenmp -lm -lpthread --fast-math -march="athlon64" -pipe -static

# waddle with gcc
CC = gcc
CFLAGS =  -fopenmp  -g -O3 -march="athlon64" -lm  
LDFLAGS = -I/ltraid3/ashao/waddle/local/include/ -L/ltraid3/ashao/waddle/local/lib -lnetcdf

# gaggle with gcc
#CC = gcc
#CFLAGS = -mfpmath=sse -flto -march=native -funroll-loops -fopenmp  -g -pipe -Ofast 
#LDFLAGS = -I/cm/shared/apps/netcdf/gcc/64/4.3.0/include -L/cm/shared/apps/netcdf/gcc/64/4.3.0/lib

#CFLAGS = -g 
#debug for gcc
#CFLAGS = -fbounds-check -g

#CFLAGS = -ffast-math -O2 -march="athlon64" -pipe
#CFLAGS = -ffast-math -O2 -pipe

#ocean
#CFLAGS = -O1 -pipe  -pg -g
#CFLAGS = -ffast-math -O3 -march="athlon64" -pipe -g

#yucatan with icc:
#CFLAGS = -O1 -march="i686"  -pipe
#yucatan gcc debugging
#CFLAGS = -O1 -march="i686"  -pipe  -pg -Wuninitialized
#CFLAGS = -g -march="i686"  -pipe -fbounds-check -Wall
#yucatan icc debugging with idb
#CFLAGS = -march="i686"  -pipe -g -lm

#CDFFLAGS = -lm -I/usr/local/include -L/usr/local/lib -E -source_listing
#CDFFLAGS = -g -lm -lnetcdf -I/usr/local/include -L/usr/local/lib 

#CDFFLAGS = -lm -lnetcdf -I/usr/local/include -L/usr/local/lib 
#CDFFLAGS = -lm -lnetcdf -I/usr/include/netcdf-3 -L/usr/lib64 

#LDFLAGS = -lm -ftrap -lnetcdf -I/usr/local/include -L/usr/local/lib
#LDFLAGS = -lm -ftrap -lnetcdf -I/usr/include -L/usr/lib
#LDFLAGS = -lm -lnetcdf -I/usr/local/include -L/usr/local/lib
SRCDIR = src

#pendragon with icc
#CC=icc
#CFLAGS= -ip -ipo -inline-level=2 -xhost -O3 -g -openmp -lpthread 
#CFLAGS= -ip -ipo -inline-level=2 -xhost -O3 -g -lpthread 

CDFFLAGS= 
# Make sure to add the Gibbs Seawater routines
GSW_DIR= $(SRCDIR)/gsw_src/
GSW_LIB= $(GSW_DIR)/libgswteos-10.so
GSW_INC= -I$(GSW_DIR)

OFFSRC = $(SRCDIR)/offtrac.c $(SRCDIR)/read.c \
	$(SRCDIR)/initialize.c $(SRCDIR)/iocdf.c $(SRCDIR)/par_IO.c \
	$(SRCDIR)/tracadv.openmp.c $(SRCDIR)/step.c \
        $(SRCDIR)/alloc_trac.c \
        $(SRCDIR)/alloc_util.c $(SRCDIR)/alloc_arrays.c \
	$(SRCDIR)/masks.c $(SRCDIR)/set_metrics.c\
	$(SRCDIR)/util.c \
        $(SRCDIR)/timekeeper.c \
        $(SRCDIR)/output_variables.c \
        $(SRCDIR)/tracer_utilities.c \
        $(SRCDIR)/gas_exchange.c \
	$(SRCDIR)/ideal_age.c \
	$(SRCDIR)/cfc11.c $(SRCDIR)/cfc12.c $(SRCDIR)/sf6.c \
	$(SRCDIR)/ttd_bp.c \
	$(SRCDIR)/n2_module.c $(SRCDIR)/ar_module.c \
	$(SRCDIR)/oxygen.c $(SRCDIR)/phosphate.c $(SRCDIR)/biotic.c \
	$(SRCDIR)/conservation_check.c

offtrac: $(OFFSRC) $(SRCDIR)/init.h Makefile 
	echo compiling $(OUTNAME)
	$(CC)  $(OFFSRC) -o $(OUTNAME) $(CDFFLAGS) $(CFLAGS) $(LDFLAGS) $(GSW_LIB)		

propre:
	rm -f *.o *.u

cleaner:
	rm -f offtrac *.o *.u 

clean:
	rm -f offtrac $(SRCDIR)/*.o *.u off

