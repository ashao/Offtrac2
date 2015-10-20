##############################
# Things you need to set
##############################

NETCDF = /usr/local/netcdf4
# Other options:
#NETCDF = /ltraid3/ashao/waddle/local/
#NETCDF = /cm/shared/apps/netcdf/gcc/64/4.3.0/
#NETCDF = /usr/local

##############################
# Things you might need to set
##############################

OPENMP = 0

DEBUG = 1
# 0 is no debugging at all
# 1 is -g
# 2 is more runtime debugging

CC = gcc

#ARCHFLAGS = -march="athlon64"
#ARCHFLAGS = -march="i686"

#CFLAGS += -lpthread
#CFLAGS += -pipe
#CFLAGS += -static
#CFLAGS += -funroll-loops

#OUTNAME = offtrac.cfcs
SRCDIR = src

##############################
# Things you shouldn't change
##############################

ifeq ($(CC),gcc)
# Detect if we're dealing with gcc 4.6 or later
CC_IS_GCC46PLUS := $(shell echo `gcc -dumpversion | cut -f1-2 -d.` \>= 4.6 | bc )
# Detect if we're compiling to a 64-bit target
CC_IS_GCCX64 := $(shell gcc -v 2>&1 | egrep -c '^Target:.*x86_64' )
endif

ifeq ($(strip $(OUTNAME)),)
OUTNAME = offtrac
endif

ifeq ($(strip $(ARCHFLAGS)),)
ARCHFLAGS = -march=native
endif

ifeq ($(CC),icc)
# icc flags
CFLAGS += -Wall -O2 -ip -ipo -inline-level=2 -xHOST
else ifeq ($(CC),gcc)
# gcc flags
CFLAGS += -Wall -Wno-unknown-pragmas -Werror
ifeq ($(CC_IS_GCC46PLUS),1)
# gcc >= 4.6 flags
CFLAGS += -Ofast -flto
else
# gcc < 4.6 flags
CFLAGS += -O3 -ffast-math
endif
ifneq ($(CC_IS_GCCX64),1)
# 32-bit gcc flags
CFLAGS += -mfpmath=sse
endif
endif

ifneq ($(DEBUG),0)
CFLAGS += -g
ifneq ($(DEBUG),1)
ifeq ($(CC),icc)
CFLAGS += -fp-trap-all=all,noinexact
else ifeq ($(CC),gcc)
CFLAGS += -ftrapv
endif
endif
endif

CFLAGS += -lm

ifneq ($(OPENMP),0)
CFLAGS += -fopenmp
endif

LDFLAGS = -I$(NETCDF)/include/ -L$(NETCDF)/lib -lnetcdf -lrt

# Make sure to add the Gibbs Seawater routines
GSW_DIR= $(SRCDIR)/gsw_src/
GSW_LIB= $(GSW_DIR)/libgswteos-10.so
GSW_INC= -I$(GSW_DIR)

OFFSRC = $(wildcard $(SRCDIR)/*.c)

$(OUTNAME): $(OFFSRC) $(SRCDIR)/init.h Makefile 
	echo compiling $(OUTNAME)
	$(CC)  $(OFFSRC) -o $(OUTNAME) $(CFLAGS) $(LDFLAGS) $(GSW_LIB)		

propre:
	rm -f *.o *.u

cleaner:
	rm -f offtrac *.o *.u 

clean:
	rm -f offtrac $(SRCDIR)/*.o *.u off

