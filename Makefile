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

#CFLAGS += -pipe
#CFLAGS += -static
#CFLAGS += -funroll-loops
#LDFLAGS += -lpthread

#OUTNAME = offtrac.cfcs
SRCDIR = src
OBJDIR = obj
DEPDIR = .d

##############################
# Things you probably shouldn't change
##############################

# Make compilation use 1 compiler per physical CPU core
# This is a slightly conservative choice, and could be changed
# up or down depending on user tolerance
CORECOUNT := $(shell lscpu -p | egrep -v ^\# | cut -d, -f2 | sort -u | wc -l )
MAKEFLAGS = -j$(CORECOUNT)

# Make the output directory and dependency directory
$(shell mkdir -p $(OBJDIR) $(DEPDIR))

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

ifneq ($(OPENMP),0)
CFLAGS += -fopenmp
endif

CFLAGS += -I$(NETCDF)/include/
LDFLAGS = -L$(NETCDF)/lib -lnetcdf

LDFLAGS += -lrt -lm

DEPFLAGS = -MMD -MP -MT $@ -MF $(DEPDIR)/$(*F).Td

# Make sure to add the Gibbs Seawater routines
GSW_DIR= $(SRCDIR)/gsw_src/
GSW_LIB= $(GSW_DIR)/libgswteos-10.so
GSW_INC= -I$(GSW_DIR)

# Compile all files in $(SRCDIR)
SRCS := $(wildcard $(SRCDIR)/*.c)
# Put the objects into $(OBJDIR) so as not to clutter up our sources
OBJS := $(subst $(SRCDIR),$(OBJDIR),$(SRCS:.c=.o))

# Make the executable by linking together all the object files.
$(OUTNAME): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(GSW_LIB)

# This is the make rule for compiling a source file into an object
# file. Each object file depends upon its corresponding source file,
# as well as the Makefile. It also creates the dependency-tracking
# info and places it in $(DEPDIR), as a separate step so compilation
# failures won't change it. The header file dependencies for each
# source file are thus tracked by the -include directive at the bottom
# of this Makefile.
$(OBJS) : $(OBJDIR)/%.o : $(SRCDIR)/%.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@ $(DEPFLAGS)
	@mv -f $(DEPDIR)/$(*F).Td $(DEPDIR)/$(*F).d

propre:
	rm -f $(OBJDIR)/*.o $(DEPDIR)/* *.u

cleaner:
	rm -f $(OUTNAME) $(OBJDIR)/*.o $(DEPDIR)/* *.u 

clean:
	rm -f $(OUTNAME) $(OBJDIR)/*.o $(DEPDIR)/* *.u

# This tells make that these rules don't actually make anything
.PHONY: clean cleaner propre

# This ensures make won't fail if the dependency file doesn't exist
$(DEPDIR)/%.d: ;

# This gives make all the dependency information, but won't cause make
# to fail if they don't exist
-include $(subst $(SRCDIR),$(DEPDIR),$(SRCS:.c=.d))
