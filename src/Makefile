# list of executable files to produce
PAR_SPLIT = ParSplit
FIX = Fix

# .o files necessary to build the executables
OBJ_PAR_SPLIT = Utils.o Floating.o StateTree.o Tree.o Fraction.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o DrawTreeGenericCut.o Piecewise.o ContinuousParsimony.o ParsimoniousSplit.o 
OBJ_FIX = Utils.o StateTree.o Tree.o FixData.o 
########### MODIFY ONLY TO CHANGE OPTIONS ############

# compiler and its options
CC = gcc
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo

# linker and its options
LD = $(CC)

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(PAR_SPLIT) $(FIX)


# build the executable
$(PAR_SPLIT): $(OBJ_PAR_SPLIT)
	$(LD) $^ -o $@ -lm  -lcairo -O3
	
# build the executable
$(FIX): $(OBJ_FIX)
	$(LD) $^ -o $@ -lm

############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release

-include Makefile.d
