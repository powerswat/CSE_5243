CC = g++
OPTFLAGS = -O3 -DNDEBUG
COPTIONS = -Wno-write-strings
INCLUDES = -I.
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
LIBS = -lsim -lm -lgsl -lgslcblas 

LDOPTIONS = -L.
LD = $(CC) $(LDOPTIONS) 

AR = ar rv
RANLIB = ranlib

SIMLIB = libsim.a

LIBOBJS = bayeslsh.o allpairs.o lsh.o io.o utils.o \
          allpairsRealCosine.o allpairsBinary.o myqsort.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c

default: allpairs bayeslsh lsh

$(SIMLIB): $(LIBOBJS) *.h
	$(AR) $@ $(LIBOBJS)
	$(RANLIB) $@

allpairs: allpairs_main.o $(SIMLIB)
	$(LD) -o $@ allpairs_main.o $(LIBS)
	chmod 755 $@

bayeslsh: bayeslsh_main.o $(SIMLIB)
	$(LD) -o $@ bayeslsh_main.o $(LIBS)
	chmod 755 $@

lsh: lsh_main.o $(SIMLIB)
	$(LD) -o $@ lsh_main.o $(LIBS)
	chmod 755 $@

sparsify: *.c *.h
	$(CC) $(CFLAGS) -c *.c
	$(LD) $(LIBS) *.o -o $@
	chmod 755 $@

clean:
	rm -f *.o libsim.a allpairs bayeslsh lsh

