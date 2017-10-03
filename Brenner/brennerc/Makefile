# My /usr/bin/gcc is version 2.91.66, which has shown buggy behavior
# when compiling past versions of this code.  Tim Freeman  5 Sep 2000.
CC=/usr/local/bin/gcc
# .DELETE_ON_ERROR means that make will delete files it was trying to
# generate, if the command that was generating them returns an error
# status.  This is always a good thing but it is incompatible with
# some historical behavior that apparently somebody wanted.  See the
# gnu-make info file and search for ".DELETE_ON_ERROR".
.PHONY: .DELETE_ON_ERROR
.DELETE_ON_ERROR:
# $(BRENNERMD) is the Fortran version, so we can compare outputs when testing.
BRENNERMD:=/home/pcm/nanocad/brenner/Subroutines/General/brennermd
ifeq ($(wildcard $(BRENNERMD)),)
   BRENNERMD:=brennermd
endif
RASAPIDIR = ../rasmolapi
INCLUDES = -I$(RASAPIDIR)
# -funroll-loops makes exactly no difference, so leave it out.
# 5.24 5.22 5.22 5.22 elapsed without -funroll-loops, 
# 5.63 5.22 5.22 5.22 elapsed with -funroll-loops.  
# 
# -mcpu=pentiumpro doesn't help either.
# 
# On my battery-powered portable, there are large variations in the
# time things take, so don't use it for timing runs.  The same run can
# take 7 or 8 or 9 or 10 seconds.  In contrast, times from a 600 Mhz
# Gateway Pentium plugged into the wall are very reproducible after
# the first run.  The portable still has large variations when it's
# plugged into the wall. 
#
# Tim Freeman 24 Jul 2000
#
C_FLAGS = -Wall -Wcast-qual -Wno-unused -W -ffast-math -fno-exceptions -Werror -Winline
ifeq ($(DEBUG),)
   C_FLAGS += -fomit-frame-pointer -DNDEBUG -O3 
else
   C_FLAGS += -g
# Comment in next line when I'm using the debugger to get information
# about performance.
#   C_FLAGS += -DNDEBUG -O3
endif
ifneq ($(INFINITE_CUBE),)
   C_FLAGS += -DINFINITE_CUBE
endif
ifneq ($(N2_ASYM),)
   C_FLAGS += -DN2_ASYM
endif

# With comparable code, if SINGLE was defined a timing run took 4.38
# seconds and otherwise it took 4.88 seconds, about a 10% speedup.

# The convention is that a "Double" in the source is a double
# precision floating point number that has to be double precision so
# things are accurate enough.  In constrast, a "double" in the source
# has to be a double to make the compiler happy.  So the code should
# compile, run, and give wobbly answers if we have -DDouble=float.
# Once I'm sure what type I want everything to have, we can forget
# about some of these distinctions.  

# Compile with DOUBLE=yes to get double precision.  The default is
# single precision, since that seems to get good results now.
ifeq ($(DOUBLE),)
   SINGLE:=yes
endif
ifeq ($(SINGLE),)
   C_FLAGS += -DFloat=double -DDouble=double
else
   C_FLAGS += -DFloat=float -DDouble=double
endif
LIBS = -L/usr/X11R6/lib -L$(RASAPIDIR) -lrasmol -lXi -lXext -lX11 -lm

OBJDIR = .
OBJS = $(OBJDIR)/AtomPairInfoState.o $(OBJDIR)/api.o \
       $(OBJDIR)/ccapi.o $(OBJDIR)/caguts.o \
       $(OBJDIR)/ljguts.o $(OBJDIR)/bren2.o $(OBJDIR)/pibond.o \
       $(OBJDIR)/sili_germ.o $(OBJDIR)/bcuint.o $(OBJDIR)/radic.o \
       $(OBJDIR)/mtable.o $(OBJDIR)/vector.o \
       $(OBJDIR)/safe_fgets.o $(OBJDIR)/spgch.o $(OBJDIR)/myassert.o \
       $(OBJDIR)/expand.o $(OBJDIR)/minimize.o $(OBJDIR)/xalloc.o \
       $(OBJDIR)/CCNeighborState.o $(OBJDIR)/vscale.o $(OBJDIR)/rangetree.o

$(OBJDIR)/%.o: %.c
	$(CC) $(INCLUDES) $(C_FLAGS) -c $*.c -o $(OBJDIR)/$*.o

.PHONY: all

# We need libbren.a to build brenfort.  If we don't have
# libbren.a, then don't try to build brenfort.
ifeq ($(wildcard libbren.a),)
   all: brennerfile TAGS
else
   all: brenfort libbrenner.a brennerfile brennerrasmol TAGS
endif

$(OBJDIR)/CCNeighborState.o: CCNeighborState.c CCNeighborState.h

.PHONY: test
new.d old.d coord.d test: brennerfile
	cp coord_safe.d coord.d
	-rm -f xmol.d new.d
	@#If you change the number of timesteps (that is, the number
	@#immediately after the -s), be sure to change the
	@#first number in input.d to be the same, otherwise the movie
	@#lengths will be different and the comparison will fail.
	time ./brennerfile -s200 -t300
	perl ./normalize_d.pl < xmol.d > new.d
	cp coord_safe.d coord.d
	-rm -f xmol.d output.d old.d
	time $(BRENNERMD)
	perl ./normalize_d.pl < xmol.d > old.d
	@# We need the larger error bound if we're doing single precision.
	perl ./diff_d.pl old.d new.d 0.004

# To pass the diff check at the end of timerun, you have to specify
# thermostat 4 and 10 steps between xmol writes in input.d.  That is,
# the first line of input.d can look like this:
# 20 1 4 10  / # steps, # steps between data, thermostat, xmol writes 

# For this test data, stopping the program with the debugger and
# looking at the stack 8 times gave 4 stops in the neighboring code,
# so we spend somewhere around half our time there.  Tim Freeman 29 Jul 2000.

EXE := brennerfile
.PHONY: timerunNoDiff
timerunNoDiff:
	cp coord_safe.d coord.d
	-rm -f xmol.d
	time ./$(EXE) -s200 -t300

.PHONY: timerun
timerun:
	-rm -f new.d
	$(MAKE) timerunNoDiff
	perl ./normalize_d.pl < xmol.d > new.d
	perl ./diff_d.pl new.d correct_timerun.d 0.003

# The diff target generates code in both single precision mode and
# double precision mode, runs the two executables on the same inputs,
# and compares the output xmol.d files and the standard output for
# significant differences.  In the text saved from standard output, if
# there's the string "CMP", we expect a floating point number
# immediately afterward and we'll compare them for approximate
# equality in diff_run.pl.
.PHONY: diff
diff:
	-rm -f new.d xmol.d
	$(MAKE) brennerfile.single DEBUG=yes
	$(MAKE) timerunNoDiff EXE=brennerfile.single> out.single
	perl ./normalize_d.pl < xmol.d > new.d
	@#Can't use make clean here, since that would delete out.single.
	-rm -f old.d xmol.d
	$(MAKE) brennerfile.double DEBUG=yes
	$(MAKE) timerunNoDiff EXE=brennerfile.double> out.double
	perl ./normalize_d.pl < xmol.d > old.d
	perl ./diff_run.pl out.single out.double
	perl ./diff_d.pl old.d new.d 

.PHONY: clean
clean:
	-rm -f libbrenner.a *.o new.d old.d xmol.d output.d coord.d \
	   eigenenergies.d eigenvectors.d dos.d brennerfile ldos.d \
	   load.d max_ke.d overwrite.d pair_energy.d TAGS \
	   out.single out.double
	-rm -rf single double brennerfile.single brennerfile.double

brennerfile: $(OBJS) bren_file.c
	$(CC) -o brennerfile $(INCLUDES) $(C_FLAGS) bren_file.c $(OBJS) -lm

# We always have to rebuild brennerfile.single and ....double so
# declare them PHONY. 
.PHONY: brennerfile.single
brennerfile.single:
	-rm -f brennerfile
	-mkdir single
	$(MAKE) brennerfile SINGLE=yes OBJDIR=single 
	mv brennerfile brennerfile.single

.PHONY: brennerfile.double
brennerfile.double:
	-rm -f brennerfile
	-mkdir double
	$(MAKE) OBJDIR=double brennerfile
	mv brennerfile brennerfile.double

brennerrasmol: $(OBJS) bren_ras.c libbrenner.a
	$(CC) -o brennerrasmol $(INCLUDES) $(C_FLAGS) bren_ras.c -L. -lbrenner $(LIBS)

libbrenner.a : $(OBJS)
	rm -f $@
	ar cq $@ $(OBJS)
	ranlib $@

brenfort.o: brenfort.c brenfort.h brenner.h

brenfort: libbren.a brenfort.o brenffile.c
	$(CC) -o brenfort $(C_FLAGS)  -Wno-missing-prototypes brenfort.o brenffile.c -L. -lbren -lf2c -lm

brenfras: libbren.a brenfort.o brenfras.c
	$(CC) -o brenfras $(C_FLAGS) $(INCLUDES) -Wno-missing-prototypes brenfort.o brenfras.c -L. $(LIBS) -lbren -lf2c -lm

tar:
	cd .. ; tar cf brennerc.tar `cat brennerc/TARGETS`; gzip -9 brennerc.tar
TAGS:
	etags *.c *.h
