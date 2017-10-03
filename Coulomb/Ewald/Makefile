# Makefile -  Ewald v2.0
#
# Ju Li <liju99@mit.edu>
#
# Before make, set environment variable SYS to `uname`.

PROG1 = driver1.$(SYS)
PROG2 = driver2.$(SYS)

# Machine options:
# DEC alpha
FCOMPILER.OSF1 = f77
CCOMPILER.OSF1 = cc
OPTFLAG.OSF1 = -O5 -om -ifo -migrate -inline speed -fp_reorder -non_shared
LINK.OSF1 = -lm -L/home/lapack/lib/axp -llapack -lblas
# IBM RS6000
FCOMPILER.AIX = xlf
CCOMPILER.AIX = cc
OPTFLAG.AIX = -O3 -Q -qtune=601
LINK.AIX = -lm -lessl
# SGI IRIX64
FCOMPILER.IRIX64 = f77
CCOMPILER.IRIX64 = cc
OPTFLAG.IRIX64 = -64 -r10000 -TARG:proc=r10000:platform=IP27 -O3 \
-LNO:opt=1:fission=2:fusion=2:cs1=32K:cs2=4M:cs2=1024M:ou_further=3:pwr2=ON\
-OPT:reorg_common=OFF:Olimit=0:div_split=ON:roundoff=3:fast_exp=ON:fast_sqrt=ON
LINK.IRIX64 = -lm
# Linux alpha
CCOMPILER.alpha = ccc
FCOMPILER.alpha = fort
OPTFLAG.alpha = -fast -O3 -fomit-frame-pointer
LINK.alpha = -lcxml -lm
MYFCPPFLAG.alpha = -traditional
# Linux PC
FCOMPILER.Linux = /usr/local/pgcc/bin/gcc
FCOMPILER.Linux = /usr/local/pgcc/bin/g77
OPTFLAG.Linux = -O6 -mpentiumpro -march=pentiumpro -mmx -ffast-math \
	-funroll-loops -inline-functions -s
LINK.Linux = -lm -L$(LIB_PATH) -llapack -lblas
MYFCPPFLAG.alpha = -traditional

MYFLAG = -D_PRINT_EWALD -D_TABULATE_ERFC=1. -D_TABULATE_SINE
#MYFLAG = -D_PRINT_EWALD
CPP = /lib/cpp -P $(MYFCPPFLAG.$(SYS))
CPPFLAG = -D_$(SYS) $(MYFLAG)

FC = $(FCOMPILER.$(SYS))
FFLAGS = $(OPTFLAG.$(SYS)) $(MACHFLAGS.$(SYS))
FHEADER = 
FOBJS1   = driver1.$(SYS).o zhpev.$(SYS).o
FOBJS2   = driver2.$(SYS).o zhpev.$(SYS).o

CC = $(CCOMPILER.$(SYS))
CFLAGS = $(OPTFLAG.$(SYS)) $(MACHFLAGS.$(SYS))
CHEADER = 
COBJS   = ewald.$(SYS).o simple_ewald.$(SYS).o

all: $(PROG1) $(PROG2)
$(PROG1): $(FOBJS1) $(COBJS)
	$(FC) $(FFLAGS) -o $@ $(FOBJS1) $(COBJS) $(LINK.$(SYS))
$(PROG2): $(FOBJS2) $(COBJS)
	$(FC) $(FFLAGS) -o $@ $(FOBJS2) $(COBJS) $(LINK.$(SYS))
run1:	$(PROG1)
	@echo \\n ================== OUTPUT OF $(PROG1) \< con ==================\\n
	@$(PROG1) < con1
run2:	$(PROG2)
	@echo \\n ================== OUTPUT OF $(PROG2) \< con ==================\\n
	@$(PROG2)

clean:
	rm -f *.$(SYS).o *.$(SYS).f *~

cleanall:
	rm -f *.o *.f *~

ewald.$(SYS).o:  ewald.c Makefile
	$(CC) $(CPPFLAG) $(CFLAGS) -o $@ -c ewald.c 

simple_ewald.$(SYS).o:  simple_ewald.c Makefile
	$(CC) $(CPPFLAG) $(CFLAGS) -o $@ -c simple_ewald.c

zhpev.$(SYS).o: zhpev.f
	$(FC) $(FFLAGS) -o $@ -c $?

driver1.$(SYS).o: driver1.$(SYS).f
	$(FC) $(FFLAGS) -o $@ -c $?

driver2.$(SYS).o: driver2.$(SYS).f
	$(FC) $(FFLAGS) -o $@ -c $?
help:
	@echo
	@echo The machine is at `uname -n`, of `uname -m`. 
	@echo The operating system is `uname` `uname -r`.
	@echo In order to make the program for this machine, type
	@echo \\n% "setenv SYS `uname`"\\n
	@echo and then 
	@echo \\n% make \\n 
	@echo you can also add your own flag such as
	@echo \\n% make MYFLAG=\"-D_TABULATE_ERFC=30 -D_TABULATE_SIN -D_ALPHA=2.5\" \\n
	@echo To run the program, type
	@echo \\n% make run1 (or run2)\\n
	@if test `hostname | grep ssibm`; \
	then \
	echo \* Ames IBM RS600 users should compile on ssibm9.ssp.ameslab.gov\\n; \
	fi
