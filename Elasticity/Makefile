#################################################
# Elasticity: -llapack -lblas -lm               #
#             -lVecMat3 -lVecMat -lScalar -lIO  #
#                                               #
# Jan 29 2001  Ju Li <liju99@mit.edu>           #
#################################################

LIB = Elasticity
HEADERS = ${LIB}.h
STATIC_LIB_FILE = lib${LIB}.a
PROGS = IsotropicStiffnessComplete cubiC hCp T2AbsoluteToObserver \
	T2ObserverToAbsolute
OBJS = ${LIB}.o Stroh.o StraightDislocation.o \
	PBCStaticStraightDislocationDipole.o Homogeneous.o # plural
SRCS = $(addsuffix .c, $(basename ${OBJS}))
USE_SYSTEM_LIB = -lm
USE_MY_LIB = -lVecMat2 -lVecMat3 -lVecMat -lIO -lScalar -lTimer -llapack -lblas
SHELL = sh
GDBTMP = /tmp/${LIB}.gdb

ifndef SYS
 SYS = Linux
endif
ifdef INC_PATH
 MY_INC_DIR = $(INC_PATH)
else
 MY_INC_DIR = $(HOME)/Co/Include
endif
ifdef LIB_PATH
 MY_LIB_DIR = $(LIB_PATH)
else
 MY_LIB_DIR = $(HOME)/Co/Lib
endif
ifdef BIN_PATH
 MY_BIN_DIR = $(BIN_PATH)
else
 MY_BIN_DIR = $(HOME)/Co/Bin
endif
ifdef BAK_PATH
 MY_BAK_DIRS = $(BAK_PATH)
else
 MY_BAK_DIRS = $(HOME)/Co/Bak
endif

CC = cc
LD = f77 # if we link to f77 compiled libraries
AR = ar -r
RANLIB = ranlib
NM = nm -p
SYSTEM_LIB_FLAG =
ifeq ($(SYS), IRIX64)
 CCFLAGS = -64 -O3 -Ofast=ip27 -OPT:IEEE_arithmetic=3
endif
ifeq ($(SYS), IRIX)
 AR = ar -rc
 RANLIB = echo There is no \"ranlib\" to be used on
 CCFLAGS = -O3
 LD = f77 -dont_warn_unused
endif
ifeq ($(SYS), SunOS)
 CCFLAGS = -xO4 -fsimple -dalign -fast
 SHELL = bash
endif
ifeq ($(SYS), alpha)
 CCFLAGS = -fast -O3 -fomit-frame-pointer
 CC = ccc
 LD = fort -nofor_main
 USE_MY_LIB = -lVecMat3 -lVecMat -lIO -lScalar -lTimer -lcxml
endif
ifeq ($(SYS), HPUX)
 CC = gcc
 LD = gcc
 CCFLAGS = -O6 -ffast-math -funroll-loops -finline-functions
 USE_SYSTEM_LIB = -lm -L/opt/fortran90/lib/ -lcl
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
endif
ifeq ($(SYS), Linux)
#  CC = gcc #
#  LD = g77 #
#  CCFLAGS = -Wall -O3 -march=i686 -mmmx -msse -ffast-math \ #
#   -funroll-all-loops -fomit-frame-pointer #
 CC = icc -w
 LD = icc -static
 CCFLAGS = -O2 
 USE_SYSTEM_LIB = -lm
 USE_MY_LIB = -lVecMat3 -lVecMat2 -lVecMat -lIO -lScalar \
	-lTimer -L/opt/intel/mkl70/lib/32/ -lmkl_lapack -lmkl_ia32 \
	-lguide -lsvml -ldl -lpthread
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
endif

CFLAGS = -D_$(SYS) ${CCFLAGS} -I${MY_INC_DIR}
LDFLAGS = ${SYSTEM_LIB_FLAG} ${USE_SYSTEM_LIB} -L${MY_LIB_DIR} ${USE_MY_LIB}

.PHONY: all  i install  t  c clean  b backup  h  r  g  p  a access

all:	install

%.o:    %.c ${HEADERS}  Makefile
	${CC} ${CFLAGS} -c $< -o $@

${STATIC_LIB_FILE}:	${OBJS}
	-@echo "***************************************************"; \
	rm -f $@; \
	${AR} $@ ${OBJS}; \
	${RANLIB} $@; \
	${NM} $@; \
	echo "***************************************************"

i:	install
install: ${STATIC_LIB_FILE} ${HEADERS}
	-@cp -f $< ${MY_LIB_DIR}; \
	chmod a+r ${MY_LIB_DIR}/$<; \
	cp -f ${HEADERS} ${MY_INC_DIR}; \
	dump=`ls -ld ${MY_LIB_DIR}/$< | sed 's/    */  /g'`; \
	month=`echo $$dump | cut -d' ' -f 6`; \
	prev=`echo "$$dump" | sed "s/$$month.*//"`; \
	let l=`echo "$$prev" | wc -c | sed 's/ //g'`-1; \
	machine_name=`hostname | cut -f 1 -d .`:"               "; \
	machine_name=`echo Current time on "$$machine_name" | cut -c1-$$l`; \
	echo "$$machine_name`date '+%b %d %H:%M' | sed 's/    */  /g'`"; \
	a=`echo $$HOME | sed 's/\//~/g'`; \
	for file in ${HEADERS}; do \
	 chmod a+r ${MY_INC_DIR}/$$file; \
	 a=`echo $$HOME | sed 's/\//~/g'`; \
	 ls -l ${MY_INC_DIR}/$$file | sed 's/    */  /g' | sed 's/\//~/g' |\
	 sed "s/$$a//" | sed -e "s/~/\//g" | sed "s/\//~\//"; \
	done; \
	ls -l ${MY_LIB_DIR}/$< | sed 's/    */  /g' | sed 's/\//~/g' | \
	sed "s/$$a//" | sed -e "s/~/\//g" | sed "s/\//~\//"

r:
	-rm -f *_TEST; \
	rm -f ${STATIC_LIB_FILE}

h:
	cp -f ${HEADERS} ${MY_INC_DIR}

c:	clean
clean:
	@echo
	-rm -f *~; \
	rm -f \#*\#; \
	rm -f *.o; \
	rm -f ${STATIC_LIB_FILE}; \
	rm -f *_TEST; \
	rm -f ${CHECKCONST_TEMPLATE}.c; \
	rm -f a.out; \
	rm -f ${PROGS}; \
	rm -f *.m *.m~; \
	rm -f Notes/*~; \
	rm -f Notes/[0-9]; \
	rm -f core
	@echo

b:	backup
backup: clean
	@tmp=/tmp/backupto.tar; \
	echo "taring the directory $(PWD) ..."; \
	dir=`echo $(PWD) | sed -e "s/.*\///g"`; \
	cd ..; \
	echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; \
	tar cvf $$tmp $$dir/* \
	`ls -d $$dir/.* | sed -e "s/$$dir\/\.\.$$//" | \
	sed -e "s/$$dir\/\.$$//"`; \
	echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; \
	echo compressing $$tmp with gzip ...; \
	gzip --verbose -c $$tmp > $$tmp.gz; \
	echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; \
	nbytes=`ls -l $$tmp.gz | sed -e 's/  */ /g' | cut -d' ' -f 5`; \
	echo $$tmp.gz = $$nbytes bytes.; \
	for backup in $(MY_BAK_DIRS); do \
	 if [ -d $$backup ]; then \
	  echo ' '; \
	  echo "copying $$tmp.gz to $$backup/$$dir.tar.gz ..."; \
	  cp -i $$tmp.gz $$backup/$$dir.tar.gz; \
	 fi; \
	done; \
	echo

g:
	@export CC=gcc; export LD=g77; export CCFLAGS=-g; \
	make --environment-overrides --no-print-directory \
	`echo "${MAKECMDGOALS}" | sed -e 's/^g //' | sed -e 's/^g$///' `; \
	exit 99

p:	${PROGS}
${PROGS}:${SRCS} ${HEADER}  ${ADD_DEPENDENCE}
	@for obj in $(OBJS); do \
	 base=`echo $$obj | sed -e 's/\.o$$//g'`; \
	 grep _$@ $$base.c > /dev/null; \
	 if [ $$? = "0" ]; then \
	  ${CC} -D_$@ ${CFLAGS} -c $$base.c -o "$$base"_$@".o"; \
	  if [ $$? = "0" ]; then \
	   objs="$$base"_$@".o $$objs"; else exit 99; fi; \
	 else \
	  make --no-print-directory $$base.o > /dev/null; \
	  objs="$$objs $$base.o"; \
	 fi; \
	done; \
	${LD} ${CFLAGS} $$objs -o $@ ${LDFLAGS}; \
	cp -f $@ ${MY_BIN_DIR}; chmod u+x ${MY_BIN_DIR}/$@; \
	ls -l ${MY_BIN_DIR}/$@ | sed -e 's/   */  /g'

%_TEST: $(SRCS) ${HEADER}
	for obj in $(OBJS); do \
	 base=`echo $$obj | sed -e 's/\.o$$//g'`; \
	 grep _$@ $$base.c > /dev/null; \
	 if [ $$? = "0" ]; then \
	  ${CC} -D_$@ ${CFLAGS} -c $$base.c -o "$$base"_$@".o"; \
	  if [ $$? = "0" ]; then \
	   objs="$$base"_$@".o $$objs"; else exit 99; fi; \
	 else \
	  make --no-print-directory $$base.o > /dev/null; \
	  objs="$$objs $$base.o"; \
	 fi; \
	done; \
	${LD} ${CFLAGS} $$objs -o $@ ${LDFLAGS}; \
	chmod u+x $@

#	Selective test of certain function, say "make Mprint"
%:	%.f
%:	%.c
%:	%.o
%:	%_TEST ${HEADER}  Makefile
	-@echo ' '; \
	echo ====== Running \"$@_TEST\" ======; \
	if [ "${CCFLAGS}" = "-g" ] ; then \
	 echo run > ${GDBTMP}; gdb ./$@_TEST -x ${GDBTMP}; \
	else \
	 ./$@_TEST; \
	fi; \
	echo ==== Done running \"$@_TEST\" ===; \
	echo ' '

.PRECIOUS:	%_TEST
Makefile: ;

a:	access
access:
	-@chmod -R 700 *; \
	chmod 744 *.h; \
	if [ -f ${STATIC_LIB_FILE} ]; then chmod 744 ${STATIC_LIB_FILE}; fi; \
	ls -la
