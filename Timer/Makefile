######################################
# libTimer:                          #
#                                    #
# Time and Schedule Manager          #
#                                    #
# Nov 19 1999 Ju Li <liju99@mit.edu> #
######################################

LIB = Timer
HEADER = ${LIB}.h # singular
STATIC_LIB_FILE = lib${LIB}.a
PROGS = datenum datevec
OBJS = Timer.o chronometer.o taskwatch.o linear_scheduler.o fools.o \
	Gregorian.o # plural
SRCS = $(addsuffix .c, $(basename ${OBJS}))
USE_SYSTEM_LIB = -lm
USE_MY_LIB =
SHELL = sh
GDBTMP = /tmp/${LIB}.gdb
EXESUFFIX =

ifndef SYS
 SYS = Linux
endif
ifdef INC_PATH
 MY_INC_DIR = $(INC_PATH)
else
 MY_INC_DIR = ~/Co/Include
endif
ifdef LIB_PATH
 MY_LIB_DIR = $(LIB_PATH)
else
 MY_LIB_DIR = ~/Co/Lib
endif
ifdef BIN_PATH
 MY_BIN_DIR = $(BIN_PATH)
else
 MY_BIN_DIR = $(HOME)/Co/Bin
endif

CC = cc
LD = f77 # if we link to f77 compiled libraries
AR = ar -r
RANLIB = ranlib
NM = nm -p
DATE =
SYSTEM_LIB_FLAG =
BACKUPDIRS = /asm/home/Moon/Bak
ifeq ($(SYS), IRIX64)
 CCFLAGS = -64 -r10000 -TARG:proc=r10000:platform=IP27 -O3 -diag_suppress 1197\
  -LNO:opt=1:fission=2:fusion=2:cs1=32K:cs2=4M:cs2=1024M:ou_further=3:pwr2=ON\
  -OPT:Olimit=0:div_split=ON:roundoff=3:fast_exp=ON:fast_sqrt=ON -s
 BACKUPDIRS = /home/Moon/Bak /usr/people/liju99/Bak
endif
ifeq ($(SYS), IRIX)
 AR = ar -rc
 RANLIB = echo There is no \"ranlib\" to be used on
 CCFLAGS = -O3
endif
ifeq ($(SYS), SunOS)
 CCFLAGS = -xO4 -fsimple -dalign -fast
 SHELL = bash
endif
ifeq ($(SYS), alpha)
 CCFLAGS = -fast -O3 -fomit-frame-pointer
 CC = ccc
 LD = fort -nofor_main
endif
ifeq ($(SYS), HPUX)
 CC = gcc
 LD = gcc
 CCFLAGS = -O6 -ffast-math -funroll-loops -finline-functions
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
endif
ifeq ($(SYS), OSF1)
 CC = cc
 LD = cc
 CCFLAGS = -O5
 USE_SYSTEM_LIB = -lm -lfor
 USE_MY_LIB = -ldb
endif
ifeq ($(SYS), CYGWIN)
 CC = gcc
 LD = g77
 CCFLAGS = -O6 -ffast-math -funroll-loops -finline-functions
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
 BACKUPDIRS = /asm/home/Moon/Bak /Home/Sun/Bak
 EXESUFFIX = .exe
endif
ifeq ($(SYS), Linux)
#  CC = gcc #
#  LD = g77 #
#  CCFLAGS = -Wall -O3 -march=i686 -mmmx -msse -ffast-math \ #
#   -funroll-all-loops -fomit-frame-pointer  #
 CC = icc -w
 LD = icc -static
 CCFLAGS = -O2 -axKWNBP -mp1
#  CCFLAGS = -O2 #
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
 BACKUPDIRS = /home/Moon/Bak /Home/Sun/Bak
endif
ifeq ($(SYS), Darwin)
 CC = gcc
 LD = f77 
 CCFLAGS = -O3 -fomit-frame-pointer -fno-schedule-insns \
	-fschedule-insns2
endif
ifeq ($(SYS), DarwinIntel386)
 CC = gcc
 LD = cc 
 CCFLAGS = -O3 -fomit-frame-pointer -fno-schedule-insns \
	-fschedule-insns2
endif


CFLAGS = -D_$(SYS) ${CCFLAGS} -I${MY_INC_DIR}
LDFLAGS = ${SYSTEM_LIB_FLAG} ${USE_SYSTEM_LIB} -L${MY_LIB_DIR} ${USE_MY_LIB}

.PHONY: all  i install  c clean  b backup  r  g  h  a access

all:	install

%.o:	%.c ${HEADER}
	${CC} ${CFLAGS} -c $< -o $@

${STATIC_LIB_FILE}:	${OBJS}
	-@echo "***************************************************"; \
	rm -f $@; \
	${AR} $@ ${OBJS}; \
	${RANLIB} $@; \
	${NM} $@; \
	echo "***************************************************"

i:	install
install:	${STATIC_LIB_FILE}
	-@cp -f $< ${MY_LIB_DIR}/; \
	chmod a+r ${MY_LIB_DIR}/$<; \
	cp -f ${HEADER} ${MY_INC_DIR}/; \
	chmod a+r ${MY_INC_DIR}/${HEADER}; \
	dump=`ls -ld ${MY_INC_DIR}/${HEADER} | sed 's/    */  /g'`; \
	month=`echo $$dump | cut -d' ' -f 6`; \
	prev=`echo "$$dump" | sed "s/$$month.*//"`; \
	ww=`echo "$$prev" | wc -c | sed 's/ //g'`; \
	l=`expr $$ww - 1`; \
	machine_name=`hostname | cut -f 1 -d .`:"               "; \
	machine_name=`echo Current time on "$$machine_name" | cut -c1-$$l`; \
	echo "$$machine_name`date '+%b %d %H:%M' | sed 's/    */  /g'`"; \
	a=`echo $$HOME | sed 's/\//~/g'`; \
	ls -l ${MY_INC_DIR}/${HEADER} | sed 's/    */  /g' | sed 's/\//~/g' |\
	sed "s/$$a//" | sed -e "s/~/\//g" | sed "s/\//~\//"; \
	ls -l ${MY_LIB_DIR}/$< | sed 's/    */  /g' | sed 's/\//~/g' | \
	sed "s/$$a//" | sed -e "s/~/\//g" | sed "s/\//~\//"

h:
	cp -f ${HEADER} ${MY_INC_DIR}/

r:	
	rm -f *_TEST${EXESUFFIX};

c:	clean
clean:
	@echo
	-rm -fr *~; \
	rm -fr \#*\#; \
	rm -fr *.o; \
	rm -fr ${STATIC_LIB_FILE}; \
	rm -fr *_TEST${EXESUFFIX}; \
	rm -fr a.out; \
	rm -fr core
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
	for backup in $(BACKUPDIRS); do \
	 if [ -d $$backup ]; then \
	  echo ' '; \
	  echo "copying $$tmp.gz to $$backup/$$dir.tar.gz ..."; \
	  cp -i $$tmp.gz $$backup/$$dir.tar.gz; \
	 fi; \
	done; \
	echo

g:
	@export CC=gcc; export LD=g77; export CCFLAGS=-g; \
	gmake --environment-overrides --no-print-directory \
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
	${LD} ${CFLAGS} $$objs -o $@${EXESUFFIX} ${LDFLAGS}; \
	cp -f $@${EXESUFFIX} ${MY_BIN_DIR}; \
	chmod u+x ${MY_BIN_DIR}/$@${EXESUFFIX}; \
	ls -l ${MY_BIN_DIR}/$@${EXESUFFIX} | sed -e 's/   */  /g'

%_TEST: $(SRCS) ${HEADER}
	for obj in $(OBJS); do \
	 base=`echo $$obj | sed -e 's/\.o$$//g'`; \
	 grep _$@ $$base.c > /dev/null; \
	 if [ $$? = "0" ]; then \
	  ${CC} -D_$@ ${CFLAGS} -c $$base.c -o "$$base"_$@".o"; \
	  if [ $$? = "0" ]; then \
	   objs="$$base"_$@".o $$objs"; else exit 99; fi; \
	 else \
	  gmake --no-print-directory $$base.o > /dev/null; \
	  objs="$$objs $$base.o"; \
	 fi; \
	done; \
	${LD} ${CFLAGS} $$objs -o $@${EXESUFFIX} ${LDFLAGS}; \
	chmod u+x $@${EXESUFFIX}

#	Selective test of certain function, say "gmake Mprint"
%:	%.f
%:	%.c
%:	%.o
%:	%_TEST
	-@echo ' '; \
	echo ============ Running \"$@_TEST${EXESUFFIX}\" ============; \
	if [ "${CCFLAGS}" = "-g" ] ; then \
	 echo run > ${GDBTMP}; gdb ./$@_TEST${EXESUFFIX} -x ${GDBTMP}; \
	else \
	 ./$@_TEST${EXESUFFIX}; \
	fi; \
	echo ========== Done running \"$@_TEST${EXESUFFIX}\" =========; \
	echo ' '

.PRECIOUS:	%_TEST${EXESUFFIX}
Makefile: ;

a:	access
access:
	-@chmod -R 700 *; \
	chmod 744 *.h; \
	if [ -f ${STATIC_LIB_FILE} ]; then chmod 744 ${STATIC_LIB_FILE}; fi; \
	ls -la
