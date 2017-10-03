#######################################
# Brenner Potential v0.03:            #
# Apr 7, 2001  Ju Li <liju99@mit.edu> #
#######################################

# http://www.rahul.net/pcm/brenner/ #
SRC_DIR = brennerc
# runtime parameters: see Dat/liju.txt #
DAT_DIR = Dat
# additional declarations to have access to more internals #
ADD_DIR = Add
# modify the maximum number of atoms #
MAX_ATOMS = 60000
COMPONENTS = AtomPairInfoState api ccapi caguts ljguts bren2 pibond \
	sili_germ bcuint radic mtable vector safe_fgets spgch myassert \
	expand minimize xalloc CCNeighborState vscale rangetree rmax
OBJS = $(addsuffix .o, $(COMPONENTS))

SHELL = sh
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
# libBrenner reads spline tables at runtime #
# ifdef RUN_PATH #
#  MY_RUN_DIR = $(RUN_PATH)/Brenner #
# else #
#  MY_RUN_DIR = ~/Co/Runtime/Brenner #
# endif #
MY_RUN_DIR = /Co/Runtime/Brenner
ifdef BAK_PATH
 MY_BAK_DIRS = $(BAK_PATH)
else
 MY_BAK_DIRS = $(HOME)/Co/Bak
endif

CC = cc
CPP = /lib/cpp -P -C
AR = ar -r
RANLIB = ranlib
NM = nm -p

# for some compilers, void type has no size #
PATCH = cp $(SRC_DIR)/rangetree.c $(SRC_DIR)/rangetree.C; \
	sed 's/void \*mem/char \*mem/g' $(SRC_DIR)/rangetree.C \
	> $(SRC_DIR)/rangetree.c; \
	cp $(SRC_DIR)/brenner.h $(SRC_DIR)/brenner.H; \
	sed 's/\#define MAX_ATOMS .*/\#define MAX_ATOMS '${MAX_ATOMS}'/g' \
	$(SRC_DIR)/brenner.H > $(SRC_DIR)/brenner.h; \
	echo 'extern Float rmax[4][4];' > $(SRC_DIR)/rmax.h; \
	echo 'Float rmax[4][4];' > $(SRC_DIR)/rmax.c; \
	cp $(SRC_DIR)/api.h $(SRC_DIR)/api.H; \
	sed 's/.*static inline.*//g' \
	$(SRC_DIR)/api.H > $(SRC_DIR)/api.h;

BREN2.NOPBC = cp $(SRC_DIR)/bren2.cc $(SRC_DIR)/bren2.c
# make_neighbor_lists_tree() of rangetree.c doesn't work under PBC #
BREN2.PBC = sed 's/make_neighbor_lists_tree(info.*/n2_neighbor_list(info);/g' \
	$(SRC_DIR)/bren2.cc > $(SRC_DIR)/bren2.c

CCFLAGS = -Dmy_exit=exit -DNDEBUG

ifeq ($(SYS), IRIX64)
 CCFLAGS += -64 -r10000 -TARG:proc=r10000:platform=IP27 -O3 \
  -LNO:opt=1:fission=2:fusion=2:cs1=32K:cs2=4M:cs2=1024M:ou_further=3:pwr2=ON \
  -OPT:Olimit=0:div_split=ON:roundoff=3:fast_exp=ON:fast_sqrt=ON -s \
  -Dinline= -D__PRETTY_FUNCTION__=__FILE__
endif
ifeq ($(SYS), IRIX)
 AR = ar -rc
 RANLIB = echo There is no \"ranlib\" to be used on
 CCFLAGS += -O3 -Dinline= -D__PRETTY_FUNCTION__=__FILE__
 LD = f77 -dont_warn_unused
endif
ifeq ($(SYS), SunOS)
 CCFLAGS += -xO4 -xCC -fsimple -dalign -fast \
	-Dinline= -D__PRETTY_FUNCTION__=__FILE__
 SHELL = bash
endif
ifeq ($(SYS), alpha)
 CCFLAGS += -fast -O3 -fomit-frame-pointer -D__PRETTY_FUNCTION__=__FILE__
 CC = ccc
 LD = fort -nofor_main
endif
ifeq ($(SYS), Linux)
 CC = gcc # rangetree.c segfaults under pgcc: compiler bug #
 LD = g77
 CCFLAGS += -Wall -O6 -march=pentiumpro -ffast-math
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
endif

.PHONY:	all  c clean  b backup  i install  patch \
	dp  dp.pbc  sp  sp.pbc
all:	dp  dp.pbc  sp  sp.pbc

# making sure that $(SRC_DIR)/bren2.cc is an untouched backup copy #
patch:	$(SRC_DIR)/bren2.cc
	${PATCH}
# nobody should touch $(SRC_DIR)/bren2.cc after it's created !! #
$(SRC_DIR)/bren2.cc:
	cp $(SRC_DIR)/bren2.c $(SRC_DIR)/bren2.cc

dp:	c patch
	@echo Making double precision, PBC-less library...;
	${BREN2.NOPBC};
	@export FLAGS="-DFloat=double -DDouble=double -DINFINITE_CUBE"; \
	$(MAKE) --no-print-directory libBrenner.$@.a Brenner.$@.h install c
dp.pbc: patch
	@echo Making double precision, PBC library...;
	${BREN2.PBC};
	@export FLAGS="-DFloat=double -DDouble=double"; \
	$(MAKE) --no-print-directory libBrenner.$@.a Brenner.$@.h install c
sp:	patch
	@echo Making single precision, PBC-less library...;
	${BREN2.NOPBC};
	@export FLAGS="-DFloat=float -DDouble=double -DINFINITE_CUBE"; \
	$(MAKE) --no-print-directory libBrenner.$@.a Brenner.$@.h install c
sp.pbc: patch
	@echo Making single precision, PBC library...;
	${BREN2.PBC};
	@export FLAGS="-DFloat=float -DDouble=double"; \
	$(MAKE) --no-print-directory libBrenner.$@.a Brenner.$@.h install c

%.a:	${OBJS}
	-@echo "***************************************************"; \
	rm -f $@; \
	${AR} $@ ${OBJS}; \
	${RANLIB} $@; \
	${NM} $@; \
	echo "***************************************************"

%.o:    $(SRC_DIR)/%.c
	cd $(SRC_DIR); ${CC} \
	${CCFLAGS} \
	${FLAGS} \
	-c ../$< -o ../$@

# OK, pick only what the end user needs to see #
Brenner.%.h: $(SRC_DIR)/brenner.h
	@echo
	@echo Assembling a neat header file...
	cd $(SRC_DIR); \
	grep typedef vector.h > vector.H;	echo ' ' >> vector.H; \
	grep typedef bool.h   > bool.H;		echo ' ' >> bool.H; \
	grep typedef NeighborState.h > NeighborState.H; \
	echo ' '					 > NeighborState.H; \
	grep typedef IntPair.h > IntPair.H;	echo ' ' > IntPair.H; \
	grep typedef brenc.h > brenc.H;		echo ' ' > brenc.H; \
	cat ../$< | \
	sed 's/^#include <.*>.*//g' | \
	sed 's/vector\.h\"/vector\.H\"/' | \
	sed 's/bool\.h\"/bool\.H\"/' | \
	sed 's/NeighborState\.h\"/NeighborState\.H\"/' | \
	sed 's/IntPair\.h\"/IntPair\.H\"/' | \
	sed 's/brenc\.h\"/brenc\.H\"/' | \
	${CPP} \
	${FLAGS} | \
	sed '/^ *$$/d' > ../$@; \
	cat ../${ADD_DIR}/*.h >> ../$@; \
	echo '#define' BRENNER_MAX_ATOMS ${MAX_ATOMS} >> ../$@; \
	echo '#define' BRENNER_DAT \"${MY_RUN_DIR}\" >> ../$@;
	@echo

i:	install
install:
	-@dump=`ls -ld ${MY_INC_DIR}/. | sed 's/    */  /g'`; \
	month=`echo $$dump | cut -d' ' -f 6`; \
	prev=`echo "$$dump" | sed "s/$$month.*//"`; \
	let l=`echo "$$prev" | wc -c | sed 's/ //g'`-1; \
	machine_name=`hostname | cut -f 1 -d .`:"               "; \
	machine_name=`echo Current time on "$$machine_name" | cut -c1-$$l`; \
	time=`date '+%b %d %H:%M' | sed 's/    */  /g' | sed 's/ 0/  /g'`; \
	echo "$$machine_name$$time"; \
	curpath=`echo $$HOME | sed 's/\//~/g'`; \
	for header in *.h; do \
	 cp -f $$header ${MY_INC_DIR}/; \
	 chmod a+r ${MY_INC_DIR}/$$header; \
	 ls -l ${MY_INC_DIR}/$$header | sed 's/    */  /g' | sed 's/\//~/g' |\
	 sed "s/$$curpath//" | sed -e "s/~/\//g" | sed "s/\//~\//"; \
	done; \
	for lib in *.a; do \
	 cp -f $$lib ${MY_LIB_DIR}/; \
	 chmod a+r ${MY_LIB_DIR}/$$lib; \
	 ls -l ${MY_LIB_DIR}/$$lib | sed 's/    */  /g' | sed 's/\//~/g' | \
	 sed "s/$$curpath//" | sed -e "s/~/\//g" | sed "s/\//~\//"; \
	done; \
	if ( ! [ -d ${MY_RUN_DIR} ] ); then \
	 mkdir ${MY_RUN_DIR}; \
	fi; \
	cd ${DAT_DIR}; \
	for run in *.d; do \
	 cp -f $$run ${MY_RUN_DIR}/; \
	 chmod a+r ${MY_RUN_DIR}/$$run; \
	 ls -l ${MY_RUN_DIR}/$$run | \
	 sed 's/    */  /g' | sed 's/\//~/g' | \
	 sed "s/$$curpath//" | sed -e "s/~/\//g" | sed "s/\//~\//"; \
	done;

c:	clean
clean:
	@echo
	-rm -f *~; \
	rm -f \#*\#; \
	rm -f .\#*; \
	rm -f *.o; \
	rm -f *.a; \
	rm -f *.h
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
