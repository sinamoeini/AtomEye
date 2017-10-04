bin_PROGRAMS = A
A_SOURCES = A.c A.h LeastSquareStrain.c geo.c icon.c info.c primitives.c rcut_patch.c scratch.c utils.c viewport.c xtal_shift.c  

AM_CFLAGS = -I$(top_srcdir)/Timer -I$(top_srcdir)/IO -I$(top_srcdir)/Scalar -I$(top_srcdir)/VecMat -I$(top_srcdir)/VecMat3 -I$(top_srcdir)/Atoms -I$(top_srcdir)/AX

A_LDADD = $(top_srcdir)/Timer/libTimer.a $(top_srcdir)/IO/libIO.a  $(top_srcdir)/Scalar/libScalar.a  $(top_srcdir)/VecMat/libVecMat.a  $(top_srcdir)/VecMat3/libVecMat3.a  $(top_srcdir)/Atoms/libAtoms.a  $(top_srcdir)/AX/libAX.a -lm -lX11 -lXpm -lXext -lpng -lz -ljpeg -lreadline -lncurses -lpthread -llapack -lblas




