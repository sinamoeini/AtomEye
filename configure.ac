AC_INIT([AtomeEye],[0.0])
AM_INIT_AUTOMAKE

AC_PROG_RANLIB
AC_LANG(C)
AC_PROG_CC
AM_PROG_CC_C_O

#m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CANONICAL_HOST

sys=Linux

case "${host_os}" in
    linux*)
    AC_DEFINE(_Linux)
        ;;
    cygwin*|mingw*)
    AC_DEFINE(_CYGWIN)
        ;;
    darwin*)
    AC_DEFINE(_Darwin,)
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac



# A
#dyanamic
#-lreadline -lncurses
#-lXpm -lXext -lX11
#-lpng -ljpeg -lz
#-lgsl -lgslcblas
#-llapack -lblas
#-lpthread

#static
#-lreadline -lncurses
#-lXpm -lXext -lX11 -lxcb -lXau
#-lpng -ljpeg -lz
#-lgsl -lgslcblas
#-llapack -lblas -lgfortran -lm
#-lpthread
#-ldl

# Utils
#daynamic
#-lgsl -lgslcblas
#-llapack

#static:
#-lgsl -lgslcblas
#-llapack -lblas -lgfortran
#-lm


# A(static)
AC_CHECK_LIB([dl],[dlsym],[],[AC_MSG_ERROR([unable to find dl library])])
# A
AC_CHECK_LIB([pthread],[pthread_mutex_lock],[],[AC_MSG_ERROR([unable to find dl library])])
# A(static) and Utils(static)
AC_CHECK_LIB([m],[cos],[],[AC_MSG_ERROR([unable to find m library])])
# A(static) and Utils(static)
AC_CHECK_LIB([gfortran],[_gfortran_compare_string],[],[AC_MSG_ERROR([unable to find gfortran library])])
# A and Utils
AC_CHECK_LIB([blas],[dgemm_],[],[AC_MSG_ERROR([unable to find blas library])])
# A and Utils
AC_CHECK_LIB([lapack],[dgetrf_],[],[AC_MSG_ERROR([unable to find lapack library])])
# A and Utils
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR([unable to find gslcblas library])])
# A and Utils
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([unable to find gsl library])])
# A
AC_CHECK_LIB([z],[gzerror],[],[AC_MSG_ERROR([unable to find z library])])
# A
AC_CHECK_LIB([jpeg],[jpeg_stdio_dest],[],[AC_MSG_ERROR([unable to find jpeg library])])
# A
AC_CHECK_LIB([png],[png_create_write_struct],[],[AC_MSG_ERROR([unable to find png library])])
# A (static)
AC_CHECK_LIB([Xau],[XauDisposeAuth],[],[AC_MSG_ERROR([unable to find Xau library])])
# A (static)
AC_CHECK_LIB([xcb],[xcb_parse_display],[],[AC_MSG_ERROR([unable to find xcb library])])
# A
AC_CHECK_LIB([X11],[XOpenDisplay],[],[AC_MSG_ERROR([unable to find X11 library])])
# A
AC_CHECK_LIB([Xext],[DPMSQueryExtension],[],[AC_MSG_ERROR([unable to find Xext library])])
# A
AC_CHECK_LIB([Xpm],[main],[],[AC_MSG_ERROR([unable to find Xpm library])])
# A
AC_CHECK_LIB([ncurses],[tgoto],[],[AC_MSG_ERROR([unable to find ncurses library])])
# A
AC_CHECK_LIB([readline],[readline],[],[AC_MSG_ERROR([unable to find readline library])])






# A Utils
AC_CHECK_HEADER(gsl/gsl_spline.h,[], [AC_MSG_ERROR([unable to find header gsl/gsl_spline.h])],[])
# A Utils
AC_CHECK_HEADER(gsl/gsl_errno.h,[], [AC_MSG_ERROR([unable to find header gsl/gsl_errno.h])],[])
# A
AC_CHECK_HEADER(X11/Xlib.h,[], [AC_MSG_ERROR([unable to find header X11/Xlib.h])],[])



# A
AC_CHECK_HEADER(X11/Xutil.h,[], [AC_MSG_ERROR([unable to find header X11/Xutil.h])],[])
# A
AC_CHECK_HEADER(X11/xpm.h,[], [AC_MSG_ERROR([unable to find header X11/xpm.h])],[])
# A
AC_CHECK_HEADER(X11/XKBlib.h,[], [AC_MSG_ERROR([unable to find header X11/XKBlib.h])],[])
# A
AC_CHECK_HEADER(zlib.h,[], [AC_MSG_ERROR([unable to find header zlib.h])],[])
# A
AC_CHECK_HEADER(png.h,[], [AC_MSG_ERROR([unable to find header png.h])],[])
# A
AC_CHECK_HEADER(jpeglib.h,[], [AC_MSG_ERROR([unable to find header jpeglib.h])],[])








AC_CONFIG_FILES([Makefile Timer/Makefile IO/Makefile Scalar/Makefile VecMat/Makefile VecMat3/Makefile Atoms/Makefile AX/Makefile A/Makefile])

AC_OUTPUT

