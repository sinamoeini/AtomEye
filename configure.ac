AC_INIT([AtomeEye],[0.0])
AM_INIT_AUTOMAKE
AC_PROG_RANLIB
AC_LANG(C)
AC_PROG_CC

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



AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[AC_MSG_ERROR([unable to find gslcblas library])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([unable to find gsl library])])
AC_CHECK_LIB([X11],[XOpenDisplay],[],[AC_MSG_ERROR([unable to find X11 library])])
AC_CHECK_LIB([Xpm],[main],[],[AC_MSG_ERROR([unable to find Xpm library])])
AC_CHECK_LIB([Xext],[DPMSQueryExtension],[],[AC_MSG_ERROR([unable to find Xext library])])
AC_CHECK_LIB([jpeg],[jpeg_stdio_dest],[],[AC_MSG_ERROR([unable to find jpeg library])])
AC_CHECK_LIB([png],[png_create_write_struct],[],[AC_MSG_ERROR([unable to find png library])])
AC_CHECK_LIB([z],[gzerror],[],[AC_MSG_ERROR([unable to find z library])])
AC_CHECK_LIB([lapack],[dgetrf_],[],[AC_MSG_ERROR([unable to find lapack library])])
AC_CHECK_LIB([blas],[dgemm_],[],[AC_MSG_ERROR([unable to find blas library])])
AC_CHECK_LIB([readline],[readline],[],[AC_MSG_ERROR([unable to find readline library])])







AC_CHECK_HEADER(gsl/gsl_spline.h,[], [AC_MSG_ERROR([unable to find header gsl/gsl_spline.h])],[])
AC_CHECK_HEADER(gsl/gsl_errno.h,[], [AC_MSG_ERROR([unable to find header gsl/gsl_errno.h])],[])
AC_CHECK_HEADER(X11/Xlib.h,[], [AC_MSG_ERROR([unable to find header X11/Xlib.h])],[])



AC_CHECK_HEADER(X11/Xutil.h,[], [AC_MSG_ERROR([unable to find header X11/Xutil.h])],[])
AC_CHECK_HEADER(X11/xpm.h,[], [AC_MSG_ERROR([unable to find header X11/xpm.h])],[])
AC_CHECK_HEADER(X11/XKBlib.h,[], [AC_MSG_ERROR([unable to find header X11/XKBlib.h])],[])
AC_CHECK_HEADER(zlib.h,[], [AC_MSG_ERROR([unable to find header zlib.h])],[])
AC_CHECK_HEADER(png.h,[], [AC_MSG_ERROR([unable to find header png.h])],[])
AC_CHECK_HEADER(jpeglib.h,[], [AC_MSG_ERROR([unable to find header jpeglib.h])],[])








AC_CONFIG_FILES([Makefile Timer/Makefile IO/Makefile Scalar/Makefile VecMat/Makefile VecMat3/Makefile Atoms/Makefile AX/Makefile A/Makefile])

AC_OUTPUT
