#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([mfold],[3.6],[zukerm@alum.mit.edu],[mfold],[www.unafold.org])
AC_CONFIG_SRCDIR([src/auxgen.f])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

AC_CANONICAL_BUILD

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_F77
AC_PROG_RANLIB
AC_PROG_CPP
AM_PROG_CC_C_O

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([gd], [gdImageCreate])
AC_CHECK_LIB([jpeg], [jpeg_CreateCompress])
AC_CHECK_LIB([png], [png_init_io])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([gd.h limits.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([bzero floor sqrt strchr strspn strstr])

AC_CONFIG_FILES([Makefile
		c++/Makefile
		src/Makefile
		dat/Makefile
		scripts/Makefile
		scripts/filter-sort
		scripts/h-num
		scripts/reformat-seq.sh
		scripts/mfold
		scripts/mfold_quik])

AC_OUTPUT
