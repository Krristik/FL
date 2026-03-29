dnl Process this file with autoconf to produce a configure script.
dnl Created by Anjuta application wizard.

AC_INIT(deepmethod, 2.3.0)

AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)
AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE

AC_ISC_POSIX
AC_PROG_CC
AM_PROG_CC_STDC
AC_HEADER_STDC




dnl ***************************************************************************
dnl Internatinalization
dnl ***************************************************************************
GETTEXT_PACKAGE=deepmethod
AC_SUBST(GETTEXT_PACKAGE)
AC_DEFINE_UNQUOTED(GETTEXT_PACKAGE,"$GETTEXT_PACKAGE", [GETTEXT package name])
AM_GLIB_GNU_GETTEXT
IT_PROG_INTLTOOL([0.35.0])

DEEPMETHOD_LOCALEDIR=[${datadir}/locale]
AC_SUBST(DEEPMETHOD_LOCALEDIR)

AM_PROG_LIBTOOL



PKG_CHECK_MODULES(DEEP, [glib-2.0 gthread-2.0 gio-2.0])

# Define a macro that is used to parse a --with-mpi parameter
# The macro is named "CHECK_MPICC"
AC_DEFUN([CHECK_MPICC],[

        AC_ARG_WITH(
                mpi,
                [  --with-mpi[=DIR]     OpenMPI installation directory],
                ,
                [with_mpi="/usr/lib64/openmpi"]
        )

        AC_MSG_CHECKING(for OpenMPI installation directory)

        if test "$with_mpi" = "no"; then
                AC_MSG_ERROR( You need to specify the OpenMPI installation directory directory using --with-mpi)
        else
		MPICC="$with_mpi/bin/mpicc"
		MPILIBS="$with_mpi/lib"
                if test -x $MPICC; then
                        MPI_CFLAGS=`$MPICC --showme:compile 2>/dev/null`
                        MPI_LDADD=`$MPICC --showme:link 2>/dev/null`
			MPI_CLIBS="-Wl,--rpath $MPILIBS $MPI_LDADD"
                        AC_MSG_RESULT(OpenMPI installation directory found!)
                else
                        AC_MSG_ERROR( $with_mpi not found.  Check the value you specified with --with-mpi)
                fi
        fi

])

dnl Let people enable/disable the mpi stuff.
AC_ARG_ENABLE(mpi-build-target,
              [  --enable-mpi-build-target  Use mpi-build-target to run on parallel machine [default=no]],
                enable_mpi_build_target="$enableval", enable_mpi_build_target=no)
if test x$enable_mpi_build_target = xyes ; then
	# Now call the CHECK_MPICC macro that was just specified
	CHECK_MPICC
	AC_SUBST(MPI_CFLAGS)
	AC_SUBST(MPI_CLIBS)
fi
AM_CONDITIONAL(ENABLE_MPI_BUILD_TARGET, test x$enable_mpi_build_target = xyes)

AS_IF([test "x$enable_mpi_build_target" == "xyes"], [AC_CONFIG_FILES(libdeep/libdeep_openmpi-1.0.pc)] )

AC_OUTPUT([
Makefile
po/Makefile.in
libdeep/Makefile
libdeep/libdeep-1.0.pc
libxmod/Makefile
libxmod/libxmod-1.0.pc
tools/Makefile
doc/Makefile
nsis/mingwbuild
nsis/DEEPMETHOD.nsis
])

