
if ENABLE_MPI_BUILD_TARGET
MPI_BUILD_LIBTARGET = libdeep_openmpi.la
MPI_INSTALL_PKGCONF = libdeep_openmpi-1.0.pc
else
MPI_BUILD_LIBTARGET =
MPI_INSTALL_PKGCONF =
endif

lib_LTLIBRARIES = \
	libdeep.la $(MPI_BUILD_LIBTARGET)

libdeep_la_SOURCES = \
	dpindivid.c \
	dpindivid.h \
	dppopulation.c \
	dppopulation.h \
	dpmutation.c \
	dpmutation.h \
	dparchive.c \
	dparchive.h \
	dprecombination.c \
	dprecombination.h \
	dpevaluation.c \
	dpevaluation.h \
	dpdeep.c \
	dpdeep.h \
	dpopt.c \
	dpopt.h \
	dptarget.c \
	dptarget.h \
	dploop.c \
	dploop.h \
	dpsettings.c \
	dpsettings.h \
	dposda.c \
	dposda.h

libdeep_la_LIBADD = \
	$(DEEP_LIBS) \
	-lm

AM_CPPFLAGS = \
	$(DEEP_CFLAGS) \
	-DG_LOG_DOMAIN=\"DeepMethod\"

libdeep_openmpi_la_SOURCES = $(libdeep_la_SOURCES)

libdeep_openmpi_la_CFLAGS = \
	-DMPIZE \
	$(libdeep_la_CFLAGS) \
	$(MPI_CFLAGS)

libdeepinclude_HEADERS = \
	dpdeep.h \
	dpevaluation.h \
	dpindivid.h \
	dploop.h \
	dpmutation.h \
	dpopt.h \
	dposda.h \
	dppopulation.h \
	dprecombination.h \
	dpsettings.h \
	dptarget.h

EXTRA_DIST = \
	$(libdeepinclude_HEADERS)\
	$(pkgconfig_DATA) 

libdeepincludedir = \
	$(includedir)/libdeep

pkgconfig_DATA = \
	libdeep-1.0.pc $(MPI_INSTALL_PKGCONF)

pkgconfigdir = \
	$(libdir)/pkgconfig

## File created by the gnome-build tools

