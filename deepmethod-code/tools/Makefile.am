
if ENABLE_MPI_BUILD_TARGET
MPI_BUILD_TARGET = deepmethod_openmpi
else
MPI_BUILD_TARGET =
endif

bin_PROGRAMS = \
	deepmethod $(MPI_BUILD_TARGET)

deepmethod_SOURCES = \
	deepmethod.c\
	deepmethod.h \
	xmtranslate.c \
	xmtranslate.h

deepmethod_CFLAGS = \
	-I$(top_builddir)/libdeep\
	-I$(top_builddir)/libxmod

deepmethod_LDFLAGS = \
	-lm

deepmethod_LDADD = \
	$(DEEP_LIBS) \
	$(top_builddir)/libdeep/libdeep.la \
	$(top_builddir)/libxmod/libxmod.la -lm

AM_CPPFLAGS = \
	$(DEEP_CFLAGS) \
	-DDEEPMETHOD_LOCALEDIR=\"${DEEPMETHOD_LOCALEDIR}\"

deepmethod_openmpi_SOURCES = $(deepmethod_SOURCES)

deepmethod_openmpi_LDADD = \
	$(MPI_CLIBS) \
	$(DEEP_LIBS) \
	$(top_builddir)/libdeep/libdeep_openmpi.la \
	$(top_builddir)/libxmod/libxmod.la

deepmethod_openmpi_CFLAGS = \
	-DMPIZE \
	$(deepmethod_CFLAGS) \
	$(MPI_CFLAGS)

## File created by the gnome-build tools
