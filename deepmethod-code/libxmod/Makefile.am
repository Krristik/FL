
lib_LTLIBRARIES = \
	libxmod.la

libxmod_la_SOURCES = \
	xmmodel.c   \
	xmmodel.h

libxmod_la_LIBADD = \
	$(DEEP_LIBS) -lm

AM_CPPFLAGS = \
	$(DEEP_CFLAGS) \
	-DG_LOG_DOMAIN=\"XmModel\"

libxmodinclude_HEADERS = \
	xmmodel.h

EXTRA_DIST = \
	$(libxmodinclude_HEADERS)\
	$(pkgconfig_DATA) 

libxmodincludedir = \
	$(includedir)/libxmod

pkgconfig_DATA = \
	libxmod-1.0.pc

pkgconfigdir = \
	$(libdir)/pkgconfig

## File created by the gnome-build tools

