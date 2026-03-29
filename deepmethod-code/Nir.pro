TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c11
CONFIG += link_pkgconfig
PKGCONFIG += glib-2.0
PKGCONFIG += gio-2.0
DEFINES += DEEPMETHOD_LOCALEDIR=0
INCLUDEPATH += libdeep
INCLUDEPATH += libxmod
INCLUDEPATH += tools

SOURCES += \
    libdeep/dparchive.c \
    libdeep/dpdeep.c \
    libdeep/dpevaluation.c \
    libdeep/dpindivid.c \
    libdeep/dploop.c \
    libdeep/dpmutation.c \
    libdeep/dpopt.c \
    libdeep/dposda.c \
    libdeep/dppopulation.c \
    libdeep/dprecombination.c \
    libdeep/dpsettings.c \
    libdeep/dptarget.c \
    libxmod/xmmodel.c \
    tools/deepmethod.c \
    tools/xmtranslate.c

DISTFILES += \
    .gitignore \
    AUTHORS \
    CMakeLists.txt \
    COPYING \
    ChangeLog \
    Makefile.am \
    NEWS \
    Nir.pro.user \
    README \
    build_mpich \
    configure.ac \
    deep.anjuta \
    deep.cbp \
    deep.layout \
    deepmethod.geany \
    doc/deepmethod.1 \
    doc/deepmethod.conf.5 \
    libdeep/libdeep-1.0.pc.in \
    libdeep/libdeep_openmpi-1.0.pc.in \
    libxmod/libxmod-1.0.pc.in \
    nsis/DEEPMETHOD.nsis.in \
    nsis/EnvVarUpdate.nsh \
    nsis/exetype.pl \
    nsis/mingwbuild.in \
    po/ChangeLog \
    po/LINGUAS \
    po/POTFILES.in \
    po/ru.po

HEADERS += \
    config.h \
    libdeep/dparchive.h \
    libdeep/dpdeep.h \
    libdeep/dpevaluation.h \
    libdeep/dpindivid.h \
    libdeep/dploop.h \
    libdeep/dpmutation.h \
    libdeep/dpopt.h \
    libdeep/dposda.h \
    libdeep/dppopulation.h \
    libdeep/dprecombination.h \
    libdeep/dpsettings.h \
    libdeep/dptarget.h \
    libxmod/xmmodel.h \
    tools/deepmethod.h \
    tools/xmtranslate.h
