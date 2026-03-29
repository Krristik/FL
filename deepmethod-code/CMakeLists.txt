cmake_minimum_required(VERSION 3.15)
project(deepmethod C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

## Required for FindIntltool module
#set(GETTEXT_PACKAGE ${PROJECT_NAME})
#set(GETTEXT_PO_DIR ${CMAKE_SOURCE_DIR}/po)

find_package(PkgConfig)
pkg_check_modules(GLIB2 REQUIRED glib-2.0>=2.36.0)
include_directories(${GLIB2_INCLUDE_DIRS})

add_compile_definitions(DEEPMETHOD_LOCALEDIR=\"${CMAKE_INSTALL_PREFIX}/share/locale\")

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

#add the Deep library
add_subdirectory(libdeep)
#add the XMod library
add_subdirectory(libxmod)

add_executable(deepmethod
	config.h
        tools/deepmethod.c
        tools/xmtranslate.c
)

target_link_libraries(deepmethod  PUBLIC 
        ${GLIB2_LIBRARIES}
        m
        gio-2.0
        gobject-2.0
        deep
        xmod
)

target_include_directories(deepmethod PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tools"
)

link_directories(${GLIB2_LIBRARY_DIRS})

install(TARGETS deepmethod DESTINATION bin)
install(FILES
        "config.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/tools/deepmethod.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/tools/xmtranslate.h"
        DESTINATION include
)