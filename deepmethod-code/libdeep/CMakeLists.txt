add_library(deep 
    dpdeep.c
    dpevaluation.c
    dpindivid.c
    dploop.c
    dpmutation.c
    dpopt.c
    dposda.c
    dppopulation.c
    dprecombination.c
    dpsettings.c
    dptarget.c
	dparchive.c
)

target_include_directories(deep
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(deep  PRIVATE
        ${GLIB2_LIBRARIES}
        m
        gio-2.0
        gobject-2.0
)

install(TARGETS deep
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
)
install(FILES
    dpdeep.h
    dpevaluation.h
    dpindivid.h
    dploop.h
    dpmutation.h
    dpopt.h
    dposda.h
    dppopulation.h
    dprecombination.h
    dpsettings.h
    dptarget.h
	dparchive.h
    DESTINATION include
)