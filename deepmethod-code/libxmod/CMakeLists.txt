add_library(xmod xmmodel.c)

target_include_directories(xmod
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
)


install(TARGETS xmod
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
)
install(FILES xmmodel.h DESTINATION include)