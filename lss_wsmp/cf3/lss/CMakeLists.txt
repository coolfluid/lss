
list( APPEND wsmp_files
  LibLSS_WSMP.cpp
  LibLSS_WSMP.hpp
  WSMP.cpp
  WSMP.hpp )


if(CF3_HAVE_WSMP)
  coolfluid3_add_library(
    TARGET  coolfluid_lss_wsmp
    LIBS    coolfluid_lss ${WSMP_LIBRARIES}
    SOURCES ${wsmp_files} )
else()
  coolfluid_mark_not_orphan( ${wsmp_files} )
endif()

