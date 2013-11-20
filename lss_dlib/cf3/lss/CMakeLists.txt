
list( APPEND dlib_files
  LibLSS_DLIB.cpp
  LibLSS_DLIB.hpp
  dlib.cpp
  dlib.hpp )


if(CF3_HAVE_DLIB)
  coolfluid3_add_library(
    TARGET    coolfluid_lss_dlib
    LIBS      coolfluid_lss ${PARDISO_LIBRARIES}
    SOURCES   ${dlib_files} )
  # LINKFLAGS -fopenmp
else()
  coolfluid_mark_not_orphan( ${dlib_files} )
endif()

