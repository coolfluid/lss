
list( APPEND pardiso_files
  LibLSS_PARDISO.cpp
  LibLSS_PARDISO.hpp
  pardiso.cpp
  pardiso.hpp )


if(CF3_HAVE_PARDISO)
  coolfluid3_add_library(
    TARGET    coolfluid_lss_pardiso
    LIBS      coolfluid_lss ${PARDISO_LIBRARIES}
    SOURCES   ${pardiso_files} )
  # LINKFLAGS -fopenmp
else()
  coolfluid_mark_not_orphan( ${pardiso_files} )
endif()

