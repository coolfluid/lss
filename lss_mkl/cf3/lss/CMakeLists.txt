
list( APPEND mkl_files
  LibLSS_MKL.cpp
  LibLSS_MKL.hpp
  dss.cpp
  dss.hpp
  iss.cpp
  iss.hpp
  pardiso.cpp
  pardiso.hpp )


if(CF3_HAVE_INTELMKL)
  coolfluid3_add_library(
    TARGET    coolfluid_lss_mkl
    LIBS      coolfluid_lss ${INTELMKL_LIBRARIES}
    INCLUDES  ${INTELMKL_INCLUDE_DIRS}
    SOURCES   ${mkl_files} )
else()
  coolfluid_mark_not_orphan( ${mkl_files} )
endif()

