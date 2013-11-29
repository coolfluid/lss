
list( APPEND petsc_files
  LibLSS_PETSC.cpp
  LibLSS_PETSC.hpp
  petsc_seq.cpp
  petsc_seq.h )


if(CF3_HAVE_PETSC)
  coolfluid3_add_library(
    TARGET    coolfluid_lss_petsc
    SOURCES   ${petsc_files}
    LIBS      coolfluid_lss ${PETSC_LIBRARIES}
    INCLUDES  ${PETSC_INCLUDE_DIRS} )
else()
  coolfluid_mark_not_orphan( ${petsc_files} )
endif()

