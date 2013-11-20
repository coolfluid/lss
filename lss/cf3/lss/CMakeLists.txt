
coolfluid3_add_library(
  TARGET coolfluid_lss
  LIBS ${LAPACK_LIBRARIES}

  SOURCES
    LibLSS.cpp
    LibLSS.hpp
    linearsystem.cpp
    linearsystem.hpp
    matrix.hpp
    utilities.cpp
    utilities.hpp

    GaussianElimination.cpp
    GaussianElimination.hpp
    GMRES.cpp
    GMRES.hpp
    LAPACK.cpp
    LAPACK.hpp
  )

coolfluid_mark_not_orphan(
  reference/a.cpp
  reference/b.cpp
  reference/c.cpp
  reference/index.hpp
  reference/linearsystem.hpp
  reference/lss_index.hpp
  reference/lss_matrix.hpp
  reference/lss_utilities.hpp
  reference/mmatrix.h
  )
