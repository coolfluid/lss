
coolfluid3_add_library(
  TARGET coolfluid_lss
  LIBS ${LAPACK_LIBRARIES}

  SOURCES
    LibLSS.cpp
    LibLSS.hpp
    detail/index.hpp
    detail/matrix.hpp
    detail/matrix.cpp
    detail/utilities.cpp
    detail/utilities.hpp
    linearsystem.hpp

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
  reference/linearsystem.hpp
  reference/lss_index.hpp
  reference/lss_matrix.hpp
  reference/lss_utilities.hpp
  reference/mmatrix.h
  )
