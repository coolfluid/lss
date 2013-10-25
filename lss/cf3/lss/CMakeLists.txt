
coolfluid3_add_library(
  TARGET coolfluid_lss
  SOURCES
    LibLSS.cpp
    LibLSS.hpp
    detail/index.hpp
    detail/linearsystem.hpp
    detail/matrix.hpp
    detail/utilities.cpp
    detail/utilities.hpp
    linearsystem.h

    GMRES.cpp
    GMRES.hpp
    GaussianElimination.cpp
    GaussianElimination.hpp
    LAPACK.cpp
    LAPACK.hpp
  )

coolfluid_mark_not_orphan(
  reference/linearsystem.hpp
  reference/lss_matrix.hpp
  reference/mlinearsystem.cpp
  reference/mlinearsystem.h
  reference/mmatrix.h
  )

