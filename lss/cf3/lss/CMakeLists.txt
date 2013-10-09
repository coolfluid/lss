
coolfluid3_add_library(
  TARGET coolfluid_lss
  SOURCES
    LibLSS.cpp
    LibLSS.hpp

    linearsystem.hpp
    lss_index.hpp
    lss_matrix.hpp
    lss_utilities.cpp
    lss_utilities.hpp

    LAPACK.cpp
    LAPACK.hpp
  )

coolfluid_mark_not_orphan(
    GaussianElimination.cpp
    GaussianElimination.hpp
  reference/GMRES.cpp
  reference/GMRES.h
  reference/linearsystem.hpp
  reference/lss_index.cpp
  reference/lss_index.hpp
  reference/lss_matrix.hpp
  reference/mlinearsystem.cpp
  reference/mlinearsystem.h
  reference/mmatrix.h
  )

