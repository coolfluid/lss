
coolfluid3_add_library(
  TARGET coolfluid_lss
  SOURCES
    LibLSS.cpp
    LibLSS.hpp

    linearsystem.hpp
    lss_index.hpp
    lss_linearsystem.hpp
    lss_matrix.hpp
    lss_utilities.cpp
    lss_utilities.hpp
  )

coolfluid_mark_not_orphan(

    GMRES.cpp
    GMRES.hpp
    GaussianElimination.cpp
    GaussianElimination.hpp
    LAPACK.cpp
    LAPACK.hpp
  reference/linearsystem.hpp
  reference/lss_index.cpp
  reference/lss_index.hpp
  reference/lss_matrix.hpp
  reference/mlinearsystem.cpp
  reference/mlinearsystem.h
  reference/mmatrix.h
  )

