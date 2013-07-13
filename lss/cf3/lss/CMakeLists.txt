
coolfluid3_add_library(
  TARGET coolfluid_lss
  SOURCES
    LibLSS.cpp
    LibLSS.hpp

    linearsystem.hpp
    lss_index.cpp
    lss_index.hpp
    lss_matrix.hpp
    dummy_linearsystem.h
    dummy_linearsystem.cpp
 )

coolfluid_mark_not_orphan(
  mlinearsystem.cpp  # for reference, will be gone with the wind
  mlinearsystem.h    #  ...
  mmatrix.h          #  ...
  GaussianElimination.cpp
  GaussianElimination.h
  GMRES.cpp
  GMRES.h
  pardiso.cpp
  pardiso.h )

