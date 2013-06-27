
coolfluid3_add_library(
  TARGET coolfluid_lss
  SOURCES
    LibLSS.cpp
    LibLSS.hpp

    linearsystem.hpp
    GaussianElimination.h
    GaussianElimination.cpp )

coolfluid_mark_not_orphan(
  mlinearsystem.cpp  # for reference, will be gone with the wind
  mlinearsystem.h    #  ...
  mmatrix.h          #  ...
  pardiso.cpp
  pardiso.h )

