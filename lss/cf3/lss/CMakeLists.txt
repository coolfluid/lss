
coolfluid3_add_library(
  TARGET coolfluid_lss
  LIBS ${LAPACK_LIBRARIES}

  SOURCES
    LibLSS.cpp
    LibLSS.hpp
    detail/index.hpp
    detail/matrix.hpp
    detail/utilities.cpp
    detail/utilities.hpp
    linearsystem.hpp

    GMRES.cpp
    GMRES.hpp
    GaussianElimination.cpp
    GaussianElimination.hpp
    LAPACK.cpp
    LAPACK.hpp
  )

coolfluid_mark_not_orphan( reference/mmatrix.h )
