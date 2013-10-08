// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "GaussianElimination.hpp"
#include "LAPACK.hpp"


namespace cf3 {
namespace lss {


#if 0
/* double precision, Gaussian elimination solvers */
common::ComponentBuilder< GaussianElimination, common::Component, LibLSS > Builder_GaussianEliminationA_P2;


/* double precision, Gaussian elimination solvers */
template<> std::string GaussianElimination< double, dense_matrix_v<  double > >::type_name() { return "GaussianEliminationV";  }
template<> std::string GaussianElimination< double, dense_matrix_vv< double > >::type_name() { return "GaussianEliminationVV"; }

common::ComponentBuilder< GaussianElimination< double, dense_matrix_v<  double > >, common::Component, LibLSS > Builder_GaussianEliminationV_P2;
common::ComponentBuilder< GaussianElimination< double, dense_matrix_vv< double > >, common::Component, LibLSS > Builder_GaussianEliminationVV_P2;


/* single precision, Gaussian elimination solvers */
// FIXME: these don't work? there seems to be a segmentation fault from python?
template<> std::string GaussianElimination< float, dense_matrix_v<  float > >::type_name() { return "GaussianEliminationV_P1";  }
template<> std::string GaussianElimination< float, dense_matrix_vv< float > >::type_name() { return "GaussianEliminationVV_P1"; }

common::ComponentBuilder< GaussianElimination< float, dense_matrix_v<  float > >, common::Component, LibLSS > Builder_GaussianEliminationV_P1;
common::ComponentBuilder< GaussianElimination< float, dense_matrix_vv< float > >, common::Component, LibLSS > Builder_GaussianEliminationVV_P1;
#endif


}  // namespace lss
}  // namespace cf3

