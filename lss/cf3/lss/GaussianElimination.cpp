// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "GaussianElimination.hpp"


namespace cf3 {
namespace lss {


template<> std::string GaussianElimination< double  >::type_name() { return "GaussianElimination_LongPrecisionReal";  }
template<> std::string GaussianElimination< zdouble >::type_name() { return "GaussianElimination_LongPrecisionComplex";  }
template<> std::string GaussianElimination< float   >::type_name() { return "GaussianElimination_ShortPrecisionReal"; }
template<> std::string GaussianElimination< zfloat  >::type_name() { return "GaussianElimination_ShortPrecisionComplex"; }
common::ComponentBuilder< GaussianElimination< double  >, common::Component, LibLSS > Builder_GaussianElimination_LongPrecisionReal;
common::ComponentBuilder< GaussianElimination< zdouble >, common::Component, LibLSS > Builder_GaussianElimination_LongPrecisionComplex;
common::ComponentBuilder< GaussianElimination< float   >, common::Component, LibLSS > Builder_GaussianElimination_ShortPrecisionReal;
common::ComponentBuilder< GaussianElimination< zfloat  >, common::Component, LibLSS > Builder_GaussianElimination_ShortPrecisionComplex;


}  // namespace lss
}  // namespace cf3

