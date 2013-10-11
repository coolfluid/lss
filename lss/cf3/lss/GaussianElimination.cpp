// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "GaussianElimination.hpp"


namespace cf3 {
namespace lss {


template<> std::string GaussianElimination< double >::type_name() { return "GaussianElimination_DoubleP"; }
template<> std::string GaussianElimination< float  >::type_name() { return "GaussianElimination_SingleP"; }
common::ComponentBuilder< GaussianElimination< double >, common::Component, LibLSS > Builder_GaussianElimination_DoubleP;
common::ComponentBuilder< GaussianElimination< float  >, common::Component, LibLSS > Builder_GaussianElimination_SingleP;


}  // namespace lss
}  // namespace cf3

