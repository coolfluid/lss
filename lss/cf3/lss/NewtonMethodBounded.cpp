// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "NewtonMethodBounded.hpp"


namespace cf3 {
namespace lss {


template<> std::string NewtonMethodBounded< double >::type_name() { return "NewtonMethodBounded";  }
template<> std::string NewtonMethodBounded< float  >::type_name() { return "NewtonMethodBounded_ShortPrecisionReal"; }
common::ComponentBuilder< NewtonMethodBounded< double >, common::Component, LibLSS > Builder_NewtonMethodBounded;
common::ComponentBuilder< NewtonMethodBounded< float  >, common::Component, LibLSS > Builder_NewtonMethodBounded_ShortPrecisionReal;


}  // namespace lss
}  // namespace cf3

