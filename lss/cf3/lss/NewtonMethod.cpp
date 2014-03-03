// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "NewtonMethod.hpp"


namespace cf3 {
namespace lss {


template<> std::string NewtonMethod< double >::type_name() { return "NewtonMethod";  }
template<> std::string NewtonMethod< float  >::type_name() { return "NewtonMethod_ShortPrecisionReal"; }
common::ComponentBuilder< NewtonMethod< double >, common::Component, LibLSS > Builder_NewtonMethod;
common::ComponentBuilder< NewtonMethod< float  >, common::Component, LibLSS > Builder_NewtonMethod_ShortPrecisionReal;


}  // namespace lss
}  // namespace cf3

