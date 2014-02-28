// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "QuasiNewtonMethod.hpp"


namespace cf3 {
namespace lss {


template<> std::string QuasiNewtonMethod< double >::type_name() { return "QuasiNewtonMethod";  }
template<> std::string QuasiNewtonMethod< float  >::type_name() { return "QuasiNewtonMethod_ShortPrecisionReal"; }
common::ComponentBuilder< QuasiNewtonMethod< double >, common::Component, LibLSS > Builder_QuasiNewtonMethod;
common::ComponentBuilder< QuasiNewtonMethod< float  >, common::Component, LibLSS > Builder_QuasiNewtonMethod_ShortPrecisionReal;


}  // namespace lss
}  // namespace cf3

