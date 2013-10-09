// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "LAPACK.hpp"


namespace cf3 {
namespace lss {


template<> std::string LAPACK< double >::type_name() { return "LAPACK_DoubleP"; }
template<> std::string LAPACK< float  >::type_name() { return "LAPACK_SingleP"; }
common::ComponentBuilder< LAPACK< double >, common::Component, LibLSS > Builder_LAPACK_DoubleP;
common::ComponentBuilder< LAPACK< float  >, common::Component, LibLSS > Builder_LAPACK_SingleP;


}  // namespace lss
}  // namespace cf3

