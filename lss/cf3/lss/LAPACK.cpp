// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "LAPACK.hpp"


namespace cf3 {
namespace lss {


template<> std::string LAPACK< double  >::type_name() { return "LAPACK"; }
template<> std::string LAPACK< zdouble >::type_name() { return "LAPACK_LongPrecisionComplex"; }
template<> std::string LAPACK< float   >::type_name() { return "LAPACK_ShortPrecisionReal"; }
template<> std::string LAPACK< zfloat  >::type_name() { return "LAPACK_ShortPrecisionComplex"; }
common::ComponentBuilder< LAPACK< double  >, common::Component, LibLSS > Builder_LAPACK_LongPrecisionReal;
common::ComponentBuilder< LAPACK< zdouble >, common::Component, LibLSS > Builder_LAPACK_LongPrecisionComplex;
common::ComponentBuilder< LAPACK< float   >, common::Component, LibLSS > Builder_LAPACK_ShortPrecisionReal;
common::ComponentBuilder< LAPACK< zfloat  >, common::Component, LibLSS > Builder_LAPACK_ShortPrecisionComplex;


}  // namespace lss
}  // namespace cf3

