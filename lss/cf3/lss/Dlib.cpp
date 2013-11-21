// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"
#include "Dlib.hpp"


namespace cf3 {
namespace lss {


template<> std::string Dlib< double >::type_name() { return "Dlib"; }
template<> std::string Dlib< float  >::type_name() { return "Dlib_SinglePrecision"; }
common::ComponentBuilder< Dlib< double >, common::Component, LibLSS > Builder_Dlib;
common::ComponentBuilder< Dlib< float  >, common::Component, LibLSS > Builder_Dlib_SinglePrecision;


}  // namespace lss
}  // namespace cf3

