// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"

#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


common::RegisterComponent< linearsystem< double >, LibLSS > Register_linearsystem;
common::RegisterComponent< linearsystem< float  >, LibLSS > Register_linearsystem_SinglePrecision;


}  // namespace lss
}  // namespace cf3

