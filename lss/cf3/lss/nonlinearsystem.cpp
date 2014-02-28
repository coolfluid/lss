// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "common/Builder.hpp"

#include "nonlinearsystem.hpp"


namespace cf3 {
namespace lss {


common::RegisterComponent< nonlinearsystem< double  >, LibLSS > Register_nonlinearsystem_LongPrecisionReal;
common::RegisterComponent< nonlinearsystem< float   >, LibLSS > Register_nonlinearsystem_ShortPrecisionReal;


}  // namespace lss
}  // namespace cf3

