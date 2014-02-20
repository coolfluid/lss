// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "cf3/common/RegistLibrary.hpp"
#include "cf3/lss/LibLSS.hpp"


namespace cf3 {
namespace lss {


cf3::common::RegistLibrary< LibLSS > LibLSS;


void LibLSS::initiate()
{
  if (m_is_initiated)
    return;

  Handle< common::Component > lib = common::Core::instance().libraries().get_child("cf3.lss");
  cf3_assert(lib);

  m_is_initiated = true;
}


}  // lss
}  // cf3

