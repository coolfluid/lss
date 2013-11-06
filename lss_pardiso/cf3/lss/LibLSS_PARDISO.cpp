// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "cf3/common/RegistLibrary.hpp"
#include "cf3/lss/LibLSS_PARDISO.hpp"


namespace cf3 {
namespace lss {


cf3::common::RegistLibrary< LibLSS_PARDISO > LibLSS_PARDISO;


void LibLSS_PARDISO::initiate()
{
  if (m_is_initiated)
    return;

  Handle< common::Component > lib = common::Core::instance().libraries().get_child_checked("cf3.lss.pardiso");

  m_is_initiated = true;
}


}  // lss
}  // cf3

