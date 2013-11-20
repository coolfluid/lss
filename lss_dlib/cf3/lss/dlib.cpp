// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


//include <cmath>
//include <algorithm>
#include <cstdio>  // for sscanf

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "dlib.hpp"


namespace cf3 {
namespace lss {


common::ComponentBuilder< dlib, common::Component, LibLSS_PARDISO > Builder_pardiso;


dlib& dlib::solve()
{
  int nrhs = static_cast< int >(m_b.size(1));

  dlib::matrix< double, 0, 1 > dlib_x = dlib::inv(m_A.a) * dlib::vector_to_matrix(m_b.a);
  for (size_t i=0; i<size(0); ++i)
    x(i) = dlib_x(i);

  return *this;
}


}  // namespace lss
}  // namespace cf3

