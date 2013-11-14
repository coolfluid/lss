// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf

#include "mkl_rci.h" 
#include "mkl_service.h"

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "iss.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< iss, common::Component, LibLSS_MKL > Builder_MKL_iss;


iss::iss(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k, const double& _value)
  : linearsystem< double >(name)
{
  const char *nthreads = getenv("OMP_NUM_THREADS");
  int nthd = 1;
  sscanf(nthreads? nthreads:"1","%d",&nthd);
  mkl_set_num_threads(nthd);

  CFinfo  << "mkl iss: OMP_NUM_THREADS: " << nthd << (nthreads? " (set)":" (not set)") << CFendl;

  linearsystem< double >::initialize(_size_i,_size_j,_size_k,_value);
}


iss& iss::solve()
{
//nrhs = static_cast< int >(m_b.size(1));
  return *this;
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3

