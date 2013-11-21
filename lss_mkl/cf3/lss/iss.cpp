// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


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
  environment_variable_t< int > nthreads("OMP_NUM_THREADS",1);
  CFinfo << "mkl iss: OMP_NUM_THREADS: " << nthreads.description() << CFendl;
  mkl_set_num_threads(nthreads.value);

  linearsystem< double >::initialize(_size_i,_size_j,_size_k,_value);
}


iss& iss::solve()
{
//nrhs = static_cast< int >(m_b.size(1));
  return *this;
}


iss& iss::copy(const iss& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  return *this;
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3

