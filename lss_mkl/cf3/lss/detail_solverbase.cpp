// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "mkl_service.h"
#include "mkl_spblas.h"

#include "common/Log.hpp"
#include "detail_solverbase.h"


namespace cf3 {
namespace lss {
namespace mkl {
namespace detail {


solverbase::solverbase(const std::string& name)
  : linearsystem< double >(name)
{
  environment_variable_t< int > nthreads("OMP_NUM_THREADS",1);
  CFinfo << "mkl: OMP_NUM_THREADS: " << nthreads.description() << CFendl;
  mkl_set_num_threads(nthreads.value);
}


solverbase& solverbase::multi(const double& _alpha, const double& _beta)
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  char
    transa = 'N',                           // not transposed,
    matdescra[6] = "G--F-";                 // general, 1-based,
  int                                       // ...
    m = static_cast< int >(this->size(0)),  // ...
    n = static_cast< int >(this->size(2)),  // ...
    k = m;                                  // square matrix (phew!)

  mkl_dcsrmm( &transa, &m, &n, &k, const_cast< double* >(&_alpha),
    &matdescra[0], &A.a[0], &A.ja[0], &A.ia[0], &A.ia[1],
    const_cast< double* >(&m_x.a[0]), &m,
    const_cast< double* >(&_beta), &m_b.a[0], &m );

  return *this;
}


solverbase& solverbase::swap(solverbase& _other)
{
  linearsystem< double >::swap(_other);
  m_A.swap(_other.m_A);
  return *this;
}


}  // namespace detail
}  // namespace mkl
}  // namespace lss
}  // namespace cf3

