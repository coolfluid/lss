// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "mkl_service.h"
#include "mkl_spblas.h"

#include "common/Log.hpp"
#include "detail_mkl_solver_base.h"


namespace cf3 {
namespace lss {
namespace mkl {
namespace detail {


mkl_solver_base::mkl_solver_base(const std::string& name)
  : linearsystem< double >(name)
{
  environment_variable_t< int > nthreads("OMP_NUM_THREADS",1);
  CFinfo << "mkl pardiso: OMP_NUM_THREADS: " << nthreads.description() << CFendl;
  mkl_set_num_threads(nthreads.value);
}


void mkl_solver_base::A___multi(const vector_t& _x, vector_t& _b)
{
  // sparse matrix - dense matrix multiplication:
  // b(m,n) = alpha A(m,k) x(k,n) + beta b(m,n)

  matrix_t::matrix_compressed_t& A = m_A.compress();
  double
    alpha = 1.,
    beta  = 0.;
  char
    transa = 'N',                           // not transposed,
    matdescra[6] = "G--F-";                 // general, 1-based,
  int                                       // ...
    m = static_cast< int >(this->size(0)),  // ...
    n = static_cast< int >(this->size(2)),  // ...
    k = m;                                  // square matrix (phew!)

  mkl_dcsrmm(
    &transa, &m, &n, &k, &alpha,
    &matdescra[0], &A.a[0], &A.ja[0], &A.ia[0], &A.ia[1],
    const_cast< double* >(&_x.a[0]), &n, &beta, &_b.a[0], &n );
}


}  // namespace detail
}  // namespace mkl
}  // namespace lss
}  // namespace cf3

