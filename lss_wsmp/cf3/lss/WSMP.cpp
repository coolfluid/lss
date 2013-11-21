// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf


#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "WSMP.hpp"


// Prototypes for external function calls
extern "C" {
  void wsetmaxthrds_(int*);
  void wgsmp_(int*, int*, int*, double*, double*, int*, int*, double*, int*, double*);
}


namespace cf3 {
namespace lss {


common::ComponentBuilder< WSMP, common::Component, LibLSS_WSMP > Builder_WSMP;


WSMP::WSMP(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k, const double& _value)
  : linearsystem< double >(name)
{
  environment_variable_t< int > nthreads("WSMP_NUM_THREADS",1);
  environment_variable_t< std::string >
    licpath    ("WSMPLICPATH"),
    wincoremem ("WINCOREMEM"),
    woocdir0   ("WOOCDIR0"),
    malloc_trh ("MALLOC_TRIM_THRESHOLD_"),
    malloc_max ("MALLOC_MMAP_MAX_");

  CFinfo  << "WSMP: WSMP_NUM_THREADS:       " << nthreads  .description() << CFendl
          << "WSMP: WSMPLICPATH:            " << licpath   .description() << CFendl;
  CFdebug << "WSMP: WINCOREMEM:             " << wincoremem.description() << CFendl
          << "WSMP: WOOCDIR0:               " << woocdir0  .description() << CFendl
          << "WSMP: MALLOC_TRIM_THRESHOLD_: " << malloc_trh.description() << CFendl
          << "WSMP: MALLOC_MMAP_MAX_:       " << malloc_max.description() << CFendl;

  // iparm/dparm defaults
  std::fill_n(&iparm[0],64,0);
  std::fill_n(&dparm[0],64,0.);
  call_wsmp(0);

  wsetmaxthrds_(&nthreads.value);
  iparm[ 4] = 0;  // + C-style numbering
  iparm[19] = 2;  // + ordering option 5

  linearsystem< double >::initialize(_size_i,_size_j,_size_k,_value);
}


WSMP& WSMP::solve()
{
  if (call_wsmp(1) ||  // analysis and reordering
      call_wsmp(2) ||  // LU factorization
      call_wsmp(3) ||  // forward and backward elimination
      call_wsmp(4))    // iterative refinement
  {
    const int &err = iparm[63];
    std::ostringstream msg;
    msg << "WSMP: task " << iparm[2] << " error " << err << ": ";

    err>    0? msg << "matrix close enough to singular, suspected pivot at i=j=" << (err-1) << " (0-based)." :
    err< -100? msg << "error in input argument iparm[" << (-err-1) << "]." :
    err==-102? msg << "failed dynamic memory allocation." :
    err==-103? msg << "probable integer overflow in large matrix." :
    err==-300? msg << "invalid operation, maybe a previous task did not finish successfuly." :
    err==-700? msg << "internal error, maybe check the input matrix structure." :
    err==-900? msg << "license is expired, invalid, or missing." :
    (iparm[35] && err==-501)? msg << "environment variable not set (WINCOREMEM)." :
    (iparm[35] && err==-502)? msg << "environment variable not set (WOOCDIR0)." :
    (iparm[35] && err==-503)? msg << "cannot write to storage (full?)." :
    (iparm[35] && err==-504)? msg << "value of WINCOREMEM insufficient." :
               msg << "unknown error.";

    throw std::runtime_error(msg.str());
  }
  m_b.swap(m_x);

  /*
   * iparm[23]: task 1 number of nonzeros in LU factors
   * dparm[23]: task 1 number of FLOPS in factorization
   * dparm[6]:  task 4 maximum relative error
   * iparm[25]: task ? summary # iterations
   * dparm[25]: task ? summary residual
   */

  return *this;
}


WSMP& WSMP::copy(const WSMP& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  for (size_t i=0; i<64; ++i) {
    dparm[i] = _other.dparm[i];
    iparm[i] = _other.iparm[i];
  }
  return *this;
}


int WSMP::call_wsmp(int _task)
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  int nrhs = static_cast< int >(m_b.size(1)),
      ldb  = static_cast< int >(m_b.size(0)),
     &fact = iparm[30],
      ldlt_pivot(fact==2 || fact==4 || fact==6 || fact==7);

  iparm[1] = iparm[2] = _task;
  wgsmp_(
    &A.nnu,&A.ia[0],&A.ja[0],&A.a[0],
    &m_b.a[0],&ldb,&nrhs,NULL,iparm,dparm);

  iparm[63] = (iparm[63]>0 && ldlt_pivot? 0 : iparm[63]);
  return iparm[63];
}


}  // namespace lss
}  // namespace cf3

