// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
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
  void wgsmatvec_(int*, int*, int*, double*, double*, double*, int*, int*);
}


namespace cf3 {
namespace lss {


common::ComponentBuilder< WSMP, common::Component, LibLSS_WSMP > Builder_WSMP;


WSMP::WSMP(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k)
  : linearsystem< double >(name)
{
  environment_variable_t< int >
    nthreads("WSMP_NUM_THREADS",1),
    malloc_trh("MALLOC_TRIM_THRESHOLD_"),
    malloc_max("MALLOC_MMAP_MAX_");
  environment_variable_t< std::string >
    licpath    ("WSMPLICPATH"),
    wincoremem ("WINCOREMEM"),
    woocdir0   ("WOOCDIR0");

  CFinfo  << "WSMP: WSMP_NUM_THREADS:       " << nthreads  .description() << CFendl
          << "WSMP: WSMPLICPATH:            " << licpath   .description() << CFendl
          << "WSMP: MALLOC_TRIM_THRESHOLD_: " << malloc_trh.description() << " (optimal -1: disable trim, suggest 524288000 [bytes])" /* 500Mb */<< CFendl
          << "WSMP: MALLOC_MMAP_MAX_:       " << malloc_max.description() << " (optimal 0)" << CFendl;
  CFdebug << "WSMP: WINCOREMEM:             " << wincoremem.description() << CFendl
          << "WSMP: WOOCDIR0:               " << woocdir0  .description() << CFendl;

  wsetmaxthrds_(&nthreads.value);

  for (size_t i=0; i<64; ++i) iparm[i] = 0;
  for (size_t i=0; i<64; ++i) dparm[i] = 0.;

  // reset pt, iparm and dparm defaults
  call_wsmp(0);
  iparm[ 3] = 0;  // CSR matrix format
  iparm[ 4] = 0;  // + C-style numbering
  iparm[19] = 2;  // + ordering option 5

  linearsystem< double >::initialize(_size_i,_size_j,_size_k);
}


WSMP& WSMP::solve()
{
  int err;
  if ( (err=call_wsmp(1)) ||  // 1: analysis and reordering
       (err=call_wsmp(2)) ||  // 2: LU factorization
       (err=call_wsmp(3)) ||  // 3: forward and backward elimination
       (err=call_wsmp(4)) )   // 4: iterative refinement
    throw std::runtime_error(err_message(err));

  /*
   * iparm[ 2]: task where last error occured
   * iparm[23]: task 1 number of nonzeros in LU factors
   * dparm[23]: task 1 number of FLOPS in factorization
   * dparm[6]:  task 4 maximum relative error
   * iparm[25]: task ? summary # iterations
   * dparm[25]: task ? summary residual
   */

  m_b.swap(m_x);
  return *this;
}


WSMP& WSMP::multi(const double& _alpha, const double& _beta)
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  vector_t b = m_b;

  int err = 0;
  for (int k=0, fmt(iparm[3]+1); k<static_cast< int >(size(2)) && !err;  ++k)
    wgsmatvec_(
      &A.nnu, &A.ia[0], &A.ja[0], &A.a[0],
      &m_x.a[A.nnu*k], &b.a[A.nnu*k], &fmt, &err );
  if (err)
    throw std::runtime_error(err_message(err));

  for (size_t i=0; i<m_b.a.size(); ++i)
    m_b.a[i] = _alpha*b.a[i] + _beta*m_b.a[i];
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


WSMP& WSMP::swap(WSMP& _other)
{
  linearsystem< double >::swap(_other);
  m_A.swap(_other.m_A);
  return *this;
}


const std::string WSMP::err_message(const int& err)
{
  std::ostringstream s;
  s << "WSMP error: " << err << ": ";
  err>    0? s << "matrix close enough to singular, suspected pivot at A(" << (err-1) << ',' << (err-1) << ") (0-based)" :
  err==   0? s << "(success)" :
  err==-102? s << "failed dynamic memory allocation" :
  err==-103? s << "probable integer overflow in large matrix" :
  err==-300? s << "invalid operation, maybe a previous task did not finish successfuly" :
  err==-700? s << "internal error, maybe check the input matrix structure" :
  err==-900? s << "license is expired, invalid, or missing" :
  err==-501? s << "environment variable not set (WINCOREMEM)" :
  err==-502? s << "environment variable not set (WOOCDIR0)" :
  err==-503? s << "cannot write to storage" :
  err==-504? s << "value of WINCOREMEM insufficient" :
             s << "(unknown error)";
  s << '.';
  return s.str();
}


int WSMP::call_wsmp(int _phase)
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  int nrhs = static_cast< int >(m_b.size(1)),
      ldb  = static_cast< int >(m_b.size(0)),
     &fact = iparm[30],
      ldlt_pivot(fact==2 || fact==4 || fact==6 || fact==7);

  iparm[1] = iparm[2] = _phase;
  wgsmp_(
    &A.nnu,&A.ia[0],&A.ja[0],&A.a[0],
    &m_b.a[0],&ldb,&nrhs,NULL,iparm,dparm);

  iparm[63] = (iparm[63]>0 && ldlt_pivot? 0 : iparm[63]);
  return iparm[63];
}


}  // namespace lss
}  // namespace cf3

