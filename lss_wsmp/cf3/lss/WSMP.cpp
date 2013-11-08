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


WSMP::WSMP(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k, const double& _value) : linearsystem_t(name) {

  char
     *nthreads   = getenv("WSMP_NUM_THREADS"),
     *licpath    = getenv("WSMPLICPATH"),
     *wincoremem = getenv("WINCOREMEM"),
     *woocdir0   = getenv("WOOCDIR0"),
     *malloc_trh = getenv("MALLOC_TRIM_THRESHOLD_"),
     *malloc_max = getenv("MALLOC_MMAP_MAX_");

  int nthd = 1;
  sscanf(nthreads? nthreads:"1","%d",&nthd);
  wsetmaxthrds_(&nthd);

  CFinfo  << "WSMP_NUM_THREADS:       " << nthd << " (" << (nthreads? "set)":"not set)") << CFendl
          << "WSMPLICPATH:            " << (licpath? licpath:"(WSMPLICPATH: not set)")   << CFendl
          << "MALLOC_TRIM_THRESHOLD_: " << (malloc_trh? "(set)":"(not set)") << CFendl
          << "MALLOC_MMAP_MAX_:       " << (malloc_max? "(set)":"(not set)") << CFendl;
  CFdebug << "WINCOREMEM:             " << (wincoremem? wincoremem : "(not set)") << CFendl
          << "WOOCDIR0:               " << (woocdir0?   woocdir0   : "(not set)") << CFendl;

  // iparm/dparm defaults
  for (int i=0; i<64; ++i)  iparm[i] = 0;
  for (int i=0; i<64; ++i)  dparm[i] = 0.;
  call_wsmp(0);

  iparm[ 4] = 0;  // + C-style numbering
  iparm[19] = 2;  // + ordering option 5

  linearsystem_t::initialize(_size_i,_size_j,_size_k,_value);
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
    msg << "WSMP: task " << task << " error " << err << ": ";

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
  b().swap(x());

  /*
   * iparm[23]: task 1 number of nonzeros in LU factors
   * dparm[23]: task 1 number of FLOPS in factorization
   * dparm[6]:  task 4 maximum relative error
   * iparm[25]: task ? summary # iterations
   * dparm[25]: task ? summary residual
   */

  return *this;
}


int WSMP::call_wsmp(int _task)
{
  int nrhs = b().size(1),
      ldb  = b().size(0),
     &fact = iparm[30],
      ldlt_pivot(fact==2 || fact==4 || fact==6 || fact==7);

  iparm[1] = iparm[2] = task = _task;
  wgsmp_(
    &m_A.idx.nnu,&m_A.idx.ia[0],&m_A.idx.ja[0],&m_A.a[0],
    &m_b.a[0],&ldb,&nrhs,NULL,iparm,dparm);

  iparm[63] = (iparm[63]>0 && ldlt_pivot? 0 : iparm[63]);
  return iparm[63];
}


}  // namespace lss
}  // namespace cf3

