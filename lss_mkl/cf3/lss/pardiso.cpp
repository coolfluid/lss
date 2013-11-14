// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf

#include "mkl_pardiso.h"
#include "mkl_service.h"

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "pardiso.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< pardiso, common::Component, LibLSS_MKL > Builder_MKL_pardiso;


pardiso::pardiso(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k, const double& _value)
  : linearsystem< double >(name)
{
  const char
    *nthreads = getenv("OMP_NUM_THREADS"),
    *ooc_path = getenv("MKL_PARDISO_OOC_PATH"),
    *ooc_maxc = getenv("MKL_PARDISO_OOC_MAX_CORE_SIZE"),
    *ooc_swap = getenv("MKL_PARDISO_OOC_MAX_SWAP_SIZE"),
    *ooc_keep = getenv("MKL_PARDISO_OOC_KEEP_FILE");

  int nthd = 1;
  sscanf(nthreads? nthreads:"1","%d",&nthd);
  mkl_set_num_threads(nthd);

  CFinfo  << "mkl pardiso: OMP_NUM_THREADS: " << nthd << (nthreads? " (set)":" (not set)") << CFendl;
  CFdebug << "mkl pardiso: MKL_PARDISO_OOC_PATH:          " << (ooc_path? ooc_path : "(not set)") << CFendl
          << "mkl pardiso: MKL_PARDISO_OOC_MAX_CORE_SIZE: " << (ooc_maxc? ooc_maxc : "(not set)") << CFendl
          << "mkl pardiso: MKL_PARDISO_OOC_MAX_SWAP_SIZE: " << (ooc_swap? ooc_swap : "(not set)") << CFendl
          << "mkl pardiso: MKL_PARDISO_OOC_KEEP_FILE:     " << (ooc_keep? ooc_keep : "(not set)") << CFendl;

  // reset pt and iparm defaults
  for (int i=0; i<64; ++i) iparm[i] = 0;
  perm .assign(m_A.nnu,0);
  PARDISOINIT(pt,&mtype,iparm);

  iparm[ 7] = 0;  // max numbers of iterative refinement steps
//iparm[31] = 0;  // [0|1] sparse direct solver or multi-recursive iterative solver
  maxfct    = 1;  // maximum number of numerical factorizations
  mnum      = 1;  // which factorization to use
  msglvl    = 1;  // message level: output statistical information
  mtype     = 1;  // real structurally symmetric matrix

  linearsystem< double >::initialize(_size_i,_size_j,_size_k,_value);
}


pardiso::~pardiso()
{
  call_pardiso(-1);  // -1: termination and release of memory
}


pardiso& pardiso::solve()
{
  nrhs = static_cast< int >(m_b.size(1));
  if
#if 0
     (call_pardiso(11) ||  // 11: reordering and symbolic factorization
      call_pardiso(22) ||  // 22: numerical factorization and
      call_pardiso(33))    // 33: back substitution and iterative refinement
#else
     (call_pardiso(13))
#endif
  {
    std::ostringstream msg;
    msg << "mkl pardiso: phase " << phase << " error " << err << ": ";
    err==  -1? msg << "input inconsistent." :
    err==  -2? msg << "not enough memory."  :
    err==  -3? msg << "reordering problem." :
    err==  -4? msg << "zero pivot, numerical factorization or iterative refinement problem." :
    err==  -5? msg << "unclassified (internal) error."     :
    err==  -6? msg << "preordering failed (matrix types 11, 13 only)." :
    err==  -7? msg << "diagonal matrix is singular."           :
    err==  -8? msg << "32-bit integer overflow problem."   :
    err==  -9? msg << "not enough memory for OOC." :
    err== -10? msg << "error opening OOC files." :
    err== -11? msg << "read/write error with OOC files." :
    err== -12? msg << "(pardiso_64 only) pardiso_64 called from 32-bit library." :
               msg << "unknown error.";
    throw std::runtime_error(msg.str());
  }
  return *this;
}


int pardiso::call_pardiso(int _phase)
{
  err = 0;
  phase = _phase;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
    &m_A.nnu, &m_A.a[0], &m_A.ia[0], &m_A.ja[0],
    &perm[0], &nrhs, iparm, &msglvl, &m_b.a[0], &m_x.a[0], &err);
  return err;
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3

