// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "pardiso.hpp"


// Prototypes for external function calls
extern "C" {
  void pardiso_printstats_(int *, int *, double *, int *, int *, int *, double *, int *);
  int pardisoinit_(void *, int *, int *, int *, double *, int *);
  int pardiso_(void *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *, double *);  
}


namespace cf3 {
namespace lss {


common::ComponentBuilder< pardiso, common::Component, LibLSS_PARDISO > Builder_pardiso;


pardiso::pardiso(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k, const double& _value)
  : linearsystem< double >(name)
{

  char* nthreads = getenv("OMP_NUM_THREADS");
  sscanf(nthreads? nthreads:"1","%d",&iparm[2]);
  CFinfo << "pardiso: OMP_NUM_THREADS: " << iparm[2] << (nthreads? " (set)":" (not set)") << CFendl;

  for (int i=0; i<64; ++i)  iparm[i] = 0;
  for (int i=0; i<64; ++i)  dparm[i] = 0.;

  iparm[ 7] = 0;  // max numbers of iterative refinement steps
  iparm[31] = 0;  // [0|1] sparse direct solver or multi-recursive iterative solver
  maxfct    = 1;  // maximum number of numerical factorizations
  mnum      = 1;  // which factorization to use
  msglvl    = 1;  // message level: output statistical information
  mtype     = 1;  // real structurally symmetric matrix

  if (call_pardiso_init()) {
    std::ostringstream msg;
    msg << "pardiso: pardisoinit error " << err << ": ";
    err== -10? msg << "no license file pardiso.lic found." :
    err== -11? msg << "license is expired."                :
    err== -12? msg << "wrong username or hostname."        :
               msg << "unknown error.";
    throw std::runtime_error(msg.str());
  }

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
     (call_pardiso_printstats() ||  // check for matrix/vector consistency
      call_pardiso(11)          ||  // 11: reordering and symbolic factorization
      call_pardiso(22)          ||  // 22: numerical factorization and
      call_pardiso(33))             // 33: back substitution and iterative refinement
#else
     (call_pardiso(13))
#endif
  {
    std::ostringstream msg;
    msg << "pardiso: phase " << phase << " error " << err << ": ";
    err==  -1? msg << "input inconsistent." :
    err==  -2? msg << "not enough memory."  :
    err==  -3? msg << "reordering problem." :
    err==  -4? msg << "zero pivot, numerical factorization or iterative refinement problem." :
    err==  -5? msg << "unclassified (internal) error."     :
    err==  -6? msg << "preordering failed (matrix types 11, 13 only)." :
    err==  -7? msg << "diagonal matrix problem."           :
    err==  -8? msg << "32-bit integer overflow problem."   :
    err==-100? msg << "reached maximum number of Krylov-subspace iteration in iterative solver."     :
    err==-101? msg << "no sufficient convergence in Krylov-subspace iteration within 25 iterations." :
    err==-102? msg << "error in Krylov-subspace iteration."      :
    err==-103? msg << "break-down in Krylov-subspace iteration." :
               msg << "unknown error.";
    throw std::runtime_error(msg.str());
  }
  return *this;
}


int pardiso::call_pardiso_printstats()
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  err = 0;
  phase = 0;
  pardiso_printstats_(
    &mtype,
    &A.nnu, &A.a[0], &A.ia[0], &A.ja[0],
    &nrhs, &m_b.a[0], &err );
  return err;
}


int pardiso::call_pardiso_init()
{
  err = 0;
  phase = 0;
  pardisoinit_(pt,&mtype,&iparm[31],iparm,dparm,&err);
  return err;
}


int pardiso::call_pardiso(int _phase)
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  err = 0;
  phase = _phase;
  pardiso_(
    pt, &maxfct, &mnum, &mtype, &phase,
    &A.nnu, &A.a[0], &A.ia[0], &A.ja[0],
    NULL, &nrhs, iparm, &msglvl, &m_b.a[0], &m_x.a[0], &err, dparm );
  return err;
}


}  // namespace lss
}  // namespace cf3

