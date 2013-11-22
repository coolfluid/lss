// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


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
  environment_variable_t< int > nthreads("OMP_NUM_THREADS",1);
  environment_variable_t< std::string >
    licspath ("PARDISO_LIC_PATH"),
    mklserial("MKL_SERIAL");

  CFinfo << "pardiso: OMP_NUM_THREADS:  " << nthreads .description() << CFendl
         << "pardiso: PARDISO_LIC_PATH: " << licspath .description() << CFendl
         << "pardiso: MKL_SERIAL:       " << mklserial.description() << " (should be set to YES)" << CFendl;

  mtype  = 1;  // real and structurally symmetric matrix
  maxfct = 1;  // maximum number of numerical factorizations
  mnum   = 1;  // which factorization to use

  for (size_t i=0; i<64; ++i) iparm[i] = 0;
  for (size_t i=0; i<64; ++i) dparm[i] = 0.;

  // reset pt, iparm and dparm defaults
  int err = 0;
  pardisoinit_(pt,&mtype,&iparm[31],iparm,dparm,&err);
  if (err)
    throw std::runtime_error(err_message(err));
  iparm[ 2] = nthreads.value;  // + set nb. threads
  iparm[ 7] = 0;               // + max numbers of iterative refinement steps
  iparm[31] = 0;               // + [0|1] sparse direct solver or multi-recursive iterative solver

  linearsystem< double >::initialize(_size_i,_size_j,_size_k,_value);
}


pardiso::~pardiso()
{
  call_pardiso(-1,0);  // -1: termination and release of memory
}


pardiso& pardiso::solve()
{
  int err;
  if ( (err=call_pardiso_printstats()) ||  // check for matrix/vector consistency
       (err=call_pardiso(11,0))        ||  // 11: reordering and symbolic factorization
       (err=call_pardiso(22,0))        ||  // 22: numerical factorization and
       (err=call_pardiso(33,0)) )          // 33: back substitution and iterative refinement
    throw std::runtime_error(err_message(err));
  return *this;
}


pardiso& pardiso::copy(const pardiso& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  std::copy(&_other.pt   [0],&_other.pt   [0]+64,&pt   [0]);
  std::copy(&_other.dparm[0],&_other.dparm[0]+64,&dparm[0]);
  std::copy(&_other.iparm[0],&_other.iparm[0]+64,&iparm[0]);
  maxfct = _other.maxfct;
  mnum   = _other.mnum;
  mtype  = _other.mtype;
  return *this;
}


std::string pardiso::err_message(const int& err)
{
  std::ostringstream s;
  s << "pardiso: error " << err << ": ";
  err==   0? s << "(success)" :
  err==  -1? s << "input inconsistent" :
  err==  -2? s << "not enough memory"  :
  err==  -3? s << "reordering problem" :
  err==  -4? s << "zero pivot, numerical factorization or iterative refinement problem" :
  err==  -5? s << "unclassified (internal) error"     :
  err==  -6? s << "preordering failed (matrix types 11, 13 only)" :
  err==  -7? s << "diagonal matrix problem"           :
  err==  -8? s << "32-bit integer overflow problem"   :
  err== -10? s << "no license file pardiso.lic found" :
  err== -11? s << "license is expired"                :
  err== -12? s << "wrong username or hostname"        :
  err==-100? s << "reached maximum number of Krylov-subspace iteration in iterative solver"     :
  err==-101? s << "no sufficient convergence in Krylov-subspace iteration within 25 iterations" :
  err==-102? s << "error in Krylov-subspace iteration"      :
  err==-103? s << "break-down in Krylov-subspace iteration" :
             s << "(unknown error)";
  s << '.';
  return s.str();
}


int pardiso::call_pardiso(int _phase, int _msglvl)
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  int nrhs = static_cast< int >(m_b.size(1));

  int err = 0;
  pardiso_(
    pt, &maxfct, &mnum, &mtype, &_phase,
    &A.nnu, &A.a[0], &A.ia[0], &A.ja[0],
    NULL, &nrhs, iparm, &_msglvl, &m_b.a[0], &m_x.a[0], &err, dparm );
  return err;
}


int pardiso::call_pardiso_printstats()
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  int nrhs = static_cast< int >(m_b.size(1));

  int err = 0;
  pardiso_printstats_(
    &mtype,
    &A.nnu, &A.a[0], &A.ia[0], &A.ja[0],
    &nrhs, &m_b.a[0], &err );
  return err;
}


}  // namespace lss
}  // namespace cf3

