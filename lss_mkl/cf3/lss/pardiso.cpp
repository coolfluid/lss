// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "mkl_pardiso.h"
#include "mkl_service.h"

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "pardiso.h"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< pardiso, common::Component, LibLSS_MKL > Builder_MKL_pardiso;


pardiso::pardiso(
    const std::string& name,
    const size_t& _size_i,
    const size_t& _size_j,
    const size_t& _size_k )
  : detail::mkl_solver_base(name)
{
  environment_variable_t< std::string >
    ooc_path("MKL_PARDISO_OOC_PATH"),
    ooc_maxc("MKL_PARDISO_OOC_MAX_CORE_SIZE"),
    ooc_swap("MKL_PARDISO_OOC_MAX_SWAP_SIZE"),
    ooc_keep("MKL_PARDISO_OOC_KEEP_FILE");

  CFdebug << "mkl pardiso: MKL_PARDISO_OOC_PATH:          " << ooc_path.description() << CFendl
          << "mkl pardiso: MKL_PARDISO_OOC_MAX_CORE_SIZE: " << ooc_maxc.description() << CFendl
          << "mkl pardiso: MKL_PARDISO_OOC_MAX_SWAP_SIZE: " << ooc_swap.description() << CFendl
          << "mkl pardiso: MKL_PARDISO_OOC_KEEP_FILE:     " << ooc_keep.description() << CFendl;

  mtype  = 1;  // real structurally symmetric matrix
  maxfct = 1;  // maximum number of numerical factorizations
  mnum   = 1;  // which factorization to use

  for (size_t i=0; i<64; ++i) iparm[i] = 0;

  // reset pt and iparm defaults
  PARDISOINIT(pt,&mtype,iparm);
  iparm[ 7] = 0;  // + max numbers of iterative refinement steps
  iparm[31] = 0;  // + [0|1] sparse direct solver or multi-recursive iterative solver

  detail::mkl_solver_base::initialize(_size_i,_size_j,_size_k);
}


pardiso::~pardiso()
{
  call_pardiso(-1,0);  // -1: termination and release of memory
}


pardiso& pardiso::solve()
{
  int err;
  if ( (err=call_pardiso(11,0)) ||  // 11: reordering and symbolic factorization
       (err=call_pardiso(22,0)) ||  // 22: numerical factorization and
       (err=call_pardiso(33,0)) )   // 33: back substitution and iterative refinement
    throw std::runtime_error(err_message(err));
  return *this;
}


pardiso& pardiso::copy(const pardiso& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  for (size_t i=0; i<64; ++i) pt   [i] = _other.pt   [i];
  for (size_t i=0; i<64; ++i) iparm[i] = _other.iparm[i];
  maxfct = _other.maxfct;
  mnum   = _other.mnum;
  mtype  = _other.mtype;
  return *this;
}


const std::string pardiso::err_message(const int& err)
{
  std::ostringstream s;
  s << "mkl pardiso error: " << err << ": ";
  err==   0? s << "(success)"          :
  err==  -1? s << "input inconsistent" :
  err==  -2? s << "not enough memory"  :
  err==  -3? s << "reordering problem" :
  err==  -4? s << "zero pivot, numerical factorization or iterative refinement problem" :
  err==  -5? s << "unclassified (internal) error"   :
  err==  -6? s << "reordering failed (matrix types 11 and 13 only)" :
  err==  -7? s << "diagonal matrix is singular"     :
  err==  -8? s << "32-bit integer overflow problem" :
  err==  -9? s << "not enough memory for OOC"       :
  err== -10? s << "error opening OOC files"         :
  err== -11? s << "read/write error with OOC files" :
  err== -12? s << "(pardiso_64 only) pardiso_64 called from 32-bit library" :
  err==-101? s << "(internal: unimplemented)" :
  err==-102? s << "(internal: null handle)"   :
  err==-103? s << "(internal: memory error)"  :
             s << "(unknown error)";
  s << '.';
  return s.str();
}


int pardiso::call_pardiso(int _phase, int _msglvl)
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  int nrhs = static_cast< int >(m_b.size(1));

  int err = 0;
  PARDISO(
    pt, &maxfct, &mnum, &mtype, &_phase,
    &A.nnu, &A.a[0], &A.ia[0], &A.ja[0],
    NULL, &nrhs, iparm, &_msglvl, &m_b.a[0], &m_x.a[0], &err );
  return err;
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3

