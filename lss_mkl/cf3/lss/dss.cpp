// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "mkl_dss.h"
#include "mkl_service.h"

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "dss.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< dss, common::Component, LibLSS_MKL > Builder_MKL_dss;


dss::dss(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k)
  : linearsystem< double >(name)
{
  environment_variable_t< int > nthreads("OMP_NUM_THREADS",1);
  CFinfo << "mkl dss: OMP_NUM_THREADS: " << nthreads.description() << CFendl;
  mkl_set_num_threads(nthreads.value);

  handle = NULL;
  for (int i=0; i<_all_phases; ++i)
    opts[i] = MKL_DSS_DEFAULTS;
  opts[ _structure ] += MKL_DSS_SYMMETRIC_STRUCTURE;
  opts[ _factor    ] += MKL_DSS_INDEFINITE;

  // dss: initialize
  int err;
  if (err=dss_create_(&handle, &opts[_create]))
    throw std::runtime_error(err_message(err));

  linearsystem< double >::initialize(_size_i,_size_j,_size_k);
}

dss::~dss()
{
  // dss: deallocate
  int err;
  if (err=dss_delete_(&handle, &opts[_delete]))
    throw std::runtime_error(err_message(err));
}


dss& dss::solve()
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  int nrhs = static_cast< int >(m_b.size(1));
  int err;
  if ((err=dss_define_structure_(&handle, &opts[_structure], &A.ia[0], &A.nnu, &A.nnu, &A.ja[0], &A.nnz))
   || (err=dss_reorder_         (&handle, &opts[_reorder],   NULL))
   || (err=dss_factor_real_     (&handle, &opts[_factor],    &A.a[0]))
   || (err=dss_solve_real_      (&handle, &opts[_solve],     &m_b.a[0], &nrhs, &m_x.a[0])) )
    throw std::runtime_error(err_message(err));
  return *this;
}


dss& dss::copy(const dss& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  handle = _other.handle;
  for (int i=0; i<_all_phases; ++i)
    opts[i] = _other.opts[i];
  return *this;
}


std::string dss::err_message(const int& err)
{
  std::ostringstream s;
  s << "mkl dss error: " << err << ": ";
  err==   0? s << "(success)"       :
  err==  -1? s << "zero pivot"      :
  err==  -2? s << "out of memory"   :
  err==  -3? s << "failure"         :
  err==  -4? s << "row err"         :
  err==  -5? s << "col err"         :
  err==  -6? s << "too few values"  :
  err==  -7? s << "too many values" :
  err==  -8? s << "not square"      :
  err==  -9? s << "state err"       :
  err== -10? s << "invalid option"  :
  err== -11? s << "option conflict" :
  err== -12? s << "msg lvl err"     :
  err== -13? s << "term lvl err"    :
  err== -14? s << "structure err"   :
  err== -15? s << "reorder err"     :
  err== -16? s << "values err"      :
  err== -17? s << "statistics invalid matrix" :
  err== -18? s << "statistics invalid state"  :
  err== -19? s << "statistics invalid string" :
  err== -20? s << "reorder1 err"    :
  err== -21? s << "preorder err"    :
  err== -22? s << "diag err"        :
  err== -23? s << "i32bit err"      :
  err== -24? s << "ooc mem err"     :
  err== -25? s << "ooc oc err"      :
  err== -26? s << "ooc rw err"      :
             s << "(unknown error)";
  return s.str();
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3

