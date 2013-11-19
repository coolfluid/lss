// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf

#include "mkl_dss.h" 
#include "mkl_service.h"

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "dss.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< dss, common::Component, LibLSS_MKL > Builder_MKL_dss;


dss::dss(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k, const double& _value)
  : linearsystem< double >(name)
{
  const char *nthreads = getenv("OMP_NUM_THREADS");
  int nthd = 1;
  sscanf(nthreads? nthreads:"1","%d",&nthd);
  mkl_set_num_threads(nthd);

  CFinfo  << "mkl dss: OMP_NUM_THREADS: " << nthd << (nthreads? " (set)":" (not set)") << CFendl;

  handle = NULL;
  opt  = MKL_DSS_DEFAULTS;
  sym  = MKL_DSS_NON_SYMMETRIC;
  type = MKL_DSS_INDEFINITE;

  linearsystem< double >::initialize(_size_i,_size_j,_size_k,_value);
}


dss& dss::solve()
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  nrhs = static_cast< int >(m_b.size(1));
  int err;

  if (/*1: initialize       */ (err=dss_create_(&handle, &opt))
    ||/*2: m. structure     */ (err=dss_define_structure_(&handle, &sym, &A.ia[0], &A.nnu, &A.nnu, &A.ja[0], &A.nnz))
    ||/*3: m. reordering    */ (err=dss_reorder_(&handle, &opt, NULL))
    ||/*4: m. factorization */ (err=dss_factor_real_(&handle, &type, &A.a[0]))
    ||/*5: solve            */ (err=dss_solve_real_(&handle, &opt, &m_b.a[0], &nrhs, &m_x.a[0]))
    ||/*6: deallocate       */ (err=dss_delete_(&handle, &opt))
  ){
    std::ostringstream msg;
    msg << "mkl dss: error " << err << ": ";
    err==  -1? msg << "zero pivot."      :
    err==  -2? msg << "out of memory."   :
    err==  -3? msg << "failure."         :
    err==  -4? msg << "row err."         :
    err==  -5? msg << "col err."         :
    err==  -6? msg << "too few values."  :
    err==  -7? msg << "too many values." :
    err==  -8? msg << "not square."      :
    err==  -9? msg << "state err."       :
    err== -10? msg << "invalid option."  :
    err== -11? msg << "option conflict." :
    err== -12? msg << "msg lvl err."     :
    err== -13? msg << "term lvl err."    :
    err== -14? msg << "structure err."   :
    err== -15? msg << "reorder err."     :
    err== -16? msg << "values err."      :
    err== -17? msg << "statistics invalid matrix." :
    err== -18? msg << "statistics invalid state."  :
    err== -19? msg << "statistics invalid string." :
    err== -20? msg << "reorder1 err."    :
    err== -21? msg << "preorder err."    :
    err== -22? msg << "diag err."        :
    err== -23? msg << "i32bit err."      :
    err== -24? msg << "ooc mem err."     :
    err== -25? msg << "ooc oc err."      :
    err== -26? msg << "ooc rw err."      :
               msg << "unknown error.";
    throw std::runtime_error(msg.str());
  }
  return *this;
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3

