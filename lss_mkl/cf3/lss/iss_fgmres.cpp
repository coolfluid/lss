// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "mkl_rci.h"
#include "mkl_service.h"
#include "mkl_spblas.h"

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "iss_fgmres.h"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< iss_fgmres, common::Component, LibLSS_MKL > Builder_MKL_iss_fgmres;


iss_fgmres::iss_fgmres(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k)
  : linearsystem< double >(name)
{
  environment_variable_t< int > nthreads("OMP_NUM_THREADS",1);
  CFinfo << "mkl iss_fgmres: OMP_NUM_THREADS: " << nthreads.description() << CFendl;
  mkl_set_num_threads(nthreads.value);

  linearsystem< double >::initialize(_size_i,_size_j,_size_k);
}


iss_fgmres& iss_fgmres::solve()
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  const size_t N = m_A.size(0);


  // allocate storage for the ?parm parameters
  std::vector< double > tmp( N*(2*N+1) + (N*(N+9))/2 + 1 );
  MKL_INT RCI_request;
  MKL_INT itercount = 0;
  char cvar = 'N'; // non-transposed matrix


  // initialize the solver
  dfgmres_init(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
  if (RCI_request)
    throw std::runtime_error(std::string("iss_fgmres: dfgmres_init failed."));


  // set configuration parameters and check their consistency
  iparm[ 8] = 1;      // do residual stopping test
  iparm[ 9] = 0;      // do not request for the user defined stopping test
  iparm[11] = 1;      // do the check of the norm of the next generated vector automatically
  dparm[ 0] = 1.E-3;  // relative tolerance (instead of default value 1.E-6)

  dfgmres_check(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
  if (RCI_request)
    throw std::runtime_error(std::string("iss_fgmres: dfgmres_init failed."));

  CFdebug << "iss_fgmres: information about the RCI FGMRES method:" << '\n'
          << "  iparm[ 7]=" << iparm[ 7] << ": automatic test for the maximal number of iterations will be " << (iparm[7]? "performed":"skipped") << '\n'
          << "  iparm[ 8]=" << iparm[ 8] << ": automatic residual test will be " << (iparm[8]? "performed":"skipped") << '\n'
          << "  iparm[ 9]=" << iparm[ 9] << ": user-defined stopping test " << (iparm[9]? "will be requested via RCI_request=2":"will not be requested, thus RCI_request will not take the value 2") << '\n'
          << "  iparm[10]=" << iparm[10] << ": FGMRES iterations " << (iparm[10]? "will be performed, thus the preconditioner action will be requested via RCI_request=3":"will not be performed, thus RCI_request will not take the value 3") << '\n'
          << "  iparm[11]=" << iparm[11] << ": automatic test for the norm of the next generated vector is not equal to zero up to rounding and computational errors " << (iparm[11]? "will be performed, thus RCI_request will not take the value 4":"will be skipped, thus the user-defined test will be requested via RCI_request=4")
          << CFendl;


  // reverse communication loop
  for (bool finished=false; !finished; ) {
    dfgmres(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
    switch (RCI_request) {

      case 1:
        // iterative step:
        // compute vector A*tmp[ipar[21]-1] into vector tmp[ipar[22]-1]
        // NOTE: iparm[21] and [22] contain FORTRAN style addresses
        mkl_dcsrgemv(&cvar, &A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &tmp[iparm[21] - 1], &tmp[iparm[22] - 1]);
        break;

      case 0:
        // finish step:
        // - get the FGMRES solution (computed_solution still contains initial guess)
        // - get the current iteration number
        dfgmres_get(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0], &itercount);
        finished = true;
        break;

      default:
        // this indicates failure, since RCI_request!=0
        finished = true;
        break;
    }
  }
  CFdebug << "iss_fgmres: " << (RCI_request? "failed":"success") << ", iterations: " << itercount << CFendl;


  // release internal memory
  MKL_Free_Buffers();
  return *this;
}


iss_fgmres& iss_fgmres::copy(const iss_fgmres& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  for (size_t i=0; i<128; ++i) iparm[i] = _other.iparm[i];
  for (size_t i=0; i<128; ++i) dparm[i] = _other.dparm[i];
  return *this;
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3
