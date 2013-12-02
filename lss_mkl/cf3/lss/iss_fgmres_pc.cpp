// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "mkl_rci.h"
#include "mkl_service.h"
#include "mkl_spblas.h"
#include "mkl_blas.h"

#include "common/Builder.hpp"
#include "common/Log.hpp"
#include "iss_fgmres_pc.h"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< iss_fgmres_pc, common::Component, LibLSS_MKL > Builder_MKL_iss_fgmres_pc;


iss_fgmres_pc::iss_fgmres_pc(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k)
  : linearsystem< double >(name),
    m_pc_type(ILU0)
{
  environment_variable_t< int > nthreads("OMP_NUM_THREADS",1);
  CFinfo << "mkl iss_fgmres_pc: OMP_NUM_THREADS: " << nthreads.description() << CFendl;
  mkl_set_num_threads(nthreads.value);

  // reset iparm and dparm
  for (size_t i=0; i<128; ++i) iparm[i] = 0;
  for (size_t i=0; i<128; ++i) dparm[i] = 0.;

  linearsystem< double >::initialize(_size_i,_size_j,_size_k);
}


iss_fgmres_pc& iss_fgmres_pc::solve()
{
  matrix_t::matrix_compressed_t& A = m_A.compress();
  const size_t N = m_A.size(0);


  // allocate storage for the ?parm parameters
  std::vector< double >
      tmp( N*(2*N+1) + (N*(N+9))/2 + 1 ),
      trvec,
      b,  // temporary right-hand side vector
      r;  // temporary residual vector
  MKL_INT RCI_request;
  MKL_INT itercount = 0;
  MKL_INT i = 1;
  char cvar1;
  char cvar2;
  char cvar3;
  double dvar;

  // ILUT-specific
  double  tol    = 1.e-6;  // tolerance threshold for preconditioner entries
  MKL_INT maxfil = 1;      // maximum fill-in (half of preconditioner bandwidth)


  // initialize the solver
  dfgmres_init(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
  if (RCI_request)
    throw std::runtime_error(std::string("iss_fgmres_pc: dfgmres_init failed."));


  // set configuration parameters
  iparm[ 1] = 6;  // output of error messages to the screen
  iparm[ 5] = 1;  // allow output of errors
  iparm[ 6] = 1;  // output warn messages if any and continue
  iparm[ 7] = 0;  // do the preconditioned iterations of FGMRES method
  iparm[10] = 1;  // do the preconditioner call
  iparm[14] = 2;  // do not do the stopping test for the maximal number of iterations
  iparm[30] = 0;  // abort preconditioner calculations if routine meets zero/small diagonal element

  dparm[ 0] = 1.e-3;  // set relative tolerance to 1.e-3 (default is 1.e-6)

//  iparm[ 8] = 1;      // do residual stopping test
//  iparm[ 9] = 0;      // do not request for the user defined stopping test
//  iparm[11] = 1;      // do the check of the norm of the next generated vector automatically


  /*
   * build preconditioner
   * Preconditioners may worsen the iterative convergence for arbitrary cases
   * (system matrix or initial guess), so they should be used skillfully. Haha.
   * NOTE: these routine use some i/dparm parameters from DFGMRES_INIT.
   */
  int ierr = 0;
  switch (m_pc_type) {

    case ILU0:
      iparm[30] = 1;       // do not abort preconditioner calculations if routine meets zero diagonal element (use dpar[31] instead)
      dparm[30] = 1.e-20;  // small value to compare a diagonal entry with it
      dparm[31] = 1.e-16;  // preconditioner diagonal target value if it is small as compared to dpar[30], change it to this rather than abort preconditioner calculation.

      m_pc.a.assign(static_cast< size_t >(A.nnu),0.);
      dcsrilu0(&A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &m_pc.a[0], &iparm[0], &dparm[0], &ierr);
      break;

    case ILUT:
      iparm[30] = 1;
      dparm[30] = 1.e-5;

      m_pc.ia.assign(A.ia.size(),0);
      m_pc.ja.assign( (2*maxfil+1)*N-maxfil*(maxfil+1)+1, 0);
      m_pc.a .assign(m_pc.ja.size(),0.);
      dcsrilut(&A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &m_pc.a[0], &m_pc.ia[0], &m_pc.ja[0], &tol, &maxfil, &iparm[0], &dparm[0], &ierr);
      break;

    default:
      throw std::runtime_error(std::string("iss_fgmres_pc: unsupported preconditioner (only ILU0 and ILUT are available)."));

  }
  if (ierr)
    throw std::runtime_error(std::string("iss_fgmres_pc: preconditioner calculation error (ierr)")); // FIXME


  // check configuration parameters consistency
  dfgmres_check(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
  if (RCI_request)
    throw std::runtime_error(std::string("iss_fgmres_pc: dfgmres_init failed."));
  CFdebug << "iss_fgmres_pc: information about the RCI FGMRES (ILU0/ILUT) method:" << '\n'
          << "  iparm[ 7]=" << iparm[ 7] << ": automatic test for the maximal number of iterations will be " << (iparm[7]? "performed":"skipped") << '\n'
          << "  iparm[ 8]=" << iparm[ 8] << ": automatic residual test will be " << (iparm[8]? "performed":"skipped") << '\n'
          << "  iparm[ 9]=" << iparm[ 9] << ": user-defined stopping test " << (iparm[9]? "will be requested via RCI_request=2":"will not be requested, thus RCI_request will not take the value 2") << '\n'
          << "  iparm[10]=" << iparm[10] << ": preconditioned FGMRES iterations " << (iparm[10]? "will be performed, thus the preconditioner action will be requested via RCI_request=3":"will not be performed, thus RCI_request will not take the value 3") << '\n'
          << "  iparm[11]=" << iparm[11] << ": automatic test for the norm of the next generated vector is not equal to zero up to rounding and computational errors " << (iparm[11]? "will be performed, thus RCI_request will not take the value 4":"will be skipped, thus the user-defined test will be requested via RCI_request=4")
          << CFendl;


#if 0
  // reverse communication loop
  for (bool finished=false; !finished;) {
    dfgmres(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
    switch (RCI_request) {

      case 0:
        // finish step:
        // - get the FGMRES solution (x still contains initial guess)
        // - get the current iteration number
        iparm[12] = 0;
        dfgmres_get(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0], &itercount);
        finished = true;
        break;

      case 1:
        // iterative step:
        // compute vector A*tmp[ipar[21]-1] into vector tmp[ipar[22]-1]
        // NOTE: iparm[21] and [22] contain FORTRAN style addresses
        cvar1 = 'N';
        mkl_dcsrgemv(&cvar1, &A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &tmp[iparm[21] - 1], &tmp[iparm[22] - 1]);
        break;

      case 2:
        // user-defined stopping test (checking the residual)

        /*
         * NOTE: from this point vector b[N] is no longer containing the right-hand
         * side of the problem! It contains the current FGMRES approximation to the
         * solution. If you need to keep the right-hand side, save it in some other
         * vector before the call to dfgmres routine. Here we saved it in vector
         * rhs[N]. The vector b is used instead of rhs to preserve the
         * original right-hand side of the problem and guarantee the proper
         * restart of FGMRES method. Vector b will be altered when computing the
         * residual stopping criterion!
         *
         * Request to the dfgmres_get routine to put the solution into b[N] via ipar[12]
         * WARNING: beware that the call to dfgmres_get routine with ipar[12]=0 at this
         * stage may destroy the convergence of the FGMRES method, therefore, only
         * advanced users should exploit this option with care
         */

          // Get current solution into b[N] and calculate current true residual
          iparm[12] = 1;
          cvar1 = 'N';
          b.resize(N);
          r.resize(N);
          dfgmres_get(&A.nnu, &m_x.a[0], &b[0], &RCI_request, iparm, dparm, &tmp[0], &itercount);
          mkl_dcsrgemv(&cvar1, &A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &b[0], &r[0]);

          dvar = -1.;
          daxpy(&A.nnu, &dvar, &m_b.a[0], &i, &r[0], &i);
          dvar = dnrm2(&A.nnu, &r[0], &i);

          finished = (dvar < 1.e-3);
          if (finished)
            RCI_request = 0;
        break;

      case 3:
        // apply the preconditioner on the vector tmp[iparm[21]-1] and put
        // result in vector tmp[iparm[22]-1]
        // NOTE: iparm[21] and [22] contain FORTRAN style addresses
        trvec.assign(N,0.);

        cvar1 = 'L';
        cvar2 = 'N';
        cvar3 = 'U';
        mkl_dcsrtrsv( &cvar1, &cvar2, &cvar3, &A.nnu, &m_pc.a[0],
          (m_pc_type==ILU0? &A.ia[0]:(m_pc_type==ILUT? &m_pc.ia[0] : NULL)),
          (m_pc_type==ILU0? &A.ja[0]:(m_pc_type==ILUT? &m_pc.ja[0] : NULL)),
          &tmp[iparm[21] - 1], &trvec[0] );

        cvar1 = 'U';
        cvar2 = 'N';
        cvar3 = 'N';
        mkl_dcsrtrsv( &cvar1, &cvar2, &cvar3, &A.nnu, &m_pc.a[0],
          (m_pc_type==ILU0? &A.ia[0]:(m_pc_type==ILUT? &m_pc.ia[0] : NULL)),
          (m_pc_type==ILU0? &A.ja[0]:(m_pc_type==ILUT? &m_pc.ja[0] : NULL)),
          &trvec[0], &tmp[iparm[22] - 1] );
        break;

      case 4:
        // check if the norm of the next generated vector (in dparm[6]) is not
        // zero up to rounding and computational errors
        finished = (dparm[6] < 1.e-12);
        if (finished)
          RCI_request = 0;
        break;

      default:
        // this indicates failure, since RCI_request!=0
        finished = true;
        break;
    }
  }
  CFdebug << "iss_fgmres_pc: " << (RCI_request? "failed":"succeded") << ", iterations: " << itercount << CFendl;
#endif


  // release internal memory
  MKL_Free_Buffers();
  return *this;
}


iss_fgmres_pc& iss_fgmres_pc::copy(const iss_fgmres_pc& _other)
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
