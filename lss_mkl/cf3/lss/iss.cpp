// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf
#include <cmath>

#include "mkl_rci.h" 
#include "mkl_blas.h"
#include "mkl_spblas.h"

#include "common/Builder.hpp"
#include "iss.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< iss, common::Component, LibLSS_MKL > Builder_MKL_iss;


iss::iss(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k, const double& _value) : linearsystem_t(name)
{
  char* nthreads = getenv("OMP_NUM_THREADS");
  int nthd = 1;
  sscanf(nthreads? nthreads:"1","%d",&nthd);
  mkl_set_num_threads(nthd);

  CFinfo  << "mkl dss: OMP_NUM_THREADS: " << nthd << (nthreads? " (set)":" (not set)") << CFendl;

  linearsystem_t::initialize(_size_i,_size_j,_size_k,_value);
}


iss& iss::solve()
{
//  nrhs = static_cast< int >(b().size(1));
  if (false) {
    std::ostringstream msg;
//    msg << "mkl iss: phase " << phase << " error " << err << ".";
    throw std::runtime_error(msg.str());
  }







  /*
   *  MKL RCI Flexible Generalized Minimal RESidual method with ILU0 Preconditioner
   *  Example program for solving non-degenerate system of equations, showing
   *  how preconditioner reduces the number of iterations.
   */

  #define N 4

  /*---------------------------------------------------------------------------
  /* Define arrays for the upper triangle of the coefficient matrix
  /* Compressed sparse row storage is used for sparse representation
  /*---------------------------------------------------------------------------*/
    int ia[5] = { 1, 4, 7, 10, 13 };
    int ja[12] = { 1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4 };
    double A[12] = { 4., -1., -1., -1., 4., -1., -1., 4., -1., -1., -1., 4. };

  /*---------------------------------------------------------------------------
  /* Allocate storage for the ?par parameters and the solution/rhs/residual vectors
  /*---------------------------------------------------------------------------*/
    double tmp[N * (2 * N + 1) + (N * (N + 9)) / 2 + 1];
    double trvec[N];
    double bilu0[12];
    double rhs[N];
    double b[N];
    double x[N];
    double residual[N];

    int matsize = 12;
    int incx = 1;
    int ref_nit = 2;
    double ref_norm2 = 7.772387E+0;
    double nrm2;

    /* Some additional variables to use with the RCI (P)FGMRES solver*/
    int itercount;
    int err = 0;
    int RCI_request;
    int ivar;
    double dvar;
    char cvar;
    char cvar1;
    char cvar2;



    // Initialize the solver
    dfgmres_init (&ivar, x, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request != 0)
      goto FAILED;

  /*---------------------------------------------------------------------------
  /* Calculate ILU0 preconditioner.
  /*                      !ATTENTION!
  /* DCSRILU0 routine uses some IPAR, DPAR set by DFGMRES_INIT routine.
  /* Important for DCSRILU0 default entries set by DFGMRES_INIT are
  /* ipar[1] = 6 - output of error messages to the screen,
  /* ipar[5] = 1 - allow output of errors,
  /* ipar[30]= 0 - abort DCSRILU0 calculations if routine meets zero diagonal element.
  /*
  /* If ILU0 is going to be used out of MKL FGMRES context, than the values
  /* of ipar[1], ipar[5], ipar[30], dpar[30], and dpar[31] should be user
  /* provided before the DCSRILU0 routine call.
  /*
  /* In this example, specific for DCSRILU0 entries are set in turn:
  /* ipar[30]= 1 - change small diagonal value to that given by dpar[31],
  /* dpar[30]= 1.E-20 instead of the default value set by DFGMRES_INIT.
  /*                  It is a small value to compare a diagonal entry with it.
  /* dpar[31]= 1.E-16 instead of the default value set by DFGMRES_INIT.
  /*                  It is the target value of the diagonal value if it is
  /*                  small as compared to dpar[30] and the routine should change
  /*                  it rather than abort DCSRILU0 calculations.
  /*---------------------------------------------------------------------------*/

    ipar[30] = 1;
    dpar[30] = 1.E-20;
    dpar[31] = 1.E-16;

    dcsrilu0(&ivar, A, ia, ja, bilu0, ipar, dpar, &err);
    nrm2 = dnrm2 (&matsize, bilu0, &incx);
    if (err != 0) {
      printf ("Preconditioner dcsrilu0 has returned the ERROR code %d\n", err);
      goto FAILED1;
    }

    /*---------------------------------------------------------------------------
    /* Set the desired parameters:
    /*
    /* LOGICAL parameters:
    /* do not do the stopping test for the maximal number of iterations
    /* do the Preconditioned iterations of FGMRES method
    /*
    /* NOTE. Preconditioner may increase the number of iterations for an
    /* arbitrary case of the system and initial guess and even ruin the
    /* convergence. It is user's responsibility to use a suitable preconditioner
    /* and to apply it skillfully.
    /*---------------------------------------------------------------------------*/
    ipar[14] = 2;         // do the restart after 2 iterations
    ipar[ 7] = 0;         // automatic test for the maximal number of iterations will be performed (!=0) or skipped (==0)
    ipar[ 8] = ipar[ 8];  // the automatic residual test will be performed (!=0) or skipped (==0)
    ipar[ 9] = ipar[ 9];  // user-defined stopping test will be requested (!=0, via RCI_request=2) or not (==0, via RCI_request=2)
    ipar[10] = 1;         // set parameter ipar[10] for preconditioner call (!=0, the preconditioner action will be requested via RCI_request=3) or not (==0, RCI_request will not take the value 3)
    ipar[11] = ipar[11];  // automatic test for the norm of the next generated vector is not equal to zero up to rounding and computational errors will be performed (!=0, RCI_request will not take the value 4) or skipped (==0, the user-defined test will be requested via RCI_request)
    dpar[ 0] = 1.0E-3;    // set the relative tolerance to 1.0D-3 instead of default value 1.0D-6

    // Check the correctness and consistency of the newly set parameters
    dfgmres_check (&ivar, x, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request != 0)
      goto FAILED;


    // Compute the solution by RCI (P)FGMRES solver with preconditioning
    // (Reverse Communication starts here)
  ONE:dfgmres (&ivar, x, rhs, &RCI_request, ipar, dpar, tmp);

    if (RCI_request == 0) {
      // If RCI_request=0, then the solution was found with the required precision
      goto COMPLETE;
    }
    if (RCI_request == 1) {
      /*---------------------------------------------------------------------------
      /* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]
      /* NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
      /* therefore, in C code it is required to subtract 1 from them to get C style addresses
      /*---------------------------------------------------------------------------*/
        mkl_dcsrgemv (&cvar, &ivar, A, ia, ja, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);
        goto ONE;
      }
    if (RCI_request == 2) {
      /*---------------------------------------------------------------------------
      /* If RCI_request=2, then do the user-defined stopping test
      /* The residual stopping test for the computed solution is performed here
      /*---------------------------------------------------------------------------
      /* NOTE: from this point vector b[N] is no longer containing the right-hand
      /* side of the problem! It contains the current FGMRES approximation to the
      /* solution. If you need to keep the right-hand side, save it in some other
      /* vector before the call to dfgmres routine. Here we saved it in vector
      /* rhs[N]. The vector b is used instead of rhs to preserve the
      /* original right-hand side of the problem and guarantee the proper
      /* restart of FGMRES method. Vector b will be altered when computing the
      /* residual stopping criterion!
      /*---------------------------------------------------------------------------*/
        /* Request to the dfgmres_get routine to put the solution into b[N] via ipar[12]
           --------------------------------------------------------------------------------
           WARNING: beware that the call to dfgmres_get routine with ipar[12]=0 at this
           stage may destroy the convergence of the FGMRES method, therefore, only
           advanced users should exploit this option with care */

        ipar[12] = 1;
        /* Get the current FGMRES solution in the vector b[N] */
        dfgmres_get (&ivar, x, b, &RCI_request, ipar, dpar, tmp, &itercount);
        /* Compute the current true residual via MKL (Sparse) BLAS routines */
        mkl_dcsrgemv (&cvar, &ivar, A, ia, ja, b, residual);

        dvar = -1.0E0;
        int i = 1;
        daxpy (&ivar, &dvar, rhs, &i, residual, &i);
        dvar = dnrm2 (&ivar, residual, &i);
        if (dvar < 1.0E-3) goto COMPLETE;
        else               goto ONE;
    }
    if (RCI_request == 3) {
        /*---------------------------------------------------------------------------
        /* If RCI_request=3, then apply the preconditioner on the vector
        /* tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]
        /*---------------------------------------------------------------------------
        /* NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
        /* therefore, in C code it is required to subtract 1 from them to get C style
        /* addresses
        /* Here is the recommended usage of the result produced by ILU0 routine
        /* via standard MKL Sparse Blas solver routine mkl_dcsrtrsv.
        /*---------------------------------------------------------------------------*/
        cvar1 = 'L';
        cvar = 'N';
        cvar2 = 'U';
        mkl_dcsrtrsv (&cvar1, &cvar, &cvar2, &ivar, bilu0, ia, ja, &tmp[ipar[21] - 1], trvec);
        cvar1 = 'U';
        cvar = 'N';
        cvar2 = 'N';
        mkl_dcsrtrsv (&cvar1, &cvar, &cvar2, &ivar, bilu0, ia, ja, trvec, &tmp[ipar[22] - 1]);
        goto ONE;
      }
    if (RCI_request == 4) {
      /*---------------------------------------------------------------------------
      /* If RCI_request=4, then check if the norm of the next generated vector is
      /* not zero up to rounding and computational errors. The norm is contained
      /* in dpar[6] parameter
      /*---------------------------------------------------------------------------*/
        if (dpar[6] < 1.0E-12) goto COMPLETE;
        else                   goto ONE;
      }
    /*---------------------------------------------------------------------------
    /* If RCI_request=anything else, then dfgmres subroutine failed to compute
    /*---------------------------------------------------------------------------*/
    else
      {
        goto FAILED;
      }


    /*---------------------------------------------------------------------------
    /* Reverse Communication ends here
    /* Get the current iteration number and the FGMRES solution (DO NOT FORGET to
    /* call dfgmres_get routine as computed_solution is still containing
    /* the initial guess!). Request to dfgmres_get to put the solution
    /* into vector computed_solution[N] via ipar[12]
    /*---------------------------------------------------------------------------*/
  COMPLETE:ipar[12] = 0;
    dfgmres_get (&ivar, x, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
    // x: contains the solution
    // itercount: number of iterations

    // Release internal memory to avoid memory leaks
    MKL_Free_Buffers();

    if (itercount == ref_nit && fabs(ref_norm2 - nrm2) < 1.e-6) {}
    else
      {  // failed: probably preconditioner is incorrect
        printf ("Either the preconditioner norm %e differs from the ", nrm2);
        printf ("expected norm %e\n", ref_norm2);
        printf ("and/or the number of iterations %d differs from the ",  itercount);
        printf ("expected number %d\n", ref_nit);
        //return 1;
      }

  FAILED:
    printf ("The solver has returned the ERROR code %d \n", RCI_request);
  FAILED1:
    printf
      ("-------------------------------------------------------------------\n");
    printf ("Unfortunately, FGMRES+ILU0 C example has FAILED\n");
    printf("-------------------------------------------------------------------\n");

    // Release internal memory to avoid memory leaks
    MKL_Free_Buffers();

  return *this;
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3

