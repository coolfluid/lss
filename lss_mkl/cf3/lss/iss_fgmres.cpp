// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
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
#include "iss_fgmres.h"


namespace cf3 {
namespace lss {
namespace mkl {


common::ComponentBuilder< iss_fgmres, common::Component, LibLSS_MKL > Builder_MKL_iss_fgmres;


iss_fgmres::iss_fgmres(const std::string& name, const size_t& _size_i, const size_t& _size_j, const size_t& _size_k)
  : detail::solverbase(name)
{
  // reset iparm and dparm defaults
  for (size_t i=0; i<128; ++i) iparm[i] = 0;
  for (size_t i=0; i<128; ++i) dparm[i] = 0.;


  // framework scripting: options level and options
  mark_basic();

  m_pc_type    = "none";
  m_pc_refresh = true;
  m_monitor    = false;
  m_resnorm    = 1.e-12;
  opt.pc_type  = previous_opt.pc_type = NONE;
  opt.tol      = 1.e-16;
  opt.maxfil   = 1;
  options().add("PCType",      m_pc_type   ).link_to(&m_pc_type   ).mark_basic().description("preconditioner type, \"none\" (default), \"ilu0\" or \"ilut\"");
  options().add("PCRefresh",   m_pc_refresh).link_to(&m_pc_refresh).mark_basic().description("if preconditioner is to be recalculated (default true)");
  options().add("monitor",     m_monitor   ).link_to(&m_monitor   ).mark_basic().description("if each iteration should be printed (default false)");
  options().add("resnorm",     m_resnorm   ).link_to(&m_resnorm   ).mark_basic().description("user-defined stopping test, maximum residual norm (default 1.e-12)");
  options().add("ilut_tol",    opt.tol     ).link_to(&opt.tol     ).mark_basic().description("ilut only: tolerance threshold for preconditioner entries (default 1.e-16)");
  options().add("ilut_maxfil", opt.maxfil  ).link_to(&opt.maxfil  ).mark_basic().description("ilut only: maximum fill-in (half of preconditioner bandwidth (default 1)");

  iparm[ 4] = opt.iparm__4 = 150;
  iparm[ 7] = opt.iparm__7 = 1;
  iparm[ 8] = opt.iparm__8 = 0;
  iparm[ 9] = opt.iparm__9 = 1;
  iparm[11] = opt.iparm_11 = 0;
  iparm[14] = opt.iparm_14 = 150;
  iparm[30] = 1;
  options().add("maxits",    opt.iparm__4).link_to(&opt.iparm__4).mark_basic().description("maximum number of iterations to perform (default 150)");
  options().add("test_its",  opt.iparm__7).link_to(&opt.iparm__7).mark_basic().description("perform maximum iterations stopping test (default 1, yes)");
  options().add("test_res",  opt.iparm__8).link_to(&opt.iparm__8).mark_basic().description("perform residual stopping test (default 0, no)");
  options().add("test_user", opt.iparm__9).link_to(&opt.iparm__9).mark_basic().description("perform user-defined stopping test, in this case the residual norm against option 'resnorm' (default 1, yes)");
  options().add("test_norm", opt.iparm_11).link_to(&opt.iparm_11).mark_basic().description("perform test for zero norm of next generated vector (default 0, no)");
  options().add("restart",   opt.iparm_14).link_to(&opt.iparm_14).mark_basic().description("number of non-restarted iterations (default 150)");
  options().add("pc_zerodiag_subs", iparm[30]).link_to(&iparm[30]).mark_basic().description("preconditioner: perform zero diagonal element substitution (default 1, yes)");

  dparm[ 0] = opt.dparm__0 = 1.e-6;
  dparm[ 1] = opt.dparm__1 = 0.;
  dparm[30] = 1.e-14;
  dparm[31] = 1.e-12;
  options().add("rtol",   opt.dparm__0).link_to(&opt.dparm__0).mark_basic().description("relative tolerance (default 1.e-6)");
  options().add("abstol", opt.dparm__1).link_to(&opt.dparm__1).mark_basic().description("absolute tolerance (default 0.)");
  options().add("pc_zerodiag_comp",  dparm[30]).link_to(&dparm[30]).mark_basic().description("preconditioner: zero diagonal element comparison value (default 1.e-14)");
  options().add("pc_zerodiag_value", dparm[31]).link_to(&dparm[31]).mark_basic().description("preconditioner: zero diagonal element substitution value (default 1.e-12)");


  // initialize linearsystem
  solverbase::initialize(_size_i,_size_j,_size_k);
}


iss_fgmres::~iss_fgmres()
{
  // release internal memory
  MKL_Free_Buffers();
}


iss_fgmres& iss_fgmres::solve()
{
  matrix_t::matrix_compressed_t& A = m_A.compress();


  // local temporary variables
  std::vector< double >
    trvec,  // vector for preconditioner application
    rhs,    // right-hand side vector
    res;    // residual vector

  char cvar[3];
  double dvar;
  int
    RCI_request = 0,
    itercount = 0,
    inc = 1;


  CFdebug << "mkl iss_fgmres: allocate temporary storage" << CFendl;
  const size_t tmp_size =
    m_A.size(0)*(2*opt.iparm_14+1) + (opt.iparm_14*(opt.iparm_14+9))/2 + 1;
  if (tmp.size()!=tmp_size)
    tmp.assign(tmp_size,0.);  // this is really really heavy!


  // initialize the solver
  dfgmres_init(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
  if (RCI_request)
    throw std::runtime_error(err_message(RCI_request,opt.pc_type));


  // update and check configuration parameters consistency
  opt.pc_type = (m_pc_type=="ilu0"? ILU0 :
                 m_pc_type=="ilut"? ILUT : NONE );
  iparm[ 4] = opt.iparm__4;
  iparm[ 7] = opt.iparm__7;
  iparm[ 8] = opt.iparm__8;
  iparm[ 9] = opt.iparm__9;
  iparm[10] = opt.pc_type? 1:0;
  iparm[11] = opt.iparm_11;
  iparm[14] = opt.iparm_14;
  dparm[ 0] = opt.dparm__0;
  dparm[ 1] = opt.dparm__1;

  dfgmres_check(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
  if (RCI_request)
    throw std::runtime_error(err_message(RCI_request,opt.pc_type));
  CFdebug << "mkl iss_fgmres: possible RCI requests:"
          <<              " 1"
          << ( iparm[ 9]? " 2":"")
          << ( iparm[10]? " 3":"")
          << (!iparm[11]? " 4":"") << CFendl;


  /*
   * (re-)build preconditioner if necessary
   * Preconditioners may worsen the iterative convergence for arbitrary cases
   * (system matrix or initial guess), so they should be used skillfully. Haha.
   */
  CFdebug << "mkl iss_fgmres: preconditioner: "
          << (opt.pc_type==ILU0? "ilu0" :
             (opt.pc_type==ILUT? "ilut" : "none" )) << CFendl;
  RCI_request = 0;
  if ( opt.pc_type==ILU0 && (m_pc_refresh
    || opt.pc_type != previous_opt.pc_type )) {

    m_pc.nnu = A.nnu;
    m_pc.nnz = A.nnz;
    m_pc.ia.clear();
    m_pc.ja.clear();
    m_pc.a .assign( m_pc.nnz, 0.);
    dcsrilu0(&A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &m_pc.a[0], &iparm[0], &dparm[0], &RCI_request);

  }
  else if ( opt.pc_type==ILUT && (m_pc_refresh
         || opt.pc_type != previous_opt.pc_type
         || opt.tol     != previous_opt.tol
         || opt.maxfil  != previous_opt.maxfil )) {

    m_pc.nnu = A.nnu;
    m_pc.nnz = (2*opt.maxfil+1)*(A.nnu) - opt.maxfil*(opt.maxfil+1) + 1;
    m_pc.ia.assign( m_pc.nnu+1, 0 );
    m_pc.ja.assign( m_pc.nnz,   0 );
    m_pc.a .assign( m_pc.nnz,   0.);
    dcsrilut(&A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &m_pc.a[0], &m_pc.ia[0], &m_pc.ja[0], &opt.tol, &opt.maxfil, &iparm[0], &dparm[0], &RCI_request);

  }
  else if (opt.pc_type && !m_pc_refresh) {}
  else {

    opt.pc_type = NONE;
    m_pc.nnu = m_pc.nnz = 0;
    m_pc.ia.clear();
    m_pc.ja.clear();
    m_pc.a .clear();

  }
  previous_opt = opt;
  if (RCI_request)
    throw std::runtime_error(err_message(RCI_request,opt.pc_type));


  // reverse communication loop
  for (bool finished=false; !finished;) {
    dfgmres(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0]);
    switch (RCI_request) {

      case 1:
        // iterative step
        // compute vector A*tmp[iparm[21]-1] into vector tmp[iparm[22]-1]
        // NOTE: iparm[21] and [22] contain FORTRAN style addresses
        cvar[0] = 'N';
        mkl_dcsrgemv(&cvar[0], &A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &tmp[iparm[21] - 1], &tmp[iparm[22] - 1]);
        if (m_monitor)
          CFinfo << "mkl iss_fgmres: iteration " << iparm[3] << CFendl;
        break;

      case 2:
        // user-defined stopping test (check the residual norm)

        // get solution into rhs (temporary) and calculate current true residual
        iparm[12] = 1;
        cvar[0] = 'N';
        rhs.resize(A.nnu);
        res.resize(A.nnu);
        dfgmres_get(&A.nnu, &m_x.a[0], &rhs[0], &RCI_request, iparm, dparm, &tmp[0], &itercount);
        mkl_dcsrgemv(&cvar[0], &A.nnu, &A.a[0], &A.ia[0], &A.ja[0], &rhs[0], &res[0]);

        dvar = -1.;
        daxpy(&A.nnu, &dvar, &m_b.a[0], &inc, &res[0], &inc);
        dvar = dnrm2(&A.nnu, &res[0], &inc);

        finished = (dvar < m_resnorm);
        if (finished)
          RCI_request = 0;
        break;

      case 3:
        // apply preconditioner on vector tmp[iparm[21]-1] into vector
        // tmp[iparm[22]-1]
        // NOTE: iparm[21] and [22] contain FORTRAN style addresses
        trvec.assign(A.nnu,0.);

        cvar[0] = 'L';
        cvar[1] = 'N';
        cvar[2] = 'U';
        mkl_dcsrtrsv( &cvar[0], &cvar[1], &cvar[2], &A.nnu, &m_pc.a[0],
          (opt.pc_type==ILU0? &A.ia[0]:(opt.pc_type==ILUT? &m_pc.ia[0] : NULL)),
          (opt.pc_type==ILU0? &A.ja[0]:(opt.pc_type==ILUT? &m_pc.ja[0] : NULL)),
          &tmp[iparm[21] - 1], &trvec[0] );

        cvar[0] = 'U';
        cvar[1] = 'N';
        cvar[2] = 'N';
        mkl_dcsrtrsv( &cvar[0], &cvar[1], &cvar[2], &A.nnu, &m_pc.a[0],
          (opt.pc_type==ILU0? &A.ia[0]:(opt.pc_type==ILUT? &m_pc.ia[0] : NULL)),
          (opt.pc_type==ILU0? &A.ja[0]:(opt.pc_type==ILUT? &m_pc.ja[0] : NULL)),
          &trvec[0], &tmp[iparm[22] - 1] );
        break;

      case 4:
        // check the norm of the generated vector is not too small
        finished = (dparm[6] < 1.e-12);
        if (finished)
          RCI_request = 0;
        break;

      default:
        // this indicates failure if RCI_request!=0
        finished = true;
        break;

    }
  }


  // get solution if successful (x still has initial guess) and iteration number
  const bool ok(!RCI_request);
  iparm[12] = ok? 0:-1;
  dfgmres_get(&A.nnu, &m_x.a[0], &m_b.a[0], &RCI_request, iparm, dparm, &tmp[0], &itercount);
  CFinfo << "mkl iss_fgmres: " << (ok? "succeded":"failed") << ", iterations: " << itercount << CFendl;


  return *this;
}


iss_fgmres& iss_fgmres::copy(const iss_fgmres& _other)
{
  linearsystem< double >::copy(_other);
  m_A          = _other.m_A;
  m_pc         = _other.m_pc;
  m_pc_type    = _other.m_pc_type;
  m_pc_refresh = _other.m_pc_refresh;
  m_monitor    = _other.m_monitor;
  m_resnorm    = _other.m_resnorm;
  opt          = _other.opt;
  previous_opt = _other.previous_opt;
  tmp          = _other.tmp;
  for (size_t i=0; i<128; ++i) iparm[i] = _other.iparm[i];
  for (size_t i=0; i<128; ++i) dparm[i] = _other.dparm[i];
  return *this;
}


const std::string iss_fgmres::err_message(const int& err, const pc_t& _pc_type)
{
  std::ostringstream s;
  s << "mkl iss_fgmres error: " << err << ": "
    << (err==-10000?               "dfgmres_init: "  :
       (err<=-1001 && err>=-1100)? "dfgmres_check: " :
        _pc_type==ILU0?            "ILU0: " :
        _pc_type==ILUT?            "ILUT: " :
                                   "" );
  err==     0? s << "(success)" :
  err== -1001? s << "warnings occurred" :
  err== -1010? s << "routine changed some parameters to make them consistent or correct" :
  err== -1011? s << "routine changed some parameters and there are some warning messages" :
  err== -1100? s << "errors occurred" :
  err==-10000? s << "errors occurred" :

  (_pc_type==ILU0 && err==-101)? s << "at least one diagonal element is omitted from the matrix in CSR format" :
  (_pc_type==ILU0 && err==-102)? s << "matrix contains a diagonal element with the value of zero" :
  (_pc_type==ILU0 && err==-103)? s << "matrix contains a diagonal element which is so small that it could cause an overflow, or that it would cause a bad approximation to ILU0" :
  (_pc_type==ILU0 && err==-104)? s << "memory is insufficient for the internal work array" :
  (_pc_type==ILU0 && err==-105)? s << "input matrix size n is less than 0" :
  (_pc_type==ILU0 && err==-106)? s << "matrix column indices ja are not in the ascending order" :

  (_pc_type==ILUT && err==-101)? s << "number of elements in some matrix row specified in the sparse format is equal to or less than 0" :
  (_pc_type==ILUT && err==-102)? s << "value of the computed diagonal element is less than the product of the given tolerance and the current matrix row norm, and it cannot be replaced as iparm[30]=0" :
  (_pc_type==ILUT && err==-103)? s << "element ia(i+1) is less than or equal to the element ia(i)" :
  (_pc_type==ILUT && err==-104)? s << "memory is insufficient for the internal work arrays" :
  (_pc_type==ILUT && err==-105)? s << "input value of maxfil is less than 0" :
  (_pc_type==ILUT && err==-106)? s << "input matrix size n is less than 0" :
  (_pc_type==ILUT && err==-107)? s << "an element of the array ja is less than 0, or greater than n" :
  (_pc_type==ILUT && err== 101)? s << "input value maxfil is greater than or equal to n" :
  (_pc_type==ILUT && err== 102)? s << "input value tol is less than 0" :
  (_pc_type==ILUT && err== 103)? s << "absolute value of tol is greater than value of dparm[30], it can result in instability" :
  (_pc_type==ILUT && err== 104)? s << "input value dparm[30] is equal to 0, it can cause calculations to fail" :

  s << "(unknown)";
  s << '.';
  return s.str();
}


}  // namespace mkl
}  // namespace lss
}  // namespace cf3
