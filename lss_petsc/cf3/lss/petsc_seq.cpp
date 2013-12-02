// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf

#include "common/Builder.hpp"
#include "petsc_seq.h"


namespace cf3 {
namespace lss {


common::ComponentBuilder< petsc_seq, common::Component, LibLSS_PETSC > Builder_petsc_seeq;


petsc_seq::petsc_seq(const std::string& name,
    const size_t& _size_i,
    const size_t& _size_j,
    const size_t& _size_k)
  : linearsystem< double >(name)
{
  // initialize solver/preconditioner
  int argc = 0;
  char **args = NULL;
  char help[] = "";
  PetscInitialize(&argc,&args,(char*)0,help);

  opt.ksptype = KSPGMRES;
  opt.pctype  = PCASM;
  opt.monitor = false;
  opt.ovl = 1;
  PetscErrorCode err = 0;
  if ( (err=KSPCreate(PETSC_COMM_SELF,&ksp)) ||
       (err=KSPGetPC(ksp,&pc)) ||
       (err=KSPSetFromOptions(ksp)) ||
       (err=KSPSetInitialGuessNonzero(ksp,PETSC_TRUE)) ||
       (err=KSPGetTolerances(ksp,&opt.rtol,&opt.abstol,&opt.dtol,&opt.maxits)) ||
       (err=KSPGMRESGetRestart(ksp,&opt.restart)) )
    throw std::runtime_error(err_message((int) err,"initialize solver/preconditioner"));


  // framework scripting: options level and options
  mark_basic();
  options().add("KSPType",opt.ksptype).link_to(&opt.ksptype).mark_basic().description("Krylov solver type (for instance \"gmres\", \"cg\", \"bicg\", ...)");
  options().add("PCType", opt.pctype ).link_to(&opt.pctype ).mark_basic().description("preconditioner type (for instance \"jacobi\", \"bjacobi\", \"ilu\", \"asm\", ...)");
  options().add("monitor",opt.monitor).link_to(&opt.monitor).mark_basic().description("if each iteration residual norm should be printed (KSPMonitorDefault)");
  options().add("restart",opt.restart).link_to(&opt.restart).mark_basic().description("number of iterations at which \"gmres\", \"fgmres\" and \"lgmres\" restarts (number of Krylov subspaces)");
  options().add("ovl",    opt.ovl    ).link_to(&opt.ovl    ).mark_basic().description("overlap between a pair of subdomains for the additive Schwarz preconditioners (\"asm\", \"gasm\")");
  options().add("maxits", opt.maxits ).link_to(&opt.maxits ).mark_basic().description("maximum number of iterations to use");
  options().add("rtol",   opt.rtol   ).link_to(&opt.rtol   ).mark_basic().description("the relative convergence tolerance (relative decrease in the residual norm)");
  options().add("abstol", opt.abstol ).link_to(&opt.abstol ).mark_basic().description("the absolute convergence tolerance (absolute size of the residual norm)");
  options().add("dtol",   opt.dtol   ).link_to(&opt.dtol   ).mark_basic().description("the divergence tolerance (amount residual can increase before KSPDefaultConverged() concludes that the method is diverging)");


  // initialize linearsystem
  linearsystem< double >::initialize(_size_i,_size_j,_size_k);
}


petsc_seq::~petsc_seq()
{
  PetscErrorCode err = 0;
  if ( (err=KSPDestroy(&ksp)) ||
       (err=PetscFinalize()) )
    throw std::runtime_error(err_message((int) err,"destruction"));
}


petsc_seq& petsc_seq::solve()
{
  if (!m_A.m_size.is_square_size())
    throw std::runtime_error("petsc_seq: system matrix must be square.");
  matrix_t::matrix_compressed_t& A = m_A.compress();


  // set matrix/vectors
  PetscInt n = static_cast< PetscInt >(m_A.size(1));
  Mat Ap;
  Vec bp;
  Vec xp;
  PetscErrorCode err = 0;
  if ( (err=MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,n,n,&A.ia[0],&A.ja[0],&A.a[0],&Ap)) ||
       (err=VecCreateSeqWithArray(PETSC_COMM_SELF,1,n,&m_b.a[0],&bp)) ||
       (err=VecCreateSeqWithArray(PETSC_COMM_SELF,1,n,&m_x.a[0],&xp)) )
    throw std::runtime_error(err_message((int) err,"set matrix/vectors"));


  // set solver/preconditioner options
  // NOTE: here the system matrix serves as preconditioning matrix
  CFdebug << "petsc_seq: solver/preconditioner options:" << '\n'
          << "  ksptype:  " << opt.ksptype << '\n'
          << "  pctype:   " << opt.pctype  << '\n'
          << "  monitor:  " << (opt.monitor?"true":"false") << '\n'
          << "  restart:  " << opt.restart << '\n'
          << "  ovl:      " << opt.ovl     << '\n'
          << "  maxits:   " << opt.maxits  << '\n'
          << "  rtol:     " << opt.rtol    << '\n'
          << "  abstol:   " << opt.abstol  << '\n'
          << "  dtol:     " << opt.dtol    << CFendl;
  if ( (err=KSPSetType(ksp,opt.ksptype.c_str())) ||
       (err=KSPSetTolerances(ksp,opt.rtol,opt.abstol,opt.dtol,opt.maxits)) ||
       (err=KSPSetOperators(ksp,Ap,Ap,DIFFERENT_NONZERO_PATTERN)) ||
       (err=PCSetType(pc,opt.pctype.c_str())) ||
       (err=PCASMSetOverlap(pc,opt.ovl)) ||
       (err=KSPGMRESSetRestart(ksp,opt.restart)) ||
       (opt.monitor && (err=KSPMonitorSet(ksp,&KSPMonitorDefault,NULL,NULL))) )
    throw std::runtime_error(err_message((int) err,"set solver/preconditioner options"));





  // solve
  if (err=KSPSolve(ksp,bp,xp))
    throw std::runtime_error(err_message((int) err,"solve"));
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  CFinfo << converged_message(reason) << CFendl;


  // unset matrix/vectors
  if ( (err=VecDestroy(&xp)) ||
       (err=VecDestroy(&bp)) ||
       (err=MatDestroy(&Ap)) )
    throw std::runtime_error(err_message((int) err,"unset matrix/vectors"));


  return *this;
}


petsc_seq& petsc_seq::copy(const petsc_seq& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  ksp = _other.ksp;
  pc  = _other.pc;
  opt = _other.opt;
  return *this;
}


const std::string petsc_seq::err_message(const int& err, const char* basemsg)
{
  std::ostringstream s;
  const char* text = "";
  PetscErrorMessage(err,&text,NULL);
  s << "petsc_seq " << basemsg << " error: " << err << ": " << text << '.';
  return s.str();
}


const std::string petsc_seq::converged_message(const KSPConvergedReason& rsn)
{
  std::ostringstream s;
  s << "petsc_seq solve: " << (rsn<0? "diverged":"converged") << ", reason: ";
  rsn==KSP_CONVERGED_ITERATING?       s << "iterating"   :
  rsn==KSP_CONVERGED_RTOL_NORMAL?     s << "rtol normal" :
  rsn==KSP_CONVERGED_ATOL_NORMAL?     s << "atol normal" :
  rsn==KSP_CONVERGED_RTOL?            s << "rtol" :
  rsn==KSP_CONVERGED_ATOL?            s << "atol" :
  rsn==KSP_CONVERGED_ITS?             s << "its"  :
  rsn==KSP_CONVERGED_CG_NEG_CURVE?    s << "cg neg curve"    :
  rsn==KSP_CONVERGED_CG_CONSTRAINED?  s << "cg constrained"  :
  rsn==KSP_CONVERGED_STEP_LENGTH?     s << "step length"     :
  rsn==KSP_CONVERGED_HAPPY_BREAKDOWN? s << "happy breakdown" :
  rsn==KSP_DIVERGED_NULL?             s << "null" :
  rsn==KSP_DIVERGED_ITS?              s << "its"  :
  rsn==KSP_DIVERGED_DTOL?             s << "dtol" :
  rsn==KSP_DIVERGED_BREAKDOWN?        s << "breakdown"      :
  rsn==KSP_DIVERGED_BREAKDOWN_BICG?   s << "breakdown bicg" :
  rsn==KSP_DIVERGED_NONSYMMETRIC?     s << "nonsymmetric"   :
  rsn==KSP_DIVERGED_INDEFINITE_PC?    s << "indefinite pc"  :
  rsn==KSP_DIVERGED_NANORINF?         s << "nanorinf"       :
  rsn==KSP_DIVERGED_INDEFINITE_MAT?   s << "indefinite mat" :
  rsn==KSP_CONVERGED_ITERATING?       s << "iterating"      :
                                      s << "(unknown)";
  s << '.';
  return s.str();
}


}  // namespace lss
}  // namespace cf3

