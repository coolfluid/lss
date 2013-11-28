// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf

#include "petscmat.h"
#include "petscerror.h"
#include "petscksp.h"

#include "common/Builder.hpp"
#include "PETSc_Seq.hpp"


namespace cf3 {
namespace lss {


common::ComponentBuilder< PETSc_Seq, common::Component, LibLSS_PETSC > Builder_PETSc_Seq;


PETSc_Seq::PETSc_Seq(const std::string& name,
    const size_t& _size_i,
    const size_t& _size_j,
    const size_t& _size_k)
  : linearsystem< double >(name)
{
  // call PetscInitialize (probably not needed in this form)
  int argc = 0;
  char **args = NULL;
  char help[] = "";
  PetscInitialize(&argc,&args,(char*)0,help);

  linearsystem< double >::initialize(_size_i,_size_j,_size_k);
}


PETSc_Seq::~PETSc_Seq()
{
  PetscErrorCode err;
  err = PetscFinalize();
}


PETSc_Seq& PETSc_Seq::solve()
{
  PetscErrorCode err = 0;
  PetscInt
    m = static_cast< PetscInt >(m_A.size(0)),
    n = static_cast< PetscInt >(m_A.size(1));
  if (!m_A.m_size.is_square_size())
    throw std::runtime_error("PETSc_Seq: system matrix must be square.");

  matrix_t::matrix_compressed_t& A = m_A.compress();
  Mat Ap;
  Vec bp;
  Vec xp;
  if ( (err=MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,m,n,&A.ia[0],&A.ja[0],&A.a[0],&Ap)) ||
       (err=VecCreateSeqWithArray(PETSC_COMM_SELF,1,n,&m_b.a[0],&bp)) ||
       (err=VecCreateSeqWithArray(PETSC_COMM_SELF,1,n,&m_x.a[0],&xp)) )
    throw std::runtime_error(err_message((int) err,"set linear system components"));

  /*
   * - create linear solver context
   * - set operators (here the system matrix serves as preconditioning matrix)
   * - set the default preconditioner for this program to be ASM
   * - set the overlap (use the default decomposition via PCASMSetOverlap)
   */
  KSP      ksp;          // linear solver context
  PC       pc;           // PC context
  PetscInt overlap = 1;  // width of subdomain overlap
  if ( (err=KSPCreate(PETSC_COMM_WORLD,&ksp)) ||
       (err=KSPSetOperators(ksp,Ap,Ap,DIFFERENT_NONZERO_PATTERN)) ||
       (err=KSPGetPC(ksp,&pc))        ||
       (err=KSPSetFromOptions(ksp))   ||
       (err=KSPSetType(ksp,KSPGMRES)) ||
       (err=PCSetType(pc,PCASM))      ||
       (err=PCASMSetOverlap(pc,overlap)) )
    throw std::runtime_error(err_message((int) err,"set solver/preconditioner"));

  // Solve
  if (err=KSPSolve(ksp,bp,xp))
    throw std::runtime_error(err_message((int) err,"solve"));

  // PETSc objects should be destroyed when no longer needed.
  if ( (err=KSPDestroy(&ksp)) ||
       (err=VecDestroy(&xp))  ||
       (err=VecDestroy(&bp))  ||
       (err=MatDestroy(&Ap)) )
    throw std::runtime_error(err_message((int) err,"deallocation"));

  return *this;
}


PETSc_Seq& PETSc_Seq::copy(const PETSc_Seq& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  return *this;
}


const std::string PETSc_Seq::err_message(const int& err, const char* basemsg)
{
  std::ostringstream s;
  const char* text = "";
  PetscErrorMessage(err,&text,NULL);
  s << "PETSc_Seq " << basemsg << " error: " << err << ": " << text << '.';
  return s.str();
}


}  // namespace lss
}  // namespace cf3

