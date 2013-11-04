// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cmath>

#include "common/Builder.hpp"

#include "LibLSS_pardiso.hpp"
#include "pardiso.hpp"


// Prototypes for external function calls
extern "C" {
  int pardisoinit_(void *, int *, int *, int *, double *, int *);
  int pardiso_(void *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *, double *);
}


namespace cf3 {
namespace lss {


std::string pardiso::type_name() { return "pardiso"; }
common::ComponentBuilder< pardiso, common::Component, LibLSS_pardiso > Builder_pardiso;


pardiso& pardiso::solve()
{
  call_pardisoinit();  // setup
  call_pardiso(11,1);  // 11: reordering and symbolic factorization
  call_pardiso(22,0);  // 22: numerical factorization and
  call_pardiso(33,0);  // 33: back substitution and iterative refinement
  call_pardiso(-1,0);  // -1: termination and release of memory
  return *this;
}


int pardiso::call_pardisoinit() {
  int err = 0;
  pardisoinit_(pt,&mtype,&iparm[31],iparm,dparm,&err);

  std::ostringstream msg(err? "pardisoinit: ":"");
  err==-10? msg << "no license file pardiso.lic found." :
  err==-11? msg << "license is expired." :
  err==-12? msg << "wrong username or hostname." :
  err?      msg << "unknown error." :
            msg;

  if (err)
    throw std::runtime_error(msg.str());
  return err;
}


int pardiso::call_pardiso(int phase, int msglvl)
{
  int err = 0;
  pardiso_( pt, &maxfct, &mnum, &mtype, &phase, &m_A.idx.nnu,
            &m_A.a[0], &m_A.idx.ia[0], &m_A.idx.ja[0], NULL, &nrhs, iparm,
            &msglvl, &m_b.a[0], &m_x.a[0], &err, dparm );

  std::ostringstream msg;
  err?       msg << "pardiso (phase/error: " << phase << '/' << err << "): " : msg;
  err==  -1? msg << "input inconsistent." :
  err==  -2? msg << "not enough memory."  :
  err==  -3? msg << "reordering problem." :
  err==  -4? msg << "zero pivot, numerical fact. or iterative refinement problem." :
  err==  -5? msg << "unclassified (internal) error."     :
  err==  -6? msg << "preordering failed (matrix types 11, 13 only)." :
  err==  -7? msg << "diagonal matrix problem."           :
  err==  -8? msg << "32-bit integer overflow problem."   :
  err== -10? msg << "no license file pardiso.lic found." :
  err== -11? msg << "license is expired."                :
  err== -12? msg << "wrong username or hostname."        :
  err==-100? msg << "reached maximum number of Krylov-subspace iteration in iterative solver."     :
  err==-101? msg << "no sufficient convergence in Krylov-subspace iteration within 25 iterations." :
  err==-102? msg << "error in Krylov-subspace iteration."      :
  err==-103? msg << "break-down in Krylov-subspace iteration." :
             msg << "unknown error.";

  if (err)
    throw std::runtime_error(msg.str());
  return err;
}


}  // namespace lss
}  // namespace cf3

