// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cmath>

#include "common/Builder.hpp"

#include "LibLSS_PARDISO.hpp"
#include "pardiso.hpp"


// Prototypes for external function calls
extern "C" {
  int pardisoinit_(void *, int *, int *, int *, double *, int *);
  int pardiso_(void *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *, double *);  
  void pardiso_chkmatrix_(int *mtype, int *n, double *a, int *ia, int *ja, int *error);
  void pardiso_chkvec_(int *n, int *nrhs, double *b, int *error);
  void pardiso_printstats_(int *mtype, int *n, double *a, int *ia, int *ja, int *nrhs, double *b, int *error);
}


namespace cf3 {
namespace lss {


common::ComponentBuilder< pardiso, common::Component, LibLSS_PARDISO > Builder_pardiso;


pardiso& pardiso::solve()
{

#if 0
const double arr_va[] = {1.,	-1., -3., -2., 5., 4., 6., 4., -4.,	2.,	7.,	8.,	-5};
const int    arr_ja[] = {1,   2,   4,   1,  2,  3,  4,  5,   1,  3,  4,  2,   5};
const int    arr_ia[] = {1, 4, 6, 9, 12, 14};
const double arr_vb[] = {0., 0., 0., 0., 0.};
const double arr_vx[] = {0., 0., 0., 0., 0.};
#else
const double arr_va[] = {-1., -1., -1.};
const int    arr_ja[] = { 1,   2,   3 };
const int    arr_ia[] = {1, 2, 3, 4};
const double arr_vb[] = {2., 2., 2.};
const double arr_vx[] = {0., 0., 0.};
#endif
va.assign(arr_va, arr_va + sizeof(arr_va)/sizeof(arr_va[0]));
ja.assign(arr_ja, arr_ja + sizeof(arr_ja)/sizeof(arr_ja[0]));
ia.assign(arr_ia, arr_ia + sizeof(arr_ia)/sizeof(arr_ia[0]));
vb.assign(arr_vb, arr_vb + sizeof(arr_vb)/sizeof(arr_vb[0]));
vx.assign(arr_vx, arr_vx + sizeof(arr_vx)/sizeof(arr_vx[0]));

std::cout << std::endl << "a:  (size " << va.size() << ") " << std::endl; std::copy(va.begin(),va.end(), std::ostream_iterator< double >(std::cout, " ")); std::cout << std::endl;
std::cout << std::endl << "ia: (size " << ia.size() << ") " << std::endl; std::copy(ia.begin(),ia.end(), std::ostream_iterator< int    >(std::cout, " ")); std::cout << std::endl;
std::cout << std::endl << "ja: (size " << ja.size() << ") " << std::endl; std::copy(ja.begin(),ja.end(), std::ostream_iterator< int    >(std::cout, " ")); std::cout << std::endl;




int nnu = ia.size()-1,
    err,
    phase,
    msglvl;

std::cout << "pardiso_chkmatrix_... " << std::endl;
err = 0;
pardiso_chkmatrix_(&mtype, &nnu, &va[0], &ia[0], &ja[0], &err);
if (err) { std::cout << "error: " << err << std::endl; }
std::cout << "pardiso_chkmatrix_. " << std::endl;

std::cout << "pardiso_chkvec_... " << std::endl;
err = 0;
pardiso_chkvec_(&nnu, &nrhs, &vb[0],&err);
if (err) { std::cout << "error: " << err << std::endl; }
std::cout << "pardiso_chkvec_. " << std::endl;

std::cout << "pardiso_printstats_... " << std::endl;
err = 0;
pardiso_printstats_(&mtype, &nnu, &va[0], &ia[0], &ja[0], &nrhs, &vb[0], &err);
if (err) { std::cout << "error: " << err << std::endl; }
std::cout << "pardiso_printstats_. " << std::endl;

std::cout << "pardisoinit_... " << std::endl;
err = 0;
pardisoinit_(pt,&mtype,&iparm[31],iparm,dparm,&err);
if (err) { std::cout << "error: " << err << std::endl; }
std::cout << "pardisoinit_. " << std::endl;

std::cout << "pardiso_... " << std::endl; call_pardiso(11,1); std::cout << "pardiso_. " << std::endl;
std::cout << "pardiso_... " << std::endl; call_pardiso(-1,1); std::cout << "pardiso_. " << std::endl;
std::cout << "pardiso_... " << std::endl; call_pardiso(12,1); std::cout << "pardiso_. " << std::endl;
std::cout << "pardiso_... " << std::endl; call_pardiso(-1,1); std::cout << "pardiso_. " << std::endl;
//std::cout << "pardiso_... " << std::endl; call_pardiso(22,1); std::cout << "pardiso_. " << std::endl;
//std::cout << "pardiso_... " << std::endl; call_pardiso(33,1); std::cout << "pardiso_. " << std::endl;



#if 0
  call_pardisoinit();  // setup
  call_pardiso(11,1);  // 11: reordering and symbolic factorization
  call_pardiso(22,0);  // 22: numerical factorization and
  call_pardiso(33,0);  // 33: back substitution and iterative refinement
  call_pardiso(-1,0);  // -1: termination and release of memory
#endif
  return *this;
}


int pardiso::call_pardisoinit() {
  int err = 0;
  pardisoinit_(pt,&mtype,&iparm[31],iparm,dparm,&err);

  std::ostringstream msg("pardisoinit: ");
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
  /*
  pardiso_( pt, &maxfct, &mnum, &mtype, &phase, &m_A.idx.nnu,
            &m_A.a[0], &m_A.idx.ia[0], &m_A.idx.ja[0], NULL, &nrhs, iparm,
            &msglvl, &m_b.a[0], &m_x.a[0], &err, dparm );
            */
  int nnu = ia.size()-1;
  pardiso_( pt, &maxfct, &mnum, &mtype, &phase, &nnu,
            &va[0], &ia[0], &ja[0], NULL, &nrhs, iparm,
            &msglvl, &vb[0], &vx[0], &err, dparm );

  std::ostringstream msg("pardiso: ");
  err?       msg << "(phase " << phase << " error " << err << ") " : msg;
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
  err?       msg << "unknown error." :
             msg;

  if (err)
    throw std::runtime_error(msg.str());
  return err;
}


}  // namespace lss
}  // namespace cf3

