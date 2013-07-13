// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_GMRES_h
#define cf3_lss_GMRES_h


#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * implementation of a serial GMRES linear system solver (double p.)
 */
class GMRES :
  public linearsystem< double >
{
  typedef linearsystem< double > linearsystem_t;

 public:
  // framework interfacing
  static std::string type_name() { return "GMRES"; } // (mandatory!)
  GMRES(const std::string& name);

 public:
  // linear system solver addressing
  const double& A(const size_t i) const { return m_A(i); }
        double& A(const size_t i)       { return m_A(i); }
  const double& A(const size_t r, const size_t c) const { return m_A(r,c); }
        double& A(const size_t r, const size_t c)       { return m_A(r,c); }

 public:
  // linear system solver interfacing
  size_t size()               const { return m_A.size(); }
  size_t size(const size_t d) const { return m_A.size(d); }

  /*
   void initialize(const std::vector< std::vector< size_t > >& nz) { m_A.initialize(nz); }
   */

  GMRES& resize(size_t Nequations, size_t Nvariables, const double& v=double());
  GMRES& zerorow(const size_t r);
  GMRES& clear();
  GMRES& solve();

  void output_A(std::ostream& out) const { m_A.output(out); }

 private:
  // internal functions
  int daxpy(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
  int iluk( int *n, double *a, int *ja, int *ia, int *lfil, double*& aluold,
            int*& jluold, int *ju, int*& levsold, int *iwk, double *w, int *jw,
            int *ierr );
  double ddot(int *n, double *dx, int *incx, double *dy, int *incy);
  double dnrm2(int *n, double *dx, int *incx);
  void amux(int *n, double *x, double *y, double *a, int *ja, int *ia);
  void lusol(int *n, double *y, double *x, double *alu, int *jlu, int *ju);
  void pgmres( int *n, int *im, double *rhs, double *sol, double *vv,
               double *eps, int *maxits, int*iout, double *aa, int *ja, int *ia,
               double *alu, int *jlu, int *ju, int *ierr );
  void trigger_IA();
  void trigger_JA();

  // members
 private:
  sparse_matrix_csr< double, 1 > m_A;
  int c__1;
};


}  // namespace lss
}  // namespace cf3


#endif
