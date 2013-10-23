// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_GMRES_h
#define cf3_lss_GMRES_h


#include "LibLSS.hpp"
#include "linearsystem.h"
#include "detail/linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * implementation of a serial GMRES linear system solver (double p.)
 */
class GMRES : public
  linearsystem,
  detail::linearsystem< double,
    detail::sparse_matrix_csr< double, 1 >,
    detail::dense_matrix_v< double > >
{
  // utility definitions
  typedef detail::sparse_matrix_csr< double, 1 > matrix_t;
  typedef detail::dense_matrix_v< double > vector_t;
  typedef detail::linearsystem< double, matrix_t, vector_t > linearsystem_t;


  // framework interfacing
 public:
  static std::string type_name();
  GMRES(const std::string& name,
        const size_t& _size_i=size_t(),
        const size_t& _size_j=size_t(),
        const size_t& _size_k=1,
        const double& _value=double() );

  /// Linear system resizing (consistently)
  GMRES& resize(
      const size_t& _size_i,
      const size_t& _size_j,
      const size_t& _size_k=1,
      const double& _value=double()) { linearsystem_t::resize(_size_i,_size_j,_size_k,_value); return *this; }

  /// Linear system initialization from file(s)
  GMRES& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) { linearsystem_t::initialize(_Afname,_bfname,_xfname); return *this; }

  /// Linear system initialization from vectors of values (lists, in right context)
  GMRES& initialize(
      const std::vector< double >& vA,
      const std::vector< double >& vb=std::vector< double >(),
      const std::vector< double >& vx=std::vector< double >()) { linearsystem_t::initialize(vA,vb,vx); return *this; }

  /// Linear system solving
  GMRES& solve();


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


 public:
        matrix_t& A()       { return m_A; }
        vector_t& b()       { return m_b; }
        vector_t& x()       { return m_x; }
  const matrix_t& A() const { return m_A; }
  const vector_t& b() const { return m_b; }
  const vector_t& x() const { return m_x; }


 protected:
  matrix_t m_A;
  vector_t m_b;
  vector_t m_x;
  int c__1;

};


}  // namespace lss
}  // namespace cf3


#endif
