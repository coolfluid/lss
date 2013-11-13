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
class lss_API GMRES : public
  linearsystem< double,
    detail::sparse_matrix< double, 1>,
    detail::dense_matrix_v< double, detail::sort_by_row > >
{
  // utility definitions
  typedef detail::sparse_matrix< double, 1> matrix_t;
  typedef detail::dense_matrix_v< double, detail::sort_by_row > vector_t;
  typedef linearsystem< double, matrix_t, vector_t > linearsystem_t;


 public:
  // framework interfacing
  static std::string type_name();

  /// Construction
  GMRES(const std::string& name,
        const size_t& _size_i=size_t(),
        const size_t& _size_j=size_t(),
        const size_t& _size_k=1,
        const double& _value=double() ) : linearsystem_t(name), c__1(1) { linearsystem_t::initialize(_size_i,_size_j,_size_k,_value); }

  /// Solve
  GMRES& solve();


 private:
  // internal functions
  int daxpy(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
  int iluk( int *n, double *a, int *ja, int *ia, int *lfil,
            std::vector< double >& aluold,
            std::vector< int >& jluold,
            int *ju, int* levsold, int *iwk, double *w, int *jw,
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
