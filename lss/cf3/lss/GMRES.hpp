// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_GMRES_hpp
#define cf3_lss_GMRES_hpp


#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * implementation of a serial GMRES linear system solver (double p.)
 */
class lss_API GMRES : public linearsystem< double >
{
  // utility definitions
  typedef sparse_matrix< double, sort_by_row, 1 > matrix_t;

 public:
  // framework interfacing
  static std::string type_name();

  /// Construction
  GMRES(const std::string& name,
        const size_t& _size_i=size_t(),
        const size_t& _size_j=size_t(),
        const size_t& _size_k=1 ) : linearsystem< double >(name), c__1(1) {
    linearsystem< double >::initialize(_size_i,_size_j,_size_k);
  }

  /// Linear system solving: x = A^-1 b
  GMRES& solve();

  /// Linear system forward multiplication: b = alpha A x + beta b
  GMRES& multi(const double& _alpha=1., const double& _beta=0.);

  /// Linear system copy
  GMRES& copy(const GMRES& _other);

  /// Linear system swap
  GMRES& swap(GMRES& _other);


 private:
  // internal functions
  int iluk(int *n, double *a, int *ja, int *ia, int *lfil, double*& aluold, int*& jluold, int *ju, int*& levsold, int *iwk, double *w, int *jw, int *ierr);
  int daxpy(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
  double ddot(int *n, double *dx, int *incx, double *dy, int *incy);
  double dnrm2(int *n, double *dx, int *incx);
  void amux(int *n, double *x, double *y, double *a, int *ja, int *ia);
  void lusol(int *n, double *y, double *x, double *alu, int *jlu, int *ju);
  void pgmres(int *n, int *im, double *rhs, double *sol, double *vv, double *eps, int *maxits, int*iout, double *aa, int *ja, int *ia, double *alu, int *jlu, int *ju, int *ierr);


 protected:
  // linear system matrix interfacing

  /// matrix indexing
  const double& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        double& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  /// matrix modifiers
  void A___initialize(const size_t& i, const size_t& j, const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >()) { m_A.initialize(i,j,_nnz); }
  void A___initialize(const std::vector< double >& _vector) { m_A.initialize(_vector); }
  void A___initialize(const std::string& _fname)            { m_A.initialize(_fname);  }
  void A___assign(const double& _value)                 { m_A = _value;   }
  void A___clear()                                      { m_A.clear();    }
  void A___zerorow(const size_t& i)                     { m_A.zerorow(i); }
  void A___sumrows(const size_t& i, const size_t& isrc) { m_A.sumrows(i,isrc); }

  /// matrix inspecting
  void   A___print(std::ostream& o, const print_t& l=print_auto) const { m_A.print(o,l); }
  size_t A___size(const size_t& d)  const { return m_A.size(d);  }


 protected:
  // storage
  matrix_t m_A;
  int c__1;

};


}  // namespace lss
}  // namespace cf3


#endif
