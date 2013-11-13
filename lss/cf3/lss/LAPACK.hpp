// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_LAPACK_h
#define cf3_lss_LAPACK_h


#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/// Prototypes for single and double precision (as per Intel MKL documentation)
extern "C"
{
  void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
  void sgesv_(int* n, int* nrhs, float*  a, int* lda, int* ipiv, float*  b, int* ldb, int* info);
}


/**
 * @brief example linear system solver, using LAPACK
 * (available in single and double precision, only works for square matrices)
 */
template< typename T >
class lss_API LAPACK : public
  linearsystem< T,
    detail::dense_matrix_v< T >,
    detail::dense_matrix_v< T > >
{
  // utility definitions
  typedef detail::dense_matrix_v< T > matrix_t;
  typedef detail::dense_matrix_v< T > vector_t;
  typedef linearsystem< T, matrix_t, vector_t > linearsystem_t;

 public:
  // framework interfacing
  static std::string type_name();

  /// Construction
  LAPACK(const std::string& name,
         const size_t& _size_i=size_t(),
         const size_t& _size_j=size_t(),
         const size_t& _size_k=1,
         const double& _value=T() ) : linearsystem_t(name) {
    linearsystem_t::initialize(_size_i,_size_j,_size_k,_value);
  }

  /// Solve
  LAPACK& solve() {
    int n    = static_cast< int >(linearsystem_t::size(0));
    int nrhs = static_cast< int >(linearsystem_t::size(2));
    int err  = 0;
    std::vector< int > ipiv(n);

    if (!m_A.size().is_square_size()) { err = -17; }
    else if (detail::type_is_equal< T, double >()) { x()=b(); dgesv_( &n, &nrhs, (double*) &m_A.a[0], &n, &ipiv[0], (double*) &m_x.a[0], &n, &err ); }
    else if (detail::type_is_equal< T, float  >()) { x()=b(); sgesv_( &n, &nrhs, (float*)  &m_A.a[0], &n, &ipiv[0], (float*)  &m_x.a[0], &n, &err ); }
    else { err = -42; }

    std::ostringstream msg;
    err==-17? msg << "LAPACK: system matrix must be square." :
    err==-42? msg << "LAPACK: precision not implemented." :
    err<0?    msg << "LAPACK: invalid " << err << "'th argument to dgesv_()/sgesv_()." :
    err>0?    msg << "LAPACK: triangular factor matrix U(" << (err-1) << ',' << (err-1) << ") is zero, so A is singular (not invertible)." :
              msg;
    if (err)
      throw std::runtime_error(msg.str());
    return *this;
  }


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

};


}  // namespace lss
}  // namespace cf3


#endif
