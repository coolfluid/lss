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
 * example linear system solver, using LAPACK
 * (available in single and double precision, only works for square matrices)
 */
template<
    typename T,
    typename INDEX=index_hierarchy_t< index_hierarchy_t_end > >
struct lss_API LAPACK :
  public linearsystem<
    T,
    INDEX,
    detail::dense_matrix_v< T, detail::column_oriented >,
    detail::dense_matrix_v< T, detail::column_oriented > >
{
  // utility definitions
  typedef detail::dense_matrix_v< T, detail::column_oriented > matrix_t;
  typedef detail::dense_matrix_v< T, detail::column_oriented > vector_t;
  typedef linearsystem< T, INDEX, matrix_t, vector_t > linearsystem_base_t;

  // framework interfacing
  static std::string type_name();
  LAPACK(const std::string& name) : linearsystem_base_t(name) {}
  LAPACK(const size_t& _size_i=size_t(),
         const size_t& _size_j=size_t(),
         const size_t& _size_k=1,
         const T& _value=T() ) :
    linearsystem_base_t(_size_i,_size_j,_size_k,_value) {}

  LAPACK& solve() {
    int n    = static_cast< int >(this->size(0));
    int nrhs = static_cast< int >(this->size(2));
    int err  = 0;
    std::vector< int > ipiv(n);

    if (!m_A.size().is_square_size()) { err = -17; }
    else if (typeid(T)==typeid(double)) { dgesv_( &n, &nrhs, (double*) &this->m_A.a[0], &n, &ipiv[0], (double*) &this->m_b.a[0], &n, &err ); }
    else if (typeid(T)==typeid(float))  { sgesv_( &n, &nrhs, (float*)  &this->m_A.a[0], &n, &ipiv[0], (float*)  &this->m_b.a[0], &n, &err ); }
    else { err = -42; }

    std::ostringstream msg;
    err==-17? msg << "LAPACK: system matrix size must be square" :
    err==-42? msg << "LAPACK: precision not implemented: " << typeid(T).name() :
    err<0?    msg << "LAPACK: invalid " << err << "'th argument to dgesv_()/sgesv_()" :
    err>0?    msg << "LAPACK: triangular factor matrix U(" << err << ',' << err << ") is zero, so A is singular (not invertible)" :
              msg;
    if (err)
      throw std::runtime_error(msg.str());

    // solution is in b, so swap with x (with A square, size b = size x)
    this->m_b.swap(this->m_x);
    return *this;
  }

  size_t size(const size_t& d) {
    //FIXME return the right size
    return 42;
  }


  // storage
 private:
  matrix_t m_A;
  vector_t m_b;
  vector_t m_x;


 public: //FIXME: permissions

        matrix_t& A()       { return m_A; }
        vector_t& b()       { return m_b; }
        vector_t& x()       { return m_x; }
  const matrix_t& A() const { return m_A; }
  const vector_t& b() const { return m_b; }
  const vector_t& x() const { return m_x; }

  /*

        T& A(const size_t& i, const size_t& j)         { return m_A(i,j); }
        T& b(const size_t& i, const size_t& j=0)       { return m_b(i,j); }
        T& x(const size_t& i, const size_t& j=0)       { return m_x(i,j); }
  const T& A(const size_t& i, const size_t& j)   const { return m_A(i,j); }
  const T& b(const size_t& i, const size_t& j=0) const { return m_b(i,j); }
  const T& x(const size_t& i, const size_t& j=0) const { return m_x(i,j); }

  size_t A_size(const size_t& _d)   const { return m_A.size(_d); }
  size_t b_size(const size_t& _d=0) const { return m_b.size(_d); }
  size_t x_size(const size_t& _d=0) const { return m_x.size(_d); }

  LAPACK& A_resize(const size_t& _size_i, const size_t& _size_j,   const T& _value=T()) { m_A.resize(_size_i,_size_j,_value); return *this; }
  LAPACK& b_resize(const size_t& _size_i, const size_t& _size_j=1, const T& _value=T()) { m_b.resize(_size_i,_size_j,_value); return *this; }
  LAPACK& x_resize(const size_t& _size_i, const size_t& _size_j=1, const T& _value=T()) { m_x.resize(_size_i,_size_j,_value); return *this; }

  LAPACK& A_initialize(const std::vector< T >& _vector) { m_A.initialize(_vector); return *this; }
  LAPACK& b_initialize(const std::vector< T >& _vector) { m_b.initialize(_vector); return *this; }
  LAPACK& x_initialize(const std::vector< T >& _vector) { m_x.initialize(_vector); return *this; }
  LAPACK& A_initialize(const std::string& _fname) { m_A.initialize(_fname); return *this; }
  LAPACK& b_initialize(const std::string& _fname) { m_b.initialize(_fname); return *this; }
  LAPACK& x_initialize(const std::string& _fname) { m_x.initialize(_fname); return *this; }

  void A_clear() { m_A.clear(); }
  void b_clear() { m_b.clear(); }
  void x_clear() { m_x.clear(); }

  void A_print(std::ostream& o) { m_A.print(o); }
  void b_print(std::ostream& o) { m_b.print(o); }
  void x_print(std::ostream& o) { m_x.print(o); }

  LAPACK& operator=(const LAPACK& _other) {
    m_A = _other.m_A;
    return *this;
  }
  */

};


}  // namespace lss
}  // namespace cf3


#endif
