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
  public linearsystem< T, INDEX >
{
  // utility definitions
  typedef lss_matrix::dense_matrix_v< T, lss_matrix::column_oriented > matrix_t;
  typedef lss_matrix::dense_matrix_v< T, lss_matrix::column_oriented > vector_t;
  typedef linearsystem< T, INDEX > linearsystem_base_t;

  // framework interfacing
  static std::string type_name() { return "LAPACK"; }
  LAPACK(const std::string& name) :
    linearsystem_base_t(name) {}


  LAPACK(const size_t& _size_i=size_t(),
             const size_t& _size_j=size_t(),
             const size_t& _size_k=1,
             const T& _value=T() ) :
    linearsystem_base_t(_size_i,_size_j,_size_k,_value) {}

  bool solve() {
    int n    = static_cast< int >(this->size(0));
    int nrhs = static_cast< int >(this->size(2));
    int err  = 0;
    std::vector< int > ipiv(n);

    if (!this->size().is_square_size()) { err = -17; }
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
    return !err;
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
};


}  // namespace lss
}  // namespace cf3


#endif
