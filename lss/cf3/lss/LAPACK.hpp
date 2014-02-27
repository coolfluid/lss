// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_LAPACK_hpp
#define cf3_lss_LAPACK_hpp


#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/// Prototypes for single/double precisions, real/complex types
/// (as per Intel MKL documentation)
extern "C"
{
  void dgesv_(int* n, int* nrhs, double*  a, int* lda, int* ipiv, double*  b, int* ldb, int* info);
  void zgesv_(int* n, int* nrhs, zdouble* a, int* lda, int* ipiv, zdouble* b, int* ldb, int* info);
  void sgesv_(int* n, int* nrhs, float*   a, int* lda, int* ipiv, float*   b, int* ldb, int* info);
  void cgesv_(int* n, int* nrhs, zfloat*  a, int* lda, int* ipiv, zfloat*  b, int* ldb, int* info);
}


/**
 * @brief example linear system solver, using LAPACK (only for square matrices)
 */
template< typename T >
class lss_API LAPACK : public linearsystem< T >
{
  // utility definitions
  typedef dense_matrix_v< T > matrix_t;

 public:
  // framework interfacing
  static std::string type_name();

  /// Construction
  LAPACK(const std::string& name,
         const size_t& _size_i=size_t(),
         const size_t& _size_j=size_t(),
         const size_t& _size_k=1 ) : linearsystem< T >(name) {
    linearsystem< T >::initialize(_size_i,_size_j,_size_k);
  }

  /// Linear system solving
  LAPACK& solve() {
    int n    = static_cast< int >(this->size(0));
    int nrhs = static_cast< int >(this->size(2));
    int err  = 0;
    std::vector< int > ipiv(n);

    if (!m_A.m_size.is_square_size()) { err = -17; }
    else if (type_is_equal< T, double  >()) { this->m_x=this->m_b; dgesv_( &n, &nrhs, (double*)  &m_A.a[0], &n, &ipiv[0], (double*)  &this->m_x.a[0], &n, &err ); }
    else if (type_is_equal< T, zdouble >()) { this->m_x=this->m_b; zgesv_( &n, &nrhs, (zdouble*) &m_A.a[0], &n, &ipiv[0], (zdouble*) &this->m_x.a[0], &n, &err ); }
    else if (type_is_equal< T, float   >()) { this->m_x=this->m_b; sgesv_( &n, &nrhs, (float*)   &m_A.a[0], &n, &ipiv[0], (float*)   &this->m_x.a[0], &n, &err ); }
    else if (type_is_equal< T, zfloat  >()) { this->m_x=this->m_b; cgesv_( &n, &nrhs, (zfloat*)  &m_A.a[0], &n, &ipiv[0], (zfloat*)  &this->m_x.a[0], &n, &err ); }
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

  /// Linear system copy
  LAPACK& copy(const LAPACK& _other) {
    linearsystem< T >::copy(_other);
    m_A = _other.m_A;
    return *this;
  }


 protected:
  // linear system matrix interfacing

  /// matrix indexing
  const T& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        T& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  /// matrix modifiers
  void A___initialize(const size_t& i, const size_t& j, const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >()) { m_A.initialize(i,j); }
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

};


}  // namespace lss
}  // namespace cf3


#endif
