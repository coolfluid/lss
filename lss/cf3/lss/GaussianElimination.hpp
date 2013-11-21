// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_GaussianElimination_hpp
#define cf3_lss_GaussianElimination_hpp


#include <cmath>
#include "boost/progress.hpp"

#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Example linear system solver, using Gaussian elimination (configurable precision)
 * @author Pedro Maciel
 */
template< typename T >
class lss_API GaussianElimination : public linearsystem< T >
{
  // utility definitions
  typedef dense_matrix_v< T, sort_by_row > matrix_t;

  // framework interfacing
 public:
  static std::string type_name();

  /// Construction
  GaussianElimination(const std::string& name,
                      const size_t& _size_i=size_t(),
                      const size_t& _size_j=size_t(),
                      const size_t& _size_k=1,
                      const T& _value=T() ) : linearsystem< T >(name) {
    linearsystem< T >::initialize(_size_i,_size_j,_size_k,_value);
  }

  /// Linear system solving
  GaussianElimination& solve() {
    const size_t
      N(linearsystem< T >::size(0)),
      K(linearsystem< T >::size(2));

    this->m_x = this->m_b;

    boost::progress_display pbar(N-1);
    for (size_t m=0; m<N-1; ++m, ++pbar) {

      // put row with highest diagonal element on top
      T C = m_A(m,m);
      for (size_t n=m+1; n<N; ++n) {
        if (std::abs(m_A(n,m)) > std::abs(C)) {
          for (size_t p=m; p<N; ++p) {
            C = m_A(m,p);
            m_A(m,p) = m_A(n,p);
            m_A(n,p) = C;
          }
          for (size_t k=0; k<K; ++k)
            std::swap( this->x(m,k), this->x(n,k) );
          C = m_A(m,m);
        }
      }

      // check if diagonal element is (close to) zero
      if (std::abs(C)<1.e-32) {
        std::ostringstream msg;
        msg << "error: matrix is singular (line:" << m << ",C:" << std::abs(C) << ").";
        throw std::runtime_error(msg.str());
      }

      // normalize row m
      for (size_t n=m+1; n<N; ++n)
        m_A(m,n) /= C;
      for (size_t k=0; k<K; ++k)
        this->x(m,k) /= C;

      // subtract row m from subsequent rows
      for (size_t n=m+1; n<N; ++n) {
        C = m_A(n,m);
        for (size_t p=m+1; p<N; ++p)
          m_A(n,p) -= C*m_A(m,p);
        for (size_t k=0; k<K; ++k)
          this->x(n,k) -= C * this->x(m,k);
      }
    }

    // solve by back substitution
    for (size_t k=0; k<K; ++k)
      this->x(N-1,k) /= m_A(N-1,N-1);
    for (size_t p=0; p<N-1; ++p) {
      size_t m = N-p-2;
      for (size_t n=m+1; n<N; ++n)
        for (size_t k=0; k<K; ++k)
          this->x(m,k) -= m_A(m,n) * this->x(n,k);
    }

    return *this;
  }

  /// Linear system copy
  GaussianElimination& copy(const GaussianElimination& _other) {
    linearsystem< T >::copy(_other);
    m_A = _other.m_A;
  }


 protected:
  // linear system matrix interfacing

  /// matrix indexing
  const T& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        T& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  /// matrix modifiers
  void A___initialize(const size_t& i, const size_t& j, const double& _value=double()) { m_A.initialize(i,j,_value); }
  void A___initialize(const std::vector< double >& _vector) { m_A.initialize(_vector); }
  void A___initialize(const std::string& _fname)            { m_A.initialize(_fname);  }
  void A___clear()                  { m_A.clear();    }
  void A___zerorow(const size_t& i) { m_A.zerorow(i); }
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
