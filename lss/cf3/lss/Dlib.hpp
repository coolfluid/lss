// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_Dlib_hpp
#define cf3_lss_Dlib_hpp


#include "dlib/assert.h"
#include "dlib/matrix.h"

#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


namespace detail {


/**
 * @brief Matrix wrapper for Dlib matrices
 * @author Pedro Maciel
 */
template< typename T >
struct dense_matrix_dlib :
    matrix< T, dense_matrix_dlib< T > >
{
  typedef matrix< T, dense_matrix_dlib< T > > matrix_base_t;
  using matrix_base_t::size;

  // initializations

  dense_matrix_dlib& initialize(const size_t& i, const size_t& j, const double& _value=double()) {
    if (idx_t(i,j).is_valid_size()) {
      matrix_base_t::m_size = idx_t(i,j);
      a.set_size(i,j);
      a = static_cast< T >(_value);
    }
    return *this;
  }

  dense_matrix_dlib& initialize(const std::vector< double >& _vector) {
    matrix_base_t::initialize(_vector);
    return *this;
  }

  dense_matrix_dlib& initialize(const std::string& _fname) {
    matrix_base_t::initialize(_fname);
    return *this;
  }

  // assignments

  dense_matrix_dlib& operator=(const dense_matrix_dlib& _other) {
    a = _other.a;
    matrix_base_t::m_size = _other.m_size;
    matrix_base_t::m_zero = _other.m_zero;
    return *this;
  }

  dense_matrix_dlib& operator=(const double& _value) {
    return initialize(size(0),size(1),_value);
  }

  void swap(dense_matrix_dlib& other) {
    a.swap(other.a);
    matrix_base_t::swap(other);
  }

  // clearing

  dense_matrix_dlib& clear() {
    a.set_size(0,0);
    matrix_base_t::m_size.clear();
    return *this;
  }

  dense_matrix_dlib& zerorow(const size_t& i) {
    if (i>=size(0))
      throw std::runtime_error("dense_matrix_dlib: row index outside bounds.");
    dlib::set_rowm(a,i) = T();
    return *this;
  }

  dense_matrix_dlib& sumrows(const size_t& i, const size_t& isrc) {
    if (std::max(i,isrc)>=size(0))
      throw std::runtime_error("dense_matrix_dlib: row index(es) outside bounds.");
    for (long j=0; j<a.nc(); ++j)
      operator()(i,j) += operator()(isrc,j);
    return *this;
  }

  // indexing
  const T& operator()(const size_t& i, const size_t& j=0) const { return a(i,j); }
        T& operator()(const size_t& i, const size_t& j=0)       { return a(i,j); }

  // storage
  dlib::matrix< T > a;
};


}  // namespace detail



/**
 * @brief Example interface to Dlib linear system solver (configurable precision)
 * @author Pedro Maciel
 */
template< typename T >
class lss_API Dlib : public linearsystem< T >
{
  // utility definitions
  typedef detail::dense_matrix_dlib< T > matrix_t;


  // framework interfacing
 public:
  static std::string type_name();

  /// Construction
  Dlib(const std::string& name,
       const size_t& _size_i=size_t(),
       const size_t& _size_j=size_t(),
       const size_t& _size_k=1,
       const T& _value=T() ) : linearsystem< T >(name) {
    linearsystem< T >::initialize(_size_i,_size_j,_size_k,_value);
  }

  /// Linear system solving
  Dlib& solve() {

    // extremely innefficient but elegant, like a pig with lipstick...
    dlib::matrix< T > rhs(linearsystem< T >::size(0),linearsystem< T >::size(2));
    for (size_t i=0; i<linearsystem< T >::size(1); ++i)
      for (size_t j=0; j<linearsystem< T >::size(2); ++j)
        rhs(i,j) = linearsystem< T >::b(i,j);

    dlib::matrix< T > sol = dlib::inv(m_A.a) * rhs;

    for (size_t i=0; i<linearsystem< T >::size(1); ++i)
      for (size_t j=0; j<linearsystem< T >::size(2); ++j)
        linearsystem< T >::x(i,j) = sol(i,j);
    return *this;
  }

  /// Linear system copy
  Dlib& copy(const Dlib& _other) {
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
  void A___assign(const double& _value)                 { m_A.operator=(_value); }
  void A___clear()                                      { m_A.clear();           }
  void A___zerorow(const size_t& i)                     { m_A.zerorow(i);        }
  void A___sumrows(const size_t& i, const size_t& isrc) { m_A.sumrows(i,isrc);   }

  /// matrix inspecting
  void   A___print(std::ostream& o, const print_t& l=print_auto) const { m_A.print(o,l); }
  size_t A___size(const size_t& d)  const { return m_A.size(d); }


 protected:
  // storage
  matrix_t m_A;

};


}  // namespace lss
}  // namespace cf3


#endif
