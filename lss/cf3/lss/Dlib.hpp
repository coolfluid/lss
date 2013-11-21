// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_Dlib_hpp
#define cf3_lss_Dlib_hpp


#include "dlib/matrix.h"

#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Example interface to Dlib linear system solver (configurable precision)
 * @author Pedro Maciel
 */
template< typename T >
class lss_API Dlib : public linearsystem< T >
{
  // utility definitions
  typedef dlib::matrix< T > matrix_t;


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
    dlib::matrix< T > rhs(linearsystem< T >::size(1),linearsystem< T >::size(2));
    for (size_t i=0; i<linearsystem< T >::size(1); ++i)
      for (size_t j=0; j<linearsystem< T >::size(2); ++j)
        rhs(i,j) = linearsystem< T >::b(i,j);

    dlib::matrix< T > sol = dlib::inv(m_A)*rhs;

    for (size_t i=0; i<linearsystem< T >::size(0); ++i)
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
  void A___initialize(const size_t& i, const size_t& j, const double& _value=double()) {
    m_A.set_size(i,j);
    m_A = _value;
  }

  void A___initialize(const std::vector< double >& _vector) { /*m_A.initialize(_vector);*/ }
  void A___initialize(const std::string& _fname)            { /*m_A.initialize(_fname);*/  }
  void A___clear()                  { /*m_A.clear();*/    }
  void A___zerorow(const size_t& i) { dlib::set_rowm(m_A,i) = T(); }
  void A___sumrows(const size_t& i, const size_t& isrc) { /*m_A.sumrows(i,isrc);*/ }

  /// matrix inspecting
  void   A___print(std::ostream& o, const print_t& l=print_auto) const {
    o << m_A;
  }
  size_t A___size(const size_t& d)  const { return 0;/*m_A.size(d);*/  }


 protected:
  // storage
  matrix_t m_A;

};


}  // namespace lss
}  // namespace cf3


#endif
