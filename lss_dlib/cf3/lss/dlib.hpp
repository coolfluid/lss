// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_dlib_hpp
#define cf3_lss_dlib_hpp


#include "Dlib.h"

#include "LibLSS_DLIB.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Interface to Dlib linear system solver.
 * @author Pedro Maciel
 */
class lss_API dlib : public linearsystem< double >
{
  // utility definitions
  typedef dlib::matrix< double, 0, 0 > matrix_t;


 public:
  // framework interfacing
  static std::string type_name() { return "dlib"; }


  /// Construction
  dlib(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1,
    const double& _value=double() )
    : linearsystem< double >(name) {
    linearsystem< double >::initialize(_size_i,_size_j,_size_k,_value);
  }

  /// Solve
  dlib& solve();


 protected:
  // linear system matrix interfacing

  const double& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        double& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  void A___initialize(const size_t& i, const size_t& j, const double& _value=double()) { m_A.set_size(i,j; m_A = _value}
  void A___initialize(const std::vector< double >& _vector) { /*m_A.initialize(_vector);*/ }
  void A___initialize(const std::string& _fname)            { /*m_A.initialize(_fname);*/  }
  void A___clear()                  { /*m_A.clear();*/    }
  void A___zerorow(const size_t& i) { dlib::set_rowm(m_A,i) = double(); }
  void A___sumrows(const size_t& i, const size_t& isrc) { m_A.sumrows(i,isrc); }

  void   A___print(std::ostream& o, const print_t& l=print_auto) const { /*m_A.print(o,l);*/ }
  size_t A___size(const size_t& d)  const { return 0;/*m_A.size(d);*/  }


 protected:
  // storage
  matrix_t m_A;

};


}  // namespace lss
}  // namespace cf3


#endif

