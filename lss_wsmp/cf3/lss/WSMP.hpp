// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_WSMP_hpp
#define cf3_lss_WSMP_hpp


#include "LibLSS_WSMP.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Interface to WSMP linear system solver (serial version).
 * @author Pedro Maciel
 */
class lss_API WSMP : public linearsystem< double >
{
  // utility definitions
  typedef sparse_matrix< double > matrix_t;


 public:
  // framework interfacing
  static std::string type_name() { return "wsmp"; }


  /// Construction
  WSMP(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1,
    const double& _value=double() );

  /// Solve
  WSMP& solve();


 private:
  // internal functions
  int call_wsmp(int _task);


 protected:
  // linear system matrix interfacing

  const double& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        double& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  void A___initialize(const size_t& i, const size_t& j, const double& _value=double()) { m_A.initialize(i,j,_value); }
  void A___initialize(const std::vector< double >& _vector) { m_A.initialize(_vector); }
  void A___initialize(const std::string& _fname)            { m_A.initialize(_fname);  }
  void A___initialize(const index_t& _index)                { m_A.initialize(_index);  }
  void A___clear()                    { m_A.clear();    }
  void A___zerorow(const size_t& i)   { m_A.zerorow(i); }
  void A___print_level(const int& _l) { m_A.m_print = print_level(_l); }

  std::ostream& A___print(std::ostream& o) const { return m_A.print(o); }
  size_t        A___size(const size_t& d)  const { return m_A.size(d);  }


 protected:
  // storage
  matrix_t m_A;
  double dparm[64];
  int    iparm[64],
         task;

};


}  // namespace lss
}  // namespace cf3


#endif

