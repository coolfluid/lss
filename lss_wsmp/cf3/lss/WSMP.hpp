// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
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
  typedef sparse_matrix< double, sort_by_row, 0 > matrix_t;


 public:
  // framework interfacing
  static std::string type_name() { return "wsmp"; }


  /// Construction
  WSMP(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1 );

  /// Linear system solving
  WSMP& solve();

  /// Linear system forward multiplication
  WSMP& multi(const double& _alpha=1., const double& _beta=0.);

  /// Linear system copy
  WSMP& copy(const WSMP& _other);


 private:
  // internal functions

  /// Verbose error message
  static const std::string err_message(const int& err);

  /// Library call
  int call_wsmp(int _phase);


 protected:
  // linear system matrix interfacing

  const double& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        double& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  void A___initialize(const size_t& i, const size_t& j, const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >()) { m_A.initialize(i,j,_nnz); }
  void A___initialize(const std::vector< double >& _vector) { m_A.initialize(_vector); }
  void A___initialize(const std::string& _fname)            { m_A.initialize(_fname);  }
  void A___assign(const double& _value) { m_A = _value;   }
  void A___clear()                      { m_A.clear();    }
  void A___zerorow(const size_t& i)     { m_A.zerorow(i); }
  void A___sumrows(const size_t& i, const size_t& isrc) { m_A.sumrows(i,isrc); }

  void   A___print(std::ostream& o, const print_t& l=print_auto) const { m_A.print(o,l); }
  size_t A___size(const size_t& d)  const { return m_A.size(d);  }


 protected:
  // storage
  matrix_t m_A;
  double dparm[64];
  int    iparm[64];

};


}  // namespace lss
}  // namespace cf3


#endif

