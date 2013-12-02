// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_mkl_dss_h
#define cf3_lss_mkl_dss_h


#include "LibLSS_MKL.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


/**
 * @brief Interface to Intel MKL direct sparse solvers.
 * @author Pedro Maciel
 */
class lss_API dss : public
  linearsystem< double >
{
  // utility definitions
  typedef sparse_matrix< double, sort_by_row, 1 > matrix_t;


 public:
  // framework interfacing
  static std::string type_name() { return "dss"; }


  /// Construction
  dss(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1 );

  /// Destruction
  ~dss();

  /// Linear system solving
  dss& solve();

  /// Linear system copy
  dss& copy(const dss& _other);


  // internal functions
 private:

  /// Verbose error message
  static std::string err_message(const int& err);


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
  // phase options
  enum phase_t { _create=0, _structure, _reorder, _factor, _solve, _delete, _all_phases };

 protected:
  // storage
  matrix_t m_A;
  int opts[_all_phases];
  void *handle;

};


}  // namespace mkl
}  // namespace lss
}  // namespace cf3


#endif

