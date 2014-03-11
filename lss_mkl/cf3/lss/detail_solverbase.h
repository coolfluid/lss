// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_mkl_detail_solverbase_h
#define cf3_lss_mkl_detail_solverbase_h


#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {
namespace mkl {
namespace detail {


/**
 * @brief MKL solvers management of sparse matrix and common operations
 * @author Pedro Maciel
 */
struct solverbase : public linearsystem< double >
{
  /// Matrix storage
  typedef sparse_matrix< double, sort_by_row, 1 > matrix_t;
  matrix_t m_A;

  /// Construction
  solverbase(const std::string& name);

  /// Linear system forward multiplication: b = alpha A x + beta b
  solverbase& multi(const double& _alpha=1., const double& _beta=0.);

  /// Linear system swap
  solverbase& swap(solverbase& _other);

  /// Matrix indexing
  const double& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        double& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  /// Matrix operations
  void A___initialize(const size_t& i, const size_t& j, const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >()) { m_A.initialize(i,j,_nnz); }
  void A___initialize(const std::vector< double >& _vector) { m_A.initialize(_vector); }
  void A___initialize(const std::string& _fname)            { m_A.initialize(_fname);  }
  void A___assign(const double& _value) { m_A = _value;   }
  void A___clear()                      { m_A.clear();    }
  void A___zerorow(const size_t& i)     { m_A.zerorow(i); }
  void A___sumrows(const size_t& i, const size_t& isrc) { m_A.sumrows(i,isrc); }

  /// Matrix utilities
  void   A___print(std::ostream& o, const print_t& l=print_auto) const { m_A.print(o,l); }
  size_t A___size(const size_t& d)  const { return m_A.size(d);  }
};


}  // namespace detail
}  // namespace mkl
}  // namespace lss
}  // namespace cf3


#endif

