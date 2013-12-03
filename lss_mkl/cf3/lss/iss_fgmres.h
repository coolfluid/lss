// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_mkl_iss_fgmres_h
#define cf3_lss_mkl_iss_fgmres_h


#include "LibLSS_MKL.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


/**
 * @brief Interface to Intel MKL iterative sparse solvers, using RCI interface
 * to implement a non-preconditioned, ILU0 or ILUT-preconditioned flexible (F)
 * GMRES solver
 * @author Pedro Maciel
 */
class lss_API iss_fgmres : public
  linearsystem< double >
{
  // utility definitions
  typedef sparse_matrix< double, sort_by_row, 1 > matrix_t;
  enum pc_t { NONE=0, ILU0, ILUT };


 public:
  // framework interfacing
  static std::string type_name() { return "iss_fgmres"; }

  /// Construction
  iss_fgmres(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1 );

  /// Destruction
  ~iss_fgmres();

  /// Linear system solving
  iss_fgmres& solve();

  /// Linear system copy
  iss_fgmres& copy(const iss_fgmres& _other);


  // internal functions
 private:

  /// Verbose error message
  static const std::string err_message(const int& err, const pc_t& _pc_type);


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
  matrix_t                      m_A;           // system matrix
  matrix_t::matrix_compressed_t m_pc;          // preconditioner matrix
  std::string                   m_pc_type;     // ... type name
  bool                          m_pc_refresh;  // ... force recalculation
  bool                          m_monitor;     // monitor iterations
  double                        m_resnorm;     // maximum residual norm (if test_user is set)
  struct {
    pc_t   pc_type;  // preconditioner type
    double tol;      // ILUT threshold
    int    maxfil;   // ILUT maximum fill-in (half-bandwidth)
    int
      iparm__4,
      iparm__7,
      iparm__8,
      iparm__9,
      iparm_11,
      iparm_14;
    double
      dparm__0,
      dparm__1;
  } opt,                      // options (current)
    previous_opt;             // options (previous, for caching)
  std::vector< double > tmp;  // space for computations
  int    iparm[128];
  double dparm[128];

};


}  // namespace mkl
}  // namespace lss
}  // namespace cf3


#endif

