// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_pardiso_hpp
#define cf3_lss_pardiso_hpp


#include "LibLSS_PARDISO.hpp"
#include "../../../lss/cf3/lss/linearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Interface to Pardiso linear system solver (U. Basel version).
 * @author Pedro Maciel
 */
class lss_API pardiso : public linearsystem< double >
{
  // utility definitions
  typedef sparse_matrix< double, sort_by_row, 1 > matrix_t;


 public:
  // framework interfacing
  static std::string type_name() { return "pardiso"; }


  /// Construction
  pardiso(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1,
    const double& _value=double() );

  /// Destruction
  ~pardiso();

  /// Linear system solving
  pardiso& solve();

  /// Linear system copy
  pardiso& copy(const pardiso& _other);


  // internal functions
 private:

  /// Verbose error message
  static std::string err_message(const int& err);

  /// Library call
  int call_pardiso(int _phase, int _msglvl);

  /// Library call (print statistics)
  int call_pardiso_printstats();


 protected:
  // linear system matrix interfacing

  const double& A(const size_t& i, const size_t& j) const { return m_A(i,j); }
        double& A(const size_t& i, const size_t& j)       { return m_A(i,j); }

  void A___initialize(const size_t& i, const size_t& j, const double& _value=double()) { m_A.initialize(i,j,_value); }
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
  void*  pt[64];  // internal memory pointer (void* for both 32/64-bit)
  double dparm[64];
  int    iparm[64],
         maxfct,
         mnum,
         mtype;

};


}  // namespace lss
}  // namespace cf3


#endif

