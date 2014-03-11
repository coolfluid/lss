// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_mkl_pardiso_h
#define cf3_lss_mkl_pardiso_h


#include "LibLSS_MKL.hpp"
#include "detail_solverbase.h"


namespace cf3 {
namespace lss {
namespace mkl {


/**
 * @brief Interface to Pardiso linear system solver (Intel MKL version).
 * @author Pedro Maciel
 */
class lss_API pardiso : public
  detail::solverbase
{
 public:
  // framework interfacing

  /// Component type name
  static std::string type_name() { return "pardiso"; }

  /// Construction
  pardiso(const std::string& name,
    const size_t& _size_i=size_t(),
    const size_t& _size_j=size_t(),
    const size_t& _size_k=1 );

  /// Destruction
  ~pardiso();

  /// Linear system solving: x = A^-1 b
  pardiso& solve();

  /// Linear system copy
  pardiso& copy(const pardiso& _other);


 private:
  // internal functions and storage

  /// Verbose error message
  static const std::string err_message(const int& err);

  /// Library call
  int call_pardiso(int _phase, int _msglvl);

  void* pt[64];  // internal memory pointer (void* for both 32/64-bit)
  int   iparm[64],
        maxfct,
        mnum,
        mtype;

};


}  // namespace mkl
}  // namespace lss
}  // namespace cf3


#endif

