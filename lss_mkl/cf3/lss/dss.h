// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_mkl_dss_h
#define cf3_lss_mkl_dss_h


#include "LibLSS_MKL.hpp"
#include "detail_mkl_solver_base.h"


namespace cf3 {
namespace lss {
namespace mkl {


/**
 * @brief Interface to Intel MKL direct sparse solvers.
 * @author Pedro Maciel
 */
class lss_API dss : public
  detail::mkl_solver_base
{
 private:
  // utility definitions
  enum phase_t { _CREATE=0, _STRUCTURE, _REORDER, _FACTOR, _SOLVE, _DELETE, _ALL_PHASES };


 public:
  // framework interfacing

  /// Component type name
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


 private:
  // internal functions and storage

  /// Verbose error message
  static std::string err_message(const int& err);

  int opts[_ALL_PHASES];
  void *handle;

};


}  // namespace mkl
}  // namespace lss
}  // namespace cf3


#endif

