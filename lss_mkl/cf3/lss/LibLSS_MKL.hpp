// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_mkl_LibLSS_MKL_hpp
#define cf3_lss_mkl_LibLSS_MKL_hpp


#define SANDBOX
//#undef SANDBOX


#include "cf3/common/Library.hpp"
#include "../../../lss/cf3/lss/LibLSS.hpp"


namespace cf3 {
namespace lss {
namespace mkl {


/**
 * @brief LibLSS_MKL class: Interface to Intel MKL iterative and direct solvers.
 * @author Pedro Maciel
 */
struct lss_API LibLSS_MKL :
  public cf3::common::Library
{
  /// Constructor
  LibLSS_MKL(const std::string& name) : cf3::common::Library(name) {}

  /// @return library namespace, name, description and class name
  static std::string library_namespace()    { return "cf3.lss.mkl"; }
  static std::string library_name()         { return "mkl"; }
  static std::string library_description()  { return "Interface to Intel MKL iterative and direct solvers."; }
  static std::string type_name()            { return "LibLSS_MKL"; }

  /// Initiate library
  void initiate();
};


}  // mkl
}  // lss
}  // cf3


#endif // cf3_lss_mkl_LibLSS_MKL_hpp

