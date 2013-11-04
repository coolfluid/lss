// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_LibLSS_pardiso_hpp
#define cf3_lss_LibLSS_pardiso_hpp


#define SANDBOX
//#undef SANDBOX


#include "cf3/common/Library.hpp"


/// Define the macro lss_API
/// @note build system defines COOLFLUID_LSS_EXPORTS when compiling lss files
#ifdef COOLFLUID_LSS_EXPORTS
#   define lss_API      CF3_EXPORT_API
#   define lss_TEMPLATE
#else
#   define lss_API      CF3_IMPORT_API
#   define lss_TEMPLATE CF3_TEMPLATE_EXTERN
#endif

////////////////////////////////////////////////////////////////////////////////

namespace cf3 {
namespace lss {


/**
 * @brief LibLSS_pardiso class
 * This library provides an interface to Pardiso linear system solver (U. Basel version).
 * @author Pedro Maciel
 */
struct lss_API LibLSS_pardiso :
  public cf3::common::Library
{
  /// Constructor
  LibLSS_pardiso(const std::string& name) : cf3::common::Library(name) {}

  /// @return string of the library namespace
  static std::string library_namespace() { return "cf3.lss.pardiso"; }

  /// @return name of the library
  static std::string library_name() { return "pardiso"; }

  /// @return description of the library
  static std::string library_description() { return "This library provides an interface to Pardiso linear system solver (U. Basel version)."; }

  /// Gets the Class name
  static std::string type_name() { return "LibLSS_pardiso"; }

  /// Initiate library
  void initiate();

};


}  // lss
}  // cf3


#endif // cf3_lss_LibLSS_pardiso_h

