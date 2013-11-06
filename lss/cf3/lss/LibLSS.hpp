// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_LibLSS_hpp
#define cf3_lss_LibLSS_hpp


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


/**
 * @brief lss namespace
 *
 * The namespace lss provides a bare-bones interface to linear system solvers,
 * of dense and sparse matrices in different precisions. A few implementations
 * are provided in the core plugin and others, more sophisticated, are
 * available as dependent plugins.
 * @author Pedro Maciel
 */
namespace lss {


struct lss_API LibLSS :
  public cf3::common::Library
{
  /// Constructor
  LibLSS(const std::string& name) : cf3::common::Library(name) {}

  /// @return library namespace, name, description and class name
  static std::string library_namespace()    { return "cf3.lss"; }
  static std::string library_name()         { return "lss"; }
  static std::string library_description()  { return "Basic interface to linear system solvers."; }
  static std::string type_name()            { return "LibLSS"; }

  /// Initiate library
  void initiate();
};


}  // lss
}  // cf3


#endif // cf3_lss_LibLSS_h

