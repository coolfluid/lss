// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_LibLSS_WSMP_hpp
#define cf3_lss_LibLSS_WSMP_hpp


#define SANDBOX
//#undef SANDBOX


#include "cf3/common/Library.hpp"
#include "../../../lss/cf3/lss/LibLSS.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief LibLSS_WSMP class: Interface to WSMP (serial version).
 * @author Pedro Maciel
 */
struct lss_API LibLSS_WSMP :
  public cf3::common::Library
{
  /// Constructor
  LibLSS_WSMP(const std::string& name) : cf3::common::Library(name) {}

  /// @return library namespace, name, description and class name
  static std::string library_namespace()    { return "cf3.lss.wsmp"; }
  static std::string library_name()         { return "wsmp"; }
  static std::string library_description()  { return "Interface to WSMP linear system solver (serial version)."; }
  static std::string type_name()            { return "LibLSS_WSMP"; }

  /// Initiate library
  void initiate();
};


}  // lss
}  // cf3


#endif // cf3_lss_LibLSS_WSMP_hpp

