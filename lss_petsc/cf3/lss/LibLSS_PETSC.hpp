// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_LibLSS_PETSC_hpp
#define cf3_lss_LibLSS_PETSC_hpp


#define SANDBOX
//#undef SANDBOX


#include "cf3/common/Library.hpp"
#include "../../../lss/cf3/lss/LibLSS.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief LibLSS_PETSC class: Interface to PETSc.
 * @author Pedro Maciel
 */
struct lss_API LibLSS_PETSC :
  public cf3::common::Library
{
  /// Constructor
  LibLSS_PETSC(const std::string& name) : cf3::common::Library(name) {}

  /// @return library namespace, name, description and class name
  static std::string library_namespace()    { return "cf3.lss.petsc"; }
  static std::string library_name()         { return "petsc"; }
  static std::string library_description()  { return "Interface to PETSc linear system solver."; }
  static std::string type_name()            { return "LibLSS_PETSC"; }

  /// Initiate library
  void initiate();
};


}  // lss
}  // cf3


#endif // cf3_lss_LibLSS_PETSC_hpp

