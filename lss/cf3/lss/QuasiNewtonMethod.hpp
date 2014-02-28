// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_QuasiNewtonMethod_hpp
#define cf3_lss_QuasiNewtonMethod_hpp


#include "LibLSS.hpp"
#include "nonlinearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Example non-linear system solver, using Gaussian elimination
 * (configurable precision)
 * @author Pedro Maciel
 */
template< typename T >
class lss_API QuasiNewtonMethod : public nonlinearsystem< T >
{
 public:

  /// Component type name (framework interfacing)
  static std::string type_name();

  /// Construction
  QuasiNewtonMethod(const std::string& name) : nonlinearsystem< T >(name) {
    //TODO
  }

  /// Non-linear system solving
  QuasiNewtonMethod& solve() {
    //TODO
    return *this;
  }

};


}  // namespace lss
}  // namespace cf3


#endif

