// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_GaussianElimination_h
#define cf3_lss_GaussianElimination_h


#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


// implementation of a linear system solver, using gaussian elimination (double p.)
class lss_API GaussianElimination :
    public linearsystem< double, matrix_dense_aa< double > >
{
  typedef linearsystem< double, matrix_dense_aa< double > > linearsystem_t;

 public:
  // framework interfacing
  static std::string type_name() { return "GaussianElimination"; } // (mandatory!)
  GaussianElimination(const std::string& name) :
    linearsystem_t(name) {}

 public:
  // linear system solver interfacing
  GaussianElimination& solve();

};


}  // namespace lss
}  // namespace cf3


#endif
