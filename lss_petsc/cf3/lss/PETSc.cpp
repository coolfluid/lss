// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf

#include "common/Builder.hpp"
#include "PETSc.hpp"


namespace cf3 {
namespace lss {


common::ComponentBuilder< PETSc, common::Component, LibLSS_PETSC > Builder_PETSc;


PETSc::PETSc(
    const std::string& name,
    const size_t& _size_i,
    const size_t& _size_j,
    const size_t& _size_k,
    const double& _value) : linearsystem_t(name)
{
  linearsystem_t::initialize(_size_i,_size_j,_size_k,_value);
}


PETSc& PETSc::solve()
{
//  nrhs = static_cast< int >(b().size(1));
  return *this;
}


}  // namespace lss
}  // namespace cf3

