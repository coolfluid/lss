// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cstdio>  // for sscanf

#include "common/Builder.hpp"
#include "PETSc_Seq.hpp"


namespace cf3 {
namespace lss {


common::ComponentBuilder< PETSc_Seq, common::Component, LibLSS_PETSC > Builder_PETSc_Seq;


PETSc_Seq::PETSc_Seq(const std::string& name,
    const size_t& _size_i,
    const size_t& _size_j,
    const size_t& _size_k)
  : linearsystem< double >(name)
{
  linearsystem< double >::initialize(_size_i,_size_j,_size_k);
}


PETSc_Seq& PETSc_Seq::solve()
{
//  nrhs = static_cast< int >(m_b.size(1));
  return *this;
}

PETSc_Seq&PETSc_Seq::copy(const PETSc_Seq& _other)
{
  linearsystem< double >::copy(_other);
  m_A = _other.m_A;
  return *this;
}


}  // namespace lss
}  // namespace cf3

