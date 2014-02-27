// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "mkl_service.h"

#include "common/Log.hpp"
#include "detail_solverbase.h"


namespace cf3 {
namespace lss {
namespace mkl {
namespace detail {


solverbase::solverbase(const std::string& name)
  : linearsystem< double >(name)
{
  environment_variable_t< int > nthreads("OMP_NUM_THREADS",1);
  CFinfo << "mkl: OMP_NUM_THREADS: " << nthreads.description() << CFendl;
  mkl_set_num_threads(nthreads.value);
}


}  // namespace detail
}  // namespace mkl
}  // namespace lss
}  // namespace cf3

