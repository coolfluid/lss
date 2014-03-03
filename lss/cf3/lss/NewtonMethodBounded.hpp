// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_NewtonMethodBounded_hpp
#define cf3_lss_NewtonMethodBounded_hpp


#include "LibLSS.hpp"
#include "nonlinearsystem.hpp"


namespace cf3 {
namespace lss {


/**
 * @brief Non-linear solver using the Newton Method for which the solution
 * update is relaxed such that resulting solution is within configurable bounds
 * (configurable precision)
 * @author Pedro Maciel
 */
template< typename T >
class lss_API NewtonMethodBounded : public nonlinearsystem< T >
{

 public:

  /// Component type name (framework interfacing)
  static std::string type_name();

  /// Construction
  NewtonMethodBounded(const std::string& name) : nonlinearsystem< T >(name) {
    //TODO
  }

  /// Non-linear system solving
  NewtonMethodBounded& solve() {
    double relax = 1.;

    //TODO
#if 0
    // update solution fields -- dry-run to check bounds
    std::vector< Handle<System> > systems = get_systems();

    for (size_t i=0; i<systems.size(); ++i) {
      std::vector< Handle<Field> > unknowns = systems[i]->get_unknowns();
      for (size_t j=0; j<unknowns.size(); ++j) {
        Field& unknown = *unknowns[j];
        std::vector<double> const& bounds = unknown.options().value<std::vector<double> >("bounds");
        if (bounds[0] == -math::Consts::inf() && bounds[1] == math::Consts::inf()) continue;
        boost_foreach(Uint const& n, systems[i]->used_nodes()) {
          double const
              &dU  = ls.x(m_idx_unknowns[n][i]+j),
              &Un  = unknown[n],
              Unp1 = Un+dU;
          relax = (Unp1<bounds[0]?std::min(std::max(0.1,(bounds[0]-Un)/dU), relax):
                  (Unp1>bounds[1]?std::min(std::max(0.1,(bounds[1]-Un)/dU), relax):
                                                                            relax));
        }
      }
    }

    cf3_always_assert_desc("Relaxation < 0! Probably the previous solution is not within the specified bounds.", relax >= 0.);

#endif

    // FIXME only one column is updated in solution vector (matrix)?
    typename linearsystem< T >::vector_t &x = this->linearsystem_get().x();
    for (size_t j=0; j<x.size(0); ++j)
      x(j) *= relax;
    return *this;
  }

};


}  // namespace lss
}  // namespace cf3


#endif

