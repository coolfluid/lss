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
 * @brief Non-linear solver using the Newton Method together with a line search
 * strategy for relaxing the solution update (configurable precision)
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

    // line search: 3-point safeguarded parabolic model
    const Uint kmax = 50;
    double
      lambdac = 1.,       // λ at k, current best guess
      lambdam = lambdac;  // ... at k-1

    // safeguarding bounds
    const double
        sigma0 = 0.1,
        sigma1 = 0.5,
        alpha  = 1.e-4;

    //TODO
    // unperturbed/perturbed residuals (F(x) and F(x+λΔx)), and their norms
    typename linearsystem< T >::vector_t
     &F0 = this->linearsystem_get().b(),
      Ft = F0; // FIXME
    double
      normF0  = F0.norm(),
      normFtc = Ft.norm(),
      normFtm(normFtc);

    // line search
    for (Uint k=0; (k<kmax) && (normFtc>=(1.-alpha*lambdac)*normF0); ++k) {
      if (!k) { lambdac *= sigma1; }
      else {
        Ft = 1.;       // FIXME
        normFtc = Ft.norm();

        const double
          F0F0  = normF0 *normF0,
          FtFtc = normFtc*normFtc,
          FtFtm = normFtm*normFtm,
          c1 = (FtFtm-F0F0)*lambdac*lambdac - (FtFtc-F0F0)*lambdam*lambdam,
          c2 = (FtFtc-F0F0)*lambdam         - (FtFtm-F0F0)*lambdac;

        lambdac = (c2>=0.? sigma1*lambdac
                : std::min(sigma1*lambdac,
                  std::max(sigma0*lambdac, -c1*0.5/c2 )) );

      }
      normFtm = normFtc;
      lambdam = lambdac;
    }

    // FIXME only one column is updated in solution vector (matrix)?
    typename linearsystem< T >::vector_t &x = this->linearsystem_get().x();
    for (size_t j=0; j<x.size(0); ++j)
      x(j) *= lambdac;
    return *this;
  }

};


}  // namespace lss
}  // namespace cf3


#endif

