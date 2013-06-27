// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cmath>
#include <algorithm>
#include "boost/progress.hpp"

#include "common/Builder.hpp"
#include "GaussianElimination.h"


namespace cf3 {
namespace lss {


common::ComponentBuilder< GaussianElimination, common::Component, LibLSS > GaussianElimination_Builder;


GaussianElimination& GaussianElimination::solve()
{
  CFinfo << "ping pong" << name() << "ping pong" << CFendl;
  return *this;

  m_x = m_b;
  const unsigned Ne = 1;
  const unsigned N = Ne;

  double C;
  boost::progress_display pbar(N-1);
  for (unsigned m=0; m<N-1; ++m, ++pbar) {

    // put row with highest diagonal element on top
    C = A(m,m);
    for (unsigned n=m+1; n<N; n++) {
      if (std::abs(A(n,m)) > std::abs(C)) {
        for (unsigned p=m; p<N; p++) {
          C = A(m,p);
          A(m,p) = A(n,p);
          A(n,p) = C;
        }
        C    = x(m);
        x(m) = x(n);
        x(n) = C;
        C    = A(m,m);
      }
    }

    // check if diagonal element is (close to) zero
    if (std::abs(C)<1.e-32) {
      std::cerr << "error: matrix is singular (line:" << m << ",C:" << std::abs(C) << ")" << std::endl;
      throw 42;
    }

    // normalize row m
    for (unsigned n=m+1; n<N; n++)
      A(m,n) /= C;
    x(m) /= C;

    // subtract row m from subsequent rows
    for (unsigned n=m+1; n<N; n++) {
      C = A(n,m);
      for (unsigned p=m+1; p<N; p++)
        A(n,p) -= C*A(m,p);
      x(n) -= C*x(m);
    }
  }

  // solve by back substitution
  x(N-1) /= A(N-1,N-1);
  for (unsigned p=0; p<N-1; p++) {
    unsigned m = N-p-2;
    for (unsigned n=m+1; n<N; n++)
      x(m) -= A(m,n)*x(n);
  }

  return *this;
}


}  // namespace lss
}  // namespace cf3

