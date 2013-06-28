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
  GaussianElimination::resize(3,3);
  A(0,0) = 1; A(0,1) = 1; A(0,2) =  1; b(0) =  6;
  A(1,0) = 0; A(1,1) = 2; A(1,2) =  5; b(1) = -4;
  A(2,0) = 2; A(2,1) = 5; A(2,2) = -1; b(2) = 27;
  std::cout << *this << std::endl;


  const size_t N(m_A.size(0));
  double C;
  m_x = m_b;

  boost::progress_display pbar(N-1);
  for (size_t m=0; m<N-1; ++m, ++pbar) {

    // put row with highest diagonal element on top
    C = A(m,m);
    for (size_t n=m+1; n<N; n++) {
      if (std::abs(A(n,m)) > std::abs(C)) {
        for (size_t p=m; p<N; p++) {
          C = A(m,p);
          A(m,p) = A(n,p);
          A(n,p) = C;
        }
        std::swap(x(m),x(n));
        C = A(m,m);
      }
    }

    // check if diagonal element is (close to) zero
    if (std::abs(C)<1.e-32) {
      std::cerr << "error: matrix is singular (line:" << m << ",C:" << std::abs(C) << ")" << std::endl;
      throw 42;
    }

    // normalize row m
    for (size_t n=m+1; n<N; n++)
      A(m,n) /= C;
    x(m) /= C;

    // subtract row m from subsequent rows
    for (size_t n=m+1; n<N; n++) {
      C = A(n,m);
      for (size_t p=m+1; p<N; p++)
        A(n,p) -= C*A(m,p);
      x(n) -= C*x(m);
    }
  }

  // solve by back substitution
  x(N-1) /= A(N-1,N-1);
  for (size_t p=0; p<N-1; p++) {
    size_t m = N-p-2;
    for (size_t n=m+1; n<N; n++)
      x(m) -= A(m,n)*x(n);
  }

  std::cout << *this << std::endl;
  return *this;
}


}  // namespace lss
}  // namespace cf3

