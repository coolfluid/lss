
#include <cmath>
#include <algorithm>
#include "boost/progress.hpp"
//include "mfactory.h"
#include "mlinearsystem.hpp"


//Register< mlinearsystem< double >,GaussianElimination > __GaussianElimination("GaussianElimination","linear system solver, using gaussian elimination (double p.)");


namespace cf3 {
namespace lss {


void GaussianElimination::solve()
{
  const unsigned N = Ne*Nb;
  m_X = m_B;

  double C;
  boost::progress_display pbar(N-1);
  for (unsigned m=0; m<N-1; ++m, ++pbar) {

    // put row with highest diagonal element on top
    C = A(m,m);
    for (unsigned n=m+1; n<N; n++) {
      if (std::fabs(A(n,m)) > std::fabs(C)) {
        for (unsigned p=m; p<N; p++) {
          C = A(m,p);
          A(m,p) = A(n,p);
          A(n,p) = C;
        }
        C    = X(m);
        X(m) = X(n);
        X(n) = C;
        C    = A(m,m);
      }
    }

    // check if diagonal element is (close to) zero
    if (std::fabs(C)<1.e-32) {
      std::cerr << "error: matrix is singular (line:" << m << ",C:" << std::fabs(C) << ")" << std::endl;
      throw 42;
    }

    // normalize row m
    for (unsigned n=m+1; n<N; n++)
      A(m,n) /= C;
    X(m) /= C;

    // subtract row m from subsequent rows
    for (unsigned n=m+1; n<N; n++) {
      C = A(n,m);
      for (unsigned p=m+1; p<N; p++)
        A(n,p) -= C*A(m,p);
      X(n) -= C*X(m);
    }
  }

  // solve by back substitution
  X(N-1) /= A(N-1,N-1);
  for (unsigned p=0; p<N-1; p++) {
    unsigned m = N-p-2;
    for (unsigned n=m+1; n<N; n++)
      X(m) -= A(m,n)*X(n);
  }
}


}  // namespace lss
}  // namespace cf3

