// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include <cmath>

#include "common/Builder.hpp"
#include "GMRES.hpp"


namespace cf3 {
namespace lss {


std::string GMRES::type_name() { return "GMRES"; }
common::ComponentBuilder< GMRES, common::Component, LibLSS > Builder_GMRES;


GMRES& GMRES::solve()
{
  int n = m_A.size(0);
  int iwk = m_A.nnz;
  int ierr;
  double eps = 1e-5;  // tolerance, process is stopped when eps>=||current residual||/||initial residual||
  int im     = 50;    // size of krylov subspace (should not exceed 50)
  int maxits = 50;    // maximum number of iterations allowed
  int newiwk = 0; newiwk=newiwk;
  int iout   = 1;
  int lfil   = 3;

  std::vector< double > alu(iwk,0.);
  std::vector< int >
    jlu(iwk+1,0),
    ju(n+1,0);

  {
    std::vector< double > w(n+1,0.);
    std::vector< int >
      levs(iwk),
      jw(n*3,0);
    /*newiwk =*/ iluk(&n,&m_A.a[0],&m_A.ja[0],&m_A.ia[0],&lfil,alu,jlu,
      &ju[0],&levs[0],&iwk,&w[0],&jw[0],&ierr);
  }
  {
    std::vector< double > vv(n*(im+1));
    x() = 0.;
    pgmres(&n,&im,&(b().a[0]),&(x().a[0]),&vv[0],&eps,&maxits,&iout,
           &m_A.a[0],&m_A.ja[0],&m_A.ia[0],&alu[0],&jlu[0],&ju[0],&ierr);
  }

  return *this;
}


int GMRES::daxpy(int *n, double *da, double *dx, int *incx, double *dy, int *incy)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
  return 0;
    }
    if (*da == 0.) {
  return 0;
    }
    if (*incx == 1 && *incy == 1) {
  goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
  ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
  iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
  dy[iy] += *da * dx[ix];
  ix += *incx;
  iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
  goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
  dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
  return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
  dy[i__] += *da * dx[i__];
  dy[i__ + 1] += *da * dx[i__ + 1];
  dy[i__ + 2] += *da * dx[i__ + 2];
  dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;
}


/*=========================================================================*
 *     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k))    *
 *=========================================================================*
 *                                                                         *
 * on entry:                                                               *
 * ==========                                                              *
 * n       = integer. The row dimension of the matrix A. The matrix        *
 * a,ja,ia = matrix stored in Compressed Sparse Row format.                *
 * lfil    = integer. The fill-in parameter. Each element whose            *
 *           leve-of-fill exceeds lfil during the ILU process is dropped.  *
 *           lfil must be .ge. 0                                           *
 * tol     = real*8. Sets the threshold for dropping small terms in the    *
 *           factorization. See below for details on dropping strategy.    *
 * iwk     = integer. The minimum length of arrays alu, jlu, and levs.     *
 *                                                                         *
 * On return:                                                              *
 * ===========                                                             *
 * alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing  *
 *           the L and U factors together. The diagonal (stored in         *
 *           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix   *
 *           contains the i-th row of L (excluding the diagonal entry=1)   *
 *           followed by the i-th row of U.                                *
 * ju      = integer array of length n containing the pointers to          *
 *           the beginning of each row of U in the matrix alu,jlu.         *
 * levs    = integer (work) array of size iwk -- which contains the        *
 *           levels of each element in alu, jlu.                           *
 * ierr    = integer. Error message with the following meaning.            *
 *           ierr  = 0    --> successful return.                           *
 *           ierr .gt. 0  --> zero pivot encountered at step number ierr.  *
 *           ierr  = -1   --> Error. input matrix may be wrong.            *
 *                            (The elimination process has generated a     *
 *                            row in L or U whose length is .gt.  n.)      *
 *           ierr  = -2   --> The matrix L overflows the array al.         *
 *           ierr  = -3   --> The matrix U overflows the array alu.        *
 *           ierr  = -4   --> Illegal value for lfil.                      *
 *           ierr  = -5   --> zero row encountered in A or U.              *
 *                                                                         *
 * work arrays:                                                            *
 * =============                                                           *
 * jw      = integer work array of length 3*n.                             *
 * w       = real work array of length n                                   *
 *                                                                         *
 * Notes/known bugs: This is not implemented efficiently storage-wise.     *
 *       For example: Only the part of the array levs(*) associated with   *
 *       the U-matrix is needed in the routine.. So some storage can       *
 *       be saved if needed. The levels of fills in the LU matrix are      *
 *       output for information only -- they are not needed by LU-solve.   *
 *                                                                         *
 *=========================================================================*
 *                                                                         *
 * w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]         *
 * jw(n+1:2n)  stores the nonzero indicator.                               *
 *                                                                         *
 * Notes:                                                                  *
 * ------                                                                  *
 * All the diagonal elements of the input matrix must be  nonzero.         *
 *                                                                         *
 *=========================================================================*/
int GMRES::iluk(int *n, double *a, int *ja, int *ia, int *lfil, std::vector<double>& aluold, std::vector<int>& jluold, int *ju, int* levsold, int *iwk, double *w, int *jw, int *ierr)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4;

  /* Local variables */
  static double fact;
  static int lenl, jlev, lenu, jpos, jrow, i__, j, k;
  static double s, t;
  static int j1, j2, n2, ii, jj, ju0;

  std::vector<double> alu((*iwk)+2,0.);
  std::vector<int> jlu((*iwk)+2,0);
  std::vector<int> levs((*iwk)+2,0);

/*     locals */
    /* Parameter adjustments */
  --jw;
  --w;
  --ju;
  --ia;
  --a;
  --ja;
//    --alu;
//    --jlu;
//    --levs;

   /* Function Body */
  if (*lfil < 0) {
    goto L998;
  }
/* ----------------------------------------------------------------------- */
/*     initialize ju0 (points to next element to be added to alu,jlu) */
/*     and pointer array. */
/* ----------------------------------------------------------------------- */
  n2 = *n + *n;
  ju0 = *n + 2;
  jlu[1] = ju0;

/*     initialize nonzero indicator array + levs array -- */

  i__1 = *n << 1;
  for (j = 1; j <= i__1; ++j) {
    jw[j] = 0;
/* L1: */
  }
/* ----------------------------------------------------------------------- */
/*     beginning of main loop. */
/* ----------------------------------------------------------------------- */
  i__1 = *n;
  for (ii = 1; ii <= i__1; ++ii) {
    j1 = ia[ii];
    j2 = ia[ii + 1] - 1;

/*     unpack L-part and U-part of row of A in arrays w */

    lenu = 1;
    lenl = 0;
    jw[ii] = ii;
    w[ii] = (double)0.;
    jw[*n + ii] = ii;

    i__2 = j2;
    for (j = j1; j <= i__2; ++j) {
      k = ja[j];
      t = a[j];
      if (t == (double)0.) {
        goto L170;
      }
      if (k < ii) {
        ++lenl;
        jw[lenl] = k;
        w[lenl] = t;
        jw[n2 + lenl] = 0;
        jw[*n + k] = lenl;
      }
      else if (k == ii) {
        w[ii] = t;
        jw[n2 + ii] = 0;
      }
      else {
        ++lenu;
        jpos = ii + lenu - 1;
        jw[jpos] = k;
        w[jpos] = t;
        jw[n2 + jpos] = 0;
        jw[*n + k] = jpos;
      }
    L170:
      ;
    }

    jj = 0;

/*     eliminate previous rows */

    L150:
      ++jj;
    if (jj > lenl) {
      goto L160;
    }
/* ----------------------------------------------------------------------- */
/*     in order to do the elimination in the correct order we must select */
/*     the smallest column index among jw(k), k=jj+1, ..., lenl. */
/* ----------------------------------------------------------------------- */
    jrow = jw[jj];
    k = jj;

/*     determine smallest column index */

    i__2 = lenl;
    for (j = jj + 1; j <= i__2; ++j) {
      if (jw[j] < jrow) {
        jrow = jw[j];
        k = j;
      }
      /* L151: */
    }

    if (k != jj) {
/*     exchange in jw */
      j = jw[jj];
      jw[jj] = jw[k];
      jw[k] = j;
/*     exchange in jw(n+  (pointers/ nonzero indicator). */
      jw[*n + jrow] = jj;
      jw[*n + j] = k;
/*     exchange in jw(n2+  (levels) */
      j = jw[n2 + jj];
      jw[n2 + jj] = jw[n2 + k];
      jw[n2 + k] = j;
/*     exchange in w */
      s = w[jj];
      w[jj] = w[k];
      w[k] = s;
    }

/*     zero out element in row by resetting jw(n+jrow) to zero. */

    jw[*n + jrow] = 0;

/*     get the multiplier for row to be eliminated (jrow) + its level */

    fact = w[jj] * alu[jrow];
    jlev = jw[n2 + jj];
    if (jlev > *lfil) {
      goto L150;
    }

/*     combine current row and row jrow */

    i__2 = jlu[jrow + 1] - 1;
    for (k = ju[jrow]; k <= i__2; ++k) {
      s = fact * alu[k];
      j = jlu[k];
      jpos = jw[*n + j];
      if (j >= ii) {

/*     dealing with upper part. */

        if (jpos == 0) {

/*     this is a fill-in element */

          ++lenu;
          if (lenu > *n) {
            goto L995;
          }
          i__ = ii + lenu - 1;
          jw[i__] = j;
          jw[*n + j] = i__;
          w[i__] = -s;
          jw[n2 + i__] = jlev + levs[k] + 1;
        }
        else {

/*     this is not a fill-in element */

          w[jpos] -= s;
/* Computing MIN */
          i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
          jw[n2 + jpos] = std::min(i__3,i__4);
        }
      }
      else {

/*     dealing with lower part. */

        if (jpos == 0) {

/*     this is a fill-in element */

          ++lenl;
          if (lenl > *n) {
            goto L995;
          }
          jw[lenl] = j;
          jw[*n + j] = lenl;
          w[lenl] = -s;
          jw[n2 + lenl] = jlev + levs[k] + 1;
        }
        else {

/*     this is not a fill-in element */

          w[jpos] -= s;
/* Computing MIN */
          i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
          jw[n2 + jpos] = std::min(i__3,i__4);
        }
      }
      /* L203: */
    }
    w[jj] = fact;
    jw[jj] = jrow;
    goto L150;
    L160:

/*     reset double-pointer to zero (U-part) */

    i__2 = lenu;
    for (k = 1; k <= i__2; ++k) {
      jw[*n + jw[ii + k - 1]] = 0;
    /* L308: */
    }

/*     update l-matrix */

    i__2 = lenl;
    for (k = 1; k <= i__2; ++k) {

      if (ju0 > *iwk) {

      //cout << "\tju0 = " << ju0 << "\tiwk = " << *iwk << "\tsize van alu = " << alu.size() << endl;

        *ierr = -2;
        //return 0;
        //goto L996;
        alu.push_back(0.);
        jlu.push_back(0);
        levs.push_back(0);
        *iwk += 1;
      }
      if (jw[n2 + k] <= *lfil) {
        alu[ju0] = w[k];
        jlu[ju0] = jw[k];
        ++ju0;
      }
      /* L204: */
    }

  //cout << "einde van for lus : ju0 = " << ju0 << endl;

  //return 0;
  //goto L996;

/*     save pointer to beginning of row ii of U */

    ju[ii] = ju0;

/*     update u-matrix */

    i__2 = ii + lenu - 1;
    for (k = ii + 1; k <= i__2; ++k) {
      if (jw[n2 + k] <= *lfil) {
        if (ju0 > *iwk) {
          alu.push_back(0.);
          jlu.push_back(0);
          levs.push_back(0);
          *iwk += 1;
        }

        jlu[ju0] = jw[k];
        alu[ju0] = w[k];
        levs[ju0] = jw[n2 + k];
        ++ju0;
      }
    /* L302: */
    }
    if (w[ii] == (double)0.) {
      goto L999;
    }

    alu[ii] = 1. / w[ii];

/*     update pointer to beginning of next row of U. */

    jlu[ii + 1] = ju0;
/* ----------------------------------------------------------------------- */
/*     end main loop */
/* ----------------------------------------------------------------------- */
/* L500: */
    }

#if 0
  delete [] aluold;
  delete [] jluold;
  //delete [] levsold;

  aluold = new double[*iwk+1];
  if (aluold == 0)
    {
           //cout << "Not enough memory to allocate lu4" << endl;
    }
  jluold = new int[*iwk+1];
  if (jluold == 0)
    {
           //cout << "Not enough memory to allocate lu5" << endl;
    }
  //levsold = new int[*iwk+1];

  for (int i=0;i<*iwk;++i)
  {
    aluold[i] = alu[i+1];
    jluold[i] = jlu[i+1];
    //levsold[i] = levs[i+1];
  }
#endif
  aluold = alu;
  jluold = jlu;

    *ierr = 0;
  return *iwk;
    //return 0;

/*     incomprehensible error. Matrix must be wrong. */

L995:
  *ierr = -1;
  return 0;

/*     insufficient storage in L. */

/* L996: */
  *ierr = -2;
  return 0;

/*     insufficient storage in U. */

/* L997: */
  *ierr = -3;
  return 0;

// illegal lfil entered.
L998:
  *ierr = -4;
  return 0;

// zero row encountered in A or U.
L999:
  *ierr = -5;
  return 0;
}


double GMRES::ddot(int *n, double *dx, int *incx, double *dy, int *incy)
{
  /* System generated locals */
  int i__1;
  double ret_val;

  /* Local variables */
  static int i__, m;
  static double dtemp;
  static int ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
  return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
  goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
  ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
  iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
  dtemp += dx[ix] * dy[iy];
  ix += *incx;
  iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
  goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
  dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
  goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
  dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
    i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ +
    4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
}


double GMRES::dnrm2(int *n, double *dx, int *incx)
{
    static double zero = 0.;
    static double one = 1.;
    static double cutlo = 8.232e-11;
    static double cuthi = 1.304e19;

    /* Format strings */
/*
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";
*/

    /* System generated locals */
    int i__1, i__2;
    double ret_val, d__1;

    /* Builtin functions */
    //double sqrt();

    /* Local variables */
    static double xmax;
    static int next, i__, j, nn;
    static double hitest, sum;

    /* Assigned format variables */
    /*static char *next_fmt;*/

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     euclidean norm of the n-vector stored in dx() with storage */
/*     increment incx . */
/*     if    n .le. 0 return with result = 0. */
/*     if n .ge. 1 then incx must be .ge. 1 */


    if (*n > 0) {
  goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    /*static char *next_fmt = fmt_30;*/
    sum = zero;
    nn = *n * *incx;
/*                                                 begin main loop */
    i__ = 1;
L20:
    switch ((int)next) {
  case 0: goto L30;
  case 1: goto L50;
  case 2: goto L70;
  case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i__], std::abs(d__1)) > cutlo) {
  goto L85;
    }
    next = 1;
    /*static char *next_fmt = fmt_50;*/
    xmax = zero;

/*                        phase 1.  sum is zero */

L50:
    if (dx[i__] == zero) {
  goto L200;
    }
    if ((d__1 = dx[i__], std::abs(d__1)) > cutlo) {
  goto L85;
    }

/*                                prepare for phase 2. */
    next = 2;
    /*static char *next_fmt = fmt_70;*/
    goto L105;

/*                                prepare for phase 4. */

L100:
    i__ = j;
    next = 3;
    /*static char *next_fmt = fmt_110;*/
    sum = sum / dx[i__] / dx[i__];
L105:
    xmax = (d__1 = dx[i__], std::abs(d__1));
    goto L115;

/*                   phase 2.  sum is small. */
/*                             scale to avoid destructive underflow. */

L70:
    if ((d__1 = dx[i__], std::abs(d__1)) > cutlo) {
  goto L75;
    }

/*                     common code for phases 2 and 4. */
/*                     in phase 4 sum is large.  scale to avoid overflow. */

L110:
    if ((d__1 = dx[i__], std::abs(d__1)) <= xmax) {
  goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i__];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i__], std::abs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i__] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  prepare for phase 3. */

L75:
    sum = sum * xmax * xmax;


/*     for real or d.p. set hitest = cuthi/n */
/*     for complex      set hitest = cuthi/(2*n) */

L85:
    hitest = cuthi / (double) (*n);

/*                   phase 3.  sum is mid-range.  no scaling. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
  if ((d__1 = dx[j], std::abs(d__1)) >= hitest) {
      goto L100;
  }
/* L95: */
/* Computing 2nd power */
  d__1 = dx[j];
  sum += d__1 * d__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i__ += *incx;
    if (i__ <= nn) {
  goto L20;
    }

/*              end of main loop. */

/*              compute square root and adjust for scaling. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
}


/*=========================================================================*
 *         A times a vector                                                *
 *=========================================================================*
 * multiplies a matrix by a vector using the dot product form              *
 * Matrix A is stored in compressed sparse row storage.                    *
 *                                                                         *
 * on entry:                                                               *
 * ----------                                                              *
 * n     = row dimension of A                                              *
 * x     = real array of length equal to the column dimension of           *
 *         the A matrix.                                                   *
 * a, ja,                                                                  *
 *    ia = input matrix in compressed sparse row format.                   *
 *                                                                         *
 * on return:                                                              *
 * -----------                                                             *
 * y     = real array of length n, containing the product y=Ax             *
 *                                                                         *
 *=========================================================================*/
void GMRES::amux(int *n, double *x, double *y, double *a, int *ja, int *ia)
{
   /* System generated locals */
   int i__1, i__2;

   /* Local variables */
   static int i__, k;
   static double t;

   /* Parameter adjustments */
   --ia;
   --ja;
   --a;
   --y;
   --x;

   /* Function Body */
   i__1 = *n;
   for (i__ = 1; i__ <= i__1; ++i__) {

      /* compute the inner product of row i with vector x */
      t = 0.;
      i__2 = ia[i__ + 1] - 1;
      for (k = ia[i__]; k <= i__2; ++k) {
         t += a[k] * x[ja[k]];
      }

      /* store result in y(i) */
      y[i__] = t;

   }
}


/*=========================================================================*
 *                                                                         *
 * This routine solves the system (LU) x = y,                              *
 * given an LU decomposition of a matrix stored in (alu, jlu, ju)          *
 * modified sparse row format                                              *
 *                                                                         *
 *=========================================================================*
 * on entry:                                                               *
 * n   = dimension of system                                               *
 * y   = the right-hand-side vector                                        *
 * alu, jlu, ju                                                            *
 *     = the LU matrix as provided from the ILU routines.                  *
 *                                                                         *
 * on return                                                               *
 * x   = solution of LU x = y.                                             *
 *=========================================================================*
 *                                                                         *
 * Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)           *
 *       will solve the system with rhs x and overwrite the result on x .  *
 *                                                                         *
 *=========================================================================*/
void GMRES::lusol(int *n, double *y, double *x, double *alu, int *jlu, int *ju)
{
   /* System generated locals */
   int i__1, i__2;

   /* Local variables */
   static int i__, k;

/* local variables */


/* forward solve */

   /* Parameter adjustments */
   --x;
   --y;
   --alu;
   --jlu;
   --ju;
   /* Function Body */
   i__1 = *n;
   for (i__ = 1; i__ <= i__1; ++i__) {
      x[i__] = y[i__];
      i__2 = ju[i__] - 1;
      for (k = jlu[i__]; k <= i__2; ++k) {
         x[i__] -= alu[k] * x[jlu[k]];
      }
   }

   /* backward solve. */
   for (i__ = *n; i__ >= 1; --i__) {
      i__1 = jlu[i__ + 1] - 1;
      for (k = ju[i__]; k <= i__1; ++k) {
         x[i__] -= alu[k] * x[jlu[k]];
      }
      x[i__] = alu[i__] * x[i__];
   }
}


/*=========================================================================*
 *                                                                         *
 *                 *** ILUT - Preconditioned GMRES ***                     *
 *                                                                         *
 *=========================================================================*
 * This is a simple version of the ILUT preconditioned GMRES algorithm.    *
 * The ILUT preconditioner uses a dual strategy for dropping elements      *
 * instead  of the usual level of-fill-in approach. See details in ILUT    *
 * subroutine documentation. PGMRES uses the L and U matrices generated    *
 * from the subroutine ILUT to precondition the GMRES algorithm.           *
 * The preconditioning is applied to the right. The stopping criterion     *
 * utilized is based simply on reducing the residual norm by epsilon.      *
 * This preconditioning is more reliable than ilu0 but requires more       *
 * storage. It seems to be much less prone to difficulties related to      *
 * strong nonsymmetries in the matrix. We recommend using a nonzero tol    *
 * (tol=.005 or .001 usually give good results) in ILUT. Use a large       *
 * lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the       *
 * more reliable the code is. Efficiency may also be much improved.        *
 * Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as    *
 * Gaussian elimination without pivoting.                                  *
 *                                                                         *
 * ILU(0) and MILU(0) are also provided for comparison purposes            *
 * USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and    *
 * then call pgmres.                                                       *
 *=========================================================================*
 * Coded by Y. Saad - This version dated May, 7, 1990.                     *
 *=========================================================================*
 * parameters                                                              *
 * -----------                                                             *
 * on entry:                                                               *
 * ==========                                                              *
 *                                                                         *
 * n     == integer. The dimension of the matrix.                          *
 * im    == size of krylov subspace:  should not exceed 50 in this         *
 *          version (can be reset by changing parameter command for        *
 *          kmax below)                                                    *
 * rhs   == real vector of length n containing the right hand side.        *
 *          Destroyed on return.                                           *
 * sol   == real vector of length n containing an initial guess to the     *
 *          solution on input. approximate solution on output              *
 * eps   == tolerance for stopping criterion. process is stopped           *
 *          as soon as ( ||.|| is the euclidean norm):                     *
 *          || current residual||/||initial residual|| <= eps              *
 * maxits== maximum number of iterations allowed                           *
 * iout  == output unit number number for printing intermediate results    *
 *          if (iout .le. 0) nothing is printed out.                       *
 *                                                                         *
 * aa, ja,                                                                 *
 * ia    == the input matrix in compressed sparse row format:              *
 *          aa(1:nnz)  = nonzero elements of A stored row-wise in order    *
 *          ja(1:nnz) = corresponding column indices.                      *
 *          ia(1:n+1) = pointer to beginning of each row in aa and ja.     *
 *          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)     *
 *                                                                         *
 * alu,jlu== A matrix stored in Modified Sparse Row format containing      *
 *           the L and U factors, as computed by subroutine ilut.          *
 *                                                                         *
 * ju     == integer array of length n containing the pointers to          *
 *           the beginning of each row of U in alu, jlu as computed        *
 *           by subroutine ILUT.                                           *
 *                                                                         *
 * on return:                                                              *
 * ==========                                                              *
 * sol   == contains an approximate solution (upon successful return).     *
 * ierr  == integer. Error message with the following meaning.             *
 *          ierr = 0 --> successful return.                                *
 *          ierr = 1 --> convergence not achieved in itmax iterations.     *
 *          ierr =-1 --> the initial guess seems to be the exact           *
 *                       solution (initial residual computed was zero)     *
 *                                                                         *
 *=========================================================================*
 *                                                                         *
 * work arrays:                                                            *
 * =============                                                           *
 * vv    == work array of length  n x (im+1) (used to store the Arnoli     *
 *          basis)                                                         *
 *=========================================================================*
 * subroutines called :                                                    *
 * amux   : SPARSKIT routine to do the matrix by vector multiplication     *
 *          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux            *
 * lusol : combined forward and backward solves (Preconditioning ope.)     *
 * BLAS1  routines.                                                        *
 *=========================================================================*
 *                                                                         *
 * arnoldi size should not exceed kmax=50 in this version..                *
 * to reset modify paramter kmax accordingly.                              *
 *=========================================================================*/
void GMRES::pgmres(int *n, int *im, double *rhs, double *sol, double *vv, double *eps, int *maxits, int*iout, double *aa, int *ja, int *ia, double *alu, int *jlu, int *ju, int *ierr)
{
   static double epsmac = 1e-16;

   /* System generated locals */
   int vv_dim1, vv_offset, i__1, i__2;
   double d__1, d__2;

   /* Local variables */
   static double c__[50];
   static int i__, j, k;
   static double s[50], t;
   static int i1, k1;
   /*static int n1;*/
   static double hh[2550]  /* was [51][50] */;
   static int ii, jj;
   static double ro, rs[51], gam;
   static int its;
   static double eps1;

   //static AnsiString Message;


   /* Parameter adjustments */
   --ju;
   --ia;
   vv_dim1 = *n;
   vv_offset = 1 + vv_dim1 * 1;
   vv -= vv_offset;
   --sol;
   --rhs;
   --aa;
   --ja;
   --alu;
   --jlu;

   /* Function Body */
   /*static int n1 = *n + 1;*/
   its = 0;

   /* compute initial residual vector */
   amux(n, &sol[1], &vv[vv_offset], &aa[1], &ja[1], &ia[1]);
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      vv[j + vv_dim1] = rhs[j] - vv[j + vv_dim1];
   }

   /* outer loop starts here.. */
L20:
   ro = dnrm2(n, &vv[vv_offset], &c__1);
   if (ro == 0.) {
     *ierr = -1;
     return;
  //goto L999;
   }
   t = 1. / ro;
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      vv[j + vv_dim1] *= t;
   }
   if (its == 0) {
      eps1 = *eps * ro;
   }

   /* initialize 1-st term  of rhs of hessenberg system.. */
   rs[0] = ro;
   i__ = 0;

L4:
   ++i__;
   ++its;
   i1 = i__ + 1;
   lusol(n, &vv[i__ * vv_dim1 + 1], &rhs[1], &alu[1], &jlu[1], &ju[1]);
   amux(n, &rhs[1], &vv[i1 * vv_dim1 + 1], &aa[1], &ja[1], &ia[1]);

   /* modified gram - schmidt... */
   i__1 = i__;
   for (j = 1; j <= i__1; ++j) {
      t = ddot(n, &vv[j * vv_dim1 + 1], &c__1, &vv[i1 * vv_dim1 + 1], &c__1);
      hh[j + i__ * 51 - 52] = t;
      d__1 = -t;
      daxpy(n, &d__1, &vv[j * vv_dim1 + 1], &c__1, &vv[i1 * vv_dim1 + 1], &c__1);
   }
   t = dnrm2(n, &vv[i1 * vv_dim1 + 1], &c__1);
   hh[i1 + i__ * 51 - 52] = t;
   if (t == 0.) {
      goto L58;
   }
   t = 1. / t;
   i__1 = *n;
   for (k = 1; k <= i__1; ++k) {
      vv[k + i1 * vv_dim1] *= t;
   }

   /* done with modified gram schimd and arnoldi step.. */
   /* now  update factorization of hh */
L58:
   if (i__ == 1) {
      goto L121;
   }

   /* perfrom previous transformations on i-th column of h */
   i__1 = i__;
   for (k = 2; k <= i__1; ++k) {
      k1 = k - 1;
      t = hh[k1 + i__ * 51 - 52];
      hh[k1 + i__ * 51 - 52] = c__[k1 - 1] * t + s[k1 - 1] * hh[k + i__ *51 - 52];
      hh[k + i__ * 51 - 52] = -s[k1 - 1] * t + c__[k1 - 1] * hh[k + i__ *51 - 52];
   }

L121:
   /* Computing 2nd power */
   d__1 = hh[i__ + i__ * 51 - 52];

   /* Computing 2nd power */
   d__2 = hh[i1 + i__ * 51 - 52];
   gam = sqrt(d__1 * d__1 + d__2 * d__2);

   /* if gamma is zero then any small value will do... */
   /* will affect only residual estimate */
   if (gam == 0.) {
      gam = epsmac;
   }

   /* get next plane rotation */
   c__[i__ - 1] = hh[i__ + i__ * 51 - 52] / gam;
   s[i__ - 1] = hh[i1 + i__ * 51 - 52] / gam;
   rs[i1 - 1] = -s[i__ - 1] * rs[i__ - 1];
   rs[i__ - 1] = c__[i__ - 1] * rs[i__ - 1];

   /* determine residual norm and test for convergence- */
   hh[i__ + i__ * 51 - 52] = c__[i__ - 1] * hh[i__ + i__ * 51 - 52]
                             + s[i__ - 1] * hh[i1 + i__ * 51 - 52];
   ro = (d__1 = rs[i1 - 1], std::abs(d__1));

   if (*iout>0)
     std::cout << "GMRES: iteration/residual: " << its << '/' << ro << std::endl;

   if (its<=1)
     goto L4;
   if (i__<*im && ro>eps1)
     goto L4;

   /* now compute solution. first solve upper triangular system. */
   rs[i__ - 1] /= hh[i__ + i__ * 51 - 52];
   i__1 = i__;
   for (ii = 2; ii <= i__1; ++ii) {
      k = i__ - ii + 1;
      k1 = k + 1;
      t = rs[k - 1];
      i__2 = i__;
      for (j = k1; j <= i__2; ++j) {
   t -= hh[k + j * 51 - 52] * rs[j - 1];
      }
      rs[k - 1] = t / hh[k + k * 51 - 52];
   }

   /* form linear combination of v(*,i)'s to get solution */
   t = rs[0];
   i__1 = *n;
   for (k = 1; k <= i__1; ++k) {
      rhs[k] = vv[k + vv_dim1] * t;
   }
   i__1 = i__;
   for (j = 2; j <= i__1; ++j) {
      t = rs[j - 1];
      i__2 = *n;
      for (k = 1; k <= i__2; ++k) {
   rhs[k] += t * vv[k + j * vv_dim1];
      }
   }

   /* call preconditioner. */
   lusol(n, &rhs[1], &rhs[1], &alu[1], &jlu[1], &ju[1]);
   i__1 = *n;
   for (k = 1; k <= i__1; ++k) {
      sol[k] += rhs[k];
   }

   /* restart outer loop  when necessary */
   if (ro <= eps1) {
      goto L990;
   }
   if (its >= *maxits) {
      goto L991;
   }

   /* else compute residual vector and continue.. */
   i__1 = i__;
   for (j = 1; j <= i__1; ++j) {
      jj = i1 - j + 1;
      rs[jj - 2] = -s[jj - 2] * rs[jj - 1];
      rs[jj - 1] = c__[jj - 2] * rs[jj - 1];
   }
   i__1 = i1;
   for (j = 1; j <= i__1; ++j) {
      t = rs[j - 1];
      if (j == 1) {
   t += -1.;
      }
      daxpy(n, &t, &vv[j * vv_dim1 + 1], &c__1, &vv[vv_offset], &c__1);
   }

   /* restart outer loop. */
    goto L20;

L990:
    *ierr = 0;

L991:
    *ierr = 1;

/* L999: */
    *ierr = -1;
}


}  // namespace lss
}  // namespace cf3
