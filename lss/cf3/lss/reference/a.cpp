
#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <vector>


#include "linearsystem.hpp"


/// Prototypes for single and double precision (as per Intel MKL)
extern "C"
{
  void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
  void sgesv_(int* n, int* nrhs, float*  a, int* lda, int* ipiv, float*  b, int* ldb, int* info);
}


/// Lapack linear system solver, single and double precision
template< typename T >
struct lss_lapack : public linearsystem< T,
    lss_matrix::dense_matrix_v< T, lss_matrix::column_oriented >,
    lss_matrix::dense_matrix_v< T, lss_matrix::column_oriented > >
{
  typedef linearsystem< T,
    lss_matrix::dense_matrix_v< T, lss_matrix::column_oriented >,
    lss_matrix::dense_matrix_v< T, lss_matrix::column_oriented > >
  lss_lapack_base_t;

  lss_lapack(const size_t& _size_i=size_t(),
             const size_t& _size_j=size_t(),
             const size_t& _size_k=1,
             const T& _value=T() ) :
    lss_lapack_base_t(_size_i,_size_j,_size_k,_value) {}

  bool solve() {
    int n    = static_cast< int >(this->size(0));
    int nrhs = static_cast< int >(this->size(2));
    int err  = 0;
    std::vector< int > ipiv(n);

    if (!this->size().is_square_size()) { err = -17; }
    else if (typeid(T)==typeid(double)) { dgesv_( &n, &nrhs, (double*) &this->m_A.a[0], &n, &ipiv[0], (double*) &this->m_b.a[0], &n, &err ); }
    else if (typeid(T)==typeid(float))  { sgesv_( &n, &nrhs, (float*)  &this->m_A.a[0], &n, &ipiv[0], (float*)  &this->m_b.a[0], &n, &err ); }
    else { err = -42; }

    std::ostringstream msg;
    err==-17? msg << "lss_lapack: system matrix size must be square" :
    err==-42? msg << "lss_lapack: precision not implemented: " << typeid(T).name() :
    err<0?    msg << "lss_lapack: invalid " << err << "'th argument to dgesv_()" :
    err>0?    msg << "lss_lapack: triangular factor matrix U(" << err << ',' << err << ") is zero, so A is singular (not invertible)" :
              msg;
    if (err)
      throw std::runtime_error(msg.str());

    // solution is in b, so swap with x (with A square, size b = size x)
    this->m_b.swap(this->m_x);
    return !err;
  }
};


int main()
{
  /*
    lss_matrix::sparse_matrix_csr< float, 0 > A;
    lss_matrix::dense_matrix_vv< float > A;
    lss_matrix::dense_matrix_v< double, lss_matrix::column_oriented > A;
    A.initialize("/home/pmaciel/Desktop/m/matrixmarket/test1.csr");
    A.initialize("/home/pmaciel/Desktop/m/matrixmarket/test.csr");
    A.initialize("/home/pmaciel/Desktop/m/matrixmarket/test_dense.mtx");
    A.initialize("/home/pmaciel/Desktop/m/matrixmarket/e05r0100.mtx");
   */

  /*
    float vfA[] = {
       6.80, -6.05, -0.45,  8.32, -9.67,
      -2.11, -3.30,  2.58,  2.71, -5.14,
       5.66,  5.36, -2.70,  4.35, -7.26,
       5.97, -4.44,  0.27, -7.17,  6.08,
       8.23,  1.08,  9.04,  2.14, -6.87, 0. };
    float vfb[] = {
       4.02, -1.56,  9.81,
       6.19,  4.00, -4.09,
      -8.22, -8.67, -4.57,
      -7.57,  1.75, -8.61,
      -3.03,  2.86,  8.99, 0. };

    lss_lapack< float > q;
    q.resize(5,5,3)
     .assign( std::vector< float >(&vfA[0], &vfA[25]), std::vector< float >(&vfb[0], &vfb[15]))
     .solve();
    cout << q << endl;
   */


  lss_lapack< double > r;
  r.initialize("/home/pmaciel/Desktop/m/testmatrices_sparskit_drivcav/e05r0500.mtx",
               "/home/pmaciel/Desktop/m/testmatrices_sparskit_drivcav/e05r0500_rhs1.mtx")
   .solve();
  std::cout << r << std::endl;


  return 0;
}

