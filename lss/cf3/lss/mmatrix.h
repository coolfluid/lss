#ifndef mmatrix_h
#define mmatrix_h

#include <vector>


namespace m {


/*
 * this is a struct describing a matrix; since it is based on parent class
 * template extension to simulate virtual functions, it is not a polymorphic
 * type and thus can't be used with factories.
 */

// description of a matrix (T: storage type, P: parent class)
template< typename T, class P >
struct mmatrix {
  // constructor
  mmatrix() : Nr(0), Nc(0), Nb(0), issparse(false) {}
  // indexing functions (absolute indexing)
  const T& operator()(const unsigned r, const unsigned c) const { return P::operator()(r,c); };
        T& operator()(const unsigned r, const unsigned c)       { return P::operator()(r,c); };
  // indexing functions (block indexing)
  const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return P::operator()(R,C,r,c); };
        T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return P::operator()(R,C,r,c); };
  // interfacing functions
  void reset(const T& v=T())                       { P::reset(v);     };
  void zerorow(const unsigned r)                   { P::zerorow(r);   };
  void zerorow(const unsigned R, const unsigned r) { P::zerorow(R,r); };
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Nr, unsigned _Nc, unsigned _Nb=1) { Nr=_Nr; Nc=_Nc; Nb=_Nb; }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {}
  // members
  unsigned Nr;  // number of (block) rows
  unsigned Nc;  // ... (block) columns
  unsigned Nb;  // ... block size
  bool issparse;
};


// implementation of a dense matrix, std::vector< std::vector< T > > based
template< typename T >
struct mmatrix_vv : mmatrix< T,mmatrix_vv< T > > {
  typedef mmatrix< T,mmatrix_vv< T > > P;
  // indexing functions (absolute indexing)
  const T& operator()(const unsigned r, const unsigned c) const { return a[r][c]; }
        T& operator()(const unsigned r, const unsigned c)       { return a[r][c]; }
  // indexing functions (block indexing)
  const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
        T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
  // interfacing functions
  void reset(const T& v=T())                       { a.assign(P::Nb*P::Nr,std::vector< T >(P::Nb*P::Nc,v)); };
  void zerorow(const unsigned r)                   { a[r].assign(P::Nb*P::Nc,T()); }
  void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); };
  // initialize method for dense variation
  void initialize(unsigned _Nr, unsigned _Nc, unsigned _Nb=1) {
    P::initialize(_Nr,_Nc,_Nb);
    reset(T());
  }
  // members
  std::vector< std::vector< T > > a;
};


// implementation of a dense matrix, T[][] based
template< typename T >
struct mmatrix_aa : mmatrix< T,mmatrix_aa< T > > {
  typedef mmatrix< T,mmatrix_aa< T > > P;
  // destructor
  ~mmatrix_aa() {
    if (P::Nr || P::Nc) {
      delete[] a[0];
      delete[] a;
    }
  }
  // indexing functions (absolute indexing)
  const T& operator()(const unsigned r, const unsigned c) const { return a[r][c]; }
        T& operator()(const unsigned r, const unsigned c)       { return a[r][c]; }
  // indexing functions (block indexing)
  const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
        T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
  // interfacing functions
  void reset(const T& v=T())                       { for (unsigned r=0; r<P::Nb*P::Nr; ++r) for (unsigned c=0; c<P::Nb*P::Nc; ++c) a[r][c] = v; };
  void zerorow(const unsigned r)                   { for (unsigned c=0; c<P::Nb*P::Nc; ++c) a[r][c] = T(); }
  void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); };
  // initialize method for dense variation
  void initialize(unsigned _Nr, unsigned _Nc, unsigned _Nb=1) {
    P::initialize(_Nr,_Nc,_Nb);
    a    = new T*[ P::Nb*P::Nr ];
    a[0] = new T [ P::Nb*P::Nr * P::Nb*P::Nc ];
    for (unsigned r=1; r<P::Nb*P::Nr; ++r)
      a[r] = a[r-1] + P::Nb*P::Nc;
    reset(T());
  }
  // members
  T **a;
};


// implementation of a sparse matrix, compressed sparse rows (CSR) format
// note: BASE={0,1}: {0,1}-based indexing (other values probably don't work)
template< typename T, int BASE >
struct mmatrix_csr : mmatrix< T,mmatrix_csr< T,BASE > > {
  typedef mmatrix< T,mmatrix_csr< T,BASE > > P;
  // cons/destructor
  mmatrix_csr() : P(), zero(T()), nnz(0), nnu(0) {
    P::issparse = true;
  }
  ~mmatrix_csr() {
    if (nnz) {
      delete[] a;
      delete[] ja;
      delete[] ia;
    }
  }
  // indexing functions (absolute indexing)
  const T& operator()(const unsigned r, const unsigned c) const { const int i=getindex(r,c); if (i<0) return zero; return a[i]; }
        T& operator()(const unsigned r, const unsigned c)       { const int i=getindex(r,c); if (i<0) return zero; return a[i]; }
  // indexing functions (block indexing)
  const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
        T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
  // interfacing functions
  void reset(const T& v=T())                       { for (int i=0; i<nnz; ++i) a[i] = v; };
  void zerorow(const unsigned r)                   { for (int k=ia[r]-BASE; k<ia[r+1]-BASE; ++k) a[k] = T(); }
  void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); };
  // initialize methods for base/sparse variation
  void initialize(unsigned _Nr, unsigned _Nc, unsigned _Nb=1) { P::initialize(_Nr,_Nc,_Nb); }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {

    // set row indices
    nnu = (int) P::Nb * (int) nz.size();
    ia = new int[nnu+1];
    ia[0] = BASE;
    int k = 0;
    for (int R=0; R<(int) nz.size(); ++R)
      for (int i=0; i<(int) P::Nb; ++i, ++k)
        ia[k+1] = ia[k] + (int) P::Nb * (int) nz[R].size();
    nnz = ia[nnu]-BASE;

    // set column indices
    ja = new int[nnz];
    for (unsigned R=0; R<(unsigned) nz.size(); ++R)
      for (int r=0; r<(int) P::Nb; ++r) {
        k = ia[R*P::Nb+r]-BASE;
        for (unsigned I=0; I<(unsigned) nz[R].size(); ++I)
          for (int c=0; c<(int) P::Nb; ++c)
            ja[k++] = (int) (P::Nb*nz[R][I]) + c + BASE;
      }

    // set entries
    a = new T[nnz];
    reset(T());
  }
  int getindex(const unsigned r, const unsigned c) const {
    for (int k=ia[r]-BASE; k<ia[r+1]-BASE; ++k)
      if (ja[k]-BASE==(int) c)
        return k;
    return -1;
  }
  // members
  T zero;
  int *ia;
  int *ja;
  T   *a;
  int nnz;
  int nnu;
};


// implementation of a sparse matrix, modified sparse row (MSR) format
template< typename T >
struct mmatrix_msr : mmatrix< T,mmatrix_msr< T > > {
  typedef mmatrix< T,mmatrix_msr< T > > P;
  // cons/destructor
  mmatrix_msr() : P(), zero(T()), nnu(0), nnz(0) {
    P::issparse = true;
  }
  ~mmatrix_msr() {
    if (nnz) {
      delete[] val;
      delete[] bindx;
    }
  }
  // indexing functions (absolute indexing)
  const T& operator()(const unsigned r, const unsigned c) const { const int i=getindex(r,c); if (i<0) return zero; return val[i]; }
        T& operator()(const unsigned r, const unsigned c)       { const int i=getindex(r,c); if (i<0) return zero; return val[i]; }
  // indexing functions (block indexing)
  const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
        T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
  // interfacing functions
  void reset(const T& v=T())                       { for (int i=0; i<bindx[nnu]; ++i) val[i] = v; };
  void zerorow(const unsigned r)                   { for (int k=bindx[r]; k<bindx[r+1]; ++k) val[k] = T(); val[r] = T(); }
  void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); };
  // initialize methods for base/sparse variation
  void initialize(unsigned _Nr, unsigned _Nc, unsigned _Nb=1) { P::initialize(_Nr,_Nc,_Nb); }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {

    // set number of rows/non-zero entries
    nnu = (int) P::Nb * (int) nz.size();
    nnz = 0;
    for (int i=0; i<(int) nz.size(); ++i)
      nnz += (int) nz[i].size();
    nnz *= (int) (P::Nb*P::Nb);

    // set row number of off-diagonal non-zeros and their addresses
    bindx = new int[nnz + 1];
    bindx[0] = nnu + 1;
    for (int R=0; R<(int) nz.size(); ++R)
      for (int r=0; r<(int) P::Nb; ++r) {
        const int i = (int) P::Nb * R + r;
        bindx[i+1] = bindx[i] + (int) P::Nb * (int) nz[R].size() - 1;
        int k = bindx[i];
        for (unsigned I=0; I<(unsigned) nz[R].size(); ++I)
          for (int c=0; c<(int) P::Nb; ++c)
            if (R!=(int) nz[R][I] || r!=c)
              bindx[k++] = (int) (P::Nb*nz[R][I]) + c;
      }

    // set entries
    val = new T[bindx[nnu]];  // or nnz+1
    reset(T());
  }
  // utilities
  int getindex(const unsigned r, const unsigned c) const {
    if (r==c)
      return r;
    for (int k=bindx[r]; k<bindx[r+1]; ++k)
      if (bindx[k]==(int) c)
        return k;
    return -1;
  }
  // members
  T   zero;
  T   *val;
  int *bindx;
  int nnu;
  int nnz;
};


// implementation of a sparse matrix, variable block compressed sparse row (VBR) format
template< typename T >
struct mmatrix_vbr : mmatrix< T,mmatrix_vbr< T > > {
  typedef mmatrix< T,mmatrix_vbr< T > > P;
  // cons/destructor and indexing operators
  mmatrix_vbr() : P(), zero(T()), bnnu(0), bnnz(0) {
    P::issparse = true;
  }
  ~mmatrix_vbr() {
    if (bnnz) {
      delete[] val;
      delete[] bpntr;
      free(cpntr);  // (set by AZ_transform)
      delete[] rpntr;
      delete[] bindx;
      delete[] indx;
    }
  }
  // indexing functions (absolute indexing)
  const T& operator()(const unsigned r, const unsigned c) const { const int b=(int) P::Nb; return operator()(r/b,c/b,r%b,c%b); }
        T& operator()(const unsigned r, const unsigned c)       { const int b=(int) P::Nb; return operator()(r/b,c/b,r%b,c%b); }
  // indexing functions (block indexing)
  const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { const int i=getindex(R,C,r,c); if (i<0) return zero; return val[i]; }
        T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { const int i=getindex(R,C,r,c); if (i<0) return zero; return val[i]; }
  // interfacing functions
  void reset(const T& v=T())                       { for (int i=0; i<indx[bpntr[bnnu]]; ++i) val[i] = v; };
  void zerorow(const unsigned r)                   { const int b=(int) P::Nb; zerorow(r/b,r%b); }
  void zerorow(const unsigned R, const unsigned r) {
    const int b = (int) P::Nb;
    const int i = getindex(R,R,r,0);
    for (int j=0; j<b*b*(bpntr[R+1]-bpntr[R]); j+=b)
      val[i+j] = T();
  }
  // initialize methods for base/sparse variation
  void initialize(unsigned _Nr, unsigned _Nc, unsigned _Nb=1) { P::initialize(_Nr,_Nc,_Nb); }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {

    // count rows and non-zero (blocks)
    bnnu = (int) nz.size();
    bnnz = 0;
    for (int i=0; i<bnnu; ++i)
      bnnz += (int) nz[i].size();

    // allocate data structure (cpntr==rpntr in a struct. symmetric matrix)
    indx  = new int[bnnz+1];
    bindx = new int[bnnz];
    rpntr = new int[bnnu+1];
    cpntr = NULL;  // (set by AZ_transform)
    bpntr = new int[bnnu+1];

    // set it up (diagonal blocks come first)
    rpntr[0] = bpntr[0] = bindx[0] = indx[0] = 0;
    for (int i=0; i<bnnu; ++i) {
      rpntr[i+1] = rpntr[i] + (int) P::Nb;
      bpntr[i+1] = bpntr[i] + (int) nz[i].size();
      int k = bpntr[i];
      bindx[k++] = i;
      for (int j=0; j<(int) nz[i].size(); ++j)
        if ((int) nz[i][j]!=i)
          bindx[k++] = (int) nz[i][j];
      k = bpntr[i];
      for (int j=0; j<(int) nz[i].size(); ++j, ++k)
        indx[k+1] = indx[k] + (int) (P::Nb*P::Nb);
    }

    // set entries
    val = new double[indx[bpntr[bnnu]]];  // or bnnz*P::Nb*P::Nb
    reset(T());
  }
  // utilities
  int getindex(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const {
    for (int i=bpntr[R]; i<bpntr[R+1]; ++i)
      if (bindx[i]==(int) C)
        return indx[i] + P::Nb*c + r;
    return -1;
  }
  // members
  T   zero;
  T   *val;
  int *indx, *bindx, *rpntr, *cpntr, *bpntr;
  int bnnu;
  int bnnz;
};


}  // namespace m


#endif

