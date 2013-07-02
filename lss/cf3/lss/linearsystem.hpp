// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_linearsystem_h
#define cf3_lss_linearsystem_h


#include <vector>
#include "boost/lexical_cast.hpp"
#include "common/Signal.hpp"
#include "common/Action.hpp"


namespace cf3 {
namespace lss {


/**
 * description of a matrix row/column indexer
 * the idea is to return a single size_t index to a (very long and boring)
 * vector of entries that populate a sparse, dense, or block-addressable matrix.
 */
struct index {
  virtual void output(std::ostream& o) = 0;
#if 0
unsigned Ne;  // number of (block) equations
unsigned Nv;  // ... (block) variables
unsigned Nb;  // ... block size
virtual const T& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const = 0;
virtual       T& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       = 0;
const T& X(const unsigned R, const unsigned r) const { return X(R*Nb+r); }
T& X(const unsigned R, const unsigned r)       { return X(R*Nb+r); }
const T& B(const unsigned C, const unsigned c) const { return B(C*Nb+c); }
T& B(const unsigned C, const unsigned c)       { return B(C*Nb+c); }
virtual void resize(const std::vector< std::vector< unsigned > >& nz) = 0;
virtual void zerorow(const unsigned R, const unsigned r) { zerorow(R*Nb+r); }
void output(std::ostream& o) {     o << "m::linearsystem::B(Ne:" << Nv << ",Nb:" << Nb << "):" << std::endl;     }
// indexing functions (block indexing)
const double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return m_A(R*Nb+r,C*Nb+c); }
double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return m_A(R*Nb+r,C*Nb+c); }
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return P::operator()(R,C,r,c); };
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return P::operator()(R,C,r,c); };

void zerorow(const unsigned R, const unsigned r) { P::zerorow(R,r); }
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); }
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }
void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); }

void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); }
void zerorow(const unsigned R, const unsigned r) { zerorow(P::Nb*R+r); }
void zerorow(const unsigned R, const unsigned r) {
const int b = (int) P::Nb;
const int i = getindex(R,R,r,0);
val[i+j] = T();
}
// indexing functions (block indexing)
const T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return operator()(P::Nb*R+r,P::Nb*C+c); }
T& operator()(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return operator()(P::Nb*R+r,P::Nb*C+c); }


  BOOST_PARAMETER_FUNCTION(
    (void),                // 1. parenthesized return type
    depth_first_search,    // 2. name of the function template
    tag,                   // 3. namespace of tag types

    (required
      (graph, *) ) // 4. one required parameter, and

    (optional              //    four optional parameters, with defaults
      (visitor,           *, boost::dfs_visitor<>())
      (root_vertex,       *, *vertices(graph).first)
      (index_map,         *, get(boost::vertex_index,graph))
      (in_out(color_map), *, default_color_map(num_vertices(graph), index_map)) )
    )
    {
      // ... body of function goes here...
      // use graph, visitor, index_map, and color_map
    }

#endif
};


/**
 * description of a matrix
 * only describes an interface, storage is provided by the derived structs.
 * (T: storage type, P: parent class)
 * it is based on parent class template extension (CRTP) to simulate virtual
 * methods (not a true polymorphic type) so it can't be used with factories.
 */
template< typename T, class P >
struct matrix
{
  // constructor
  matrix() : Nr(0), Nc(0) {}

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return P::operator()(r,c); }
        T& operator()(const size_t r, const size_t c)       { return P::operator()(r,c); }
  const T& operator()(const size_t i) const { return P::operator()(i); }
        T& operator()(const size_t i)       { return P::operator()(i); }

  size_t size()                const { return Nr*Nc; }
  size_t size(const size_t& d) const { return (d==0? Nr : (d==1? Nc : 0)); }

  matrix& resize  (size_t r, size_t c, const T& v=T()) { Nr=r; Nc=c; return *this; }
  matrix& zerorow (const size_t r)     { return P::zerorow(r); }
  matrix& assign  (const T& v=T())     { return P::assign(v); }
  matrix& clear() { Nr = Nc = 0; }

  void output(std::ostream& out) const { P::output(out); }
  matrix& operator=(const T& v) { return P::assign(v); }

  // members
  size_t Nr;  // number of rows
  size_t Nc;  // ... columns
};


/**
 * dense matrix implementation, std::vector< std::vector< T > > based
 */
template< typename T >
struct matrix_dense_vv :
    matrix< T,matrix_dense_vv< T > > {
  typedef matrix< T,matrix_dense_vv< T > > P;

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return a[r][c]; }
        T& operator()(const size_t r, const size_t c)       { return a[r][c]; }
  const T& operator()(const size_t i) const { const size_t b(P::size(0)); return a[i/b][i%b]; }
        T& operator()(const size_t i)       { const size_t b(P::size(0)); return a[i/b][i%b]; }

  matrix_dense_vv& resize  (const size_t r, const size_t c, const T& v=T()) { P::resize(r,c); return assign(v); }
  matrix_dense_vv& zerorow (const size_t r) { if (P::size()) a[r].assign(P::Nc,T());                    return *this; }
  matrix_dense_vv& assign  (const T& v=T()) { if (P::size()) a.assign(P::Nr,std::vector< T >(P::Nc,v)); return *this; }
  matrix_dense_vv& clear() { P::clear(); a.clear(); return *this; }

  void output(std::ostream& out) const {
    out << "[ ";
    BOOST_FOREACH(const std::vector< T >& r, a) {
      std::copy(r.begin(),r.end(), std::ostream_iterator< T >(out,", "));
      out << "\n  ";
    }
    out << ']';
  }

  // members
  std::vector< std::vector< T > > a;
};


/**
 * dense matrix implementation, T[][] array based
 */
template< typename T >
struct matrix_dense_aa :
    matrix< T,matrix_dense_aa< T > > {
  typedef matrix< T,matrix_dense_aa< T > > P;

  // destructor
  ~matrix_dense_aa() { matrix_dense_aa::clear(); }

  // matrix interfacing
  const T& operator()(const size_t r, const size_t c) const { return a[r][c]; }
        T& operator()(const size_t r, const size_t c)       { return a[r][c]; }
  const T& operator()(const size_t i) const { const size_t b(P::size(1)); return a[i/b][i%b]; }
        T& operator()(const size_t i)       { const size_t b(P::size(1)); return a[i/b][i%b]; }

  matrix_dense_aa& resize(const size_t r, const size_t c, const T& v=T()) {
    if (r && c) {
      P::resize(r,c);
      a    = new T*[ P::Nr ];
      a[0] = new T [ P::Nr * P::Nc ];
      for (size_t r=1; r<P::Nr; ++r)
        a[r] = a[r-1] + P::Nc;
      assign(v);
    }
  }

  matrix_dense_aa& zerorow(const size_t r) {
    for (size_t c=0; c<P::Nc; ++c)
      a[r][c] = T();
    return *this;
  }

  matrix_dense_aa& assign(const T& v=T())  {
    for (size_t r=0; r<P::Nr; ++r)
      for (size_t c=0; c<P::Nc; ++c)
        a[r][c] = v;
    return *this;
  }

  matrix_dense_aa& clear() {
    if (P::size()) {
      delete[] a[0];
      delete[] a;
    }
    P::clear();
  }

  void output(std::ostream& out) const {
    out << "[ ";
    for (size_t r=0; r<P::Nr; ++r) {
      std::copy(&(a[r][0]),&(a[r][0])+P::Nc, std::ostream_iterator< T >(out,", "));
      out << "\n  ";
    }
    out << ']';
  }

  // members
  T **a;
};


/**
 * sparse matrix implementation, in compressed sparse rows (CSR) format
 * note: BASE={0,1}: {0,1}-based indexing (other values probably don't work)
 */
template< typename T, int BASE >
struct matrix_sparse_csr : matrix< T,matrix_sparse_csr< T,BASE > > {
  typedef matrix< T,matrix_sparse_csr< T,BASE > > P;

  // cons/destructor
  matrix_sparse_csr() : P(), zero(T()), nnz(0), nnu(0) {}
  ~matrix_sparse_csr() {
    if (nnz) {
      delete[] a;
      delete[] ja;
      delete[] ia;
    }
  }

  // interfacing
  const T& operator()(const size_t r, const size_t c) const { const int i=getindex(r,c); if (i<0) return zero; return a[i]; }
        T& operator()(const size_t r, const size_t c)       { const int i=getindex(r,c); if (i<0) return zero; return a[i]; }

  void clear(const T& v=T())   { for (int i=0; i<nnz; ++i) a[i] = v; }
  void zerorow(const size_t r) { for (int k=ia[r]-BASE; k<ia[r+1]-BASE; ++k) a[k] = T(); }

  void resize(size_t r, size_t c) { P::resize(r,c); }
  void resize(const std::vector< std::vector< size_t > >& nz) {

    // set row indices
    nnu = (int) nz.size();
    ia = new int[nnu+1];
    ia[0] = BASE;
    int k = 0;
    for (size_t r=0; r<nz.size(); ++r)
      ia[k+1] = ia[k] + (int) nz[r].size();
    nnz = ia[nnu]-BASE;

    // set column indices
    ja = new int[nnz];
    for (size_t r=0; r<nz.size(); ++r) {
      k = ia[r]-BASE;
      for (size_t i=0; i<nz[r].size(); ++i)
        ja[k++] = (int) (nz[r][i]) + BASE;
    }

    // set entries
    a = new T[nnz];
    clear();
  }
  int getindex(const size_t r, const size_t c) const {
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


// NOTE: missing Aztec MSR and VBR (maybe this will never happen)


/**
 * description of a linear system
 * solution and right-hand side vectors are included, matrix should be
 * included (by aggregation) in derived linearsystem classes
 */
template< typename T > class linearsystem;
template< typename T > std::ostream& operator<< (std::ostream&, const linearsystem< T >&);
template< typename T >
class linearsystem :
  public common::Action
{
  typedef linearsystem< T > linearsystem_t;

 public:
  // framework interfacing
  linearsystem(const std::string& name) :
    common::Action(name),
    m_zero(T())
  {

    // set component options level
    mark_basic();

    // configure signals
    regist_signal("resize")
      .connect(   boost::bind( &linearsystem::signal_resize,    this, _1 ))
      .signature( boost::bind( &linearsystem::signature_ask_rc, this, _1 ));
    regist_signal("zerorow")
      .connect(   boost::bind( &linearsystem::signal_zerorow,  this, _1 ))
      .signature( boost::bind( &linearsystem::signature_ask_r, this, _1 ));
    regist_signal("clear") .connect( boost::bind( &linearsystem::signal_clear,  this ));
    regist_signal("solve") .connect( boost::bind( &linearsystem::signal_solve,  this ));
    regist_signal("output").connect( boost::bind( &linearsystem::signal_output, this ));

    // configure options
    options().add("A",std::vector< T >()).mark_basic().link_to(&m_swap)
      .attach_trigger(boost::bind( &linearsystem::trigger_A, this ));
    options().add("b",std::vector< T >()).mark_basic().link_to(&m_swap)
      .attach_trigger(boost::bind( &linearsystem::trigger_b, this ));
    options().add("x",std::vector< T >()).mark_basic().link_to(&m_swap)
      .attach_trigger(boost::bind( &linearsystem::trigger_x, this ));

  }
  virtual ~linearsystem() {}
  void execute() { solve(); }

 private:
  // framework scripting
  void signal_resize  (common::SignalArgs& args) { common::XML::SignalOptions opts(args); resize((size_t) opts.value< unsigned int >("r"), (size_t) opts.value< unsigned int >("c"), T()); }
  void signal_zerorow (common::SignalArgs& args) { common::XML::SignalOptions opts(args); zerorow(opts.value< unsigned int >("r")); }
  void signal_clear  () { clear(); }
  void signal_solve  () { solve(); }
  void signal_output () { operator<<(std::cout,*this); }

  void signature_ask_rc    (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("r",(unsigned int) size_t()); opts.add("c",(unsigned int) size_t());  }
  void signature_ask_r     (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("r",(unsigned int) size_t()); }
  void signature_ask_value (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("value",T());  }

  virtual void trigger_A() {
    if (m_swap.size()!=size())
      throw common::BadValue(FromHere(), "linearsystem is not of the same size as given matrix ("
        + boost::lexical_cast< std::string >(size(0)) + "*"
        + boost::lexical_cast< std::string >(size(1)) + "!="
        + boost::lexical_cast< std::string >(m_swap.size()) + ")");
    for (size_t r=0; r<size(0); ++r)
      for (size_t c=0; c<size(1); ++c)
        A(r,c) = m_swap[r*size(1)+c];
    m_swap.clear();
  }
  void trigger_b() { swap_lss_vector(m_b); }
  void trigger_x() { swap_lss_vector(m_x); }
  void swap_lss_vector (std::vector< T >& v) {
    if (m_swap.size()!=v.size())
      throw common::BadValue(FromHere(), "linearsystem is not of the same size as given vector ("
        + boost::lexical_cast< std::string >(v.size()) + "!="
        + boost::lexical_cast< std::string >(m_swap.size()) + ")");
    m_swap.swap(v);
    m_swap.clear();
  }

 public:
  // linear system solver addressing
  virtual const T& A(const size_t r, const size_t c) const = 0;
  virtual       T& A(const size_t r, const size_t c)       = 0;
  virtual const T& A(const size_t i) const = 0;
  virtual       T& A(const size_t i)       = 0;
  const T& b(const size_t r) const { return m_b[r]; }
        T& b(const size_t r)       { return m_b[r]; }
  const T& x(const size_t r) const { return m_x[r]; }
        T& x(const size_t r)       { return m_x[r]; }

 public:
  // linear system solver interfacing
  virtual size_t size()             const = 0;
  virtual size_t size(const size_t) const = 0;
  virtual linearsystem& resize (const size_t Nequations, const size_t Nvariables, const T&) = 0;
  virtual linearsystem& zerorow(const size_t) = 0;
  virtual linearsystem& clear() = 0;
  virtual linearsystem& solve() = 0;

 private:
  // output
  virtual void output_A(std::ostream& out) const { out << "[ (unavailable) ]"; }
  virtual void output_b(std::ostream& out) const { out << "[ "; std::copy(m_b.begin(),m_b.end(),std::ostream_iterator< T >(out,", ")); out << ']'; }
  virtual void output_x(std::ostream& out) const { out << "[ "; std::copy(m_x.begin(),m_x.end(),std::ostream_iterator< T >(out,", ")); out << ']'; }
  friend std::ostream& operator<< < T > (std::ostream&, const linearsystem&);

 protected:
  // variables
  std::vector< T > m_b;
  std::vector< T > m_x;
  std::vector< T > m_swap;
  T m_zero;
};


template< typename T >
std::ostream& operator<<(std::ostream& out, const linearsystem< T >& lss) {
  out << "linearsystem::A("
    << "Nequations=" << lss.size(0) << ","
    << "Nvariables=" << lss.size(1) << ","
    << "Nentries="   << lss.size()  << "):\n";
  lss.output_A(out);
  out << '\n'
      << "linearsystem::b(Nequations=" << lss.m_b.size() << "):\n";
  lss.output_b(out);
  out << '\n'
      << "linearsystem::x(Nvariables=" << lss.m_x.size() << "):\n";
  lss.output_x(out);
  return out << '\n';
}


}  // namespace lss
}  // namespace cf3


#endif

