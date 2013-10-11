// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_linearsystem_h
#define cf3_lss_linearsystem_h


#include <vector>
#include <limits>
//#include <cstdarg>  // NOTE: yes, here because it's a "really ugly" header
//#include "boost/lexical_cast.hpp"
//#include <boost/foreach.hpp>


#include "common/Signal.hpp"
#include "common/Action.hpp"


#include "lss_index.hpp"
#include "lss_matrix.hpp"


namespace cf3 {
namespace lss {


// helper forward declarations
template< typename T, typename MATRIX > class linearsystem;
template< typename T, typename MATRIX > std::ostream& operator<< (std::ostream&, const linearsystem< T, MATRIX >&);


/**
 * @brief Description of a linear system (suitable for sparse matrix solvers.)
 * solution and right-hand side vectors are STL containers.
 */
template< typename T, typename MATRIX >
class linearsystem :
  public common::Action
{
  typedef linearsystem< T, MATRIX > linearsystem_t;


  /* ------------------------------------------------------------------------ *
   * framework interfacing                                                    *
   * ------------------------------------------------------------------------ */
 public:

  linearsystem(const std::string& name) :
    common::Action(name),
    m_dummy_value(std::numeric_limits< T >::quiet_NaN())
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


  /* ------------------------------------------------------------------------ *
   * framework scripting                                                      *
   * ------------------------------------------------------------------------ */
 private:

  void signal_resize  (common::SignalArgs& args) { common::XML::SignalOptions opts(args); this->resize((size_t) opts.value< unsigned int >("r"), (size_t) opts.value< unsigned int >("c"), T()); }
  void signal_zerorow (common::SignalArgs& args) { common::XML::SignalOptions opts(args); this->zerorow(opts.value< unsigned int >("r")); }
  void signal_clear  () { this->clear(); }
  void signal_solve  () { this->solve(); }
  void signal_output () { operator<<(std::cout,*this); }

  void signature_ask_rc    (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("r",(unsigned int) size_t()); opts.add("c",(unsigned int) size_t());  }
  void signature_ask_r     (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("r",(unsigned int) size_t()); }
  void signature_ask_value (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("value",T());  }

  void trigger_b() { swap_lss_vector(m_b); }
  void trigger_x() { swap_lss_vector(m_x); }
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

  void swap_lss_vector (std::vector< T >& v) {
    if (m_swap.size()!=v.size())
      throw common::BadValue(FromHere(), "linearsystem is not of the same size as given vector ("
        + boost::lexical_cast< std::string >(v.size()) + "!="
        + boost::lexical_cast< std::string >(m_swap.size()) + ")");
    m_swap.swap(v);
    m_swap.clear();
  }


  /* ------------------------------------------------------------------------ *
   * linear system addressing                                                 *
   * ------------------------------------------------------------------------ */
 public:

  const T& b(const size_t r) const { return m_b[r]; }
        T& b(const size_t r)       { return m_b[r]; }
  const T& x(const size_t r) const { return m_x[r]; }
        T& x(const size_t r)       { return m_x[r]; }
  const T& A(const size_t r, const size_t c) const { return m_dummy_value; /*m_A(index_dereference(r,c));*/ }
        T& A(const size_t r, const size_t c)       { return m_dummy_value; /*m_A(index_dereference(r,c));*/ }


  /* ------------------------------------------------------------------------ *
   * linear system interfacing                                                *
   * ------------------------------------------------------------------------ */
 public:

  virtual linearsystem& solve() = 0;

  size_t size()               const { return m_A.size(); }
  size_t size(const size_t d) const { return m_A.size(d); }

  linearsystem& clear() {
    m_A.clear();
    m_b.clear();
    m_x.clear();
    return *this;
  }

  linearsystem& resize(const size_t Nequations, const size_t Nvariables, const T& v=T()) {
    m_A.resize(Nequations,Nvariables,v);
    m_b.assign(Nequations,v);
    m_x.assign(Nvariables,v);
    return *this;
  }

  linearsystem& zerorow(const size_t r) {
    m_A.zerorow(r);
    m_b[r] = T();
    return *this;
  }

 private:

  virtual void output_A(std::ostream& out) const { out << "[ (unavailable) ]"; }
  virtual void output_b(std::ostream& out) const { out << "[ "; std::copy(m_b.begin(),m_b.end(),std::ostream_iterator< T >(out,", ")); out << ']'; }
  virtual void output_x(std::ostream& out) const { out << "[ "; std::copy(m_x.begin(),m_x.end(),std::ostream_iterator< T >(out,", ")); out << ']'; }
  friend std::ostream& operator<< < T, MATRIX >(std::ostream&, const linearsystem&);


  /* ------------------------------------------------------------------------ *
   * indexing manipulation                                                    *
   * - implemented in stack style, forcing indexing hierarchy on client side  *
   *   (simulate a FILO container, std::stack is a piece of crap, really!)    *
   * - methods don't return anything, to force deliberate statements          *
   * ------------------------------------------------------------------------ */
 public:

  /// Registration of the indexing hierarchy
  void index(const size_t n, lss_index::index_t idx, ...) {
    if (n==1) m_indexer.push_back(idx);
    if (n<2) return;
    va_list ap;
    va_start(ap,idx);
    for (size_t i=1; i<n; ++i)
      m_indexer.push_back( va_arg(ap,lss_index::index_t) );
    va_end(ap);
  }

#if 0
  /// Indexing insertion
  void index_push(lss_index::index_t& _idx) { m_indexer.push_back(_idx); }

  /// Indexing removal
  void index_pop() { if (!m_indexer.empty()) m_indexer.pop_back(); }

  /// Indexing configuration (top-most only, to force hierarchy in client side)
  void index(size_t& _i, size_t& _j=size_t()) { m_indexer.back().index(_i,_j); }
#endif

 private:  // internal ugly stuff

  /// Indexing dereference
  const lss_index::ij index_dereference(const size_t r, const size_t c) {
    lss_index::ij ij(r,c);
    BOOST_REVERSE_FOREACH(const lss_index::index_t& idx, m_indexer)
      idx.index(ij);
    return ij;
  }


  /* ------------------------------------------------------------------------ *
   * internal variables                                                       *
   * ------------------------------------------------------------------------ */
 protected:

  /// Linear system matrix, right and left-hand sides vectors, and dummy value
  MATRIX m_A;
  std::vector< T > m_b;
  std::vector< T > m_x;
  T m_dummy_value;  // should keep NaN throughout

  /// Indexing stack
  std::vector< lss_index::index_t > m_indexer;

 private:

  /// Scripting temporary storage (to be swapped with internal containers)
  std::vector< T > m_swap;

};


/**
 * printing of a linear system to a given stream
 */
template< typename T, typename MATRIX >
std::ostream& operator<<(std::ostream& out, const linearsystem< T, MATRIX >& lss) {
  const bool
    output_A(false),
    output_b(false),
    output_x(true);
  if (output_A) {
    out << "linearsystem::A("
        << "Nequations=" << lss.size(0) << ","
        << "Nvariables=" << lss.size(1) << ","
        << "Nentries="   << lss.size()  << "):\n";
    lss.output_A(out);
    out << '\n';
  }
  if (output_b) {
    out << "linearsystem::b(Nequations=" << lss.m_b.size() << "):\n";
    lss.output_b(out);
    out << '\n';
  }
  if (output_x) {
    out << "linearsystem::x(Nvariables=" << lss.m_x.size() << "):\n";
    lss.output_x(out);
    out << '\n';
  }
  return out;
}


}  // namespace lss
}  // namespace cf3


#endif

