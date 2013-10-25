// Copyright (C) 2013 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_linearsystem_hpp
#define cf3_lss_linearsystem_hpp


#include <limits>
#include <sstream>
#include <stdexcept>

#include "common/Signal.hpp"
#include "common/Action.hpp"

#include "detail/index.hpp"
#include "detail/matrix.hpp"


namespace cf3 {
namespace lss {


/* -- linear system --------------------------------------------------------- */

/**
 * @brief Description of a linear system (suitable for sparse matrix solvers.)
 * solution and right-hand side vectors are STL containers.
 * @note the class name forces the distinction of (this) acessible plugin and
 * the detail::linearsystem building class, to avoid 'using namespace's.
 */
class linearsystem : public common::Action
{
  // -- Construction and destruction
 public:

  /// Construct the linear system
  linearsystem(const std::string& name) :
    common::Action(name),
    m_dummy_value(std::numeric_limits< double >::quiet_NaN())
  {
    // set component options level
    mark_basic();

#if 0
    // configure signals
    regist_signal("initialize").connect(   boost::bind( &linearsystem::signal_initialize, this, _1 )).signature( boost::bind( &linearsystem::signature_ask_rc, this, _1 ));
    regist_signal("zerorow")   .connect(   boost::bind( &linearsystem::signal_zerorow,    this, _1 )).signature( boost::bind( &linearsystem::signature_ask_r,  this, _1 ));
    regist_signal("clear") .connect( boost::bind( &linearsystem::signal_clear,  this ));
    regist_signal("solve") .connect( boost::bind( &linearsystem::signal_solve,  this ));
    regist_signal("output").connect( boost::bind( &linearsystem::signal_output, this ));

    // configure options
    options().add("A",std::vector< T >()).mark_basic().link_to(&m_swap).attach_trigger(boost::bind( &linearsystem::trigger_A, this ));
    options().add("b",std::vector< T >()).mark_basic().link_to(&m_swap).attach_trigger(boost::bind( &linearsystem::trigger_b, this ));
    options().add("x",std::vector< T >()).mark_basic().link_to(&m_swap).attach_trigger(boost::bind( &linearsystem::trigger_x, this ));
#endif
  }

  /// Destruct the linear system
  virtual ~linearsystem() {}


  // -- Framework scripting
 private:

#if 0
  void signal_initialize (common::SignalArgs& args) { common::XML::SignalOptions opts(args); this->initialize((size_t) opts.value< unsigned int >("r"), (size_t) opts.value< unsigned int >("c"), T()); }
  void signal_zerorow    (common::SignalArgs& args) { common::XML::SignalOptions opts(args); this->zerorow(opts.value< unsigned int >("r")); }
  void signal_clear  () { this->clear(); }
  void signal_solve  () { this->solve(); }
  void signal_output () { operator<<(std::cout,*this); }

  void signature_ask_rc    (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("r",(unsigned int) size_t()); opts.add("c",(unsigned int) size_t());  }
  void signature_ask_r     (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("r",(unsigned int) size_t()); }
  void signature_ask_value (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("value",T());  }

  void trigger_b() { swap_lss_vector(m_b); }
  void trigger_x() { swap_lss_vector(m_x); }
  void trigger_A() {
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
#endif


  // -- Interfacing (pure)
 public:

  /// Initialize the linear system (resizing consistently)
  virtual linearsystem& initialize(
      const size_t& _size_i,
      const size_t& _size_j,
      const size_t& _size_k=1,
      const double& _value=double()) = 0;

  /// Initialize linear system from file(s)
  virtual linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) = 0;

  /// Initialize linear system from vectors of values (lists, in the right context)
  virtual linearsystem& initialize(
      const std::vector< double >& vA,
      const std::vector< double >& vb=std::vector< double >(),
      const std::vector< double >& vx=std::vector< double >()) = 0;

  /// Linear system solving
  virtual linearsystem& solve() = 0;

  /// Linear system solving, aliased from execute
  void execute() { solve(); }

  // -- Storage
 protected:

  double m_dummy_value;          // should keep NaN throughout
  std::vector< double > m_swap;  // scripting temporary storage (to swap with internal storage)

};



#if 0
#include <vector>
//#include "boost/lexical_cast.hpp"
//#include <boost/foreach.hpp>
#include "lss_index.hpp"
#include "lss_matrix.hpp"


// helper forward declarations
template< typename T, typename MATRIX > class linearsystem;
template< typename T, typename MATRIX > std::ostream& operator<< (std::ostream&, const linearsystem< T, MATRIX >&);


class xpto
{
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

  size_t size()               const { return m_A.size(); }
  size_t size(const size_t d) const { return m_A.size(d); }

 private:

  virtual void output_A(std::ostream& out) const { out << "[ (unavailable) ]"; }
  virtual void output_b(std::ostream& out) const { out << "[ "; std::copy(m_b.begin(),m_b.end(),std::ostream_iterator< T >(out,", ")); out << ']'; }
  virtual void output_x(std::ostream& out) const { out << "[ "; std::copy(m_x.begin(),m_x.end(),std::ostream_iterator< T >(out,", ")); out << ']'; }
  friend std::ostream& operator<< < T, MATRIX >(std::ostream&, const linearsystem&);

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

#endif


}  // namespace lss
}  // namespace cf3


#endif
