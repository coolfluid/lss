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

#include "detail/utilities.hpp"
#include "detail/index.hpp"
#include "detail/matrix.hpp"


namespace cf3 {
namespace lss {


/* -- linear system --------------------------------------------------------- */


// helper forward declarations
template< typename T, typename MATRIX, typename VECTOR > class linearsystem;
template< typename T, typename MATRIX, typename VECTOR > std::ostream& operator<< (std::ostream&, const linearsystem< T, MATRIX, VECTOR >&);


/**
 * @brief Description of a linear system (suitable for sparse matrix solvers.)
 */
template<
    typename T,
    typename MATRIX=detail::dense_matrix_v< T >,
    typename VECTOR=detail::dense_matrix_v< T > >
class linearsystem : public common::Action
{
  // -- Construction and destruction
 public:

  /// Construct the linear system
  linearsystem(const std::string& name,
               const size_t& _size_i=size_t(),
               const size_t& _size_j=size_t(),
               const size_t& _size_k=1,
               const double& _value=T() ) :
    common::Action(name),
    m_dummy_value(std::numeric_limits< T >::quiet_NaN())
  {
    // framework scripting: options level, signals and options
    mark_basic();

    regist_signal("initialize").connect( boost::bind( &linearsystem::signal_initialize, this, _1 )).signature( boost::bind( &linearsystem::ask_ijk, this, _1 ));
    regist_signal("zerorow")   .connect( boost::bind( &linearsystem::signal_zerorow,    this, _1 )).signature( boost::bind( &linearsystem::ask_i,   this, _1 ));
    regist_signal("clear")     .connect( boost::bind( &linearsystem::signal_clear,      this ));
    regist_signal("solve")     .connect( boost::bind( &linearsystem::signal_solve,      this ));
    regist_signal("output")    .connect( boost::bind( &linearsystem::signal_output,     this ));

    options().add("A",std::vector< T >()).mark_basic().link_to(&m_swap).attach_trigger(boost::bind( &linearsystem::trigger_A, this ));
    options().add("b",std::vector< T >()).mark_basic().link_to(&m_swap).attach_trigger(boost::bind( &linearsystem::trigger_b, this ));
    options().add("x",std::vector< T >()).mark_basic().link_to(&m_swap).attach_trigger(boost::bind( &linearsystem::trigger_x, this ));

    /*
     * note: you cannot call any pure methods (or any method that would call a
     * pure one) in the constructor! specifically, initialize() can't be called
     * here, so the child constructor has to call initialize() itself.
     */
  }

  /// Destruct the linear system
  virtual ~linearsystem() {}


  // -- Framework scripting
 private:

  void signal_initialize (common::SignalArgs& args)
  {
    common::XML::SignalOptions opts(args);
    initialize(
      opts.value< size_t >("i"),
      opts.value< size_t >("j"),
      opts.check("k")? opts.value< size_t >("k") : 1,
      T());
  }

  void signal_zerorow    (common::SignalArgs& args) { common::XML::SignalOptions opts(args); zerorow(opts.value< unsigned int >("i")); }
  void signal_clear  () { clear(); }
  void signal_solve  () { solve(); }
  void signal_output () { operator<<(std::cout,*this); }

  void ask_ijk   (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("i",size_t(1)); opts.add("j",size_t(1)); opts.add("k",size_t(1));  }
  void ask_i     (common::SignalArgs& args) { common::XML::SignalOptions opts(args); opts.add("i",size_t(1)); }

  void trigger_A() { try { A().assign(m_swap); } catch(std::logic_error& e) { std::cout << "problem avoided with A" << std::endl; } }
  void trigger_b() { try { b().assign(m_swap); } catch(std::logic_error& e) { std::cout << "problem avoided with b" << std::endl; } }
  void trigger_x() { try { x().assign(m_swap); } catch(std::logic_error& e) { std::cout << "problem avoided with x" << std::endl; } }

#if 0
  void swap_lss_vector(std::vector< T >& v) {
      if (m_swap.size()!=size())
        throw common::BadValue(FromHere(), "linearsystem is not of the same size as given matrix ("
          + boost::lexical_cast< std::string >(size(0)) + "*"
          + boost::lexical_cast< std::string >(size(1)) + "!="
          + boost::lexical_cast< std::string >(m_swap.size()) + ")");
      for (size_t r=0; r<size(0); ++r)
        for (size_t c=0; c<size(1); ++c)
          A(r,c) = m_swap[r*size(1)+c];
      m_swap.clear();

    if (m_swap.size()!=v.size())
      throw common::BadValue(FromHere(), "linearsystem is not of the same size as given vector ("
        + boost::lexical_cast< std::string >(v.size()) + "!="
        + boost::lexical_cast< std::string >(m_swap.size()) + ")");
    m_swap.swap(v);
    m_swap.clear();
  }
#endif


  // -- Basic functionality
 public:

  /// Linear system solving, aliased from execute
  void execute() { solve(); }

  /// Initialize the linear system
  virtual linearsystem& initialize(
      const size_t& _size_i=size_t(),
      const size_t& _size_j=size_t(),
      const size_t& _size_k=1,
      const double& _value=double())
  {
    const T value(static_cast< T >(_value));
    A().initialize(_size_i,_size_j,value);
    b().initialize(_size_i,_size_k,value);
    x().initialize(_size_j,_size_k,value);

    consistent(A().size(0),A().size(1),b().size(0),b().size(1),x().size(0),x().size(1));
    return *this;
  }

  /// Initialize linear system from file(s)
  virtual linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) {
    if (_Afname.length()) A().initialize(_Afname);
    if (_bfname.length()) b().initialize(_bfname); else b().initialize(size(0),1);
    if (_xfname.length()) x().initialize(_xfname); else x().initialize(size(1),size(2));

    consistent(A().size(0),A().size(1),b().size(0),b().size(1),x().size(0),x().size(1));
    return *this;
  }

  /// Initialize linear system from vectors of values (lists, in the right context)
  virtual linearsystem& initialize(
      const std::vector< double >& vA,
      const std::vector< double >& vb=std::vector< double >(),
      const std::vector< double >& vx=std::vector< double >()) {

    // input parameter type insulation from the templated base class
    const bool conversion_needed(typeid(T)!=typeid(double));
    std::vector< T >
      another_A,
      another_b,
      another_x;
    if (conversion_needed) {
      another_A.resize(vA.size()),
      another_b.resize(vb.size()),
      another_x.resize(vx.size());
      std::transform( vA.begin(),vA.end(),another_A.begin(),detail::storage_conversion_t< double, T >() );
      std::transform( vb.begin(),vb.end(),another_b.begin(),detail::storage_conversion_t< double, T >() );
      std::transform( vx.begin(),vx.end(),another_x.begin(),detail::storage_conversion_t< double, T >() );
    }

    // propper initialization of linear system components
    std::vector< T >&
      correct_A(conversion_needed? another_A : (std::vector< T >&) vA),
      correct_b(conversion_needed? another_b : (std::vector< T >&) vb),
      correct_x(conversion_needed? another_x : (std::vector< T >&) vx);
    if (vA.size()) A().initialize(correct_A);
    if (vb.size()) b().initialize(correct_b); else b().initialize(size(0),1);
    if (vx.size()) x().initialize(correct_x); else x().initialize(size(1),size(2));

    consistent(A().size(0),A().size(1),b().size(0),b().size(1),x().size(0),x().size(1));
    return *this;
  }

  /// Clear contents in all system components
  virtual linearsystem& clear() {
    A().clear();
    b().clear();
    x().clear();
    return *this;
  }

  /// Zero row in all system components
  virtual linearsystem& zerorow(const size_t r) {
    A().zerorow(r);
    b().zerorow(r);
    x().zerorow(r);
    return *this;
  }

  /// Returns the specific dimension of the system
  size_t size(const size_t& d) const {
    return (d< 2? A().size(d) :
           (d==2? b().size(1) : 0));
  }

  /// Checks whether the linear system matrix is empty
  bool empty() { return !(size(0)*size(1)*size(2)); }



  // -- Internal functionality
 private:

  /// Checks whether the matrix/vectors sizes are consistent in the system
  bool consistent(const size_t& Ai, const size_t& Aj,
                  const size_t& bi, const size_t& bj,
                  const size_t& xi, const size_t& xj) const {
    if (Ai!=bi || Aj!=xi || bj!=xj) {
      std::ostringstream msg;
      msg << "linearsystem: size is not consistent: "
          << "A(" << Ai << 'x' << Aj << ") "
          << "x(" << xi << 'x' << xj << ") = "
          << "b(" << bi << 'x' << bj << ") ";
      throw std::runtime_error(msg.str());
      return false;
    }
    return true;
  }

  /// Output
  template< typename aT, typename aMATRIX, typename aVECTOR >
  friend std::ostream& operator<< (std::ostream&, const linearsystem< aT, aMATRIX, aVECTOR >&);
  // (a private generic friend? how promiscuous!)


  // -- Storage
 protected:

  T m_dummy_value;          // should keep NaN throughout
  std::vector< T > m_swap;  // scripting temporary storage (to swap with internal storage)


  // -- Pure methods
 public:

  /// Linear system solving and components access
  virtual linearsystem& solve() = 0;
  virtual       MATRIX& A()       = 0;
  virtual       VECTOR& b()       = 0;
  virtual       VECTOR& x()       = 0;
  virtual const MATRIX& A() const = 0;
  virtual const VECTOR& b() const = 0;
  virtual const VECTOR& x() const = 0;

#if 0
  /// Linear system copy
  linearsystem& operator=(const linearsystem& _other) {
    A() = _other.A();
    b() = _other.b();
    x() = _other.x();
    return *this;
  }

  /// Value assignment (operator)
  linearsystem& operator=(const T& _value) { return initialize(size(0),size(1),size(2),_value); }

  /// Value assignment (method)
  linearsystem& assign(const T& _value=T()) { return operator=(_value); }
#endif


};


#if 0
class xpto
{
  virtual void output_A(std::ostream& out) const { out << "[ (unavailable) ]"; }
  virtual void output_b(std::ostream& out) const { out << "[ "; std::copy(m_b.begin(),m_b.end(),std::ostream_iterator< T >(out,", ")); out << ']'; }
  virtual void output_x(std::ostream& out) const { out << "[ "; std::copy(m_x.begin(),m_x.end(),std::ostream_iterator< T >(out,", ")); out << ']'; }
  friend std::ostream& operator<< < T, MATRIX, VECTOR >(std::ostream&, const linearsystem&);
};
#endif


/// Output to given std::ostream (non-member version)
template< typename T, typename MATRIX, typename VECTOR >
std::ostream& operator<< (std::ostream& o, const linearsystem< T, MATRIX, VECTOR >& lss)
{
#if 0
  o << "linearsystem: A: "; lss.output_A(o); o << std::endl;
  o << "linearsystem: b: "; lss.output_b(o); o << std::endl;
  o << "linearsystem: x: "; lss.output_x(o); o << std::endl;
#endif
  return o;
}


}  // namespace lss
}  // namespace cf3


#endif
