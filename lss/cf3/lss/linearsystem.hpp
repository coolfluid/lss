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

#include "common/BasicExceptions.hpp"
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
  linearsystem(const std::string& name) :
    common::Action(name),
    m_dummy_value(std::numeric_limits< T >::quiet_NaN())
  {
    // framework scripting: options level, signals and options
    mark_basic();

    regist_signal("initialize")
        .description("Initialize with given size (i,j,k) and value (value), or from file (A, b, and/or x)")
        .connect   ( boost::bind( &linearsystem::signal_initialize, this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_ijkvalue,   this, _1 ));

    regist_signal("zerorow")
        .description("Erase given row (i) in all components")
        .connect   ( boost::bind( &linearsystem::signal_zerorow,  this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_ijkvalue, this, _1 ));

    regist_signal("clear") .connect( boost::bind( &linearsystem::signal_clear,  this )).description("Empty linear system components");
    regist_signal("solve") .connect( boost::bind( &linearsystem::signal_solve,  this )).description("Solve linear system, returning solution in x()");
    regist_signal("output").connect( boost::bind( &linearsystem::signal_output, this )).description("Print a pretty linear system");

    regist_signal("A").connect(boost::bind( &linearsystem::signal_A, this, _1 )).signature(boost::bind( &linearsystem::signat_ijkvalue, this, _1 )).description("Set entry in matrix A, by given index (i,j) and value (value)");
    regist_signal("b").connect(boost::bind( &linearsystem::signal_b, this, _1 )).signature(boost::bind( &linearsystem::signat_ijkvalue, this, _1 )).description("Set entry in vector b, by given index (i,k) and value (value)");
    regist_signal("x").connect(boost::bind( &linearsystem::signal_x, this, _1 )).signature(boost::bind( &linearsystem::signat_ijkvalue, this, _1 )).description("Set entry in vector x, by given index (j,k) and value (value)");

    options().add("A",std::vector< double >())
        .link_to(&m_dummy_vector).mark_basic()
        .attach_trigger(boost::bind( &linearsystem::trigger_A, this ));

    options().add("b",std::vector< double >())
        .link_to(&m_dummy_vector).mark_basic()
        .attach_trigger(boost::bind( &linearsystem::trigger_b, this ));

    options().add("x",std::vector< double >())
        .link_to(&m_dummy_vector).mark_basic()
        .attach_trigger(boost::bind( &linearsystem::trigger_x, this ));

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


  void signat_ijkvalue(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    opts.add< unsigned >("i");
    opts.add< unsigned >("j");
    opts.add< unsigned >("k",1);
    opts.add< std::string >("A");
    opts.add< std::string >("b");
    opts.add< std::string >("x");
    opts.add< double >("value");
  }

  void signal_initialize(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    const std::string
        Afname(opts.value< std::string >("A")),
        bfname(opts.value< std::string >("b")),
        xfname(opts.value< std::string >("x"));
    const double value(opts.value< double >("value"));
      if (Afname.length() || xfname.length() || bfname.length()) {
        if (Afname.length()) component_initialize_with_file(A(),"A",Afname);
        if (bfname.length()) component_initialize_with_file(b(),"b",bfname); else b().initialize(size(0),1);
        if (xfname.length()) component_initialize_with_file(x(),"x",xfname); else x().initialize(size(1),size(2));
//        consistent(A().size(0),A().size(1),b().size(0),b().size(1),x().size(0),x().size(1));
      }
      else {
        const unsigned
            i(opts.value< unsigned >("i")),
            j(opts.value< unsigned >("j")),
            k(opts.value< unsigned >("k"));
        A().initialize(i,j,value);
        b().initialize(i,k,value);
        x().initialize(j,k,value);
      }
  }


  void signal_zerorow(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    zerorow(opts.value< unsigned >("i"));
  }

  void signal_clear () { clear(); }
  void signal_solve () { solve(); }
  void signal_output() { operator<<(std::cout,*this); }

  void signal_A(common::SignalArgs& args) { common::XML::SignalOptions opts(args); A().operator()(opts.value< unsigned >("i"),opts.value< unsigned >("j")) = opts.value< double   >("value"); }
  void signal_b(common::SignalArgs& args) { common::XML::SignalOptions opts(args); b().operator()(opts.value< unsigned >("i"),opts.value< unsigned >("k")) = opts.value< unsigned >("value"); }
  void signal_x(common::SignalArgs& args) { common::XML::SignalOptions opts(args); x().operator()(opts.value< unsigned >("j"),opts.value< unsigned >("k")) = opts.value< unsigned >("value"); }

  void trigger_A() { component_initialize_with_vector< MATRIX >(A(),"A"); }
  void trigger_b() { component_initialize_with_vector< VECTOR >(b(),"b"); }
  void trigger_x() { component_initialize_with_vector< VECTOR >(x(),"x"); }


  template< typename COMP >
  void component_initialize_with_vector(COMP& c, const std::string& name) {
    try { c.initialize(m_dummy_vector); }
    catch (const std::runtime_error& e) {
      std::cout << "linearsystem: " << name << ": " << e.what() << std::endl;
    }
    m_dummy_vector.clear();
  }

  template< typename COMP >
  void component_initialize_with_file(COMP& c, const std::string& name, const std::string& fname) {
    try { c.initialize(fname); }
    catch (const std::runtime_error& e) {
      std::cout << "linearsystem: " << name << ": " << e.what() << std::endl;
    }
    m_dummy_vector.clear();
  }


  // -- Basic functionality
 public:

  /// Linear system solving, aliased from execute
  void execute() {
    try { solve(); }
    catch (const std::runtime_error& e) {
      std::cout << "linearsystem: " << e.what() << std::endl;
    }
  }

  /// Initialize the linear system
  virtual linearsystem& initialize(
      const size_t& i=size_t(),
      const size_t& j=size_t(),
      const size_t& k=1,
      const double& _value=double())
  {
    A().initialize(i,j,_value);
    b().initialize(i,k,_value);
    x().initialize(j,k,_value);
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
    if (vA.size()) A().initialize(vA);
    if (vb.size()) b().initialize(vb); else b().initialize(size(0),1);
    if (vx.size()) x().initialize(vx); else x().initialize(size(1),size(2));
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
  virtual linearsystem& zerorow(const size_t& i) {
    A().zerorow(i);
    b().zerorow(i);
    x().zerorow(i);
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
          << "b(" << bi << 'x' << bj << ").";
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

  T m_dummy_value;               // should keep NaN throughout
  std::vector< double > m_dummy_vector;  // scripting temporary storage (to swap with internal storage)


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


/// Output to given std::ostream (non-member version)
template< typename T, typename MATRIX, typename VECTOR >
std::ostream& operator<< (std::ostream& o, const linearsystem< T, MATRIX, VECTOR >& lss)
{
  o << "linearsystem: A: "; lss.A().print(o); o << std::endl;
  o << "linearsystem: b: "; lss.b().print(o); o << std::endl;
  o << "linearsystem: x: "; lss.x().print(o); o << std::endl;
  return o;
}


}  // namespace lss
}  // namespace cf3


#endif
