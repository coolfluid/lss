// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_linearsystem_hpp
#define cf3_lss_linearsystem_hpp


#include <iostream>

#include "common/BasicExceptions.hpp"
#include "common/Signal.hpp"
#include "common/Action.hpp"

#include "utilities.hpp"
#include "matrix.hpp"


namespace cf3 {
namespace lss {


/* -- linear system --------------------------------------------------------- */


// helper forward declarations
template< typename T > class linearsystem;
template< typename T > std::ostream& operator<< (std::ostream&, const linearsystem< T >&);


/**
 * @brief Description of a linear system (suitable for sparse matrix solvers.)
 */
template< typename T >
class linearsystem : public common::Action
{
  // -- Construction and destruction
 public:

  /// Definition of linear system vector (as a column-oriented matrix)
  typedef dense_matrix_v< T, sort_by_column > vector_t;

  /// Construct the linear system
  linearsystem(const std::string& name) :
    common::Action(name),
    m_dummy_value(std::numeric_limits< T >::signaling_NaN())
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

    regist_signal("output")
        .description("Print a pretty linear system, at print level per component where 0:auto (default), 1:size, 2:signs, and 3:full")
        .connect   ( boost::bind( &linearsystem::signal_output,  this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_abcfile, this, _1 ));

    regist_signal("clear")
        .description("Empty linear system components")
        .connect( boost::bind( &linearsystem::signal_clear, this ));

    regist_signal("solve")
        .description("Linear system solving, x = A^-1 b")
        .connect( boost::bind( &linearsystem::signal_solve, this ));

    regist_signal("multi")
        .description("Linear system forward multiplication, b = alpha A x + beta b")
        .connect   ( boost::bind( &linearsystem::signal_multi, this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_multi, this, _1 ));

    regist_signal("A")
        .description("Set entry in matrix A, by given index (i,j) and value (value)")
        .connect   ( boost::bind( &linearsystem::signal_A,        this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_ijkvalue, this, _1 ));

    regist_signal("b")
        .description("Set entry in vector b, by given index (i,k) and value (value)")
        .connect   ( boost::bind( &linearsystem::signal_b,        this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_ijkvalue, this, _1 ));

    regist_signal("x")
        .description("Set entry in vector x, by given index (j,k) and value (value)")
        .connect   ( boost::bind( &linearsystem::signal_x,        this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_ijkvalue, this, _1 ));

    regist_signal("bnorm")
        .description("Calculate vector b p-norm, given column index (j) and Lp space dimension (p)")
        .connect   ( boost::bind( &linearsystem::signal_bnorm, this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_jp,    this, _1 ));

    regist_signal("xnorm")
        .description("Calculate vector x p-norm, given column index (j) and Lp space dimension (p)")
        .connect   ( boost::bind( &linearsystem::signal_xnorm, this, _1 ))
        .signature ( boost::bind( &linearsystem::signat_jp,    this, _1 ));

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

  void signat_abcfile(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    opts.add< int >("A",(int) print_auto);
    opts.add< int >("b",(int) print_auto);
    opts.add< int >("x",(int) print_auto);
    opts.add< std::string >("file","");
  }

  void signat_multi(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    opts.add< double >("alpha",1.);
    opts.add< double >("beta", 0.);
  }

  void signat_jp(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    opts.add< unsigned >("j",(unsigned) 0);
    opts.add< double >("p",2.);
  }

  void signal_initialize(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    const double value(opts.value< double >("value"));
    const std::string
        Afname(opts.value< std::string >("A")),
        bfname(opts.value< std::string >("b")),
        xfname(opts.value< std::string >("x"));
    if (Afname.length() || bfname.length() || xfname.length()) {
      if (Afname.length()) A___initialize(Afname);
      if (bfname.length() && !component_initialize_with_file(m_b,"b",bfname) || !bfname.length()) m_b.initialize(size(0),1      );
      if (xfname.length() && !component_initialize_with_file(m_x,"x",xfname) || !xfname.length()) m_x.initialize(size(1),size(2));
      consistent(A___size(0),A___size(1),m_b.size(0),m_b.size(1),m_x.size(0),m_x.size(1));
    }
    else {
      const unsigned
          i(opts.value< unsigned >("i")),
          j(opts.value< unsigned >("j")),
          k(opts.value< unsigned >("k"));
      A___initialize(i,j);
      m_b.initialize(i,k);
      m_x.initialize(j,k);
    }
  }

  void signal_zerorow(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    zerorow(opts.value< unsigned >("i"));
  }

  void signal_output(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    using namespace std;

    m_print[0] = static_cast< print_t >(opts.value< int >("A"));
    m_print[1] = static_cast< print_t >(opts.value< int >("b"));
    m_print[2] = static_cast< print_t >(opts.value< int >("x"));
    string bname = opts.value< string >("file");
    if (bname.length()) {
      try {
        struct fhelper {
          ofstream f;
          fhelper(const bool& _print, const string& _fname)
            : f(_print? _fname.c_str():"") {
            if (_print && !f)
              throw runtime_error("cannot write to file \""+_fname+"\"");
          }
        } hA(m_print[0],bname+"_A.mtx"),
          hb(m_print[1],bname+"_b.mtx"),
          hx(m_print[2],bname+"_x.mtx");
        if (m_print[0]) A___print(hA.f,print_file);
        if (m_print[1]) m_b.print(hb.f,print_file);
        if (m_print[2]) m_x.print(hx.f,print_file);
      }
      catch (const std::runtime_error& e) {
        CFwarn << "linearsystem: " << e.what() << CFendl;
      }
    }
    else {
      operator<<(cout,*this);
    }
  }

  void signal_clear() { clear(); }

  void signal_solve() { execute(); }

  void signal_multi(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    multi(opts.value< double >("alpha"),opts.value< double >("beta"));
  }

  void signal_A(common::SignalArgs& args) {
    common::XML::SignalFrame reply(args.create_reply(uri()));
    common::XML::SignalOptions opts(args), repl(reply);
    repl.add("return_value",A(opts.value< unsigned >("i"),opts.value< unsigned >("j")) = opts.value< double >("value"));
  }

  void signal_b(common::SignalArgs& args) {
    common::XML::SignalFrame reply(args.create_reply(uri()));
    common::XML::SignalOptions opts(args), repl(reply);
    repl.add("return_value",m_b(opts.value< unsigned >("i"),opts.value< unsigned >("k")) = opts.value< double >("value"));
  }

  void signal_x(common::SignalArgs& args) {
    common::XML::SignalFrame reply(args.create_reply(uri()));
    common::XML::SignalOptions opts(args), repl(reply);
    repl.add("return_value",m_x(opts.value< unsigned >("i"),opts.value< unsigned >("k")) = opts.value< double >("value"));
  }

  void signal_bnorm(common::SignalArgs& args) {
    common::XML::SignalFrame reply(args.create_reply(uri()));
    common::XML::SignalOptions
      opts(args),
      repl(reply);
    repl.add("return_value",m_b.norm(opts.value< unsigned >("j"),opts.value< double >("p")));
  }

  void signal_xnorm(common::SignalArgs& args) {
    common::XML::SignalFrame reply(args.create_reply(uri()));
    common::XML::SignalOptions
      opts(args),
      repl(reply);
    repl.add("return_value",m_x.norm(opts.value< unsigned >("j"),opts.value< double >("p")));
  }

  void trigger_A() { try { A___initialize(m_dummy_vector); } catch (const std::runtime_error& e) { CFwarn << "linearsystem: A: " << e.what() << CFendl; } m_dummy_vector.clear(); }
  void trigger_b() { try { m_b.initialize(m_dummy_vector); } catch (const std::runtime_error& e) { CFwarn << "linearsystem: b: " << e.what() << CFendl; } m_dummy_vector.clear(); }
  void trigger_x() { try { m_x.initialize(m_dummy_vector); } catch (const std::runtime_error& e) { CFwarn << "linearsystem: x: " << e.what() << CFendl; } m_dummy_vector.clear(); }

  bool component_initialize_with_file(dense_matrix_v< T >& c, const std::string& name, const std::string& fname) {
    try { c.initialize(fname); }
    catch (const std::runtime_error& e) {
      CFwarn << "linearsystem: " << name << ": " << e.what() << CFendl;
      return false;
    }
    m_dummy_vector.clear();
    return true;
  }


  // -- Basic functionality
 public:

  /// Linear system solving, aliased from execute
  void execute() {
    try { solve(); }
    catch (const std::runtime_error& e) {
      CFwarn << "linearsystem: " << e.what() << CFendl;
    }
  }

  /// Linear system forward multiplication
  linearsystem& multi(const vector_t& _x, vector_t& _b) {
    try {
      consistent(A___size(0),A___size(1),_b.size(0),_b.size(1),_x.size(0),_x.size(1));
      A___multi(_x,_b);
    }
    catch (const std::runtime_error& e) {
      CFwarn << "linearsystem: " << e.what() << CFendl;
    }
    return *this;
  }

  /// Linear system initialization
  linearsystem& initialize(
      const size_t& i=size_t(),
      const size_t& j=size_t(),
      const size_t& k=1,
      const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >() )
  {
    A___initialize(i,j,_nnz);
    m_b.initialize(i,k);
    m_x.initialize(j,k);
    return *this;
  }

  /// Linear system initialization, from file(s)
  linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) {
    if (_Afname.length()) A___initialize(_Afname);
    if (_bfname.length()) m_b.initialize(_bfname); else m_b.initialize(size(0),1);
    if (_xfname.length()) m_x.initialize(_xfname); else m_x.initialize(size(1),size(2));
    consistent(A___size(0),A___size(1),m_b.size(0),m_b.size(1),m_x.size(0),m_x.size(1));
    return *this;
  }

  /// Linear system initialization, from vectors of values (lists, in the right context)
  linearsystem& initialize(
      const std::vector< double >& vA,
      const std::vector< double >& vb=std::vector< double >(),
      const std::vector< double >& vx=std::vector< double >()) {
    if (vA.size()) A___initialize(vA);
    if (vb.size()) m_b.initialize(vb); else m_b.initialize(size(0),1);
    if (vx.size()) m_x.initialize(vx); else m_x.initialize(size(1),size(2));
    consistent(A___size(0),A___size(1),m_b.size(0),m_b.size(1),m_x.size(0),m_x.size(1));
    return *this;
  }

  /// Clear contents in all system components
  linearsystem& clear() {
    A___clear();
    m_b.clear();
    m_x.clear();
    return *this;
  }

  /// Zero row in all system components
  linearsystem& zerorow(const size_t& i) {
    A___zerorow(i);
    m_b.zerorow(i);
    m_x.zerorow(i);
    return *this;
  }

  /// Sum entries into row, from another given row
  linearsystem& sumrows(const size_t& i, const size_t& isrc) {
    A___sumrows(i,isrc);
    m_b.sumrows(i,isrc);
    m_x.sumrows(i,isrc);
    return *this;
  }

  /// Value assignment (method)
  linearsystem& assign(const double& _value=double()) {
    A___assign(_value);
    m_b = _value;
    m_x = _value;
    return *this;
  }

  /// Operators assign value and copy
  linearsystem& operator=(const T& _value)            { return assign(_value); }
  linearsystem& operator=(const linearsystem& _other) { return copy(_other);   }

  /// Output system components
  linearsystem& output(
      std::ostream& _sA, const print_t& _lA,
      std::ostream& _sb, const print_t& _lb,
      std::ostream& _sx, const print_t& _lx ) {
    A___print(_sA,(m_print[0]=_lA));
    m_b.print(_sb,(m_print[1]=_lb));
    m_x.print(_sx,(m_print[2]=_lx));
    return *this;
  }

  /// Returns the specific dimension of the system
  size_t size(const size_t& d) const {
    return (d< 2? A___size(d) : (d==2? m_b.size(1) : 0));
  }

  /// Checks whether the linear system matrix is empty
  bool empty() { return !(size(0)*size(1)*size(2)); }

  /// Linear system b vector access
  vector_t& b() { return m_b; }

  /// Linear system x vector access
  vector_t& x() { return m_x; }

  /// Linar system componets indexing (absolute)
  virtual const T& A(const size_t& i, const size_t& j)   const = 0;
  virtual       T& A(const size_t& i, const size_t& j)         = 0;
          const T& b(const size_t& i, const size_t& j=0) const { return m_b(i,j); }
          const T& x(const size_t& i, const size_t& j=0) const { return m_x(i,j); }
                T& b(const size_t& i, const size_t& j=0)       { return m_b(i,j); }
                T& x(const size_t& i, const size_t& j=0)       { return m_x(i,j); }


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
  template< typename aT >
  friend std::ostream& operator<< (std::ostream&, const linearsystem< aT >&);
  // (a private generic friend? how promiscuous!)


  // -- Storage
 protected:

  vector_t m_b;  // linear system b vector
  vector_t m_x;  // linear system x vector

  T                     m_dummy_value;   // (scripting temp.) value
  std::vector< double > m_dummy_vector;  // (scripting temp.) vector
  print_t m_print[3];                    // (scripting temp.) print levels


  // -- Interfacing (public)
 public:

  /// Linear system solving: x = A^-1 b
  /// @note: might destroy system matrix contents (structure or non-zero values)
  virtual linearsystem& solve() = 0;

  /// Linear system forward multiplication: b = alpha A x + beta b
  /// @note: might should not destroy system matrix contents (structure or non-zero values)
  virtual linearsystem& multi(const double& _alpha, const double& _beta) = 0;

  /// Linear system copy
  /// @note (to be specialized in concrete derived class to account for matrix)
  virtual linearsystem& copy(const linearsystem& _other) {
    m_b = _other.m_b;
    m_x = _other.m_x;
    std::copy(_other.m_print,_other.m_print+3,m_print);
    m_dummy_value  = _other.m_dummy_value;
    m_dummy_vector = _other.m_dummy_vector;
    return *this;
  }

  /// Linear system swap
  /// @note (to be specialized in concrete derived class to account for matrix)
  virtual linearsystem& swap(linearsystem& _other) {
    m_b.swap(_other.m_b);
    m_x.swap(_other.m_x);
    print_t *pa1(m_print), *pa2(_other.m_print);
    std::swap(pa1,pa2);
    std::swap(m_dummy_value,_other.m_dummy_value);
    m_dummy_vector.swap(_other.m_dummy_vector);
    return *this;
  }

  // -- Interfacing (protected)
 protected:

  /// Linear system matrix modifiers
  virtual void A___initialize(const size_t& i, const size_t& j, const std::vector< std::vector< size_t > >& _nnz=std::vector< std::vector< size_t > >()) = 0;
  virtual void A___initialize(const std::vector< double >& _vector) = 0;
  virtual void A___initialize(const std::string& _fname)            = 0;
  virtual void A___assign(const double& _value)                 = 0;
  virtual void A___clear()                                      = 0;
  virtual void A___zerorow(const size_t& i)                     = 0;
  virtual void A___sumrows(const size_t& i, const size_t& isrc) = 0;

  /// Linear system matrix inspecting
  virtual void   A___print(std::ostream& o, const print_t &l=print_auto) const = 0;
  virtual size_t A___size(const size_t& d ) const = 0;

};


/// Output to given std::ostream (non-member version)
template< typename T >
std::ostream& operator<< (std::ostream& o, const linearsystem< T >& lss)
{
  o << "linearsystem: A: "; lss.A___print(o,lss.m_print[0]); o << std::endl;
  o << "linearsystem: b: "; lss.m_b.print(o,lss.m_print[1]); o << std::endl;
  o << "linearsystem: x: "; lss.m_x.print(o,lss.m_print[2]); o << std::endl;
  return o;
}


}  // namespace lss
}  // namespace cf3


#endif
