// Copyright (C) 2014 Vrije Universiteit Brussel, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#ifndef cf3_lss_nonlinearsystem_hpp
#define cf3_lss_nonlinearsystem_hpp


#include "LibLSS.hpp"
#include "linearsystem.hpp"


namespace cf3 {
namespace lss {


/* -- non-linear system ----------------------------------------------------- */


/**
 * @brief Description of a non-linear system, suitable for trust region and line
 * search strategies
 */
template< typename T >
class nonlinearsystem : public common::Action
{
  // -- Construction and destruction
 public:

  /// Construct the non-linear system
  nonlinearsystem(const std::string& name) :
    common::Action(name)
  {
    // framework scripting: options level, signals and options
    mark_basic();

    options().add("linearsystem", h_linearsystem)
      .description("Linear system to use")
      .link_to(&h_linearsystem).mark_basic()
      .attach_trigger(boost::bind( &nonlinearsystem::linearsystem_recreate, this ));

    regist_signal("clear")
      .description("Non-linear system components clearing")
      .connect( boost::bind( &nonlinearsystem::signal_clear, this ));

    regist_signal("solve")
      .description("Non-linear system solving")
      .connect( boost::bind( &nonlinearsystem::signal_solve, this ));

    regist_signal("setup")
      .description("Setup which linear system to use")
      .connect(   boost::bind( &nonlinearsystem::signal_setup, this, _1 ))
      .signature( boost::bind( &nonlinearsystem::signat_setup, this, _1 ));

    regist_signal("linearsystem")
      .description("Linear system: access currently active linearsystem")
      .connect(   boost::bind( &nonlinearsystem::signal_linearsystem, this, _1 ))
      .signature( boost::bind( &nonlinearsystem::signat_none,         this, _1 ));

    regist_signal("swap")
      .description("Linear system: swap active linear system contents, creating internal copy if necessary (if this has not occurred before)")
      .connect(   boost::bind( &nonlinearsystem::signal_swap, this, _1 ))
      .signature( boost::bind( &nonlinearsystem::signat_none, this, _1 ));
  }

  /// Destruct the non-linear system solver
  virtual ~nonlinearsystem() {}


  // -- Framework scripting
 private:

  void linearsystem_recreate() {
    try {
      linearsystem< T >& ls = linearsystem_get();
      h_linearsystem_pert = Handle< linearsystem< T > >(create_component(
          ls.name()+"_nonlinearsystem_pert",
          ls.derived_type_name() ));
      if (is_null(h_linearsystem_pert))
        throw std::runtime_error("nonlinearsystem: cannot recreate component");
    }
    catch (const std::runtime_error& e) {
      CFwarn << "nonlinearsystem: linearsystem cannot be recreated (" << e.what() << ')' << CFendl;
    }
  }

  void signat_setup(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    opts.add("linearsystem",Handle< linearsystem< T > >());
  }

  void signat_none(common::SignalArgs& args) {}

  void signal_clear() { clear(); }

  void signal_solve() { execute(); }

  void signal_setup(common::SignalArgs& args) {
    common::XML::SignalOptions opts(args);
    if (opts.check("linearsystem")) {
      h_linearsystem = opts.value< Handle< linearsystem< T > > >("linearsystem");
      linearsystem_get();
    }
  }

  void signal_linearsystem(common::SignalArgs& args) {
    common::XML::SignalFrame reply(args.create_reply(uri()));
    common::XML::SignalOptions repl(reply);
    try {
      linearsystem< T >& ls = linearsystem_get();
      repl.add("created_component",ls.uri());
    }
    catch (const std::runtime_error& e) {
      CFwarn << "nonlinearsystem: " << e.what() << CFendl;
    }
  }

  void signal_swap(common::SignalArgs& args) {
    common::XML::SignalFrame reply(args.create_reply(uri()));
    common::XML::SignalOptions repl(reply);
    try {
      linearsystem_swap();
      repl.add("created_component",uri());
    }
    catch (const std::runtime_error& e) {
      CFwarn << "nonlinearsystem: " << e.what() << CFendl;
    }
  }


  // -- Basic functionality
 public:

  /// Non-linear system components clearing
  nonlinearsystem& clear() {
    if (is_not_null(h_linearsystem))       h_linearsystem->clear();
    if (is_not_null(h_linearsystem_pert))  h_linearsystem_pert->clear();
    return *this;
  }

  /// Non-linear system solving, aliased from execute
  void execute() {
    try { solve(); }
    catch (const std::runtime_error& e) {
      CFwarn << "nonlinearsystem: " << e.what() << CFendl;
    }
  }

  /// Linear system: access base, "unperturbed" linearsystem
  linearsystem< T >& linearsystem_get() {
    if (is_null(h_linearsystem))
      throw std::runtime_error("nonlinearsystem: linearsystem not set");
    return *h_linearsystem;
  }

  /// Linear system: access internal, "perturbed" linearsystem
  linearsystem< T >& linearsystem_pert_get() {
    if (is_null(h_linearsystem_pert))
      throw std::runtime_error("nonlinearsystem: linearsystem not recreated");
    return *h_linearsystem_pert;
  }

  /// Linear system: swap contents of "unperturbed" with "perturbed" linear sys.
  linearsystem< T >& linearsystem_swap() {
    return linearsystem_get().swap(linearsystem_pert_get());
  }

  /// Linear system: copy contents of "unperturbed" to "perturbed" linear sys.
  linearsystem< T >& linearsystem_copy() {
    return linearsystem_pert_get().copy(linearsystem_get());
  }


  // -- Storage
 protected:

  /// Handle to linear system, in unperturbed state (provided externally)
  Handle< linearsystem< T > > h_linearsystem;

  /// Handle to linear system, in perturbed state (provided internally)
  Handle< linearsystem< T > > h_linearsystem_pert;


  // -- Interfacing (public)
 public:

  /// Non-linear system solving
  /// @note: relies on the matrix non-destructive solving (structure or non-zero
  /// values) since copying and swapping components is involved
  virtual nonlinearsystem& solve() = 0;

};


}  // namespace lss
}  // namespace cf3


#endif

