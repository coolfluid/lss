
#ifndef cf3_lss_dummy_linearsystem_h
#define cf3_lss_dummy_linearsystem_h


#include "linearsystem.hpp"
#include "lss_matrix.hpp"


namespace cf3 {
namespace lss {


struct dummy_linearsystem :
  public linearsystem<
    double,
    lss_matrix::dense_matrix_v< double > >
{
  typedef linearsystem< double, lss_matrix::dense_matrix_v< double > > linearsystem_t;

  // cons/destructor
  dummy_linearsystem(const std::string& name) : linearsystem_t(name) {}
  ~dummy_linearsystem() {}

  // solve
  dummy_linearsystem& solve();

};


}  // namespace lss
}  // namespace cf3


#endif
