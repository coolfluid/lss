#ifndef lss_linearsystem_hpp
#define lss_linearsystem_hpp


#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>


#include "lss_utilities.hpp"
#include "lss_index.hpp"
#include "lss_matrix.hpp"


/* -- linear system --------------------------------------------------------- */

template<
    typename T,
    typename MATRIX=lss_matrix::dense_matrix_vv< T >,
    typename VECTOR=lss_matrix::dense_matrix_vv< T, lss_matrix::column_oriented >,
    typename INDEX_HIERARCHY=lss_index::index_hierarchy_t< lss_index::index_hierarchy_t_end > >
class linearsystem
{ public:  //FIXME update permissions

  // utility definitions

  typedef MATRIX matrix_t;
  typedef VECTOR vector_t;
  typedef INDEX_HIERARCHY index_hierarchy_t;

  // construction, destruction and initialization

  /// Construct the linear system, by direct initialization
  linearsystem(
      const size_t& _size_i=size_t(),
      const size_t& _size_j=size_t(),
      const size_t& _size_k=1,
      const T& _value=T() ) { resize(_size_i,_size_j,_size_k,_value); }

  /// Construct the linear system, by copy
  linearsystem(const linearsystem& _other) { operator=(_other); }

  /// Destructs the linear system
  virtual ~linearsystem() {}

  /// Initialize linear system from file(s)
  linearsystem& initialize(
      const std::string& _Afname,
      const std::string& _bfname="",
      const std::string& _xfname="" ) {
    if (_Afname.length()) m_A.initialize(_Afname);
    if (_bfname.length()) m_b.initialize(_bfname); else m_b.resize(size(0),1);
    if (_xfname.length()) m_x.initialize(_xfname); else m_x.resize(size(1),size(2));
    if (size(0)!=m_b.size(0) || size(1)!=m_x.size(0)
     || size(2)!=m_b.size(1) || size(2)!=m_x.size(1)) {
      std::ostringstream msg;
      msg << "linearsystem: inconsistent size: "
          << "A(" <<     size(0) << 'x' <<     size(1) << ") "
          << "x(" << m_x.size(0) << 'x' << m_x.size(1) << ") = "
          << "b(" << m_b.size(0) << 'x' << m_b.size(1) << ") ";
      throw std::runtime_error(msg.str());
    }
    return *this;
  }

  // accessing

        T& A(const size_t& i, const size_t& j)         { return m_A(i,j); }
        T& b(const size_t& i, const size_t& j=0)       { return m_b(i,j); }
        T& x(const size_t& i, const size_t& j=0)       { return m_x(i,j); }
  const T& A(const size_t& i, const size_t& j)   const { return m_A(i,j); }
  const T& b(const size_t& i, const size_t& j=0) const { return m_b(i,j); }
  const T& x(const size_t& i, const size_t& j=0) const { return m_x(i,j); }

  // interfacing

  /// Assign values to the linear system
  linearsystem& operator=(const T& _value) { return resize(size(0),size(1),size(2),_value); }

  /// Assign values to the linear system
  linearsystem& operator=(const linearsystem& _other) {
    m_A = _other.m_A;
    m_b = _other.m_b;
    m_x = _other.m_x;
    return *this;
  }

  /// assigns values to the linear system, specifying an initialization value
  linearsystem& assign(const T& _value=T()) { return resize(m_b.size(),m_x.size(),_value); }

  /// assign fom vectors of values (lists, in the right context)
  linearsystem& assign(
      const std::vector< T >& vA,
      const std::vector< T >& vb=std::vector< T >(),
      const std::vector< T >& vx=std::vector< T >()) {
    m_A.assign(vA);
    if (vb.size()) m_b.assign(vb); else m_b.resize(size(0),1);
    if (vx.size()) m_x.assign(vx); else m_x.resize(size(1),size(2));
    if (size(0)!=m_b.size(0) || size(1)!=m_x.size(0)
     || size(2)!=m_b.size(1) || size(2)!=m_x.size(1)) {
      std::ostringstream msg;
      msg << "linearsystem: inconsistent size: "
          << "A(" <<     size(0) << 'x' <<     size(1) << ") "
          << "x(" << m_x.size(0) << 'x' << m_x.size(1) << ") = "
          << "b(" << m_b.size(0) << 'x' << m_b.size(1) << ") ";
      throw std::runtime_error(msg.str());
    }
    return *this;
  }

  /// Changes the number of elements stored
  linearsystem& resize(const size_t& _size_i, const size_t& _size_j, const size_t& _size_k=1, const T& _value=T()) {
    if (idx_t(_size_i,_size_j)!=m_A.m_size) m_A.clear();
    if (idx_t(_size_i,_size_k)!=m_b.m_size) m_b.clear();
    if (idx_t(_size_j,_size_k)!=m_x.m_size) m_x.clear();
    if (idx_t(_size_i,_size_j).is_valid_size()) m_A.resize(_size_i,_size_j,_value);
    if (idx_t(_size_i,_size_k).is_valid_size()) m_b.resize(_size_i,_size_k,_value);
    if (idx_t(_size_j,_size_k).is_valid_size()) m_x.resize(_size_j,_size_k,_value);
    return *this;
  }

  /// Clears the contents
  linearsystem& clear() {
    m_A.clear();
    m_b.clear();
    m_x.clear();
    return *this;
  }

  /// Solve (this is THE instruction everyone is waiting for!)
  virtual bool solve() = 0;

  // -- Linear system checks

  /// Checks whether the linear system matrix is empty
  bool empty() { return !m_A.size(); }

  /// Returns the generic size of the system matrix
  const idx_t& size() const { return m_A.m_size; }

  /// Returns the specific dimension size of the system matrix
  size_t size(const size_t& d) const {
      return (d==0? m_A.m_size.i :
             (d==1? m_A.m_size.j :
             (d==2? m_b.size(1) : std::numeric_limits< size_t >::max())));
    }

  // a private friend, how promiscuous of you
 private:
  template< typename aT, typename aMATRIX, typename aVECTOR, typename aINDEX_HIERARCHY >
  friend std::ostream& operator<< (std::ostream&, const linearsystem< aT, aMATRIX, aVECTOR, aINDEX_HIERARCHY >&);

  // storage

 protected:
  /// linear system matrix and vectors (containers)
  index_hierarchy_t m_idx;
  matrix_t m_A;
  vector_t m_b;
  vector_t m_x;
};


/* -- linear system output -------------------------------------------------- */

template< typename T, typename MATRIX, typename VECTOR, typename INDEX_HIERARCHY >
std::ostream& operator<< (std::ostream& o, const linearsystem< T, MATRIX, VECTOR, INDEX_HIERARCHY >& lss)
{
  o << "linearsystem: A: "; lss.m_A.print(o); o << std::endl;
  o << "linearsystem: b: "; lss.m_b.print(o); o << std::endl;
  o << "linearsystem: x: "; lss.m_x.print(o); o << std::endl;
  return o;//.print(o);
}


#endif

