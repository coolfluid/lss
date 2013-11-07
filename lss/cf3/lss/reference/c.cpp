
#include <iostream>
#include <limits>
#include <vector>


#if 0
struct idx_t
{
  size_t i, j;
  idx_t(const size_t& _i=std::numeric_limits< size_t >::max(),
        const size_t& _j=std::numeric_limits< size_t >::max()) : i(_i), j(_j) {}

  bool operator<  (const idx_t& _other) const { return (i<_other.i? true : i>_other.i? false : (j<_other.j)); }
  bool operator>  (const idx_t& _other) const { return idx_t(_other)<*this; }
  bool operator== (const idx_t& _other) const { return i==_other.i && j==_other.j; }
  bool operator!= (const idx_t& _other) const { return i!=_other.i || j!=_other.j; }

  idx_t& invalidate() { return (*this = idx_t()); }
  bool is_valid_size()  const { return operator>(idx_t(0,0)) && operator<(idx_t()); }
  bool is_square_size() const { return i==j; }
  bool is_diagonal()    const { return is_square_size(); }
};
#endif


void log(const std::string& p, const int& i) {
  std::cout << p << i << std::endl;
}


struct index_conversion_t
{
  virtual int& apply(int&) const = 0; 
}; 


struct index_hierarchy_t_end : index_conversion_t {                                                                   void setup(const index_hierarchy_t_end& t) { log(".index_hierarchy_t_end.setup.",0); } int& apply(int& i) const { log(".index_hierarchy_t_end.apply.",i); return i; } }; 
struct test_mul_t            : index_conversion_t { test_mul_t(const int& _in) { log(".test_mul_t.construct.",_in); } void setup(const test_mul_t& t)            { log(".test_mul_t.setup.",0);            } int& apply(int& i) const { log(".test_mul_t.apply.",i);            return i*=x; } test_mul_t() : x(1) {} int x; }; 
struct test_sum_t            : index_conversion_t { test_sum_t(const int& _in) { log(".test_sum_t.construct.",_in); } void setup(const test_sum_t& t)            { log(".test_sum_t.setup.",0);            } int& apply(int& i) const { log(".test_sum_t.apply.",i);            return i+=x; } test_sum_t() : x(0) {} int x; };


template< typename Ta, typename Tb=index_hierarchy_t_end >
struct index_hierarchy_t : Tb
{
  // copy constructor 
  index_hierarchy_t() {}  // (necessary to initialize a linearsystem without a setup indices)



  Tb& operator()(const Ta _Oa) {
    std::cout << "s -> " << Oa.x << std::endl;
    return Ob;
  }


  index_hierarchy_t(const Ta& _Oa, const Tb& _Ob=index_hierarchy_t_end()) {}

  void setup(const index_hierarchy_t& t) { Ob.setup(t.Ob); Oa.setup(t.Oa); }  
  int& apply(int& i) const         { return Oa.apply(Ob.apply(i)); }

  Ta Oa;
  Tb Ob;
};


struct index_multidomain_singlepointcloud_t : index_conversion_t {
};




template<
  typename T,
  typename INDEX=index_hierarchy_t< index_hierarchy_t_end >
  >
struct linearsystem
{
  typedef INDEX index_t;
  std::vector< T > a;
  index_t idx;

  void setup(index_t& _other) { idx.setup(_other); }
};


main() {

  int i = 42;

  typedef linearsystem< double,
            index_hierarchy_t< test_mul_t,
            index_hierarchy_t< test_mul_t,
            index_hierarchy_t< test_sum_t,
            index_hierarchy_t< test_mul_t,
            index_hierarchy_t< test_sum_t,
            index_hierarchy_t< test_mul_t > > > > > >
              > lss_t;
  lss_t lss;

  lss.idx
        (test_mul_t(1))
        (test_mul_t(1))
        (test_sum_t(0))
        (test_mul_t(1))
        (test_sum_t(10))
        (test_mul_t(1));


  lss_t::index_t t;

  std::cout << "X1" << std::endl;
  //lss.setup(t("DomainA","DomainB")(10,5));
  //lss.setup(index_hierarchy_t< test_mul_t >("DomainA","DomainB"));
  std::cout << "X2" << std::endl;

  std::cout << "hello, " << lss.idx.apply(i) << "!" << std::endl;

}

