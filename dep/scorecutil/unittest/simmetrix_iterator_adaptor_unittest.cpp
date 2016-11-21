
struct LatticeModel;
// Ugly hack should go away on 1/9/2007 when Ottmar updates us.
struct LatticeModelIter {
  LatticeModelIter(LatticeModel*);
  LatticeModelIter(LatticeModelIter const&);
  unsigned char buff[1024]; // Maybe that's enough?
};
struct LatticeIter {
  LatticeIter(LatticeIter const&);
  unsigned char buff[1024]; // Maybe that's enough?
};


LatticeModelIter* LMIter_clone(LatticeModelIter* i) {
  return new LatticeModelIter(*i);
};
LatticeIter* LIter_clone(LatticeIter* i) {
  //  return new LatticeIter(*i);
};


#include <boost/test/minimal.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <vector>
#include <set>
class LatticeDomain {}; // Not sure why this is needed.
#include "beta/simmetrix_iterator_adaptor.hpp"
#include "mPoint.h"




using std::vector;
using std::set;
using std::cout;
using std::endl;
using SCOREC::Util::mPoint;
using SCOREC::Util::mVector;

namespace {
  mPoint const position(LIter i) { 
    mPoint result;
    LIter_position(i, data(result));
    return result;
  }

  const std::string name(LIter i) {
    return std::string(LIter_name(i));
  }
  
  unsigned int atomic_number(LIter i) {
    return LIter_AtomicNumber(i);
  }

  using SCOREC::Util::Simmetrix::lattice_iterator;
  using SCOREC::Util::Simmetrix::lattice_model_iterator;

  class scoped_lattice_model {
  public:
    explicit scoped_lattice_model(pLatticeModel lm) 
      : lattice_model_(lm) {}
    ~scoped_lattice_model() {
      LatticeModel_release(lattice_model_);
    }
    operator pLatticeModel() { return lattice_model_; }
  private:
    pLatticeModel lattice_model_;
  };

  // A RAII class for SimLattice 
  class SimLattice_instance {
  public:
    SimLattice_instance() {
      SimLattice_start();
    }
    ~SimLattice_instance() {
      SimLattice_stop();
    }
  };

  bool region_contains(pLattice l, const mPoint& p) {
    return (Lattice_isPointInDomain(l, const_cast<double*>(data(p))) 
	    ? 
	    true : false);
  }


  std::ostream& operator<<(std::ostream& os, LIter i) {
    os << name(i) << " (" << atomic_number(i) << ")\t" 
       << LIter_LatticeIndex1(i) << "/"
       << LIter_LatticeIndex2(i) << "/"
       << LIter_LatticeIndex3(i) << "/"
       << LIter_LatticeIndex4(i) << "\t"
       << position(i) << "\t" << static_cast<void*>(i);
    return os;
  }


}

int old_style_test();

int iterator_test();


int test_main(int, char**) {
  SimLattice_instance theSimLatticeInstance;

  old_style_test();
  cout << "=============\n";
  return iterator_test();
}


int old_style_test() {
    //**************************************************
  //                    Cubic Perovskite Lattice
  //**************************************************  
    pLatticeModel latModel = LatticeModel_new();
    pLatticeDomain d = LatticeDomain_newRect(0,0,0,3.8,3.8,3.8);
    pLattice lat = LatticeModel_newLattice(latModel,1,0,0,
                                   0,1,0,
                       3.795,0,0,
                       0,3.795,0,
                       0,0,3.795,d);
    Lattice_addBasisVector(lat,0,0,0,"Ca");
    Lattice_addBasisVector(lat,1.8975,1.8975,1.8975,"Ti");
    Lattice_addBasisVector(lat,0.,1.8975,1.8975,"O");
    Lattice_addBasisVector(lat,1.8975,0.,1.8975,"O");
    Lattice_addBasisVector(lat,1.8975,1.8975,0.,"O");
    double min[3], max[3];
    LatticeModel_bounds(latModel,min,max);


    LMIter iter = LatticeModel_Iter(latModel);
    while(lat = LMIter_next(iter))
    {
      cout << "New lattice. (" << min[0] << "," << min[1] << "," << min[2] << ") (" << max[0] << "," << max[1] << "," << max[2] << ")" << endl;
      LIter latIter = Lattice_Iter(lat, min, max);
      while (LIter_next(latIter))
      {
	double xyz[3];
	LIter_position(latIter, xyz);
	printf("%s %d/%d/%d/%d  %e %e %e\n",LIter_name(latIter),
	       LIter_LatticeIndex1(latIter),
	       LIter_LatticeIndex2(latIter),
	       LIter_LatticeIndex3(latIter),
	       LIter_LatticeIndex4(latIter),
	       xyz[0],xyz[1],xyz[2]);     
      }
      LIter_delete(latIter);
    }
    LMIter_delete(iter);
    LatticeModel_release(latModel);

  return 0;
}

template <typename ForwardIterator>
void cout_range(ForwardIterator i, ForwardIterator end) {
  for (; i != end; ++i) cout << *i << endl;
}



int iterator_test()
{
/*
  using  ::Util::Simmetrix;

  int error = 0;
    
  //**************************************************
  //                    Cubic Perovskite Lattice
  //**************************************************  
  scoped_lattice_model latModel(LatticeModel_new());
  //  pLatticeModel latModel = LatticeModel_new();
  pLattice lattice 
    = LatticeModel_newLattice(latModel,
			      1,0,0, // Orientation 1
			      0,1,0, // Orientation 2
			      3.795,    0,    0, // lat1
			          0,3.795,    0, // lat2
			          0,    0,3.795, // lat3
			      LatticeDomain_newRect(0,0,0,3.8,3.8,3.8));
  Lattice_addBasisVector(lattice,      0,      0,      0, "Ca");
  Lattice_addBasisVector(lattice, 1.8975, 1.8975, 1.8975, "Ti");
  Lattice_addBasisVector(lattice,    0.0, 1.8975, 1.8975, "O");
  Lattice_addBasisVector(lattice, 1.8975,    0.0, 1.8975, "O");
  Lattice_addBasisVector(lattice, 1.8975, 1.8975,    0.0, "O");
  double min[3], max[3];
  LatticeModel_bounds(latModel, min, max);
  cout << " Model bounds: (" << min[0] << "," << min[1] << "," << min[2] << ") (" 
       << max[0] << "," << max[1] << "," << max[2] << ")" << endl;  
  //  using SCOREC::Util::Simmetrix::lattice_range_getter;

  const lattice_model_iterator modelEnd = end(latModel);
  for (lattice_model_iterator lat = begin(latModel);
       lat != modelEnd;
       ++lat) {
       
    const lattice_iterator latEnd = end(*lat, min, max);
    for (lattice_iterator latIter = begin(*lat, min, max);
	 latIter != latEnd;
	 ++latIter) {
      cout << *latIter << endl;
    }
  }
  cout << "Iterating over all atoms.\n";
  all_lattice_iterator i(begin(latModel), end(latModel));
  const all_lattice_iterator all_end(end(latModel), end(latModel));
  for (; i != all_end; ++i) {
    cout << *i << endl;

    typedef lattice_model_nested_traits::position_type pos_type;
    mPoint nearMinPoint = position(*i) - 2.0*mVector(1.0, 1.0, 1.0);
    mPoint nearMaxPoint = position(*i) + 2.0*mVector(1.0, 1.0, 1.0); 
    pos_type nearMin, nearMax;
    std::copy(data(nearMinPoint), data(nearMinPoint) + 3, nearMin.begin());
    std::copy(data(nearMaxPoint), data(nearMaxPoint) + 3, nearMax.begin());
    const lattice_model_nested_traits range(nearMin, nearMax);

    all_lattice_iterator ineigh(begin(latModel), end(latModel), 
				range);
    all_lattice_iterator ineigh_end(end(latModel), end(latModel), 
				    range);
    for(; ineigh != ineigh_end; ++ineigh) {
      cout << "  near: " << *ineigh << endl;
      if (ineigh == i) cout << "\t\tIt's me, it's me!\n";
    }
    // Try boost filter itertor.

    LatticeModelIter testing(latModel);
    std::auto_ptr<LatticeModelIter> testing2(new LatticeModelIter(testing));




    all_lattice_iterator ineigh2(begin(latModel), end(latModel), 
				 range);
    using boost::cref;
    using boost::ref;
    using boost::lambda::_1;
    using boost::make_filter_iterator;
    cout << "=====\n";
    vector<int> asdf(3, 42);
    asdf[1] = 12;
    int test = 12;
    cout_range(make_filter_iterator(_1 != cref(test), asdf.begin(), asdf.end()),
	       make_filter_iterator(_1 != cref(test), asdf.end(),   asdf.end()));
	       //    cout_range(make_filter_iterator(_1 != cref(i), ineigh2,    ineigh_end),
	  make_filter_iterator(_1 != cref(i), ineigh_end, ineigh_end));
    cout << "=====\n";

  }
  return error;
*/
return 0;
}

