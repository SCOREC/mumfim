#ifndef BIO_INIT_GUESS_H_
#define BIO_INIT_GUESS_H_

namespace bio
{
  // basically just shape function values...
  void tricubicInterpCoefs(const apf::Vector3 & crd,
			   const apf::Vector3 & dms,
			   const apf::Vector3 & trns,
			   double (&cfs)[8])
  {
    double x = crd[0]*dms[0] + trns[0];
    double y = crd[1]*dms[1] + trns[1];
    double z = crd[2]*dms[2] + trns[2];
    
    double mx = dms[0] - x;
    double my = dms[1] - y;
    double mz = dms[2] - z;
    
    cfs[0] = mx * my * mz;
    cfs[1] =  x * my * mz;
    cfs[2] = mx *  y * mz;
    cfs[3] =  x *  y * mz;
    cfs[4] = mx * my *  z;
    cfs[5] =  x * my *  z;
    cfs[6] = mx *  y *  z;
    cfs[7] =  x *  y *  z;
  }
};

#endif
