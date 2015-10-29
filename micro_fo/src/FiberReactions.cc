#include "FiberReactions.h"

namespace bio
{
  double NonLinearReaction::f(const Element & e)
  {
    double result = 0.0;
    double lngth = 1.0; // calcFiberLength(e);
    double olngth = e.orig_len;
    double length_ratio = lnght / olngth;
    double tension = E * fiber_are / B;
    double grn_strn = 0.5 * (length_ratio * length_ratio - 1.0);
    double expBeps = exp(B * grn_strn);
    result = tension * (expBeps - 1.0);
    return result;
  }

  double NonLinearReaction::df_dl(const Element & e)
  {
    double result = 0.0;
    double lngth = 1.0; // calcFiberLength(e);
    double olngth = e.orig_len;
    double lngth_rtio = lngth / olngth;
    double tension = E * fbr_ar;
    double grn_strn = 0.5 * (lngth_rtio * lngth_rtio - 1.0);
    double expBeps = exp(B * grn_strn);
    result = ((fbr_ar * E * lngth) / (olngth * olngth)) * expBeps;
    return result;
  }

  double LinearReaction::f(const Element & e)
  {
    double result = 0.0;
    double lngth = 1.0; // calcFiberLength(e);
    double olngth = e.orig_len;
    double lngth_rtio = lngth / olngth;
    result = (E * fbr_ar) / (lngth_rtio - 1.0);
    return result;
  }

  double LinearReaction::df_dl(const Element & e)
  {
    double result = 0.0;
    double lngth = 1.0; // calcFiberLength(e);
    double olngth = e.orig_len;
    double lngth_rtio = lngth / olngth;
    result = (E * fbr_ar) / olngth;
    return result;
  }
}
