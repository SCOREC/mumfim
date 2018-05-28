#include "io.h"
#include <cassert>
#include <cstring>
#include <sstream>
#include <vector>
std::istream & operator>>(std::istream & in,
                          apf::DynamicMatrix & m)
{
  std::vector<std::vector<double> > tmp;
  std::string ln;
  int ln_nm = 0;
  int ln_lng = 0;
  int trms = 0;
  while(std::getline(in,ln))
  {
    trms = 0;
    std::vector<double> vls;
    std::stringstream sin;
    sin.str(ln);
    double d = 0.0;
    while(sin >> d)
    {
      vls.push_back(d);
      trms++;
    }
    if(ln_nm == 0)
      ln_lng = trms;
    else
      assert(trms == ln_lng);
    tmp.push_back(vls);
    ln_nm++;
  }
  m.setSize(ln_nm,trms);
  for(int rw = 0; rw < ln_nm; ++rw)
    for(int cl = 0; cl < trms; ++cl)
      m(rw,cl) = tmp[rw][cl];
  return in;
}
std::istream & operator>>(std::istream & in,
                           apf::DynamicVector & v)
{
  std::vector<double> tmp;
  double d = 0.0;
  while(in >> d)
    tmp.push_back(d);
  unsigned sz = tmp.size();
  v.setSize(sz);
  memcpy(&v[0],&tmp[0],sizeof(double)*sz);
  return in;
}
bool close(double a, double b, double eps = 1e-12)
{
  double df = fabs(a-b);
  double fa = fabs(a);
  double fb = fabs(b);
  double ref = ( (fa < fb ? fb : fa ) * eps );
  bool cls = df <= ref;
  return cls;
}
bool operator==(const apf::DynamicMatrix & a,
                const apf::DynamicMatrix & b)
{
  bool eq = true;
  int a_rw_cnt = a.getRows();
  int a_cl_cnt = a.getColumns();
  int b_rw_cnt = b.getRows();
  int b_cl_cnt = b.getColumns();
  assert(a_rw_cnt == b_rw_cnt && a_cl_cnt == b_cl_cnt);
  for(int rw = 0; rw < a_rw_cnt; ++rw)
    for(int cl = 0; cl < a_cl_cnt; ++cl)
      if(!close(a(rw,cl),b(rw,cl),1e-4))
      {
        eq = false;
        goto finish;
      }
finish:
  return eq;
}
void rowPermute(const apf::DynamicMatrix & a,
                int * prmt,
                apf::DynamicMatrix & b)
{
  int rws = a.getRows();
  b.setSize(rws,a.getColumns());
  for(int rw = 0; rw < rws; ++rw)
  {
    apf::DynamicVector arw;
    a.getRow(rw,arw);
    b.setRow(prmt[rw],arw);
  }
}
void colPermute(const apf::DynamicMatrix & a,
                int * prmt,
                apf::DynamicMatrix & b)
{
  int cls = a.getColumns();
  b.setSize(a.getRows(),cls);
  for(int cl = 0; cl < cls; ++cl)
  {
    apf::DynamicVector acl;
    a.getColumn(cl,acl);
    b.setColumn(prmt[cl],acl);
  }
}
