#include "CSR.h"
#include <apfField.h>
#include <cassert>
namespace bio
{
  bool needNonzero(int * rws, int rw, int * cls, int cl)
  {
    bool result = true;
    for(int ii = rws[rw]; ii < rws[rw+1]; ii++)
    {
      if(cls[ii-1] == (cl+1))
      {
	result = false;
	break;
      }
    }
    return result;
  }
  int findLocation(int * rws, int rw, int * cls, int cl)
  {
    int loc;
    for(loc = rws[rw]; (loc < rws[rw+1]) && (cls[loc-1] < cl + 1) ; loc++) {}
    return loc;
  }
  bool addNonzero(int * rws, int rw, int * cls, int cl, int neq)
  {
    bool result = false;
    if((result = needNonzero(rws,rw,cls,cl)))
    {
      int ii = findLocation(rws,rw,cls,cl) - 1;
      for(int jj = rw + 1; jj <= neq; jj++)
	rws[jj]++;
      for(int jj = rws[neq] - 2; jj > ii; jj--)
	cls[jj] = cls[jj - 1];
      cls[ii] = cl + 1;
    }
    return result;
  }
  class CSRBuilder : public apf::FieldOp
  {
  protected:
    apf::Numbering * nm;
    apf::MeshEntity * ment;
    int nedofs;
    int ndofs;
    int nnz;
    int * rws;
    int * cls;
  public:
    CSRBuilder(apf::Numbering * n, int nd)
      : nm(n)
      , ment(NULL)
      , nedofs(0)
      , ndofs(nd)
      , nnz(0)
      , rws(NULL)
      , cls(NULL)
    {
      assert(nm);
      rws = new int[ndofs]();
      cls = new int[ndofs*ndofs]();
    }
    ~CSRBuilder()
    {
      // todo (h) Bill : memory leak? but double-free error if these are here...
      //delete [] rws;
      //delete [] cls;
    }
    bool inEntity(apf::MeshEntity * me)
    {
      ment = me;
    }
    void outEntity() {}
    void atNode(int nde)
    {
      apf::NewArray<int> dofs;
      int nedofs = apf::getElementNumbers(nm,ment,dofs);
      for(int ii = 0; ii < nedofs; ii++)
	for(int jj = 0; jj < nedofs; jj++)
	  if(addNonzero(rws,dofs[ii],cls,dofs[jj],ndofs))
	    nnz++;
    }
    CSR * finalize()
    {
      return new CSR(ndofs,nnz,rws,cls);
    }
  };

  CSR * createCSR(apf::Numbering * num, int ndofs)
  {
    CSRBuilder bldr(num,ndofs);
    bldr.apply(apf::getField(num));
    return bldr.finalize();
  }
}
