#include "lasCSR.h"
#include <apfFieldOp.h>
#include <cassert>
namespace las
{
  bool needNonzero(int * rws, int rw, int * cls, int cl)
  {
    bool result = true;
    // negative rw and cl correspond to fixed dofs and won't add a nonzero
    if(rw < 0 || cl < 0)
      result = false;
    else
    {
      for(int ii = rws[rw]; ii < rws[rw+1]; ii++)
      {
        if(cls[ii-1] == (cl+1))
        {
          result = false;
          break;
        }
      }
    }
    return result;
  }
  int findLocation(int * rws, int rw, int * cls, int cl, int neq)
  {
    int loc = 0;
    assert(rw < neq);
    for(loc = rws[rw]; (loc < rws[rw+1]) && (cls[loc-1] < cl+1); ++loc) {}
    return loc;
  }
  bool addNonzero(int * rws, int rw, int * cls, int cl, int neq)
  {
    bool result = false;
    if((result = needNonzero(rws,rw,cls,cl)))
    {
      int ii = findLocation(rws,rw,cls,cl,neq) - 1;
      for(int jj = rw+1; jj < neq+1; ++jj)
        rws[jj]++;
      for(int jj = rws[neq] - 2; jj > ii; --jj)
        cls[jj] = cls[jj-1];
      cls[ii] = cl+1;
    }
    return result;
  }
  class CSRBuilder
  {
  protected:
    apf::Numbering * nm;
    apf::MeshEntity * ment;
    int nedofs;
    int ndofs;
    int nnz;
    std::vector<int> rws;
    std::vector<int> cls;
  public:
    CSRBuilder(apf::Numbering * n, int nd)
      : nm(n)
      , ment(NULL)
      , nedofs(0)
      , ndofs(nd)
      , nnz(0)
      , rws(ndofs+1)
      , cls(ndofs*ndofs)
    {
      assert(nm);
      rws.assign(ndofs+1,1);
      cls.assign(ndofs*ndofs,0);
    }
    void apply()
    {
      apf::Mesh * msh = apf::getMesh(apf::getField(nm));
      // find the highest dimension with mesh entities
      int dim = -1;
      for(dim = 3; dim >= 0; --dim)
        if(msh->count(dim) != 0)
          break;
      apf::MeshEntity * ent = NULL;
      apf::MeshIterator * it = NULL;
      for(it = msh->begin(dim); (ent = msh->iterate(it));)
      {
        apf::NewArray<int> dofs;
        int nedofs = apf::getElementNumbers(nm, ent, dofs);
        for(int ii = 0; ii < nedofs; ii++)
          for(int jj = 0; jj < nedofs; jj++)
            if(addNonzero(&rws[0],dofs[ii],&cls[0],dofs[jj],ndofs))
              nnz++;
      }
    }
    CSR * finalize()
    {
      return new CSR(ndofs,nnz,&rws[0],&cls[0]);
    }
  };
  CSR * createCSR(apf::Numbering * num, int ndofs)
  {
    CSRBuilder bldr(num,ndofs);
    bldr.apply();
    return bldr.finalize();
  }
  CSR * csrFromFull(double * mat, int rws, int cls)
  {
    int nnz = 0;
    std::vector<int> rwb(rws+1,1);
    std::vector<int> clb(rws*rws,0.0);
    for(int ii = 0; ii < rws; ii++)
      for(int jj = 0; jj < cls; jj++)
        if(mat[ii*cls + jj] != 0.0 && addNonzero(&rwb[0],ii,&clb[0],jj,rws))
          nnz++;
    rwb[rws+1] = nnz;
    CSR * rslt = new CSR(rws,nnz,&rwb[0],&clb[0]);
    return rslt;
  }
}
