#include "lasSparskitExterns.h"
#ifdef BGQ
void ilut_(int *n,double a[],int ja[],int ia[],int *lfil,double *droptol,double *alu,int *jlu,int *ju,int *iwk,double *w,int *jw,int *ierr)
{
  ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr);
}
void lusol_(int *n,double a[],double *x,double *alu,int *jlu,int *ju)
{
  lusol(n,a,x,alu,jlu,ju);
}
void pgmres_(int *n,int *imk,double *rhs,double sol[],double *vv,double *eps,int *maxits,int *iout,double aa[],int *ja,int *ia,double *alu,int *jlu,int *ju,int *ierr)
{
  pgmres(n,imk,rhs,sol,vv,eps,maxits,iout,aa,ja,ia,alu,jlu,ju,ierr);
}
void amux_(int *n,double x[],double y[],double a[],int ja[],int ia[])
{
  amux(n,x,y,a,ja,ia);
}
#endif
namespace las
{
  SparskitBuffers::SparskitBuffers(int num_dofs) :
    heuristic_length(num_dofs * sqrt(num_dofs) * 100),
    int_work_array(2*num_dofs),
    rows(num_dofs),
    cols(heuristic_length),
    double_work_array(num_dofs),
    matrix(heuristic_length)
  { }
  void SparskitBuffers::zero()
  {
    int_work_array.assign(int_work_array.size(),0.0);
    rows.assign(rows.size(),0.0);
    cols.assign(heuristic_length,0.0);
    double_work_array.assign(double_work_array.size(),0.0);
    matrix.assign(heuristic_length,0.0);
  }
};
