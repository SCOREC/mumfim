#ifndef BIO_SOLVER_H_
#define BIO_SOLVER_H_

namespace bio
{

  void RVEanalysis(RVE * rve, FiberNetwork * fn, double epsilon)
  {
    SparskitLinearSystem ls; // would be better to pass this in...
    SparskitLinearSolver solver;
    LSOps * ops = ls.getOps();

    int iter = 0;
    double ni = 0.0;
    double nim1 = 0.0;
    Vector * sol = ls.getSolution();
    Vector * dmp = NULL;
    do
    {
      bool damping = false;
      do
      {
	assembleLinearSystem(fn->getMesh(),
			     fn->getNumbering(),
			     1,
			     fn,
			     &ls);
	applyRVEForceBC(&ls,rve,fn);
	ni = ops->norm(ls.getVector());
	if(iter > 0 && (ni / nim1) > 100.0 )
	{
	  damping = true;
	  //dmp = AXPY(-0.5,sol,NULL);
	  //Applysolution(fn->getNumbering(),dmp).apply(fn->getIncrementalDispField());
	}
	else
	  damping = false;
      } while( !damping )

      solver.solve(ls);     
      sol = ls.getSolution();
      ApplySolution(fn->getNumbering(),sol).apply(fn->getIncrementalDispField());
      nim1 = ni;
      iter++;
    } while ( ni > epsilon )
  }
  
  void multiscaleRVEAnalysis(MacroInfo * macro,
			     RVE * rve,
			     FiberNetwork * fn)
  {
    // use macroscale information to setup initial displacements

    RVEAnalysis(rve,fn,1e-12);

    // compute the macroscale stress terms
  }
}

#endif
