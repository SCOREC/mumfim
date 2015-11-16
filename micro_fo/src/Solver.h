#ifndef BIO_SOLVER_H_
#define BIO_SOLVER_H_

namespace bio
{

  void RVEanalysis(RVE * rve, FiberNetwork * fn)
  {
    SparskitLinearSystem ls; // would be better to pass this in...
    SparskitLinearSolver solver;
    LSOps * ops = ls.getOps();

    double ni = 0.0;
    double nim1 = 0.0;
    assembleLinearSystem(fn->getMesh(),
			 fn->getNumbering(),
			 1,
			 fn,
			 &ls);
    applyRVEForceBC(&ls,rve,fn);
    ni = ops->norm(ls.getVector());
    solver.solve(ls);
    Vector * sol = ls.getSolution();
    ApplySolution(fn->getNumbering(),sol).apply(fn->getIncrementalDispField());
  }
  
  void multiscaleRVEAnalysis(MacroInfo * macro,
			     RVE * rve,
			     FiberNetwork * fn)
  {
    // use macroscale information to setup initial displacements

    RVEAnalysis(rve,fn);

    // compute the macroscale stress terms
  }
}

#endif
