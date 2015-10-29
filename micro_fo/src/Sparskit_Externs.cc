#include "Sparskit_Externs.h"

namespace Biotissue
{
  SparskitBuffers::SparskitBuffers(int num_dofs) :
    heuristic_length(num_dofs * sqrt(num_dofs)),
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
