#include "RepresentVolElem.h"
#include "globals.h"
#include "Sparskit_Externs.h"
#include "RVE_Util.h"
#include "Util.h"
#include <cstring> //memset
#include <amsiUtil.h>
namespace bio
{
  int MicroFO::Solver()
  {
    int result = 0;
    int num_dofs = fiber_network->numDofs();
    int num_elements = fiber_network->numElements();
    //int num_nodes = fiber_network->numNodes();
    int step = 0;
    int cont = 0;
    double tol = 1e-6; // tolerance for dropping relatively small off-diagonal terms during ilu factorization
    double tolerance = 1e-8; // convergence tolerance
    double current_norm = std::numeric_limits<double>::max();
    double previous_norm = std::numeric_limits<double>::max();
    double relative_norm = std::numeric_limits<double>::max();
    relative_norm = 0.0; previous_norm = 1.0;
    int ierr = 0;
    double * solution = new double[num_dofs]();
    double * fib_str = new double[num_elements]();
    double * dfdl = new double[num_elements]();
    // debug vars
    std::vector<double> mat;
    // Calculates the fiber length
    std::vector<double> len(num_elements);
    calcFiberLengths(*fiber_network,len);
    // Calculates the fiber force
    calc_fiber_str(&len[0],
                   &fib_str[0],
                   &dfdl[0]);
    // Assembles the Jacobian matrix
    calc_precond(&matrix[0],
                 &len[0],
                 &fib_str[0],
                 &dfdl[0],
                 false);
    // Calculates the x,y,z conponents of the fiber force
    calc_force_vector(&len[0],
                      &fib_str[0]); // resets the force vector to 0.0 before accumulating values...
    // Calculates boundary conditions of the microscopic problem
    force_vector_bcs();
    // Calculates the norm of the microscopic residuals
    current_norm = calc_norm(&force_vector[0],
                             num_dofs);
    buffers->zero();
    // begin timing solver
    double start_time = MPI_Wtime();
    // Newton iterations begin
    do // while (norm > tolerance)
    {
      //std::cout<<"number of iterations in solver_gmres_df.cc"<<std::endl;
      // Preconditioner
      int buffer_length = buffers->matrixLength();
      ilut_(&num_dofs,
            &matrix[0],
            sparse_structure->getCols(), // might need to switch these...
            sparse_structure->getRows(),
            &num_dofs,
            &tol,
            buffers->matrixBuffer(),
            buffers->colsBuffer(),
            buffers->rowsBuffer(),
            &buffer_length,
            buffers->doubleWorkBuffer(),
            buffers->intWorkBuffer(),
            &ierr);
      if(ierr != 0)
      {
        std::cerr << "ERROR: ilut_ returned error code " << ierr << std::endl;
        result++;
      }
      // Solver
      lusol_(&num_dofs,
             &force_vector[0],
             &solution[0],
             buffers->matrixBuffer(),
             buffers->colsBuffer(),
             buffers->rowsBuffer());
      for(int ii = 0; ii < num_dofs; ii++)
      {
        /*
          if(fabs(solution[ii]) > 1.0)
          std::cerr << "solution at position " << ii << " is " << solution[ii] << std::endl;
        */
        coordinate_vector[ii] += solution[ii];
      }
      double df = 1.0;
      int ncv = 0;
      do
      {
        cont = 0;
        update_nodes(); //updates fiber_network->nodes with values in coordinate_vector
        // fiber length determined from fiber_network->nodes.
        calcFiberLengths(*fiber_network,len);
        // fiber force determined from len calculated in calc_fiber_len function.
        calc_fiber_str(&len[0],
                       &fib_str[0],
                       &dfdl[0]);
        //force vector determined from combination of fiber_network->nodes and fib_str determined from calc_fiber_str function.
        calc_force_vector(&len[0],
                          &fib_str[0]);
        // bcs imposed by setting terms in force_vector that are on boundary to 0.
        force_vector_bcs();
        // calc_norm determines L2 norm of force vector.
        current_norm = calc_norm(&force_vector[0],num_dofs);
        // when the norm increases drastically; when norm is inf
        // or nan; when it keeps oscillating even after a lot of iterations
        // (Damping)
        // checking for NaN (f != f is only true if f=NaN)
        if( current_norm != current_norm || ( (current_norm/previous_norm) > 100.0) || (ncv == 1) )
        {
          std::cout << "damping!" << std::endl;
          df = df - df/2;
          for(int ii = 0; ii < num_dofs; ii++)
          {
            solution[ii] = solution[ii] / 2;
            coordinate_vector[ii] -= solution[ii];
          }
          cont = 1;
          ncv = 0;
        }
      } while(cont == 1);
      //relative_norm = current_norm;
      if (current_norm != 0)
	relative_norm = fabs(previous_norm - current_norm);
      else
	relative_norm = current_norm;
      
      /*
        AMSI_DEBUG
        (
        std::cout << "Micro relative norm: " << relative_norm << " < " << tolerance << std::endl;
        std::cout << "Micro current norm: " << current_norm << std::endl;
        std::cout << "Micro previous norm: " << previous_norm << std::endl;
        )
      */
      previous_norm = current_norm;
      calc_precond(&matrix[0],
                   &len[0],
                   &fib_str[0],
                   &dfdl[0],
                   false);
      step++;
      if((step == 50)||(step == 100)||(step == 150)||(step == 200))
        ncv = 1; //Damping is required when ncv = 1.
      if(step == 60)
      {
        tolerance = tolerance * 10.0;
        AMSI_DEBUG(std::cout<<"step = "<<step<<", relative_norm = "<<relative_norm<<std::endl);
      }
      else if(step == 80)
      {
        tolerance = tolerance * 10.0;
        AMSI_DEBUG(std::cout<<"step = "<<step<<", relative_norm = "<<relative_norm<<std::endl);
      }
      else if(step == 100)
      {
        tolerance = tolerance * 10.0;
        std::cout<<"step = "<<step<<", relative_norm = "<<relative_norm<<std::endl;
      }
      else if(step == 120)
      {
        tolerance = tolerance * 10.0;
        AMSI_DEBUG(std::cout<<"step = "<<step<<", relative_norm = "<<relative_norm<<std::endl);
      }
      else if(step == 130)
      {
	std::cout<<"step = "<<step<<", relative_norm = "<<relative_norm<<std::endl;
        std::cerr << "Warning: unusual number of newton iterations in micro_fo! step = "<<step<<", relative_norm = "<<relative_norm <<std::endl;
        break;
      }
    } while (relative_norm > tolerance);
    double end_time = MPI_Wtime();
    // Newton iteration ended. The network got a new equilibrium position and the fiber length and fiber force are updated
    calcFiberLengths(*fiber_network,len); //added in to match Minnesota code.
    calc_fiber_str(&len[0],
                   &fib_str[0],
                   &dfdl[0]);
    calc_force_vector(&len[0],
                      &fib_str[0]);
//    calc_mean_fiber_stretch_omega();
    // Calculation of the derivative dSdy that is used in function calc_femjacob_newmethod for the calculation of the derivative dSdx
    calc_precond(&matrix[0],
                 &len[0],
                 &fib_str[0],
                 &dfdl[0],
                 true);
    // update_nodes();
    rve_iterations.push_back(step);
    rve_timing.push_back(end_time - start_time);
    delete [] solution;
    delete [] fib_str;
    delete [] dfdl;
    return result;
  }
  void MicroFO::force_vector_bcs()
  {
    int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::ALL);
    for(int ii = 0; ii < num_boundary_nodes; ii++)
    {
      //  Boundary conditions for the microscopic problem
      int node = fiber_network->boundaryNode(FiberNetwork::ALL,ii);
      force_vector[node*3] = 0.0;
      force_vector[node*3+1] = 0.0;
      force_vector[node*3+2] = 0.0;
    }
    SupportFiberNetwork * spfn = dynamic_cast<SupportFiberNetwork*>(fiber_network);
    if (spfn)
    {
      int num_support_nodes = spfn->numSupportNodes();
      for(int ii = 0; ii < num_support_nodes; ii++)
      {
        int node = spfn->supportNode(ii);
        force_vector[node*3] = 0.0;
        force_vector[node*3+1] = 0.0;
        force_vector[node*3+2] = 0.0;
      }
    }
  }
  void MicroFO::calc_precond(double * matrix,
                             double * lengths,
                             double * fib_str,
                             double * dfdl,
                             bool calc_dSdy)
  {
    // Assembly of the Jacobian matrix for the microscopic problem
    // Also the derivative dSdy is also calculated and will be used in the function calc_femjacob_newmethod
    int num_dofs = fiber_network->numDofs();
    int num_elements = fiber_network->numElements();
    memset(matrix,0,sparse_structure->numNonzeros()*sizeof(double));
    for(int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fiber_network->element(ii);
      double dfcdx2 = gfunc(ii,1,lengths,fib_str,dfdl);
      double dfcdy2 = gfunc(ii,2,lengths,fib_str,dfdl);
      double dfcdz2 = gfunc(ii,5,lengths,fib_str,dfdl);
      double dfsdx2 = dfcdy2;
      double dfsdy2 = gfunc(ii,4,lengths,fib_str,dfdl);
      double dfsdz2 = gfunc(ii,6,lengths,fib_str,dfdl);
      double dfzdx2 = dfcdz2;
      double dfzdy2 = dfsdz2;
      double dfzdz2 = gfunc(ii,7,lengths,fib_str,dfdl);
      // Assembly of the Jacobian Matrix
      matrix[e.Ke[0]]  += -dfcdx2;
      matrix[e.Ke[1]]  += -dfcdy2;
      matrix[e.Ke[2]]  += -dfcdz2;
      matrix[e.Ke[3]]  +=  dfcdx2;
      matrix[e.Ke[4]]  +=  dfcdy2;
      matrix[e.Ke[5]]  +=  dfcdz2;
      matrix[e.Ke[6]]  += -dfsdx2;
      matrix[e.Ke[7]]  += -dfsdy2;
      matrix[e.Ke[8]]  += -dfsdz2;
      matrix[e.Ke[9]]  +=  dfsdx2;
      matrix[e.Ke[10]] +=  dfsdy2;
      matrix[e.Ke[11]] +=  dfsdz2;
      matrix[e.Ke[12]] += -dfzdx2;
      matrix[e.Ke[13]] += -dfzdy2;
      matrix[e.Ke[14]] += -dfzdz2;
      matrix[e.Ke[15]] +=  dfzdx2;
      matrix[e.Ke[16]] +=  dfzdy2;
      matrix[e.Ke[17]] +=  dfzdz2;
      matrix[e.Ke[18]] +=  dfcdx2;
      matrix[e.Ke[19]] +=  dfcdy2;
      matrix[e.Ke[20]] +=  dfcdz2;
      matrix[e.Ke[21]] += -dfcdx2;
      matrix[e.Ke[22]] += -dfcdy2;
      matrix[e.Ke[23]] += -dfcdz2;
      matrix[e.Ke[24]] +=  dfsdx2;
      matrix[e.Ke[25]] +=  dfsdy2;
      matrix[e.Ke[26]] +=  dfsdz2;
      matrix[e.Ke[27]] += -dfsdx2;
      matrix[e.Ke[28]] += -dfsdy2;
      matrix[e.Ke[29]] += -dfsdz2;
      matrix[e.Ke[30]] +=  dfzdx2;
      matrix[e.Ke[31]] +=  dfzdy2;
      matrix[e.Ke[32]] +=  dfzdz2;
      matrix[e.Ke[33]] += -dfzdx2;
      matrix[e.Ke[34]] += -dfzdy2;
      matrix[e.Ke[35]] += -dfzdz2;
    }
    std::vector<int> con;
    if(!calc_dSdy)
    {
      int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::ALL);
      // Boundary Conditions for the microscopic problem
      for(int ii = 0; ii < num_boundary_nodes; ii++)
      {
        int node = fiber_network->boundaryNode(FiberNetwork::ALL,ii);
        con.push_back(node);
        //con[0] = node;
        //int s = 1;
        // collect all nodes which are attached to the current boundary node
        for(int jj = 0; jj < num_elements; jj++)
        {
          const Element & e = fiber_network->element(jj);
          if(e.node1_id == node)
            con.push_back(e.node2_id);
          if(e.node2_id == node)
            con.push_back(e.node1_id);
        }
        for(std::vector<int>::iterator iter = con.begin(), iterend = con.end();
            iter != iterend; iter++)
        {
          for(int nn = 0; nn < 3; nn++)
          {
            int jj = 3 * (*iter) + nn;
            int location = sparse_structure->sparseLocation(jj,node*3);
            matrix[location] = 0.0;
            location = sparse_structure->sparseLocation(jj,node*3+1);
            matrix[location] = 0.0;
            location = sparse_structure->sparseLocation(jj,node*3+2);
            matrix[location] = 0.0;
          }
        }
        matrix[sparse_structure->sparseLocation(node*3  ,node*3  )] = 1.0;
        matrix[sparse_structure->sparseLocation(node*3+1,node*3+1)] = 1.0;
        matrix[sparse_structure->sparseLocation(node*3+2,node*3+2)] = 1.0;
        con.erase(con.begin(),con.end());
      }
      /* Affinely deform supported fibers*/
      SupportFiberNetwork * spfn = dynamic_cast<SupportFiberNetwork*>(fiber_network);
      if (spfn)
      {
        int num_support_nodes = spfn->numSupportNodes();
        for (int ii = 0; ii < num_support_nodes; ii++)
        {
          //int node = supportNodes[ii];
          int node = spfn->supportNode(ii);
//      std::cout<<"support node is "<<node<<std::endl;
          con.push_back(node);
          // collect all nodes which are attached to a support node
          for(int jj = 0; jj < num_elements; jj++)
          {
            const Element & e = fiber_network->element(jj);
            if(e.node1_id == node)
              con.push_back(e.node2_id);
            if(e.node2_id == node)
              con.push_back(e.node1_id);
          }
          for(std::vector<int>::iterator iter = con.begin(), iterend = con.end();
              iter != iterend; iter++)
          {
            for(int nn = 0; nn < 3; nn++)
            {
              int jj = 3 * (*iter) + nn;
              int location = sparse_structure->sparseLocation(jj,node*3);
              matrix[location] = 0.0;
              location = sparse_structure->sparseLocation(jj,node*3+1);
              matrix[location] = 0.0;
              location = sparse_structure->sparseLocation(jj,node*3+2);
              matrix[location] = 0.0;
            }
          }
          matrix[sparse_structure->sparseLocation(node*3  ,node*3  )] = 1.0;
          matrix[sparse_structure->sparseLocation(node*3+1,node*3+1)] = 1.0;
          matrix[sparse_structure->sparseLocation(node*3+2,node*3+2)] = 1.0;
          con.erase(con.begin(),con.end());
        } // End num_support_nodes loop.
      } // End If (SPECIFY_FIBER_TYPE) loop.
    }
    else
    {
      ttdSdy.assign(num_dofs*6,0);
      // Calculation of the derivative dSdy which is used later for the calculation
      // of the derivative dSdx in the function calc_femjacob_newmethod
      int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::ALL);
      for(int ii = 0; ii < num_boundary_nodes; ii++)
      {
        int node = fiber_network->boundaryNode(FiberNetwork::ALL,ii);
        Node n = fiber_network->node(node);
        double xyz[3];
        xyz[0] = n.x;
        xyz[1] = n.y;
        xyz[2] = n.z;
        con.push_back(node);
        for(int jj = 0; jj < num_elements; jj++)
        {
          const Element & e = fiber_network->element(jj);
          if(e.node1_id == node)
            con.push_back(e.node2_id);
          if(e.node2_id == node)
            con.push_back(e.node1_id);
        }
        for(std::vector<int>::iterator iter = con.begin(), iterend = con.end();
            iter != iterend; iter++)
        {
          for(int nn = 0; nn < 3; nn++)
          {
            int kk = 3 * (*iter) + nn;
            int location = sparse_structure->sparseLocation(kk,node * 3);
            ttdSdy[               kk] -= matrix[location] * xyz[0];
            ttdSdy[    num_dofs + kk] -= matrix[location] * xyz[1] / 2;
            ttdSdy[2 * num_dofs + kk] -= matrix[location] * xyz[2] / 2;
            if(kk == node * 3)
              ttdSdy[               kk] += force_vector[node * 3];
            if(kk == node * 3 + 1)
              ttdSdy[    num_dofs + kk] += force_vector[node * 3] / 2;
            if(kk == node * 3 + 2)
              ttdSdy[2 * num_dofs + kk] += force_vector[node * 3] / 2;
            matrix[location] = 0.0;
            location = sparse_structure->sparseLocation(kk,node * 3 + 1);
            ttdSdy[    num_dofs + kk] -= matrix[location] * xyz[0] / 2;
            ttdSdy[3 * num_dofs + kk] -= matrix[location] * xyz[1];
            ttdSdy[4 * num_dofs + kk] -= matrix[location] * xyz[2] / 2;
            if(kk == node * 3)
              ttdSdy[    num_dofs + kk] += force_vector[node * 3 + 1] / 2;
            if(kk == node * 3 + 1)
              ttdSdy[3 * num_dofs + kk] += force_vector[node * 3 + 1];
            if(kk == node * 3 + 2)
              ttdSdy[4 * num_dofs + kk] += force_vector[node * 3 + 1] / 2;
            matrix[location] = 0.0;
            location = sparse_structure->sparseLocation(kk,node * 3 + 2);
            ttdSdy[2 * num_dofs + kk] -= matrix[location] * xyz[0] / 2;
            ttdSdy[4 * num_dofs + kk] -= matrix[location] * xyz[1] / 2;
            ttdSdy[5 * num_dofs + kk] -= matrix[location] * xyz[2];
            if(kk == node * 3)
              ttdSdy[2 * num_dofs + kk] += force_vector[node * 3 + 2] / 2;
            if(kk == node * 3 + 1)
              ttdSdy[4 * num_dofs + kk] += force_vector[node * 3 + 2] / 2;
            if(kk == node * 3 + 2)
              ttdSdy[5 * num_dofs + kk] += force_vector[node * 3 + 2];
            matrix[location] = 0.0;
          }
        }
        matrix[sparse_structure->sparseLocation(node * 3    ,node * 3    )] = 1.0;
        matrix[sparse_structure->sparseLocation(node * 3 + 1,node * 3 + 1)] = 1.0;
        matrix[sparse_structure->sparseLocation(node * 3 + 2,node * 3 + 2)] = 1.0;
        con.erase(con.begin(),con.end());
      }
    }
  }
}
