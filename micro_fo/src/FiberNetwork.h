#ifndef BIO_FIBER_NETWORK_H_
#define BIO_FIBER_NETWORK_H_

#include "apfUtil.h"
#include "FiberReactions.h"
#include "RVE.h"
#include "SparseStructure.h"
#include "Sparskit_Externs.h"

#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfField.h>

#include <cassert>

namespace bio
{
  /**
   * Responsible for managing the internal state of a single fiber-network
   *  quasistatics simulation. 
   */
  class FiberNetwork
  {
  private:
    apf::Mesh * fn;
    apf::Field * fn_u;
    apf::Numbering * fn_dof;
    
    int dim;
  protected:
  public:

    /**
     * Construct a FiberRVE object.
     * @param f A pointer to a fiber network mesh (contains only vertices and edges)
     *          typically loaded using the NetworkLoader classes 
     */
    FiberNetwork(apf::Mesh * f);

    /**
     *  Gives the dimensionality of the managed fiber network
     *  @return the dimensionality of the fiber network (2 or 3)
     */
    int getDim() const { return fn->getDimension(); }

    /**
     * Calculate the current length of each fiber in the network
     * @param lngths A vector containing the lengths of all fibers in the network;
     * @note Possibly implement local caching of the answer instead of recalculating
     */
    void calcFiberLengths(std::vector<double> & lngths) const { calcDimMeasures(fn,1,lngths); }

    void assembleJacobian();

    const apf::Mesh * getNetworkMesh() {return fn;}
    const apf::Field * getDisplacementField() {return fn_u;}
  };

/*  
struct RVE_Info
{
  int mType;
  int order;
  double * derivS;
};

class MicroFO
{
private:
  SparskitBuffers * buffers;
  FiberReactionAssignment * fiber_reactions;

public:
  MicroFO ();
  
  MicroFO(int elemType,
	  double * gpts,
	  int networkID,
	  FiberNetwork * fn,
	  SparseMatrix * st,
	  SparskitBuffers * buffers,
	  double * init_coords,
	  int coor_size);
  
  ~MicroFO();

  void SetRveCorners();
  void UpdateNodeArray();
  double calc_stiffness(); // Temporary for size effect test

  // setup the current coordinates and the disp - need to be updated incrementally
  void SetDisplacement(double * displacement);
  void SetDeformationGradient(int gauss_pt, double * grad);  // get deformation gradient at gauss point i
  void SetDeformationGradients(double * grads);
  
  void eval_derivstress( double *dSdx);
  void getRVEs(double *ss);
  void setResults(double * ss);

  // Function to precompute the stress information for each RVE which is
  // a wrapper of the the main_solver of the original code
  int run(); 
  
  // There is no contribution for the integration order from the microscopic problem at this moment.
  int getIntegrationOrder () const { return 0; }

  //output the RVE
  void output(std::ofstream &);
  void output(const std::string & filename);
 
  void printD(char *, double *, int);
  void printI(char *, int *, int);
  void outputFiber(const std:: string & fileame) const;
  void outputFiberToVTK(const std::string & filename);
  void outputFiberToVectFile(const std::string & filename);

  
  // Migration - need to pull into AMSI - almost entirely
  void collectMigrationData();
  void clearMigrationData();
  void constructRVEFromMigrationData(FiberNetwork ** networks,
				     SparseMatrix ** matrices,
				     SparskitBuffers * b);
  std::vector<double> doubleMigrationData;
  std::vector<int> intMigrationData;
  void getMigrationData(std::vector<char> & rveData);
  void setMigrationData(std::vector<char> & rveData);

  // hard-coded weight of 1.0 for migration
  double getWeight();
  double getTiming();
  
  void resetWeight() {rve_iterations.erase(rve_iterations.begin(),rve_iterations.end());}
  void resetTiming() {rve_timing.erase(rve_timing.begin(),rve_timing.end());}
private:
  int network_id;
  FiberNetwork * fiber_network;
  SparseMatrix * sparse_structure;
  const FiberNetwork * init_fiber_network;

  // HACKY refactor to simplify this sort of thing
  bool firstTimeThrough;
  bool post_migration;

  // RVE parameters
  double fiber_volume_fraction;
  double scale_conversion;
  double HALF_RVE_DIM;

  //double fiber_radius; // unit: mm ??
  //double fiber_area;
  double fiber_E;
  double fiber_B;
  double lexp;

  // Fiber parameters for beam implementation
  double fiber_G;
  double fiber_J;
  double fiber_I1;
  double fiber_I2;

  // initial total fiber length
  double total_fiber_length; 

  // The coordinates and displacements
  double * coords;
  std::vector<double> initial_coords;
  double * displacement;

  // Macroscale element information
  //apf::Mesh2 * element_mesh;
  //apf::MeshEntity * macro_entity;
  //apf::MeshElement * macro_element;
  
  int element_type;
  int num_field_comp;
  int num_element_nodes;
  int num_rve_doubles;
  double gauss_pt[3];

  double F[3][3];  //deformation gradient for the rve
  double FItp[8][3][3];  //deformation gradient for 8 vertices of RVE

  // The output information for macroscopic problem
  double * rve_info;
  
  // rve_iterations is a vector that stores the number of Newton iterations
  std::vector<double> rve_iterations; // stores number of iterations for each run, after convergence
  std::vector<double> rve_timing; // stores time to solve RVE problem.
  // The following are the variables related to microscopic problem
  int num_gauss_pts;        // number of Gauss Points employed

  std::vector<double> tdydxr;
  std::vector<double> ttdSdy;

  std::vector<double> matrix;
  std::vector<double> matrix_axial; // for beam solver.
  std::vector<double> force_vector;
  std::vector<double> force_vector_axial; //for beam solver.
  std::vector<double> coordinate_vector;

  apf::Vector3 rve[8];
  apf::Vector3 u[8];

  //The following are the functions that must access the variables
  
  int main_solver(double * coords_loc,
		   double * fedisp,
		   double & local_S11, double & local_S12, double & local_S13,
		                       double & local_S22, double & local_S23,
                                                           double & local_S33,
		   double * loc_dSdx,
		   double & loc_vastrx,
		   double & loc_vastry,
		   double & loc_vastrz,
		   double & fem_res_norm);

  // map node dof values to solution vector
  void update_coordinate_vector();

  // map solution vector values to node dof values
  void update_nodes();

  // read the fiber network in from the specified file, only done once
  //void read_init_net(const std::string & filename);

  void create_element_shape_vars(double & vol,
				 double * dvol, 
				 double * fedisp);

  void make_dRVEdFE(double * dRVEdFE,double * lcoords);
  
  int Solver();
  int Solver_Beam();

  void calc_precond(double * matrix,
		    double * lengths, 
		    double * fib_str, 
		    double * dfdE,
		    bool calc_dSdy);
  
  double gfunc(int element,
	       int flag,
	       double * lengths,
	       double * fib_str,
	       double * dfdE);
  
  void calc_quantities(int side,
		       double * dridx,
		       int node,
		       int node1,
		       int node2,
		       int node3,
		       int node4);
  
  void calc_dridx(double * dridx);

  void calc_stress(double * stress);

  void avg_vol_stress(double * stress,
		      double & loc_vastrx, 
		      double & loc_vastry,
		      double & loc_vastrz,
		      double vol,
		      double fem_res_norm);

  void calc_tdydxr();

  void calc_femjacob_newmethod(double * dSdx,
                               double vol,
                               double * dvol,
                               double * stress);

  void fvec_bcs(double fvec[],int pas);
  void computeRVEboundary(double rvedisp[]);

  void computeInterpolatedRVEboundary(double rvedisp[]);
  
  double compute_total_fiber_length();
  double calc_fiber_len_avg();
  void calc_fiber_str(double * lengths, double * fib_str, double * dfdE);
  void calc_force_vector(double * lengths, double * fib_str);
  void force_vector_bcs();
  void calc_mean_fiber_stretch_omega();

  // Beam specific functions
  void calc_current_disps_beam(double * solution);
  void calc_matrix_beam(double * matrix,double * lengths);
  void calc_matrix_beam(double * matrix, double * matrix_axial, double * lengths);
  void calc_force_vector_beam(double * x);
  void matrix_bcs_beam();
  void matrix_periodic_bcs_helper(const std::vector<PBCRelation> & bcs, int ddof, int * pdofs);
  void make_material_stiffness_matrix(double * matrix,int ielem,double L);
  void make_geometry_stiffness_matrix(double * matrix,int ielem,double L);
  double calc_norm_beam(double * x);
  void force_vector_bcs_beam();
  void force_vector_periodic_bcs_helper(const std::vector<PBCRelation> & bcs, int ddof, int * pdofs);
  void calc_dsdy_beam();

  // constants carried over from previous implementation
  // moved to functions so we can have fiber values differ between individual fibers
  double fiber_get_E(int fiber) {return fiber_E;}
  double fiber_get_B(int fiber) {return fiber_B;}
  //double fiber_get_A(int fiber) {return fiber_area;}

  double fiber_get_G(int fiber) {return fiber_G;}
  double fiber_get_J(int fiber) {return fiber_J;}
  double fiber_get_I1(int fiber) {return fiber_I1;}
  double fiber_get_I2(int fiber) {return fiber_I2;}

  //PHI FUNCTION TO COMPUTE THE BASE FUNCTION FOR DIFFERENT TOPOLOGY ELEMENT
  double PHI(int, double, double, double);

};
*/
 }

#endif
