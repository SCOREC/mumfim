#ifndef H_RepresentVolElem
#define H_RepresentVolElem
#include "MicroFOMultiscaleTypes.h"
#include "FiberNetwork.h"
#include "SparseMatrix.h"
#include "Sparskit_Externs.h"
#include <apf.h>
#include <apfMesh2.h>
#include <iomanip>
#include <fstream>
#include <vector>
namespace bio
{
  struct RVE_Info
  {
    int mType;
    int order;
    double * derivS;
  };
  // the structure glo_node is used for bookkeeping the position of the edges of the RVEs.
  struct glo_node
  {
    double x[8];     // x-coordinate of the edges
    double y[8];     // y-coordinate of the edges
    double z[8];     // z-coordinate of the edges
  };
  class MicroFO
  {
  public:
    MicroFO ();
    MicroFO(int * hdr,
            double * gpts,
            int rnd_id,
            FiberNetwork * fn,
            SparseMatrix * st,
            SparskitBuffers * buffers,
            double * prms,
            double * init_coords,
            int coor_size);
    ~MicroFO();
    void init();
    void SetRveCorners();
    void UpdateNodeArray();
    double calc_stiffness(); // Temporary for size effect test
    // setup the current coordinates and the disp - need to be updated incrementally
    void setDisplacement(double * displacement);
    void setDeformationGradient(double * grad) { deformation_gradient = grad; }  // set deformation gradient at gauss point i
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
    void constructRVEFromMigrationData(FiberNetwork *** networks,
                                       SparseMatrix *** matrices,
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
    int rve_tp;
    int rnd_id;
    SparskitBuffers * buffers;
    FiberNetwork * fiber_network;
    SparseMatrix * sparse_structure;
    const FiberNetwork * init_fiber_network;
    // HACKY refactor to simplify this sort of thing
    bool firstTimeThrough;
    bool post_migration;
    // RVE parameters
    double fiber_volume_fraction;
    double fiber_radius;
    // computed from RVE parameters and FiberNetwork
    double total_fiber_length;    // initial
    double fiber_area;
    double rve_dim;
    double half_rve_dim;
    double init_rve_dim;
    double scale_conversion;
    std::vector<FiberReaction*> fiber_types;
    // Fiber parameters for beam implementation
    //double fiber_G;
    //double fiber_J;
    //double fiber_I1;
    //double fiber_I2;
    // The coordinates and displacements
    double * coords;
    std::vector<double> initial_coords;
    double * displacement;
    // The deformation gradient
    double * deformation_gradient;
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
    int main_solver(double * deformation_grad,
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
    void create_element_shape_vars(double & vol,
                                   double * dvol,
                                   double * deform_grad);
    void make_dRVEdFE(double * dRVEdFE,double * lcoords);
    void getRVECornerDisp(const double F[], double rvedisp[]);
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
    //PHI FUNCTION TO COMPUTE THE BASE FUNCTION FOR DIFFERENT TOPOLOGY ELEMENT
    double PHI(int, double, double, double);
  };
}
#endif
