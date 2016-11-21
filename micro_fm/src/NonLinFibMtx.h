// todo: maintain a list of embedded fibers rather than accumulating them in every damn assembly function
#ifndef H_NONLINFIBMTX
#define H_NONLINFIBMTX
#include "NonFiniteElement.h"
// Analysis
//#include <SimFEA.h>
#include <apfSimFEA.h>
// Fields
#include <apf.h>
// Mesh/Geom Library
#include <MeshSim.h>
// Standard Libs
#include <memory.h>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <list>
#include <map>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <cstring>
// not exported from simmetrix libs // pull up into femanalysis with conditional compilation
extern const apf::Matrix3x3 Identity;
namespace Biotissue
{
  using namespace amsi;
  class NonLinFibMtx : public apfSimFEA
  {
  protected:
    // Nonlinear Multi-step FEM Solver State Variables
    int nonlinear_iteration;
    int load_step;
    int integration_order;
    // FEM Tensor fields and dof information
    int field_components;
    apf::FieldShape * disp_field_shape;
    int num_local_unknowns;
    int global_unknown_offset;
    int num_global_unknowns;
    apf::Field * accumulated_disp;
    // Parallel variables
    int local_rank;
    int local_size;
    // Material parameters
    double A;
    double B;
    double shear_modulus;
    double poisson_ratio;
    double standard_volume;
    // Derived mechanical values
    double total_embedded_length;
    double average_stretch_ratio;
    // Micro-to Macro multiscale linking terms
    apf::Matrix3x3 total_stress;
    apf::Matrix3x3 matrix_stress;
    apf::Matrix3x3 fiber_stress;
    apf::Vector3 unbalanced_force;
    apf::Vector3 matrix_term;
    apf::Vector3 fiber_term;
    // Macro-to-micro multiscale linking terms
    apf::Matrix3x3 deformation_gradient;   // deformation gradient
    apf::Matrix3x3 displacement_gradient;  // spatial displacement gradient
    // convenience lists to retain commonly-queried data
    std::list<pGEntity> embedded_fibers;
    std::vector<pGEntity> boundary_geom_faces;
    //std::vector<apf::Vector3> boundary_face_normals;
  public:
    NonLinFibMtx(pGModel in_model,
                 pParMesh in_mesh,
                 pACase ipd,
                 double A, double B,
                 double E, double V,
                 MPI_Comm cm = AMSI_COMM_SCALE);
    ~NonLinFibMtx();
    // initialize member variables that require more complicated initializations then can be handled in an initialized list
    void Initialize();
    void DirectionAndStretch(apf::Vector3 & net_direction,
                             double & mean_stretch);
    void BoundaryFiberEdges(std::list<pEdge> & boundary_edges);
    void GetAverageFiberStress(LAS * las, apf::Matrix3x3 & fiber_stress);
    void GetAverageMatrixStress(apf::Matrix3x3 & matrix_stress);
    void Output_TopSurfaceDisplacement(const std::string & filename1,
                                       const std::string & filename2);
    void Run(LAS * las, int max_load_step);
    void Output_FiberNetwork(const std::string & filename);
    void Output_FiberMatrix(const std::string & filename);
    void Output_MatrixElementStrain(const std::string & filename);
    void Output_MatrixElementSmallStrain(const std::string & filename);
    void Output_MatrixElementPrincipleStrain(const std::string & filename);
    void Output_FiberForce(const std::string & filename);
    void Output_FiberCrosslinkDisplacement(const std::string & filename);
    void Output_FiberMatrixPrincipleStress(const std::string & filename);
    void Output_FiberStretch(const std::string & filename);
    void Output_FiberInTension(const std::string & filename);
    // interface functions and also functions used to calculate scale linking info
    void GetAverageStress(apf::Matrix3x3 & stress) {stress = total_stress;}
    void GetUnbalancedForce(apf::Vector3 & Q) {Q = unbalanced_force;}
    void AverageStress(apf::Matrix3x3 & stress);
    void UnbalancedForce(apf::Vector3 & term);
    void AverageMatrixStress(apf::Matrix3x3 & stress);
    void AverageFiberStress(apf::Matrix3x3 & stress);
    void UnbalancedMatrixForce(apf::Vector3 & term);
    void UnbalancedFiberForce(apf::Vector3 & term);
    // todo: shift to apf::MeshEntity and add asserts in simmetrix-specific stuff
    void Mesh_BoundaryFace_GetRegion(pFace face,
                                     pRegion & region);
    void LocalSpatialDispGrad(pRegion region,
                              const apf::Vector3 & pos,
                              apf::Matrix3x3 &disp_grad);
    void FaceNormByGrad(pGFace gface,
                        apf::Vector3 & normal);
    void CauchyStress(pEntity entity,
                      apf::Matrix3x3 & stress);
    void PK2ToCauchy(const apf::Matrix3x3 & def_grad,
                     const apf::Matrix3x3 & pk2_stress,
                     apf::Matrix3x3 & cauchy_stress);
    void SpatialDispGrad(apf::Matrix3x3 & grad);
  protected:
    void InitialGuess();
    void AccumulateDisplacement();
    void SetBoundaryGeomFaces();
    void FiberForceVector(LAS * las);
    void MoveMesh();
    void Output_WriteMesh();
    void Config_CheckRegionValidity();
    void ComputeVector(LAS * las, double & residual_norm);
    void ComputeMatrix(LAS * las);
    //void ComputeForceResidual(double & resid);
    virtual void Assemble(LAS * las) {};
    //void get3DMtrlTensorC(double C[6][6]);
    void MatrixStiffnessMatrix(LAS * las);
    void MatrixForceVector(LAS * las);
    void FiberStiffnessMatrix(LAS * las);
    void MatrixElementStiffnessMatrix(apf::MeshEntity * me,
                                      apf::DynamicMatrix & Ke);
    void MatrixElementForceVector(apf::MeshEntity * me,
                                  apf::DynamicVector & fe);
/*
  void FiberAxialForce(double StchRat,
  double initLen,
  double & DfrcDlen,
  double & FibFrc);
*/
    void AssembleFiberNetwork(LAS * las,
                              pGEntity gentity);
    void TrussElementForce(apf::MeshEntity * edge,
                           apf::DynamicVector & f);
    void TrussElementStiffness(apf::MeshEntity * edge,
                               apf::DynamicMatrix & Ke);
/*
  void DeformationGradient(apf::MeshEntity * me,
  const apf::Vector3 & p,
  apf::Matrix3x3 & grad,
  double & grad_det);
  void LeftCauchy(const apf::Matrix3x3 & F,
  apf::Matrix3x3 & B);
  void RightCauchy(const apf::Matrix3x3 & F,
  apf::Matrix3x3 & C);
  void LinearizedConstitutive(const apf::Matrix3x3 & C,
  double J,
  apf::Matrix<6,6> & hookes);
  void PK2Stress(const apf::Matrix3x3 & C,
  double J,
  apf::Matrix3x3 & PK2_stress);
*/
    void SurfItg(pFace face, apf::Vector3 & values);
    void ElementStress(apf::MeshEntity * me,
                       apf::Matrix3x3 & pk2);
    // todo: pull out into seperate file (no class)
    // Compute the prinicple stresses
    void PrincipleStress(const apf::Matrix3x3 & stress,
                         apf::Vector3 principle_stress);
    void FiberForce(pEntity entity, double & fiber_force);
    // todo: pull out into seperate file (no class)
    // use the deformation gradient to compute the small strain
    void SmallStrain(const apf::Matrix3x3 & F,
                     apf::Matrix3x3 & strain);
    // compute the elemental small strain tensor
    void ElementSmallStrain(pEntity entity,
                            apf::Matrix3x3 & strain);
  };
}
#endif
