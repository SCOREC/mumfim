#include "NonLinFibMtx.h"

#include <PetscLAS.h>

// Simmetrix
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimParasolidKrnl.h> 
#include <SimModel.h>

// PETSc
#include <petsc.h>

// C/C++
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include <list>
#include <memory.h>
#include <fenv.h>

void display_help_string()
{
  std::cout << "Usage: MICRO_FM [OPTIONS]\n"
	    << "  [-h, --help]                              display this help text\n"
	    << "  [-l, --license path_to_license]           the path to the folder that contains\n"
	    << "                                            the Simmetrix license\n"
            << "  [-s, --model model_file]                  the model file (.smd)\n"
	    << "  [-m, --mesh mesh_file]                    the mesh file (.sms)\n";
}

std::string model_filename;
std::string mesh_filename;
std::string license_filename;

bool parse_options(int & argc, char ** & argv)
{
  bool result = true;

  while(1)
  {
    bool quit_loop = false;
    static struct option long_options[] =
      {
	/* These options don't set a flag */
	{"help",        no_argument,        0, 'h'},
	{"license",     required_argument,  0, 'l'},
	{"model",       required_argument,  0, 's'},
	{"mesh",        required_argument,  0, 'm'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int option = getopt_long(argc, argv, "h:l:m:s:", long_options, &option_index);

    switch (option)
    {
    case 0:
      // do we need this? for flags if we add some
      break;
    case 'h':
      display_help_string();
      result = false;
      break;
    case 'l':
      license_filename = optarg;
      break;
    case 's':
      model_filename = optarg;
      break;
    case 'm':
      mesh_filename = optarg;
      break;
    case '?':
      /* error, but getopt_long should alrady have printed an error message */
      display_help_string();
      result = false;
      break;
    case -1:
      // end of options
      quit_loop = true;
      break;
    default:
      break;
    }
    if (quit_loop)
      break;
  }
  return result;
}

int main(int argc, char** argv)
{
  using namespace Biotissue;
  using namespace amsi;
  using namespace Analysis;
  
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
  MPI_Init(&argc,&argv);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  parse_options(argc,argv);

  // Initialize simmetrix
  SimPartitionedMesh_start(&argc,&argv);

  //todo : change this to be passed in
  Sim_readLicenseFile(license_filename.c_str());

  MPI_Comm comm = MPI_COMM_WORLD;

  // micro_fm specific
  
  // Initialize petsc
  PETSC_COMM_WORLD = comm;
  PetscInitialize(&argc,&argv,"PetscOpt",PETSC_NULL);
  SimModel_start();
  
  // Load the model from the first filename passed in
  pProgress progress = Progress_new();
  pGModel model = GM_load(model_filename.c_str(), 0, progress);
/*
    GM_createFromNativeModel(
    ParasolidNM_createFromFile(model_file_name.c_str(), 0),
    0);
*/
  assert(model);
  
  pParMesh mesh = PM_load(mesh_filename.c_str(), sthreadNone, model, progress);
  
  NonLinFibMtx * micro_fm = new NonLinFibMtx("micro_fibermatrix",
					     model,
					     mesh,
					     0.65824e-9,3.8, // A and B parameters for nonlinear fibers
					     110, 0.10); // Young's modulus and Poission's ratio
  LAS * las = new PetscLAS(0,0);
  micro_fm->Run(las,40);
    
  GM_release(model);
  M_release(mesh);
  SimParasolid_stop(1);
  SimPartitionedMesh_stop();
  PetscFinalize();
  
  return EXIT_SUCCESS;
}
