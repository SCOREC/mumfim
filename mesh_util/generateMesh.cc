#include "MeshSim.h"
#include "SimParasolidKrnl.h"
#include "SimMessages.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include "SimMeshingErrorCodes.h"
#include <stdio.h>
#include <stdlib.h>

// Usage : generateMesh [model_in] [model_out] [mesh_out] [mesh_size]

void messageHandler(int type, const char *msg);
void progressHandler(const char *what, int level, int startVal,
                     int endVal, int currentVal, void *);

double mesh_size = 0.1;

int main(int argc, char *argv[])
{

  pParasolidNativeModel nModel;
  pGModel model;
  pMesh mesh;

  mesh_size = atof(argv[4]);

  // You will want to place a try/catch around all Simmetrix calls
  // as errors are thrown.
  try {
    // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
    // pass in the location of a file containing your keys.  For a release
    // product, use Sim_registerKey()
    Sim_readLicenseFile(0);

    Sim_logOn("generate_mesh.log");
    MS_init();
    SimParasolid_start(1);

    Sim_setMessageHandler(messageHandler);
    pProgress progress = Progress_new();
    Progress_setCallback(progress, progressHandler);

    nModel = ParasolidNM_createFromFile(argv[1],0);
    model = GM_createFromNativeModel(nModel,progress);
    NM_release(nModel);

    // check that the input model is topologically valid before meshing
    pPList modelErrors = PList_new();
    if(!GM_isValid(model,0,modelErrors)) {
      printf("Input model not valid\n");
      printf("Number of errors returned: %d\n",PList_size(modelErrors));
      GM_release(model);
      return 1;
    }
    PList_delete(modelErrors);

    pModelItem modelDomain = GM_domain(model);
    mesh = M_new(0,model);
    pACase meshCase = MS_newMeshCase(model); // create a meshing case
    MS_setMeshSize(meshCase,modelDomain,2,mesh_size,0); // set a relative mesh size on the model

    MS_setSurfaceMeshType(meshCase,modelDomain,3); // only generate quads via surface meshing

    // create the meshers and run
    pSurfaceMesher surfaceMesher = SurfaceMesher_new(meshCase, mesh);
    SurfaceMesher_execute(surfaceMesher,progress);
    pVolumeMesher volumeMesher = VolumeMesher_new(meshCase, mesh);
    VolumeMesher_execute(volumeMesher,progress);

    GM_write(model,argv[2],0,progress); // write out the model before the mesh!
    M_write(mesh,argv[3], 0,progress); // writes out the mesh to Simmetrix format

    // clean up meshers, meshing case, and mesh
    SurfaceMesher_delete(surfaceMesher);
    VolumeMesher_delete(volumeMesher);
    MS_deleteMeshCase(meshCase);
    M_release(mesh);

    // cleanup
    GM_release(model);
    Progress_delete(progress);
    SimParasolid_stop(1);
    MS_exit();
    Sim_logOff();
    Sim_unregisterAllKeys();

  } catch (pSimError err) {
    printf("Simmetrix error caught:\n");
    printf("  Error code: %d\n",SimError_code(err));
    printf("  Error string: %s\n",SimError_toString(err));
    SimError_delete(err);
    return 1;
  } catch (...) {
    printf("Unhandled exception caught\n");
    return 1;
  }
  return 0;
}

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    printf("Info: %s\n",msg);
    break;
  case Sim_DebugMsg:
    printf("Debug: %s\n",msg);
    break;
  case Sim_WarningMsg:
    printf("Warning: %s\n",msg);
    break;
  case Sim_ErrorMsg:
    printf("Error: %s\n",msg);
    break;
  }
  return;
}

void progressHandler(const char *what, int level, int startVal,
                     int endVal, int currentVal, void *)
{
  printf("Progress: %s, level: %d, startVal: %d, endVal: %d, currentVal: %d\n",
         what,level,startVal,endVal,currentVal);
  return;
}

