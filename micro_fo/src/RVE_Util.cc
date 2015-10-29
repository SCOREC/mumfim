/******************************************************************************

  (c) 2014 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#include "RVE_Util.h"
#include "Util.h"

#include "globals.h"

// Fiber-only 
#include "RepresentVolElem.h"
#include "MicroFOMultiscaleTypes.h"
#include "SparseMatrix.h"
#include "Sparskit_Externs.h"

// for integration point information (parametric coords)
#include <apf.h>
#include <apfMesh.h>

// Multiscale control stuff (AMSI)
#include <ControlService.h>
#include <Migration.h>
#include <amsiReporter.h>

// Simmetrix libs
#include <MeshSim.h>

// Standard libs
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <list>
#include <sstream>
#include <string>

  namespace Biotissue {

    int rve_load_balancing = -2;
    bool lb_per_iteration = true;

    int P_computeRVEs(const std::string & fiber_network_filename, int num_fiber_files)
    {
      int result = 0;
      srand(8675309);
      
      // AMSI communication variables
      amsi::ControlService * cs = amsi::ControlService::Instance();
      amsi::Task * local_task = amsi::getLocal();
      amsi::DataDistribution * dd_rve = local_task->getDD("macro_fo_data");
      int local_rank = local_task->localRank();

      amsi::TaskManager * tm = cs->GetTaskManager();
      int macro_size = taskSize(tm->getTask("macro"));
      int micro_size = taskSize(tm->getTask("micro_fo"));

#     ifdef LOGRUN
      amsi::Log * micro_fo_efficiency = amsi::makeLog("micro_fo_efficiency");
      amsi::Log * micro_fo_weights = amsi::makeLog("micro_fo_weights");
      amsi::Log * micro_fo_timing = amsi::makeLog("micro_fo_timing");
      
      amsi::log(micro_fo_efficiency) << "MACRO_STEP, MACRO_ITER, TIME, STATUS, DESCRIPTION" << std::endl;
      amsi::log(micro_fo_efficiency) << "0, 0, 0.0, ACTIVE, START" << std::endl;

      amsi::log(micro_fo_weights) << "MACRO_STEP, WEIGHT" << std::endl;
      amsi::log(micro_fo_timing) << "MACRO_STEP, TIME, NUM_MACRO_ITER" << std::endl;
#     endif
      
      // Read in all the fiber networks
      std::vector<FiberNetwork*> networks;
      std::vector<SparseMatrix*> sparse_structs;
      networks.reserve(num_fiber_files);
      std::stringstream current_file;
      int dof_max = -1;
      for(int ii = 0; ii < num_fiber_files; ii++)
      {
	current_file << fiber_network_filename << ii+1 << ".txt"; // Files count from 1

	FiberNetwork * fn = NULL;

	if (SPECIFY_FIBER_TYPE) //SPECIFY_FIBER_TYPE is true
	  fn = new SupportFiberNetwork();
	else
	  fn = new FiberNetwork();

	fn->readFromFile(current_file.str());
	networks.push_back(fn);
	current_file.str("");

	sparse_structs.push_back(Make_Structure(networks.back()));

	int dofs = networks.back()->numDofs();
	dof_max = dofs > dof_max ? dofs : dof_max;
      }

      SparskitBuffers * buffers = new SparskitBuffers(dof_max);

      size_t r_id = cs->CommRelation_GetID("macro","micro_fo");;
      size_t recv_pattern_id = -1;
      size_t recv_pattern_init_id = -1;
      size_t send_pattern_id = -1;
      
      InitializeScaleCoupling(recv_pattern_id,send_pattern_id);
      
      // Declare and commit (to mpi) multiscale data types
      MicroFOMultiscaleDataTypes mdt;
      mdt.MultiscaleDataTypesMPICommit();
      MPI_Datatype micro_fo_init_data_type = mdt.micro_fo_init_data_type;    
      MPI_Datatype micro_fo_data_type = mdt.micro_fo_data_type;
      MPI_Datatype micro_fo_result_type = mdt.micro_fo_result_type;
      
      // load step and iteration counters
      int step_iter[2];
      step_iter[0] = step_iter[1] = 0; // Load step and iteration step, respectively
      int timestep_converged = false;
      int iteration_computed = false;

      // List of current RVEs
      std::vector<MicroFO*> rves;

      // Multiscale data
      std::vector<micro_fo_init_data> initial_data;
      std::vector<micro_fo_data> data;
      std::vector<micro_fo_result> results;

      double local_weight = 0.0;
      double local_time = 0.0;
      
      while(!timestep_converged)
      {
        iteration_computed = false;
        step_iter[1] = 0;
	local_weight = 0.0;
	local_time = 0.0;
	
	// Get RVEs to be removed
	std::vector<int> to_delete;
	cs->RemoveData(recv_pattern_id,to_delete);
	std::sort(to_delete.begin(),to_delete.end());
	
	// iterate over indices to be removed in reverse order so that latter indices are not altered by the removal of earlier ones
	for(std::vector<int>::reverse_iterator iter = to_delete.rbegin(), iterend = to_delete.rend(); iter != iterend; iter++)
	  erase_ptr_at(rves,*iter);
	
	// Get RVEs to be added
	std::vector<int> to_add;
	std::vector<int> dummy;
	recv_pattern_init_id = cs->AddData(recv_pattern_id,dummy,to_add);
	
	// Add RVEs to list
	int ii = 0;
	int jj = 0;
	std::vector<MicroFO*>::iterator currentRVE = rves.begin();
	// the *old* indices of these shouldn't matter... just shove them on the end...
	rves.resize(rves.size()+to_add.size());
	
	// Get initial data for those RVEs
	cs->Communicate(recv_pattern_init_id,
			initial_data,
			micro_fo_init_data_type);
	  
	// Update new RVEs with initial data
	ii = 0;
	for(currentRVE = rves.begin(); currentRVE != rves.end(); ++currentRVE)
	{
	  if(*currentRVE == NULL)
	  {
	    // Construct new RVE
	    double pt[3];
	    micro_fo_init_data * cid = &initial_data[ii];
	    
	    apf::Vector3 gauss_pt;
	    apf::getGaussPoint(cid->element_type,1,cid->gauss_pt_id,gauss_pt);
	    gauss_pt.toArray(pt);
	    
	    int num_nodes = NumElementNodes(cid->element_type);

	    // select one of the preloaded fiber-networks to associate with this RVE
	    int rand_network = rand() % networks.size();
	    (*currentRVE) = new MicroFO(cid->element_type,
					pt,
					rand_network,
					networks[rand_network],
					sparse_structs[rand_network],
					buffers,
					&cid->init_data[0],
					num_nodes);
	    ii++;
	  }
	}

	//  Update macro->micro communication
	
	// Update micro->macro communication
	/*
	cs->CommPattern_UpdateInverted(recv_pattern_id,send_pattern_id);
	cs->CommPattern_Assemble(send_pattern_id);
	cs->CommPattern_Reconcile(send_pattern_id);
	*/
	
	// ***************************************
	// Migration of RVEs
	// ***************************************
	/*
	if(step_iter[0] > 0)
	{
#         ifdef LOGRUN
	  double pre_migration = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 << pre_migration << ", IDLE, PRE_MIGRATION" << std::endl;
#         endif
	  
	  std::vector< int > migration_indeces;
	  std::vector< std::vector<char> > rve_data;
	  
//
	  amsi::CommPattern * scale_coupling = cs->getCouplingManager()->CommPattern_Get(recv_pattern_id);
	  amsi::Migration * mg = new amsi::ScaleSensitiveMigration(scale_coupling,
								   AMSI_COMM_LOCAL,
								   amsi::Migration::ZOLTAN_ALGO);

	  std::vector<double> weights;
	  dd_rve->getWeights(local_rank,weights);
	  mg->plan(migration_indeces,micro_size,&weights[0]);
//


	  cs->planMigration(migration_indeces,recv_pattern_id,rve_load_balancing);
	  
	  // Collect data of RVEs to migrate elsewhere
	  //   and delete and remove those RVEs
	  if(migration_indeces.size() > 0)
	  {
	    rve_data.resize( migration_indeces.size() );
	    std::sort(migration_indeces.begin(),migration_indeces.end(),std::greater<int>());
	    int jj = 0;
	    for(std::vector<int>::iterator index = migration_indeces.begin();
		index != migration_indeces.end(); index++)
	    {
	      rves[*index]->collectMigrationData();
	      rves[*index]->getMigrationData(rve_data[jj]);
	      jj++;
	      erase_ptr_at(rves,*index); // removes the rve from the vector, and deallocates the MicroFO object
	    }
	    std::reverse(rve_data.begin(),rve_data.end()); // THIS IS HUGELY IMPORTANT TO MAINTAIN THE CORRECT ORDERING SINCE THESE ARE ASSUMED TO BE SENT IN ASCENDING ORDER ... i think...
	  }
	  
	  // AMSI call to do migration communication... needs to be called regardless of whether there are local RVES to be migrated
	  cs->migration(rve_data,recv_pattern_id,rves);

	  //mg->execute(rve_data);
	  
	  // After migration use data to reconstruct RVEs that have migrated here
	  //   NULL values have been added to RVE list, just need to fill
#         ifdef LOGRUN
	  double pre_construct_new = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 << pre_construct_new << ", ACTIVE, PRE_CONSTRUCT_NEW" << std::endl;
#         endif	  
	  if(rve_data.size() > 0)
	  {
	    int jj = 0;
	    for(std::vector<MicroFO*>::iterator rve = rves.begin(); rve != rves.end(); ++rve)
	    {
	      if((*rve) == NULL)
	      {
		(*rve) = new MicroFO();
		(*rve)->setMigrationData(rve_data[jj]);
		(*rve)->constructRVEFromMigrationData(&networks[0],&sparse_structs[0],buffers);
		jj++;
	      }
	    }
	  }
#         ifdef LOGRUN
	  double post_construct_new = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 << post_construct_new << ", ACTIVE, POST_CONSTRUCT_NEW" << std::endl;

#         endif 
	}
	
	// Share new comm pattern with Macroscale (macroscale will need to rearrange element list)
	cs->shareMigration(recv_pattern_id,dummy);
*/	
	// Update micro->macro communication
	cs->CommPattern_UpdateInverted(recv_pattern_id,send_pattern_id);
	cs->CommPattern_Assemble(send_pattern_id);
	cs->CommPattern_Reconcile(send_pattern_id);
/*
#       ifdef LOGRUN
	double post_migration = amsi::getElapsedTime(micro_fo_efficiency);
	amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
				       << post_migration << ", ACTIVE, POST_MIGRATION" << std::endl;
p*/
#       ifdef LOGRUN
	double pre_iterations = amsi::getElapsedTime(micro_fo_efficiency);
	amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 << pre_iterations << ", ACTIVE, PRE_ITERATIONS" << std::endl;
#       endif
	
        while(!iteration_computed)
        {
	  if(lb_per_iteration && step_iter[1] > 0) // start load balancing after first iteration
	  {
#           ifdef LOGRUN
	    double pre_migration = amsi::getElapsedTime(micro_fo_efficiency);
	    amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					   << pre_migration << ", ACTIVE, PRE_MIGRATION" << std::endl;
#           endif
	    
	    std::vector<int> migration_indices;
	    std::vector< std::vector<char> > rve_data;

	    cs->planMigration(migration_indices,recv_pattern_id,rve_load_balancing);

	    if(migration_indices.size() > 0)
	    {
	      rve_data.resize( migration_indices.size() );
	      std::sort(migration_indices.begin(),migration_indices.end(),std::greater<int>());
	      int jj = 0;
	      for(std::vector<int>::iterator index = migration_indices.begin();
		  index != migration_indices.end(); index++)
	      {
		rves[*index]->collectMigrationData();
		rves[*index]->getMigrationData(rve_data[jj]);
		jj++;
		erase_ptr_at(rves,*index); // removes the rve from the vector, and deallocates the MicroFO object
	      }
	      std::reverse(rve_data.begin(),rve_data.end()); // THIS IS HUGELY IMPORTANT TO MAINTAIN THE CORRECT ORDERING SINCE THESE ARE ASSUMED TO BE SENT IN ASCENDING ORDER ... i think...
	    }
	    
	    cs->migration(rve_data,recv_pattern_id,rves);

	    if(rve_data.size() > 0)
	    {
	      int jj = 0;
	      for(std::vector<MicroFO*>::iterator rve = rves.begin(); rve != rves.end(); ++rve)
	      {
		if((*rve) == NULL)
		{
		  (*rve) = new MicroFO(/*rve_data[jj]*/);
		  (*rve)->setMigrationData(rve_data[jj]);
		  (*rve)->constructRVEFromMigrationData(&networks[0],&sparse_structs[0],buffers);
		  jj++;
		}
	      }
	    }
	    
	    std::vector<int> dummy;
	    cs->shareMigration(recv_pattern_id,dummy);

	    cs->CommPattern_UpdateInverted(recv_pattern_id,send_pattern_id);
	    cs->CommPattern_Assemble(send_pattern_id);
	    cs->CommPattern_Reconcile(send_pattern_id);
#           ifdef LOGRUN
	    double post_migration = amsi::getElapsedTime(micro_fo_efficiency);
	    amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					   << post_migration << ", ACTIVE, POST_MIGRATION" << std::endl;
#           endif
	  }

#         ifdef LOGRUN
          double pre_recv = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 << pre_recv << ", IDLE, PRE_RECV" << std::endl;
#         endif

          // Get macroscale data (current coordinates and incremental displacements)
	  cs->Communicate(recv_pattern_id,data,micro_fo_data_type);

#         ifdef LOGRUN
	  double post_recv = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					   << post_recv << ", ACTIVE, POST_RECV" << std::endl;
#         endif
	  results.resize(data.size()); // shouldn't change except between load steps
#         ifdef LOGRUN
	  double pre_compute = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 << pre_compute << ", ACTIVE, PRE_SOLVE" << std::endl; 
#         endif
          int ii = 0;
          for(std::vector<MicroFO*>::iterator rve = rves.begin(); rve != rves.end(); ++rve)
          {
	    result += P_computeFiberOnlyRVE((*rve),&data[ii].data[0],&results[ii].data[0]);
	    double w = (*rve)->getWeight();
	    //std::cout << w << std::endl;
	    local_weight += w;
	    dd_rve->setWeight(local_rank,ii,w);
	    (*rve)->resetWeight();
            ii++;
	  }
#         ifdef LOGRUN
          double post_compute = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					   << post_compute << ", ACTIVE, POST_SOLVE" << std::endl;
	  
	  double pre_send = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					   << pre_send << ", IDLE, PRE_SEND" << std::endl;
#         endif
	  cs->Communicate(send_pattern_id,results,micro_fo_result_type);
#         ifdef LOGRUN
	  double post_send = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 << post_send << ", ACTIVE, POST_SEND" << std::endl;
#         endif

          // Erase data to prepare for next step
	  data.erase(data.begin(),data.end());
	  results.erase(results.begin(),results.end());

#         ifdef LOGRUN
	  double pre_continue_iters = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					   << pre_continue_iters << ", IDLE, PRE_CONTINUE_ITERS" << std::endl;
#         endif
	  
          cs->Relation_Broadcast(r_id,iteration_computed,MPI_INTEGER);
	  
#         ifdef LOGRUN
	  double post_continue_iters = amsi::getElapsedTime(micro_fo_efficiency);
	  amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 << post_continue_iters<< ", ACTIVE, POST_CONTINUE_ITERS" << std::endl;
#         endif

#         if defined(SCOREC) && defined(LOGRUN)
	  AMSI_DEBUG (
	  if (step_iter[0] == 2 || step_iter[0] == 35)
	  {
	    std::stringstream output_file;
	    output_file<<"MicroFiberNetwork_rank"<<local_rank<<"_"<<"ldstp_"<<step_iter[0]+1<<"_"<<step_iter[1]+1<<".txt"; //+1 to step_iter for post processing purposes, Matlab reads from index 1 and not 0.
	    rves[0]->outputFiber(output_file.str());
	  }
	    )
#         endif
	  step_iter[1]++;
        } // load_step computed
	step_iter[1]--;

#       ifdef LOGRUN
	double post_iterations = amsi::getElapsedTime(micro_fo_efficiency);
	amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
                                         << post_iterations << ", ACTIVE, POST_ITERATION" << std::endl; 
#       endif
	
	ii = 0;
	for(std::vector<MicroFO*>::iterator rve = rves.begin(); rve != rves.end(); ++rve)
	{
	  double w = (*rve)->getWeight();
	  local_weight += w;
	  dd_rve->setWeight(local_rank,ii++,w);
	  (*rve)->resetWeight();

#         ifdef LOGRUN
	  double time = (*rve)->getTiming();
	  //std::cout << "RVE took " << time << " seconds to converge" << std::endl;
	  local_time += time;
	  (*rve)->resetTiming();
#         endif
	}

#       if defined(LOGRUN) && defined(SCOREC)
	std::cout << "[" << local_rank << "] Total local weight: " << local_weight << std::endl;
	std::cout << "[" << local_rank << "] Total local time (sec): " << local_time << std::endl; 
	std::cout << "[" << local_rank << "] number of macro iterations: " << step_iter[1] << std::endl;
	
	AMSI_DEBUG (
	  std::stringstream output_file;	
	  output_file<<"FiberNetwork_rank"<<local_rank<<"_"<<step_iter[0]+1<<".txt"; //step_iter + 1 for post processing, Matlab reads from index 1 and not 0.
	  rves[0]->outputFiber(output_file.str());
	  )
#       endif
	
#       ifdef LOGRUN
        double pre_continue_steps = amsi::getElapsedTime(micro_fo_efficiency);
	amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
				       << pre_continue_steps << ", IDLE, PRE_CONTINUE_STEPS" << std::endl;
#       endif

        cs->Relation_Broadcast(r_id,timestep_converged,MPI_INTEGER);

#       ifdef LOGRUN
        double post_continue_steps = amsi::getElapsedTime(micro_fo_efficiency);
	amsi::log(micro_fo_efficiency) << step_iter[0] << ", " << step_iter[1] << ", "
					 <<post_continue_steps << ", ACTIVE, POST_CONTINUE_STEPS" << std::endl;
	amsi::log(micro_fo_weights) << step_iter[0] << ", " << local_weight << std::endl;
	amsi::log(micro_fo_timing) << step_iter[0] << ", " << local_time << ", " << step_iter[1] << std::endl;
#       endif
	
	step_iter[0]++;
	step_iter[1]=0;

#       ifdef LOGRUN
	std::stringstream filename;
	filename << "micro_fo_efficiency." << local_rank << ".log";
	amsi::flush2FStream(micro_fo_efficiency, filename.str());
	
	filename.str("");
	filename << "micro_fo_weights." << local_rank << ".log";
	amsi::flush2FStream(micro_fo_weights, filename.str());

	filename.str("");
	filename << "micro_fo_timing." << local_rank << ".log";
	amsi::flush2FStream(micro_fo_timing, filename.str());
#       endif	
	
      } // done

#     ifdef LOGRUN
      amsi::deleteLog(micro_fo_efficiency);
      amsi::deleteLog(micro_fo_weights);
      amsi::deleteLog(micro_fo_timing);
#     endif
      
      return result;
    }
    
    int P_computeFiberOnlyRVE(MicroFO * rve,
			      double * disp_coords,
			      double * results)
    {
      int result = 0;
      // Setup the displaced coords
      rve->SetDisplacement(disp_coords);
      
      // Set results array (sets pointer to rve info to results, 
      // so fills directly instead of copying, this replaces getRVEs)
      rve->setResults(results);
      
      // Run the fiber-only simulation
      result += rve->run();
      return result;
    }

    void InitializeScaleCoupling(size_t & recv_pattern_id, size_t & send_pattern_id)
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      size_t r_id = cs->CommRelation_GetID("macro","micro_fo");
      recv_pattern_id = cs->RecvCommPattern("micro_fo_data",
					    "macro",
					    "macro_fo_data",
					    "micro_fo");
      cs->CommPattern_Reconcile(recv_pattern_id);

      send_pattern_id = cs->CommPattern_UseInverted(recv_pattern_id,
						    "macro_fo_data",
						    "micro_fo",
						    "macro");
      cs->CommPattern_Assemble(send_pattern_id);
      cs->CommPattern_Reconcile(send_pattern_id);
    }

    int NumElementNodes(int element_type)
    {
      int result = -1;
      switch(element_type) {
      case apf::Mesh::HEX: // Hex
	result = 8;
	break;
      case apf::Mesh::PRISM: // Prism
	result = 6;
	break;
      case apf::Mesh::TET: // Tet
	result = 4;
	break;
      case apf::Mesh::PYRAMID: // Pyramid
	result = 5;
	break;
      default:
	std::cerr << "Unsupported element topology!" << std::endl;
	break;
      }
      return result;
    }


    // basis functions
    void tsfun(double * phi, 
	       double * phic, 
	       double * phie, 
	       double * phis, 
	       double x, 
	       double y, 
	       double z)
    {
      // Shape functions for bilinear quadrilateral.
      phi[0] = (1. / 8.) * (1 - x) * (1 - y) * (1 + z);
      phi[1] = (1. / 8.) * (1 + x) * (1 - y) * (1 + z);
      phi[2] = (1. / 8.) * (1 - x) * (1 - y) * (1 - z);
      phi[3] = (1. / 8.) * (1 + x) * (1 - y) * (1 - z);
      phi[4] = (1. / 8.) * (1 - x) * (1 + y) * (1 + z);
      phi[5] = (1. / 8.) * (1 + x) * (1 + y) * (1 + z);
      phi[6] = (1. / 8.) * (1 - x) * (1 + y) * (1 - z);
      phi[7] = (1. / 8.) * (1 + x) * (1 + y) * (1 - z);
// Derivative of phi w.r.t x (N_A,x) 
      phic[0] = -(1. / 8.) * (1 - y) * (1 + z);
      phic[1] = (1. / 8.) * (1 - y) * (1 + z);
      phic[2] = -(1. / 8.) * (1 - y) * (1 - z);
      phic[3] = (1. / 8.) * (1 - y) * (1 - z);
      phic[4] = -(1. / 8.) * (1 + y) * (1 + z);
      phic[5] = (1. / 8.) * (1 + y) * (1 + z);
      phic[6] = -(1. / 8.) * (1 + y) * (1 - z);
      phic[7] = (1. / 8.) * (1 + y) * (1 - z);
// Derivative of phi w.r.t y (N_A,y)
      phie[0] = -(1. / 8.) * (1 - x) * (1 + z);
      phie[1] = -(1. / 8.) * (1 + x) * (1 + z);
      phie[2] = -(1. / 8.) * (1 - x) * (1 - z);
      phie[3] = -(1. / 8.) * (1 + x) * (1 - z);
      phie[4] = (1. / 8.) * (1 - x) * (1 + z);
      phie[5] = (1. / 8.) * (1 + x) * (1 + z);
      phie[6] = (1. / 8.) * (1 - x) * (1 - z);
      phie[7] = (1. / 8.) * (1 + x) * (1 - z);
// Derivative of phi w.r.t z (N_A,z)
      phis[0] = (1. / 8.) * (1 - x) * (1 - y);
      phis[1] = (1. / 8.) * (1 + x) * (1 - y);
      phis[2] = -(1. / 8.) * (1 - x) * (1 - y);
      phis[3] = -(1. / 8.) * (1 + x) * (1 - y);
      phis[4] = (1. / 8.) * (1 - x) * (1 + y);
      phis[5] = (1. / 8.) * (1 + x) * (1 + y);
      phis[6] = -(1. / 8.) * (1 - x) * (1 + y);
      phis[7] = -(1. / 8.) * (1 + x) * (1 + y);
    }
    
    // bilinear basis functions
    void tsfun_2d(double phi[], 
		  double phic[], 
		  double phie[], 
		  double x, 
		  double y)
    {
      phi[0] = (1. / 4.) * (1 - x) * (1 - y);
      phi[1] = (1. / 4.) * (1 + x) * (1 - y);
      phi[2] = (1. / 4.) * (1 - x) * (1 + y);
      phi[3] = (1. / 4.) * (1 + x) * (1 + y);
      
      phic[0] = -(1. / 4.) * (1 - y);
      phic[1] = (1. / 4.) * (1 - y);
      phic[2] = -(1. / 4.) * (1 + y);
      phic[3] = (1. / 4.) * (1 + y);
      
      phie[0] = -(1. / 4.) * (1 - x);
      phie[1] = -(1. / 4.) * (1 + x);
      phie[2] = (1. / 4.) * (1 - x);
      phie[3] = (1. / 4.) * (1 + x);
    }

    // to multiply two matrices given the matrices and their dimensions. C=A X B
    void matrix_multiply(double * amat, 
			 int arow, 
			 int acol, 
			 double * bmat, 
			 int brow, 
			 int bcol, 
			 double * cmat)
    {

      if (acol != brow)
	std::cerr << "Matrix dimensions do not match, cannot multiply" << std::endl;

      for (int ii = 0; ii < arow * bcol; ii++)
        cmat[ii] = 0.0;

      for (int ii = 0; ii < arow; ii++)
      {
        for (int jj = 0; jj < bcol; jj++)
        {
	  cmat[ii * bcol + jj] = 0.0;

	  for (int kk = 0; kk < acol; kk++)
	  {
	    cmat[ii * bcol + jj] += amat[ii * acol + kk] * bmat[kk * bcol + jj];
	  }
        }
      }
    }

    // to multiply two matrices given the matrices and their dimensions. C=A^T X B
    void matrix_multiply_ATranspose(double * amat, 
			            int arow, 
			            int acol, 
			            double * bmat, 
			            int brow, 
			            int bcol, 
			            double * cmat)
    {
      if (arow != brow)
	std::cerr << "Matrix dimensions do not match, cannot multiply" << std::endl;

      for (int ii = 0; ii < acol * bcol; ii++)
        cmat[ii] = 0.0;

      for (int ii = 0; ii < acol; ii++)
      {
        for (int jj = 0; jj < bcol; jj++)
        {
	  cmat[ii * bcol + jj] = 0.0;

	  for (int kk = 0; kk < arow; kk++)
	  {
	    cmat[ii * bcol + jj] += amat[kk * acol + ii] * bmat[kk * bcol + jj];
	  }
        }
      }
    }

    
    // This function performs Gauss Elimination 
    // I do not use it any more
    void solve_matrix_system(double * A, 
			     double * B, 
			     double * soln, 
			     int order)
    {
      int i, j, k;//, p, u;
      double ratio, sum, maxe, *a, *b;
      
      for (i = 0; i <= order - 1; i++)
      {
        *(soln + i) = 0;
      }
      
      a = new double[order * order];
      b = new double[order];
      
      for (i = 0; i < order; i++)
      {
        for (j = 0; j < order; j++)
        {
	  a[i * order + j] = A[i * order + j];
        }
	
        b[i] = B[i];
      }
      
      for (i = 0; i <= order - 2; i++)
      {
        maxe = *(a + i * order + i);
	
        for (j = i + 1; j < order - 1; j++)
        {
	  if (fabs(*(a + j * order + i)) > fabs(maxe))
	  {
	    maxe = *(a + j * order + i);
	  }
        }
	
        for (j = i + 1; j <= order - 1; j++)
        {
	  ratio = -(*(a + j * order + i) / (*(a + i * order + i)));
	  
	  for (k = i; k <= order - 1; k++)
	  {
	    *(a + j * order + k) += ratio * (*(a + i * order + k));
	  }
	  
	  *(b + j) += ratio * (*(b + i));
        }
      }
      
      for (i = order - 1; i >= 0; i--)
      {
        sum = 0;
	
        for (j = i + 1; j <= order - 1; j++)
        {
	  sum += *(a + i * order + j)* *(soln + j);
        }
	
        *(soln + i) = (*(b + i) - sum) / *(a + i * order + i);
      }
      
      delete [] a;
      delete [] b;
    }


    void calc_ksi_eta_zeta(double xr, 
			   double yr, 
			   double zr, 
			   double u[], 
			   double init_coords_loc[], 
			   double gpx, 
			   double gpy, 
			   double gpz)
{
    int i, iter;
    double F[3], dfdu[9], du[3];    /* F is the vector of the functions, dfdu is the Jacombian matrix, du is the error */
    double x[8], y[8], z[8];        /* u is the solution, x,y,z are the coordinates of the FE in the physical domain */
    double a[3], b[3], c[3], d[3], e[3], f[3], g[3], h[3];
    double eps, norm;              /* tolerance and norm of residuals */

    for (i = 0; i < 8; i++)
    {
        x[i] = init_coords_loc[3 * i];
        y[i] = init_coords_loc[3 * i + 1];
        z[i] = init_coords_loc[3 * i + 2];
    }

    a[0] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7];
    a[1] = y[0] + y[1] + y[2] + y[3] + y[4] + y[5] + y[6] + y[7];
    a[2] = z[0] + z[1] + z[2] + z[3] + z[4] + z[5] + z[6] + z[7];
    b[0] = -x[0] + x[1] - x[2] + x[3] - x[4] + x[5] - x[6] + x[7];
    b[1] = -y[0] + y[1] - y[2] + y[3] - y[4] + y[5] - y[6] + y[7];
    b[2] = -z[0] + z[1] - z[2] + z[3] - z[4] + z[5] - z[6] + z[7];
    c[0] = -x[0] - x[1] - x[2] - x[3] + x[4] + x[5] + x[6] + x[7];
    c[1] = -y[0] - y[1] - y[2] - y[3] + y[4] + y[5] + y[6] + y[7];
    c[2] = -z[0] - z[1] - z[2] - z[3] + z[4] + z[5] + z[6] + z[7];
    d[0] = x[0] - x[1] + x[2] - x[3] - x[4] + x[5] - x[6] + x[7];
    d[1] = y[0] - y[1] + y[2] - y[3] - y[4] + y[5] - y[6] + y[7];
    d[2] = z[0] - z[1] + z[2] - z[3] - z[4] + z[5] - z[6] + z[7];
    e[0] = x[0] + x[1] - x[2] - x[3] + x[4] + x[5] - x[6] - x[7];
    e[1] = y[0] + y[1] - y[2] - y[3] + y[4] + y[5] - y[6] - y[7];
    e[2] = z[0] + z[1] - z[2] - z[3] + z[4] + z[5] - z[6] - z[7];
    f[0] = -x[0] + x[1] + x[2] - x[3] - x[4] + x[5] + x[6] - x[7];
    f[1] = -y[0] + y[1] + y[2] - y[3] - y[4] + y[5] + y[6] - y[7];
    f[2] = -z[0] + z[1] + z[2] - z[3] - z[4] + z[5] + z[6] - z[7];
    g[0] = -x[0] - x[1] + x[2] + x[3] + x[4] + x[5] - x[6] - x[7];
    g[1] = -y[0] - y[1] + y[2] + y[3] + y[4] + y[5] - y[6] - y[7];
    g[2] = -z[0] - z[1] + z[2] + z[3] + z[4] + z[5] - z[6] - z[7];
    h[0] = x[0] - x[1] - x[2] + x[3] - x[4] + x[5] + x[6] - x[7];
    h[1] = y[0] - y[1] - y[2] + y[3] - y[4] + y[5] + y[6] - y[7];
    h[2] = z[0] - z[1] - z[2] + z[3] - z[4] + z[5] + z[6] - z[7];

    /* Newton-Raphson to solve the system */
    /* initial guess */
    u[0] = gpx;
    u[1] = gpy;
    u[2] = gpz;
    norm = 12;
    eps = 0.0000001;
    iter = 0;

    while (norm > eps)
    {
        iter++;
        /*
             if(iter>1000)
        {
          printf("too many iterations in calc_ksi_eta_zeta %d\n",iter);
        }
          */
        /* Calculate Residuals */
        F[0] = -(-8 * xr + a[0] + u[0] * b[0] + u[1] * c[0] + u[0] * u[1] * d[0] + u[2] * e[0] + u[0] * u[2] * f[0] + u[1] * u[2] * g[0] + u[0] * u[1] * u[2] * h[0]);
        F[1] = -(-8 * yr + a[1] + u[0] * b[1] + u[1] * c[1] + u[0] * u[1] * d[1] + u[2] * e[1] + u[0] * u[2] * f[1] + u[1] * u[2] * g[1] + u[0] * u[1] * u[2] * h[1]);
        F[2] = -(-8 * zr + a[2] + u[0] * b[2] + u[1] * c[2] + u[0] * u[1] * d[2] + u[2] * e[2] + u[0] * u[2] * f[2] + u[1] * u[2] * g[2] + u[0] * u[1] * u[2] * h[2]);

        /* Calculate Jacobian */
        dfdu[0] = b[0] + u[1] * d[0] + u[2] * f[0] + u[1] * u[2] * h[0];
        dfdu[1] = c[0] + u[0] * d[0] + u[2] * g[0] + u[0] * u[2] * h[0];
        dfdu[2] = e[0] + u[0] * f[0] + u[1] * g[0] + u[0] * u[1] * h[0];
        dfdu[3] = b[1] + u[1] * d[1] + u[2] * f[1] + u[1] * u[2] * h[1];
        dfdu[4] = c[1] + u[0] * d[1] + u[2] * g[1] + u[0] * u[2] * h[1];
        dfdu[5] = e[1] + u[0] * f[1] + u[1] * g[1] + u[0] * u[1] * h[1];
        dfdu[6] = b[2] + u[1] * d[2] + u[2] * f[2] + u[1] * u[2] * h[2];
        dfdu[7] = c[2] + u[0] * d[2] + u[2] * g[2] + u[0] * u[2] * h[2];
        dfdu[8] = e[2] + u[0] * f[2] + u[1] * g[2] + u[0] * u[1] * h[2];

        /* Solve Matrix System */
        solve_matrix_system(dfdu, F, du, 3);

        /* Update the solution */
        u[0] += du[0];
        u[1] += du[1];
        u[2] += du[2];

        /* calculate the norm */
        norm = sqrt(du[0] * du[0] + du[1] * du[1] + du[2] * du[2]);
    }

    /*  printf("%lf, %lf, %lf,%d\n",u[0],u[1],u[2],iter); */
}


void calc_dos(apf::Vector3 rve[8],
	      double phic[],
	      double phie[],
	      double
	      phis[],
	      double *dett,
	      double ddett[])
{
  int in, n, kk;
  double x1, x2, x3, y1, y2, y3, z1, z2, z3;
  double dx1[24], dx2[24], dx3[24], dy1[24], dy2[24], dy3[24], dz1[24], dz2[24], dz3[24];
  /* isoparametric transformation */
  x1 = 0.0;
  x2 = 0.0;
  x3 = 0.0;
  y1 = 0.0;
  y2 = 0.0;
  y3 = 0.0;
  z1 = 0.0;
  z2 = 0.0;
  z3 = 0.0;
  
  for (in = 0; in < 24; in++)
  {
    dx1[in] = 0;
    dx2[in] = 0;
    dx3[in] = 0;
    dy1[in] = 0;
    dy2[in] = 0;
    dy3[in] = 0;
    dz1[in] = 0;
    dz2[in] = 0;
    dz3[in] = 0;
    ddett[in] = 0;
  }
  
  *dett = 0;
  
  for (n = 0; n < 8; n++)
  {
    double x = rve[n][0];
    double y = rve[n][1];
    double z = rve[n][2];
    x1 = x1 + x * phic[n];
    x2 = x2 + x * phie[n];
    x3 = x3 + x * phis[n];
    y1 = y1 + y * phic[n];
    y2 = y2 + y * phie[n];
    y3 = y3 + y * phis[n];
    z1 = z1 + z * phic[n];
    z2 = z2 + z * phie[n];
    z3 = z3 + z * phis[n];
  }
  
  for (n = 0; n < 8; n++)
  {
    kk = 3 * n;
    dx1[kk] = phic[n];
    dx2[kk] = phie[n];
    dx3[kk] = phis[n];
  }
  
  for (n = 0; n < 8; n++)
  {
    kk = 3 * n + 1;
    dy1[kk] = phic[n];
    dy2[kk] = phie[n];
    dy3[kk] = phis[n];
  }
  
  for (n = 0; n < 8; n++)
  {
    kk = 3 * n + 2;
    dz1[kk] = phic[n];
    dz2[kk] = phie[n];
    dz3[kk] = phis[n];
  }
  
  *dett = x1 * (y2 * z3 - y3 * z2) - y1 * (x2 * z3 - x3 * z2) + z1 * (x2 * y3 - y2 * x3);
  
  for (in = 0; in < 24; in = in + 3)
  {
    ddett[in] = dx1[in] * (y2 * z3 - y3 * z2) - y1 * (dx2[in] * z3 - dx3[in] * z2) + z1 * (dx2[in] * y3 - y2 * dx3[in]);
    ddett[in + 1] = x1 * (dy2[in + 1] * z3 - dy3[in + 1] * z2) - dy1[in + 1] * (x2 * z3 - x3 * z2) + z1 * (x2 * dy3[in + 1] - dy2[in + 1] * x3);
    ddett[in + 2] = x1 * (y2 * dz3[in + 2] - y3 * dz2[in + 2]) - y1 * (x2 * dz3[in + 2] - x3 * dz2[in + 2]) + dz1[in + 2] * (x2 * y3 - y2 * x3);
  }
}

void calc_deriv(double x[], double y[], double z[], double ux[], double uy[], double uz[], double phic[], double phie[], double phis[], double *dudx, double *dudy, double *dudz, double *dvdx, double *dvdy, double *dvdz, double *dgdx, double *dgdy, double *dgdz, double *dett)
{
    double x1(0), x2(0), x3(0), y1(0), y2(0), y3(0), z1(0), z2(0), z3(0);
    double ttky[8], ttkx[8], ttkz[8];
    /* isoparametric transformation */
    *dudx = 0;
    *dudy = 0;
    *dudz = 0;
    *dvdx = 0;
    *dvdy = 0;
    *dvdz = 0;
    *dgdx = 0;
    *dgdy = 0;
    *dgdz = 0;

    for (unsigned short in = 0; in < 8; ++in)
    {
        ttkx[in] = 0;
        ttky[in] = 0;
        ttkz[in] = 0;
    }

    for (unsigned short n = 0; n < 8; n++)
    {
        x1 += x[n] * phic[n];
        x2 += x[n] * phie[n];
        x3 += x[n] * phis[n];
        y1 += y[n] * phic[n];
        y2 += y[n] * phie[n];
        y3 += y[n] * phis[n];
        z1 += z[n] * phic[n];
        z2 += z[n] * phie[n];
        z3 += z[n] * phis[n];
    }

    *dett = x1 * (y2 * z3 - y3 * z2) - y1 * (x2 * z3 - x3 * z2) + z1 * (x2 * y3 - y2 * x3);

    for (unsigned short in = 0; in < 8; ++in)
    {
        ttkx[in] = ((y2 * z3 - z2 * y3) * phic[in] + (z1 * y3 - y1 * z3) * phie[in] + (y1 * z2 - z1 * y2) * phis[in]) / (*dett);
        ttky[in] = ((z2 * x3 - x2 * z3) * phic[in] + (x1 * z3 - z1 * x3) * phie[in] + (z1 * x2 - x1 * z2) * phis[in]) / (*dett);
        ttkz[in] = ((x2 * y3 - y2 * x3) * phic[in] + (y1 * x3 - x1 * y3) * phie[in] + (x1 * y2 - y1 * x2) * phis[in]) / (*dett);
    }

    for (unsigned short in = 0; in < 8; ++in)
    {
        *dudx += ux[in] * ttkx[in];
        *dudy += ux[in] * ttky[in];
        *dudz += ux[in] * ttkz[in];
        *dvdx += uy[in] * ttkx[in];
        *dvdy += uy[in] * ttky[in];
        *dvdz += uy[in] * ttkz[in];
        *dgdx += uz[in] * ttkx[in];
        *dgdy += uz[in] * ttky[in];
        *dgdz += uz[in] * ttkz[in];
    }
}


  double MicroFO::calc_fiber_len_avg()
  {
    // Calculation of the length of the fibers
    int num_elements = fiber_network->numElements();
    double avg = 0.0;
    double len;
    for (int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fiber_network->element(ii);
      int node1 = e.node1_id;
      int node2 = e.node2_id;

      const Node & n1 = fiber_network->node(node1);
      const Node & n2 = fiber_network->node(node2);
      
      len = pow(n2.x - n1.x, 2) + pow(n2.y - n1.y, 2) + pow(n2.z - n1.z, 2);
      len = sqrt(len);

      avg += len;
    }

    avg /= num_elements;
    return avg;
  }

    // T. Stylianopolous, V.H. Barocas, Volume-averaging theory for the study of the mechanics of collangen networks (2006), equation (1)
    void MicroFO::calc_fiber_str(double * len,
				 double * fib_str,
				 double * dfdl)
  {
    int num_elements = fiber_network->numElements();
    for(int ii = 0; ii < num_elements; ii++)
    {
      const Element & fiber = fiber_network->element(ii);
      const FiberReaction * fiber_reaction = fiber_network->getFiberParams(ii);
      
      std::pair<double,double> force_reaction = fiber_reaction->forceReaction(fiber.orig_len,
									      len[ii]);
      fib_str[ii] = force_reaction.first;
      dfdl[ii] = force_reaction.second;
    }
  }

    void MicroFO::calc_force_vector(double * lengths,
				    double * fib_str)
    {
      int num_dofs = fiber_network->numDofs();
      force_vector.assign(num_dofs,0.0);
      
      int num_elements = fiber_network->numElements();
      for (int ii = 0; ii < num_elements; ii++)
      {
	const Element & e = fiber_network->element(ii);
	int node1 = e.node1_id;
        int node2 = e.node2_id;

	int dof0 = node1 * 3;
	int dof1 = node1 * 3 + 1;
	int dof2 = node1 * 3 + 2;

	int dof3 = node2 * 3;
	int dof4 = node2 * 3 + 1;
	int dof5 = node2 * 3 + 2;

	const Node & n1 = fiber_network->node(node1);
	const Node & n2 = fiber_network->node(node2);
	
        double cosalpha = (n2.x - n1.x) / lengths[ii];
        double cosbeta  = (n2.y - n1.y) / lengths[ii];
        double cosgamma = (n2.z - n1.z) / lengths[ii];

	force_vector[dof0] -= fib_str[ii] * cosalpha;
	force_vector[dof1] -= fib_str[ii] * cosbeta;
	force_vector[dof2] -= fib_str[ii] * cosgamma;

	force_vector[dof3] += fib_str[ii] * cosalpha;
	force_vector[dof4] += fib_str[ii] * cosbeta;
	force_vector[dof5] += fib_str[ii] * cosgamma;
      }
    }

    // calculates the norm of a vector
    double calc_norm(double * vector, int size)
    {
      double sum = 0.0;
      for(int ii = 0; ii < size; ii++)
	sum += vector[ii] * vector[ii];
      
      return sqrt(sum);
    }

    // calculation of the entries of the microscopic Jacobian this function is calced by the calc_precond
    double MicroFO::gfunc(int ii, 
			  int flag, 
			  double * lengths, 
			  double * fib_str, 
			  double * dFdL)
    {
      double & dfdl = dFdL[ii];

      const Element & e = fiber_network->element(ii);
      
      int node1 = e.node1_id;
      int node2 = e.node2_id;
      
      double fl = fib_str[ii] / lengths[ii];

      const Node & n1 = fiber_network->node(node1);
      const Node & n2 = fiber_network->node(node2);
      
      double xl = (n2.x - n1.x) / lengths[ii];
      double yl = (n2.y - n1.y) / lengths[ii];
      double z1 = (n2.z - n1.z) / lengths[ii];

      double result = 0.0;

      switch (flag)
      {
      case 1:
        result = (dfdl - fl) * (xl * xl) + fl;
        break;
      case 2:
        result = (dfdl - fl) * (xl * yl);
        break;
      case 4:
        result = (dfdl - fl) * (yl * yl) + fl;
        break;
      case 5:
        result = (dfdl - fl) * (xl * z1);
        break;
      case 6:
        result = (dfdl - fl) * (yl * z1);
        break;
      case 7:
	result = (dfdl - fl) * (z1 * z1) + fl;
        break;
      default:
	std::cerr << "ERROR: STRANGE VALUE ENTERED IN DOPHI_DOETA" << std::endl;
	result = -1.0;
      }

      return result;
    }

    // todo: pull into FiberNetwork, pass in the vector
    void MicroFO::update_nodes()
    {
      int num_nodes = fiber_network->numNodes();
      for(int ii = 0; ii < num_nodes; ii++)
      {
	Node n;
	n.x = coordinate_vector[ii*3    ];
	n.y = coordinate_vector[ii*3 + 1];
	n.z = coordinate_vector[ii*3 + 2];
	if(!TRUSS)
	{
	  n.rx = coordinate_vector[ii*6 + 3];
	  n.ry = coordinate_vector[ii*6 + 4];
	  n.rz = coordinate_vector[ii*6 + 5];	  
	}
	fiber_network->setNode(ii,n);
      }
    }

    // todo: pull into FiberNetwork, pass in/out a vector
    void MicroFO::update_coordinate_vector()
    {
      // sets all values in the vector to 0
      coordinate_vector.assign(coordinate_vector.size(),0.0);
      int num_nodes = fiber_network->numNodes();
      
      for(int ii = 0; ii < num_nodes ; ii++)
      {
	const Node & n = fiber_network->node(ii);
	if(TRUSS)
	{
	  int dof_ = ii * 3;
	  coordinate_vector[dof_]     = n.x;
	  coordinate_vector[dof_ + 1] = n.y;
	  coordinate_vector[dof_ + 2] = n.z;
	}
	else
	{
	  int dof_ = ii * 6;
	  coordinate_vector[dof_]     = n.x;
	  coordinate_vector[dof_ + 1] = n.y;
	  coordinate_vector[dof_ + 2] = n.z;
	  coordinate_vector[dof_ + 3] = n.rx;
	  coordinate_vector[dof_ + 4] = n.ry;
	  coordinate_vector[dof_ + 5] = n.rz;
	}
      }
    }

    void MicroFO::calc_mean_fiber_stretch_omega()
    {
      double mean_fiber_stretch = 0.0;
      double omega_11 = 0.0;
      double total_length = 0.0;
      int num_elements = fiber_network->numElements();
      for(int ii = 0; ii < num_elements; ii++)
      {
	const Element & e = fiber_network->element(ii);
	
	int node1 = e.node1_id;
	int node2 = e.node2_id;

	const Node & n1 = fiber_network->node(node1);
	const Node & n2 = fiber_network->node(node2);
	
	double len = pow((n2.x - n1.x),2) + 
                     pow((n2.y - n1.y),2) + 
                     pow((n2.z - n1.z),2);
	
	len = sqrt(len);
	
	mean_fiber_stretch += len / e.orig_len;
	
	omega_11 += pow((n2.x - n1.x),2) / len;
	total_length += len;
      }

      mean_fiber_stretch /= num_elements;
      omega_11 /= total_length;
  
      AMSI_DEBUG(std::cerr << "mean fiber stretch is " << std::setprecision(16) << mean_fiber_stretch << std::endl);
      AMSI_DEBUG(std::cerr << "omega 11 is " << std::setprecision(16) << omega_11 << std::endl);
}

 } // end of namespace Biotissue
