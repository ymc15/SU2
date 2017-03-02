/*!
 * \file iteration_structure.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/iteration_structure.hpp"

CIteration::CIteration(CConfig *config) { }
CIteration::~CIteration(void) { }

void CIteration::SetGrid_Movement(CGeometry ***geometry_container, 
          CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
          CFreeFormDefBox ***FFDBox,
                            CSolver ****solver_container,
          CConfig **config_container,
                            unsigned short val_iZone,
          unsigned long IntIter,
          unsigned long ExtIter)   {

  unsigned short iDim, iMGlevel, nMGlevels = config_container[val_iZone]->GetnMGLevels();
  unsigned short Kind_Grid_Movement = config_container[val_iZone]->GetKind_GridMovement(val_iZone);
  unsigned long nIterMesh;
  unsigned long iPoint;
  bool stat_mesh = true;
  bool adjoint = config_container[val_iZone]->GetContinuous_Adjoint();
  bool harmonic_balance = (config_container[val_iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE);

  /*--- For a harmonic balance case, set "iteration number" to the zone number,
   so that the meshes are positioned correctly for each instance. ---*/
  if (harmonic_balance) {
    ExtIter = val_iZone;
    Kind_Grid_Movement = config_container[val_iZone]->GetKind_GridMovement(ZONE_0);
  }

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Perform mesh movement depending on specified type ---*/
  switch (Kind_Grid_Movement) {

    case MOVING_WALL:

      /*--- Fixed wall velocities: set the grid velocities only one time
       before the first iteration flow solver. ---*/

      if (ExtIter == 0) {

        if (rank == MASTER_NODE)
          cout << endl << " Setting the moving wall velocities." << endl;

        surface_movement[val_iZone]->Moving_Walls(geometry_container[val_iZone][MESH_0],
                                       config_container[val_iZone], val_iZone, ExtIter);

        /*--- Update the grid velocities on the coarser multigrid levels after
         setting the moving wall velocities for the finest mesh. ---*/

        grid_movement[val_iZone]->UpdateMultiGrid(geometry_container[val_iZone], config_container[val_iZone]);

      }

      break;


    case ROTATING_FRAME:

      /*--- Steadily rotating frame: set the grid velocities just once
       before the first iteration flow solver. ---*/

      if (ExtIter == 0) {

        if (rank == MASTER_NODE) {
          cout << endl << " Setting rotating frame grid velocities";
          cout << " for zone " << val_iZone << "." << endl;
        }

        /*--- Set the grid velocities on all multigrid levels for a steadily
         rotating reference frame. ---*/

        for (iMGlevel = 0; iMGlevel <= nMGlevels; iMGlevel++)
          geometry_container[val_iZone][iMGlevel]->SetRotationalVelocity(config_container[val_iZone], val_iZone);

      }

      break;

    case STEADY_TRANSLATION:

      /*--- Set the translational velocity and hold the grid fixed during
       the calculation (similar to rotating frame, but there is no extra
       source term for translation). ---*/

      if (ExtIter == 0) {

        if (rank == MASTER_NODE)
          cout << endl << " Setting translational grid velocities." << endl;

        /*--- Set the translational velocity on all grid levels. ---*/

        for (iMGlevel = 0; iMGlevel <= nMGlevels; iMGlevel++)
          geometry_container[val_iZone][iMGlevel]->SetTranslationalVelocity(config_container[val_iZone]);

      }

      break;

    case RIGID_MOTION:

      if (rank == MASTER_NODE) {
        cout << endl << " Performing rigid mesh transformation." << endl;
      }

      /*--- Move each node in the volume mesh using the specified type
       of rigid mesh motion. These routines also compute analytic grid
       velocities for the fine mesh. ---*/

      grid_movement[val_iZone]->Rigid_Translation(geometry_container[val_iZone][MESH_0],
                                       config_container[val_iZone], val_iZone, ExtIter);
      grid_movement[val_iZone]->Rigid_Plunging(geometry_container[val_iZone][MESH_0],
                                    config_container[val_iZone], val_iZone, ExtIter);
      grid_movement[val_iZone]->Rigid_Pitching(geometry_container[val_iZone][MESH_0],
                                    config_container[val_iZone], val_iZone, ExtIter);
      grid_movement[val_iZone]->Rigid_Rotation(geometry_container[val_iZone][MESH_0],
                                    config_container[val_iZone], val_iZone, ExtIter);

      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

      grid_movement[val_iZone]->UpdateMultiGrid(geometry_container[val_iZone], config_container[val_iZone]);

      break;

    case DEFORMING:

      if (rank == MASTER_NODE)
        cout << endl << " Updating surface positions." << endl;

      /*--- Translating ---*/

      /*--- Compute the new node locations for moving markers ---*/

      surface_movement[val_iZone]->Surface_Translating(geometry_container[val_iZone][MESH_0],
                                            config_container[val_iZone], ExtIter, val_iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/

      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement[val_iZone]->SetVolume_Deformation(geometry_container[val_iZone][MESH_0],
                                           config_container[val_iZone], true);

      /*--- Plunging ---*/

      /*--- Compute the new node locations for moving markers ---*/

      surface_movement[val_iZone]->Surface_Plunging(geometry_container[val_iZone][MESH_0],
                                         config_container[val_iZone], ExtIter, val_iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/

      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement[val_iZone]->SetVolume_Deformation(geometry_container[val_iZone][MESH_0],
                                           config_container[val_iZone], true);

      /*--- Pitching ---*/

      /*--- Compute the new node locations for moving markers ---*/

      surface_movement[val_iZone]->Surface_Pitching(geometry_container[val_iZone][MESH_0],
                                         config_container[val_iZone], ExtIter, val_iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/

      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement[val_iZone]->SetVolume_Deformation(geometry_container[val_iZone][MESH_0],
                                           config_container[val_iZone], true);

      /*--- Rotating ---*/

      /*--- Compute the new node locations for moving markers ---*/

      surface_movement[val_iZone]->Surface_Rotating(geometry_container[val_iZone][MESH_0],
                                         config_container[val_iZone], ExtIter, val_iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/

      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement[val_iZone]->SetVolume_Deformation(geometry_container[val_iZone][MESH_0],
                                           config_container[val_iZone], true);

      /*--- Update the grid velocities on the fine mesh using finite
       differencing based on node coordinates at previous times. ---*/

      if (!adjoint) {
        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[val_iZone][MESH_0]->SetGridVelocity(config_container[val_iZone], ExtIter);
      }

      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

      grid_movement[val_iZone]->UpdateMultiGrid(geometry_container[val_iZone], config_container[val_iZone]);

      break;

    case EXTERNAL: case EXTERNAL_ROTATION:

      /*--- Apply rigid rotation to entire grid first, if necessary ---*/

      if (Kind_Grid_Movement == EXTERNAL_ROTATION) {
        if (rank == MASTER_NODE)
          cout << " Updating node locations by rigid rotation." << endl;
        grid_movement[val_iZone]->Rigid_Rotation(geometry_container[val_iZone][MESH_0],
                                      config_container[val_iZone], val_iZone, ExtIter);
      }

      /*--- Load new surface node locations from external files ---*/

      if (rank == MASTER_NODE)
        cout << " Updating surface locations from file." << endl;
      surface_movement[val_iZone]->SetExternal_Deformation(geometry_container[val_iZone][MESH_0],
                                                config_container[val_iZone], val_iZone, ExtIter);

      /*--- Deform the volume grid around the new boundary locations ---*/

      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement[val_iZone]->SetVolume_Deformation(geometry_container[val_iZone][MESH_0],
                                           config_container[val_iZone], true);

      /*--- Update the grid velocities on the fine mesh using finite
       differencing based on node coordinates at previous times. ---*/

      if (!adjoint) {
        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[val_iZone][MESH_0]->SetGridVelocity(config_container[val_iZone], ExtIter);
      }

      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

      grid_movement[val_iZone]->UpdateMultiGrid(geometry_container[val_iZone], config_container[val_iZone]);

      break;

    case AEROELASTIC: case AEROELASTIC_RIGID_MOTION:

      /*--- Apply rigid mesh transformation to entire grid first, if necessary ---*/
      if (IntIter == 0) {
        if (Kind_Grid_Movement == AEROELASTIC_RIGID_MOTION) {

          if (rank == MASTER_NODE) {
            cout << endl << " Performing rigid mesh transformation." << endl;
          }

          /*--- Move each node in the volume mesh using the specified type
           of rigid mesh motion. These routines also compute analytic grid
           velocities for the fine mesh. ---*/

          grid_movement[val_iZone]->Rigid_Translation(geometry_container[val_iZone][MESH_0],
                                           config_container[val_iZone], val_iZone, ExtIter);
          grid_movement[val_iZone]->Rigid_Plunging(geometry_container[val_iZone][MESH_0],
                                        config_container[val_iZone], val_iZone, ExtIter);
          grid_movement[val_iZone]->Rigid_Pitching(geometry_container[val_iZone][MESH_0],
                                        config_container[val_iZone], val_iZone, ExtIter);
          grid_movement[val_iZone]->Rigid_Rotation(geometry_container[val_iZone][MESH_0],
                                        config_container[val_iZone], val_iZone, ExtIter);

          /*--- Update the multigrid structure after moving the finest grid,
           including computing the grid velocities on the coarser levels. ---*/

          grid_movement[val_iZone]->UpdateMultiGrid(geometry_container[val_iZone], config_container[val_iZone]);
        }

      }

      /*--- Use the if statement to move the grid only at selected dual time step iterations. ---*/
      else if (IntIter % config_container[val_iZone]->GetAeroelasticIter() == 0) {

        if (rank == MASTER_NODE)
          cout << endl << " Solving aeroelastic equations and updating surface positions." << endl;

        /*--- Solve the aeroelastic equations for the new node locations of the moving markers(surfaces) ---*/

        solver_container[val_iZone][MESH_0][FLOW_SOL]->Aeroelastic(surface_movement[val_iZone], geometry_container[val_iZone][MESH_0], config_container[val_iZone], ExtIter);

        /*--- Deform the volume grid around the new boundary locations ---*/

        if (rank == MASTER_NODE)
          cout << " Deforming the volume grid due to the aeroelastic movement." << endl;
        grid_movement[val_iZone]->SetVolume_Deformation(geometry_container[val_iZone][MESH_0],
                                             config_container[val_iZone], true);

        /*--- Update the grid velocities on the fine mesh using finite
         differencing based on node coordinates at previous times. ---*/

        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[val_iZone][MESH_0]->SetGridVelocity(config_container[val_iZone], ExtIter);

        /*--- Update the multigrid structure after moving the finest grid,
         including computing the grid velocities on the coarser levels. ---*/

        grid_movement[val_iZone]->UpdateMultiGrid(geometry_container[val_iZone], config_container[val_iZone]);
      }

      break;

    case ELASTICITY:

      if (ExtIter != 0) {

        if (rank == MASTER_NODE)
          cout << " Deforming the grid using the Linear Elasticity solution." << endl;

        /*--- Update the coordinates of the grid using the linear elasticity solution. ---*/
        for (iPoint = 0; iPoint < geometry_container[val_iZone][MESH_0]->GetnPoint(); iPoint++) {

          su2double *U_time_nM1 = solver_container[val_iZone][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_time_n1();
          su2double *U_time_n   = solver_container[val_iZone][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_time_n();

          for (iDim = 0; iDim < geometry_container[val_iZone][MESH_0]->GetnDim(); iDim++)
            geometry_container[val_iZone][MESH_0]->node[iPoint]->AddCoord(iDim, U_time_n[iDim] - U_time_nM1[iDim]);

        }

      }

      break;

    case FLUID_STRUCTURE:

      if (rank == MASTER_NODE)
        cout << endl << "Deforming the grid for Fluid-Structure Interaction applications." << endl;

      /*--- Deform the volume grid around the new boundary locations ---*/

      if (rank == MASTER_NODE)
        cout << "Deforming the volume grid." << endl;
      grid_movement[val_iZone]->SetVolume_Deformation(geometry_container[val_iZone][MESH_0],
                                           config_container[val_iZone], true);

      nIterMesh = grid_movement[val_iZone]->Get_nIterMesh();
      stat_mesh = (nIterMesh == 0);

      if (!adjoint && !stat_mesh) {
        if (rank == MASTER_NODE)
          cout << "Computing grid velocities by finite differencing." << endl;
        geometry_container[val_iZone][MESH_0]->SetGridVelocity(config_container[val_iZone], ExtIter);
      }
      else if (stat_mesh) {
          if (rank == MASTER_NODE)
            cout << "The mesh is up-to-date. Using previously stored grid velocities." << endl;
      }

      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

      grid_movement[val_iZone]->UpdateMultiGrid(geometry_container[val_iZone], config_container[val_iZone]);

      break;

    case NO_MOVEMENT: case GUST: default:

      /*--- There is no mesh motion specified for this zone. ---*/
      if (rank == MASTER_NODE)
        cout << "No mesh motion specified." << endl;

      break;
  }

}

void CIteration::Preprocess(COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox,
                            unsigned short val_iZone) { }
void CIteration::Iterate(COutput *output,
                         CIntegration ***integration_container,
                         CGeometry ***geometry_container,
                         CSolver ****solver_container,
                         CNumerics *****numerics_container,
                         CConfig **config_container,
                         CSurfaceMovement **surface_movement,
                         CVolumetricMovement **grid_movement,
                         CFreeFormDefBox*** FFDBox,
                         unsigned short val_iZone) { }
void CIteration::Update(COutput *output,
                        CIntegration ***integration_container,
                        CGeometry ***geometry_container,
                        CSolver ****solver_container,
                        CNumerics *****numerics_container,
                        CConfig **config_container,
                        CSurfaceMovement **surface_movement,
                        CVolumetricMovement **grid_movement,
                        CFreeFormDefBox*** FFDBox,
                        unsigned short val_iZone)      { }
void CIteration::Monitor()     { }
void CIteration::Output()      { }
void CIteration::Postprocess() { }



CMeanFlowIteration::CMeanFlowIteration(CConfig *config) : CIteration(config) { }
CMeanFlowIteration::~CMeanFlowIteration(void) { }

void CMeanFlowIteration::Preprocess(COutput *output,
                                    CIntegration ***integration_container,
                                    CGeometry ***geometry_container,
                                    CSolver ****solver_container,
                                    CNumerics *****numerics_container,
                                    CConfig **config_container,
                                    CSurfaceMovement **surface_movement,
                                    CVolumetricMovement **grid_movement,
                                    CFreeFormDefBox*** FFDBox,
                                    unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[val_iZone]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[val_iZone]->GetExtIter();
  
  bool fsi = config_container[val_iZone]->GetFSI_Simulation();
  unsigned long FSIIter = config_container[val_iZone]->GetFSIIter();

  
  /*--- Set the initial condition for FSI problems with subiterations ---*/
  /*--- This must be done only in the first subiteration ---*/
  if( fsi  && ( FSIIter == 0 ) )
   solver_container[val_iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], ExtIter);

  
  /*--- Apply a Wind Gust ---*/
  
  if (config_container[val_iZone]->GetWind_Gust()) {
    SetWind_GustField(config_container[val_iZone], geometry_container[val_iZone], solver_container[val_iZone]);
  }
  
  
  /*--- Calculate and set Mixing Plane averaged quantities at interfaces ---*/
  
  if(config_container[val_iZone]->GetBoolMixingPlane())
    SetMixingPlane(geometry_container, solver_container, config_container, val_iZone);
  
  /*--- Compute turboperformance ---*/
  
  if(config_container[val_iZone]->GetBoolTurboPerf())
    SetTurboPerformance(geometry_container, solver_container, config_container, output, val_iZone);
}

void CMeanFlowIteration::Iterate(COutput *output,
                                 CIntegration ***integration_container,
                                 CGeometry ***geometry_container,
                                 CSolver ****solver_container,
                                 CNumerics *****numerics_container,
                                 CConfig **config_container,
                                 CSurfaceMovement **surface_movement,
                                 CVolumetricMovement **grid_movement,
                                 CFreeFormDefBox*** FFDBox,
                                 unsigned short val_iZone) {
  
  unsigned long IntIter, ExtIter;
  
  bool unsteady = (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  
  ExtIter = config_container[val_iZone]->GetExtIter();
  
  /* --- Setting up iteration values depending on if this is a
   steady or an unsteady simulaiton */
  
  if ( !unsteady ) IntIter = ExtIter;
  else IntIter = config_container[val_iZone]->GetIntIter();
  
  /*--- Update global parameters ---*/
  
  switch( config_container[val_iZone]->GetKind_Solver() ) {
      
    case EULER: case DISC_ADJ_EULER:
      config_container[val_iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter); break;
      
    case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES:
      config_container[val_iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter); break;
      
    case RANS: case DISC_ADJ_RANS:
      config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter); break;
      
  }
  
  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
  
  integration_container[val_iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                  config_container, RUNTIME_FLOW_SYS, IntIter, val_iZone);
  
  
  if ((config_container[val_iZone]->GetKind_Solver() == RANS) ||
      (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
    
    /*--- Solve the turbulence model ---*/
    
    config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
    integration_container[val_iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                     config_container, RUNTIME_TURB_SYS, IntIter, val_iZone);
    
    /*--- Solve transition model ---*/
    
    if (config_container[val_iZone]->GetKind_Trans_Model() == LM) {
      config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
      integration_container[val_iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                        config_container, RUNTIME_TRANS_SYS, IntIter, val_iZone);
    }
    
  }
  
  /*--- Call Dynamic mesh update if AEROELASTIC motion was specified ---*/
  
  if ((config_container[val_iZone]->GetGrid_Movement()) && (config_container[val_iZone]->GetAeroelastic_Simulation()) && unsteady) {
      
    SetGrid_Movement(geometry_container, surface_movement, grid_movement, FFDBox, solver_container, config_container, val_iZone, IntIter, ExtIter);
    
    /*--- Apply a Wind Gust ---*/
    
    if (config_container[val_iZone]->GetWind_Gust()) {
      if (IntIter % config_container[val_iZone]->GetAeroelasticIter() == 0 && IntIter != 0)
        SetWind_GustField(config_container[val_iZone], geometry_container[val_iZone], solver_container[val_iZone]);
    }
    
  }
  
  
  if ( unsteady && !config_container[val_iZone]->GetDiscrete_Adjoint() )
    
  /*--- Write the convergence history (only screen output) ---*/
    
    output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);
  
}

void CMeanFlowIteration::Update(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                unsigned short val_iZone)      {
  
  unsigned short iMesh;
  su2double Physical_dt, Physical_t;
  unsigned long ExtIter = config_container[val_iZone]->GetExtIter();

  /*--- Dual time stepping strategy ---*/
  
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    /*--- Update dual time solver on all mesh levels ---*/
    
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      integration_container[val_iZone][FLOW_SOL]->SetDualTime_Solver(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh][FLOW_SOL], config_container[val_iZone], iMesh);
      integration_container[val_iZone][FLOW_SOL]->SetConvergence(false);
    }
    
    /*--- Update dual time solver for the turbulence model ---*/
    
    if ((config_container[val_iZone]->GetKind_Solver() == RANS) ||
        (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
      integration_container[val_iZone][TURB_SOL]->SetDualTime_Solver(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0][TURB_SOL], config_container[val_iZone], MESH_0);
      integration_container[val_iZone][TURB_SOL]->SetConvergence(false);
    }
    
    /*--- Update dual time solver for the transition model ---*/
    
    if (config_container[val_iZone]->GetKind_Trans_Model() == LM) {
      integration_container[val_iZone][TRANS_SOL]->SetDualTime_Solver(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0][TRANS_SOL], config_container[val_iZone], MESH_0);
      integration_container[val_iZone][TRANS_SOL]->SetConvergence(false);
    }
    
    /*--- Verify convergence criteria (based on total time) ---*/
    
    Physical_dt = config_container[val_iZone]->GetDelta_UnstTime();
    Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_UnstTime())
      integration_container[val_iZone][FLOW_SOL]->SetConvergence(true);
    
  }
  
}

void CMeanFlowIteration::Monitor()     { }
void CMeanFlowIteration::Output()      { }
void CMeanFlowIteration::Postprocess() { }

void CMeanFlowIteration::SetWind_GustField(CConfig *config_container, CGeometry **geometry_container, CSolver ***solver_container) {
  // The gust is imposed on the flow field via the grid velocities. This method called the Field Velocity Method is described in the
  // NASA TM–2012-217771 - Development, Verification and Use of Gust Modeling in the NASA Computational Fluid Dynamics Code FUN3D
  // the desired gust is prescribed as the negative of the grid velocity.
  
  // If a source term is included to account for the gust field, the method is described by Jones et al. as the Split Velocity Method in
  // Simulation of Airfoil Gust Responses Using Prescribed Velocities.
  // In this routine the gust derivatives needed for the source term are calculated when applicable.
  // If the gust derivatives are zero the source term is also zero.
  // The source term itself is implemented in the class CSourceWindGust
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (rank == MASTER_NODE)
    cout << endl << "Running simulation with a Wind Gust." << endl;
  unsigned short iDim, nDim = geometry_container[MESH_0]->GetnDim(); //We assume nDim = 2
  if (nDim != 2) {
    if (rank == MASTER_NODE) {
      cout << endl << "WARNING - Wind Gust capability is only verified for 2 dimensional simulations." << endl;
    }
  }
  
  /*--- Gust Parameters from config ---*/
  unsigned short Gust_Type = config_container->GetGust_Type();
  su2double xbegin = config_container->GetGust_Begin_Loc();    // Location at which the gust begins.
  su2double L = config_container->GetGust_WaveLength();        // Gust size
  su2double tbegin = config_container->GetGust_Begin_Time();   // Physical time at which the gust begins.
  su2double gust_amp = config_container->GetGust_Ampl();       // Gust amplitude
  su2double n = config_container->GetGust_Periods();           // Number of gust periods
  unsigned short GustDir = config_container->GetGust_Dir(); // Gust direction
  
  /*--- Variables needed to compute the gust ---*/
  unsigned short Kind_Grid_Movement = config_container->GetKind_GridMovement(ZONE_0);
  unsigned long iPoint;
  unsigned short iMGlevel, nMGlevel = config_container->GetnMGLevels();
  
  su2double x, y, x_gust, dgust_dx, dgust_dy, dgust_dt;
  su2double *Gust, *GridVel, *NewGridVel, *GustDer;
  
  su2double Physical_dt = config_container->GetDelta_UnstTime();
  unsigned long ExtIter = config_container->GetExtIter();
  su2double Physical_t = ExtIter*Physical_dt;
  
  su2double Uinf = solver_container[MESH_0][FLOW_SOL]->GetVelocity_Inf(0); // Assumption gust moves at infinity velocity
  
  Gust = new su2double [nDim];
  NewGridVel = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Gust[iDim] = 0.0;
    NewGridVel[iDim] = 0.0;
  }
  
  GustDer = new su2double [3];
  for (unsigned short i = 0; i < 3; i++) {
    GustDer[i] = 0.0;
  }
  
  // Vortex variables
  unsigned long nVortex = 0;
  vector<su2double> x0, y0, vort_strenth, r_core; //vortex is positive in clockwise direction.
  if (Gust_Type == VORTEX) {
    InitializeVortexDistribution(nVortex, x0, y0, vort_strenth, r_core);
  }
  
  /*--- Check to make sure gust lenght is not zero or negative (vortex gust doesn't use this). ---*/
  if (L <= 0.0 && Gust_Type != VORTEX) {
    if (rank == MASTER_NODE) cout << "ERROR: The gust length needs to be positive" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  /*--- Loop over all multigrid levels ---*/
  
  for (iMGlevel = 0; iMGlevel <= nMGlevel; iMGlevel++) {
    
    /*--- Loop over each node in the volume mesh ---*/
    
    for (iPoint = 0; iPoint < geometry_container[iMGlevel]->GetnPoint(); iPoint++) {
      
      /*--- Reset the Grid Velocity to zero if there is no grid movement ---*/
      if (Kind_Grid_Movement == GUST) {
        for (iDim = 0; iDim < nDim; iDim++)
          geometry_container[iMGlevel]->node[iPoint]->SetGridVel(iDim, 0.0);
      }
      
      /*--- initialize the gust and derivatives to zero everywhere ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {Gust[iDim]=0.0;}
      dgust_dx = 0.0; dgust_dy = 0.0; dgust_dt = 0.0;
      
      /*--- Begin applying the gust ---*/
      
      if (Physical_t >= tbegin) {
        
        x = geometry_container[iMGlevel]->node[iPoint]->GetCoord()[0]; // x-location of the node.
        y = geometry_container[iMGlevel]->node[iPoint]->GetCoord()[1]; // y-location of the node.
        
        // Gust coordinate
        x_gust = (x - xbegin - Uinf*(Physical_t-tbegin))/L;
        
        /*--- Calculate the specified gust ---*/
        switch (Gust_Type) {
            
          case TOP_HAT:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp;
              // Still need to put the gust derivatives. Think about this.
            }
            break;
            
          case SINE:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp*(sin(2*PI_NUMBER*x_gust));
              
              // Gust derivatives
              //dgust_dx = gust_amp*2*PI_NUMBER*(cos(2*PI_NUMBER*x_gust))/L;
              //dgust_dy = 0;
              //dgust_dt = gust_amp*2*PI_NUMBER*(cos(2*PI_NUMBER*x_gust))*(-Uinf)/L;
            }
            break;
            
          case ONE_M_COSINE:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp*(1-cos(2*PI_NUMBER*x_gust));
              
              // Gust derivatives
              //dgust_dx = gust_amp*2*PI_NUMBER*(sin(2*PI_NUMBER*x_gust))/L;
              //dgust_dy = 0;
              //dgust_dt = gust_amp*2*PI_NUMBER*(sin(2*PI_NUMBER*x_gust))*(-Uinf)/L;
            }
            break;
            
          case EOG:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = -0.37*gust_amp*sin(3*PI_NUMBER*x_gust)*(1-cos(2*PI_NUMBER*x_gust));
            }
            break;
            
          case VORTEX:
            
            /*--- Use vortex distribution ---*/
            // Algebraic vortex equation.
            for (unsigned long i=0; i<nVortex; i++) {
              su2double r2 = pow(x-(x0[i]+Uinf*(Physical_t-tbegin)), 2) + pow(y-y0[i], 2);
              su2double r = sqrt(r2);
              su2double v_theta = vort_strenth[i]/(2*PI_NUMBER) * r/(r2+pow(r_core[i],2));
              Gust[0] = Gust[0] + v_theta*(y-y0[i])/r;
              Gust[1] = Gust[1] - v_theta*(x-(x0[i]+Uinf*(Physical_t-tbegin)))/r;
            }
            break;
            
          case NONE: default:
            
            /*--- There is no wind gust specified. ---*/
            if (rank == MASTER_NODE) {
              cout << "No wind gust specified." << endl;
            }
            break;
            
        }
      }
      
      /*--- Set the Wind Gust, Wind Gust Derivatives and the Grid Velocities ---*/
      
      GustDer[0] = dgust_dx;
      GustDer[1] = dgust_dy;
      GustDer[2] = dgust_dt;
      
      solver_container[iMGlevel][FLOW_SOL]->node[iPoint]->SetWindGust(Gust);
      solver_container[iMGlevel][FLOW_SOL]->node[iPoint]->SetWindGustDer(GustDer);
      
      GridVel = geometry_container[iMGlevel]->node[iPoint]->GetGridVel();
      
      /*--- Store new grid velocity ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {
        NewGridVel[iDim] = GridVel[iDim] - Gust[iDim];
        geometry_container[iMGlevel]->node[iPoint]->SetGridVel(iDim, NewGridVel[iDim]);
      }
      
    }
  }
  
  delete [] Gust;
  delete [] GustDer;
  delete [] NewGridVel;
  
}

void CMeanFlowIteration::InitializeVortexDistribution(unsigned long &nVortex, vector<su2double>& x0, vector<su2double>& y0, vector<su2double>& vort_strength, vector<su2double>& r_core) {
  /*--- Read in Vortex Distribution ---*/
  std::string line;
  std::ifstream file;
  su2double x_temp, y_temp, vort_strength_temp, r_core_temp;
  file.open("vortex_distribution.txt");
  /*--- In case there is no vortex file ---*/
  if (file.fail()) {
    cout << "There is no vortex data file!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get(); exit(EXIT_FAILURE);
  }
  
  // Ignore line containing the header
  getline(file, line);
  // Read in the information of the vortices (xloc, yloc, lambda(strength), eta(size, gradient))
  while (file.good())
  {
    getline(file, line);
    std::stringstream ss(line);
    if (line.size() != 0) { //ignore blank lines if they exist.
      ss >> x_temp;
      ss >> y_temp;
      ss >> vort_strength_temp;
      ss >> r_core_temp;
      x0.push_back(x_temp);
      y0.push_back(y_temp);
      vort_strength.push_back(vort_strength_temp);
      r_core.push_back(r_core_temp);
    }
  }
  file.close();
  // number of vortices
  nVortex = x0.size();
  
}

void CMeanFlowIteration::SetMixingPlane(CGeometry ***geometry_container, CSolver ****solver_container, CConfig **config_container, unsigned short iZone) {
  
  unsigned short jZone;
  unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  int intMarker, extMarker, intMarkerMix;
  string intMarker_Tag, extMarker_Tag;
  
  /*-- Loop on all the boundary to find MIXING_PLANE boundary --*/
  for (intMarker = 0; intMarker < config_container[iZone]->GetnMarker_All(); intMarker++) {
    for (intMarkerMix=0; intMarkerMix < config_container[iZone]->Get_nMarkerMixingPlane(); intMarkerMix++)
      if (config_container[iZone]->GetMarker_All_TagBound(intMarker) == config_container[iZone]->GetMarker_MixingPlane_Bound(intMarkerMix) ) {
        solver_container[iZone][MESH_0][FLOW_SOL]->Mixing_Process(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], config_container[iZone], intMarker);
        extMarker_Tag = config_container[iZone]->GetMarker_MixingPlane_Donor(intMarkerMix);
        for (jZone = 0; jZone < nZone; jZone++) {
          for (extMarker = 0; extMarker < config_container[jZone]->GetnMarker_All(); extMarker++)
            if (config_container[jZone]->GetMarker_All_TagBound(extMarker) == extMarker_Tag) {
              solver_container[jZone][MESH_0][FLOW_SOL]->SetExtAveragedValue(solver_container[iZone][MESH_0][FLOW_SOL], intMarker, extMarker);
            }
        }
        
      }
    
  }
  
}

void CMeanFlowIteration::SetTurboPerformance(CGeometry ***geometry_container, CSolver ****solver_container, CConfig **config_container, COutput *output, unsigned short iZone) {
  
  unsigned short  jZone, inMarker, outMarker, inMarkerTP, Kind_TurboPerf;
  unsigned short nZone = geometry_container[iZone][MESH_0]->GetnZone();
  string inMarker_Tag, outMarker_Tag;
  
  
  /*-- Loop on all the boundary to find MIXING_PLANE boundary --*/
  for (inMarker = 0; inMarker < config_container[iZone]->GetnMarker_All(); inMarker++)
    for (inMarkerTP=0; inMarkerTP < config_container[iZone]->Get_nMarkerTurboPerf(); inMarkerTP++)
      if (config_container[iZone]->GetMarker_All_TagBound(inMarker) == config_container[iZone]->GetMarker_TurboPerf_BoundIn(inMarkerTP) ) {
        outMarker_Tag =  config_container[iZone]->GetMarker_TurboPerf_BoundOut(inMarkerTP);
        Kind_TurboPerf = config_container[iZone]->GetKind_TurboPerf(inMarkerTP);
        for (jZone = 0; jZone < nZone; jZone++)
          for (outMarker = 0; outMarker < config_container[jZone]->GetnMarker_All(); outMarker++)
            if (config_container[jZone]->GetMarker_All_TagBound(outMarker) == outMarker_Tag) {
              solver_container[iZone][MESH_0][FLOW_SOL]->Mixing_Process(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], config_container[iZone], inMarker);
              solver_container[jZone][MESH_0][FLOW_SOL]->Mixing_Process(geometry_container[jZone][MESH_0], solver_container[jZone][MESH_0], config_container[jZone], outMarker);
              solver_container[iZone][MESH_0][FLOW_SOL]->TurboPerformance(solver_container[jZone][MESH_0][FLOW_SOL], config_container[iZone], inMarker, outMarker, Kind_TurboPerf, inMarkerTP);
              solver_container[ZONE_0][MESH_0][FLOW_SOL]->StoreTurboPerformance(solver_container[iZone][MESH_0][FLOW_SOL], inMarkerTP);
            }
      }
}


CWaveIteration::CWaveIteration(CConfig *config) : CIteration(config) { }
CWaveIteration::~CWaveIteration(void) { }
void CWaveIteration::Preprocess(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                unsigned short val_iZone) { }
void CWaveIteration::Iterate(COutput *output,
                             CIntegration ***integration_container,
                             CGeometry ***geometry_container,
                             CSolver ****solver_container,
                             CNumerics *****numerics_container,
                             CConfig **config_container,
                             CSurfaceMovement **surface_movement,
                             CVolumetricMovement **grid_movement,
                             CFreeFormDefBox*** FFDBox,
                             unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Set the value of the internal iteration ---*/
  IntIter = ExtIter;
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
  
  /*--- Wave equations ---*/
  config_container[val_iZone]->SetGlobalParam(WAVE_EQUATION, RUNTIME_WAVE_SYS, ExtIter);
  integration_container[val_iZone][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                   config_container, RUNTIME_WAVE_SYS, IntIter, val_iZone);
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    for (IntIter = 1; IntIter < config_container[val_iZone]->GetUnst_nIntIter(); IntIter++) {
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);
      config_container[val_iZone]->SetIntIter(IntIter);
      integration_container[val_iZone][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                       config_container, RUNTIME_WAVE_SYS, IntIter, val_iZone);
      if (integration_container[val_iZone][WAVE_SOL]->GetConvergence()) break;
    }
    
  }
  
}
void CWaveIteration::Update(COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox,
                            unsigned short val_iZone)      {
  
  unsigned short iMesh;
  su2double Physical_dt, Physical_t;
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    /*--- Update dual time solver ---*/
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      integration_container[val_iZone][WAVE_SOL]->SetDualTime_Solver(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh][WAVE_SOL], config_container[val_iZone], iMesh);
      integration_container[val_iZone][WAVE_SOL]->SetConvergence(false);
    }
    
    Physical_dt = config_container[val_iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_UnstTime()) integration_container[val_iZone][WAVE_SOL]->SetConvergence(true);
  }
}

void CWaveIteration::Monitor()     { }
void CWaveIteration::Output()      { }
void CWaveIteration::Postprocess() { }


CHeatIteration::CHeatIteration(CConfig *config) : CIteration(config) { }
CHeatIteration::~CHeatIteration(void) { }
void CHeatIteration::Preprocess(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                unsigned short val_iZone) { }
void CHeatIteration::Iterate(COutput *output,
                             CIntegration ***integration_container,
                             CGeometry ***geometry_container,
                             CSolver ****solver_container,
                             CNumerics *****numerics_container,
                             CConfig **config_container,
                             CSurfaceMovement **surface_movement,
                             CVolumetricMovement **grid_movement,
                             CFreeFormDefBox*** FFDBox,
                             unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Set the value of the internal iteration ---*/
  IntIter = ExtIter;
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
  
  /*--- Heat equation ---*/
  config_container[val_iZone]->SetGlobalParam(HEAT_EQUATION, RUNTIME_HEAT_SYS, ExtIter);
  integration_container[val_iZone][HEAT_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                   config_container, RUNTIME_HEAT_SYS, IntIter, val_iZone);
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    for (IntIter = 1; IntIter < config_container[val_iZone]->GetUnst_nIntIter(); IntIter++) {
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);
      config_container[val_iZone]->SetIntIter(IntIter);
      integration_container[val_iZone][HEAT_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                       config_container, RUNTIME_HEAT_SYS, IntIter, val_iZone);
      if (integration_container[val_iZone][HEAT_SOL]->GetConvergence()) break;
    }
  }
}

void CHeatIteration::Update(COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox,
                            unsigned short val_iZone)      {
  
  unsigned short iMesh;
  su2double Physical_dt, Physical_t;
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    /*--- Update dual time solver ---*/
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      integration_container[val_iZone][HEAT_SOL]->SetDualTime_Solver(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh][HEAT_SOL], config_container[val_iZone], iMesh);
      integration_container[val_iZone][HEAT_SOL]->SetConvergence(false);
    }
    
    Physical_dt = config_container[val_iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_UnstTime()) integration_container[val_iZone][HEAT_SOL]->SetConvergence(true);
  }
}
void CHeatIteration::Monitor()     { }
void CHeatIteration::Output()      { }
void CHeatIteration::Postprocess() { }


CPoissonIteration::CPoissonIteration(CConfig *config) : CIteration(config) { }
CPoissonIteration::~CPoissonIteration(void) { }
void CPoissonIteration::Preprocess(COutput *output,
                                   CIntegration ***integration_container,
                                   CGeometry ***geometry_container,
                                   CSolver ****solver_container,
                                   CNumerics *****numerics_container,
                                   CConfig **config_container,
                                   CSurfaceMovement **surface_movement,
                                   CVolumetricMovement **grid_movement,
                                   CFreeFormDefBox*** FFDBox,
                                   unsigned short val_iZone) { }
void CPoissonIteration::Iterate(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Set the value of the internal iteration ---*/
  IntIter = ExtIter;
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
  
  /*--- Poisson equation ---*/
  config_container[val_iZone]->SetGlobalParam(POISSON_EQUATION, RUNTIME_POISSON_SYS, ExtIter);
  integration_container[val_iZone][POISSON_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                      config_container, RUNTIME_POISSON_SYS, IntIter, val_iZone);
  
  
}
void CPoissonIteration::Update(COutput *output,
                               CIntegration ***integration_container,
                               CGeometry ***geometry_container,
                               CSolver ****solver_container,
                               CNumerics *****numerics_container,
                               CConfig **config_container,
                               CSurfaceMovement **surface_movement,
                               CVolumetricMovement **grid_movement,
                               CFreeFormDefBox*** FFDBox,
                               unsigned short val_iZone)      { }
void CPoissonIteration::Monitor()     { }
void CPoissonIteration::Output()      { }
void CPoissonIteration::Postprocess() { }


CFEM_StructuralAnalysis::CFEM_StructuralAnalysis(CConfig *config) : CIteration(config) { }
CFEM_StructuralAnalysis::~CFEM_StructuralAnalysis(void) { }
void CFEM_StructuralAnalysis::Preprocess() { }
void CFEM_StructuralAnalysis::Iterate(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                  unsigned short val_iZone
                                ) {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  su2double loadIncrement;
  unsigned long IntIter = 0; config_container[val_iZone]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[val_iZone]->GetExtIter();

  bool fsi = config_container[val_iZone]->GetFSI_Simulation();

  unsigned long iIncrement;
  unsigned long nIncrements = config_container[val_iZone]->GetNumberIncrements();

  bool nonlinear = (config_container[val_iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Geometrically non-linear problems
  bool linear = (config_container[val_iZone]->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Geometrically non-linear problems

  bool initial_calc = config_container[val_iZone]->GetExtIter() == 0;        // Checks if it is the first calculation.
  bool first_iter = config_container[val_iZone]->GetIntIter() == 0;        // Checks if it is the first iteration
  bool restart = config_container[val_iZone]->GetRestart();                        // Restart analysis
  bool initial_calc_restart = (SU2_TYPE::Int(config_container[val_iZone]->GetExtIter()) == config_container[val_iZone]->GetDyn_RestartIter()); // Initial calculation for restart

  su2double CurrentTime = config_container[val_iZone]->GetCurrent_DynTime();
  su2double Static_Time = config_container[val_iZone]->GetStatic_Time();

  bool statTime = (CurrentTime <= Static_Time);

  bool incremental_load = config_container[val_iZone]->GetIncrementalLoad();              // If an incremental load is applied

  /*--- This is to prevent problems when running a linear solver ---*/
  if (!nonlinear) incremental_load = false;

  /*--- Set the convergence monitor to false, to prevent the solver to stop in intermediate FSI subiterations ---*/
  integration_container[val_iZone][FEA_SOL]->SetConvergence(false);

  if (linear) {

    /*--- Set the value of the internal iteration ---*/

    IntIter = ExtIter;

    /*--- FEA equations ---*/

    config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

    /*--- Run the iteration ---*/

    integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
        config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);

  }
  /*--- If the structure is held static and the solver is nonlinear, we don't need to solve for static time, but we need to compute Mass Matrix and Integration constants ---*/
  else if ((nonlinear) && ((!statTime) || (!fsi))) {

    /*--- THIS IS THE DIRECT APPROACH (NO INCREMENTAL LOAD APPLIED) ---*/

    if (!incremental_load) {

      /*--- Set the value of the internal iteration ---*/

      IntIter = 0;

      /*--- FEA equations ---*/

      config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

      /*--- Run the iteration ---*/

      integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
          config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);


      /*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

      for (IntIter = 1; IntIter < config_container[val_iZone]->GetDyn_nIntIter(); IntIter++) {

        /*--- Write the convergence history (only screen output) ---*/

        output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

        config_container[val_iZone]->SetIntIter(IntIter);

        integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
            config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);

        if (integration_container[val_iZone][FEA_SOL]->GetConvergence()) break;

      }

    }
    /*--- The incremental load is only used in nonlinear cases ---*/
    else if (incremental_load) {

      /*--- Set the initial condition: store the current solution as Solution_Old ---*/

      solver_container[val_iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], ExtIter);

      /*--- The load increment is 1.0 ---*/
      loadIncrement = 1.0;
      solver_container[val_iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);

      /*--- Set the value of the internal iteration ---*/

      IntIter = 0;

      /*--- FEA equations ---*/

      config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

      /*--- Run the first iteration ---*/

      integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
          config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);


      /*--- Write the convergence history (only screen output) ---*/

      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

      /*--- Run the second iteration ---*/

      IntIter = 1;

      config_container[val_iZone]->SetIntIter(IntIter);

      integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
          config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);


      bool meetCriteria;
      su2double Residual_UTOL, Residual_RTOL, Residual_ETOL;
      su2double Criteria_UTOL, Criteria_RTOL, Criteria_ETOL;

      Criteria_UTOL = config_container[val_iZone]->GetIncLoad_Criteria(0);
      Criteria_RTOL = config_container[val_iZone]->GetIncLoad_Criteria(1);
      Criteria_ETOL = config_container[val_iZone]->GetIncLoad_Criteria(2);

      Residual_UTOL = log10(solver_container[val_iZone][MESH_0][FEA_SOL]->GetRes_FEM(0));
      Residual_RTOL = log10(solver_container[val_iZone][MESH_0][FEA_SOL]->GetRes_FEM(1));
      Residual_ETOL = log10(solver_container[val_iZone][MESH_0][FEA_SOL]->GetRes_FEM(2));

      meetCriteria = ( ( Residual_UTOL <  Criteria_UTOL ) &&
          ( Residual_RTOL <  Criteria_RTOL ) &&
          ( Residual_ETOL <  Criteria_ETOL ) );

      /*--- If the criteria is met and the load is not "too big", do the regular calculation ---*/
      if (meetCriteria) {

        for (IntIter = 2; IntIter < config_container[val_iZone]->GetDyn_nIntIter(); IntIter++) {

          /*--- Write the convergence history (only screen output) ---*/

          output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

          config_container[val_iZone]->SetIntIter(IntIter);

          integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
              config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);

          if (integration_container[val_iZone][FEA_SOL]->GetConvergence()) break;

        }

      }

      /*--- If the criteria is not met, a whole set of subiterations for the different loads must be done ---*/

      else {

        /*--- Here we have to restart the solution to the original one of the iteration ---*/
        /*--- Retrieve the Solution_Old as the current solution before subiterating ---*/

        solver_container[val_iZone][MESH_0][FEA_SOL]->ResetInitialCondition(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], ExtIter);

        /*--- For the number of increments ---*/
        for (iIncrement = 0; iIncrement < nIncrements; iIncrement++) {

          loadIncrement = (iIncrement + 1.0) * (1.0 / nIncrements);

          /*--- Set the load increment and the initial condition, and output the parameters of UTOL, RTOL, ETOL for the previous iteration ---*/

          /*--- Set the convergence monitor to false, to force se solver to converge every subiteration ---*/
          integration_container[val_iZone][FEA_SOL]->SetConvergence(false);

          output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

          /*--- FEA equations ---*/

          config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);


          solver_container[val_iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);

          if (rank == MASTER_NODE) {
            cout << endl;
            cout << "-- Incremental load: increment " << iIncrement + 1 << " ------------------------------------------" << endl;
          }

          /*--- Set the value of the internal iteration ---*/
          IntIter = 0;
          config_container[val_iZone]->SetIntIter(IntIter);

          /*--- FEA equations ---*/

          config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

          /*--- Run the iteration ---*/

          integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
              config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);


          /*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

          for (IntIter = 1; IntIter < config_container[val_iZone]->GetDyn_nIntIter(); IntIter++) {

            /*--- Write the convergence history (only screen output) ---*/

            output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

            config_container[val_iZone]->SetIntIter(IntIter);

            integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);

            if (integration_container[val_iZone][FEA_SOL]->GetConvergence()) break;

          }

        }

      }

    }


  }
  else if (
      (nonlinear && statTime) &&
      ((first_iter && initial_calc) || (restart && initial_calc_restart))
  ) {

    /*--- We need to do the preprocessing to compute the Mass Matrix and integration constants ---*/
    solver_container[val_iZone][MESH_0][FEA_SOL]->Preprocessing(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0],
        config_container[val_iZone], numerics_container[val_iZone][MESH_0][FEA_SOL], MESH_0, 0, RUNTIME_FEA_SYS, false);

  }

}

void CFEM_StructuralAnalysis::Update(COutput *output,
       CIntegration ***integration_container,
       CGeometry ***geometry_container,
       CSolver ****solver_container,
       CNumerics *****numerics_container,
       CConfig **config_container,
       CSurfaceMovement **surface_movement,
       CVolumetricMovement **grid_movement,
       CFreeFormDefBox*** FFDBox,
       unsigned short val_iZone) {

  su2double Physical_dt, Physical_t;
    unsigned long ExtIter = config_container[val_iZone]->GetExtIter();
  bool dynamic = (config_container[val_iZone]->GetDynamic_Analysis() == DYNAMIC);          // Dynamic problems

  /*----------------- Compute averaged nodal stress and reactions ------------------------*/

  solver_container[val_iZone][MESH_0][FEA_SOL]->Compute_NodalStress(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0], numerics_container[val_iZone][MESH_0][FEA_SOL], config_container[val_iZone]);

  /*----------------- Update structural solver ----------------------*/

  if (dynamic) {
    integration_container[val_iZone][FEA_SOL]->SetFEM_StructuralSolver(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0], config_container[val_iZone], MESH_0);
    integration_container[val_iZone][FEA_SOL]->SetConvergence(false);

      /*--- Verify convergence criteria (based on total time) ---*/

    Physical_dt = config_container[val_iZone]->GetDelta_DynTime();
    Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_DynTime())
      integration_container[val_iZone][FEA_SOL]->SetConvergence(true);
  }

}
void CFEM_StructuralAnalysis::Monitor()     { }
void CFEM_StructuralAnalysis::Output()      { }
void CFEM_StructuralAnalysis::Postprocess() { }


CAdjMeanFlowIteration::CAdjMeanFlowIteration(CConfig *config) : CIteration(config) { }
CAdjMeanFlowIteration::~CAdjMeanFlowIteration(void) { }
void CAdjMeanFlowIteration::Preprocess(COutput *output,
                                       CIntegration ***integration_container,
                                       CGeometry ***geometry_container,
                                       CSolver ****solver_container,
                                       CNumerics *****numerics_container,
                                       CConfig **config_container,
                                       CSurfaceMovement **surface_movement,
                                       CVolumetricMovement **grid_movement,
                                       CFreeFormDefBox*** FFDBox,
                                       unsigned short val_iZone) {
  
  unsigned short iMesh;
  bool harmonic_balance = (config_container[ZONE_0]->GetUnsteady_Simulation() == HARMONIC_BALANCE);
  bool dynamic_mesh = config_container[ZONE_0]->GetGrid_Movement();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- For the unsteady adjoint, load a new direct solution from a restart file. ---*/
  
  if (((dynamic_mesh && ExtIter == 0) || config_container[val_iZone]->GetUnsteady_Simulation()) && !harmonic_balance) {
    int Direct_Iter = SU2_TYPE::Int(config_container[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 1;
    if (rank == MASTER_NODE && val_iZone == ZONE_0 && config_container[val_iZone]->GetUnsteady_Simulation())
      cout << endl << " Loading flow solution from direct iteration " << Direct_Iter << "." << endl;
    solver_container[val_iZone][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], Direct_Iter);
  }
  
  /*--- Continuous adjoint Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations ---*/
  
  if ((ExtIter == 0) || config_container[val_iZone]->GetUnsteady_Simulation()) {
    
    if (config_container[val_iZone]->GetKind_Solver() == ADJ_EULER)
      config_container[val_iZone]->SetGlobalParam(ADJ_EULER, RUNTIME_FLOW_SYS, ExtIter);
    if (config_container[val_iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES)
      config_container[val_iZone]->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
    if (config_container[val_iZone]->GetKind_Solver() == ADJ_RANS)
      config_container[val_iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_FLOW_SYS, ExtIter);
    
    /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
    
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "Begin direct solver to store flow data (single iteration)." << endl;
    
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
    
    integration_container[val_iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                    config_container, RUNTIME_FLOW_SYS, 0, val_iZone);
    
    if (config_container[val_iZone]->GetKind_Solver() == ADJ_RANS) {
      
      /*--- Solve the turbulence model ---*/
      
      config_container[val_iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_TURB_SYS, ExtIter);
      integration_container[val_iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                       config_container, RUNTIME_TURB_SYS, IntIter, val_iZone);
      
      /*--- Solve transition model ---*/
      
      if (config_container[val_iZone]->GetKind_Trans_Model() == LM) {
        config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
        integration_container[val_iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                          config_container, RUNTIME_TRANS_SYS, IntIter, val_iZone);
      }
      
    }
    
    /*--- Output the residual (visualization purpouses to identify if
     the direct solution is converged)---*/
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "log10[Maximum residual]: " << log10(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(0))
      <<", located at point "<< solver_container[val_iZone][MESH_0][FLOW_SOL]->GetPoint_Max(0) << "." << endl;
    
    /*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
     note that in the direct Euler problem we are not computing the gradients of the primitive variables ---*/
    
    if (config_container[val_iZone]->GetKind_Gradient_Method() == GREEN_GAUSS)
      solver_container[val_iZone][MESH_0][FLOW_SOL]->SetPrimitive_Gradient_GG(geometry_container[val_iZone][MESH_0], config_container[val_iZone]);
    if (config_container[val_iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
      solver_container[val_iZone][MESH_0][FLOW_SOL]->SetPrimitive_Gradient_LS(geometry_container[val_iZone][MESH_0], config_container[val_iZone]);
    
    /*--- Set contribution from cost function for boundary conditions ---*/
    
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      
      /*--- Set the value of the non-dimensional coefficients in the coarse levels, using the fine level solution ---*/
      
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetTotal_CD(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CD());
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetTotal_CL(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CL());
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetTotal_CT(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CT());
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetTotal_CQ(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CQ());
      
      /*--- Compute the adjoint boundary condition on Euler walls ---*/
      
      solver_container[val_iZone][iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh], config_container[val_iZone]);
      
      /*--- Set the internal boundary condition on nearfield surfaces ---*/
      
      if ((config_container[val_iZone]->GetKind_ObjFunc() == EQUIVALENT_AREA) ||
          (config_container[val_iZone]->GetKind_ObjFunc() == NEARFIELD_PRESSURE))
        solver_container[val_iZone][iMesh][ADJFLOW_SOL]->SetIntBoundary_Jump(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh], config_container[val_iZone]);
      
    }
    
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "End direct solver, begin adjoint problem." << endl;
    
  }
  
}
void CAdjMeanFlowIteration::Iterate(COutput *output,
                                    CIntegration ***integration_container,
                                    CGeometry ***geometry_container,
                                    CSolver ****solver_container,
                                    CNumerics *****numerics_container,
                                    CConfig **config_container,
                                    CSurfaceMovement **surface_movement,
                                    CVolumetricMovement **grid_movement,
                                    CFreeFormDefBox*** FFDBox,
                                    unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  bool unsteady = (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  
  /*--- Set the value of the internal iteration ---*/
  
  ExtIter = config_container[val_iZone]->GetExtIter();
  
  /* --- Setting up iteration values depending on if this is a 
  steady or an unsteady simulaiton */

  if ( !unsteady ) 
    IntIter = ExtIter;
  else
    IntIter = config_container[val_iZone]->GetIntIter();
    
    
  switch( config_container[val_iZone]->GetKind_Solver() ) {

  case ADJ_EULER:
    config_container[val_iZone]->SetGlobalParam(ADJ_EULER, RUNTIME_ADJFLOW_SYS, ExtIter); break;

  case ADJ_NAVIER_STOKES:
    config_container[val_iZone]->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_ADJFLOW_SYS, ExtIter); break;

  case ADJ_RANS:
    config_container[val_iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_ADJFLOW_SYS, ExtIter); break;          
  }
    
  /*--- Iteration of the flow adjoint problem ---*/
  
  integration_container[val_iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                     config_container, RUNTIME_ADJFLOW_SYS, IntIter, val_iZone);
  
  /*--- Iteration of the turbulence model adjoint ---*/
  
  if ((config_container[val_iZone]->GetKind_Solver() == ADJ_RANS) && (!config_container[val_iZone]->GetFrozen_Visc())) {
    
    /*--- Adjoint turbulence model solution ---*/
    
    config_container[val_iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_ADJTURB_SYS, ExtIter);
    integration_container[val_iZone][ADJTURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                        config_container, RUNTIME_ADJTURB_SYS, IntIter, val_iZone);
    
  }
  
  /*--- Dual time stepping strategy ---*/
  
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    for (IntIter = 1; IntIter < config_container[val_iZone]->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Write the convergence history (only screen output) ---*/
      
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);
      
      /*--- Set the value of the internal iteration ---*/
      
      config_container[val_iZone]->SetIntIter(IntIter);
      
      /*--- All zones must be advanced and coupled with each pseudo timestep ---*/
      
      integration_container[val_iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                         config_container, RUNTIME_ADJFLOW_SYS, IntIter, val_iZone);
      
      /*--- Check to see if the convergence criteria has been met ---*/
      
      if (integration_container[val_iZone][ADJFLOW_SOL]->GetConvergence()) break;
    }
    
  }
  
}
void CAdjMeanFlowIteration::Update(COutput *output,
                                   CIntegration ***integration_container,
                                   CGeometry ***geometry_container,
                                   CSolver ****solver_container,
                                   CNumerics *****numerics_container,
                                   CConfig **config_container,
                                   CSurfaceMovement **surface_movement,
                                   CVolumetricMovement **grid_movement,
                                   CFreeFormDefBox*** FFDBox,
                                   unsigned short val_iZone)      {
  
  su2double Physical_dt, Physical_t;
  unsigned short iMesh;
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Dual time stepping strategy ---*/
  
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    /*--- Update dual time solver ---*/
    
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      integration_container[val_iZone][ADJFLOW_SOL]->SetDualTime_Solver(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh][ADJFLOW_SOL], config_container[val_iZone], iMesh);
      integration_container[val_iZone][ADJFLOW_SOL]->SetConvergence(false);
    }
    
    Physical_dt = config_container[val_iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_UnstTime()) integration_container[val_iZone][ADJFLOW_SOL]->SetConvergence(true);
    
  }
}

void CAdjMeanFlowIteration::Monitor()     { }
void CAdjMeanFlowIteration::Output()      { }
void CAdjMeanFlowIteration::Postprocess() { }

CDiscAdjMeanFlowIteration::CDiscAdjMeanFlowIteration(CConfig *config) : CIteration(config), CurrentRecording(NONE) {
  
  meanflow_iteration = new CMeanFlowIteration(config);
  
  turbulent = config->GetKind_Solver() == DISC_ADJ_RANS;
  
}

CDiscAdjMeanFlowIteration::~CDiscAdjMeanFlowIteration(void) { }
void CDiscAdjMeanFlowIteration::Preprocess(COutput *output,
                                           CIntegration ***integration_container,
                                           CGeometry ***geometry_container,
                                           CSolver ****solver_container,
                                           CNumerics *****numerics_container,
                                           CConfig **config_container,
                                           CSurfaceMovement **surface_movement,
                                           CVolumetricMovement **grid_movement,
                                           CFreeFormDefBox*** FFDBox,
                                           unsigned short val_iZone) {
  
  unsigned long IntIter = 0, iPoint;
  config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned short ExtIter = config_container[val_iZone]->GetExtIter();
  bool unsteady = config_container[val_iZone]->GetUnsteady_Simulation() != NONE;
  bool dual_time_1st = (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  unsigned short iMesh;
  int Direct_Iter;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- For the unsteady adjoint, load direct solutions from restart files. ---*/

  if (config_container[val_iZone]->GetUnsteady_Simulation()) {

    Direct_Iter = SU2_TYPE::Int(config_container[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 2;

    /*--- For dual-time stepping we want to load the already converged solution at timestep n ---*/

    if (dual_time) {
      Direct_Iter += 1;
    }

    if (dual_time_2nd) {

      /*--- Load solution at timestep n-2 ---*/

      LoadUnsteady_Solution(geometry_container, solver_container,config_container, val_iZone, Direct_Iter-2);

      /*--- Push solution back to correct array ---*/

      for (iMesh=0; iMesh<=config_container[val_iZone]->GetnMGLevels();iMesh++) {
        for(iPoint=0; iPoint<geometry_container[val_iZone][iMesh]->GetnPoint();iPoint++) {
          solver_container[val_iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          solver_container[val_iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
          if (turbulent) {
            solver_container[val_iZone][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
            solver_container[val_iZone][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
          }
        }
      }
    }
    if (dual_time) {

      /*--- Load solution at timestep n-1 ---*/

      LoadUnsteady_Solution(geometry_container, solver_container,config_container, val_iZone, Direct_Iter-1);

      /*--- Push solution back to correct array ---*/

      for (iMesh=0; iMesh<=config_container[val_iZone]->GetnMGLevels();iMesh++) {
        for(iPoint=0; iPoint<geometry_container[val_iZone][iMesh]->GetnPoint();iPoint++) {
          solver_container[val_iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          if (turbulent) {
            solver_container[val_iZone][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          }
        }
      }
    }

    /*--- Load solution timestep n ---*/

    LoadUnsteady_Solution(geometry_container, solver_container,config_container, val_iZone, Direct_Iter);


    /*--- Store flow solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry_container[val_iZone][MESH_0]->GetnPoint(); iPoint++) {
      solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->node[iPoint]->SetSolution_Direct(solver_container[val_iZone][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution());
    }
    if (turbulent) {
      for (iPoint = 0; iPoint < geometry_container[val_iZone][MESH_0]->GetnPoint(); iPoint++) {
        solver_container[val_iZone][MESH_0][ADJTURB_SOL]->node[iPoint]->SetSolution_Direct(solver_container[val_iZone][MESH_0][TURB_SOL]->node[iPoint]->GetSolution());
      }
    }
  }

  solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->Preprocessing(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0],  config_container[val_iZone] , MESH_0, 0, RUNTIME_ADJFLOW_SYS, false);
  if (turbulent) {
    solver_container[val_iZone][MESH_0][ADJTURB_SOL]->Preprocessing(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0],  config_container[val_iZone] , MESH_0, 0, RUNTIME_ADJTURB_SYS, false);
  }

  if (CurrentRecording != FLOW_VARIABLES || unsteady) {
    
    if (rank == MASTER_NODE && ((ExtIter == 0) || unsteady )) {
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << "Direct iteration to store computational graph." << endl;
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
      cout << "-------------------------------------------------------------------------" << endl << endl;
    }
    
    /*--- Record one mean flow iteration with flow variables as input ---*/
    
    SetRecording(output, integration_container, geometry_container, solver_container, numerics_container,
                 config_container, surface_movement, grid_movement, FFDBox, val_iZone, FLOW_VARIABLES);
    
    /*--- Print residuals in the first iteration ---*/
    
    if (rank == MASTER_NODE && ((ExtIter == 0) || unsteady )) {
      cout << "log10[RMS Density]: "<< log10(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_RMS(0))
           <<", CD: " <<solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CD()
          <<", CL: " << solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CL() << "." << endl;

      if (turbulent) {
        cout << "log10[RMS k]: " << log10(solver_container[val_iZone][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
      }
    }
    
  }
  
}



void CDiscAdjMeanFlowIteration::LoadUnsteady_Solution(CGeometry ***geometry_container,
                                           CSolver ****solver_container,
                                           CConfig **config_container,
                                           unsigned short val_iZone, int val_DirectIter) {
  unsigned short iMesh;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (val_DirectIter >= 0) {
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << " Loading flow solution from direct iteration " << val_DirectIter  << "." << endl;
    solver_container[val_iZone][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], val_DirectIter);
    solver_container[val_iZone][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[val_iZone][MESH_0],solver_container[val_iZone][MESH_0], config_container[val_iZone], MESH_0, val_DirectIter, RUNTIME_FLOW_SYS, false);
    if (turbulent) {
      solver_container[val_iZone][MESH_0][TURB_SOL]->LoadRestart(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], val_DirectIter);
      solver_container[val_iZone][MESH_0][TURB_SOL]->Postprocessing(geometry_container[val_iZone][MESH_0],solver_container[val_iZone][MESH_0], config_container[val_iZone], MESH_0);
    }
  } else {
    /*--- If there is no solution file we set the freestream condition ---*/
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << " Setting freestream conditions at direct iteration " << val_DirectIter << "." << endl;
    for (iMesh=0; iMesh<=config_container[val_iZone]->GetnMGLevels();iMesh++) {
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetFreeStream_Solution(config_container[val_iZone]);
      solver_container[val_iZone][iMesh][FLOW_SOL]->Preprocessing(geometry_container[val_iZone][iMesh],solver_container[val_iZone][iMesh], config_container[val_iZone], iMesh, val_DirectIter, RUNTIME_FLOW_SYS, false);
      if (turbulent) {
        solver_container[val_iZone][iMesh][TURB_SOL]->SetFreeStream_Solution(config_container[val_iZone]);
        solver_container[val_iZone][iMesh][TURB_SOL]->Postprocessing(geometry_container[val_iZone][iMesh],solver_container[val_iZone][iMesh], config_container[val_iZone], iMesh);
      }
    }
  }
}


void CDiscAdjMeanFlowIteration::Iterate(COutput *output,
                                        CIntegration ***integration_container,
                                        CGeometry ***geometry_container,
                                        CSolver ****solver_container,
                                        CNumerics *****numerics_container,
                                        CConfig **config_container,
                                        CSurfaceMovement **surface_movement,
                                        CVolumetricMovement **volume_grid_movement,
                                        CFreeFormDefBox*** FFDBox,
                                        unsigned short val_iZone) {
  
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  unsigned long IntIter=0, nIntIter = 1;
  bool dual_time_1st = (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);

  config_container[val_iZone]->SetIntIter(IntIter);

  if(dual_time)
    nIntIter = config_container[val_iZone]->GetUnst_nIntIter();

  for(IntIter=0; IntIter< nIntIter; IntIter++) {

    /*--- Set the internal iteration ---*/

    config_container[val_iZone]->SetIntIter(IntIter);

    /*--- Set the adjoint values of the flow and objective function ---*/

    InitializeAdjoint(solver_container, geometry_container, config_container, val_iZone);

    /*--- Run the adjoint computation ---*/

    AD::ComputeAdjoint();

    /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

    solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry_container[val_iZone][MESH_0],
                                                                              config_container[val_iZone]);

    solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry_container[val_iZone][MESH_0],
                                                                               config_container[val_iZone]);

    if (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS) {
      solver_container[val_iZone][MESH_0][ADJTURB_SOL]->ExtractAdjoint_Solution(geometry_container[val_iZone][MESH_0],
                                                                                config_container[val_iZone]);
    }

    /*--- Clear all adjoints to re-use the stored computational graph in the next iteration ---*/

    AD::ClearAdjoints();

    /*--- Set the convergence criteria (only residual possible) ---*/

    integration_container[val_iZone][ADJFLOW_SOL]->Convergence_Monitoring(geometry_container[val_iZone][MESH_0],config_container[val_iZone],
                                                                          IntIter,log10(solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0)), MESH_0);

    if(integration_container[val_iZone][ADJFLOW_SOL]->GetConvergence()) {
      break;
    }

    /*--- Write the convergence history (only screen output) ---*/

    if(dual_time && (IntIter != nIntIter-1))
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

  }

  if (dual_time) {
    integration_container[val_iZone][ADJFLOW_SOL]->SetConvergence(false);
  }
  
  if (((ExtIter+1 >= config_container[val_iZone]->GetnExtIter()) || (integration_container[val_iZone][ADJFLOW_SOL]->GetConvergence()) ||
      ((ExtIter % config_container[val_iZone]->GetWrt_Sol_Freq() == 0))) || (dual_time)) {
    
    /*--- Record one mean flow iteration with geometry variables as input ---*/
    
    SetRecording(output, integration_container, geometry_container, solver_container, numerics_container,
                 config_container, surface_movement, volume_grid_movement, FFDBox, val_iZone, GEOMETRY_VARIABLES);
    
    /*--- Set the adjoint values of the flow and objective function ---*/
    
    InitializeAdjoint(solver_container, geometry_container, config_container, val_iZone);
    
    /*--- Run the adjoint computation ---*/
    
    AD::ComputeAdjoint();
    
    /*--- Extract the sensitivities (adjoint of node coordinates) ---*/
    
    solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[val_iZone][MESH_0],config_container[val_iZone]);
    
  }
  
}

void CDiscAdjMeanFlowIteration::SetRecording(COutput *output,
                                             CIntegration ***integration_container,
                                             CGeometry ***geometry_container,
                                             CSolver ****solver_container,
                                             CNumerics *****numerics_container,
                                             CConfig **config_container,
                                             CSurfaceMovement **surface_movement,
                                             CVolumetricMovement **grid_movement,
                                             CFreeFormDefBox*** FFDBox,
                                             unsigned short val_iZone,
                                             unsigned short kind_recording)      {
  
  unsigned long IntIter = config_container[ZONE_0]->GetIntIter();
  unsigned long ExtIter = config_container[val_iZone]->GetExtIter(), DirectExtIter;
  bool unsteady = config_container[val_iZone]->GetUnsteady_Simulation() != NONE;
  unsigned short iMesh;

  DirectExtIter = 0;
  if (unsteady) {
    DirectExtIter = SU2_TYPE::Int(config_container[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 1;
  }

  /*--- Reset the tape ---*/

  AD::Reset();

  /*--- We only need to reset the indices if the current recording is different from the recording we want to have ---*/

  if (CurrentRecording != kind_recording && (CurrentRecording != NONE) ) {

    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      solver_container[val_iZone][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[val_iZone][MESH_0], config_container[val_iZone], kind_recording);
    }

    if (turbulent) {
      solver_container[val_iZone][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[val_iZone][MESH_0], config_container[val_iZone], kind_recording);
    }

    /*--- Clear indices of coupling variables ---*/

    SetDependencies(solver_container, geometry_container, config_container, val_iZone, ALL_VARIABLES);

    /*--- Run one iteration while tape is passive - this clears all indices ---*/

    meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                                config_container,surface_movement,grid_movement,FFDBox,val_iZone);

  }

    /*--- Prepare for recording ---*/

    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      solver_container[val_iZone][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[val_iZone][MESH_0], config_container[val_iZone], kind_recording);
    }

    if (turbulent) {
      solver_container[val_iZone][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[val_iZone][MESH_0], config_container[val_iZone], kind_recording);
    }

  /*--- Start the recording of all operations ---*/
  
  AD::StartRecording();
  
  /*--- Register flow variables ---*/
  
  RegisterInput(solver_container, geometry_container, config_container, val_iZone, kind_recording);
  
  /*--- Compute coupling or update the geometry ---*/

  SetDependencies(solver_container, geometry_container, config_container, val_iZone, kind_recording);
  
  /*--- Set the correct direct iteration number ---*/

  if (unsteady) {
    config_container[val_iZone]->SetExtIter(DirectExtIter);
  }

  /*--- Run the direct iteration ---*/

  config_container[val_iZone]->SetIntIter(0);
  meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                              config_container,surface_movement,grid_movement,FFDBox, val_iZone);

  config_container[val_iZone]->SetExtIter(ExtIter);

  /*--- Register flow variables and objective function as output ---*/
  
  /*--- For flux-avg or area-avg objective functions the 1D values must be calculated first ---*/
  for (unsigned short iObj=0; iObj<config_container[val_iZone]->GetnObj(); iObj++) {
    if (config_container[val_iZone]->GetKind_ObjFunc(iObj)==AVG_OUTLET_PRESSURE ||
        config_container[val_iZone]->GetKind_ObjFunc(iObj)==AVG_TOTAL_PRESSURE ||
        config_container[val_iZone]->GetKind_ObjFunc(iObj)==MASS_FLOW_RATE) {
      output->OneDimensionalOutput(solver_container[val_iZone][MESH_0][FLOW_SOL],
                                   geometry_container[val_iZone][MESH_0], config_container[val_iZone]);
      break;
    }
  }

  /*--- For a combined objective function, the total should be computed and stored ---*/
  if (config_container[val_iZone]->GetnObj()>0) {
    solver_container[val_iZone][MESH_0][FLOW_SOL]->Compute_ComboObj(config_container[val_iZone]);
  }

  RegisterOutput(solver_container, geometry_container, config_container, val_iZone);
  
  /*--- Stop the recording ---*/
  
  AD::StopRecording();
  
  /*--- Set the recording status ---*/
  
  CurrentRecording = kind_recording;

  /*--- Reset the number of the internal iterations---*/

  config_container[ZONE_0]->SetIntIter(IntIter);

}


void CDiscAdjMeanFlowIteration::RegisterInput(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone, unsigned short kind_recording) {
  
  
  if (kind_recording == FLOW_VARIABLES) {
    
    /*--- Register flow and turbulent variables as input ---*/
    
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0], config_container[iZone]);
    
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterVariables(geometry_container[iZone][MESH_0], config_container[iZone]);
    
    if (turbulent) {
      solver_container[iZone][MESH_0][ADJTURB_SOL]->RegisterSolution(geometry_container[iZone][MESH_0], config_container[iZone]);
    }
  }
  if (kind_recording == GEOMETRY_VARIABLES) {
    
    /*--- Register node coordinates as input ---*/
    
    geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);
    
  }

}

void CDiscAdjMeanFlowIteration::SetDependencies(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone, unsigned short kind_recording) {

  if ((kind_recording == GEOMETRY_VARIABLES) || (kind_recording == ALL_VARIABLES)) {

    /*--- Update geometry to get the influence on other geometry variables (normals, volume etc) ---*/
    
    geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone], config_container[iZone]);
    
  }

  /*--- Compute coupling between flow and turbulent equations ---*/

  if (turbulent) {
    solver_container[iZone][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[iZone][MESH_0],solver_container[iZone][MESH_0], config_container[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

    solver_container[iZone][MESH_0][TURB_SOL]->Postprocessing(geometry_container[iZone][MESH_0],solver_container[iZone][MESH_0], config_container[iZone], MESH_0);
  }

}

void CDiscAdjMeanFlowIteration::RegisterOutput(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone) {
  
  /*--- Register objective function as output of the iteration ---*/
  
  solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
  
  /*--- Register conservative variables as output of the iteration ---*/
  
  solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],config_container[iZone]);
  
  if (turbulent) {
    solver_container[iZone][MESH_0][ADJTURB_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                 config_container[iZone]);
  }
}

void CDiscAdjMeanFlowIteration::InitializeAdjoint(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone) {
  
  /*--- Initialize the adjoint of the objective function (typically with 1.0) ---*/
  
  solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone]);
  
  /*--- Initialize the adjoints the conservative variables ---*/
  
  solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                  config_container[iZone]);
  
  if (turbulent) {
    solver_container[iZone][MESH_0][ADJTURB_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                    config_container[iZone]);
  }
}
void CDiscAdjMeanFlowIteration::Update(COutput *output,
                                       CIntegration ***integration_container,
                                       CGeometry ***geometry_container,
                                       CSolver ****solver_container,
                                       CNumerics *****numerics_container,
                                       CConfig **config_container,
                                       CSurfaceMovement **surface_movement,
                                       CVolumetricMovement **grid_movement,
                                       CFreeFormDefBox*** FFDBox,
                                       unsigned short val_iZone)      { }
void CDiscAdjMeanFlowIteration::Monitor()     { }
void CDiscAdjMeanFlowIteration::Output()      { }
void CDiscAdjMeanFlowIteration::Postprocess() { }


void FEM_StructuralIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                                 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                                 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {

  su2double Physical_dt, Physical_t;
  su2double loadIncrement;
  unsigned short iZone;
  unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
    unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

    unsigned long iIncrement;
    unsigned long nIncrements = config_container[ZONE_0]->GetNumberIncrements();

  bool dynamic = (config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC);          // Dynamic problems
  bool nonlinear = (config_container[ZONE_0]->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Geometrically non-linear problems

  bool incremental_load = config_container[ZONE_0]->GetIncrementalLoad();              // If an incremental load is applied

  /*--- This is to prevent problems when running a linear solver ---*/
  if (!nonlinear) incremental_load = false;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


  /*--- THIS IS THE DIRECT APPROACH (NO INCREMENTAL LOAD APPLIED) ---*/

  if (!incremental_load) {

    /*--- Set the initial condition ---*/

//    for (iZone = 0; iZone < nZone; iZone++)
//      solver_container[iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

    for (iZone = 0; iZone < nZone; iZone++) {

      /*--- Set the value of the internal iteration ---*/

      IntIter = ExtIter;
      if (nonlinear) IntIter = 0;

      /*--- FEA equations ---*/

      config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

      /*--- Run the iteration ---*/

      integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                                                                      config_container, RUNTIME_FEA_SYS, IntIter, iZone);



    }

    /*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

    if (nonlinear) {
      for (IntIter = 1; IntIter < config_container[ZONE_0]->GetDyn_nIntIter(); IntIter++) {

        for (iZone = 0; iZone < nZone; iZone++) {

          /*--- Write the convergence history (only screen output) ---*/

          output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

          config_container[iZone]->SetIntIter(IntIter);

          integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                                                                          config_container, RUNTIME_FEA_SYS, IntIter, iZone);

        }

        if (integration_container[ZONE_0][FEA_SOL]->GetConvergence()) break;

      }

    }

  }
  /*--- The incremental load is only used in nonlinear cases ---*/
  else if (incremental_load) {

    /*--- Set the initial condition: store the current solution as Solution_Old ---*/

    for (iZone = 0; iZone < nZone; iZone++)
      solver_container[iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

    for (iZone = 0; iZone < nZone; iZone++) {

        /*--- The load increment is 1.0 ---*/
        loadIncrement = 1.0;
        solver_container[iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);

        /*--- Set the value of the internal iteration ---*/

        IntIter = 0;

        /*--- FEA equations ---*/

        config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

        /*--- Run the first iteration ---*/

        integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                                                                        config_container, RUNTIME_FEA_SYS, IntIter, iZone);


        /*--- Write the convergence history (only screen output) ---*/

        output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

        /*--- Run the second iteration ---*/

        IntIter = 1;

        config_container[iZone]->SetIntIter(IntIter);

        integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                                                                        config_container, RUNTIME_FEA_SYS, IntIter, iZone);

    }

    bool meetCriteria;
    su2double Residual_UTOL, Residual_RTOL, Residual_ETOL;
    su2double Criteria_UTOL, Criteria_RTOL, Criteria_ETOL;

    Criteria_UTOL = config_container[ZONE_0]->GetIncLoad_Criteria(0);
    Criteria_RTOL = config_container[ZONE_0]->GetIncLoad_Criteria(1);
    Criteria_ETOL = config_container[ZONE_0]->GetIncLoad_Criteria(2);

    Residual_UTOL = log10(solver_container[ZONE_0][MESH_0][FEA_SOL]->GetRes_FEM(0));
    Residual_RTOL = log10(solver_container[ZONE_0][MESH_0][FEA_SOL]->GetRes_FEM(1));
    Residual_ETOL = log10(solver_container[ZONE_0][MESH_0][FEA_SOL]->GetRes_FEM(2));

    meetCriteria = ( ( Residual_UTOL <  Criteria_UTOL ) &&
               ( Residual_RTOL <  Criteria_RTOL ) &&
             ( Residual_ETOL <  Criteria_ETOL ) );

    /*--- If the criteria is met and the load is not "too big", do the regular calculation ---*/
    if (meetCriteria) {

      for (IntIter = 2; IntIter < config_container[ZONE_0]->GetDyn_nIntIter(); IntIter++) {

        for (iZone = 0; iZone < nZone; iZone++) {

        /*--- Write the convergence history (only screen output) ---*/

        output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

        config_container[iZone]->SetIntIter(IntIter);

        integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                                        config_container, RUNTIME_FEA_SYS, IntIter, iZone);

        }

        if (integration_container[ZONE_0][FEA_SOL]->GetConvergence()) break;

      }

    }

    /*--- If the criteria is not met, a whole set of subiterations for the different loads must be done ---*/

    else {

      /*--- Here we have to restart the solution to the original one of the iteration ---*/
      /*--- Retrieve the Solution_Old as the current solution before subiterating ---*/

      for (iZone = 0; iZone < nZone; iZone++)
        solver_container[iZone][MESH_0][FEA_SOL]->ResetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

      /*--- For the number of increments ---*/
      for (iIncrement = 0; iIncrement < nIncrements; iIncrement++) {

        loadIncrement = (iIncrement + 1.0) * (1.0 / nIncrements);

        /*--- Set the load increment and the initial condition, and output the parameters of UTOL, RTOL, ETOL for the previous iteration ---*/

        for (iZone = 0; iZone < nZone; iZone++) {

          /*--- Set the convergence monitor to false, to force se solver to converge every subiteration ---*/
          integration_container[iZone][FEA_SOL]->SetConvergence(false);

          output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

          /*--- FEA equations ---*/

          config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);


          solver_container[iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);
        }

        if (rank == MASTER_NODE) {
          cout << endl;
          cout << "-- Incremental load: increment " << iIncrement + 1 << " ------------------------------------------" << endl;
        }

        for (iZone = 0; iZone < nZone; iZone++) {

          /*--- Set the value of the internal iteration ---*/
          IntIter = 0;
          config_container[iZone]->SetIntIter(IntIter);

          /*--- FEA equations ---*/

          config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

          /*--- Run the iteration ---*/

          integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                                                                          config_container, RUNTIME_FEA_SYS, IntIter, iZone);



        }

        /*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

        for (IntIter = 1; IntIter < config_container[ZONE_0]->GetDyn_nIntIter(); IntIter++) {

          for (iZone = 0; iZone < nZone; iZone++) {

            /*--- Write the convergence history (only screen output) ---*/

            output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

            config_container[iZone]->SetIntIter(IntIter);

            integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                                                                              config_container, RUNTIME_FEA_SYS, IntIter, iZone);

          }

          if (integration_container[ZONE_0][FEA_SOL]->GetConvergence()) break;

        }

      }

    }

  }



  /*----------------- Compute averaged nodal stress and reactions ------------------------*/

  for (iZone = 0; iZone < nZone; iZone++)
    solver_container[iZone][MESH_0][FEA_SOL]->Compute_NodalStress(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], numerics_container[iZone][MESH_0][FEA_SOL], config_container[iZone]);

  /*----------------- Update structural solver ----------------------*/

  if (dynamic) {
    for (iZone = 0; iZone < nZone; iZone++) {
      integration_container[iZone][FEA_SOL]->SetFEM_StructuralSolver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], config_container[iZone], MESH_0);
      integration_container[iZone][FEA_SOL]->SetConvergence(false);
    }

      /*--- Verify convergence criteria (based on total time) ---*/

    Physical_dt = config_container[ZONE_0]->GetDelta_DynTime();
    Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[ZONE_0]->GetTotal_DynTime())
      integration_container[ZONE_0][FEA_SOL]->SetConvergence(true);
  }


}
