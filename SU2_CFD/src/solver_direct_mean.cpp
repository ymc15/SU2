/*!
 * \file solution_direct_mean.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solver_structure.hpp"

CEulerSolver::CEulerSolver(void) : CSolver() {
  
  /*--- Basic array initialization ---*/
  
  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  
  CD_Mnt = NULL; CL_Mnt = NULL; CSF_Mnt = NULL;  CEff_Mnt = NULL;
  CMx_Mnt = NULL; CMy_Mnt = NULL; CMz_Mnt = NULL;
  CFx_Mnt = NULL; CFy_Mnt = NULL; CFz_Mnt = NULL;

  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  ForceMomentum = NULL; MomentMomentum = NULL;

  /*--- Surface based array initialization ---*/
  
  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;
  
  Surface_CL_Mnt = NULL; Surface_CD_Mnt = NULL; Surface_CSF_Mnt = NULL; Surface_CEff_Mnt = NULL;
  Surface_CFx_Mnt = NULL; Surface_CFy_Mnt = NULL; Surface_CFz_Mnt = NULL;
  Surface_CMx_Mnt = NULL; Surface_CMy_Mnt = NULL; Surface_CMz_Mnt = NULL;
  
  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;
  
  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;
  
  CMerit_Mnt = NULL;  CT_Mnt = NULL;  CQ_Mnt = NULL;

  /*--- Supersonic simulation array initialization ---*/
  
  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;
  
  /*--- Engine simulation array initialization ---*/
  
  Inflow_MassFlow = NULL;   Inflow_Pressure = NULL;
  Inflow_Mach = NULL;       Inflow_Area = NULL;
  Exhaust_Pressure = NULL;  Exhaust_Temperature = NULL;
  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;
  
  /*--- Numerical methods array initialization ---*/
  
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  DonorPrimVar = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT = NULL;

  Smatrix = NULL; Cvector = NULL;
 
  Secondary = NULL; Secondary_i = NULL; Secondary_j = NULL;

  /*--- Fixed CL mode initialization (cauchy criteria) ---*/
  
  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;
  

  SlidingState = NULL;
  
  FluidModel = NULL;

  AveragedVelocity = NULL;
  AveragedNormal   = NULL;
  AveragedGridVel  = NULL;
  AveragedFlux     = NULL;
  TotalFlux        = NULL;

  AveragedNormalVelocity     = NULL;
  AveragedTangVelocity       = NULL;
  ExtAveragedNormalVelocity  = NULL;
  ExtAveragedTangVelocity    = NULL;
  MassFlow                   = NULL;
  FlowAngle                  = NULL;
  AveragedEnthalpy           = NULL;
  AveragedPressure           = NULL;
  AveragedTotPressure        = NULL;
  AveragedTotTemperature     = NULL;
  ExtAveragedTotPressure     = NULL;
  ExtAveragedTotTemperature  = NULL;
  AveragedDensity            = NULL;
  ExtAveragedPressure        = NULL;
  ExtAveragedDensity         = NULL;
  AveragedSoundSpeed         = NULL;
  AveragedEntropy            = NULL;
  AveragedTangGridVelocity   = NULL;
  AveragedMach               = NULL;
  AveragedNormalMach         = NULL;
  AveragedTangMach           = NULL;

  TotalStaticEfficiency = NULL;
  TotalTotalEfficiency  = NULL;
  KineticEnergyLoss     = NULL;
  TotalPressureLoss     = NULL;
  MassFlowIn            = NULL;
  MassFlowOut           = NULL; 
  FlowAngleIn           = NULL;
  FlowAngleOut          = NULL;
  EulerianWork          = NULL;
  TotalEnthalpyIn       = NULL;
  PressureRatio         = NULL;
  PressureOut           = NULL;
  EnthalpyOut           = NULL;
  MachIn                = NULL;
  MachOut               = NULL;
  NormalMachIn          = NULL;
  NormalMachOut         = NULL;
  VelocityOutIs         = NULL;
 
}

CEulerSolver::CEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  unsigned long iPoint, index, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  su2double StaticEnergy, Density, Velocity2, Pressure, Temperature, dull_val;
  int Unst_RestartIter;
  ifstream restart_file;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
    bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  su2double AoA_, AoS_, BCThrust_;
  string filename = config->GetSolution_FlowFileName();
  string filename_ = config->GetSolution_FlowFileName();
  string::size_type position;
  unsigned long ExtIter_;
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));
  unsigned short direct_diff = config->GetDirectDiff();
  unsigned short nMarkerTurboPerf = config->Get_nMarkerTurboPerf();
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/
  
  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {
    
    /*--- Multizone problems require the number of the zone to be appended. ---*/
    
    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone);
    
    /*--- Modify file name for a dual-time unsteady restart ---*/
    
    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }
    
    /*--- Modify file name for a time stepping unsteady restart ---*/
    
    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Open the restart file, throw an error if this fails. ---*/
    
    restart_file.open(filename_.data(), ios::in);
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file!! " << filename_.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }
    
    unsigned long iPoint_Global = 0;
    string text_line;
    
    /*--- The first line is the header (General description) ---*/
    
    getline (restart_file, text_line);
    
    /*--- Space for the solution ---*/
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
    }
    
    /*--- Space for extra info (if any) ---*/
    
    while (getline (restart_file, text_line)) {
      
      /*--- Angle of attack ---*/
      
      position = text_line.find ("AOA=",0);
      if (position != string::npos) {
        text_line.erase (0,4); AoA_ = atof(text_line.c_str());
        if (config->GetDiscard_InFiles() == false) {
          if ((config->GetAoA() != AoA_) &&  (rank == MASTER_NODE)) {
            cout.precision(6);
            cout << fixed <<"WARNING: AoA in the solution file (" << AoA_ << " deg.) +" << endl;
            cout << "         AoA offset in mesh file (" << config->GetAoA_Offset() << " deg.) = " << AoA_ + config->GetAoA_Offset() << " deg." << endl;
          }
          config->SetAoA(AoA_ + config->GetAoA_Offset());
        }
        else {
          if ((config->GetAoA() != AoA_) &&  (rank == MASTER_NODE))
            cout <<"WARNING: Discarding the AoA in the solution file." << endl;
        }
      }
      
      /*--- Sideslip angle ---*/
      
      position = text_line.find ("SIDESLIP_ANGLE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); AoS_ = atof(text_line.c_str());
        if (config->GetDiscard_InFiles() == false) {
          if ((config->GetAoS() != AoS_) &&  (rank == MASTER_NODE)) {
            cout.precision(6);
            cout << fixed <<"WARNING: AoS in the solution file (" << AoS_ << " deg.) +" << endl;
            cout << "         AoS offset in mesh file (" << config->GetAoS_Offset() << " deg.) = " << AoS_ + config->GetAoS_Offset() << " deg." << endl;
          }
          config->SetAoS(AoS_ + config->GetAoS_Offset());
        }
        else {
          if ((config->GetAoS() != AoS_) &&  (rank == MASTER_NODE))
            cout <<"WARNING: Discarding the AoS in the solution file." << endl;
        }
      }
      
      /*--- BCThrust angle ---*/
      
      position = text_line.find ("INITIAL_BCTHRUST=",0);
      if (position != string::npos) {
        text_line.erase (0,17); BCThrust_ = atof(text_line.c_str());
        if (config->GetDiscard_InFiles() == false) {
          if ((config->GetInitial_BCThrust() != BCThrust_) &&  (rank == MASTER_NODE))
            cout <<"WARNING: ACDC will use the initial BC Thrust provided in the solution file: " << BCThrust_ << " lbs." << endl;
          config->SetInitial_BCThrust(BCThrust_);
        }
        else {
          if ((config->GetInitial_BCThrust() != BCThrust_) &&  (rank == MASTER_NODE))
            cout <<"WARNING: Discarding the BC Thrust in the solution file." << endl;
        }
      }
      
      /*--- External iteration ---*/
      
      position = text_line.find ("EXT_ITER=",0);
      if (position != string::npos) {
        text_line.erase (0,9); ExtIter_ = atoi(text_line.c_str());
        if (!config->GetContinuous_Adjoint() && !config->GetDiscrete_Adjoint())
          config->SetExtIter_OffSet(ExtIter_);
      }
      
    }
    
    /*--- Close the restart file... we will open this file again... ---*/
    
    restart_file.close();
    
  }

  /*--- Array initialization ---*/
  
  /*--- Basic array initialization ---*/

  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;

  CD_Mnt= NULL; CL_Mnt= NULL; CSF_Mnt= NULL; CEff_Mnt= NULL;
  CMx_Mnt= NULL;   CMy_Mnt= NULL;   CMz_Mnt= NULL;
  CFx_Mnt= NULL;   CFy_Mnt= NULL;   CFz_Mnt= NULL;

  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  ForceMomentum = NULL;  MomentMomentum = NULL;

  /*--- Surface based array initialization ---*/
  
  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;
  
  Surface_CL_Mnt= NULL; Surface_CD_Mnt= NULL; Surface_CSF_Mnt= NULL; Surface_CEff_Mnt= NULL;
  Surface_CFx_Mnt= NULL;   Surface_CFy_Mnt= NULL;   Surface_CFz_Mnt= NULL;
  Surface_CMx_Mnt= NULL;   Surface_CMy_Mnt= NULL;   Surface_CMz_Mnt = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Rotorcraft simulation array initialization ---*/

  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;

  CMerit_Mnt = NULL; CT_Mnt = NULL; CQ_Mnt = NULL;

  /*--- Supersonic simulation array initialization ---*/
  
  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;
  
  /*--- Engine simulation array initialization ---*/
  
  Inflow_MassFlow = NULL;   Inflow_Pressure = NULL;
  Inflow_Mach = NULL;       Inflow_Area = NULL;
  Exhaust_Pressure = NULL;  Exhaust_Temperature = NULL;
  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;
  
  /*--- Numerical methods array initialization ---*/
  
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  DonorPrimVar = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT = NULL;

  Smatrix = NULL; Cvector = NULL;

  Secondary=NULL; Secondary_i=NULL; Secondary_j=NULL;

  /*--- Fixed CL mode initialization (cauchy criteria) ---*/

  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;
  
  /*--- Fluid model pointer initialization ---*/

  FluidModel = NULL;

  /*--- Turbo array initialization ---*/

  AveragedVelocity = NULL;
  AveragedNormal   = NULL;
  AveragedGridVel  = NULL;
  AveragedFlux     = NULL;
  TotalFlux        = NULL;

  AveragedNormalVelocity     = NULL;
  AveragedTangVelocity       = NULL;
  ExtAveragedNormalVelocity  = NULL;
  ExtAveragedTangVelocity    = NULL;
  MassFlow                   = NULL;
  FlowAngle                  = NULL;
  AveragedEnthalpy           = NULL;
  AveragedPressure           = NULL;
  AveragedTotPressure        = NULL;
  AveragedTotTemperature     = NULL;
  ExtAveragedTotPressure     = NULL;
  ExtAveragedTotTemperature  = NULL;
  AveragedDensity            = NULL;
  ExtAveragedPressure        = NULL;
  ExtAveragedDensity         = NULL;
  AveragedSoundSpeed         = NULL;
  AveragedEntropy            = NULL;
  AveragedTangGridVelocity   = NULL;
  AveragedMach               = NULL;
  AveragedNormalMach         = NULL;
  AveragedTangMach           = NULL;

  TotalStaticEfficiency = NULL;
  TotalTotalEfficiency  = NULL;
  KineticEnergyLoss     = NULL;
  TotalPressureLoss     = NULL;
  MassFlowIn            = NULL; 
  MassFlowOut           = NULL;
  FlowAngleIn           = NULL;
  FlowAngleOut          = NULL;
  EulerianWork          = NULL;
  TotalEnthalpyIn       = NULL;
  PressureRatio         = NULL;
  PressureOut           = NULL;
  EnthalpyOut           = NULL;
  MachIn                = NULL;
  MachOut               = NULL;
  NormalMachIn          = NULL;
  NormalMachOut         = NULL;
  VelocityOutIs         = NULL;

  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables (T, vx, vy, vz, P, rho, h, c, lamMu, EddyMu, ThCond, Cp).
   ---*/
  
  nDim = geometry->GetnDim();

  nVar = nDim+2;
  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  nSecondaryVar = 2; nSecondaryVarGrad = 2;

  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nPrimVarGrad;
  
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
 
  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) 
    nVertex[iMarker] = geometry->nVertex[iMarker];
 
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  SetNondimensionalization(geometry, config, iMesh);
  
  /*--- Allocate the node variables ---*/
  
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  
  Primitive   = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the Secondary solution ---*/
  
  Secondary   = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
  Secondary_i = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
  Secondary_j = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  
  if (roe_turkel) {
    LowMach_Precontioner = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new su2double[nVar];
  }
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  }
  
  else {
    if (rank == MASTER_NODE) cout << "Explicit scheme. No Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Define some auxiliary vectors for computing flow variable
   gradients by least squares, S matrix := inv(R)*traspose(inv(R)),
   c vector := transpose(WA)*(Wb) ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    Cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      Cvector[iVar] = new su2double [nDim];
    
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  CharacPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/
  
  DonorPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      if (rans) {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar+2];
        for (iVar = 0; iVar < nPrimVar + 2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
        for (iVar = 0; iVar < nPrimVar ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
    }
  }
  
  /*--- Store the value of the characteristic primitive variables index at the boundaries ---*/
  
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Delta P at the Actuator Disk ---*/
  
  ActDisk_DeltaP = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaP[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaP[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Delta T at the Actuator Disk ---*/
  
  ActDisk_DeltaT = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaT[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaT[iMarker][iVertex] = 0;
    }
  }

  /*--- Force definition and coefficient arrays for all of the markers ---*/
  
  CPressure = new su2double* [nMarker];
  CPressureTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new su2double [geometry->nVertex[iMarker]];
    CPressureTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Non-dimensional coefficients ---*/
  
  ForceInviscid     = new su2double[nDim];
  MomentInviscid    = new su2double[3];
  CD_Inv         = new su2double[nMarker];
  CL_Inv         = new su2double[nMarker];
  CSF_Inv        = new su2double[nMarker];
  CMx_Inv           = new su2double[nMarker];
  CMy_Inv           = new su2double[nMarker];
  CMz_Inv           = new su2double[nMarker];
  CEff_Inv          = new su2double[nMarker];
  CFx_Inv           = new su2double[nMarker];
  CFy_Inv           = new su2double[nMarker];
  CFz_Inv           = new su2double[nMarker];
  
  ForceMomentum     = new su2double[nDim];
  MomentMomentum    = new su2double[3];
  CD_Mnt        = new su2double[nMarker];
  CL_Mnt        = new su2double[nMarker];
  CSF_Mnt       = new su2double[nMarker];
  CMx_Mnt          = new su2double[nMarker];
  CMy_Mnt          = new su2double[nMarker];
  CMz_Mnt          = new su2double[nMarker];
  CEff_Mnt         = new su2double[nMarker];
  CFx_Mnt          = new su2double[nMarker];
  CFy_Mnt          = new su2double[nMarker];
  CFz_Mnt          = new su2double[nMarker];
  
  Surface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Inv     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  
  Surface_CL_Mnt     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Mnt     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Mnt= new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Rotorcraft coefficients ---*/
  
  CT_Inv           = new su2double[nMarker];
  CQ_Inv           = new su2double[nMarker];
  CMerit_Inv       = new su2double[nMarker];
  
  CT_Mnt           = new su2double[nMarker];
  CQ_Mnt           = new su2double[nMarker];
  CMerit_Mnt       = new su2double[nMarker];

  /*--- Supersonic coefficients ---*/
  
  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];
  
  /*--- Engine simulation ---*/
  
  Inflow_MassFlow     = new su2double[nMarker];
  Inflow_Pressure     = new su2double[nMarker];
  Inflow_Mach         = new su2double[nMarker];
  Inflow_Area         = new su2double[nMarker];
  
  Exhaust_MassFlow    = new su2double[nMarker];
  Exhaust_Pressure    = new su2double[nMarker];
  Exhaust_Temperature = new su2double[nMarker];
  Exhaust_Area        = new su2double[nMarker];
  
  /*--- Init total coefficients ---*/
  
  Total_CD      = 0.0;    Total_CL           = 0.0;    Total_CSF          = 0.0;
  Total_CMx     = 0.0;    Total_CMy          = 0.0;    Total_CMz          = 0.0;
  Total_CEff    = 0.0;    Total_CEquivArea   = 0.0;    Total_CNearFieldOF = 0.0;
  Total_CFx     = 0.0;    Total_CFy          = 0.0;    Total_CFz          = 0.0;
  Total_CT      = 0.0;    Total_CQ           = 0.0;    Total_CMerit       = 0.0;
  Total_MaxHeat = 0.0;    Total_Heat         = 0.0;    Total_ComboObj     = 0.0;
  Total_CpDiff  = 0.0;    Total_HeatFluxDiff = 0.0;
  Total_NetCThrust = 0.0; Total_NetCThrust_Prev = 0.0; Total_BCThrust_Prev = 0.0;
  Total_Power = 0.0;      AoA_Prev           = 0.0;
  Total_CL_Prev = 0.0;    Total_CD_Prev      = 0.0;
  Total_AeroCD = 0.0;     Total_RadialDistortion   = 0.0;    Total_CircumferentialDistortion   = 0.0;

  /*--- Read farfield conditions ---*/
  
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Mach_Inf        = config->GetMach();
  
  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  
  switch(direct_diff) {
    case NO_DERIVATIVE:
      /*--- Default ---*/
      break;
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
      break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
      break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
      break;
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Already done in postprocessing of config ---*/
      break;
    default:
      break;
  }
  
  
  /*--- Initializate fan face pressure, fan face mach number, and mass flow rate ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;
    
    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;
  }
  
  /*--- Initializate quantities for SlidingMesh Interface ---*/
  
  SlidingState = new su2double** [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    SlidingState[iMarker] = NULL;
      
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker] = new su2double* [geometry->GetnVertex(iMarker)];

      for (iPoint = 0; iPoint < geometry->nVertex[iMarker]; iPoint++) {
        SlidingState[iMarker][iPoint] = new su2double[nPrimVar];
      for (iVar = 0; iVar < nVar; iVar++)
        SlidingState[iMarker][iPoint][iVar] = -1;
      }
    }
    else
      SlidingState[iMarker] = NULL;
  }
 
  /*--- Initializate quantities for the mixing process ---*/
  
  AveragedVelocity = new su2double* [nMarker];
  AveragedNormal   = new su2double* [nMarker];
  AveragedGridVel  = new su2double* [nMarker];
  AveragedFlux     = new su2double* [nMarker];
  TotalFlux        = new su2double* [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AveragedVelocity[iMarker] = new su2double [nDim];
    AveragedNormal[iMarker]   = new su2double [nDim];
    AveragedGridVel[iMarker]  = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      AveragedVelocity[iMarker][iDim] = 0.0;
      AveragedNormal[iMarker][iDim]   = 0.0;
      AveragedGridVel [iMarker][iDim] = 0.0;
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AveragedFlux[iMarker] = new su2double [nVar];
    TotalFlux[iMarker]    = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      AveragedFlux[iMarker][iVar] = 0.0;
      TotalFlux[iMarker][iVar]    = 0.0;
    }
  }
  
  AveragedNormalVelocity     = new su2double[nMarker];
  AveragedTangVelocity       = new su2double[nMarker];
  ExtAveragedNormalVelocity  = new su2double[nMarker];
  ExtAveragedTangVelocity    = new su2double[nMarker];
  MassFlow                   = new su2double[nMarker];
  FlowAngle                  = new su2double[nMarker];
  AveragedEnthalpy           = new su2double[nMarker];
  AveragedPressure           = new su2double[nMarker];
  AveragedTotPressure        = new su2double[nMarker];
  AveragedTotTemperature     = new su2double[nMarker];
  ExtAveragedTotPressure     = new su2double[nMarker];
  ExtAveragedTotTemperature  = new su2double[nMarker];
  AveragedDensity            = new su2double[nMarker];
  ExtAveragedPressure        = new su2double[nMarker];
  ExtAveragedDensity         = new su2double[nMarker];
  AveragedSoundSpeed         = new su2double[nMarker];
  AveragedEntropy            = new su2double[nMarker];
  AveragedTangGridVelocity   = new su2double[nMarker];
  AveragedMach               = new su2double[nMarker];
  AveragedNormalMach         = new su2double[nMarker];
  AveragedTangMach           = new su2double[nMarker];
  
  /*--- Initializate quantities for turboperformace ---*/
  
  TotalStaticEfficiency = new su2double[nMarkerTurboPerf];
  TotalTotalEfficiency  = new su2double[nMarkerTurboPerf];
  KineticEnergyLoss     = new su2double[nMarkerTurboPerf];
  TotalPressureLoss     = new su2double[nMarkerTurboPerf];
  MassFlowIn            = new su2double[nMarkerTurboPerf];
  MassFlowOut           = new su2double[nMarkerTurboPerf];
  FlowAngleIn           = new su2double[nMarkerTurboPerf];
  FlowAngleOut          = new su2double[nMarkerTurboPerf];
  EulerianWork          = new su2double[nMarkerTurboPerf];
  TotalEnthalpyIn       = new su2double[nMarkerTurboPerf];
  PressureRatio         = new su2double[nMarkerTurboPerf];
  PressureOut           = new su2double[nMarkerTurboPerf];
  EnthalpyOut           = new su2double[nMarkerTurboPerf];
  MachIn                = new su2double[nMarkerTurboPerf];
  MachOut               = new su2double[nMarkerTurboPerf];
  NormalMachIn          = new su2double[nMarkerTurboPerf];
  NormalMachOut         = new su2double[nMarkerTurboPerf];
  VelocityOutIs         = new su2double[nMarkerTurboPerf];
  
  for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++) {
    TotalStaticEfficiency[iMarker] = 0.0;
    TotalTotalEfficiency[iMarker]  = 0.0;
    KineticEnergyLoss[iMarker]     = 0.0;
    TotalPressureLoss[iMarker]     = 0.0;
    MassFlowIn[iMarker]            = 0.0;
    MassFlowOut[iMarker]           = 0.0;
    FlowAngleIn[iMarker]           = 0.0;
    FlowAngleOut[iMarker]          = 0.0;
    EulerianWork[iMarker]          = 0.0;
    TotalEnthalpyIn[iMarker]       = 0.0;
    PressureRatio[iMarker]         = 0.0;
    PressureOut[iMarker]           = 0.0;
    EnthalpyOut[iMarker]           = 0.0;
    MachIn[iMarker]                = 0.0;
    MachOut[iMarker]               = 0.0;
    NormalMachIn[iMarker]          = 0.0;
    NormalMachOut[iMarker]         = 0.0;
    VelocityOutIs[iMarker]         = 0.0;
  }
  
  
  /*--- Initialize the cauchy critera array for fixed CL mode ---*/
  
  if (config->GetFixed_CL_Mode())
    
    Cauchy_Serie = new su2double [config->GetCauchy_Elems()+1];
  
  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  
  if (!restart || (iMesh != MESH_0)) {
    
    /*--- Restart the solution from the free-stream state ---*/
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
    
  } else {
        
    /*--- Multizone problems require the number of the zone to be appended. ---*/
    
    if (nZone > 1) filename = config->GetMultizone_FileName(filename, iZone);
    
    /*--- Modify file name for a dual-time unsteady restart ---*/
    
    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }
    
    /*--- Modify file name for a time stepping unsteady restart ---*/
    
    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }

        /*--- Open the restart file, throw an error if this fails. ---*/
    
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    
    map<unsigned long,unsigned long> Global2Local;
    map<unsigned long,unsigned long>::const_iterator MI;
    
    /*--- Now fill array with the transform values only for local points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    
    /*--- Read all lines in the restart file ---*/
    
    long iPoint_Local;
    unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
    
    /*--- The first line is the header ---*/
    
    getline (restart_file, text_line);
    
    /*--- Solution ---*/
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
      istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/
      
      MI = Global2Local.find(iPoint_Global);
      if (MI != Global2Local.end()) {
      
      iPoint_Local = Global2Local[iPoint_Global];
      
      /*--- Load the solution for this node. Note that the first entry
       on the restart file line is the global index, followed by the
       node coordinates, and then the conservative variables. ---*/
      
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];

        node[iPoint_Local] = new CEulerVariable(Solution, nDim, nVar, config);
        iPoint_Global_Local++;
      }

    }
    
    /*--- Detect a wrong solution file ---*/
    
    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
    
#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
        cout << "It could be empty lines at the end of the file." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      node[iPoint] = new CEulerVariable(Solution, nDim, nVar, config);
    
    /*--- Close the restart file ---*/
    
    restart_file.close();
    
  }
  
  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
    
  counter_local = 0;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    Density = node[iPoint]->GetSolution(0);

    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);

    StaticEnergy= node[iPoint]->GetSolution(nDim+1)/Density - 0.5*Velocity2;

    FluidModel->SetTDState_rhoe(Density, StaticEnergy);
    Pressure= FluidModel->GetPressure();
    Temperature= FluidModel->GetTemperature();

    /*--- Use the values at the infinity ---*/

    if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);
      counter_local++;
    }

  }

  /*--- Warning message about non-physical points ---*/

  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }


  
  /*--- Define solver parameters needed for execution of destructor ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
}

CEulerSolver::~CEulerSolver(void) {

  unsigned short iVar, iMarker;
  unsigned long iVertex;

  /*--- Array deallocation ---*/

  if (CD_Inv != NULL)         delete [] CD_Inv;
  if (CL_Inv != NULL)         delete [] CL_Inv;
  if (CSF_Inv != NULL)    delete [] CSF_Inv;
  if (CMx_Inv != NULL)           delete [] CMx_Inv;
  if (CMy_Inv != NULL)           delete [] CMy_Inv;
  if (CMz_Inv != NULL)           delete [] CMz_Inv;
  if (CFx_Inv != NULL)           delete [] CFx_Inv;
  if (CFy_Inv != NULL)           delete [] CFy_Inv;
  if (CFz_Inv != NULL)           delete [] CFz_Inv;
  if (Surface_CL_Inv != NULL) delete[] Surface_CL_Inv;
  if (Surface_CD_Inv != NULL) delete[] Surface_CD_Inv;
  if (Surface_CSF_Inv != NULL) delete[] Surface_CSF_Inv;
  if (Surface_CEff_Inv != NULL) delete[] Surface_CEff_Inv;
  if (Surface_CFx_Inv != NULL)  delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv != NULL)  delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv != NULL)  delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;
  
  if (CD_Mnt != NULL)         delete [] CD_Mnt;
  if (CL_Mnt != NULL)         delete [] CL_Mnt;
  if (CSF_Mnt != NULL)    delete [] CSF_Mnt;
  if (CMx_Mnt != NULL)           delete [] CMx_Mnt;
  if (CMy_Mnt != NULL)           delete [] CMy_Mnt;
  if (CMz_Mnt != NULL)           delete [] CMz_Mnt;
  if (CFx_Mnt != NULL)           delete [] CFx_Mnt;
  if (CFy_Mnt != NULL)           delete [] CFy_Mnt;
  if (CFz_Mnt != NULL)           delete [] CFz_Mnt;
  if (Surface_CL_Mnt != NULL) delete[] Surface_CL_Mnt;
  if (Surface_CD_Mnt != NULL) delete[] Surface_CD_Mnt;
  if (Surface_CSF_Mnt != NULL) delete[] Surface_CSF_Mnt;
  if (Surface_CEff_Mnt != NULL) delete[] Surface_CEff_Mnt;
  if (Surface_CFx_Mnt != NULL)  delete [] Surface_CFx_Mnt;
  if (Surface_CFy_Mnt != NULL)  delete [] Surface_CFy_Mnt;
  if (Surface_CFz_Mnt != NULL)  delete [] Surface_CFz_Mnt;
  if (Surface_CMx_Mnt != NULL)  delete [] Surface_CMx_Mnt;
  if (Surface_CMy_Mnt != NULL)  delete [] Surface_CMy_Mnt;
  if (Surface_CMz_Mnt != NULL)  delete [] Surface_CMz_Mnt;

  if (Surface_CL != NULL)    delete [] Surface_CL;
  if (Surface_CD != NULL)    delete [] Surface_CD;
  if (Surface_CSF != NULL) delete [] Surface_CSF;
  if (Surface_CEff != NULL) delete [] Surface_CEff;
  if (Surface_CFx != NULL)      delete [] Surface_CFx;
  if (Surface_CFy != NULL)      delete [] Surface_CFy;
  if (Surface_CFz != NULL)      delete [] Surface_CFz;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)          delete [] CEff_Inv;
  if (CMerit_Inv != NULL)        delete [] CMerit_Inv;
  if (CT_Inv != NULL)            delete [] CT_Inv;
  if (CQ_Inv != NULL)            delete [] CQ_Inv;
  if (CEquivArea_Inv != NULL)    delete [] CEquivArea_Inv;
  if (CNearFieldOF_Inv != NULL)  delete [] CNearFieldOF_Inv;
  
  if (CEff_Mnt != NULL)          delete [] CEff_Mnt;
  if (CMerit_Mnt != NULL)        delete [] CMerit_Mnt;
  if (CT_Mnt != NULL)            delete [] CT_Mnt;
  if (CQ_Mnt != NULL)            delete [] CQ_Mnt;

  if (ForceInviscid != NULL)     delete [] ForceInviscid;
  if (MomentInviscid != NULL)    delete [] MomentInviscid;
  if (ForceMomentum != NULL)     delete [] ForceMomentum;
  if (MomentMomentum != NULL)    delete [] MomentMomentum;
  if (Inflow_MassFlow != NULL)  delete [] Inflow_MassFlow;
  if (Exhaust_MassFlow != NULL)  delete [] Exhaust_MassFlow;
  if (Exhaust_Area != NULL)      delete [] Exhaust_Area;
  if (Inflow_Pressure != NULL)  delete [] Inflow_Pressure;
  if (Inflow_Mach != NULL)      delete [] Inflow_Mach;
  if (Inflow_Area != NULL)      delete [] Inflow_Area;

  if (Exhaust_Pressure != NULL)  delete [] Exhaust_Pressure;
  if (Exhaust_Temperature != NULL)      delete [] Exhaust_Temperature;

  if (iPoint_UndLapl != NULL)       delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL)       delete [] jPoint_UndLapl;

  if (Primitive != NULL)        delete [] Primitive;
  if (Primitive_i != NULL)      delete [] Primitive_i;
  if (Primitive_j != NULL)      delete [] Primitive_j;

  if (Secondary != NULL)        delete [] Secondary;
  if (Secondary_i != NULL)      delete [] Secondary_i;
  if (Secondary_j != NULL)      delete [] Secondary_j;

  if (LowMach_Precontioner != NULL) {
    for (iVar = 0; iVar < nVar; iVar ++)
      delete [] LowMach_Precontioner[iVar];
    delete [] LowMach_Precontioner;
  }
  
  if (CPressure != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CPressure[iMarker];
    delete [] CPressure;
  }
  
  if (CPressureTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CPressureTarget[iMarker];
    delete [] CPressureTarget;
  }
  
  if (CharacPrimVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] CharacPrimVar[iMarker][iVertex];
      delete [] CharacPrimVar[iMarker];
    }
    delete [] CharacPrimVar;
  }
  
  if ( SlidingState != NULL ) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          delete [] SlidingState[iMarker][iVertex];
          
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }
  
  if (DonorPrimVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] DonorPrimVar[iMarker][iVertex];
      delete [] DonorPrimVar[iMarker];
    }
    delete [] DonorPrimVar;
  }
  
  if (DonorGlobalIndex != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] DonorGlobalIndex[iMarker];
    delete [] DonorGlobalIndex;
  }

  if (ActDisk_DeltaP != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ActDisk_DeltaP[iMarker];
    delete [] ActDisk_DeltaP;
  }

  if (ActDisk_DeltaT != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] ActDisk_DeltaT[iMarker];
    delete [] ActDisk_DeltaT;
  }

  if (nVertex != NULL)  delete [] nVertex;

  if (HeatFlux != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] HeatFlux[iMarker];
    }
    delete [] HeatFlux;
  }
  
  if (HeatFluxTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] HeatFluxTarget[iMarker];
    }
    delete [] HeatFluxTarget;
  }
  
  if (YPlus != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] YPlus[iMarker];
    }
    delete [] YPlus;
  }
  
  if (Cauchy_Serie != NULL)  delete [] Cauchy_Serie;
  
  if (FluidModel != NULL) delete FluidModel;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (AveragedVelocity[iMarker] != NULL) delete [] AveragedVelocity[iMarker];
    if (AveragedNormal[iMarker]   != NULL) delete [] AveragedNormal[iMarker];
    if (AveragedGridVel[iMarker]  != NULL) delete [] AveragedGridVel[iMarker];
    if (AveragedFlux[iMarker]     != NULL) delete [] AveragedFlux[iMarker];
    if (TotalFlux[iMarker]        != NULL) delete [] TotalFlux[iMarker];
  }
  if (AveragedVelocity != NULL) delete [] AveragedVelocity;
  if (AveragedNormal   != NULL) delete [] AveragedNormal;
  if (AveragedGridVel  != NULL) delete [] AveragedGridVel;
  if (AveragedFlux     != NULL) delete [] AveragedFlux;
  if (TotalFlux        != NULL) delete [] TotalFlux;

  if (AveragedNormalVelocity    != NULL) delete [] AveragedNormalVelocity;
  if (AveragedTangVelocity      != NULL) delete [] AveragedTangVelocity;
  if (ExtAveragedNormalVelocity != NULL) delete [] ExtAveragedNormalVelocity;
  if (ExtAveragedTangVelocity   != NULL) delete [] ExtAveragedTangVelocity;
  if (MassFlow                  != NULL) delete [] MassFlow;
  if (FlowAngle                 != NULL) delete [] FlowAngle;
  if (AveragedEnthalpy          != NULL) delete [] AveragedEnthalpy;
  if (AveragedPressure          != NULL) delete [] AveragedPressure;
  if (AveragedTotPressure       != NULL) delete [] AveragedTotPressure;
  if (AveragedTotTemperature    != NULL) delete [] AveragedTotTemperature;
  if (ExtAveragedTotPressure    != NULL) delete [] ExtAveragedTotPressure;
  if (ExtAveragedTotTemperature != NULL) delete [] ExtAveragedTotTemperature;
  if (AveragedDensity           != NULL) delete [] AveragedDensity;
  if (ExtAveragedPressure       != NULL) delete [] ExtAveragedPressure;
  if (ExtAveragedDensity        != NULL) delete [] ExtAveragedDensity;
  if (AveragedSoundSpeed        != NULL) delete [] AveragedSoundSpeed;
  if (AveragedEntropy           != NULL) delete [] AveragedEntropy;
  if (AveragedTangGridVelocity  != NULL) delete [] AveragedTangGridVelocity;
  if (AveragedMach              != NULL) delete [] AveragedMach;
  if (AveragedNormalMach        != NULL) delete [] AveragedNormalMach;
  if (AveragedTangMach          != NULL) delete [] AveragedTangMach;

  if (TotalStaticEfficiency != NULL) delete [] TotalStaticEfficiency;
  if (TotalTotalEfficiency  != NULL) delete [] TotalTotalEfficiency;
  if (KineticEnergyLoss     != NULL) delete [] KineticEnergyLoss;
  if (TotalPressureLoss     != NULL) delete [] TotalPressureLoss;
  if (MassFlowIn            != NULL) delete [] MassFlowIn;
  if (MassFlowOut           != NULL) delete [] MassFlowOut;
  if (FlowAngleIn           != NULL) delete [] FlowAngleIn;
  if (FlowAngleOut          != NULL) delete [] FlowAngleOut;
  if (EulerianWork          != NULL) delete [] EulerianWork;
  if (TotalEnthalpyIn       != NULL) delete [] TotalEnthalpyIn;
  if (PressureRatio         != NULL) delete [] PressureRatio;
  if (PressureOut           != NULL) delete [] PressureOut;
  if (EnthalpyOut           != NULL) delete [] EnthalpyOut;
  if (MachIn                != NULL) delete [] MachIn;
  if (MachOut               != NULL) delete [] MachOut;
  if (NormalMachIn          != NULL) delete [] NormalMachIn;
  if (NormalMachOut         != NULL) delete [] NormalMachOut;
  if (VelocityOutIs         != NULL) delete [] VelocityOutIs;
  
}

void CEulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CEulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new su2double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex] = Buffer_Send_Undivided_Laplacian[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Undivided_Laplacian;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetUndivided_Laplacian(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Undivided_Laplacian;
      
    }
    
  }
  
}

void CEulerSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR, *Buffer_Receive_Neighbor = NULL, *Buffer_Send_Neighbor = NULL;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];
      Buffer_Receive_Neighbor = new unsigned short [nBufferR_Vector];
      Buffer_Send_Neighbor = new unsigned short[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetLambda();
        Buffer_Send_Neighbor[iVertex] = geometry->node[iPoint]->GetnPoint();
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_Neighbor, nBufferS_Vector, MPI_UNSIGNED_SHORT, send_to, 1,
                        Buffer_Receive_Neighbor, nBufferR_Vector, MPI_UNSIGNED_SHORT, receive_from, 1, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
        Buffer_Receive_Neighbor[iVertex] = Buffer_Send_Neighbor[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      delete [] Buffer_Send_Neighbor;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetLambda(Buffer_Receive_Lambda[iVertex]);
        geometry->node[iPoint]->SetnNeighbor(Buffer_Receive_Neighbor[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      delete [] Buffer_Receive_Neighbor;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetSensor();
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetSensor(Buffer_Receive_Lambda[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CEulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nVar];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Limiter[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}

void CEulerSolver::Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad*nDim;        nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Primitive(iVar, iDim);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient_Primitive(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
  }
  
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CEulerSolver::Set_MPI_Primitive_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nPrimVarGrad];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad;        nBufferR_Vector = nVertexR*nPrimVarGrad;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter_Primitive(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          node[iPoint]->SetLimiter_Primitive(iVar, Limiter[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}

void CEulerSolver::Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0;
  long iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  //bool ActDisk_Perimeter;
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));
  
  unsigned short nPrimVar_ = nPrimVar;
  if (rans) nPrimVar_ += 2; // Add two extra variables for the turbulence.
  
#ifdef HAVE_MPI
  
  /*--- MPI initialization ---*/
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  MPI_Status status;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double *Buffer_Send_PrimVar = NULL;
  long      *Buffer_Send_Data    = NULL;
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  su2double *iPrimVar = new su2double [nPrimVar_];
  
  unsigned long Buffer_Size_PrimVar = 0;
  unsigned long Buffer_Size_Data    = 0;
  
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double *Buffer_Receive_PrimVar = NULL;
  long      *Buffer_Receive_Data    = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          //ActDisk_Perimeter = geometry->vertex[iMarker][iVertex]->GetActDisk_Perimeter();
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
//          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain()) && (!ActDisk_Perimeter)) {
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar += nPointTotal_s[iDomain]*(nPrimVar_);
    Buffer_Size_Data += nPointTotal_s[iDomain]*(3);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar = new su2double[Buffer_Size_PrimVar];
  Buffer_Send_Data    = new long[Buffer_Size_Data];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          //ActDisk_Perimeter = geometry->vertex[iMarker][iVertex]->GetActDisk_Perimeter();
          
//          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain()) && (!ActDisk_Perimeter)) {
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetPrimitive(iVar);
            }
            if (rans) {
              Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+nPrimVar] = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
              Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = 0.0;
            }
            
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            
            Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(0)]  = iGlobalIndex;
            Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(1)] = jVertex;
            Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(2)]  = jMarker;
            
            iPointTotal++;
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar_)],
                     nPointTotal_s[iDomain]*(nPrimVar_), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
      SU2_MPI::Bsend(&Buffer_Send_Data[PointTotal_Counter*(3)],
                     nPointTotal_s[iDomain]*(3), MPI_LONG, iDomain,
                     iDomain+nDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar_)];
      Buffer_Receive_Data               = new long[nPointTotal_s[iDomain]*(3)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar_); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar_)+iter];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(3); iter++)
        Buffer_Receive_Data[iter] = Buffer_Send_Data[PointTotal_Counter*(3)+iter];
      
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        for (iVar = 0; iVar < nPrimVar_; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar_)+iVar];
        
        iGlobal       =  Buffer_Receive_Data[iPoint*(3)+(0)];
        iVertex      =  Buffer_Receive_Data[iPoint*(3)+(1)];
        iMarker      = Buffer_Receive_Data[iPoint*(3)+(2)];
        
        for (iVar = 0; iVar < nPrimVar_; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      delete [] Buffer_Receive_Data;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar_)];
      Buffer_Receive_Data               = new long [nPointTotal_r[iDomain]*(3)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar_) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status);
      
      SU2_MPI::Recv(Buffer_Receive_Data, nPointTotal_r[iDomain]*(3) , MPI_LONG,
                    iDomain, rank+nDomain, MPI_COMM_WORLD, &status);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = Buffer_Receive_Data[iPoint*(3)+(0)];
        iVertex      = Buffer_Receive_Data[iPoint*(3)+(1)];
        iMarker      = Buffer_Receive_Data[iPoint*(3)+(2)];
        
        for (iVar = 0; iVar < nPrimVar_; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar_)+iVar];
        
        for (iVar = 0; iVar < nPrimVar_; iVar++) {
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        }
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      delete [] Buffer_Receive_Data;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  delete[] Buffer_Send_Data;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
}

void CEulerSolver::Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0;
  long iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  
  /*--- MPI initialization ---*/
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  MPI_Status status, status_;
  

#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  su2double        *iPrimVar          = new su2double [nPrimVar];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(nPrimVar+3);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetPrimitive(iVar);
            }
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+0)]  = su2double(iGlobalIndex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = su2double(jVertex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+2)]  = su2double(jMarker);
            
            iPointTotal++;
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)],
                     nPointTotal_s[iDomain]*(nPrimVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar+3)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar+3); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal       =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar+3)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar+3) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar,  iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
}

void CEulerSolver::Set_MPI_Interface(CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker;
  long nDomain = 0, iDomain, jDomain;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  
  /*--- MPI initialization ---*/
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  MPI_Status status, status_;
  

#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [nPrimVar];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            Buffer_Send_nPointTotal++;
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(nPrimVar+3);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
          
          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
            
            iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
            jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
            jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetPrimitive(iVar);
            }
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+0)]  = su2double(iGlobalIndex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = su2double(jVertex);
            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+2)]  = su2double(jMarker);
            
            iPointTotal++;
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)],
                     nPointTotal_s[iDomain]*(nPrimVar+3), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar+3)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar+3); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal       =  SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar+3)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar+3) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
}

//void CEulerSolver::Set_MPI_Secondary_Gradient(CGeometry *geometry, CConfig *config) {
//  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
//  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
//  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
//  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
//  int send_to, receive_from;
//
//  su2double **Gradient = new su2double* [nSecondaryVarGrad];
//  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//    Gradient[iVar] = new su2double[nDim];
//
//#ifdef HAVE_MPI
//  MPI_Status status;
//#endif
//
//  for (iMarker = 0; iMarker < nMarker; iMarker++) {
//
//    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
//        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
//
//      MarkerS = iMarker;  MarkerR = iMarker+1;
//
//      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
//      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
//
//      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
//      nBufferS_Vector = nVertexS*nSecondaryVarGrad*nDim;        nBufferR_Vector = nVertexR*nSecondaryVarGrad*nDim;
//
//      /*--- Allocate Receive and send buffers  ---*/
//      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
//      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
//
//      /*--- Copy the solution old that should be sended ---*/
//      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
//        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Buffer_Send_Gradient[iDim*nSecondaryVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Secondary(iVar, iDim);
//      }
//
//#ifdef HAVE_MPI
//
//      /*--- Send/Receive information using Sendrecv ---*/
//      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
//                   Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
//
//#else
//
//      /*--- Receive information without MPI ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Buffer_Receive_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//      }
//
//#endif
//
//      /*--- Deallocate send buffer ---*/
//      delete [] Buffer_Send_Gradient;
//
//      /*--- Do the coordinate transformation ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//
//        /*--- Find point and its type of transformation ---*/
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
//
//        /*--- Retrieve the supplied periodic information. ---*/
//        angles = config->GetPeriodicRotation(iPeriodic_Index);
//
//        /*--- Store angles separately for clarity. ---*/
//        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
//        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
//        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
//
//        /*--- Compute the rotation matrix. Note that the implicit
//         ordering is rotation about the x-axis, y-axis,
//         then z-axis. Note that this is the transpose of the matrix
//         used during the preprocessing stage. ---*/
//        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
//        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
//        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
//
//        /*--- Copy conserved variables before performing transformation. ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//
//        /*--- Need to rotate the gradients for all conserved variables. ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//          if (nDim == 2) {
//            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//          }
//          else {
//            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//          }
//        }
//
//        /*--- Store the received information ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            node[iPoint]->SetGradient_Secondary(iVar, iDim, Gradient[iVar][iDim]);
//
//      }
//
//      /*--- Deallocate receive buffer ---*/
//      delete [] Buffer_Receive_Gradient;
//
//    }
//
//  }
//
//  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//    delete [] Gradient[iVar];
//  delete [] Gradient;
//
//}

//void CEulerSolver::Set_MPI_Secondary_Limiter(CGeometry *geometry, CConfig *config) {
//  unsigned short iVar, iMarker, MarkerS, MarkerR;
//  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
//  su2double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
//  int send_to, receive_from;
//
//  su2double *Limiter = new su2double [nSecondaryVarGrad];
//
//#ifdef HAVE_MPI
//  MPI_Status status;
//#endif
//
//  for (iMarker = 0; iMarker < nMarker; iMarker++) {
//
//    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
//        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
//
//      MarkerS = iMarker;  MarkerR = iMarker+1;
//
//      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
//      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
//
//      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
//      nBufferS_Vector = nVertexS*nSecondaryVarGrad;        nBufferR_Vector = nVertexR*nSecondaryVarGrad;
//
//      /*--- Allocate Receive and send buffers  ---*/
//      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
//      Buffer_Send_Limit = new su2double[nBufferS_Vector];
//
//      /*--- Copy the solution old that should be sended ---*/
//      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
//        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter_Secondary(iVar);
//      }
//
//#ifdef HAVE_MPI
//
//      /*--- Send/Receive information using Sendrecv ---*/
//      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
//                   Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
//
//#else
//
//      /*--- Receive information without MPI ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
//      }
//
//#endif
//
//      /*--- Deallocate send buffer ---*/
//      delete [] Buffer_Send_Limit;
//
//      /*--- Do the coordinate transformation ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//
//        /*--- Find point and its type of transformation ---*/
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//
//        /*--- Copy conserved variables before performing transformation. ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
//
//        /*--- Copy transformed conserved variables back into buffer. ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          node[iPoint]->SetLimiter_Secondary(iVar, Limiter[iVar]);
//
//      }
//
//      /*--- Deallocate receive buffer ---*/
//      delete [] Buffer_Receive_Limit;
//
//    }
//
//  }
//
//  delete [] Limiter;
//
//}

void CEulerSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) {
  
  su2double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0,
  Energy_FreeStream = 0.0, ModVel_FreeStreamND = 0.0, Velocity_Reynolds = 0.0,
  Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
  Temperature_Ref = 0.0, Time_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Energy_Ref= 0.0,
  Froude = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0, TgammaR = 0.0;
  su2double Mu_RefND, Mu_Temperature_RefND, Mu_SND;

  unsigned short iDim;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Local variables ---*/
  
  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double Reynolds         = config->GetReynolds();
  bool unsteady           = (config->GetUnsteady_Simulation() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool gravity            = config->GetGravityForce();
  bool turbulent          = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
  bool tkeNeeded          = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool standard_air       = (config->GetKind_FluidModel() == STANDARD_AIR);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);
  bool aeroelastic        = config->GetAeroelastic_Simulation();
  
  /*--- Set temperature via the flutter speed index ---*/
  if (aeroelastic) {
    su2double vf             = config->GetAeroelastic_Flutter_Speed_Index();
    su2double w_alpha        = config->GetAeroelastic_Frequency_Pitch();
    su2double b              = config->GetLength_Reynolds()/2.0; // airfoil semichord, Reynolds length is by defaul 1.0
    su2double mu             = config->GetAeroelastic_Airfoil_Mass_Ratio();
    // The temperature times gamma times the gas constant. Depending on the FluidModel temp is calculated below.
    TgammaR = ((vf*vf)*(b*b)*(w_alpha*w_alpha)*mu) / (Mach*Mach);
  }
  
  /*--- Compressible non dimensionalization ---*/

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);

      FluidModel = new CIdealGas(1.4, config->GetGas_Constant());
      if (free_stream_temp) {
        if (aeroelastic) {
          Temperature_FreeStream = TgammaR / (config->GetGas_Constant()*1.4);
          config->SetTemperature_FreeStream(Temperature_FreeStream);
        }
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case IDEAL_GAS:

      FluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case VW_GAS:

      FluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                       config->GetPressure_Critical(), config->GetTemperature_Critical());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case PR_GAS:

      FluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                     config->GetTemperature_Critical(), config->GetAcentric_Factor());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

  }

  Mach2Vel_FreeStream = FluidModel->GetSoundSpeed();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }

  /*--- Compute the modulus of the free stream velocity ---*/

  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Viscous initialization ---*/

  if (viscous) {

    /*--- Reynolds based initialization ---*/

    if (reynolds_init) {

      /*--- First, check if there is mesh motion. If yes, use the Mach
         number relative to the body to initialize the flow. ---*/

      if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
      else Velocity_Reynolds = ModVel_FreeStream;

      /*--- Change of measurement system, hard coded value working only with STANDAR AIR model ---*/

      if (standard_air) {
        Mu_RefND = 1.716E-5;
        Mu_Temperature_RefND = 273.15; // 273.0
        Mu_SND = 110.4;                // 111.0
        if (config->GetSystemMeasurements() == SI) {
          config->SetMu_RefND(Mu_RefND);
          config->SetMu_Temperature_RefND(Mu_Temperature_RefND);
          config->SetMu_SND(Mu_SND);
        }
        if (config->GetSystemMeasurements() == US) {
          Mu_RefND = Mu_RefND/47.88025898;
          Mu_Temperature_RefND = (Mu_Temperature_RefND - 273.15) * 1.8 + 491.67;
          Mu_SND = (Mu_SND - 273.15) * 1.8 + 491.67;
          config->SetMu_RefND(Mu_RefND);
          config->SetMu_Temperature_RefND(Mu_Temperature_RefND);
          config->SetMu_SND(Mu_SND);
        }
      }

      /*--- For viscous flows, pressure will be computed from a density
         that is found from the Reynolds number. The viscosity is computed
         from the dimensional version of Sutherland's law ---*/

      FluidModel->SetLaminarViscosityModel(config);

      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);

      Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
      config->SetDensity_FreeStream(Density_FreeStream);
      FluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
      Pressure_FreeStream = FluidModel->GetPressure();
      config->SetPressure_FreeStream(Pressure_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Thermodynamics quantities based initialization ---*/

    else {

      FluidModel->SetLaminarViscosityModel(config);
      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Turbulence kinetic energy ---*/

    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }
  else {

    /*--- For inviscid flow, energy is calculated from the specified
       FreeStream quantities using the proper gas law. ---*/

    Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

  }

  /*-- Compute the freestream energy. ---*/

  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute non dimensional quantities. By definition,
     Lref is one because we have converted the grid to meters. ---*/

  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
    Density_Ref       = 1.0;
    Temperature_Ref   = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
    Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;                        config->SetForce_Ref(Force_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);         config->SetFroude(Froude);
  
  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  
  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }
  
  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  
  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);
  
  
  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);
  
  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);
  
  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);
  
  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);
  
  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/(Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStream(Omega_FreeStream);
  
  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/(Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);
  
  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/
 
  /*--- Delete the original (dimensional) FluidModel object before replacing. ---*/
  
  delete FluidModel;
 
  switch (config->GetKind_FluidModel()) {
      
    case STANDARD_AIR:
      FluidModel = new CIdealGas(1.4, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case IDEAL_GAS:
      FluidModel = new CIdealGas(Gamma, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case VW_GAS:
      FluidModel = new CVanDerWaalsGas(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                       config->GetTemperature_Critical()/config->GetTemperature_Ref());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case PR_GAS:
      FluidModel = new CPengRobinson(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                     config->GetTemperature_Critical()/config->GetTemperature_Ref(), config->GetAcentric_Factor());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
  }
  
  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
  
  if (viscous) {
    
    /*--- Constant viscosity model ---*/
    config->SetMu_ConstantND(config->GetMu_ConstantND()/Viscosity_Ref);
    
    /*--- Sutherland's model ---*/
    
    config->SetMu_RefND(config->GetMu_RefND()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_SND()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_RefND()/config->GetTemperature_Ref());
    
    /* constant thermal conductivity model */
    config->SetKt_ConstantND(config->GetKt_ConstantND()/Conductivity_Ref);
    
    FluidModel->SetLaminarViscosityModel(config);
    FluidModel->SetThermalConductivityModel(config);
    
  }
  
  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);
  
  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);
  
  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);
  
  /*--- Write output to the console if this is the master node and first domain ---*/
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {
    
    cout.precision(6);
    
    if (viscous) {
      cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
      cout << "based on the free-stream temperature and a density computed" << endl;
      cout << "from the Reynolds number." << endl;
    } else {
      cout << "Inviscid flow: Computing density based on free-stream" << endl;
      cout << "temperature and pressure using the ideal gas law." << endl;
    }
    
    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;
    
    cout <<"-- Input conditions:"<< endl;
    
    switch (config->GetKind_FluidModel()) {

      case STANDARD_AIR:
        cout << "Fluid Model: STANDARD_AIR "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant();
        if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;

      case IDEAL_GAS:
        cout << "Fluid Model: IDEAL_GAS "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        break;

      case VW_GAS:
        cout << "Fluid Model: Van der Waals "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;

      case PR_GAS:
        cout << "Fluid Model: Peng-Robinson "<< endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
        cout << "Specific Heat Ratio: "<< Gamma << endl;
        cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
        cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
        cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
        cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
        break;

    }
    if (viscous) {
      switch (config->GetKind_ViscosityModel()) {

        case CONSTANT_VISCOSITY:
          cout << "Viscosity Model: CONSTANT_VISCOSITY  "<< endl;
          cout << "Laminar Viscosity: " << config->GetMu_ConstantND()*Viscosity_Ref;
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          break;

        case SUTHERLAND:
          cout << "Viscosity Model: SUTHERLAND "<< endl;
          cout << "Ref. Laminar Viscosity: " << config->GetMu_RefND()*Viscosity_Ref;
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Ref. Temperature: " << config->GetMu_Temperature_RefND()*config->GetTemperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Sutherland Constant: "<< config->GetMu_SND()*config->GetTemperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          cout << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND()<< endl;
          cout << "Sutherland constant (non-dim): "<< config->GetMu_SND()<< endl;
          break;

      }
      switch (config->GetKind_ConductivityModel()) {

        case CONSTANT_PRANDTL:
          cout << "Conductivity Model: CONSTANT_PRANDTL  "<< endl;
          cout << "Prandtl: " << config->GetPrandtl_Lam()<< endl;
          break;

        case CONSTANT_CONDUCTIVITY:
          cout << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< endl;
          cout << "Molecular Conductivity: " << config->GetKt_ConstantND()*Conductivity_Ref<< " W/m^2.K." << endl;
          cout << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< endl;
          break;

      }
    }

    
    cout << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Free-stream total pressure: " << config->GetPressure_FreeStream() * pow( 1.0+Mach*Mach*0.5*(Gamma-1.0), Gamma/(Gamma-1.0) );
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Free-stream temperature: " << config->GetTemperature_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    
    cout << "Free-stream density: " << config->GetDensity_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    if (nDim == 2) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
    }
    if (config->GetSystemMeasurements() == SI) cout << " m/s. ";
    else if (config->GetSystemMeasurements() == US) cout << " ft/s. ";
    
    cout << "Magnitude: "   << config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m/s (" << config->GetModVel_FreeStream()*1.94384 << " KTS)." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s (" << config->GetModVel_FreeStream()*0.592484 << " KTS)." << endl;
    
    cout << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
    
    if (viscous) {
      cout << "Free-stream viscosity: " << config->GetViscosity_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
        cout << "Free-stream specific dissipation: " << config->GetOmega_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " 1/s." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " 1/s." << endl;
      }
    }
    
    if (unsteady) { cout << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << endl; }
    
    /*--- Print out reference values. ---*/
    
    cout <<"-- Reference values:"<< endl;
    
    cout << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;

    cout << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Reference temperature: " << config->GetTemperature_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    
    cout << "Reference density: " << config->GetDensity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    cout << "Reference velocity: " << config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;
    
    cout << "Reference energy per unit mass: " << config->GetEnergy_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
    
    if (viscous) {
      cout << "Reference viscosity: " << config->GetViscosity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      cout << "Reference conductivity: " << config->GetConductivity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " W/m^2.K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf/ft.s.R." << endl;
    }
    
    
    if (unsteady) cout << "Reference time: " << config->GetTime_Ref() <<" s." << endl;
    
    /*--- Print out resulting non-dim values here. ---*/
    
    cout << "-- Resulting non-dimensional state:" << endl;
    cout << "Mach number (non-dim): " << config->GetMach() << endl;
    if (viscous) {
      cout << "Reynolds number (non-dim): " << config->GetReynolds() <<". Re length: " << config->GetLength_Reynolds();
      if (config->GetSystemMeasurements() == SI) cout << " m." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " ft." << endl;
    }
    if (gravity) {
      cout << "Froude number (non-dim): " << Froude << endl;
      cout << "Lenght of the baseline wave (non-dim): " << 2.0*PI_NUMBER*Froude*Froude << endl;
    }
    
    cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
    cout << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << endl;

    cout << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << endl;
    
    cout << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << endl;
    
    if (nDim == 2) {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    cout << "Magnitude: "    << config->GetModVel_FreeStreamND() << endl;
    
    cout << "Free-stream total energy per unit mass (non-dim): " << config->GetEnergy_FreeStreamND() << endl;
    
    if (viscous) {
      cout << "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << endl;
        cout << "Free-stream specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << endl;
      }
    }
    
    if (unsteady) {
      cout << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << endl;
      cout << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << endl;
    }
    
    cout << endl;
    
  }
  
}

void CEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar, iDim;
  su2double Area_Children, Area_Parent,
  *Solution_Fine, *Solution;
  su2double X0[3] = {0.0,0.0,0.0}, X1[3] = {0.0,0.0,0.0}, X2[3] = {0.0,0.0,0.0},
  X1_X0[3] = {0.0,0.0,0.0}, X2_X0[3] = {0.0,0.0,0.0}, X2_X1[3] = {0.0,0.0,0.0},
  CP[3] = {0.0,0.0,0.0}, Distance, DotCheck, Radius;
  
  unsigned short nDim = geometry[MESH_0]->GetnDim();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == RANS) ||
               (config->GetKind_Solver() == ADJ_RANS) ||
               (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool SubsonicEngine = config->GetSubsonicEngine();
  

  /*--- Set subsonic initial condition for engine intakes ---*/
  
  if (SubsonicEngine) {
    
    /*--- Set initial boundary condition at iteration 0 ---*/
    
    if ((ExtIter == 0) && (!restart)) {
      
      su2double Velocity_Cyl[3] = {0.0, 0.0, 0.0}, Velocity_CylND[3] = {0.0, 0.0, 0.0}, Viscosity_Cyl,
      Density_Cyl, Density_CylND, Pressure_CylND, ModVel_Cyl, ModVel_CylND, Energy_CylND,
      T_ref = 0.0, S = 0.0, Mu_ref = 0.0, *Coord, *SubsonicEngine_Cyl, *SubsonicEngine_Values;
      
      SubsonicEngine_Values = config->GetSubsonicEngine_Values();
      su2double Mach_Cyl        = SubsonicEngine_Values[0];
      su2double Alpha_Cyl       = SubsonicEngine_Values[1];
      su2double Beta_Cyl        = SubsonicEngine_Values[2];
      su2double Pressure_Cyl    = SubsonicEngine_Values[3];
      su2double Temperature_Cyl = SubsonicEngine_Values[4];
      
      su2double Alpha = Alpha_Cyl*PI_NUMBER/180.0;
      su2double Beta  = Beta_Cyl*PI_NUMBER/180.0;
      
      su2double Gamma_Minus_One = Gamma - 1.0;
      su2double Gas_Constant = config->GetGas_Constant();
      
      su2double Mach2Vel_Cyl = sqrt(Gamma*Gas_Constant*Temperature_Cyl);
      
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          
          Velocity_Cyl[0] = cos(Alpha)*cos(Beta)*Mach_Cyl*Mach2Vel_Cyl;
          Velocity_Cyl[1] = sin(Beta)*Mach_Cyl*Mach2Vel_Cyl;
          Velocity_Cyl[2] = sin(Alpha)*cos(Beta)*Mach_Cyl*Mach2Vel_Cyl;
          
          ModVel_Cyl = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            ModVel_Cyl += Velocity_Cyl[iDim]*Velocity_Cyl[iDim];
          }
          ModVel_Cyl = sqrt(ModVel_Cyl);
          
          if (config->GetViscous()) {
            if (config->GetSystemMeasurements() == SI) { T_ref = 273.15; S = 110.4; Mu_ref = 1.716E-5; }
            if (config->GetSystemMeasurements() == US) {
              T_ref = (273.15 - 273.15) * 1.8 + 491.67;
              S = (110.4 - 273.15) * 1.8 + 491.67;
              Mu_ref = 1.716E-5/47.88025898;
            }
            Viscosity_Cyl = Mu_ref*(pow(Temperature_Cyl/T_ref, 1.5) * (T_ref+S)/(Temperature_Cyl+S));
            Density_Cyl   = config->GetReynolds()*Viscosity_Cyl/(ModVel_Cyl*config->GetLength_Reynolds());
            Pressure_Cyl  = Density_Cyl*Gas_Constant*Temperature_Cyl;
          }
          else {
            Density_Cyl = Pressure_Cyl/(Gas_Constant*Temperature_Cyl);
          }
          
          Density_CylND  = Density_Cyl/config->GetDensity_Ref();
          Pressure_CylND = Pressure_Cyl/config->GetPressure_Ref();
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_CylND[iDim] = Velocity_Cyl[iDim]/config->GetVelocity_Ref();
          }
          
          ModVel_CylND = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            ModVel_CylND += Velocity_CylND[iDim]*Velocity_CylND[iDim];
          }
          ModVel_CylND = sqrt(ModVel_CylND);
          
          Energy_CylND = Pressure_CylND/(Density_CylND*Gamma_Minus_One)+0.5*ModVel_CylND*ModVel_CylND;
          
          Coord = geometry[iMesh]->node[iPoint]->GetCoord();
          
          SubsonicEngine_Cyl = config->GetSubsonicEngine_Cyl();
          
          X0[0] = Coord[0];                             X0[1] = Coord[1];                           X0[2] = Coord[2];
          X1[0] = SubsonicEngine_Cyl[0]; X1[1] = SubsonicEngine_Cyl[1]; X1[2] = SubsonicEngine_Cyl[2];
          X2[0] = SubsonicEngine_Cyl[3]; X2[1] = SubsonicEngine_Cyl[4]; X2[2] = SubsonicEngine_Cyl[5];
          Radius = SubsonicEngine_Cyl[6];
          
          for (iDim = 0; iDim < nDim; iDim++) {
            X2_X1[iDim]= X1[iDim] - X2[iDim];
            X1_X0[iDim]= X0[iDim] - X1[iDim];
            X2_X0[iDim]= X0[iDim] - X2[iDim];
          }
          
          CP[0] = (X2_X1[1]*X1_X0[2] - X2_X1[2]*X1_X0[1]);
          CP[1] = (X2_X1[2]*X1_X0[0] - X2_X1[0]*X1_X0[2]);
          CP[2] = (X2_X1[0]*X1_X0[1] - X2_X1[1]*X1_X0[0]);
          
          Distance = sqrt((CP[0]*CP[0]+CP[1]*CP[1]+CP[2]*CP[2])/(X2_X1[0]*X2_X1[0]+X2_X1[1]*X2_X1[1]+X2_X1[2]*X2_X1[2]));
          
          DotCheck = -(X1_X0[0]*X2_X1[0]+X1_X0[1]*X2_X1[1]+X1_X0[2]*X2_X1[2]);
          if (DotCheck < 0.0) Distance = sqrt(X1_X0[0]*X1_X0[0]+X1_X0[1]*X1_X0[1]+X1_X0[2]*X1_X0[2]);
          
          DotCheck = (X2_X0[0]*X2_X1[0]+X2_X0[1]*X2_X1[1]+X2_X0[2]*X2_X1[2]);
          if (DotCheck < 0.0) Distance = sqrt(X2_X0[0]*X2_X0[0]+X2_X0[1]*X2_X0[1]+X2_X0[2]*X2_X0[2]);
          
          if (Distance < Radius) {
            
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(0, Density_CylND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(iDim+1, Density_CylND*Velocity_CylND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(nVar-1, Density_CylND*Energy_CylND);
            
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(0, Density_CylND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(iDim+1, Density_CylND*Velocity_CylND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(nVar-1, Density_CylND*Energy_CylND);
            
          }
          
        }
        
        /*--- Set the MPI communication ---*/
        
        solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
        solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution_Old(geometry[iMesh], config);
        
      }
      
    }
    
  }
  
  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  
  if (restart && (ExtIter == 0)) {
    
    Solution = new su2double[nVar];
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          Solution_Fine = solver_container[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
      }
      solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    }
    delete [] Solution;
    
    /*--- Interpolate the turblence variable also, if needed ---*/
    
    if (rans) {
      
      unsigned short nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
      Solution = new su2double[nVar_Turb];
      for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
          for (iVar = 0; iVar < nVar_Turb; iVar++) Solution[iVar] = 0.0;
          for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
            Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
            Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
            Solution_Fine = solver_container[iMesh-1][TURB_SOL]->node[Point_Fine]->GetSolution();
            for (iVar = 0; iVar < nVar_Turb; iVar++) {
              Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
            }
          }
          solver_container[iMesh][TURB_SOL]->node[iPoint]->SetSolution(Solution);
        }
        solver_container[iMesh][TURB_SOL]->Set_MPI_Solution(geometry[iMesh], config);
        solver_container[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver_container[iMesh], config, iMesh);
      }
      delete [] Solution;
    }
    
  }
  
  /*--- The value of the solution for the first iteration of the dual time ---*/
  
  if (dual_time && (ExtIter == 0 || (restart && (long)ExtIter == config->GetUnst_RestartIter()))) {
    
    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/
    
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
        if (rans) {
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
        }
      }
    }
    
    if ((restart && (long)ExtIter == config->GetUnst_RestartIter()) &&
        (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
      
      /*--- Load an additional restart file for a 2nd-order restart ---*/
      
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1));
      
      /*--- Load an additional restart file for the turbulence model ---*/
      if (rans)
        solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1));
      
      /*--- Push back this new solution to time level N. ---*/
      
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          if (rans) {
            solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          }
        }
      }
    }
  }
}

void CEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long ErrorCounter = 0;
  
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned long ExtIter = config->GetExtIter();
  bool adjoint          = config->GetContinuous_Adjoint();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool low_fidelity     = (config->GetLowFidelitySim() && (iMesh == MESH_1));
  bool second_order     = ((config->GetSpatialOrder_Flow() == SECOND_ORDER) || (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
  bool limiter          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (!low_fidelity) && (ExtIter <= config->GetLimiterIter()));
  bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst       = center && (config->GetKind_Centered_Flow() == JST);
  bool engine           = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk    = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool nearfield        = (config->GetnMarker_NearFieldBound() != 0);
  bool interface        = (config->GetnMarker_InterfaceBound() != 0);
  bool marker_analyze   = (config->GetnMarker_Analyze() != 0);
  bool fixed_cl         = config->GetFixed_CL_Mode();

  /*--- Update the angle of attack at the far-field for fixed CL calculations. ---*/
  
  if (fixed_cl) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }

  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
 
  /*--- Set the primitive variables ---*/

  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Compute the engine properties ---*/

  if (engine) { GetPower_Properties(geometry, config, iMesh, Output); }

  /*--- Compute the control volume properties ---*/

  if (marker_analyze) {
    GetSurface_Properties(geometry, NULL, NULL, config, iMesh, Output);
    GetSurface_Distortion(geometry, config, iMesh, Output);
  }

  /*--- Compute the actuator disk properties and distortion levels ---*/

  if (actuator_disk) {
    Set_MPI_ActDisk(solver_container, geometry, config);
    GetPower_Properties(geometry, config, iMesh, Output);
    SetActDisk_BCThrust(geometry, solver_container, config, iMesh, Output);
  }

  /*--- Compute Interface MPI ---*/

  if (interface) { Set_MPI_Interface(geometry, config); }

  /*--- Compute NearField MPI ---*/

  if (nearfield) { Set_MPI_Nearfield(geometry, config); }

 
  /*--- Upwind second order reconstruction ---*/
  
  if ((second_order && !center) && ((iMesh == MESH_0) || low_fidelity) && !Output) {
    
    /*--- Gradient computation ---*/
    
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config);
      //        if (compressible && !ideal_gas) SetSecondary_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config);
      //        if (compressible && !ideal_gas) SetSecondary_Gradient_LS(geometry, config);
    }
    
    
    /*--- Limiter computation ---*/
    
    if ((limiter) && (iMesh == MESH_0) && !Output) {
      SetPrimitive_Limiter(geometry, config);
      //        if (compressible && !ideal_gas) SetSecondary_Limiter(geometry, config);
    }
    
  }
  
  /*--- Artificial dissipation ---*/
  
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && ((iMesh == MESH_0) || low_fidelity)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !config->GetDiscrete_Adjoint()) Jacobian.SetValZero();

  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CEulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) { }

unsigned long CEulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
    
    RightSol = node[iPoint]->SetPrimVar(FluidModel);
    node[iPoint]->SetSecondaryVar(FluidModel);

    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
}
void CEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {
  
  su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
    bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetMax_Lambda_Inv(0.0);
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
      
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
    
  }
  
  
  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;
    
    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }
  
  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  
  if (time_steping) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Sets the regular CFL equal to the unsteady CFL ---*/
            config->SetCFL(iMesh,config->GetUnst_CFL());
            
            /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
             it computes the time step based on the unsteady CFL ---*/
            if (config->GetCFL(iMesh) == 0.0) {
                node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
            } else {
                node[iPoint]->SetDelta_Time(Global_Delta_Time);
            }
        }
  }
  
  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
    
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }
  
  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool second_order = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));
  bool low_fidelity = (config->GetLowFidelitySim() && (iMesh == MESH_1));
  bool grid_movement = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, set normal vectors, and number of neighbors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
    /*--- Set primitive variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    
    /*--- Set the largest convective eigenvalue ---*/
    
    numerics->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());
    
    /*--- Set undivided laplacian an pressure based sensor ---*/
    
    if ((second_order || low_fidelity)) {
      numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
      numerics->SetSensor(node[iPoint]->GetSensor(), node[jPoint]->GetSensor());
    }
    
    /*--- Grid movement ---*/
    
    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- Compute residuals, and Jacobians ---*/
    
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    
    /*--- Update convective and artificial dissipation residuals ---*/
    
    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Set implicit computation ---*/
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
  }
  
}

void CEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  su2double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, RoeVelocity[3] = {0.0,0.0,0.0}, R, sq_vel, RoeEnthalpy,
  *V_i, *V_j, *S_i, *S_j, *Limiter_i = NULL, *Limiter_j = NULL, sqvel, Non_Physical = 1.0;
  
  su2double z, velocity2_i, velocity2_j, mach_i, mach_j, vel_i_corr[3], vel_j_corr[3];
  
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;
  
  bool neg_density_i = false, neg_density_j = false, neg_pressure_i = false, neg_pressure_j = false, neg_sound_speed = false;
  
  
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool low_fidelity     = (config->GetLowFidelitySim() && (iMesh == MESH_1));
  bool second_order     = (((config->GetSpatialOrder_Flow() == SECOND_ORDER) || (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER)) && ((iMesh == MESH_0) || low_fidelity));
  bool limiter          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && !low_fidelity);
  bool grid_movement    = config->GetGrid_Movement();
  bool roe_turkel       = (config->GetKind_Upwind_Flow() == TURKEL);
  bool ideal_gas        = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );
  bool low_mach_corr    = config->Low_Mach_Correction();

  /*--- Loop over all the edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Roe Turkel preconditioning ---*/
    
    if (roe_turkel) {
      sqvel = 0.0;
      for (iDim = 0; iDim < nDim; iDim ++)
        sqvel += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
      numerics->SetVelocity2_Inf(sqvel);
    }
    
    /*--- Grid movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    /*--- Get primitive variables ---*/
    
    V_i = node[iPoint]->GetPrimitive(); V_j = node[jPoint]->GetPrimitive();
    S_i = node[iPoint]->GetSecondary(); S_j = node[jPoint]->GetSecondary();

    /*--- High order reconstruction using MUSCL strategy ---*/
    
    if (second_order) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      if (limiter) {
        Limiter_i = node[iPoint]->GetLimiter_Primitive();
        Limiter_j = node[jPoint]->GetLimiter_Primitive();
      }
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        Non_Physical = node[iPoint]->GetNon_Physical()*node[jPoint]->GetNon_Physical();
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim]*Non_Physical;
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim]*Non_Physical;
        }
        if (limiter) {
          Primitive_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Primitive_i[iVar] = V_i[iVar] + Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }

      /*--- Recompute the extrapolated quantities in a
       thermodynamic consistent way  ---*/

      if (!ideal_gas || low_mach_corr) { ComputeConsExtrapolation(config); }

      /*--- Low-Mach number correction ---*/

      if (low_mach_corr) {

        velocity2_i = 0.0;
        velocity2_j = 0.0;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          velocity2_i += Primitive_i[iDim+1]*Primitive_i[iDim+1];
          velocity2_j += Primitive_j[iDim+1]*Primitive_j[iDim+1];
        }
        mach_i = sqrt(velocity2_i)/Primitive_i[nDim+4];
        mach_j = sqrt(velocity2_j)/Primitive_j[nDim+4];

        z = min(max(mach_i,mach_j),1.0);
        velocity2_i = 0.0;
        velocity2_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
            vel_i_corr[iDim] = ( Primitive_i[iDim+1] + Primitive_j[iDim+1] )/2.0 \
                    + z * ( Primitive_i[iDim+1] - Primitive_j[iDim+1] )/2.0;
            vel_j_corr[iDim] = ( Primitive_i[iDim+1] + Primitive_j[iDim+1] )/2.0 \
                    + z * ( Primitive_j[iDim+1] - Primitive_i[iDim+1] )/2.0;

            velocity2_j += vel_j_corr[iDim]*vel_j_corr[iDim];
            velocity2_i += vel_i_corr[iDim]*vel_i_corr[iDim];

            Primitive_i[iDim+1] = vel_i_corr[iDim];
            Primitive_j[iDim+1] = vel_j_corr[iDim];
        }

        FluidModel->SetEnergy_Prho(Primitive_i[nDim+1],Primitive_i[nDim+2]);
        Primitive_i[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_i[nDim+1]/Primitive_i[nDim+2] + 0.5*velocity2_i;
        FluidModel->SetEnergy_Prho(Primitive_j[nDim+1],Primitive_j[nDim+2]);
        Primitive_j[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_j[nDim+1]/Primitive_j[nDim+2] + 0.5*velocity2_j;
      }
      
      /*--- Check for non-physical solutions after reconstruction. If found,
       use the cell-average value of the solution. This results in a locally
       first-order approximation, but this is typically only active
       during the start-up of a calculation. If non-physical, use the 
       cell-averaged state. ---*/
      
      neg_pressure_i = (Primitive_i[nDim+1] < 0.0); neg_pressure_j = (Primitive_j[nDim+1] < 0.0);
      neg_density_i  = (Primitive_i[nDim+2] < 0.0); neg_density_j  = (Primitive_j[nDim+2] < 0.0);

      R = sqrt(fabs(Primitive_j[nDim+2]/Primitive_i[nDim+2]));
      sq_vel = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        RoeVelocity[iDim] = (R*Primitive_j[iDim+1]+Primitive_i[iDim+1])/(R+1);
        sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
      }
      RoeEnthalpy = (R*Primitive_j[nDim+3]+Primitive_i[nDim+3])/(R+1);
      neg_sound_speed = ((Gamma-1)*(RoeEnthalpy-0.5*sq_vel) < 0.0);
      
      if (neg_sound_speed) {
        for (iVar = 0; iVar < nPrimVar; iVar++) {
          Primitive_i[iVar] = V_i[iVar];
          Primitive_j[iVar] = V_j[iVar]; }
        Secondary_i[0] = S_i[0]; Secondary_i[1] = S_i[1];
        Secondary_j[0] = S_i[0]; Secondary_j[1] = S_i[1];
        counter_local++;
      }
      
      if (neg_density_i || neg_pressure_i) {
        for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = V_i[iVar];
        Secondary_i[0] = S_i[0]; Secondary_i[1] = S_i[1];
        counter_local++;
      }
      
      if (neg_density_j || neg_pressure_j) {
        for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = V_j[iVar];
        Secondary_j[0] = S_j[0]; Secondary_j[1] = S_j[1];
        counter_local++;
      }

      numerics->SetPrimitive(Primitive_i, Primitive_j);
      numerics->SetSecondary(Secondary_i, Secondary_j);
      
    }
    else {
      
      /*--- Set conservative variables without reconstruction ---*/
      
      numerics->SetPrimitive(V_i, V_j);
      numerics->SetSecondary(S_i, S_j);
      
    }
    
    /*--- Compute the residual ---*/
    
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

    /*--- Update residual value ---*/
    
    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Set implicit Jacobians ---*/
    
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
    
    /*--- Roe Turkel preconditioning, set the value of beta ---*/
    
    if (roe_turkel) {
      node[iPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
      node[jPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
    }
    
  }
  
  /*--- Warning message about non-physical reconstructions ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Reconstr(counter_global);
  }
  
}

void CEulerSolver::ComputeConsExtrapolation(CConfig *config) {
  
  unsigned short iDim;
  
  su2double density_i = Primitive_i[nDim+2];
  su2double pressure_i = Primitive_i[nDim+1];
  su2double velocity2_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    velocity2_i += Primitive_i[iDim+1]*Primitive_i[iDim+1];
  }
  
  FluidModel->SetTDState_Prho(pressure_i, density_i);
  
  Primitive_i[0]= FluidModel->GetTemperature();
  Primitive_i[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_i[nDim+1]/Primitive_i[nDim+2] + 0.5*velocity2_i;
  Primitive_i[nDim+4]= FluidModel->GetSoundSpeed();
  Secondary_i[0]=FluidModel->GetdPdrho_e();
  Secondary_i[1]=FluidModel->GetdPde_rho();
  
  
  su2double density_j = Primitive_j[nDim+2];
  su2double pressure_j = Primitive_j[nDim+1];
  su2double velocity2_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    velocity2_j += Primitive_j[iDim+1]*Primitive_j[iDim+1];
  }
  
  FluidModel->SetTDState_Prho(pressure_j, density_j);
  
  Primitive_j[0]= FluidModel->GetTemperature();
  Primitive_j[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_j[nDim+1]/Primitive_j[nDim+2] + 0.5*velocity2_j;
  Primitive_j[nDim+4]=FluidModel->GetSoundSpeed();
  Secondary_j[0]=FluidModel->GetdPdrho_e();
  Secondary_j[1]=FluidModel->GetdPde_rho();
  
}

void CEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
  bool axisymmetric   = config->GetAxisymmetric();
  bool gravity        = (config->GetGravityForce() == YES);
  bool harmonic_balance  = (config->GetUnsteady_Simulation() == HARMONIC_BALANCE);
  bool windgust       = config->GetWind_Gust();
  
  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  
  if (rotating_frame) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the conservative variables ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());
      
      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Set solution  ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());

      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set y coordinate ---*/
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      
      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }
  
  if (gravity) {
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Set solution  ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, config);
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
    
  }
  
  if (harmonic_balance) {
    
    su2double Volume, Source;
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Get control volume ---*/
      Volume = geometry->node[iPoint]->GetVolume();
      
      /*--- Get stored time spectral source term ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Source = node[iPoint]->GetHarmonicBalance_Source(iVar);        
        Residual[iVar] = Source*Volume;
      }
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }
  
  if (windgust) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the wind gust ---*/
      numerics->SetWindGust(node[iPoint]->GetWindGust(), node[iPoint]->GetWindGust());
      
      /*--- Load the wind gust derivatives ---*/
      numerics->SetWindGustDer(node[iPoint]->GetWindGustDer(), node[iPoint]->GetWindGustDer());
      
      /*--- Load the primitive variables ---*/
      numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[iPoint]->GetPrimitive());
      
      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
}

void CEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  /* This method should be used to call any new source terms for a particular problem*/
  /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */
  
  /* Next we describe how to get access to some important quanties for this method */
  /* Access to all points in the current geometric mesh by saying: nPointDomain */
  /* Get the vector of conservative variables at some point iPoint = node[iPoint]->GetSolution() */
  /* Get the volume (or area in 2D) associated with iPoint = node[iPoint]->GetVolume() */
  /* Get the vector of geometric coordinates of point iPoint = node[iPoint]->GetCoord() */
  
}

void CEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
  
  su2double *Normal, Area, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda,
  ProjVel, ProjVel_i, ProjVel_j, *GridVel, *GridVel_i, *GridVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  bool grid_movement = config->GetGrid_Movement();
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetLambda(0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;

    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddLambda(Lambda);
      }
      
    }
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_MaxEigenvalue(geometry, config);
  
}

void CEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, jPoint, iEdge;
  su2double Pressure_i = 0, Pressure_j = 0, *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;
    
  Diff = new su2double[nVar];
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetUnd_LaplZero();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Solution differences ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
    
    /*--- Correction for compressible flows which use the enthalpy ---*/
    
    Pressure_i = node[iPoint]->GetPressure();
    Pressure_j = node[jPoint]->GetPressure();
    Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + Pressure_i) - (node[jPoint]->GetSolution(nVar-1) + Pressure_j);
    
    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both in the boundary ---*/
    
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_Undivided_Laplacian(geometry, config);
  
  delete [] Diff;
  
}

void CEulerSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  su2double Pressure_i = 0.0, Pressure_j = 0.0;
  bool boundary_i, boundary_j;
  
  /*--- Reset variables to store the undivided pressure ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    iPoint_UndLapl[iPoint] = 0.0;
    jPoint_UndLapl[iPoint] = 0.0;
  }
  
  /*--- Evaluate the pressure sensor ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Pressure_i = node[iPoint]->GetPressure();
    Pressure_j = node[jPoint]->GetPressure();
    
    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both on the boundary ---*/
    
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    
  }
  
  /*--- Set pressure switch for each point ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetSensor(fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_Dissipation_Switch(geometry, config);
  
}

void CEulerSolver::Pressure_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double Pressure = 0.0, *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord,
  factor, NFPressOF, RefVel2, RefTemp, RefDensity, RefPressure, Mach2Vel, Mach_Motion,
  Force[3] = {0.0,0.0,0.0};
  string Marker_Tag, Monitoring_Tag;
  su2double AxiFactor;
  
#ifdef HAVE_MPI
  su2double MyAllBound_CD_Inv, MyAllBound_CL_Inv, MyAllBound_CSF_Inv, MyAllBound_CMx_Inv, MyAllBound_CMy_Inv, MyAllBound_CMz_Inv, MyAllBound_CFx_Inv, MyAllBound_CFy_Inv, MyAllBound_CFz_Inv, MyAllBound_CT_Inv, MyAllBound_CQ_Inv, MyAllBound_CNearFieldOF_Inv, *MySurface_CL_Inv = NULL, *MySurface_CD_Inv = NULL, *MySurface_CSF_Inv = NULL, *MySurface_CEff_Inv = NULL, *MySurface_CFx_Inv = NULL, *MySurface_CFy_Inv = NULL, *MySurface_CFz_Inv = NULL, *MySurface_CMx_Inv = NULL, *MySurface_CMy_Inv = NULL, *MySurface_CMz_Inv = NULL;
#endif
  
  su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefAreaCoeff    = config->GetRefAreaCoeff();
  su2double RefLengthMoment = config->GetRefLengthMoment();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double *Origin         = config->GetRefOriginMoment(0);
  bool grid_movement        = config->GetGrid_Movement();
  bool axisymmetric         = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream
   values, which is the standard convention. ---*/
  
  RefTemp     = Temperature_Inf;
  RefDensity  = Density_Inf;
  RefPressure = Pressure_Inf;
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  /*-- Variables initialization ---*/
  
  Total_CD = 0.0;           Total_CL = 0.0;    Total_CSF = 0.0;     Total_CEff = 0.0;
  Total_CMx = 0.0;          Total_CMy = 0.0;   Total_CMz = 0.0;
  Total_CFx = 0.0;          Total_CFy = 0.0;   Total_CFz = 0.0;
  Total_CT = 0.0;           Total_CQ = 0.0;    Total_CMerit = 0.0;
  Total_CNearFieldOF = 0.0; Total_Heat = 0.0;  Total_MaxHeat = 0.0;
  
  AllBound_CD_Inv = 0.0;        AllBound_CL_Inv = 0.0; AllBound_CSF_Inv = 0.0;
  AllBound_CMx_Inv = 0.0;          AllBound_CMy_Inv = 0.0;   AllBound_CMz_Inv = 0.0;
  AllBound_CFx_Inv = 0.0;          AllBound_CFy_Inv = 0.0;   AllBound_CFz_Inv = 0.0;
  AllBound_CT_Inv = 0.0;           AllBound_CQ_Inv = 0.0;    AllBound_CMerit_Inv = 0.0;
  AllBound_CNearFieldOF_Inv = 0.0; AllBound_CEff_Inv = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Inv[iMarker_Monitoring]      = 0.0; Surface_CD_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring] = 0.0; Surface_CEff_Inv[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0; Surface_CFy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0; Surface_CMx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0; Surface_CMz_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CL[iMarker_Monitoring]          = 0.0; Surface_CD[iMarker_Monitoring]          = 0.0;
    Surface_CSF[iMarker_Monitoring]     = 0.0; Surface_CEff[iMarker_Monitoring]           = 0.0;
    Surface_CFx[iMarker_Monitoring]            = 0.0; Surface_CFy[iMarker_Monitoring]            = 0.0;
    Surface_CFz[iMarker_Monitoring]            = 0.0; Surface_CMx[iMarker_Monitoring]            = 0.0;
    Surface_CMy[iMarker_Monitoring]            = 0.0; Surface_CMz[iMarker_Monitoring]            = 0.0;
  }
  
  /*--- Loop over the Euler and Navier-Stokes markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
        (Boundary == ISOTHERMAL) || (Boundary == NEARFIELD_BOUNDARY) ||
        (Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) ||
        (Boundary == ACTDISK_INLET) || (Boundary == ACTDISK_OUTLET)||
        (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CD_Inv[iMarker] = 0.0;        CL_Inv[iMarker] = 0.0; CSF_Inv[iMarker] = 0.0;
      CMx_Inv[iMarker] = 0.0;          CMy_Inv[iMarker] = 0.0;   CMz_Inv[iMarker] = 0.0;
      CFx_Inv[iMarker] = 0.0;          CFy_Inv[iMarker] = 0.0;   CFz_Inv[iMarker] = 0.0;
      CT_Inv[iMarker] = 0.0;           CQ_Inv[iMarker] = 0.0;    CMerit_Inv[iMarker] = 0.0;
      CNearFieldOF_Inv[iMarker] = 0.0; CEff_Inv[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;
      MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
      NFPressOF = 0.0;
      
      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        Pressure = node[iPoint]->GetPressure();
        
        CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefAreaCoeff;
        
        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          
          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/
          
          NFPressOF += 0.5*(Pressure - Pressure_Inf)*(Pressure - Pressure_Inf)*Normal[nDim-1];
          
          for (iDim = 0; iDim < nDim; iDim++) {
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }
          
          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf) * Normal[iDim] * factor * AxiFactor;
            ForceInviscid[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (nDim == 3) {
            MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLengthMoment;
            MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
          }
          MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLengthMoment;
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if (Monitoring == YES) {
        
        if (Boundary != NEARFIELD_BOUNDARY) {
          if (nDim == 2) {
            CD_Inv[iMarker]  =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
            CL_Inv[iMarker]  = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
            CEff_Inv[iMarker]   = CL_Inv[iMarker] / (CD_Inv[iMarker]+EPS);
            CMz_Inv[iMarker]    = MomentInviscid[2];
            CFx_Inv[iMarker]    = ForceInviscid[0];
            CFy_Inv[iMarker]    = ForceInviscid[1];
            CT_Inv[iMarker]     = -CFx_Inv[iMarker];
            CQ_Inv[iMarker]     = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker] = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
          }
          if (nDim == 3) {
            CD_Inv[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
            CL_Inv[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
            CSF_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
            CEff_Inv[iMarker]       = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
            CMx_Inv[iMarker]        = MomentInviscid[0];
            CMy_Inv[iMarker]        = MomentInviscid[1];
            CMz_Inv[iMarker]        = MomentInviscid[2];
            CFx_Inv[iMarker]        = ForceInviscid[0];
            CFy_Inv[iMarker]        = ForceInviscid[1];
            CFz_Inv[iMarker]        = ForceInviscid[2];
            CT_Inv[iMarker]         = -CFz_Inv[iMarker];
            CQ_Inv[iMarker]         = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker]     = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
          }
          
          AllBound_CD_Inv        += CD_Inv[iMarker];
          AllBound_CL_Inv        += CL_Inv[iMarker];
          AllBound_CSF_Inv   += CSF_Inv[iMarker];
          AllBound_CEff_Inv          = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
          AllBound_CMx_Inv          += CMx_Inv[iMarker];
          AllBound_CMy_Inv          += CMy_Inv[iMarker];
          AllBound_CMz_Inv          += CMz_Inv[iMarker];
          AllBound_CFx_Inv          += CFx_Inv[iMarker];
          AllBound_CFy_Inv          += CFy_Inv[iMarker];
          AllBound_CFz_Inv          += CFz_Inv[iMarker];
          AllBound_CT_Inv           += CT_Inv[iMarker];
          AllBound_CQ_Inv           += CQ_Inv[iMarker];
          AllBound_CMerit_Inv        = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
          
          /*--- Compute the coefficients per surface ---*/
          
          for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
            Marker_Tag = config->GetMarker_All_TagBound(iMarker);
            if (Marker_Tag == Monitoring_Tag) {
              Surface_CL_Inv[iMarker_Monitoring]      += CL_Inv[iMarker];
              Surface_CD_Inv[iMarker_Monitoring]      += CD_Inv[iMarker];
              Surface_CSF_Inv[iMarker_Monitoring] += CSF_Inv[iMarker];
              Surface_CEff_Inv[iMarker_Monitoring]        = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
              Surface_CFx_Inv[iMarker_Monitoring]        += CFx_Inv[iMarker];
              Surface_CFy_Inv[iMarker_Monitoring]        += CFy_Inv[iMarker];
              Surface_CFz_Inv[iMarker_Monitoring]        += CFz_Inv[iMarker];
              Surface_CMx_Inv[iMarker_Monitoring]        += CMx_Inv[iMarker];
              Surface_CMy_Inv[iMarker_Monitoring]        += CMy_Inv[iMarker];
              Surface_CMz_Inv[iMarker_Monitoring]        += CMz_Inv[iMarker];
            }
          }
          
        }
        
        /*--- At the Nearfield SU2 only cares about the pressure coeffient ---*/
        
        else {
          CNearFieldOF_Inv[iMarker] = NFPressOF;
          AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];
        }
        
      }
      
      
    }
  }
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  MyAllBound_CD_Inv        = AllBound_CD_Inv;        AllBound_CD_Inv = 0.0;
  MyAllBound_CL_Inv        = AllBound_CL_Inv;        AllBound_CL_Inv = 0.0;
  MyAllBound_CSF_Inv   = AllBound_CSF_Inv;   AllBound_CSF_Inv = 0.0;
  AllBound_CEff_Inv = 0.0;
  MyAllBound_CMx_Inv          = AllBound_CMx_Inv;          AllBound_CMx_Inv = 0.0;
  MyAllBound_CMy_Inv          = AllBound_CMy_Inv;          AllBound_CMy_Inv = 0.0;
  MyAllBound_CMz_Inv          = AllBound_CMz_Inv;          AllBound_CMz_Inv = 0.0;
  MyAllBound_CFx_Inv          = AllBound_CFx_Inv;          AllBound_CFx_Inv = 0.0;
  MyAllBound_CFy_Inv          = AllBound_CFy_Inv;          AllBound_CFy_Inv = 0.0;
  MyAllBound_CFz_Inv          = AllBound_CFz_Inv;          AllBound_CFz_Inv = 0.0;
  MyAllBound_CT_Inv           = AllBound_CT_Inv;           AllBound_CT_Inv = 0.0;
  MyAllBound_CQ_Inv           = AllBound_CQ_Inv;           AllBound_CQ_Inv = 0.0;
  AllBound_CMerit_Inv = 0.0;
  MyAllBound_CNearFieldOF_Inv = AllBound_CNearFieldOF_Inv; AllBound_CNearFieldOF_Inv = 0.0;
  
  SU2_MPI::Allreduce(&MyAllBound_CD_Inv, &AllBound_CD_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CL_Inv, &AllBound_CL_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CSF_Inv, &AllBound_CSF_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Inv = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CMx_Inv, &AllBound_CMx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMy_Inv, &AllBound_CMy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMz_Inv, &AllBound_CMz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFx_Inv, &AllBound_CFx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFy_Inv, &AllBound_CFy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFz_Inv, &AllBound_CFz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CT_Inv, &AllBound_CT_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CQ_Inv, &AllBound_CQ_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Inv = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CNearFieldOF_Inv, &AllBound_CNearFieldOF_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  /*--- Add the forces on the surfaces using all the nodes ---*/
  
  MySurface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Inv = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CL_Inv[iMarker_Monitoring]      = Surface_CL_Inv[iMarker_Monitoring];
    MySurface_CD_Inv[iMarker_Monitoring]      = Surface_CD_Inv[iMarker_Monitoring];
    MySurface_CSF_Inv[iMarker_Monitoring] = Surface_CSF_Inv[iMarker_Monitoring];
    MySurface_CEff_Inv[iMarker_Monitoring]       = Surface_CEff_Inv[iMarker_Monitoring];
    MySurface_CFx_Inv[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    MySurface_CFy_Inv[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    MySurface_CFz_Inv[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    MySurface_CMx_Inv[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    MySurface_CMy_Inv[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    MySurface_CMz_Inv[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];
    
    Surface_CL_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CD_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring] = 0.0;
    Surface_CEff_Inv[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Inv[iMarker_Monitoring]        = 0.0;
  }
  
  SU2_MPI::Allreduce(MySurface_CL_Inv, Surface_CL_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CD_Inv, Surface_CD_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CSF_Inv, Surface_CSF_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
    Surface_CEff_Inv[iMarker_Monitoring] = Surface_CL_Inv[iMarker_Monitoring] / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
  SU2_MPI::Allreduce(MySurface_CFx_Inv, Surface_CFx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFy_Inv, Surface_CFy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFz_Inv, Surface_CFz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMx_Inv, Surface_CMx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMy_Inv, Surface_CMy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMz_Inv, Surface_CMz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  delete [] MySurface_CL_Inv; delete [] MySurface_CD_Inv; delete [] MySurface_CSF_Inv;
  delete [] MySurface_CEff_Inv;  delete [] MySurface_CFx_Inv;   delete [] MySurface_CFy_Inv;
  delete [] MySurface_CFz_Inv;   delete [] MySurface_CMx_Inv;   delete [] MySurface_CMy_Inv;
  delete [] MySurface_CMz_Inv;
  
#endif
  
  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
  
  Total_CD            = AllBound_CD_Inv;
  Total_CL            = AllBound_CL_Inv;
  Total_CSF           = AllBound_CSF_Inv;
  Total_CEff          = Total_CL / (Total_CD + EPS);
  Total_CMx           = AllBound_CMx_Inv;
  Total_CMy           = AllBound_CMy_Inv;
  Total_CMz           = AllBound_CMz_Inv;
  Total_CFx           = AllBound_CFx_Inv;
  Total_CFy           = AllBound_CFy_Inv;
  Total_CFz           = AllBound_CFz_Inv;
  Total_CT            = AllBound_CT_Inv;
  Total_CQ            = AllBound_CQ_Inv;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);
  Total_CNearFieldOF  = AllBound_CNearFieldOF_Inv;
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]      = Surface_CL_Inv[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]      = Surface_CD_Inv[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring] = Surface_CSF_Inv[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CL_Inv[iMarker_Monitoring] / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];
  }
  
}

void CEulerSolver::Momentum_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord, Area,
  factor, RefVel2, RefTemp, RefDensity,  Mach2Vel, Mach_Motion,
  Force[3] = {0.0,0.0,0.0}, Velocity[3], MassFlow, Density;
  string Marker_Tag, Monitoring_Tag;
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0}, MomentZ_Force[3] = {0.0,0.0,0.0};
  su2double AxiFactor;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Mnt, MyAllBound_CL_Mnt, MyAllBound_CSF_Mnt,
MyAllBound_CMx_Mnt, MyAllBound_CMy_Mnt, MyAllBound_CMz_Mnt,
  MyAllBound_CFx_Mnt, MyAllBound_CFy_Mnt, MyAllBound_CFz_Mnt, MyAllBound_CT_Mnt,
  MyAllBound_CQ_Mnt,
  *MySurface_CL_Mnt = NULL, *MySurface_CD_Mnt = NULL, *MySurface_CSF_Mnt = NULL,
  *MySurface_CEff_Mnt = NULL, *MySurface_CFx_Mnt = NULL, *MySurface_CFy_Mnt = NULL,
  *MySurface_CFz_Mnt = NULL,
  *MySurface_CMx_Mnt = NULL, *MySurface_CMy_Mnt = NULL,  *MySurface_CMz_Mnt = NULL;
#endif
  
  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefAreaCoeff     = config->GetRefAreaCoeff();
  su2double RefLengthMoment  = config->GetRefLengthMoment();
  su2double Gas_Constant     = config->GetGas_ConstantND();
  su2double *Origin          = config->GetRefOriginMoment(0);
  bool grid_movement         = config->GetGrid_Movement();
  bool axisymmetric          = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  RefTemp     = Temperature_Inf;
  RefDensity  = Density_Inf;
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  /*-- Variables initialization ---*/
  
  AllBound_CD_Mnt = 0.0;        AllBound_CL_Mnt = 0.0; AllBound_CSF_Mnt = 0.0;
  AllBound_CMx_Mnt = 0.0;          AllBound_CMy_Mnt = 0.0;   AllBound_CMz_Mnt = 0.0;
  AllBound_CFx_Mnt = 0.0;          AllBound_CFy_Mnt = 0.0;   AllBound_CFz_Mnt = 0.0;
  AllBound_CT_Mnt = 0.0;           AllBound_CQ_Mnt = 0.0;    AllBound_CMerit_Mnt = 0.0;
  AllBound_CEff_Mnt = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Mnt[iMarker_Monitoring]      = 0.0; Surface_CD_Mnt[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Mnt[iMarker_Monitoring] = 0.0; Surface_CEff_Mnt[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Mnt[iMarker_Monitoring]        = 0.0; Surface_CFy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Mnt[iMarker_Monitoring]        = 0.0; Surface_CMy_Mnt[iMarker_Monitoring]        = 0.0; Surface_CMz_Mnt[iMarker_Monitoring]        = 0.0;
  }
  
  /*--- Loop over the Inlet -Outlet Markers  ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) ||
        (Boundary == ACTDISK_INLET) || (Boundary == ACTDISK_OUTLET)||
        (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CD_Mnt[iMarker] = 0.0;        CL_Mnt[iMarker] = 0.0; CSF_Mnt[iMarker] = 0.0;
      CMx_Mnt[iMarker] = 0.0;          CMy_Mnt[iMarker] = 0.0;   CMz_Mnt[iMarker] = 0.0;
      CFx_Mnt[iMarker] = 0.0;          CFy_Mnt[iMarker] = 0.0;   CFz_Mnt[iMarker] = 0.0;
      CT_Mnt[iMarker] = 0.0;           CQ_Mnt[iMarker] = 0.0;    CMerit_Mnt[iMarker] = 0.0;
      CEff_Mnt[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceMomentum[iDim] = 0.0;
      MomentMomentum[0] = 0.0; MomentMomentum[1] = 0.0; MomentMomentum[2] = 0.0;
      MomentX_Force[0] = 0.0; MomentX_Force[1] = 0.0; MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0; MomentY_Force[1] = 0.0; MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0; MomentZ_Force[1] = 0.0; MomentZ_Force[2] = 0.0;

      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          Density   = node[iPoint]->GetDensity();
          
          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/
          
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
          
          MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim]  = node[iPoint]->GetVelocity(iDim);
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
            MassFlow -= Normal[iDim]*Velocity[iDim]*Density;
          }
          
          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = MassFlow * Velocity[iDim] * factor * AxiFactor;
            ForceMomentum[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (iDim == 3) {
            MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLengthMoment;
                    MomentX_Force[1] += (-Force[1]*MomentDist[2])/RefLengthMoment;
                    MomentX_Force[2] += (Force[2]*MomentDist[1])/RefLengthMoment;

            MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
            MomentY_Force[0] += (Force[0]*MomentDist[2])/RefLengthMoment;
            MomentY_Force[2] += (-Force[2]*MomentDist[0])/RefLengthMoment;
          }
          MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLengthMoment;
          MomentZ_Force[0] += (-Force[0]*MomentDist[1])/RefLengthMoment;
                    MomentZ_Force[1] += (Force[1]*MomentDist[0])/RefLengthMoment;
          
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if (Monitoring == YES) {
        
        if (nDim == 2) {
          CD_Mnt[iMarker]  =  ForceMomentum[0]*cos(Alpha) + ForceMomentum[1]*sin(Alpha);
          CL_Mnt[iMarker]  = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[1]*cos(Alpha);
          CEff_Mnt[iMarker]   = CL_Mnt[iMarker] / (CD_Mnt[iMarker]+EPS);
          CMz_Mnt[iMarker]    = MomentInviscid[2];
          CFx_Mnt[iMarker]    = ForceMomentum[0];
          CFy_Mnt[iMarker]    = ForceMomentum[1];
          CT_Mnt[iMarker]     = -CFx_Mnt[iMarker];
          CQ_Mnt[iMarker]     = -CMz_Mnt[iMarker];
          CMerit_Mnt[iMarker] = CT_Mnt[iMarker] / (CQ_Mnt[iMarker] + EPS);
        }
        if (nDim == 3) {
          CD_Mnt[iMarker]      =  ForceMomentum[0]*cos(Alpha)*cos(Beta) + ForceMomentum[1]*sin(Beta) + ForceMomentum[2]*sin(Alpha)*cos(Beta);
          CL_Mnt[iMarker]      = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[2]*cos(Alpha);
          CSF_Mnt[iMarker] = -ForceMomentum[0]*sin(Beta)*cos(Alpha) + ForceMomentum[1]*cos(Beta) - ForceMomentum[2]*sin(Beta)*sin(Alpha);
          CEff_Mnt[iMarker]       = CL_Mnt[iMarker] / (CD_Mnt[iMarker] + EPS);
          CMx_Mnt[iMarker]        = MomentInviscid[0];
          CMy_Mnt[iMarker]        = MomentInviscid[1];
          CMz_Mnt[iMarker]        = MomentInviscid[2];
          CFx_Mnt[iMarker]        = ForceMomentum[0];
          CFy_Mnt[iMarker]        = ForceMomentum[1];
          CFz_Mnt[iMarker]        = ForceMomentum[2];
          CT_Mnt[iMarker]         = -CFz_Mnt[iMarker];
          CQ_Mnt[iMarker]         = -CMz_Mnt[iMarker];
          CMerit_Mnt[iMarker]     = CT_Mnt[iMarker] / (CQ_Mnt[iMarker] + EPS);
        }
        
        AllBound_CD_Mnt        += CD_Mnt[iMarker];
        AllBound_CL_Mnt        += CL_Mnt[iMarker];
        AllBound_CSF_Mnt   += CSF_Mnt[iMarker];
        AllBound_CEff_Mnt          = AllBound_CL_Mnt / (AllBound_CD_Mnt + EPS);
        AllBound_CMx_Mnt          += CMx_Mnt[iMarker];
        AllBound_CMy_Mnt          += CMy_Mnt[iMarker];
        AllBound_CMz_Mnt          += CMz_Mnt[iMarker];
        AllBound_CFx_Mnt          += CFx_Mnt[iMarker];
        AllBound_CFy_Mnt          += CFy_Mnt[iMarker];
        AllBound_CFz_Mnt          += CFz_Mnt[iMarker];
        AllBound_CT_Mnt           += CT_Mnt[iMarker];
        AllBound_CQ_Mnt           += CQ_Mnt[iMarker];
        AllBound_CMerit_Mnt        += AllBound_CT_Mnt / (AllBound_CQ_Mnt + EPS);
        
        /*--- Compute the coefficients per surface ---*/
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Mnt[iMarker_Monitoring]      += CL_Mnt[iMarker];
            Surface_CD_Mnt[iMarker_Monitoring]      += CD_Mnt[iMarker];
            Surface_CSF_Mnt[iMarker_Monitoring] += CSF_Mnt[iMarker];
            Surface_CEff_Mnt[iMarker_Monitoring]        = CL_Mnt[iMarker] / (CD_Mnt[iMarker] + EPS);
            Surface_CFx_Mnt[iMarker_Monitoring]        += CFx_Mnt[iMarker];
            Surface_CFy_Mnt[iMarker_Monitoring]        += CFy_Mnt[iMarker];
            Surface_CFz_Mnt[iMarker_Monitoring]        += CFz_Mnt[iMarker];
            Surface_CMx_Mnt[iMarker_Monitoring]        += CMx_Mnt[iMarker];
            Surface_CMy_Mnt[iMarker_Monitoring]        += CMy_Mnt[iMarker];
            Surface_CMz_Mnt[iMarker_Monitoring]        += CMz_Mnt[iMarker];
          }
        }
        
      }
      
      
    }
  }
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  MyAllBound_CD_Mnt        = AllBound_CD_Mnt;        AllBound_CD_Mnt = 0.0;
  MyAllBound_CL_Mnt        = AllBound_CL_Mnt;        AllBound_CL_Mnt = 0.0;
  MyAllBound_CSF_Mnt   = AllBound_CSF_Mnt;   AllBound_CSF_Mnt = 0.0;
  MyAllBound_CMx_Mnt          = AllBound_CMx_Mnt;          AllBound_CMx_Mnt = 0.0;
  MyAllBound_CMy_Mnt          = AllBound_CMy_Mnt;          AllBound_CMy_Mnt = 0.0;
  MyAllBound_CMz_Mnt          = AllBound_CMz_Mnt;          AllBound_CMz_Mnt = 0.0;
  MyAllBound_CFx_Mnt          = AllBound_CFx_Mnt;          AllBound_CFx_Mnt = 0.0;
  MyAllBound_CFy_Mnt          = AllBound_CFy_Mnt;          AllBound_CFy_Mnt = 0.0;
  MyAllBound_CFz_Mnt          = AllBound_CFz_Mnt;          AllBound_CFz_Mnt = 0.0;
  MyAllBound_CT_Mnt           = AllBound_CT_Mnt;           AllBound_CT_Mnt = 0.0;
  MyAllBound_CQ_Mnt           = AllBound_CQ_Mnt;           AllBound_CQ_Mnt = 0.0;
  
  SU2_MPI::Allreduce(&MyAllBound_CD_Mnt, &AllBound_CD_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CL_Mnt, &AllBound_CL_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CSF_Mnt, &AllBound_CSF_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Mnt = AllBound_CL_Mnt / (AllBound_CD_Mnt + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CMx_Mnt, &AllBound_CMx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMy_Mnt, &AllBound_CMy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMz_Mnt, &AllBound_CMz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFx_Mnt, &AllBound_CFx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFy_Mnt, &AllBound_CFy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFz_Mnt, &AllBound_CFz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CT_Mnt, &AllBound_CT_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CQ_Mnt, &AllBound_CQ_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Mnt = AllBound_CT_Mnt / (AllBound_CQ_Mnt + EPS);
  
  /*--- Add the forces on the surfaces using all the nodes ---*/
  
  MySurface_CL_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Mnt = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CL_Mnt[iMarker_Monitoring]      = Surface_CL_Mnt[iMarker_Monitoring];
    MySurface_CD_Mnt[iMarker_Monitoring]      = Surface_CD_Mnt[iMarker_Monitoring];
    MySurface_CSF_Mnt[iMarker_Monitoring] = Surface_CSF_Mnt[iMarker_Monitoring];
    MySurface_CEff_Mnt[iMarker_Monitoring]       = Surface_CEff_Mnt[iMarker_Monitoring];
    MySurface_CFx_Mnt[iMarker_Monitoring]        = Surface_CFx_Mnt[iMarker_Monitoring];
    MySurface_CFy_Mnt[iMarker_Monitoring]        = Surface_CFy_Mnt[iMarker_Monitoring];
    MySurface_CFz_Mnt[iMarker_Monitoring]        = Surface_CFz_Mnt[iMarker_Monitoring];
    MySurface_CMx_Mnt[iMarker_Monitoring]        = Surface_CMx_Mnt[iMarker_Monitoring];
    MySurface_CMy_Mnt[iMarker_Monitoring]        = Surface_CMy_Mnt[iMarker_Monitoring];
    MySurface_CMz_Mnt[iMarker_Monitoring]        = Surface_CMz_Mnt[iMarker_Monitoring];

    Surface_CL_Mnt[iMarker_Monitoring]      = 0.0;
    Surface_CD_Mnt[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Mnt[iMarker_Monitoring] = 0.0;
    Surface_CEff_Mnt[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Mnt[iMarker_Monitoring]        = 0.0;
  }
  
  SU2_MPI::Allreduce(MySurface_CL_Mnt, Surface_CL_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CD_Mnt, Surface_CD_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CSF_Mnt, Surface_CSF_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
    Surface_CEff_Mnt[iMarker_Monitoring] = Surface_CL_Mnt[iMarker_Monitoring] / (Surface_CD_Mnt[iMarker_Monitoring] + EPS);
  SU2_MPI::Allreduce(MySurface_CFx_Mnt, Surface_CFx_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFy_Mnt, Surface_CFy_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFz_Mnt, Surface_CFz_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMx_Mnt, Surface_CMx_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMy_Mnt, Surface_CMy_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMz_Mnt, Surface_CMz_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete [] MySurface_CL_Mnt; delete [] MySurface_CD_Mnt; delete [] MySurface_CSF_Mnt;
  delete [] MySurface_CEff_Mnt;  delete [] MySurface_CFx_Mnt;   delete [] MySurface_CFy_Mnt;
  delete [] MySurface_CFz_Mnt;
 delete [] MySurface_CMx_Mnt;   delete [] MySurface_CMy_Mnt;  delete [] MySurface_CMz_Mnt;

#endif
  
  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
  
  Total_CD         += AllBound_CD_Mnt;
  Total_CL        += AllBound_CL_Mnt;
  Total_CSF    += AllBound_CSF_Mnt;
  Total_CEff          = Total_CL / (Total_CD + EPS);
  Total_CMx           += AllBound_CMx_Mnt;
  Total_CMy           += AllBound_CMy_Mnt;
  Total_CMz          += AllBound_CMz_Mnt;
  Total_CFx           += AllBound_CFx_Mnt;
  Total_CFy           += AllBound_CFy_Mnt;
  Total_CFz           += AllBound_CFz_Mnt;
  Total_CT            += AllBound_CT_Mnt;
  Total_CQ            += AllBound_CQ_Mnt;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]      += Surface_CL_Mnt[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]      += Surface_CD_Mnt[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring] += Surface_CSF_Mnt[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       += Surface_CL_Mnt[iMarker_Monitoring] / (Surface_CD_Mnt[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        += Surface_CFx_Mnt[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        += Surface_CFy_Mnt[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        += Surface_CFz_Mnt[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        += Surface_CMx_Mnt[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        += Surface_CMy_Mnt[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]       += Surface_CMz_Mnt[iMarker_Monitoring];
  }
  
}

void CEulerSolver::TurboPerformance(CSolver *solver, CConfig *config, unsigned short inMarker,  unsigned short outMarker, unsigned short Kind_TurboPerf, unsigned short inMarkerTP ) {
  
  su2double  avgVel2In, avgVel2Out,avgVelRel2In, avgVelRel2Out, avgGridVel2In, avgGridVel2Out, avgTotalEnthalpyIn= 0.0,avgTotalRothalpyIn,
  avgTotalEnthalpyOut, avgTotalRothalpyOut, avgTotalEnthalpyOutIs, avgEnthalpyOut, avgEnthalpyOutIs,
  avgPressureOut, avgTotalRelPressureIn, avgTotalRelPressureOut, avgEntropyIn, avgEntropyOut;
  unsigned short iDim;
  
  
  /*--- compute or retrieve inlet information ---*/
  avgVelRel2In= 0.0;
  avgGridVel2In= 0.0;
  avgVel2In= 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    avgVelRel2In +=( AveragedVelocity[inMarker][iDim] - AveragedGridVel[inMarker][iDim])*( AveragedVelocity[inMarker][iDim] - AveragedGridVel[inMarker][iDim]);
    avgGridVel2In += AveragedGridVel[inMarker][iDim]*AveragedGridVel[inMarker][iDim];
    avgVel2In += AveragedVelocity[inMarker][iDim]*AveragedVelocity[inMarker][iDim];
  }
  
  avgTotalRothalpyIn = AveragedEnthalpy[inMarker] + 0.5*avgVelRel2In - 0.5*avgGridVel2In;
  avgTotalEnthalpyIn = AveragedEnthalpy[inMarker] + 0.5*avgVel2In;
  avgEntropyIn = AveragedEntropy[inMarker];
  FluidModel->SetTDState_hs(avgTotalRothalpyIn, avgEntropyIn);
  avgTotalRelPressureIn  = FluidModel->GetPressure();
  
  
  
  /*--- compute or retrieve outlet information ---*/
  avgVelRel2Out = 0.0;
  avgGridVel2Out = 0.0;
  avgVel2Out = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    avgVelRel2Out += (solver->GetAveragedVelocity(outMarker)[iDim]- solver->GetAveragedGridVelocity(outMarker)[iDim])*(solver->GetAveragedVelocity(outMarker)[iDim]- solver->GetAveragedGridVelocity(outMarker)[iDim]);
    avgGridVel2Out += solver->GetAveragedGridVelocity(outMarker)[iDim]*solver->GetAveragedGridVelocity(outMarker)[iDim];
    avgVel2Out += solver->GetAveragedVelocity(outMarker)[iDim]*solver->GetAveragedVelocity(outMarker)[iDim];
  }
  avgTotalRothalpyOut = solver->GetAveragedEnthalpy(outMarker) + 0.5*avgVelRel2Out - 0.5*avgGridVel2Out;
  avgTotalEnthalpyOut = solver->GetAveragedEnthalpy(outMarker) + 0.5*avgVel2Out;
  avgEntropyOut = solver->GetAveragedEntropy(outMarker);
  avgEnthalpyOut = solver->GetAveragedEnthalpy(outMarker);
  FluidModel->SetTDState_hs(avgTotalRothalpyOut, avgEntropyOut);
  avgTotalRelPressureOut  =  FluidModel->GetPressure();
  avgPressureOut= solver->GetAveragedPressure(outMarker);
  
  /*--- compute outlet isoentropic conditions ---*/
  FluidModel->SetTDState_Ps(avgPressureOut, avgEntropyIn);
  avgEnthalpyOutIs = FluidModel->GetStaticEnergy() + avgPressureOut/FluidModel->GetDensity();
  avgTotalEnthalpyOutIs = avgEnthalpyOutIs + 0.5*avgVel2Out;
  
  /*--- store turboperformance informations ---*/
  PressureOut[inMarkerTP] = avgPressureOut;
  PressureRatio[inMarkerTP] = avgTotalRelPressureIn/avgPressureOut;
  
  switch(Kind_TurboPerf) {
    case BLADE:
      
      TotalPressureLoss[inMarkerTP] = (avgTotalRelPressureIn - avgTotalRelPressureOut)/(avgTotalRelPressureOut - avgPressureOut) ;
      KineticEnergyLoss[inMarkerTP] = (avgEnthalpyOut - avgEnthalpyOutIs)/(avgTotalRothalpyIn - avgEnthalpyOut + 0.5*avgGridVel2Out);
      EulerianWork[inMarkerTP] = avgTotalEnthalpyIn - avgTotalEnthalpyOut;
      TotalEnthalpyIn[inMarkerTP] = avgTotalRothalpyIn;
      FlowAngleIn[inMarkerTP]= FlowAngle[inMarker];
      FlowAngleOut[inMarkerTP]= solver->GetFlowAngle(outMarker);
      MassFlowIn[inMarkerTP]= MassFlow[inMarker];
      MassFlowOut[inMarkerTP]= solver->GetMassFlow(outMarker);
      MachIn[inMarkerTP]= AveragedMach[inMarker];
      MachOut[inMarkerTP]= solver->GetAveragedMach(outMarker);
      NormalMachIn[inMarkerTP]= AveragedNormalMach[inMarker];
      NormalMachOut[inMarkerTP]= solver->GetAveragedNormalMach(outMarker);
      EnthalpyOut[inMarkerTP]= avgEnthalpyOut;
      VelocityOutIs[inMarkerTP]=sqrt(2.0*(avgTotalRothalpyIn - avgEnthalpyOut + 0.5*avgGridVel2Out));
      break;
      
    case STAGE: case TURBINE:
      
      TotalTotalEfficiency[inMarkerTP] = (avgTotalEnthalpyIn - avgTotalEnthalpyOut)/(avgTotalEnthalpyIn - avgTotalEnthalpyOutIs);
      TotalStaticEfficiency[inMarkerTP] = (avgTotalEnthalpyIn - avgTotalEnthalpyOut)/(avgTotalEnthalpyIn - avgEnthalpyOutIs);
      TotalEnthalpyIn[inMarkerTP]= avgTotalEnthalpyIn;
      EnthalpyOut[inMarkerTP] = avgTotalEnthalpyOut;
      break;
      
    default:
      cout << "Warning! Invalid TurboPerformance option!" << endl;
      exit(EXIT_FAILURE);
      break;
  }
  
  
  
  
  
}

void CEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  bool adjoint = config->GetContinuous_Adjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);
    
    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
    
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
  
}

void CEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  bool adjoint = config->GetContinuous_Adjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);
    
    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = local_Residual[iVar] + local_Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
    
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar, jVar;
  unsigned long iPoint, total_index, IterLinSol = 0;
  su2double Delta, *local_Res_TruncError, Vol;
  
  bool adjoint = config->GetContinuous_Adjoint();
  bool roe_turkel = config->GetKind_Upwind_Flow() == TURKEL;
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the residual ---*/
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    
    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    
    if (node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / node[iPoint]->GetDelta_Time();
      if (roe_turkel) {
        SetPreconditioner(config, iPoint);
        for (iVar = 0; iVar < nVar; iVar ++ )
          for (jVar = 0; jVar < nVar; jVar ++ )
            LowMach_Precontioner[iVar][jVar] = Delta*LowMach_Precontioner[iVar][jVar];
        Jacobian.AddBlock(iPoint, iPoint, LowMach_Precontioner);
      }
      else {
        Jacobian.AddVal2Diag(iPoint, Delta);
      }
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
  
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- The the number of iterations of the linear solver ---*/
  
  SetIterLinSolver(IterLinSol);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        node[iPoint]->AddSolution(iVar, config->GetRelaxation_Factor_Flow()*LinSysSol[iPoint*nVar+iVar]);
      }
    }
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CEulerSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  su2double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
  Partial_Gradient, Partial_Res, *Normal;
  
  /*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho) ---*/

  PrimVar_Vertex = new su2double [nPrimVarGrad];
  PrimVar_i = new su2double [nPrimVarGrad];
  PrimVar_j = new su2double [nPrimVarGrad];
  
  /*--- Set Gradient_Primitive to zero ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);

  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
      PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);
    }
    
    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = PrimVar_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
      }
    }
  }

  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()) {
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          PrimVar_Vertex[iVar] = node[iPoint]->GetPrimitive(iVar);
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++) {
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
            node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
          }
      }
    }
  }

  /*--- Update gradient value ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar, iDim) / (geometry->node[iPoint]->GetVolume());
        node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
      }
    }
  }

  delete [] PrimVar_Vertex;
  delete [] PrimVar_i;
  delete [] PrimVar_j;

  Set_MPI_Primitive_Gradient(geometry, config);

}

void CEulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular;
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Set the value of the singular ---*/
    singular = false;
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    
    /*--- Get primitives from CVariable ---*/
    
    PrimVar_i = node[iPoint]->GetPrimitive();
    
    /*--- Inizialization of variables ---*/
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
    
    AD::StartPreacc();
    AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    AD::SetPreaccIn(Coord_i, nDim);
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      PrimVar_j = node[jPoint]->GetPrimitive();
      
      AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (weight != 0.0) {
        
        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
        
        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
        
      }
      
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
    
    if (nDim == 3) {
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }
    
    /*--- Computation of the gradient: S*c ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        }
        
        node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
      }
    }
    
    AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
    AD::EndPreacc();
  }
  
  Set_MPI_Primitive_Gradient(geometry, config);
  
}

void CEulerSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j, *Primitive_i, *Primitive_j,
  dave, LimK, eps2, eps1, dm, dp, du, y, limiter;
  
  /*--- Initialize solution max and solution min and the limiter in the entire domain --*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      node[iPoint]->SetSolution_Max(iVar, -EPS);
      node[iPoint]->SetSolution_Min(iVar, EPS);
      node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
    }
  }
  
  /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Get the primitive variables ---*/
    
    Primitive_i = node[iPoint]->GetPrimitive();
    Primitive_j = node[jPoint]->GetPrimitive();
    
    /*--- Compute the maximum, and minimum values for nodes i & j ---*/
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      du = (Primitive_j[iVar] - Primitive_i[iVar]);
      node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
      node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
      node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
      node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
    }
    
  }
  
  
  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  
  if (config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
          else dp = node[iPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }
        
        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar))
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
          else dp = node[jPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }
        
        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar))
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
        
      }
      
    }
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        y =  node[iPoint]->GetLimiter_Primitive(iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
      }
    }
    
  }
  
  /*--- Venkatakrishnan limiter ---*/
  
  if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) {
    
    /*-- Get limiter parameters from the configuration file ---*/
    
    dave = config->GetRefElemLength();
    LimK = config->GetLimiterCoeff();
    eps1 = LimK*dave;
    eps2 = eps1*eps1*eps1;
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      

      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);


      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
        else dp = node[iPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }
      }

      AD::EndPreacc();
      
    }
    
  }
  
  /*--- Limiter MPI ---*/
  
  Set_MPI_Primitive_Limiter(geometry, config);
  
}

//void CEulerSolver::SetSecondary_Gradient_GG(CGeometry *geometry, CConfig *config) {
//  unsigned long iPoint, jPoint, iEdge, iVertex;
//  unsigned short iDim, iVar, iMarker;
//  su2double *SecondaryVar_Vertex, *SecondaryVar_i, *SecondaryVar_j, SecondaryVar_Average,
//  Partial_Gradient, Partial_Res, *Normal;
//
//  /*--- Gradient Secondary variables compressible (temp, vx, vy, vz, P, rho)
//   Gradient Secondary variables incompressible (rho, vx, vy, vz, beta) ---*/
//  SecondaryVar_Vertex = new su2double [nSecondaryVarGrad];
//  SecondaryVar_i = new su2double [nSecondaryVarGrad];
//  SecondaryVar_j = new su2double [nSecondaryVarGrad];
//
//  /*--- Set Gradient_Secondary to zero ---*/
//  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
//    node[iPoint]->SetGradient_SecondaryZero(nSecondaryVarGrad);
//
//  /*--- Loop interior edges ---*/
//  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//    iPoint = geometry->edge[iEdge]->GetNode(0);
//    jPoint = geometry->edge[iEdge]->GetNode(1);
//
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      SecondaryVar_i[iVar] = node[iPoint]->GetSecondary(iVar);
//      SecondaryVar_j[iVar] = node[jPoint]->GetSecondary(iVar);
//    }
//
//    Normal = geometry->edge[iEdge]->GetNormal();
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      SecondaryVar_Average =  0.5 * ( SecondaryVar_i[iVar] + SecondaryVar_j[iVar] );
//      for (iDim = 0; iDim < nDim; iDim++) {
//        Partial_Res = SecondaryVar_Average*Normal[iDim];
//        if (geometry->node[iPoint]->GetDomain())
//          node[iPoint]->AddGradient_Secondary(iVar, iDim, Partial_Res);
//        if (geometry->node[jPoint]->GetDomain())
//          node[jPoint]->SubtractGradient_Secondary(iVar, iDim, Partial_Res);
//      }
//    }
//  }
//
//  /*--- Loop boundary edges ---*/
//  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
//    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
//      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
//      if (geometry->node[iPoint]->GetDomain()) {
//
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          SecondaryVar_Vertex[iVar] = node[iPoint]->GetSecondary(iVar);
//
//        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++) {
//            Partial_Res = SecondaryVar_Vertex[iVar]*Normal[iDim];
//            node[iPoint]->SubtractGradient_Secondary(iVar, iDim, Partial_Res);
//          }
//      }
//    }
//  }
//
//  /*--- Update gradient value ---*/
//  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        Partial_Gradient = node[iPoint]->GetGradient_Secondary(iVar, iDim) / (geometry->node[iPoint]->GetVolume());
//        node[iPoint]->SetGradient_Secondary(iVar, iDim, Partial_Gradient);
//      }
//    }
//  }
//
//  delete [] SecondaryVar_Vertex;
//  delete [] SecondaryVar_i;
//  delete [] SecondaryVar_j;
//
//  Set_MPI_Secondary_Gradient(geometry, config);
//
//}

//void CEulerSolver::SetSecondary_Gradient_LS(CGeometry *geometry, CConfig *config) {
//
//  unsigned short iVar, iDim, jDim, iNeigh;
//  unsigned long iPoint, jPoint;
//  su2double *SecondaryVar_i, *SecondaryVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
//  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
//  bool singular;
//
//  /*--- Loop over points of the grid ---*/
//
//  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
//
//    /*--- Set the value of the singular ---*/
//    singular = false;
//
//    /*--- Get coordinates ---*/
//
//    Coord_i = geometry->node[iPoint]->GetCoord();
//
//    /*--- Get Secondarys from CVariable ---*/
//
//    SecondaryVar_i = node[iPoint]->GetSecondary();
//
//    /*--- Inizialization of variables ---*/
//
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//      for (iDim = 0; iDim < nDim; iDim++)
//        Cvector[iVar][iDim] = 0.0;
//
//    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
//    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0; detR2 = 0.0;
//
//    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
//      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
//      Coord_j = geometry->node[jPoint]->GetCoord();
//
//      SecondaryVar_j = node[jPoint]->GetSecondary();
//
//      weight = 0.0;
//      for (iDim = 0; iDim < nDim; iDim++)
//        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
//
//      /*--- Sumations for entries of upper triangular matrix R ---*/
//
//      if (weight != 0.0) {
//
//        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
//        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
//        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
//
//        if (nDim == 3) {
//          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
//          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
//          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
//          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
//        }
//
//        /*--- Entries of c:= transpose(A)*b ---*/
//
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(SecondaryVar_j[iVar]-SecondaryVar_i[iVar])/weight;
//
//      }
//
//    }
//
//    /*--- Entries of upper triangular matrix R ---*/
//
//    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
//    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
//    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
//
//    if (nDim == 3) {
//      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
//      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
//      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
//    }
//
//    /*--- Compute determinant ---*/
//
//    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
//    else detR2 = (r11*r22*r33)*(r11*r22*r33);
//
//    /*--- Detect singular matrices ---*/
//
//    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
//
//    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
//
//    if (singular) {
//      for (iDim = 0; iDim < nDim; iDim++)
//        for (jDim = 0; jDim < nDim; jDim++)
//          Smatrix[iDim][jDim] = 0.0;
//    }
//    else {
//      if (nDim == 2) {
//        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
//        Smatrix[0][1] = -r11*r12/detR2;
//        Smatrix[1][0] = Smatrix[0][1];
//        Smatrix[1][1] = r11*r11/detR2;
//      }
//      else {
//        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
//        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
//        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
//        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
//        Smatrix[0][2] = (z13*z33)/detR2;
//        Smatrix[1][0] = Smatrix[0][1];
//        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
//        Smatrix[1][2] = (z23*z33)/detR2;
//        Smatrix[2][0] = Smatrix[0][2];
//        Smatrix[2][1] = Smatrix[1][2];
//        Smatrix[2][2] = (z33*z33)/detR2;
//      }
//    }
//
//    /*--- Computation of the gradient: S*c ---*/
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        product = 0.0;
//        for (jDim = 0; jDim < nDim; jDim++) {
//          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
//        }
//
//        node[iPoint]->SetGradient_Secondary(iVar, iDim, product);
//      }
//    }
//
//  }
//
//  Set_MPI_Secondary_Gradient(geometry, config);
//
//}

//void CEulerSolver::SetSecondary_Limiter(CGeometry *geometry, CConfig *config) {
//
//  unsigned long iEdge, iPoint, jPoint;
//  unsigned short iVar, iDim;
//  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j, *Secondary_i, *Secondary_j,
//  dave, LimK, eps2, dm, dp, du, limiter;
//
//  /*--- Initialize solution max and solution min in the entire domain --*/
//  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      node[iPoint]->SetSolution_Max(iVar, -EPS);
//      node[iPoint]->SetSolution_Min(iVar, EPS);
//    }
//  }
//
//  /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
//  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//
//    /*--- Point identification, Normal vector and area ---*/
//    iPoint = geometry->edge[iEdge]->GetNode(0);
//    jPoint = geometry->edge[iEdge]->GetNode(1);
//
//    /*--- Get the conserved variables ---*/
//    Secondary_i = node[iPoint]->GetSecondary();
//    Secondary_j = node[jPoint]->GetSecondary();
//
//    /*--- Compute the maximum, and minimum values for nodes i & j ---*/
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      du = (Secondary_j[iVar] - Secondary_i[iVar]);
//      node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
//      node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
//      node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
//      node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
//    }
//  }
//
//  /*--- Initialize the limiter --*/
//  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      node[iPoint]->SetLimiter_Secondary(iVar, 2.0);
//    }
//  }
//
//  /*--- Venkatakrishnan limiter ---*/
//
//  if (config->GetKind_SlopeLimit() == VENKATAKRISHNAN) {
//
//    /*-- Get limiter parameters from the configuration file ---*/
//    dave = config->GetRefElemLength();
//    LimK = config->GetLimiterCoeff();
//    eps2 = pow((LimK*dave), 3.0);
//
//    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//
//      iPoint     = geometry->edge[iEdge]->GetNode(0);
//      jPoint     = geometry->edge[iEdge]->GetNode(1);
//      Gradient_i = node[iPoint]->GetGradient_Secondary();
//      Gradient_j = node[jPoint]->GetGradient_Secondary();
//      Coord_i    = geometry->node[iPoint]->GetCoord();
//      Coord_j    = geometry->node[jPoint]->GetCoord();
//
//      for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//
//        /*--- Calculate the interface left gradient, delta- (dm) ---*/
//        dm = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
//
//        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
//        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
//        else dp = node[iPoint]->GetSolution_Min(iVar);
//
//        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
//
//        if (limiter < node[iPoint]->GetLimiter_Secondary(iVar))
//          if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SetLimiter_Secondary(iVar, limiter);
//
//        /*-- Repeat for point j on the edge ---*/
//        dm = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
//
//        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
//        else dp = node[jPoint]->GetSolution_Min(iVar);
//
//        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
//
//        if (limiter < node[jPoint]->GetLimiter_Secondary(iVar))
//          if (geometry->node[jPoint]->GetDomain()) node[jPoint]->SetLimiter_Secondary(iVar, limiter);
//      }
//    }
//  }
//
//  /*--- Limiter MPI ---*/
//  Set_MPI_Secondary_Limiter(geometry, config);
//
//}

void CEulerSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) {
  unsigned short iDim, jDim, iVar, jVar;
  su2double Beta, local_Mach, Beta2, rho, enthalpy, soundspeed, sq_vel;
  su2double *U_i = NULL;
  su2double Beta_min = config->GetminTurkelBeta();
  su2double Beta_max = config->GetmaxTurkelBeta();
  
  
  /*--- Variables to calculate the preconditioner parameter Beta ---*/
  local_Mach = sqrt(node[iPoint]->GetVelocity2())/node[iPoint]->GetSoundSpeed();
  Beta          = max(Beta_min, min(local_Mach, Beta_max));
  Beta2             = Beta*Beta;
  
  U_i = node[iPoint]->GetSolution();
  
  rho = U_i[0];
  enthalpy = node[iPoint]->GetEnthalpy();
  soundspeed = node[iPoint]->GetSoundSpeed();
  sq_vel = node[iPoint]->GetVelocity2();
  
  /*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
  LowMach_Precontioner[0][0] = 0.5*sq_vel;
  LowMach_Precontioner[0][nVar-1] = 1.0;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[0][1+iDim] = -1.0*U_i[iDim+1]/rho;
  
  for (iDim = 0; iDim < nDim; iDim ++) {
    LowMach_Precontioner[iDim+1][0] = 0.5*sq_vel*U_i[iDim+1]/rho;
    LowMach_Precontioner[iDim+1][nVar-1] = U_i[iDim+1]/rho;
    for (jDim = 0; jDim < nDim; jDim ++) {
      LowMach_Precontioner[iDim+1][1+jDim] = -1.0*U_i[jDim+1]/rho*U_i[iDim+1]/rho;
    }
  }
  
  LowMach_Precontioner[nVar-1][0] = 0.5*sq_vel*enthalpy;
  LowMach_Precontioner[nVar-1][nVar-1] = enthalpy;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[nVar-1][1+iDim] = -1.0*U_i[iDim+1]/rho*enthalpy;
  
  
  for (iVar = 0; iVar < nVar; iVar ++ ) {
    for (jVar = 0; jVar < nVar; jVar ++ ) {
      LowMach_Precontioner[iVar][jVar] = (1.0/(Beta2+EPS) - 1.0) * (Gamma-1.0)/(soundspeed*soundspeed)*LowMach_Precontioner[iVar][jVar];
      if (iVar == jVar)
        LowMach_Precontioner[iVar][iVar] += 1.0;
    }
  }
  
}

void CEulerSolver::GetSurface_Properties(CGeometry *geometry, CNumerics *conv_numerics,
                                               CNumerics *visc_numerics, CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned short iDim, iMarker, iMarker_Analyze;
  unsigned long iVertex, iPoint;
  su2double Mach, Pressure, Temperature, TotalPressure, TotalTemperature, Velocity[3], Velocity2, MassFlow, Density, Energy, Area, AxiFactor = 1.0, SoundSpeed;
  
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Gamma        = config->GetGamma();

  bool axisymmetric               = config->GetAxisymmetric();
  unsigned short nMarker_Analyze    = config->GetnMarker_Analyze();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  su2double  *Surface_MassFlow    = new su2double[nMarker];
  su2double  *Surface_Mach        = new su2double[nMarker];
  su2double  *Surface_Pressure    = new su2double[nMarker];
  su2double  *Surface_Temperature = new su2double[nMarker];
  su2double  *Surface_TotalPressure    = new su2double[nMarker];
  su2double  *Surface_TotalTemperature = new su2double[nMarker];
  su2double  *Surface_Area        = new su2double[nMarker];
  
  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    Surface_MassFlow[iMarker] = 0.0;
    Surface_Mach[iMarker] = 0.0;
    Surface_Pressure[iMarker] = 0.0;
    Surface_Temperature[iMarker] = 0.0;
    Surface_TotalPressure[iMarker] = 0.0;
    Surface_TotalTemperature[iMarker] = 0.0;
    Surface_Area[iMarker] = 0.0;
    
    if (config->GetMarker_All_Analyze(iMarker) == YES) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          Density = node[iPoint]->GetSolution(0);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
            Velocity[iDim] = node[iPoint]->GetSolution(iDim+1) / Density;
            Velocity2 += Velocity[iDim] * Velocity[iDim];
            MassFlow += Vector[iDim] * AxiFactor * node[iPoint]->GetSolution(iDim+1);
          }
          
          Area              = sqrt (Area);
          Energy            = node[iPoint]->GetSolution(nVar-1)/Density;
          Pressure          = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
          SoundSpeed        = sqrt(Gamma*Pressure/Density);
          Mach              = sqrt(Velocity2)/SoundSpeed;
          Temperature       = Pressure / (Gas_Constant * Density);
          TotalPressure     = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma    / (Gamma - 1.0));
          TotalTemperature  = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));

          /*--- Compute the mass Surface_MassFlow ---*/
          
          Surface_MassFlow[iMarker] += MassFlow;
          Surface_Mach[iMarker] += Mach*MassFlow;
          Surface_Pressure[iMarker] += Pressure*MassFlow;
          Surface_Temperature[iMarker] += Temperature*MassFlow;
          Surface_TotalPressure[iMarker] += TotalPressure*MassFlow;
          Surface_TotalTemperature[iMarker] += TotalTemperature*MassFlow;
          Surface_Area[iMarker] += Area;
          
        }
      }
      
    }

  }
  
  /*--- Copy to the appropriate structure ---*/
  
  su2double *Surface_MassFlow_Local = new su2double [nMarker_Analyze];
  su2double *Surface_Mach_Local = new su2double [nMarker_Analyze];
  su2double *Surface_Temperature_Local = new su2double [nMarker_Analyze];
  su2double *Surface_Pressure_Local = new su2double [nMarker_Analyze];
  su2double *Surface_TotalTemperature_Local = new su2double [nMarker_Analyze];
  su2double *Surface_TotalPressure_Local = new su2double [nMarker_Analyze];
  su2double *Surface_Area_Local = new su2double [nMarker_Analyze];
  
  su2double *Surface_MassFlow_Total = new su2double [nMarker_Analyze];
  su2double *Surface_Mach_Total = new su2double [nMarker_Analyze];
  su2double *Surface_Temperature_Total = new su2double [nMarker_Analyze];
  su2double *Surface_Pressure_Total = new su2double [nMarker_Analyze];
  su2double *Surface_TotalTemperature_Total = new su2double [nMarker_Analyze];
  su2double *Surface_TotalPressure_Total = new su2double [nMarker_Analyze];
  su2double *Surface_Area_Total = new su2double [nMarker_Analyze];
  
  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    Surface_MassFlow_Local[iMarker_Analyze]    = 0.0;
    Surface_Mach_Local[iMarker_Analyze]        = 0.0;
    Surface_Temperature_Local[iMarker_Analyze] = 0.0;
    Surface_Pressure_Local[iMarker_Analyze]    = 0.0;
    Surface_TotalTemperature_Local[iMarker_Analyze] = 0.0;
    Surface_TotalPressure_Local[iMarker_Analyze]    = 0.0;
    Surface_Area_Local[iMarker_Analyze]        = 0.0;
    
    Surface_MassFlow_Total[iMarker_Analyze]    = 0.0;
    Surface_Mach_Total[iMarker_Analyze]        = 0.0;
    Surface_Temperature_Total[iMarker_Analyze] = 0.0;
    Surface_Pressure_Total[iMarker_Analyze]    = 0.0;
    Surface_TotalTemperature_Total[iMarker_Analyze] = 0.0;
    Surface_TotalPressure_Total[iMarker_Analyze]    = 0.0;
    Surface_Area_Total[iMarker_Analyze]        = 0.0;
  }
  
  /*--- Compute the numerical fan face Mach number, mach number, temperature and the total area ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_Analyze(iMarker) == YES)  {
      
      for (iMarker_Analyze= 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
        
        /*--- Add the Surface_MassFlow, and Surface_Area to the particular boundary ---*/
        
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_Analyze_TagBound(iMarker_Analyze)) {
          Surface_MassFlow_Local[iMarker_Analyze]    += Surface_MassFlow[iMarker];
          Surface_Mach_Local[iMarker_Analyze]        += Surface_Mach[iMarker];
          Surface_Temperature_Local[iMarker_Analyze] += Surface_Temperature[iMarker];
          Surface_Pressure_Local[iMarker_Analyze]    += Surface_Pressure[iMarker];
          Surface_TotalTemperature_Local[iMarker_Analyze] += Surface_TotalTemperature[iMarker];
          Surface_TotalPressure_Local[iMarker_Analyze]    += Surface_TotalPressure[iMarker];
          Surface_Area_Local[iMarker_Analyze]        += Surface_Area[iMarker];
        }
        
      }
      
    }
    
  }
  
#ifdef HAVE_MPI
  
  SU2_MPI::Allreduce(Surface_MassFlow_Local, Surface_MassFlow_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Mach_Local, Surface_Mach_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Temperature_Local, Surface_Temperature_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Pressure_Local, Surface_Pressure_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_TotalTemperature_Local, Surface_TotalTemperature_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_TotalPressure_Local, Surface_TotalPressure_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Area_Local, Surface_Area_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#else
  
  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    Surface_MassFlow_Total[iMarker_Analyze]    = Surface_MassFlow_Local[iMarker_Analyze];
    Surface_Mach_Total[iMarker_Analyze]        = Surface_Mach_Local[iMarker_Analyze];
    Surface_Temperature_Total[iMarker_Analyze] = Surface_Temperature_Local[iMarker_Analyze];
    Surface_Pressure_Total[iMarker_Analyze]    = Surface_Pressure_Local[iMarker_Analyze];
    Surface_TotalTemperature_Total[iMarker_Analyze] = Surface_TotalTemperature_Local[iMarker_Analyze];
    Surface_TotalPressure_Total[iMarker_Analyze]    = Surface_TotalPressure_Local[iMarker_Analyze];
    Surface_Area_Total[iMarker_Analyze]        = Surface_Area_Local[iMarker_Analyze];
  }
  
#endif
  
  /*--- Compute the value of Surface_Area_Total, and Surface_Pressure_Total, and
   set the value in the config structure for future use ---*/
  
  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    if (Surface_MassFlow_Total[iMarker_Analyze] != 0.0) Surface_Mach_Total[iMarker_Analyze] /= Surface_MassFlow_Total[iMarker_Analyze];
    else Surface_Mach_Total[iMarker_Analyze] = 0.0;
    if (Surface_MassFlow_Total[iMarker_Analyze] != 0.0) Surface_Temperature_Total[iMarker_Analyze] /= Surface_MassFlow_Total[iMarker_Analyze];
    else Surface_Temperature_Total[iMarker_Analyze] = 0.0;
    if (Surface_MassFlow_Total[iMarker_Analyze] != 0.0) Surface_Pressure_Total[iMarker_Analyze] /= Surface_MassFlow_Total[iMarker_Analyze];
    else Surface_Pressure_Total[iMarker_Analyze] = 0.0;
    if (Surface_MassFlow_Total[iMarker_Analyze] != 0.0) Surface_TotalTemperature_Total[iMarker_Analyze] /= Surface_MassFlow_Total[iMarker_Analyze];
    else Surface_TotalTemperature_Total[iMarker_Analyze] = 0.0;
    if (Surface_MassFlow_Total[iMarker_Analyze] != 0.0) Surface_TotalPressure_Total[iMarker_Analyze] /= Surface_MassFlow_Total[iMarker_Analyze];
    else Surface_TotalPressure_Total[iMarker_Analyze] = 0.0;

  }
  
  bool write_heads = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output && !config->GetDiscrete_Adjoint()) {
    
    cout.precision(4);
    cout.setf(ios::fixed, ios::floatfield);
    
    cout << endl << "--------------------------- Surface Analysis -------------------------" << endl;

    for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
      cout << "Surface ("<< config->GetMarker_Analyze_TagBound(iMarker_Analyze);

      su2double MassFlow = Surface_MassFlow_Total[iMarker_Analyze] * config->GetDensity_Ref() * config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == US) MassFlow *= 32.174;

      if (config->GetSystemMeasurements() == SI) cout << "): Mass flow (kg/s): " << MassFlow;
      else if (config->GetSystemMeasurements() == US) cout << "): Mass flow (lbs/s): " << MassFlow;
      config->SetSurface_MassFlow(iMarker_Analyze, MassFlow);

      cout << ", Mach: " << Surface_Mach_Total[iMarker_Analyze];
      if (config->GetSystemMeasurements() == SI) cout << ", T(K): ";
      else if (config->GetSystemMeasurements() == US) cout << ", T(R): ";
      cout << Surface_Temperature_Total[iMarker_Analyze] * config->GetTemperature_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", P(Pa): ";
      else if (config->GetSystemMeasurements() == US) cout << ", P(psf): ";
      cout << Surface_Pressure_Total[iMarker_Analyze] * config->GetPressure_Ref();

      if (config->GetSystemMeasurements() == SI) cout << ", TT(K): ";
       else if (config->GetSystemMeasurements() == US) cout << ", TT(R): ";
       cout << Surface_TotalTemperature_Total[iMarker_Analyze] * config->GetTemperature_Ref();
       if (config->GetSystemMeasurements() == SI) cout << ", PT(Pa): ";
       else if (config->GetSystemMeasurements() == US) cout << ", PT(psf): ";
       cout << Surface_TotalPressure_Total[iMarker_Analyze] * config->GetPressure_Ref();

      if (config->GetSystemMeasurements() == SI) cout << ", Area (m^2): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Area (ft^2): ";
      cout << Surface_Area_Total[iMarker_Analyze] <<"."<< endl;
    }
    
  }
  
  delete [] Surface_MassFlow_Local;
  delete [] Surface_Mach_Local;
  delete [] Surface_Temperature_Local;
  delete [] Surface_Pressure_Local;
  delete [] Surface_TotalTemperature_Local;
  delete [] Surface_TotalPressure_Local;
  delete [] Surface_Area_Local;
  
  delete [] Surface_MassFlow_Total;
  delete [] Surface_Mach_Total;
  delete [] Surface_Temperature_Total;
  delete [] Surface_Pressure_Total;
  delete [] Surface_TotalTemperature_Total;
  delete [] Surface_TotalPressure_Total;
  delete [] Surface_Area_Total;
  
  delete [] Surface_MassFlow;
  delete [] Surface_Mach;
  delete [] Surface_Temperature;
  delete [] Surface_Pressure;
  delete [] Surface_TotalTemperature;
  delete [] Surface_TotalPressure;
  delete [] Surface_Area;
  
}

void CEulerSolver::GetSurface_Distortion(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned short iMarker, iDim;
  unsigned long iPoint, iVertex;
  su2double xCoord = 0.0, yCoord = 0.0, zCoord = 0.0, q = 0.0, PT = 0.0, Area = 0.0, *Vector, TotalArea = 0.0;
  su2double xCoord_CG = 0.0, yCoord_CG = 0.0, zCoord_CG = 0.0, TipRadius, HubRadius, Distance = 0.0, Distance_Mirror = 0.0;
  su2double *r, MinDistance, xCoord_ = 0.0, yCoord_ = 0.0, zCoord_ = 0, dx = 0.0, dy = 0.0, dz = 0.0, dx_ = 0.0, dy_ = 0.0, dz_ = 0.0;
  unsigned short iStation, iAngle, nAngle, Theta, nStation;
  su2double *** ProbeArray, *Mach_Station, *Mach_Station_Min, *PT_Sector, *PT_Station, *PT_Station_Min,
  PT_Sector_Min, PT_Mean, Mach_Mean, q_Mean, UpVector[3], radians, RotatedVector[3],
//  RefDensity, RefVel, Factor,
  DC60, IDR, IDC, IDC_Mach;
  su2double Pressure, Density, SoundSpeed, Velocity2, Mach,  Gamma, TotalPressure, Mach_Inf, TotalPressure_Inf, Pressure_Inf;
//  su2double dMach_dVel_x = 0.0, dMach_dVel_y = 0.0, dMach_dVel_z = 0.0, dMach_dT = 0.0
//  su2double dMach_dx = 0.0, dMach_dy = 0.0, dMach_dz = 0.0, dPT_dP = 0.0, dPT_dMach = 0.0, Aux = 0.0;
  
  unsigned short iMarker_Analyze;
  int rank, iProcessor, nProcessor;
  rank = MASTER_NODE;
  nProcessor = SINGLE_NODE;
  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
  unsigned long nVertex_Surface, nLocalVertex_Surface, MaxLocalVertex_Surface;
  unsigned long Total_Index;
  unsigned short nDim = geometry->GetnDim();
  unsigned short Theta_DC60 = 60, nStation_DC60 = 5;
  
  bool Evaluate_Distortion = ((((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0)) || (config->GetExtIter() == 1) || (config->GetDiscrete_Adjoint()));
  bool write_heads         = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));
  bool Engine_HalfModel    = config->GetEngine_HalfModel();
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
  
  if (Evaluate_Distortion) {
    
    for (iMarker_Analyze = 0; iMarker_Analyze < config->GetnMarker_Analyze(); iMarker_Analyze++) {
      
      string Analyze_TagBound = config->GetMarker_Analyze_TagBound(iMarker_Analyze);
      
      nVertex_Surface = 0, nLocalVertex_Surface = 0; MaxLocalVertex_Surface = 0;
      
      /*--- Find the max number of surface vertices among all
       partitions and set up buffers. The master node will handle the
       writing of the CSV file after gathering all of the data. ---*/
      
      nLocalVertex_Surface = 0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        string Marker_TagBound = config->GetMarker_All_TagBound(iMarker);
        if (Marker_TagBound == Analyze_TagBound) {
          for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
          }
        }
      }
      
      /*--- Communicate the number of local vertices on each partition
       to the master node ---*/
      
      Buffer_Send_nVertex[0] = nLocalVertex_Surface;
      if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];
      
#ifdef HAVE_MPI
      SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
      SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
      MaxLocalVertex_Surface = nLocalVertex_Surface;
      Buffer_Recv_nVertex[MASTER_NODE] = Buffer_Send_nVertex[MASTER_NODE];
#endif
      
      /*--- Send and Recv buffers ---*/
      
      su2double *Buffer_Send_Coord_x = NULL, *Buffer_Recv_Coord_x = NULL;
      Buffer_Send_Coord_x = new su2double [MaxLocalVertex_Surface];
      
      su2double *Buffer_Send_Coord_y = NULL, *Buffer_Recv_Coord_y = NULL;
      Buffer_Send_Coord_y = new su2double [MaxLocalVertex_Surface];
      
      su2double *Buffer_Send_Coord_z = NULL, *Buffer_Recv_Coord_z = NULL;
      if (nDim == 3)  Buffer_Send_Coord_z = new su2double [MaxLocalVertex_Surface];
      
      su2double *Buffer_Send_PT = NULL, *Buffer_Recv_PT = NULL;
      Buffer_Send_PT = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dPT_dx = NULL, *Buffer_Recv_dPT_dx = NULL;
      //      Buffer_Send_dPT_dx = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dPT_dy = NULL, *Buffer_Recv_dPT_dy = NULL;
      //      Buffer_Send_dPT_dy = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dPT_dz = NULL, *Buffer_Recv_dPT_dz = NULL;
      //      if (nDim == 3) Buffer_Send_dPT_dz = new su2double [MaxLocalVertex_Surface];
      
      su2double *Buffer_Send_Mach = NULL, *Buffer_Recv_Mach = NULL;
      Buffer_Send_Mach = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dMach_dx = NULL, *Buffer_Recv_dMach_dx = NULL;
      //      Buffer_Send_dMach_dx = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dMach_dy = NULL, *Buffer_Recv_dMach_dy = NULL;
      //      Buffer_Send_dMach_dy = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dMach_dz = NULL, *Buffer_Recv_dMach_dz = NULL;
      //      if (nDim == 3) Buffer_Send_dMach_dz = new su2double [MaxLocalVertex_Surface];
      
      su2double *Buffer_Send_q = NULL, *Buffer_Recv_q = NULL;
      Buffer_Send_q = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dq_dx = NULL, *Buffer_Recv_dq_dx = NULL;
      //      Buffer_Send_dq_dx = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dq_dy = NULL, *Buffer_Recv_dq_dy = NULL;
      //      Buffer_Send_dq_dy = new su2double [MaxLocalVertex_Surface];
      
      //      su2double *Buffer_Send_dq_dz = NULL, *Buffer_Recv_dq_dz = NULL;
      //      if (nDim == 3) Buffer_Send_dq_dz = new su2double [MaxLocalVertex_Surface];
      
      su2double *Buffer_Send_Area = NULL, *Buffer_Recv_Area = NULL;
      Buffer_Send_Area = new su2double [MaxLocalVertex_Surface];

      /*--- Prepare the receive buffers on the master node only. ---*/
      
      if (rank == MASTER_NODE) {
        Buffer_Recv_Coord_x = new su2double [nProcessor*MaxLocalVertex_Surface];
        Buffer_Recv_Coord_y = new su2double [nProcessor*MaxLocalVertex_Surface];
        if (nDim == 3) Buffer_Recv_Coord_z = new su2double [nProcessor*MaxLocalVertex_Surface];
        Buffer_Recv_PT = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        Buffer_Recv_dPT_dx = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        Buffer_Recv_dPT_dy = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        if (nDim == 3) Buffer_Recv_dPT_dz = new su2double [nProcessor*MaxLocalVertex_Surface];
        Buffer_Recv_Mach = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        Buffer_Recv_dMach_dx = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        Buffer_Recv_dMach_dy = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        if (nDim == 3) Buffer_Recv_dMach_dz = new su2double [nProcessor*MaxLocalVertex_Surface];
        Buffer_Recv_q = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        Buffer_Recv_dq_dx = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        Buffer_Recv_dq_dy = new su2double [nProcessor*MaxLocalVertex_Surface];
        //        if (nDim == 3) Buffer_Recv_dq_dz = new su2double [nProcessor*MaxLocalVertex_Surface];
        Buffer_Recv_Area = new su2double [nProcessor*MaxLocalVertex_Surface];
      }
      
      /*--- Loop over all vertices in this partition and load the
       data of the specified type into the buffer to be sent to
       the master node. ---*/
      
      nVertex_Surface = 0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        string Marker_TagBound = config->GetMarker_All_TagBound(iMarker);
        if (Marker_TagBound == Analyze_TagBound) {
          
          for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            
            if (geometry->node[iPoint]->GetDomain()) {
              
              Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
              Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
              if (nDim == 3) { Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2); }
              
              Pressure     = node[iPoint]->GetPressure();
              Density      = node[iPoint]->GetDensity();
              SoundSpeed   = node[iPoint]->GetSoundSpeed();
              Velocity2    = node[iPoint]->GetVelocity2();
              Mach         = sqrt(Velocity2)/SoundSpeed;
              Gamma        = config->GetGamma();
              Mach_Inf     = config->GetMach();
              Pressure_Inf = config->GetPressure_FreeStreamND();
              
              Buffer_Send_Mach[nVertex_Surface] = Mach;
              //              dMach_dVel_x = 0.0; dMach_dVel_y = 0.0; dMach_dVel_z = 0.0;
              //              if ((Velocity2 != 0.0) && (Mach != 0.0)) {
              //                  dMach_dVel_x = node[iPoint]->GetVelocity(0) / (Mach * sqrt(Velocity2));
              //                  dMach_dVel_y = node[iPoint]->GetVelocity(1) / (Mach * sqrt(Velocity2));
              //                  if (nDim == 3) { dMach_dVel_z = node[iPoint]->GetVelocity(2) / (Mach * sqrt(Velocity2)); }
              //              }
              //              Aux = Gas_Constant*Temperature;
              //              dMach_dT = - Gas_Constant * sqrt(Velocity2) / (2.0 * sqrt(Gamma) * pow(Aux, 1.5));
              
              //              dMach_dx = dMach_dT*node[iPoint]->GetGradient_Primitive(0, 0) +
              //              dMach_dVel_x*node[iPoint]->GetGradient_Primitive(1, 0) + dMach_dVel_y*node[iPoint]->GetGradient_Primitive(2, 0);
              //              if (nDim == 3) { Buffer_Send_dMach_dx[nVertex_Surface] += dMach_dVel_z*node[iPoint]->GetGradient_Primitive(3, 0); }
              //              Buffer_Send_dMach_dx[nVertex_Surface] = dMach_dx;
              
              //              dMach_dy = dMach_dT*node[iPoint]->GetGradient_Primitive(0, 1) +
              //              dMach_dVel_x*node[iPoint]->GetGradient_Primitive(1, 1) + dMach_dVel_y*node[iPoint]->GetGradient_Primitive(2, 1);
              //              if (nDim == 3) { Buffer_Send_dMach_dx[nVertex_Surface] += dMach_dVel_z*node[iPoint]->GetGradient_Primitive(3, 1); }
              //              Buffer_Send_dMach_dy[nVertex_Surface] = dMach_dy;
              
              //              if (nDim == 3) {
              //                dMach_dz = dMach_dT*node[iPoint]->GetGradient_Primitive(0, 2) +
              //                dMach_dVel_x*node[iPoint]->GetGradient_Primitive(1, 2) + dMach_dVel_y*node[iPoint]->GetGradient_Primitive(2, 2) +
              //                dMach_dVel_z*node[iPoint]->GetGradient_Primitive(3, 2);
              //                Buffer_Send_dMach_dz[nVertex_Surface] = dMach_dz;
              //              }
              
              TotalPressure      = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
              TotalPressure_Inf  = Pressure_Inf * pow( 1.0 + Mach_Inf * Mach_Inf * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
              
              Buffer_Send_PT[nVertex_Surface] = TotalPressure / TotalPressure_Inf;
              //              dPT_dP = (pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0))) / TotalPressure_Inf;
              //              dPT_dMach = (Gamma * Mach * Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), 1.0 / (Gamma - 1.0))) / TotalPressure_Inf;
              //              Buffer_Send_dPT_dx[nVertex_Surface] = dPT_dP*node[iPoint]->GetGradient_Primitive(nDim+1, 0) + dPT_dMach * dMach_dx;
              //              Buffer_Send_dPT_dy[nVertex_Surface] = dPT_dP*node[iPoint]->GetGradient_Primitive(nDim+1, 1) + dPT_dMach * dMach_dy;
              //              if (nDim == 3) { Buffer_Send_dPT_dz[nVertex_Surface] = dPT_dP*node[iPoint]->GetGradient_Primitive(nDim+1, 2) + dPT_dMach * dMach_dz; }
              
              Buffer_Send_q[nVertex_Surface] = 0.5*Density*Velocity2;
              //              Buffer_Send_dq_dx[nVertex_Surface] = 0.0;
              //              Buffer_Send_dq_dy[nVertex_Surface] = 0.0;
              //              if (nDim == 3) { Buffer_Send_dq_dz[nVertex_Surface] = 0.0; }
              
              Vector = geometry->vertex[iMarker][iVertex]->GetNormal();
              Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) { Area += Vector[iDim]*Vector[iDim]; } Area = sqrt(Area);
              Buffer_Send_Area[nVertex_Surface] = Area;
              
              /*--- If US system, the output should be in inches ---*/
              
              if (config->GetSystemMeasurements() == US) {
                Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
                Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
                if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
                Buffer_Send_Area[nVertex_Surface] *= 144.0;
                
                //                Buffer_Send_dPT_dx[nVertex_Surface] /= 12.0;
                //                Buffer_Send_dPT_dy[nVertex_Surface] /= 12.0;
                //                if (nDim == 3) Buffer_Send_dPT_dz[nVertex_Surface] /= 12.0;
                
                //                Buffer_Send_dMach_dx[nVertex_Surface] /= 12.0;
                //                Buffer_Send_dMach_dy[nVertex_Surface] /= 12.0;
                //                if (nDim == 3) Buffer_Send_dMach_dz[nVertex_Surface] /= 12.0;
                
                //                Buffer_Send_dq_dx[nVertex_Surface] /= 12.0;
                //                Buffer_Send_dq_dy[nVertex_Surface] /= 12.0;
                //                if (nDim == 3) Buffer_Send_dq_dz[nVertex_Surface] /= 12.0;
              }
                            
              nVertex_Surface++;
              
            }
          }
        }
      }
      
      /*--- Send the information to the master node ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_PT, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_PT, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      SU2_MPI::Gather(Buffer_Send_dPT_dx, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dPT_dx, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      SU2_MPI::Gather(Buffer_Send_dPT_dy, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dPT_dy, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      if (nDim == 3) SU2_MPI::Gather(Buffer_Send_dPT_dz, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dPT_dz, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      SU2_MPI::Gather(Buffer_Send_dMach_dx, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dMach_dx, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      SU2_MPI::Gather(Buffer_Send_dMach_dy, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dMach_dy, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      if (nDim == 3) SU2_MPI::Gather(Buffer_Send_dMach_dz, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dMach_dz, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_q, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_q, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      SU2_MPI::Gather(Buffer_Send_dq_dx, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dq_dx, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      SU2_MPI::Gather(Buffer_Send_dq_dy, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dq_dy, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      //      if (nDim == 3) SU2_MPI::Gather(Buffer_Send_dq_dz, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_dq_dz, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Area, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Area, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iVertex = 0; iVertex < MaxLocalVertex_Surface; iVertex++) {
        Buffer_Recv_Coord_x[iVertex] = Buffer_Send_Coord_x[iVertex];
        Buffer_Recv_Coord_y[iVertex] = Buffer_Send_Coord_y[iVertex];
        if (nDim == 3) Buffer_Recv_Coord_z[iVertex] = Buffer_Send_Coord_z[iVertex];
        Buffer_Recv_PT[iVertex] = Buffer_Send_PT[iVertex];
        //        Buffer_Recv_dPT_dx[iVertex] = Buffer_Send_dPT_dx[iVertex];
        //        Buffer_Recv_dPT_dy[iVertex] = Buffer_Send_dPT_dy[iVertex];
        //        if (nDim == 3) Buffer_Recv_dPT_dz[iVertex] = Buffer_Send_dPT_dz[iVertex];
        Buffer_Recv_Mach[iVertex] = Buffer_Send_Mach[iVertex];
        //        Buffer_Recv_dMach_dx[iVertex] = Buffer_Send_dMach_dx[iVertex];
        //        Buffer_Recv_dMach_dy[iVertex] = Buffer_Send_dMach_dy[iVertex];
        //        if (nDim == 3) Buffer_Recv_dMach_dz[iVertex] = Buffer_Send_dMach_dz[iVertex];
        Buffer_Recv_q[iVertex] = Buffer_Send_q[iVertex];
        //        Buffer_Recv_dq_dx[iVertex] = Buffer_Send_dq_dx[iVertex];
        //        Buffer_Recv_dq_dy[iVertex] = Buffer_Send_dq_dy[iVertex];
        //        if (nDim == 3) Buffer_Recv_dq_dz[iVertex] = Buffer_Send_dq_dz[iVertex];
        Buffer_Recv_Area[iVertex] = Buffer_Send_Area[iVertex];
      }
#endif
      
      if (rank == MASTER_NODE) {
        
        /*--- Compute center of gravity ---*/
        
        TotalArea = 0.0; xCoord_CG = 0.0; yCoord_CG = 0.0; zCoord_CG = 0.0; PT_Mean = 0.0; Mach_Mean = 0.0;  q_Mean = 0.0;
        
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {
            
            /*--- Current index position and global index ---*/
            
            Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;
            
            /*--- Retrieve the merged data for this node ---*/
            
            xCoord = Buffer_Recv_Coord_x[Total_Index];
            yCoord = Buffer_Recv_Coord_y[Total_Index];
            if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
            PT   = Buffer_Recv_PT[Total_Index];
            Mach = Buffer_Recv_Mach[Total_Index];
            q    = Buffer_Recv_q[Total_Index];
            
            Area       = Buffer_Recv_Area[Total_Index];
            TotalArea += Area;
            xCoord_CG += xCoord*Area;
            yCoord_CG += yCoord*Area;
            zCoord_CG += zCoord*Area;
            PT_Mean   += PT*Area;
            Mach_Mean += PT*Area;
            q_Mean    += q*Area;
            
          }
        }
        
        /*--- Evaluate the area averaged pressure and CG ---*/
        
        xCoord_CG /= TotalArea;
        yCoord_CG /= TotalArea;
        zCoord_CG /= TotalArea;
        PT_Mean   /= TotalArea;
        Mach_Mean /= TotalArea;
        q_Mean    /=  TotalArea;
        
        /*--- If it is a half model, CGy = 0 ---*/
        
        if (Engine_HalfModel) { yCoord_CG = 0.0; }
        
        /*--- Compute hub and tip radius ---*/
        
        TipRadius = 1E-6; HubRadius = 1E6;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {
            
            /*--- Current index position and global index ---*/
            
            Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;
            
            /*--- Retrieve the merged data for this node ---*/
            
            xCoord = Buffer_Recv_Coord_x[Total_Index];
            yCoord = Buffer_Recv_Coord_y[Total_Index];
            if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
            
            if (nDim == 2)
              Distance = sqrt((xCoord_CG-xCoord)*(xCoord_CG-xCoord) +
                              (yCoord_CG-yCoord)*(yCoord_CG-yCoord));
            
            if (nDim == 3)
              Distance = sqrt((xCoord_CG-xCoord)*(xCoord_CG-xCoord) +
                              (yCoord_CG-yCoord)*(yCoord_CG-yCoord) +
                              (zCoord_CG-zCoord)*(zCoord_CG-zCoord));
            
            if (Distance > TipRadius) TipRadius = Distance;
            if (Distance < HubRadius) HubRadius = Distance;
            
          }
        }
        
        if (HubRadius/TipRadius < 0.05) HubRadius = 0.0;
        
        /*--- Evaluate the DC60 parameter ---*/
        
        Theta = Theta_DC60;
        nStation = nStation_DC60;
        
        nAngle = SU2_TYPE::Int(360/float(Theta));
        r = new su2double [nStation+1];
        
        /*--- Allocate memory ---*/
        
        PT_Sector = new su2double [nAngle];
        ProbeArray = new su2double ** [nAngle];
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          ProbeArray[iAngle] = new su2double * [nStation];
          for (iStation = 0; iStation < nStation; iStation++) {
            ProbeArray[iAngle][iStation] = new su2double [5];
          }
        }
        
        /*--- Define the radius for each probe ---*/
        
        r[0] = HubRadius;
        r[nStation] = TipRadius;
        
        for (iStation = 1; iStation < nStation; iStation++) {
          r[iStation] = sqrt(  r[iStation-1]*r[iStation-1] + (r[nStation]*r[nStation] - r[0]*r[0])/float(nStation) );
        }
        
        /*--- Define the probe rack ---*/
        
        UpVector[0] = 0.0; UpVector[1] = 0.0; UpVector[2] = 1.0;
        
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          
          radians = -iAngle*Theta*2.0*PI_NUMBER/360;
          RotatedVector[0] =  UpVector[0];
          RotatedVector[1] =  UpVector[1] * cos(radians) - UpVector[2] * sin(radians);
          RotatedVector[2] =  UpVector[1] * sin(radians) + UpVector[2] * cos(radians);
          
          for (iStation = 1; iStation <= nStation; iStation++) {
            ProbeArray[iAngle][iStation-1][0] = xCoord_CG+RotatedVector[0]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
            ProbeArray[iAngle][iStation-1][1] = yCoord_CG+RotatedVector[1]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
            ProbeArray[iAngle][iStation-1][2] = zCoord_CG+RotatedVector[2]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
          }
          
        }
        
        /*--- Compute the Total pressure at each probe, closes grid point to the location ---*/
        
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          for (iStation = 0; iStation < nStation; iStation++) {
            xCoord_ = ProbeArray[iAngle][iStation][0];
            yCoord_ = ProbeArray[iAngle][iStation][1];
            zCoord_ = ProbeArray[iAngle][iStation][2];
            
            MinDistance = 1E6;
            
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
              for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {
                
                /*--- Current index position and global index ---*/
                
                Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;
                
                /*--- Retrieve the merged data for this node ---*/
                
                xCoord = Buffer_Recv_Coord_x[Total_Index];
                yCoord = Buffer_Recv_Coord_y[Total_Index];
                if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
                
                dx = (xCoord_ - xCoord);
                dy = (yCoord_ - yCoord);
                if (nDim == 3) dz = (zCoord_ - zCoord);
                
                Distance = dx*dx + dy*dy;
                if (nDim == 3) Distance += dz*dz;
                Distance = sqrt(Distance);
                                
                if (Engine_HalfModel) {
                  
                  yCoord = -yCoord;
                  
                  dx_ = (xCoord_ - xCoord);
                  dy_ = (yCoord_ - yCoord);
                  if (nDim == 3) dz_ = (zCoord_ - zCoord);
                  
                  Distance_Mirror = dx_*dx_ + dy_*dy_;
                  if (nDim == 3) Distance_Mirror += dz_*dz_;
                  Distance_Mirror = sqrt(Distance_Mirror);
                  
                  if (Distance_Mirror < Distance) {
                    Distance = Distance_Mirror;
                    dx = dx_; dy = dy_;
                    if (nDim == 3) dz = dz_;
                  }
                  
                }
                
                if (Distance <= MinDistance) {
                  MinDistance = Distance;
                  ProbeArray[iAngle][iStation][3] = Buffer_Recv_PT[Total_Index]; //+ Buffer_Recv_dPT_dx[Total_Index]*dx + Buffer_Recv_dPT_dy[Total_Index]*dy*SignFlip;
                  //                  if (nDim == 3) ProbeArray[iAngle][iStation][3] += Buffer_Recv_dPT_dz[Total_Index]*dz;
                  
                  ProbeArray[iAngle][iStation][4] = Buffer_Recv_q[Total_Index]; //  + Buffer_Recv_dq_dx[Total_Index]*dx + Buffer_Recv_dq_dy[Total_Index]*dy*SignFlip;
                  //                  if (nDim == 3) ProbeArray[iAngle][iStation][4] += Buffer_Recv_dq_dz[Total_Index]*dz;
                  
                }
                
              }
            }
            
          }
          
        }
        
        /*--- Evaluate the average pressure at each sector ---*/
        
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          PT_Sector[iAngle] = 0.0;
          for (iStation = 0; iStation < nStation; iStation++) {
            PT_Sector[iAngle] += ProbeArray[iAngle][iStation][3]/float(nStation);
          }
        }
        
        /*--- Compute the total average pressure at the fan face ---*/
        
        PT_Mean = 0.0;
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          for (iStation = 0; iStation < nStation; iStation++) {
            PT_Mean += ProbeArray[iAngle][iStation][3]/float(nStation*nAngle);
          }
        }
        
        /*--- Compute the total average dynamic pressure at the fan face ---*/
        
        q_Mean = 0.0;
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          for (iStation = 0; iStation < nStation; iStation++) {
            q_Mean += ProbeArray[iAngle][iStation][4]/float(nStation*nAngle);
          }
        }
        
        /*--- Compute the min value of the averaged pressure at each sector ---*/
        
        PT_Sector_Min = PT_Sector[0];
        for (iAngle = 1; iAngle < nAngle; iAngle++) {
          if (PT_Sector[iAngle] <= PT_Sector_Min) PT_Sector_Min = PT_Sector[iAngle];
        }

        /*--- Set the value of the distortion, it only works for one surface ---*/
        
        Mach_Inf           = config->GetMach();
        Gamma              = config->GetGamma();
        TotalPressure_Inf  = config->GetPressure_FreeStreamND() * pow( 1.0 + Mach_Inf * Mach_Inf *
                                                                      0.5 * (Gamma - 1.0), Gamma    / (Gamma - 1.0));
        
        if (q_Mean != 0.0) DC60 = ((PT_Mean - PT_Sector_Min)*TotalPressure_Inf)/q_Mean;
        else DC60 = 0.0;
        
        config->SetSurface_DC60(iMarker_Analyze, DC60);
        
        SetTotal_DC60(DC60);
        
        /*--- Deallocate the memory ---*/
        
        delete[] r;
        
        delete [] PT_Sector;
        
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          for (iStation = 0; iStation < nStation; iStation++) {
            delete[] ProbeArray[iAngle][iStation];
          }
        }
        delete[] ProbeArray;
        
        
        /*--- Evaluate the IDC, and IDR parameters ---*/
        
        nStation = SU2_TYPE::Int(config->GetDistortionRack()[0]);
        Theta = SU2_TYPE::Int(config->GetDistortionRack()[1]);
        
        nAngle = SU2_TYPE::Int(360/float(Theta));
        
        
        /*--- Allocate memory ---*/
        
        r = new su2double [nStation+1];
        ProbeArray = new su2double ** [nAngle];
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          ProbeArray[iAngle] = new su2double * [nStation];
          for (iStation = 0; iStation < nStation; iStation++) {
            ProbeArray[iAngle][iStation] = new su2double [4];
          }
        }
        
        /*--- Define the radius for each probe ---*/
        
        r[0] = HubRadius;
        r[nStation] = TipRadius;
        
        for (iStation = 1; iStation < nStation; iStation++) {
          r[iStation] = sqrt(  r[iStation-1]*r[iStation-1] + (r[nStation]*r[nStation] - r[0]*r[0])/float(nStation) );
        }
        
        /*--- Define the probe rack ---*/
        
        UpVector[0] = 0.0; UpVector[1] = 0.0; UpVector[2] = 1.0;
        
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          
          radians = -iAngle*Theta*2.0*PI_NUMBER/360;
          RotatedVector[0] =  UpVector[0];
          RotatedVector[1] =  UpVector[1] * cos(radians) - UpVector[2] * sin(radians);
          RotatedVector[2] =  UpVector[1] * sin(radians) + UpVector[2] * cos(radians);
          
          for (iStation = 1; iStation <= nStation; iStation++) {
            ProbeArray[iAngle][iStation-1][0] = xCoord_CG+RotatedVector[0]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
            ProbeArray[iAngle][iStation-1][1] = yCoord_CG+RotatedVector[1]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
            ProbeArray[iAngle][iStation-1][2] = zCoord_CG+RotatedVector[2]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
          }
          
        }
        
        /*--- Compute the Total pressure at each probe, closes grid point to the location ---*/
        
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          for (iStation = 0; iStation < nStation; iStation++) {
            xCoord_ = ProbeArray[iAngle][iStation][0];
            yCoord_ = ProbeArray[iAngle][iStation][1];
            zCoord_ = ProbeArray[iAngle][iStation][2];
            
            MinDistance = 1E6;
            
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
              for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {
                
                /*--- Current index position and global index ---*/
                
                Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;
                
                /*--- Retrieve the merged data for this node ---*/
                
                xCoord = Buffer_Recv_Coord_x[Total_Index];
                yCoord = Buffer_Recv_Coord_y[Total_Index];
                if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
                
                dx = (xCoord_ - xCoord);
                dy = (yCoord_ - yCoord);
                if (nDim == 3) dz = (zCoord_ - zCoord);
                
                Distance = dx*dx + dy*dy;
                if (nDim == 3) Distance += dz*dz;
                Distance = sqrt(Distance);
                                
                if (Engine_HalfModel) {
                  
                  yCoord = -yCoord;
                  
                  dx_ = (xCoord_ - xCoord);
                  dy_ = (yCoord_ - yCoord);
                  if (nDim == 3) dz_ = (zCoord_ - zCoord);
                  
                  Distance_Mirror = dx_*dx_ + dy_*dy_;
                  if (nDim == 3) Distance_Mirror += dz_*dz_;
                  Distance_Mirror = sqrt(Distance_Mirror);
                  
                  if (Distance_Mirror < Distance) {
                    Distance = Distance_Mirror;
                    dx = dx_; dy = dy_;
                    if (nDim == 3) dz = dz_;
                  }
                  
                }
                
                if (Distance <= MinDistance) {
                  MinDistance = Distance;
                  ProbeArray[iAngle][iStation][3] = Buffer_Recv_PT[Total_Index]; // + Buffer_Recv_dPT_dx[Total_Index]*dx + Buffer_Recv_dPT_dy[Total_Index]*dy*SignFlip;
                  //                  if (nDim == 3) ProbeArray[iAngle][iStation][3] += Buffer_Recv_dPT_dz[Total_Index]*dz;
                }
                
              }
            }
            
          }
          
        }
        
        /*--- Evaluate the average and min. pressure at each station/radius  ---*/
        
        PT_Station = new su2double [nStation];
        PT_Station_Min = new su2double [nStation];
        
        for (iStation = 0; iStation < nStation; iStation++) {
          PT_Station[iStation] = 0.0;
          PT_Station_Min[iStation] = ProbeArray[0][iStation][3];
          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            PT = ProbeArray[iAngle][iStation][3];
            PT_Station[iStation] += PT / float(nAngle);
            if (PT <= PT_Station_Min[iStation] ) PT_Station_Min[iStation] = PT;
          }
        }
        
        /*--- Compute the total average pressure at the fan face ---*/
        
        PT_Mean = 0.0;
        for (iStation = 0; iStation < nStation; iStation++) {
          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            PT_Mean += ProbeArray[iAngle][iStation][3]/float(nStation*nAngle);
          }
        }
        
        /*--- Set the value of the distortion, it only works for one surface ---*/
        
        IDC = 0.0;
        for (iStation = 0; iStation < nStation-1; iStation++) {
          IDC = max (IDC, 0.5*((PT_Station[iStation] -PT_Station_Min[iStation])/PT_Mean
                               + (PT_Station[iStation+1] -PT_Station_Min[iStation+1])/PT_Mean)   );
          
        }
        
        config->SetSurface_IDC(iMarker_Analyze, IDC);
        
        SetTotal_IDC(IDC);
        SetTotal_CircumferentialDistortion(IDC);
        
        IDR = 0.0;
        for (iStation = 0; iStation < nStation; iStation++) {
          IDR = max (IDR, (PT_Mean-PT_Station[iStation])/PT_Mean);
        }
        
        config->SetSurface_IDR(iMarker_Analyze, IDR);
        
        SetTotal_IDR(IDR);
        SetTotal_RadialDistortion(IDR);
        
        /*--- Release IDX parameters ---*/
        
        delete [] PT_Station_Min;
        delete [] PT_Station;
        
        /*--- Evaluate the IDC Mach parameter ---*/
        
        /*--- Compute the Mach number at each probe, closes grid point to the location ---*/
        
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          for (iStation = 0; iStation < nStation; iStation++) {
            xCoord_ = ProbeArray[iAngle][iStation][0];
            yCoord_ = ProbeArray[iAngle][iStation][1];
            zCoord_ = ProbeArray[iAngle][iStation][2];
            
            MinDistance = 1E6;
            
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
              for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {
                
                /*--- Current index position and global index ---*/
                
                Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;
                
                /*--- Retrieve the merged data for this node ---*/
                
                xCoord = Buffer_Recv_Coord_x[Total_Index];
                yCoord = Buffer_Recv_Coord_y[Total_Index];
                if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
                
                dx = (xCoord_ - xCoord);
                dy = (yCoord_ - yCoord);
                if (nDim == 3) dz = (zCoord_ - zCoord);
                
                Distance = dx*dx + dy*dy;
                if (nDim == 3) Distance += dz*dz;
                Distance = sqrt(Distance);
                                
                if (Engine_HalfModel) {
                  
                  yCoord = -yCoord;
                  
                  dx_ = (xCoord_ - xCoord);
                  dy_ = (yCoord_ - yCoord);
                  if (nDim == 3) dz_ = (zCoord_ - zCoord);
                  
                  Distance_Mirror = dx_*dx_ + dy_*dy_;
                  if (nDim == 3) Distance_Mirror += dz_*dz_;
                  Distance_Mirror = sqrt(Distance_Mirror);
                  
                  if (Distance_Mirror < Distance) {
                    Distance = Distance_Mirror;
                    dx = dx_; dy = dy_;
                    if (nDim == 3) dz = dz_;
                  }
                  
                }
                
                if (Distance <= MinDistance) {
                  MinDistance = Distance;
                  ProbeArray[iAngle][iStation][3] = Buffer_Recv_Mach[Total_Index]; // + Buffer_Recv_dMach_dx[Total_Index]*dx + Buffer_Recv_dMach_dy[Total_Index]*dy*SignFlip;
                  //                  if (nDim == 3) ProbeArray[iAngle][iStation][3] += Buffer_Recv_dMach_dz[Total_Index]*dz;
                }
                
              }
            }
            
          }
          
        }
        
        /*--- Evaluate the average and min. pressure at each station/radius  ---*/
        
        Mach_Station = new su2double [nStation];
        Mach_Station_Min = new su2double [nStation];
        
        for (iStation = 0; iStation < nStation; iStation++) {
          Mach_Station[iStation] = 0.0;
          Mach_Station_Min[iStation] = ProbeArray[0][iStation][3];
          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            Mach = ProbeArray[iAngle][iStation][3];
            Mach_Station[iStation] += Mach / float(nAngle);
            if (Mach <= Mach_Station_Min[iStation] ) Mach_Station_Min[iStation] = Mach;
          }
        }
        
        /*--- Compute the total average pressure at the fan face ---*/
        
        Mach_Mean = 0.0;
        for (iStation = 0; iStation < nStation; iStation++) {
          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            Mach_Mean += ProbeArray[iAngle][iStation][3]/float(nStation*nAngle);
          }
        }
        
        /*--- Set the value of the distortion, it only works for one surface ---*/
        
        IDC_Mach = 0.0;
        for (iStation = 0; iStation < nStation-1; iStation++) {
          if (Mach_Mean!=0)
            IDC_Mach = max (IDC_Mach, 0.5*((Mach_Station[iStation] - Mach_Station_Min[iStation])/Mach_Mean
                                           + (Mach_Station[iStation+1] - Mach_Station_Min[iStation+1])/Mach_Mean)   );
          
        }
        
        config->SetSurface_IDC_Mach(iMarker_Analyze, IDC_Mach);
        
        SetTotal_IDC_Mach(IDC_Mach);
        
        delete [] Mach_Station_Min;
        delete [] Mach_Station;
        
        /*--- Release distortion parameters ---*/
        
        delete[] r;
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          for (iStation = 0; iStation < nStation; iStation++) {
            delete[] ProbeArray[iAngle][iStation];
          }
        }
        delete[] ProbeArray;
        
        /*--- Release the recv buffers on the master node ---*/
        
        delete [] Buffer_Recv_Coord_x;
        delete [] Buffer_Recv_Coord_y;
        if (nDim == 3) delete [] Buffer_Recv_Coord_z;
        
        delete [] Buffer_Recv_PT;
        //        delete [] Buffer_Recv_dPT_dx;
        //        delete [] Buffer_Recv_dPT_dy;
        //        if (nDim == 3) delete [] Buffer_Recv_dPT_dz;
        
        delete [] Buffer_Recv_Mach;
        //        delete [] Buffer_Recv_dMach_dx;
        //        delete [] Buffer_Recv_dMach_dy;
        //        if (nDim == 3) delete [] Buffer_Recv_dMach_dz;
        
        delete [] Buffer_Recv_q;
        //        delete [] Buffer_Recv_dq_dx;
        //        delete [] Buffer_Recv_dq_dy;
        //        if (nDim == 3) delete [] Buffer_Recv_dq_dz;
        
        delete [] Buffer_Recv_Area;

        delete [] Buffer_Recv_nVertex;
        
        
      }
      
      if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output && !config->GetDiscrete_Adjoint()) {
        
        cout << "Surface ("<< Analyze_TagBound << "): ";
        cout.precision(4);
        cout.setf(ios::fixed, ios::floatfield);
        cout << setprecision(1) << "Dist. coeff.: IDC " << 100*config->GetSurface_IDC(iMarker_Analyze)
        << "%. IDC Mach " << 100*config->GetSurface_IDC_Mach(iMarker_Analyze)
        << "%. IDR " << 100*config->GetSurface_IDR(iMarker_Analyze)
        << "%. DC60 " << config->GetSurface_DC60(iMarker_Analyze) << "." << endl;
        
      }
      
      /*--- Release the memory for the remaining buffers and exit ---*/
      
      delete [] Buffer_Send_Coord_x;
      delete [] Buffer_Send_Coord_y;
      if (nDim == 3) delete [] Buffer_Send_Coord_z;
      
      delete [] Buffer_Send_PT;
      //      delete [] Buffer_Send_dPT_dx;
      //      delete [] Buffer_Send_dPT_dy;
      //      if (nDim == 3) delete [] Buffer_Send_dPT_dz;
      
      delete [] Buffer_Send_q;
      //      delete [] Buffer_Send_dq_dx;
      //      delete [] Buffer_Send_dq_dy;
      //      if (nDim == 3) delete [] Buffer_Send_dq_dz;
      
      delete [] Buffer_Send_Mach;
      //      delete [] Buffer_Send_dMach_dx;
      //      delete [] Buffer_Send_dMach_dy;
      //      if (nDim == 3) delete [] Buffer_Send_dMach_dz;
      
      delete [] Buffer_Send_Area;
            
    }
    
  }
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output && !config->GetDiscrete_Adjoint()) {
    
    cout << "-------------------------------------------------------------------------" << endl;
    
  }
  
}

void CEulerSolver::GetPower_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned short iDim, iMarker, jMarker;
  unsigned long iVertex, iPoint;
  su2double  *V_inlet = NULL, *V_outlet = NULL, Pressure, Temperature, Velocity[3], Vn,
  Velocity2, Density, Area, SoundSpeed, TotalPressure, Vel_Infty2, RamDrag,
  TotalTemperature, VelocityJet,
  Vel_Infty, MaxPressure, MinPressure, MFR, InfVel2;
  unsigned short iMarker_Inlet, iMarker_Outlet, nMarker_Inlet, nMarker_Outlet;
  string Inlet_TagBound, Outlet_TagBound;
  su2double DeltaPress = 0.0, DeltaTemp = 0.0, TotalPressRatio = 0.0, TotalTempRatio = 0.0, StaticPressRatio = 0.0, StaticTempRatio = 0.0,
  NetThrust = 0.0, GrossThrust = 0.0, Power = 0.0, MassFlow = 0.0, Mach = 0.0, Force = 0.0;
  bool ReverseFlow, Engine = false, Pair = true;
  
  su2double Gas_Constant                         = config->GetGas_ConstantND();
  su2double Cp                                               = Gas_Constant*Gamma / (Gamma-1.0);
  su2double Alpha                                         = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta                                           = config->GetAoS()*PI_NUMBER/180.0;
  bool write_heads = ((((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0) && (config->GetExtIter()!= 0))
                      || (config->GetExtIter() == 1));
  bool Evaluate_BC = ((((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0))
                      || (config->GetExtIter() == 1) || (config->GetDiscrete_Adjoint()));
  
  if ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0)) Engine = true;
  if ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0)) Engine = false;
  if ((config->GetnMarker_EngineInflow()) != (config->GetnMarker_EngineExhaust())) Pair = false;
  
  
  if (Engine) { nMarker_Inlet  = config->GetnMarker_EngineInflow(); nMarker_Outlet = config->GetnMarker_EngineExhaust(); }
  else  { nMarker_Inlet   = config->GetnMarker_ActDiskInlet(); nMarker_Outlet  = config->GetnMarker_ActDiskOutlet(); }
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Evaluate the MPI for the actuator disk IO ---*/
  
  if (Evaluate_BC) {
    
    /*--- Allocate memory ---*/
    
    su2double *Inlet_MassFlow         = new su2double [config->GetnMarker_All()];
    su2double *Inlet_ReverseMassFlow  = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Pressure         = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Mach             = new su2double [config->GetnMarker_All()];
    su2double *Inlet_MaxPressure      = new su2double [config->GetnMarker_All()];
    su2double *Inlet_MinPressure      = new su2double [config->GetnMarker_All()];
    su2double *Inlet_TotalPressure    = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Temperature      = new su2double [config->GetnMarker_All()];
    su2double *Inlet_TotalTemperature   = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Area             = new su2double [config->GetnMarker_All()];
    su2double *Inlet_RamDrag          = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Force            = new su2double [config->GetnMarker_All()];
    su2double *Inlet_Power            = new su2double [config->GetnMarker_All()];
    su2double *Inlet_XCG                        = new su2double [config->GetnMarker_All()];
    su2double *Inlet_YCG                        = new su2double [config->GetnMarker_All()];
    su2double *Inlet_ZCG                        = new su2double [config->GetnMarker_All()];
    
    su2double *Outlet_MassFlow           = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Pressure           = new su2double [config->GetnMarker_All()];
    su2double *Outlet_TotalPressure    = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Temperature        = new su2double [config->GetnMarker_All()];
    su2double *Outlet_TotalTemperature = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Area                   = new su2double [config->GetnMarker_All()];
    su2double *Outlet_GrossThrust        = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Force                    = new su2double [config->GetnMarker_All()];
    su2double *Outlet_Power                    = new su2double [config->GetnMarker_All()];
    
    /*--- Comute MassFlow, average temp, press, etc. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      Inlet_MassFlow[iMarker] = 0.0;      Inlet_ReverseMassFlow[iMarker] = 0.0;  Inlet_MinPressure[iMarker] = 0.0;
      Inlet_Pressure[iMarker] = 0.0;      Inlet_Mach[iMarker] = 0.0;             Inlet_Temperature[iMarker] = 0.0;
      Inlet_MinPressure[iMarker] = +1E10; Inlet_MaxPressure[iMarker] = -1E10;    Inlet_Power[iMarker] = 0.0;
      Inlet_TotalPressure[iMarker] = 0.0; Inlet_TotalTemperature[iMarker] = 0.0;
      Inlet_Area[iMarker] = 0.0;
      Inlet_RamDrag[iMarker] = 0.0; Inlet_Force[iMarker]    = 0.0;
      Inlet_XCG[iMarker] = 0.0;  Inlet_YCG[iMarker]    = 0.0; Inlet_ZCG[iMarker]    = 0.0;
      
      Outlet_MassFlow[iMarker] = 0.0;
      Outlet_Pressure[iMarker] = 0.0; Outlet_Temperature[iMarker] = 0.0;
      Outlet_TotalPressure[iMarker] = 0.0; Outlet_TotalTemperature[iMarker] = 0.0;
      Outlet_Area[iMarker] = 0.0;
      Outlet_GrossThrust[iMarker]    = 0.0; Outlet_Force[iMarker]    = 0.0; Outlet_Power[iMarker]    = 0.0;
      
      MinPressure = +1E10; MaxPressure = -1E10;
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            V_inlet = node[iPoint]->GetPrimitive();
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
            
            Temperature     = V_inlet[0];
            Pressure                = V_inlet[nDim+1];
            
            Density                         = V_inlet[nDim+2];
            SoundSpeed  = sqrt(Gamma*Pressure/Density);
            
            Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vel_Infty2 =0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Area += Vector[iDim]*Vector[iDim];
              Velocity[iDim] = V_inlet[iDim+1];
              Velocity2 += Velocity[iDim]*Velocity[iDim];
              Vel_Infty2 += GetVelocity_Inf(iDim)*GetVelocity_Inf(iDim);
              MassFlow -= Vector[iDim]*Velocity[iDim]*Density;
            }
            
            Vn = 0.0; ReverseFlow = false;
            for (iDim = 0; iDim < nDim; iDim++) {  Vn -= Velocity[iDim]*Vector[iDim]/Area; }
            if (Vn < 0.0) { ReverseFlow = true; }
            
            Vel_Infty                   = sqrt (Vel_Infty2);
            Area                                = sqrt (Area);
            Mach                                = sqrt(Velocity2)/SoundSpeed;
            TotalPressure           = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma    / (Gamma - 1.0));
            TotalTemperature  = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));
            MinPressure               = min(MinPressure, TotalPressure);
            MaxPressure               = max(MaxPressure, TotalPressure);
            
            RamDrag     = MassFlow * Vel_Infty;
            
            Inlet_MassFlow[iMarker]         += MassFlow;
            Inlet_Pressure[iMarker]         += Pressure*MassFlow;
            Inlet_Mach[iMarker]             += Mach*MassFlow;
            Inlet_MinPressure[iMarker]       = min (MinPressure, Inlet_MinPressure[iMarker]);
            Inlet_MaxPressure[iMarker]       = max(MaxPressure, Inlet_MaxPressure[iMarker]);
            Inlet_TotalPressure[iMarker]    += TotalPressure*MassFlow;
            Inlet_Temperature[iMarker]      += Temperature*MassFlow;
            Inlet_TotalTemperature[iMarker] += TotalTemperature*MassFlow;
            Inlet_Area[iMarker]             += Area;
            Inlet_RamDrag[iMarker]          += RamDrag;
            Inlet_Power[iMarker] += MassFlow*Cp*TotalTemperature;
            if (ReverseFlow) Inlet_ReverseMassFlow[iMarker]  += MassFlow;
            
            su2double Inlet_ForceX =  -(Pressure - Pressure_Inf)*Vector[0] + MassFlow*Velocity[0];
            su2double Inlet_ForceY = -(Pressure - Pressure_Inf)*Vector[1] + MassFlow*Velocity[1];
            su2double Inlet_ForceZ = 0.0;
            if (nDim == 3) Inlet_ForceZ = -(Pressure - Pressure_Inf)*Vector[2] + MassFlow*Velocity[2];
            Inlet_Force[iMarker] +=  Inlet_ForceX*cos(Alpha)*cos(Beta) + Inlet_ForceY*sin(Beta) +Inlet_ForceZ*sin(Alpha)*cos(Beta);
            
            Inlet_XCG[iMarker]                        += geometry->node[iPoint]->GetCoord(0)*Area;
            Inlet_YCG[iMarker]                        += geometry->node[iPoint]->GetCoord(1)*Area;
            if (nDim == 3) Inlet_ZCG[iMarker] += geometry->node[iPoint]->GetCoord(2)*Area;
            
          }
        }
        
      }
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            V_outlet = node[iPoint]->GetPrimitive();
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
            
            Temperature         = V_outlet[0];
            Pressure                    = V_outlet[nDim+1];
            
            Density                         = V_outlet[nDim+2];
            SoundSpeed  = sqrt(Gamma*Pressure/Density);
            
            Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vel_Infty2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Area += Vector[iDim]*Vector[iDim];
              Velocity[iDim] = V_outlet[iDim+1];
              Velocity2 += Velocity[iDim]*Velocity[iDim];
              Vel_Infty2 += GetVelocity_Inf(iDim)*GetVelocity_Inf(iDim);
              MassFlow += Vector[iDim]*Velocity[iDim]*Density;
            }
            
            Vel_Infty                   = sqrt (Vel_Infty2);
            Area                                            = sqrt (Area);
            Mach                                = sqrt(Velocity2)/SoundSpeed;
            TotalPressure           = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma    / (Gamma - 1.0));
            TotalTemperature    = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));
            VelocityJet         = sqrt(Velocity2) ;

            GrossThrust     = MassFlow * VelocityJet;

            Outlet_MassFlow[iMarker]                            += MassFlow;
            Outlet_Pressure[iMarker]                                += Pressure*MassFlow;
            Outlet_TotalPressure[iMarker]               += TotalPressure*MassFlow;
            Outlet_Temperature[iMarker]                     += Temperature*MassFlow;
            Outlet_TotalTemperature[iMarker]    += TotalTemperature*MassFlow;
            Outlet_Area[iMarker]                                            += Area;
            Outlet_GrossThrust[iMarker]             += GrossThrust;
            Outlet_Power[iMarker] += MassFlow*Cp*TotalTemperature;
            
            su2double Outlet_ForceX = -(Pressure - Pressure_Inf)*Vector[0] -MassFlow*Velocity[0];
            su2double Outlet_ForceY =  -(Pressure - Pressure_Inf)*Vector[1] -MassFlow*Velocity[1];
            su2double Outlet_ForceZ = 0.0;
            if (nDim == 3) Outlet_ForceZ = -(Pressure - Pressure_Inf)*Vector[2] -MassFlow*Velocity[2];
            
            if (nDim == 2) Outlet_Force[iMarker] +=  Outlet_ForceX*cos(Alpha) + Outlet_ForceY*sin(Alpha);
            if (nDim == 3) Outlet_Force[iMarker] +=  Outlet_ForceX*cos(Alpha)*cos(Beta) + Outlet_ForceY*sin(Beta) + Outlet_ForceZ*sin(Alpha)*cos(Beta);
            
          }
        }
        
      }
      
    }
    
    /*--- Copy to the appropriate structure ---*/
    
    su2double *Inlet_MassFlow_Local             = new su2double [nMarker_Inlet];
    su2double *Inlet_ReverseMassFlow_Local              = new su2double [nMarker_Inlet];
    su2double *Inlet_Temperature_Local          = new su2double [nMarker_Inlet];
    su2double *Inlet_TotalTemperature_Local = new su2double [nMarker_Inlet];
    su2double *Inlet_Pressure_Local             = new su2double [nMarker_Inlet];
    su2double *Inlet_Mach_Local             = new su2double [nMarker_Inlet];
    su2double *Inlet_MinPressure_Local      = new su2double [nMarker_Inlet];
    su2double *Inlet_MaxPressure_Local      = new su2double [nMarker_Inlet];
    su2double *Inlet_Power_Local        = new su2double [nMarker_Inlet];
    su2double *Inlet_TotalPressure_Local        = new su2double [nMarker_Inlet];
    su2double *Inlet_RamDrag_Local              = new su2double [nMarker_Inlet];
    su2double *Inlet_Force_Local                            = new su2double [nMarker_Inlet];
    su2double *Inlet_Area_Local              = new su2double [nMarker_Inlet];
    su2double *Inlet_XCG_Local                      = new su2double [nMarker_Inlet];
    su2double *Inlet_YCG_Local                      = new su2double [nMarker_Inlet];
    su2double *Inlet_ZCG_Local                      = new su2double [nMarker_Inlet];
    
    su2double *Inlet_MassFlow_Total                         = new su2double [nMarker_Inlet];
    su2double *Inlet_ReverseMassFlow_Total              = new su2double [nMarker_Inlet];
    su2double *Inlet_Pressure_Total                             = new su2double [nMarker_Inlet];
    su2double *Inlet_Mach_Total                             = new su2double [nMarker_Inlet];
    su2double *Inlet_MinPressure_Total            = new su2double [nMarker_Inlet];
    su2double *Inlet_MaxPressure_Total         = new su2double [nMarker_Inlet];
    su2double *Inlet_Power_Total                            = new su2double [nMarker_Inlet];
    su2double *Inlet_TotalPressure_Total            = new su2double [nMarker_Inlet];
    su2double *Inlet_Temperature_Total                  = new su2double [nMarker_Inlet];
    su2double *Inlet_TotalTemperature_Total     = new su2double [nMarker_Inlet];
    su2double *Inlet_RamDrag_Total                          = new su2double [nMarker_Inlet];
    su2double *Inlet_Force_Total                                = new su2double [nMarker_Inlet];
    su2double *Inlet_Area_Total                                         = new su2double [nMarker_Inlet];
    su2double *Inlet_XCG_Total                                          = new su2double [nMarker_Inlet];
    su2double *Inlet_YCG_Total                                          = new su2double [nMarker_Inlet];
    su2double *Inlet_ZCG_Total                                          = new su2double [nMarker_Inlet];
    
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      Inlet_MassFlow_Local[iMarker_Inlet]                       = 0.0;
      Inlet_ReverseMassFlow_Local[iMarker_Inlet]            = 0.0;
      Inlet_Pressure_Local[iMarker_Inlet]                           = 0.0;
      Inlet_Mach_Local[iMarker_Inlet]                           = 0.0;
      Inlet_MinPressure_Local[iMarker_Inlet]            = 0.0;
      Inlet_MaxPressure_Local[iMarker_Inlet]            = 0.0;
      Inlet_TotalPressure_Local[iMarker_Inlet]          = 0.0;
      Inlet_Temperature_Local[iMarker_Inlet]                    = 0.0;
      Inlet_TotalTemperature_Local[iMarker_Inlet] = 0.0;
      Inlet_RamDrag_Local[iMarker_Inlet]                    = 0.0;
      Inlet_Force_Local[iMarker_Inlet]                              = 0.0;
      Inlet_Power_Local[iMarker_Inlet]                          = 0.0;
      Inlet_Area_Local[iMarker_Inlet]                                           = 0.0;
      Inlet_XCG_Local[iMarker_Inlet]                                            = 0.0;
      Inlet_YCG_Local[iMarker_Inlet]                                            = 0.0;
      Inlet_ZCG_Local[iMarker_Inlet]                                            = 0.0;
      
      Inlet_MassFlow_Total[iMarker_Inlet]                           = 0.0;
      Inlet_ReverseMassFlow_Total[iMarker_Inlet]                = 0.0;
      Inlet_Pressure_Total[iMarker_Inlet]                               = 0.0;
      Inlet_Mach_Total[iMarker_Inlet]                               = 0.0;
      Inlet_MinPressure_Total[iMarker_Inlet]                = 0.0;
      Inlet_MaxPressure_Total[iMarker_Inlet]                = 0.0;
      Inlet_TotalPressure_Total[iMarker_Inlet]              = 0.0;
      Inlet_Temperature_Total[iMarker_Inlet]                    = 0.0;
      Inlet_TotalTemperature_Total[iMarker_Inlet]   = 0.0;
      Inlet_RamDrag_Total[iMarker_Inlet]                        = 0.0;
      Inlet_Force_Total[iMarker_Inlet]                                  = 0.0;
      Inlet_Power_Total[iMarker_Inlet]              = 0.0;
      Inlet_Area_Total[iMarker_Inlet]                                           = 0.0;
      Inlet_XCG_Total[iMarker_Inlet]                                            = 0.0;
      Inlet_YCG_Total[iMarker_Inlet]                                            = 0.0;
      Inlet_ZCG_Total[iMarker_Inlet]                                            = 0.0;
    }
    
    su2double *Outlet_MassFlow_Local                            = new su2double [nMarker_Outlet];
    su2double *Outlet_Pressure_Local                                = new su2double [nMarker_Outlet];
    su2double *Outlet_TotalPressure_Local           = new su2double [nMarker_Outlet];
    su2double *Outlet_Temperature_Local                     = new su2double [nMarker_Outlet];
    su2double *Outlet_TotalTemperature_Local    = new su2double [nMarker_Outlet];
    su2double *Outlet_GrossThrust_Local             = new su2double [nMarker_Outlet];
    su2double *Outlet_Force_Local               = new su2double [nMarker_Outlet];
    su2double *Outlet_Power_Local               = new su2double [nMarker_Outlet];
    su2double *Outlet_Area_Local                                            = new su2double [nMarker_Outlet];
    
    su2double *Outlet_MassFlow_Total                        = new su2double [nMarker_Outlet];
    su2double *Outlet_Pressure_Total                            = new su2double [nMarker_Outlet];
    su2double *Outlet_TotalPressure_Total           = new su2double [nMarker_Outlet];
    su2double *Outlet_Temperature_Total                     = new su2double [nMarker_Outlet];
    su2double *Outlet_TotalTemperature_Total    = new su2double [nMarker_Outlet];
    su2double *Outlet_GrossThrust_Total             = new su2double [nMarker_Outlet];
    su2double *Outlet_Force_Total                             = new su2double [nMarker_Outlet];
    su2double *Outlet_Power_Total                             = new su2double [nMarker_Outlet];
    su2double *Outlet_Area_Total                                            = new su2double [nMarker_Outlet];
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      Outlet_MassFlow_Local[iMarker_Outlet]                     = 0.0;
      Outlet_Pressure_Local[iMarker_Outlet]                         = 0.0;
      Outlet_TotalPressure_Local[iMarker_Outlet]            = 0.0;
      Outlet_Temperature_Local[iMarker_Outlet]                  = 0.0;
      Outlet_TotalTemperature_Local[iMarker_Outlet] = 0.0;
      Outlet_GrossThrust_Local[iMarker_Outlet]              = 0.0;
      Outlet_Force_Local[iMarker_Outlet]                        = 0.0;
      Outlet_Power_Local[iMarker_Outlet]                        = 0.0;
      Outlet_Area_Local[iMarker_Outlet]                                         = 0.0;
      
      Outlet_MassFlow_Total[iMarker_Outlet]                             = 0.0;
      Outlet_Pressure_Total[iMarker_Outlet]                             = 0.0;
      Outlet_TotalPressure_Total[iMarker_Outlet]                = 0.0;
      Outlet_Temperature_Total[iMarker_Outlet]                  = 0.0;
      Outlet_TotalTemperature_Total[iMarker_Outlet] = 0.0;
      Outlet_GrossThrust_Total[iMarker_Outlet]              = 0.0;
      Outlet_Force_Total[iMarker_Outlet]                        = 0.0;
      Outlet_Power_Total[iMarker_Outlet]                        = 0.0;
      Outlet_Area_Total[iMarker_Outlet]                                         = 0.0;
    }
    
    /*--- Copy the values to the local array for MPI ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||  (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) {
        for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
          
          if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
          if (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW) Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
          
          if (config->GetMarker_All_TagBound(iMarker) == Inlet_TagBound) {
            Inlet_MassFlow_Local[iMarker_Inlet]                             += Inlet_MassFlow[iMarker];
            Inlet_ReverseMassFlow_Local[iMarker_Inlet]              += Inlet_ReverseMassFlow[iMarker];
            Inlet_Pressure_Local[iMarker_Inlet]                                 += Inlet_Pressure[iMarker];
            Inlet_Mach_Local[iMarker_Inlet]                                     += Inlet_Mach[iMarker];
            Inlet_MinPressure_Local[iMarker_Inlet]              += Inlet_MinPressure[iMarker];
            Inlet_MaxPressure_Local[iMarker_Inlet]          += Inlet_MaxPressure[iMarker];
            Inlet_TotalPressure_Local[iMarker_Inlet]            += Inlet_TotalPressure[iMarker];
            Inlet_Temperature_Local[iMarker_Inlet]                  += Inlet_Temperature[iMarker];
            Inlet_TotalTemperature_Local[iMarker_Inlet] += Inlet_TotalTemperature[iMarker];
            Inlet_RamDrag_Local[iMarker_Inlet]                          += Inlet_RamDrag[iMarker];
            Inlet_Force_Local[iMarker_Inlet]                                    += Inlet_Force[iMarker];
            Inlet_Power_Local[iMarker_Inlet]                          += Inlet_Power[iMarker];
            Inlet_Area_Local[iMarker_Inlet]                                                 += Inlet_Area[iMarker];
            Inlet_XCG_Local[iMarker_Inlet]                                              += Inlet_XCG[iMarker];
            Inlet_YCG_Local[iMarker_Inlet]                                              += Inlet_YCG[iMarker];
            if (nDim == 3) Inlet_ZCG_Local[iMarker_Inlet]                                               += Inlet_ZCG[iMarker];
          }
          
        }
      }
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) || (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST)) {
        for (iMarker_Outlet= 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
          
          if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) Outlet_TagBound = config->GetMarker_ActDiskOutlet_TagBound(iMarker_Outlet);
          if (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST) Outlet_TagBound = config->GetMarker_EngineExhaust_TagBound(iMarker_Outlet);
          
          if (config->GetMarker_All_TagBound(iMarker) == Outlet_TagBound) {
            Outlet_MassFlow_Local[iMarker_Outlet]                           += Outlet_MassFlow[iMarker];
            Outlet_Pressure_Local[iMarker_Outlet]                               += Outlet_Pressure[iMarker];
            Outlet_TotalPressure_Local[iMarker_Outlet]          += Outlet_TotalPressure[iMarker];
            Outlet_Temperature_Local[iMarker_Outlet]                    += Outlet_Temperature[iMarker];
            Outlet_TotalTemperature_Local[iMarker_Outlet]   += Outlet_TotalTemperature[iMarker];
            Outlet_GrossThrust_Local[iMarker_Outlet]                    += Outlet_GrossThrust[iMarker];
            Outlet_Force_Local[iMarker_Outlet]                              += Outlet_Force[iMarker];
            Outlet_Power_Local[iMarker_Outlet]                              += Outlet_Power[iMarker];
            Outlet_Area_Local[iMarker_Outlet]                                               += Outlet_Area[iMarker];
          }
          
        }
      }
      
    }
    
    /*--- Correct the min max values for the MPI ---*/
    
    bool ActDisk  = false;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) { ActDisk  = true; break; }
    }
    
    if (!ActDisk) {
      for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
        Inlet_MinPressure_Local[iMarker_Inlet]          = 1E10;
        Inlet_MaxPressure_Local[iMarker_Inlet]      = -1E10;
      }
    }
    
    /*--- All the ranks to compute the total value ---*/
    
#ifdef HAVE_MPI
    
    SU2_MPI::Allreduce(Inlet_MassFlow_Local, Inlet_MassFlow_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_ReverseMassFlow_Local, Inlet_ReverseMassFlow_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Pressure_Local, Inlet_Pressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Mach_Local, Inlet_Mach_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_MinPressure_Local, Inlet_MinPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_MaxPressure_Local, Inlet_MaxPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_TotalPressure_Local, Inlet_TotalPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Temperature_Local, Inlet_Temperature_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_TotalTemperature_Local, Inlet_TotalTemperature_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_RamDrag_Local, Inlet_RamDrag_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Force_Local, Inlet_Force_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Power_Local, Inlet_Power_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_Area_Local, Inlet_Area_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_XCG_Local, Inlet_XCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Inlet_YCG_Local, Inlet_YCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (nDim == 3) SU2_MPI::Allreduce(Inlet_ZCG_Local, Inlet_ZCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    SU2_MPI::Allreduce(Outlet_MassFlow_Local, Outlet_MassFlow_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Pressure_Local, Outlet_Pressure_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_TotalPressure_Local, Outlet_TotalPressure_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Temperature_Local, Outlet_Temperature_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_TotalTemperature_Local, Outlet_TotalTemperature_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_GrossThrust_Local, Outlet_GrossThrust_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Force_Local, Outlet_Force_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Power_Local, Outlet_Power_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Area_Local, Outlet_Area_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
#else
    
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      Inlet_MassFlow_Total[iMarker_Inlet]                            = Inlet_MassFlow_Local[iMarker_Inlet];
      Inlet_ReverseMassFlow_Total[iMarker_Inlet] = Inlet_ReverseMassFlow_Local[iMarker_Inlet];
      Inlet_Pressure_Total[iMarker_Inlet]                              = Inlet_Pressure_Local[iMarker_Inlet];
      Inlet_Mach_Total[iMarker_Inlet]                                        = Inlet_Mach_Local[iMarker_Inlet];
      Inlet_MinPressure_Total[iMarker_Inlet]                 = Inlet_MinPressure_Local[iMarker_Inlet];
      Inlet_MaxPressure_Total[iMarker_Inlet]             = Inlet_MaxPressure_Local[iMarker_Inlet];
      Inlet_TotalPressure_Total[iMarker_Inlet]               = Inlet_TotalPressure_Local[iMarker_Inlet];
      Inlet_Temperature_Total[iMarker_Inlet]                     = Inlet_Temperature_Local[iMarker_Inlet];
      Inlet_TotalTemperature_Total[iMarker_Inlet]    = Inlet_TotalTemperature_Local[iMarker_Inlet];
      Inlet_RamDrag_Total[iMarker_Inlet]                         = Inlet_RamDrag_Local[iMarker_Inlet];
      Inlet_Force_Total[iMarker_Inlet]                               = Inlet_Force_Local[iMarker_Inlet];
      Inlet_Power_Total[iMarker_Inlet]                           = Inlet_Power_Local[iMarker_Inlet];
      Inlet_Area_Total[iMarker_Inlet]                                            = Inlet_Area_Local[iMarker_Inlet];
      Inlet_XCG_Total[iMarker_Inlet]                                             = Inlet_XCG_Local[iMarker_Inlet];
      Inlet_YCG_Total[iMarker_Inlet]                                             = Inlet_YCG_Local[iMarker_Inlet];
      if (nDim == 3) Inlet_ZCG_Total[iMarker_Inlet]                                          = Inlet_ZCG_Local[iMarker_Inlet];
    }
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      Outlet_MassFlow_Total[iMarker_Outlet]                         = Outlet_MassFlow_Local[iMarker_Outlet];
      Outlet_Pressure_Total[iMarker_Outlet]                             = Outlet_Pressure_Local[iMarker_Outlet];
      Outlet_TotalPressure_Total[iMarker_Outlet]                = Outlet_TotalPressure_Local[iMarker_Outlet];
      Outlet_Temperature_Total[iMarker_Outlet]                  = Outlet_Temperature_Local[iMarker_Outlet];
      Outlet_TotalTemperature_Total[iMarker_Outlet] = Outlet_TotalTemperature_Local[iMarker_Outlet];
      Outlet_GrossThrust_Total[iMarker_Outlet]              = Outlet_GrossThrust_Local[iMarker_Outlet];
      Outlet_Force_Total[iMarker_Outlet]                          = Outlet_Force_Local[iMarker_Outlet];
      Outlet_Power_Total[iMarker_Outlet]                          = Outlet_Power_Local[iMarker_Outlet];
      Outlet_Area_Total[iMarker_Outlet]                                         = Outlet_Area_Local[iMarker_Outlet];
    }
    
#endif
    
    /*--- Compute the value of the average surface temperature and pressure and
     set the value in the config structure for future use ---*/
    
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      if (Inlet_Area_Total[iMarker_Inlet] != 0.0) {
        Inlet_Pressure_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_Mach_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_TotalPressure_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_Temperature_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_TotalTemperature_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
        Inlet_XCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
        Inlet_YCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
        if (nDim == 3) Inlet_ZCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
      }
      else {
        Inlet_Pressure_Total[iMarker_Inlet] = 0.0;
        Inlet_Mach_Total[iMarker_Inlet] = 0.0;
        Inlet_TotalPressure_Total[iMarker_Inlet] = 0.0;
        Inlet_Temperature_Total[iMarker_Inlet] = 0.0;
        Inlet_TotalTemperature_Total[iMarker_Inlet] = 0.0;
        Inlet_XCG_Total[iMarker_Inlet] = 0.0;
        Inlet_YCG_Total[iMarker_Inlet] = 0.0;
        if (nDim == 3) Inlet_ZCG_Total[iMarker_Inlet] = 0.0;
      }
      
      if (iMesh == MESH_0) {
        
        if (Engine) {
          config->SetInflow_MassFlow(iMarker_Inlet, Inlet_MassFlow_Total[iMarker_Inlet]);
          config->SetInflow_ReverseMassFlow(iMarker_Inlet, Inlet_ReverseMassFlow_Total[iMarker_Inlet]);
          config->SetInflow_Pressure(iMarker_Inlet, Inlet_Pressure_Total[iMarker_Inlet]);
          config->SetInflow_TotalPressure(iMarker_Inlet, Inlet_TotalPressure_Total[iMarker_Inlet]);
          config->SetInflow_Temperature(iMarker_Inlet, Inlet_Temperature_Total[iMarker_Inlet]);
          config->SetInflow_TotalTemperature(iMarker_Inlet, Inlet_TotalTemperature_Total[iMarker_Inlet]);
          config->SetInflow_RamDrag(iMarker_Inlet, Inlet_RamDrag_Total[iMarker_Inlet]);
          config->SetInflow_Force(iMarker_Inlet, Inlet_Force_Total[iMarker_Inlet]);
          config->SetInflow_Power(iMarker_Inlet, Inlet_Power_Total[iMarker_Inlet]);
        }
        else {
          config->SetActDiskInlet_MassFlow(iMarker_Inlet, Inlet_MassFlow_Total[iMarker_Inlet]);
          config->SetActDiskInlet_ReverseMassFlow(iMarker_Inlet, Inlet_ReverseMassFlow_Total[iMarker_Inlet]);
          config->SetActDiskInlet_Pressure(iMarker_Inlet, Inlet_Pressure_Total[iMarker_Inlet]);
          config->SetActDiskInlet_TotalPressure(iMarker_Inlet, Inlet_TotalPressure_Total[iMarker_Inlet]);
          config->SetActDiskInlet_Temperature(iMarker_Inlet, Inlet_Temperature_Total[iMarker_Inlet]);
          config->SetActDiskInlet_TotalTemperature(iMarker_Inlet, Inlet_TotalTemperature_Total[iMarker_Inlet]);
          config->SetActDiskInlet_RamDrag(iMarker_Inlet, Inlet_RamDrag_Total[iMarker_Inlet]);
          config->SetActDiskInlet_Force(iMarker_Inlet, Inlet_Force_Total[iMarker_Inlet]);
          config->SetActDiskInlet_Power(iMarker_Inlet, Inlet_Power_Total[iMarker_Inlet]);
        }
        
      }
      
    }
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      if (Outlet_Area_Total[iMarker_Outlet] != 0.0) {
        Outlet_Pressure_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
        Outlet_TotalPressure_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
        Outlet_Temperature_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
        Outlet_TotalTemperature_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
      }
      else {
        Outlet_Pressure_Total[iMarker_Outlet] = 0.0;
        Outlet_TotalPressure_Total[iMarker_Outlet] = 0.0;
        Outlet_Temperature_Total[iMarker_Outlet] = 0.0;
        Outlet_TotalTemperature_Total[iMarker_Outlet] = 0.0;
      }
      
      if (iMesh == MESH_0) {
        
        if (Engine) {
          config->SetExhaust_MassFlow(iMarker_Outlet, Outlet_MassFlow_Total[iMarker_Outlet]);
          config->SetExhaust_Pressure(iMarker_Outlet, Outlet_Pressure_Total[iMarker_Outlet]);
          config->SetExhaust_TotalPressure(iMarker_Outlet, Outlet_TotalPressure_Total[iMarker_Outlet]);
          config->SetExhaust_Temperature(iMarker_Outlet, Outlet_Temperature_Total[iMarker_Outlet]);
          config->SetExhaust_TotalTemperature(iMarker_Outlet, Outlet_TotalTemperature_Total[iMarker_Outlet]);
          config->SetExhaust_GrossThrust(iMarker_Outlet, Outlet_GrossThrust_Total[iMarker_Outlet]);
          config->SetExhaust_Force(iMarker_Outlet, Outlet_Force_Total[iMarker_Outlet]);
          config->SetExhaust_Power(iMarker_Outlet, Outlet_Power_Total[iMarker_Outlet]);
        }
        else {
          config->SetActDiskOutlet_MassFlow(iMarker_Outlet, Outlet_MassFlow_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_Pressure(iMarker_Outlet, Outlet_Pressure_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_TotalPressure(iMarker_Outlet, Outlet_TotalPressure_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_Temperature(iMarker_Outlet, Outlet_Temperature_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_TotalTemperature(iMarker_Outlet, Outlet_TotalTemperature_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_GrossThrust(iMarker_Outlet, Outlet_GrossThrust_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_Force(iMarker_Outlet, Outlet_Force_Total[iMarker_Outlet]);
          config->SetActDiskOutlet_Power(iMarker_Outlet, Outlet_Power_Total[iMarker_Outlet]);
        }
        
      }
      
    }
    
    
    if (Pair) {
      
      /*--- Store delta pressure, temperature, thrust, and area ---*/
      
      for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
        
        if (Engine) {
          Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
          jMarker = config->GetMarker_CfgFile_EngineExhaust(Inlet_TagBound);
          Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
        }
        else {
          Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
          jMarker = config->GetMarker_CfgFile_ActDiskOutlet(Inlet_TagBound);
          Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
        }
        
        
        su2double DeltaPress = 0.0, DeltaTemp = 0.0, NetThrust = 0.0, GrossThrust = 0.0, TotalPressRatio = 0.0, TotalTempRatio = 0.0, StaticPressRatio = 0.0, StaticTempRatio = 0.0;
        
        if (Engine) {
          DeltaPress   = config->GetExhaust_Pressure(Outlet_TagBound) - config->GetInflow_Pressure(Inlet_TagBound);
          DeltaTemp         = config->GetExhaust_Temperature(Outlet_TagBound) - config->GetInflow_Temperature(Inlet_TagBound);
          NetThrust    = config->GetExhaust_GrossThrust(Outlet_TagBound) - config->GetInflow_RamDrag(Inlet_TagBound);
          GrossThrust   = config->GetExhaust_GrossThrust(Outlet_TagBound);
          TotalPressRatio   = config->GetExhaust_TotalPressure(Outlet_TagBound)/config->GetInflow_TotalPressure(Inlet_TagBound);
          TotalTempRatio    = config->GetExhaust_TotalTemperature(Outlet_TagBound)/config->GetInflow_TotalTemperature(Inlet_TagBound);
          StaticPressRatio   = config->GetExhaust_Pressure(Outlet_TagBound)/config->GetInflow_Pressure(Inlet_TagBound);
          StaticTempRatio    = config->GetExhaust_Temperature(Outlet_TagBound)/config->GetInflow_Temperature(Inlet_TagBound);
          Force = config->GetInflow_Force(Inlet_TagBound) + config->GetExhaust_Force(Outlet_TagBound);
          Power = config->GetExhaust_Power(Outlet_TagBound) - config->GetInflow_Power(Inlet_TagBound);
        }
        else {
          DeltaPress   = config->GetActDiskOutlet_Pressure(Outlet_TagBound) - config->GetActDiskInlet_Pressure(Inlet_TagBound);
          DeltaTemp         = config->GetActDiskOutlet_Temperature(Outlet_TagBound) - config->GetActDiskInlet_Temperature(Inlet_TagBound);
          NetThrust    = config->GetActDiskOutlet_GrossThrust(Outlet_TagBound) - config->GetActDiskInlet_RamDrag(Inlet_TagBound);
          GrossThrust   = config->GetActDiskOutlet_GrossThrust(Outlet_TagBound);
          TotalPressRatio   = config->GetActDiskOutlet_TotalPressure(Outlet_TagBound)/config->GetActDiskInlet_TotalPressure(Inlet_TagBound);
          TotalTempRatio    = config->GetActDiskOutlet_TotalTemperature(Outlet_TagBound)/config->GetActDiskInlet_TotalTemperature(Inlet_TagBound);
          StaticPressRatio   = config->GetActDiskOutlet_Pressure(Outlet_TagBound)/config->GetActDiskInlet_Pressure(Inlet_TagBound);
          StaticTempRatio    = config->GetActDiskOutlet_Temperature(Outlet_TagBound)/config->GetActDiskInlet_Temperature(Inlet_TagBound);
          Force = config->GetActDiskInlet_Force(Inlet_TagBound) + config->GetActDiskOutlet_Force(Outlet_TagBound);
          Power =  config->GetActDiskOutlet_Power(Outlet_TagBound) - config->GetActDiskInlet_Power(Inlet_TagBound);
          MassFlow =  config->GetActDiskInlet_MassFlow(Inlet_TagBound);
        }
        
        Mach        = Inlet_Mach_Total[iMarker_Inlet];
        Area              = Inlet_Area_Total[iMarker_Inlet];
        
        if (Engine) {
          config->SetEngine_Mach(iMarker_Inlet, Mach);
          config->SetEngine_Force(iMarker_Inlet, Force);
          config->SetEngine_Power(iMarker_Inlet, Power);
          config->SetEngine_NetThrust(iMarker_Inlet, NetThrust);
          config->SetEngine_GrossThrust(iMarker_Inlet, GrossThrust);
          config->SetEngine_Area(iMarker_Inlet, Area);
        }
        else {
          config->SetActDisk_DeltaPress(iMarker_Inlet, DeltaPress);
          config->SetActDisk_DeltaTemp(iMarker_Inlet, DeltaTemp);
          config->SetActDisk_Mach(iMarker_Inlet, Mach);
          config->SetActDisk_Force(iMarker_Inlet, Force);
          config->SetActDisk_Power(iMarker_Inlet, Power);
          config->SetActDisk_MassFlow(iMarker_Inlet, MassFlow);
          config->SetActDisk_TotalPressRatio(iMarker_Inlet, TotalPressRatio);
          config->SetActDisk_TotalTempRatio(iMarker_Inlet, TotalTempRatio);
          config->SetActDisk_StaticPressRatio(iMarker_Inlet, StaticPressRatio);
          config->SetActDisk_StaticTempRatio(iMarker_Inlet, StaticTempRatio);
          config->SetActDisk_NetThrust(iMarker_Inlet, NetThrust);
          config->SetActDisk_GrossThrust(iMarker_Inlet, GrossThrust);
          config->SetActDisk_Area(iMarker_Inlet, Area);
        }
        
      }
      
      /*--- Screen output using the values already stored in the config container ---*/
      
      if ((rank == MASTER_NODE) && (iMesh == MESH_0) ) {
        
        cout.precision(5);
        cout.setf(ios::fixed, ios::floatfield);
        
        if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
          if (Engine) cout << endl   << "---------------------------- Engine properties --------------------------" << endl;
          else cout << endl   << "------------------------ Actuator Disk properties -----------------------" << endl;
        }
        
        for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
          
          if (Engine) {
            Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
            jMarker = config->GetMarker_CfgFile_EngineExhaust(Inlet_TagBound);
            Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
          }
          else {
            Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
            jMarker = config->GetMarker_CfgFile_ActDiskOutlet(Inlet_TagBound);
            Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
          }
          
          
          if (Engine) {
            NetThrust             =  config->GetEngine_NetThrust(iMarker_Inlet);
            GrossThrust   = config->GetEngine_GrossThrust(iMarker_Inlet);
            Power               = config->GetEngine_Power(iMarker_Inlet);
            Mach                  = config->GetEngine_Mach(iMarker_Inlet);
            Force               = config->GetEngine_Force(iMarker_Inlet);
          }
          else {
            DeltaPress      = config->GetActDisk_DeltaPress(iMarker_Inlet);
            DeltaTemp         = config->GetActDisk_DeltaTemp(iMarker_Inlet);
            TotalPressRatio       = config->GetActDisk_TotalPressRatio(iMarker_Inlet);
            TotalTempRatio        = config->GetActDisk_TotalTempRatio(iMarker_Inlet);
            StaticPressRatio      = config->GetActDisk_StaticPressRatio(iMarker_Inlet);
            StaticTempRatio           = config->GetActDisk_StaticTempRatio(iMarker_Inlet);
            NetThrust             =  config->GetActDisk_NetThrust(iMarker_Inlet);
            GrossThrust   = config->GetActDisk_GrossThrust(iMarker_Inlet);
            Power               = config->GetActDisk_Power(iMarker_Inlet);
            Mach                  = config->GetActDisk_Mach(iMarker_Inlet);
            Force               = config->GetActDisk_Force(iMarker_Inlet);
          }
          
          su2double Mach_Inf                  = config->GetMach();
          su2double Pressure_Inf  = config->GetPressure_FreeStreamND();
          
          su2double TotalPressure_Inf  = Pressure_Inf * pow( 1.0 + Mach_Inf * Mach_Inf * 0.5 * (Gamma - 1.0), Gamma     / (Gamma - 1.0));
          
          su2double MinPressure = Inlet_MinPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
          su2double MaxPressure = Inlet_MaxPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
          su2double AvePressure = Inlet_TotalPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
          
          su2double RefDensity  = Density_Inf;
          su2double RefAreaCoeff     = config->GetRefAreaCoeff();
          su2double RefVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
          
          su2double Factor = (0.5*RefDensity*RefAreaCoeff*RefVel2);
          su2double Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
          su2double DmT = GetTotal_CD() * Factor;
          
//          su2double ModDmT = 0.0;
//          if (nDim == 2) ModDmT = sqrt(GetTotal_CFx()*GetTotal_CFx() +
//                                       GetTotal_CFy()*GetTotal_CFy());
//
//          if (nDim == 3) ModDmT = sqrt(GetTotal_CFx()*GetTotal_CFx() +
//                                       GetTotal_CFy()*GetTotal_CFy() +
//                                       GetTotal_CFz()*GetTotal_CFz());
//
//          DmTVector[0] = GetTotal_CFx()/ModDmT;
//          DmTVector[1] = GetTotal_CFy()/ModDmT;
//          if (nDim == 3)  DmTVector[2] = GetTotal_CFz()/ModDmT;
          
          /*--- Set the aero drag ---*/

          su2double Aero_Drag = DmT - Force;
          su2double Aero_CD = Aero_Drag / Factor;

          SetTotal_AeroCD(Aero_CD);

          /*--- Set the solid surface drag ---*/
          
          su2double SolidSurf_Drag = DmT - Force;
          su2double SolidSurf_CD = SolidSurf_Drag / Factor;

          SetTotal_CD_SolidSurf(SolidSurf_CD);
          
          /*--- Set the net thrust value---*/
          
          su2double CT = NetThrust / Factor;

          SetTotal_NetCThrust(CT);
          
          /*--- Set the total power ---*/
          
          su2double PowerHP = Power * Ref *  config->GetVelocity_Ref() / 550.0;
          
          SetTotal_Power(PowerHP);
          
          /*--- Set the total ReverseFlow ---*/
          
          su2double ReverseFlow;
          if (Engine) ReverseFlow = fabs(config->GetInflow_ReverseMassFlow(iMarker_Inlet)  / config->GetInflow_MassFlow(Inlet_TagBound));
          else ReverseFlow = fabs(config->GetActDisk_ReverseMassFlow(iMarker_Inlet)  / config->GetActDiskInlet_MassFlow(Inlet_TagBound));
          
          SetTotal_ReverseFlow(ReverseFlow);
          
          /*--- Set the total mass flow ratio ---*/
          
          InfVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) InfVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
          if (Engine) MFR =fabs(config->GetInflow_MassFlow(Inlet_TagBound)) / (Density_Inf * sqrt(InfVel2) * config->GetHighlite_Area());
          else MFR = fabs(config->GetActDiskInlet_MassFlow(Inlet_TagBound)) / (Density_Inf * sqrt(InfVel2) * config->GetHighlite_Area());
          SetTotal_MFR(MFR);
          
          /*--- Evaluate shaft power and adiabatic efficiency (average) ---*/
          
          su2double Pstatic1, P1, P2, T1, T2;
          if (Engine) {
            Pstatic1 = config->GetInflow_Pressure(Inlet_TagBound);
            P1 = config->GetInflow_TotalPressure(Inlet_TagBound);
            P2 = config->GetExhaust_TotalPressure(Outlet_TagBound);
            T1 = config->GetInflow_TotalTemperature(Inlet_TagBound);
            T2 = config->GetExhaust_TotalTemperature(Outlet_TagBound);
          }
          else {
            Pstatic1 = config->GetActDiskInlet_Pressure(Inlet_TagBound);
            P1 = config->GetActDiskInlet_TotalPressure(Inlet_TagBound);
            P2 = config->GetActDiskOutlet_TotalPressure(Outlet_TagBound);
            T1 = config->GetActDiskInlet_TotalTemperature(Inlet_TagBound);
            T2 = config->GetActDiskOutlet_TotalTemperature(Outlet_TagBound);
          }
          
          /*-- Set the propulsive efficiency ---*/
          
          su2double mu_prop = fabs(DmT)*sqrt(RefVel2)/Power;
          SetTotal_Prop_Eff(mu_prop);
          
          /*-- Set the bypass propulsive efficiency ---*/
          
          su2double mu_bypass_prop = NetThrust*sqrt(RefVel2)/Power;
          SetTotal_ByPassProp_Eff(mu_bypass_prop);
          
          /*-- Set the fan adiabatic efficiency ---*/
          
          su2double mu_isentropic = 0.0;
          if ((P2/P1) > 0.0) mu_isentropic =    (T1/(T2-T1))*(pow((P2/P1),(Gamma-1.0)/Gamma)-1.0);
          SetTotal_Adiab_Eff(mu_isentropic);
          
          /*-- Set the polytropic efficiency ---*/
          
          su2double poly_coeff = 1.0/(1.0-log(T2/T1)/log(P2/P1));
          su2double mu_polytropic = ((Gamma-1.0)/Gamma)/((poly_coeff-1.0)/poly_coeff);
          SetTotal_Poly_Eff(mu_polytropic);
          
          if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
            
            if (iMarker_Inlet > 0) cout << endl;
            
            /*--- Geometry defintion ---*/
            
            if (Engine) cout <<"Engine surfaces: " << Inlet_TagBound << ", " << Outlet_TagBound << "." << endl;
            else cout <<"Actuator disk surfaces: " << Inlet_TagBound << ", " << Outlet_TagBound << "." << endl;
            
            if (nDim == 2) {
              if (config->GetSystemMeasurements() == SI)
                cout <<"CG (m): (" << Inlet_XCG_Total[iMarker_Inlet] <<", " << Inlet_YCG_Total[iMarker_Inlet] << "). Length (m): " << Inlet_Area_Total[iMarker_Inlet] << "." << endl;
              else if (config->GetSystemMeasurements() == US)
                cout <<"CG (in): (" << Inlet_XCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_YCG_Total[iMarker_Inlet]*12.0 << "). Length (in): " << Inlet_Area_Total[iMarker_Inlet]*12.0 << "." << endl;
              cout << endl;
            }
            
            if (nDim ==3) {
              if (config->GetSystemMeasurements() == SI)
                cout <<"CG (m): (" << Inlet_XCG_Total[iMarker_Inlet] <<", " << Inlet_YCG_Total[iMarker_Inlet] <<", " << Inlet_ZCG_Total[iMarker_Inlet] << "). Area (m^2): " << Inlet_Area_Total[iMarker_Inlet] << ". Radius (m): " << sqrt(Inlet_Area_Total[iMarker_Inlet]/PI_NUMBER) << "." << endl;
              else if (config->GetSystemMeasurements() == US)
                cout <<"CG (in): (" << Inlet_XCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_YCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_ZCG_Total[iMarker_Inlet]*12.0 << "). Area (in^2): " << Inlet_Area_Total[iMarker_Inlet]*12.0*12.0 << "." << endl;
              cout << endl;
            }
            
            
            /*--- Flow field descritption  ---*/
            
            if (config->GetSystemMeasurements() == SI) {
              cout << setprecision(2) << "Inlet Ave. P (Pa): " << Pstatic1*config->GetPressure_Ref() << setprecision(3) <<  ". Inlet Ave. Mach: " << Mach << "." << endl;
              cout << setprecision(2) << "Outlet Ave. PT (Pa): " << P2*config->GetPressure_Ref() << ". Outlet Ave. TT (K): " << T2*config->GetTemperature_Ref() << "." << endl;
            }
            else if (config->GetSystemMeasurements() == US) {
              cout << setprecision(2) << "Inlet Ave. P (psf): " << Pstatic1*config->GetPressure_Ref() << setprecision(3) <<  ". Inlet Ave. Mach: " << Mach << "." << endl;
              cout << setprecision(2) << "Outlet Ave. PT (psf): " << P2*config->GetPressure_Ref() << ". Outlet Ave. TT (R): " << T2*config->GetTemperature_Ref() << "." << endl;
            }
            
            cout << "Inlet min. PT/PTinf: " << MinPressure << ". Inlet max. PT/PTinf: " << MaxPressure << ". Inlet Ave. PT/PTinf: " << AvePressure << endl;
            
            su2double InfVel2, Inlet_MassFlow, Outlet_MassFlow;
            
            if (Engine) Inlet_MassFlow = fabs(config->GetInflow_MassFlow(Inlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
            else Inlet_MassFlow = fabs(config->GetActDiskInlet_MassFlow(Inlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
            
            if (config->GetSystemMeasurements() == SI) { cout << "Inlet mass flow (kg/s): "; cout << setprecision(2) << Inlet_MassFlow; }
            else if (config->GetSystemMeasurements() == US) { cout << "Inlet mass flow (lbs/s): "; cout << setprecision(2) << Inlet_MassFlow * 32.174; }
            
            if (Engine) Outlet_MassFlow = fabs(config->GetExhaust_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
            else Outlet_MassFlow = fabs(config->GetActDiskOutlet_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
            
            //          if (config->GetSystemMeasurements() == SI) { cout << ". Outlet mass flow (kg/s): "; cout << setprecision(2) << Outlet_MassFlow; }
            //          else if (config->GetSystemMeasurements() == US) { cout << ". Outlet mass flow (lbs/s): "; cout << setprecision(2) << Outlet_MassFlow * 32.174; }
            
            if (Inlet_MassFlow > Outlet_MassFlow) cout << ". I/O diff.: " << setprecision(2) << 100.0*fabs(1.0-(Outlet_MassFlow/Inlet_MassFlow)) << "%";
            else cout << ". I/O diff.: " << setprecision(2) << -100.0*fabs(1.0-(Inlet_MassFlow/Outlet_MassFlow)) << "%";
            
            InfVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) InfVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
            cout << setprecision(2) << ". MFR: " << MFR << "." << endl;
            
            if (!Engine) {
              
              cout << setprecision(3) << "PT in/out ratio: " << TotalPressRatio << ". TT in/out ratio: " << TotalTempRatio  <<  "." << endl;
              
              if (config->GetActDisk_Jump() == VARIABLES_JUMP) {
                if (config->GetSystemMeasurements() == SI) cout << setprecision(3) << "P in/out jump (Pa): ";
                else if (config->GetSystemMeasurements() == US) cout << setprecision(3) << "P in/out jump (psf): ";
                cout << setprecision(3) << DeltaPress * config->GetPressure_Ref();
                if (config->GetSystemMeasurements() == SI) cout << setprecision(3) << ". T in/out jump (K): ";
                else if (config->GetSystemMeasurements() == US) cout << setprecision(3) << ". T in/out jump (R): ";
                cout << setprecision(3) << DeltaTemp * config->GetTemperature_Ref() <<"."<< endl;
              }
              else  if (config->GetActDisk_Jump() == RATIO) {
                cout << setprecision(3) << "P in/out ratio: ";
                cout << setprecision(3) << StaticPressRatio;
                cout  << setprecision(3) <<". T in/out ratio: ";
                cout << setprecision(3) << StaticTempRatio <<"."<< endl;
              }
            }
            
            cout << setprecision(1) << "\nProp. eff. (D-T.V/Shaft P): " << 100*mu_prop << "%. By-pass prop. eff. (NetT.V/Shaft P): " << 100*mu_bypass_prop <<   "%." << endl;
            cout << setprecision(1) << "Fan adiabatic eff.: " << 100*mu_isentropic  << "%. Fan poly. eff.: " << 100*mu_polytropic << "%. Poly coeff. (n): " << setprecision(4) << poly_coeff << "." << endl;
            
            cout << endl;
            
            
            /*--- Forces descritption  ---*/
            
            if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << "Ram Drag (N): ";
            else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << "Ram Drag (lbf): ";
            cout << (GrossThrust-NetThrust) * Ref;
            
            if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << ". Gross Thrust (N): ";
            else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << ". Gross Thrust (lbf): ";
            cout << -GrossThrust * Ref  << "." << endl;
            
            if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << "Open surfaces Thurst (N): ";
            else if (config->GetSystemMeasurements() == US) cout  << setprecision(1) << "Open surfaces Thrust (lbf): ";
            cout<< setprecision(1) << Force * Ref << ". Open surfaces CT: " << setprecision(5) << -Force / Factor << "." << endl;
            
            if (config->GetSystemMeasurements() == SI) cout << "Solid surfaces Drag (N): ";
            else if (config->GetSystemMeasurements() == US) cout << "Solid surfaces Drag (lbf): ";
            cout << setprecision(1) << SolidSurf_Drag * Ref << ". Solid surfaces CD: " << setprecision(5) << SolidSurf_CD << "." << endl;

            if (config->GetSystemMeasurements() == SI) cout << "Aero Drag (N): ";
            else if (config->GetSystemMeasurements() == US) cout << "Aero Drag (lbf): ";
            cout << setprecision(1) << Aero_Drag * Ref << ". Aero CD: " << setprecision(5) << Aero_CD << "." << endl;

            if (config->GetSystemMeasurements() == SI) cout << setprecision(1) <<"Net Thrust (N): ";
            else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << "Net Thrust (lbf): ";
            cout << setprecision(5) << -NetThrust * Ref  << ". Net CT: " << CT;
            
            if (config->GetSystemMeasurements() == SI) {
              cout << ". Power (W): ";
              cout << setprecision(1) << Power * Ref *  config->GetVelocity_Ref()  << "." << endl;
            }
            else if (config->GetSystemMeasurements() == US) {
              cout << ". Power (HP): ";
              cout << setprecision(1) << Power * Ref *  config->GetVelocity_Ref() / 550.0 << "." << endl;
            }
            
          }
          
        }
        
        if (write_heads && Output && !config->GetDiscrete_Adjoint()) cout << "-------------------------------------------------------------------------" << endl << endl;
        
      }
      
    }
    
    
    delete [] Outlet_MassFlow_Local;
    delete [] Outlet_Temperature_Local;
    delete [] Outlet_TotalTemperature_Local;
    delete [] Outlet_Pressure_Local;
    delete [] Outlet_TotalPressure_Local;
    delete [] Outlet_Area_Local;
    delete [] Outlet_GrossThrust_Local;
    delete [] Outlet_Force_Local;
    delete [] Outlet_Power_Local;
    
    delete [] Outlet_MassFlow_Total;
    delete [] Outlet_Temperature_Total;
    delete [] Outlet_TotalTemperature_Total;
    delete [] Outlet_Pressure_Total;
    delete [] Outlet_TotalPressure_Total;
    delete [] Outlet_Area_Total;
    delete [] Outlet_GrossThrust_Total;
    delete [] Outlet_Force_Total;
    delete [] Outlet_Power_Total;
    
    delete [] Inlet_MassFlow_Local;
    delete [] Inlet_ReverseMassFlow_Local;
    delete [] Inlet_Temperature_Local;
    delete [] Inlet_TotalTemperature_Local;
    delete [] Inlet_Pressure_Local;
    delete [] Inlet_Mach_Local;
    delete [] Inlet_MinPressure_Local;
    delete [] Inlet_MaxPressure_Local;
    delete [] Inlet_TotalPressure_Local;
    delete [] Inlet_Area_Local;
    delete [] Inlet_RamDrag_Local;
    delete [] Inlet_Force_Local;
    delete [] Inlet_Power_Local;
    delete [] Inlet_XCG_Local;
    delete [] Inlet_YCG_Local;
    delete [] Inlet_ZCG_Local;
    
    delete [] Inlet_MassFlow_Total;
    delete [] Inlet_ReverseMassFlow_Total;
    delete [] Inlet_Temperature_Total;
    delete [] Inlet_TotalTemperature_Total;
    delete [] Inlet_Pressure_Total;
    delete [] Inlet_Mach_Total;
    delete [] Inlet_MinPressure_Total;
    delete [] Inlet_MaxPressure_Total;
    delete [] Inlet_TotalPressure_Total;
    delete [] Inlet_Area_Total;
    delete [] Inlet_RamDrag_Total;
    delete [] Inlet_Force_Total;
    delete [] Inlet_Power_Total;
    delete [] Inlet_XCG_Total;
    delete [] Inlet_YCG_Total;
    delete [] Inlet_ZCG_Total;
    
    delete [] Inlet_MassFlow;
    delete [] Inlet_Mach;
    delete [] Inlet_MinPressure;
    delete [] Inlet_MaxPressure;
    delete [] Inlet_ReverseMassFlow;
    delete [] Inlet_Pressure;
    delete [] Inlet_TotalPressure;
    delete [] Inlet_Temperature;
    delete [] Inlet_TotalTemperature;
    delete [] Inlet_Area;
    delete [] Inlet_RamDrag;
    delete [] Inlet_Force;
    delete [] Inlet_Power;
    delete [] Inlet_XCG;
    delete [] Inlet_YCG;
    delete [] Inlet_ZCG;
    
    delete [] Outlet_MassFlow;
    delete [] Outlet_Pressure;
    delete [] Outlet_TotalPressure;
    delete [] Outlet_Temperature;
    delete [] Outlet_TotalTemperature;
    delete [] Outlet_Area;
    delete [] Outlet_GrossThrust;
    delete [] Outlet_Force;
    delete [] Outlet_Power;
    
  }
  
}

void CEulerSolver::SetActDisk_BCThrust(CGeometry *geometry, CSolver **solver_container,
                                       CConfig *config, unsigned short iMesh, bool Output) {
  
  su2double Massflow = 0.0 , Target_Massflow = 0.0, DragMinusThrust = 0.0 , Target_DragMinusThrust = 0.0, Target_NetThrust = 0.0, BCThrust = 0.0, BCThrust_inc = 0.0;
  unsigned short iDim, iMarker;
  unsigned long iVertex, iPoint;
  su2double  *V_inlet = NULL, Pressure,
  Density, T0_Ti, ATerm, BTerm, LHS, RHS, RHS_PDelta, RHS_MDelta, F, DF_DLa, CTerm_, DTerm_,
  ETerm, La, La_old, TotalArea, To_Ti, DeltaT, Po_Pi, DeltaP, Area, Velocity_Normal,
  SoundSpeed2, Force_Normal,
  RefDensity, RefAreaCoeff, RefVel2, Factor, Ref;
  unsigned short iter;
  string Marker_Tag;
  su2double Target_Force, Force, Target_Power, Power, NetThrust, BCThrust_old, Initial_BCThrust;
  bool ActDisk_Info;
  su2double MyBCThrust, BCThrust_Init;
  
  su2double dNetThrust_dBCThrust        = config->GetdNetThrust_dBCThrust();
  unsigned short Kind_ActDisk           = config->GetKind_ActDisk();
  bool ratio                            = (config->GetActDisk_Jump() == RATIO);
  unsigned long Update_BCThrust         = config->GetUpdate_BCThrust();
  unsigned long Iter_Fixed_NetThrust    = config->GetIter_Fixed_NetThrust();
  unsigned long ExtIter                 = config->GetExtIter();
  bool Update_BCThrust_Bool             = false;
  bool restart                          = (config->GetRestart() || config->GetRestart_Flow());
  su2double Fan_Poly_Eff                = config->GetFan_Poly_Eff();
  su2double PolyCoeff                   = 1.0/(1.0-((Gamma-1.0)/Gamma)/Fan_Poly_Eff);
  
  RefDensity   = Density_Inf;
  RefAreaCoeff = config->GetRefAreaCoeff();
  RefVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  
  Factor = (0.5*RefDensity*RefAreaCoeff*RefVel2);
  Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
  
  /*--- Delta P and delta T are inputs ---*/
  
  if (Kind_ActDisk == VARIABLES_JUMP) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        
        if (ratio) {
          if (config->GetMach()  < 0.5) {
            DeltaP       = config->GetActDisk_PressJump(Marker_Tag, 0);
            DeltaT       = config->GetActDisk_TempJump(Marker_Tag, 0);
          }
          else {
            DeltaP       = config->GetActDisk_PressJump(Marker_Tag, 1);
            DeltaT       = config->GetActDisk_TempJump(Marker_Tag, 1);
          }
        }
        else {
          if (config->GetMach()  < 0.5) {
            DeltaP       = max(0.0, config->GetActDisk_PressJump(Marker_Tag, 0) / config->GetPressure_Ref());
            DeltaT       = max(0.0, config->GetActDisk_TempJump(Marker_Tag, 0) / config->GetTemperature_Ref());
          }
          else {
            DeltaP       = max(0.0, config->GetActDisk_PressJump(Marker_Tag, 1) / config->GetPressure_Ref());
            DeltaT       = max(0.0, config->GetActDisk_TempJump(Marker_Tag, 1) / config->GetTemperature_Ref());
          }
        }
        
        /*--- Set the Delta P, Delta T values at each discrete point (uniform distribution)  ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          SetActDisk_DeltaP(iMarker, iVertex, DeltaP);
          SetActDisk_DeltaT(iMarker, iVertex, DeltaT);
        }
        
      }
    }
  }
  
  /*--- Iteration using BCThrust ---*/
  
  else {
    
    if (ExtIter == 0) BCThrust_Counter = 0;
    
    /*--- Only the fine mesh level should check the convergence criteria ---*/
    
    if ((iMesh == MESH_0) && Output) {
      
      /*--- Initialize the update flag to false ---*/
      
      Update_BCThrust_Bool = false;
      
      /*--- Reevaluate BCThrust at a fix number of iterations ---*/
      
      if ((ExtIter % Iter_Fixed_NetThrust == 0) && (ExtIter != 0)) {
        BCThrust_Counter++;
        if ((BCThrust_Counter != 0) &&
            (BCThrust_Counter != 1) &&
            (BCThrust_Counter != Update_BCThrust) &&
            (BCThrust_Counter != Update_BCThrust + 2) &&
            (BCThrust_Counter != Update_BCThrust + 4) ) Update_BCThrust_Bool = true;
        else Update_BCThrust_Bool = false;
      }
      
      /*--- Store the update boolean for use on other mesh levels in the MG ---*/
      
      config->SetUpdate_BCThrust_Bool(Update_BCThrust_Bool);
      
    }
    
    else {
      Update_BCThrust_Bool = config->GetUpdate_BCThrust_Bool();
    }
    
    
    /*--- If it is the first iteration, set the BCThrust to a meaning full target value,
     * this can be done at an initialization level, for the time being it is OK here ---*/
    
    if (ExtIter == 0) {
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          
          if (Kind_ActDisk == NET_THRUST) {
            if (restart)
              Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            else {
              if (config->GetMach() < 0.5) Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
              else Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
            }
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
          if (Kind_ActDisk == BC_THRUST) {
            if (restart)
              Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            else {
              if (config->GetMach() < 0.5) Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
              else Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
            }
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
          if (Kind_ActDisk == POWER) {
            Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
          if (Kind_ActDisk == DRAG_MINUS_THRUST) {
            Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
          if (Kind_ActDisk == MASSFLOW) {
            Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
            config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
            config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
          }
          
        }
      }
    }
    
    /*--- Typical iteration to set the value of BC Thrust at each actuator disk ---*/
    
    if (Update_BCThrust_Bool && Output) {
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
          
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          
          if (Kind_ActDisk == NET_THRUST) {
            
            if (config->GetMach() < 0.5) {
              Target_NetThrust    = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
            }
            else {
              Target_NetThrust    = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
            }
            NetThrust    = config->GetActDisk_NetThrust(Marker_Tag);
            BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_NetThrust - NetThrust);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
          }
          
          if (Kind_ActDisk == BC_THRUST) {
            
            if (config->GetMach() < 0.5) {
              Target_Force =  fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
            }
            else {
              Target_Force =  fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
            }
            
            Force        = -config->GetActDisk_Force(Marker_Tag);
            BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_Force - Force);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
            
          }
          
          if (Kind_ActDisk == POWER) {
            
            if (config->GetMach() < 0.5) {
              Target_Power =  fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / (Ref * config->GetVelocity_Ref() /  550.0));
            }
            else {
              Target_Power =  fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / (Ref * config->GetVelocity_Ref() /  550.0));
            }
            
            Power        = config->GetActDisk_Power(Marker_Tag);
            BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_Power - Power);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
          }
          
          if (Kind_ActDisk == DRAG_MINUS_THRUST) {
            
            if (config->GetMach() < 0.5) {
              Target_DragMinusThrust  =  -fabs(config->GetActDisk_PressJump(Marker_Tag, 0)) * Factor;
            }
            else {
              Target_DragMinusThrust  =  -fabs(config->GetActDisk_PressJump(Marker_Tag, 1)) * Factor;
            }
            
            DragMinusThrust = GetTotal_CD() * Factor;
            BCThrust_old    = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc    = -(1.0/dNetThrust_dBCThrust)*(Target_DragMinusThrust - DragMinusThrust);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
            
          }
          
          if (Kind_ActDisk == MASSFLOW) {
            
            if (config->GetMach() < 0.5) {
              Target_Massflow  =  fabs(config->GetActDisk_PressJump(Marker_Tag, 0) / (config->GetDensity_Ref() * config->GetVelocity_Ref()));
              if (config->GetSystemMeasurements() == US) Target_Massflow /= 32.174;
            }
            else {
              Target_Massflow  =  fabs(config->GetActDisk_PressJump(Marker_Tag, 1) / (config->GetDensity_Ref() * config->GetVelocity_Ref()));
              if (config->GetSystemMeasurements() == US) Target_Massflow /= 32.174;
            }
            
            Massflow = config->GetActDisk_MassFlow(Marker_Tag);
            BCThrust_old    = config->GetActDisk_BCThrust_Old(Marker_Tag);
            BCThrust_inc    = (1.0/dNetThrust_dBCThrust)*(Target_Massflow - Massflow);
            if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
            else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            if (iMesh == MESH_0) {
              config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
              BCThrust_Init = BCThrust*Ref;
              config->SetInitial_BCThrust(BCThrust_Init);
            }
            
          }
          
        }
        
      }
      
      /*--- After a complete update of BC_Thrust
       update the value of BC Thrust (old) for future iterations ---*/
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if ((Kind_ActDisk == NET_THRUST) || (Kind_ActDisk == BC_THRUST) ||
              (Kind_ActDisk == POWER) || (Kind_ActDisk == DRAG_MINUS_THRUST) ||
              (Kind_ActDisk == MASSFLOW)) {
            BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
            config->SetActDisk_BCThrust_Old(Marker_Tag, BCThrust);
          }
        }
      }
      
    }
    
    /*--- Evaluate the pressure jump at each node using the total thrust ---*/
    
    if ((Update_BCThrust_Bool && Output) || (ExtIter == 0)) {
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
          
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          RefDensity  = Density_Inf;
          RefAreaCoeff = config->GetRefAreaCoeff();
          RefVel2 = 0.0; for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
          
          Factor = (0.5*RefDensity*RefAreaCoeff*RefVel2);
          Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
          BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            
            if (geometry->node[iPoint]->GetDomain()) {
              
              geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
              
              if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) {
                for (iDim = 0; iDim < nDim; iDim++) { Vector[iDim] = -Vector[iDim]; }
              }
              
              Area = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) { Area += Vector[iDim]*Vector[iDim]; }
              Area = sqrt (Area);
              
              /*--- Use the inlet state to compute the Pressure and Temperature jumps ---*/
              
              if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET)
                V_inlet = node[iPoint]->GetPrimitive();
              if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)
                V_inlet = GetDonorPrimVar(iMarker, iVertex);
              
              Density       = V_inlet[nDim+2];
              Pressure      = V_inlet[nDim+1];
              SoundSpeed2   = Pressure*Gamma/Density;
              TotalArea     = config->GetActDisk_Area(Marker_Tag);
              Force_Normal  = Area*(BCThrust/TotalArea);
              
              
              Velocity_Normal = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                Velocity_Normal += V_inlet[iDim+1]*Vector[iDim]/Area;
              }
              
              
              if (Velocity_Normal > EPS) {
                
                /*--- Ratio of the total temperature to the temperature at the inflow ---*/
                
                T0_Ti = 1.0 + ((Gamma-1.0)/SoundSpeed2)*(0.5*Velocity_Normal*Velocity_Normal + Force_Normal/(Density*Area));
                
                
                ATerm = 2.0*T0_Ti/(Gamma+1.0);
                BTerm = 0.5*(Gamma+1.0)/(Gamma-1.0);
                LHS = fabs(Velocity_Normal)/(sqrt(SoundSpeed2)*pow(ATerm,BTerm));
                
                CTerm_ = (PolyCoeff-1.0)/(PolyCoeff+1.0);
                DTerm_ = 1.0/(PolyCoeff-1.0);
                
                La = EPS, La_old = EPS;
                
                for (iter = 0; iter < 100; iter++) {
                  
                  ETerm = ((1.0-CTerm_*La*La)/(1.0-CTerm_+EPS));
                  
                  RHS = La*pow(ETerm, DTerm_);
                  
                  ETerm = ((1.0-CTerm_*(La+1E-6)*(La+1E-6))/(1.0-CTerm_+EPS));
                  RHS_PDelta = (La+1E-6)*pow(ETerm, DTerm_);
                  
                  ETerm = ((1.0-CTerm_*(La-1E-6)*(La-1E-6))/(1.0-CTerm_+EPS));
                  RHS_MDelta = (La-1E-6)*pow(ETerm, DTerm_);
                  
                  /*--- Objective function and finitte differences derivative ---*/
                  
                  F = RHS - LHS;
                  DF_DLa = (RHS_PDelta - RHS_MDelta)/2E-6;
                  
                  /*--- Newton's step ---*/
                  
                  La_old = La;
                  La = La_old - 0.75*(F/DF_DLa);
                  
                  if (fabs(F) < 1E-10) break;
                  
                }
                
                if (iter == 99) cout << "The laval number evaluation is not converging." << endl;
                
                /*--- Laval is bounded ---*/
                
                La = min(La, sqrt(6.0));  La = max(La, 0.0);
                
                To_Ti = max(1.0, T0_Ti*(1.0-CTerm_*La*La));
                SetActDisk_DeltaT(iMarker, iVertex, To_Ti);
                
                Po_Pi = max(1.0, pow(To_Ti, PolyCoeff*DTerm_));
                SetActDisk_DeltaP(iMarker, iVertex, Po_Pi);
                
              }
              else {
                SetActDisk_DeltaT(iMarker, iVertex, 1.0);
                SetActDisk_DeltaP(iMarker, iVertex, 1.0);
              }
              
            }
            
          }
        }
        
      }
    }
  }
  
  /*--- Broadcast some information to the master node ---*/
  
  ActDisk_Info = false;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
        (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
      ActDisk_Info = true;
    }
  }
  if (!ActDisk_Info) config->SetInitial_BCThrust(0.0);
  
  MyBCThrust = config->GetInitial_BCThrust();
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&MyBCThrust, &BCThrust, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  BCThrust = MyBCThrust;
#endif
  config->SetInitial_BCThrust(BCThrust);
  
}

void CEulerSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, bool Output) {
  
  su2double Target_CL = 0.0, AoA = 0.0, Vel_Infty[3], AoA_inc = 0.0, Vel_Infty_Mag, Delta_AoA, Old_AoA,
  dCL_dAlpha_, dCD_dCL_;
  
  unsigned short iDim;
  
  unsigned long Iter_Fixed_CL = config->GetIter_Fixed_CL();
  unsigned long Update_Alpha = config->GetUpdate_Alpha();
  
  unsigned long ExtIter       = config->GetExtIter();
  bool write_heads = ((ExtIter % Iter_Fixed_CL == 0) && (ExtIter != 0));
  su2double Beta                 = config->GetAoS()*PI_NUMBER/180.0;
  su2double dCL_dAlpha           = config->GetdCL_dAlpha()*180.0/PI_NUMBER;
  bool Update_AoA             = false;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (ExtIter == 0) AoA_Counter = 0;
  
  /*--- Only the fine mesh level should check the convergence criteria ---*/
  
  if ((iMesh == MESH_0) && Output) {
    
    /*--- Initialize the update flag to false ---*/
    
    Update_AoA = false;
    
    /*--- Reevaluate Angle of Attack at a fix number of iterations ---*/

    if ((ExtIter % Iter_Fixed_CL == 0) && (ExtIter != 0)) {
      AoA_Counter++;
      if ((AoA_Counter != 0) &&
          (AoA_Counter != 1) &&
          (AoA_Counter != Update_Alpha) &&
          (AoA_Counter != Update_Alpha + 2) &&
          (AoA_Counter != Update_Alpha + 4) ) Update_AoA = true;
      else Update_AoA = false;
    }
    
    /*--- Store the update boolean for use on other mesh levels in the MG ---*/
    
    config->SetUpdate_AoA(Update_AoA);
    
  }
  
  else {
    Update_AoA = config->GetUpdate_AoA();
  }
  
  if (Update_AoA && Output) {
    
    /*--- Retrieve the specified target CL value. ---*/
    
    Target_CL = config->GetTarget_CL();
    
    /*--- Retrieve the old AoA (radians) ---*/
    
    AoA_old = config->GetAoA()*PI_NUMBER/180.0;
    
    /*--- Estimate the increment in AoA based on dCL_dAlpha (radians) ---*/
    
    AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - Total_CL);
    
    /*--- Compute a new value for AoA on the fine mesh only (radians)---*/
    
    if (iMesh == MESH_0) AoA = AoA_old + AoA_inc;
    else { AoA = config->GetAoA()*PI_NUMBER/180.0; }
    
    /*--- Only the fine mesh stores the updated values for AoA in config ---*/
    
    if (iMesh == MESH_0) {
      config->SetAoA(AoA*180.0/PI_NUMBER);
    }
    
    /*--- Update the freestream velocity vector at the farfield ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty[iDim] = GetVelocity_Inf(iDim);
    
    /*--- Compute the magnitude of the free stream velocity ---*/
    
    Vel_Infty_Mag = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty_Mag += Vel_Infty[iDim]*Vel_Infty[iDim];
    Vel_Infty_Mag = sqrt(Vel_Infty_Mag);
    
    /*--- Compute the new freestream velocity with the updated AoA ---*/
    
    if (nDim == 2) {
      Vel_Infty[0] = cos(AoA)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(AoA)*Vel_Infty_Mag;
    }
    if (nDim == 3) {
      Vel_Infty[0] = cos(AoA)*cos(Beta)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(Beta)*Vel_Infty_Mag;
      Vel_Infty[2] = sin(AoA)*cos(Beta)*Vel_Infty_Mag;
    }
    
    /*--- Store the new freestream velocity vector for the next iteration ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_Inf[iDim] = Vel_Infty[iDim];
    }
    
    /*--- Only the fine mesh stores the updated values for velocity in config ---*/
    
    if (iMesh == MESH_0) {
      for (iDim = 0; iDim < nDim; iDim++)
        config->SetVelocity_FreeStreamND(Vel_Infty[iDim], iDim);
    }
    
    /*--- Output some information to the console with the headers ---*/
    
    if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && !config->GetDiscrete_Adjoint()) {
      Old_AoA = config->GetAoA() - AoA_inc*(180.0/PI_NUMBER);
      Delta_AoA = Old_AoA - AoA_Prev;
      
      cout.precision(7);
      cout.setf(ios::fixed, ios::floatfield);
      cout << endl << "----------------------------- Fixed CL Mode -----------------------------" << endl;
      cout << "CL: " << Total_CL;
      cout << " (target: " << config->GetTarget_CL() <<")." << endl;
      cout.precision(4);
      cout << "Previous AoA: " << Old_AoA << " deg";
      cout << ", new AoA: " << config->GetAoA() << " deg." << endl;
      
      cout.precision(7);
      if ((fabs(Delta_AoA) > EPS) && (AoA_Counter != 2))  {
        
        dCL_dAlpha_ = (Total_CL-Total_CL_Prev)/Delta_AoA;
        dCD_dCL_ = (Total_CD-Total_CD_Prev)/(Total_CL-Total_CL_Prev);
        
        cout << "Approx. Delta CL / Delta AoA: " << dCL_dAlpha_ << " (1/deg)." << endl;
        cout << "Approx. Delta CD / Delta CL: " << dCD_dCL_ << ". " << endl;
        
      }
      
      Total_CD_Prev = Total_CD;
      Total_CL_Prev = Total_CL;
      AoA_Prev =config->GetAoA() - AoA_inc*(180.0/PI_NUMBER);
      
      cout << "-------------------------------------------------------------------------" << endl << endl;
    }
    
  }
  
}

void CEulerSolver::Compute_ComboObj(CConfig *config) {
  unsigned short iMarker_Monitoring;
  su2double Weight_ObjFunc;

  /*--- Loop over all monitored markers, add to the 'combo' objective ---*/
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
    switch(config->GetKind_ObjFunc(iMarker_Monitoring))
    {
    case DRAG_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*(Surface_CD[iMarker_Monitoring]);
      break;
    case LIFT_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*(Surface_CL[iMarker_Monitoring]);
      break;
    case SIDEFORCE_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*(Surface_CSF[iMarker_Monitoring]);
      break;
    case EFFICIENCY:
      Total_ComboObj+=Weight_ObjFunc*(Surface_CEff[iMarker_Monitoring]);
      break;
    case MOMENT_X_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*(Surface_CMx[iMarker_Monitoring]);
      break;
    case MOMENT_Y_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*(Surface_CMy[iMarker_Monitoring]);
      break;
    case MOMENT_Z_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*(Surface_CMz[iMarker_Monitoring]);
      break;
    case FORCE_X_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*Surface_CFx[iMarker_Monitoring];
      break;
    case FORCE_Y_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*Surface_CFy[iMarker_Monitoring];
      break;
    case FORCE_Z_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*Surface_CFz[iMarker_Monitoring];
      break;
    case TOTAL_HEATFLUX:
      Total_ComboObj+=Weight_ObjFunc*Surface_HF_Visc[iMarker_Monitoring];
      break;
    case MAXIMUM_HEATFLUX:
      Total_ComboObj+=Weight_ObjFunc*Surface_MaxHF_Visc[iMarker_Monitoring];
      break;
    /*--- The following are not per-surface, and as a result will be
     * double-counted iff multiple surfaces are specified as well as multi-objective
     * TODO: print a warning to the user about that possibility. ---*/
    case EQUIVALENT_AREA:
      Total_ComboObj+=Weight_ObjFunc*Total_CEquivArea;
      break;
    case AERO_DRAG_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*Total_AeroCD;
      break;
    case RADIAL_DISTORTION:
      Total_ComboObj+=Weight_ObjFunc*Total_RadialDistortion;
      break;
    case CIRCUMFERENTIAL_DISTORTION:
      Total_ComboObj+=Weight_ObjFunc*Total_CircumferentialDistortion;
      break;
    case NEARFIELD_PRESSURE:
      Total_ComboObj+=Weight_ObjFunc*Total_CNearFieldOF;
      break;
    case INVERSE_DESIGN_PRESSURE:
      Total_ComboObj+=Weight_ObjFunc*Total_CpDiff;
      break;
    case INVERSE_DESIGN_HEATFLUX:
      Total_ComboObj+=Weight_ObjFunc*Total_HeatFluxDiff;
      break;
    case THRUST_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*Total_CT;
      break;
    case TORQUE_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*Total_CQ;
      break;
    case FIGURE_OF_MERIT:
      Total_ComboObj+=Weight_ObjFunc*Total_CMerit;
      break;
    case AVG_TOTAL_PRESSURE:
      Total_ComboObj+=Weight_ObjFunc*OneD_TotalPress;
      break;
    case AVG_OUTLET_PRESSURE:
      Total_ComboObj+=Weight_ObjFunc*OneD_PressureRef;
      break;
    case MASS_FLOW_RATE:
      Total_ComboObj+=Weight_ObjFunc*OneD_MassFlowRate;
      break;
    default:
      break;
    }
  }

}

void CEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar, kVar, jDim;
  unsigned long iPoint, iVertex;
  su2double *Normal = NULL, *GridVel = NULL, Area, UnitNormal[3], *NormalArea,
  ProjGridVel = 0.0, turb_ke;
  su2double Density_b, StaticEnergy_b, Enthalpy_b, *Velocity_b, Kappa_b, Chi_b, Energy_b, VelMagnitude2_b, Pressure_b;
  su2double Density_i, *Velocity_i, ProjVelocity_i = 0.0, Energy_i, VelMagnitude2_i;
  su2double **Jacobian_b, **DubDu;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  
  Normal = new su2double[nDim];
  NormalArea = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_i = new su2double[nDim];
  Jacobian_b = new su2double*[nVar];
  DubDu = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_b[iVar] = new su2double[nVar];
    DubDu[iVar] = new su2double[nVar];
  }
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negative for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++) {
        NormalArea[iDim] = -Normal[iDim];
        UnitNormal[iDim] = -Normal[iDim]/Area;
      }
      
      /*--- Get the state i ---*/

      VelMagnitude2_i = 0.0; ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
        VelMagnitude2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      Density_i = node[iPoint]->GetDensity();
      Energy_i = node[iPoint]->GetEnergy();

      /*--- Compute the boundary state b ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        Velocity_b[iDim] = Velocity_i[iDim] - ProjVelocity_i * UnitNormal[iDim]; //Force the velocity to be tangential to the surface.

      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        for (iDim = 0; iDim < nDim; iDim++) Velocity_b[iDim] += GridVel[iDim] - ProjGridVel * UnitNormal[iDim];
      }

      VelMagnitude2_b = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        VelMagnitude2_b += Velocity_b[iDim] * Velocity_b[iDim];

      /*--- Compute the residual ---*/

      turb_ke = 0.0;
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);

      Density_b = Density_i;
      StaticEnergy_b = Energy_i - 0.5 * VelMagnitude2_i - turb_ke;
      Energy_b = StaticEnergy_b + 0.5 * VelMagnitude2_b + turb_ke;

      FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);
      Kappa_b = FluidModel->GetdPde_rho() / Density_b;
      Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
      Pressure_b = FluidModel->GetPressure();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;

      numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, NormalArea, Residual);

      /*--- Grid velocity correction to the energy term ---*/
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
        Residual[nVar-1] += Pressure_b*ProjGridVel*Area;
      }

      /*--- Add the Reynolds stress tensor contribution ---*/

      if (tkeNeeded) {
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] += (2.0/3.0)*Density_b*turb_ke*NormalArea[iDim];
      }
      
      /*--- Add value to the residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Form Jacobians for implicit computations ---*/
      
      if (implicit) {
        
        /*--- Initialize Jacobian ---*/
        
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
        
        /*--- Compute DubDu ---*/

        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0.0;
          DubDu[iVar][iVar]= 1.0;
        }

        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim<nDim; jDim++)
            DubDu[iDim+1][jDim+1] -= UnitNormal[iDim]*UnitNormal[jDim];
        DubDu[nVar-1][0] += 0.5*ProjVelocity_i*ProjVelocity_i;
        for (iDim = 0; iDim < nDim; iDim++) {
          DubDu[nVar-1][iDim+1] -= ProjVelocity_i*UnitNormal[iDim];
        }

        /*--- Compute flux Jacobian in state b ---*/

        numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, NormalArea, 1, Jacobian_b);

        // Check for grid movement, should be already considered since Jacobian b is computed from u_b
        // if (grid_movement) {
        // Jacobian_b[nVar-1][0] += 0.5*ProjGridVel*ProjGridVel;
        // for (iDim = 0; iDim < nDim; iDim++)
        // Jacobian_b[nVar-1][iDim+1] -= ProjGridVel * UnitNormal[iDim];
        // }

        /*--- Compute numerical flux Jacobian at node i ---*/

        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            for (kVar = 0; kVar < nVar; kVar++)
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];

        /*--- Add the Jacobian to the sparse matrix ---*/

        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      }
    }
  }
  
  delete [] Normal;
  delete [] NormalArea;
  delete [] Velocity_b;
  delete [] Velocity_i;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_b[iVar];
    delete [] DubDu[iVar];
  }
  delete [] Jacobian_b;
  delete [] DubDu;
  
}

void CEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  
  su2double *GridVel;
  su2double Area, UnitNormal[3] = {0.0,0.0,0.0};
  su2double Density, Pressure, Energy,  Velocity[3] = {0.0,0.0,0.0};
  su2double Density_Bound, Pressure_Bound, Vel_Bound[3] = {0.0,0.0,0.0};
  su2double Density_Infty, Pressure_Infty, Vel_Infty[3] = {0.0,0.0,0.0};
  su2double SoundSpeed, Entropy, Velocity2, Vn;
  su2double SoundSpeed_Bound, Entropy_Bound, Vel2_Bound, Vn_Bound;
  su2double SoundSpeed_Infty, Entropy_Infty, Vel2_Infty, Vn_Infty, Qn_Infty;
  su2double RiemannPlus, RiemannMinus;
  su2double *V_infty, *V_domain;
  
  su2double Gas_Constant     = config->GetGas_ConstantND();
  
  bool implicit       = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool grid_movement  = config->GetGrid_Movement();
  bool viscous        = config->GetViscous();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS ) ||
                     (config->GetKind_Solver() == DISC_ADJ_RANS))
                    && (config->GetKind_Turb_Model() == SST));
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Allocate the value at the infinity ---*/
    V_infty = GetCharacPrimVar(val_marker, iVertex);
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Construct solution state at infinity for compressible flow by
         using Riemann invariants, and then impose a weak boundary condition
         by computing the flux using this new state for U. See CFD texts by
         Hirsch or Blazek for more detail. Adapted from an original
         implementation in the Stanford University multi-block (SUmb) solver
         in the routine bcFarfield.f90 written by Edwin van der Weide,
         last modified 06-12-2005. First, compute the unit normal at the
         boundary nodes. ---*/

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Store primitive variables (density, velocities, velocity squared,
         energy, pressure, and sound speed) at the boundary node, and set some
         other quantities for clarity. Project the current flow velocity vector
         at this boundary node into the local normal direction, i.e. compute
         v_bound.n.  ---*/

      Density_Bound = V_domain[nDim+2];
      Vel2_Bound = 0.0; Vn_Bound = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Bound[iDim] = V_domain[iDim+1];
        Vel2_Bound     += Vel_Bound[iDim]*Vel_Bound[iDim];
        Vn_Bound       += Vel_Bound[iDim]*UnitNormal[iDim];
      }
      Pressure_Bound   = node[iPoint]->GetPressure();
      SoundSpeed_Bound = sqrt(Gamma*Pressure_Bound/Density_Bound);
      Entropy_Bound    = pow(Density_Bound, Gamma)/Pressure_Bound;

      /*--- Store the primitive variable state for the freestream. Project
         the freestream velocity vector into the local normal direction,
         i.e. compute v_infty.n. ---*/

      Density_Infty = GetDensity_Inf();
      Vel2_Infty = 0.0; Vn_Infty = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Infty[iDim] = GetVelocity_Inf(iDim);
        Vel2_Infty     += Vel_Infty[iDim]*Vel_Infty[iDim];
        Vn_Infty       += Vel_Infty[iDim]*UnitNormal[iDim];
      }
      Pressure_Infty   = GetPressure_Inf();
      SoundSpeed_Infty = sqrt(Gamma*Pressure_Infty/Density_Infty);
      Entropy_Infty    = pow(Density_Infty, Gamma)/Pressure_Infty;

      /*--- Adjust the normal freestream velocity for grid movement ---*/

      Qn_Infty = Vn_Infty;
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          Qn_Infty -= GridVel[iDim]*UnitNormal[iDim];
      }

      /*--- Compute acoustic Riemann invariants: R = u.n +/- 2c/(gamma-1).
         These correspond with the eigenvalues (u+c) and (u-c), respectively,
         which represent the acoustic waves. Positive characteristics are
         incoming, and a physical boundary condition is imposed (freestream
         state). This occurs when either (u.n+c) > 0 or (u.n-c) > 0. Negative
         characteristics are leaving the domain, and numerical boundary
         conditions are required by extrapolating from the interior state
         using the Riemann invariants. This occurs when (u.n+c) < 0 or
         (u.n-c) < 0. Note that grid movement is taken into account when
         checking the sign of the eigenvalue. ---*/

      /*--- Check whether (u.n+c) is greater or less than zero ---*/

      if (Qn_Infty > -SoundSpeed_Infty) {
        /*--- Subsonic inflow or outflow ---*/
        RiemannPlus = Vn_Bound + 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Supersonic inflow ---*/
        RiemannPlus = Vn_Infty + 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }

      /*--- Check whether (u.n-c) is greater or less than zero ---*/

      if (Qn_Infty > SoundSpeed_Infty) {
        /*--- Supersonic outflow ---*/
        RiemannMinus = Vn_Bound - 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Subsonic outflow ---*/
        RiemannMinus = Vn_Infty - 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }

      /*--- Compute a new value for the local normal velocity and speed of
         sound from the Riemann invariants. ---*/

      Vn = 0.5 * (RiemannPlus + RiemannMinus);
      SoundSpeed = 0.25 * (RiemannPlus - RiemannMinus)*Gamma_Minus_One;

      /*--- Construct the primitive variable state at the boundary for
         computing the flux for the weak boundary condition. The values
         that we choose to construct the solution (boundary or freestream)
         depend on whether we are at an inflow or outflow. At an outflow, we
         choose boundary information (at most one characteristic is incoming),
         while at an inflow, we choose infinity values (at most one
         characteristic is outgoing). ---*/

      if (Qn_Infty > 0.0)   {
        /*--- Outflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Bound[iDim] + (Vn-Vn_Bound)*UnitNormal[iDim];
        Entropy = Entropy_Bound;
      } else  {
        /*--- Inflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Infty[iDim] + (Vn-Vn_Infty)*UnitNormal[iDim];
        Entropy = Entropy_Infty;
      }

      /*--- Recompute the primitive variables. ---*/

      Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
      Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      if (tkeNeeded) Energy += GetTke_Inf();

      /*--- Store new primitive state for computing the flux. ---*/

      V_infty[0] = Pressure/(Gas_Constant*Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_infty[iDim+1] = Velocity[iDim];
      V_infty[nDim+1] = Pressure;
      V_infty[nDim+2] = Density;
      V_infty[nDim+3] = Energy + Pressure/Density;


      
      /*--- Set various quantities in the numerics class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      if (grid_movement) {
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      }
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Convective Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous residual contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_infty[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_infty[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_infty);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                          node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                              solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update viscous residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Viscous Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container,
                              CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jVar, kVar;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double P_Total, T_Total, P_static, T_static, Rho_static, *Mach, *Flow_Dir, Area, UnitNormal[3];
  su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b, Temperature_b;
  su2double *Velocity_e, Velocity2_e, VelMag_e, Enthalpy_e, Entropy_e, Energy_e = 0.0, StaticEnthalpy_e, StaticEnergy_e, Density_e = 0.0, Pressure_e;
  su2double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
  su2double ProjVelocity_i;
  su2double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **DubDu, *dw, *u_e, *u_i, *u_b;
  su2double *gridVel;
  su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  
  su2double *Normal, *FlowDirMix, TangVelocity, NormalVelocity;
  Normal = new su2double[nDim];
  su2double ext_flow_angle;
  
  Velocity_i = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_e = new su2double[nDim];
  FlowDirMix = new su2double[nDim];
  Lambda_i = new su2double[nVar];
  u_i = new su2double[nVar];
  u_e = new su2double[nVar];
  u_b = new su2double[nVar];
  dw = new su2double[nVar];
  
  S_boundary = new su2double[8];
  
  P_Tensor = new su2double*[nVar];
  invP_Tensor = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    P_Tensor[iVar] = new su2double[nVar];
    invP_Tensor[iVar] = new su2double[nVar];
  }
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    V_boundary= GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Compute the internal state u_i ---*/
      Velocity2_i = 0;
      for (iDim=0; iDim < nDim; iDim++)
      {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      
      
      Density_i = node[iPoint]->GetDensity();
      
      Energy_i = node[iPoint]->GetEnergy();
      StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
      
      FluidModel->SetTDState_rhoe(Density_i, StaticEnergy_i);
      
      Pressure_i = FluidModel->GetPressure();
      Enthalpy_i = Energy_i + Pressure_i/Density_i;
      
      SoundSpeed_i = FluidModel->GetSoundSpeed();
      
      Kappa_i = FluidModel->GetdPde_rho() / Density_i;
      Chi_i = FluidModel->GetdPdrho_e() - Kappa_i * StaticEnergy_i;
      
      ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      
      /*--- Build the external state u_e from boundary data and internal node ---*/
      
      switch(config->GetKind_Data_Riemann(Marker_Tag))
      {
          //TODO(turbo), generilize for 3D case
          //TODO(turbo), generilize for Inlet and Outlet in for backflow treatment
          //TODO(turbo), implement not uniform inlet and radial equilibrium for the outlet
        case TOTAL_CONDITIONS_PT:
          
          /*--- Retrieve the specified total conditions for this boundary. ---*/
          if (gravity) P_Total = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
          else P_Total  = config->GetRiemann_Var1(Marker_Tag);
          T_Total  = config->GetRiemann_Var2(Marker_Tag);
          Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          P_Total /= config->GetPressure_Ref();
          T_Total /= config->GetTemperature_Ref();
          
          /*--- Computes the total state ---*/
          FluidModel->SetTDState_PT(P_Total, T_Total);
          Enthalpy_e = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
          Entropy_e = FluidModel->GetEntropy();
          
          /*--- Compute the boundary state u_e ---*/
          Velocity2_e = Velocity2_i;
          if (nDim == 2) {
            NormalVelocity= -sqrt(Velocity2_e)*Flow_Dir[0];
            TangVelocity= -sqrt(Velocity2_e)*Flow_Dir[1];
            Velocity_e[0]= UnitNormal[0]*NormalVelocity - UnitNormal[1]*TangVelocity;
            Velocity_e[1]= UnitNormal[1]*NormalVelocity + UnitNormal[0]*TangVelocity;
          }else {
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity_e[iDim] = sqrt(Velocity2_e)*Flow_Dir[iDim];
          }
          StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
          FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
          Density_e = FluidModel->GetDensity();
          StaticEnergy_e = FluidModel->GetStaticEnergy();
          Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
          if (tkeNeeded) Energy_e += GetTke_Inf();
          break;
          
        case STATIC_SUPERSONIC_INFLOW_PT:
          
          /*--- Retrieve the specified total conditions for this boundary. ---*/
          if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
          else P_static  = config->GetRiemann_Var1(Marker_Tag);
          T_static  = config->GetRiemann_Var2(Marker_Tag);
          Mach = config->GetRiemann_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          P_static /= config->GetPressure_Ref();
          T_static /= config->GetTemperature_Ref();
          
          /*--- Computes the total state ---*/
          FluidModel->SetTDState_PT(P_static, T_static);
          
          /*--- Compute the boundary state u_e ---*/
          Velocity2_e = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_e[iDim] = Mach[iDim]*FluidModel->GetSoundSpeed();
            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
          }
          Density_e = FluidModel->GetDensity();
          StaticEnergy_e = FluidModel->GetStaticEnergy();
          Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
          if (tkeNeeded) Energy_e += GetTke_Inf();
          break;
          
        case STATIC_SUPERSONIC_INFLOW_PD:
          
          /*--- Retrieve the specified total conditions for this boundary. ---*/
          
          if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
          else P_static  = config->GetRiemann_Var1(Marker_Tag);
          Rho_static  = config->GetRiemann_Var2(Marker_Tag);
          Mach = config->GetRiemann_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          P_static /= config->GetPressure_Ref();
          Rho_static /= config->GetDensity_Ref();
          
          /*--- Computes the total state ---*/
          FluidModel->SetTDState_Prho(P_static, Rho_static);
          
          /*--- Compute the boundary state u_e ---*/
          Velocity2_e = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_e[iDim] = Mach[iDim]*FluidModel->GetSoundSpeed();
            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
          }
          Density_e = FluidModel->GetDensity();
          StaticEnergy_e = FluidModel->GetStaticEnergy();
          Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
          if (tkeNeeded) Energy_e += GetTke_Inf();
          break;
          
        case MIXING_IN:
          
          /*--- Retrieve the specified total conditions for this boundary. ---*/
          P_Total = ExtAveragedTotPressure[val_marker];
          T_Total = ExtAveragedTotTemperature[val_marker];
          ext_flow_angle = atan(ExtAveragedTangVelocity[val_marker]/ExtAveragedNormalVelocity[val_marker]);
          FlowDirMix[0] = cos(ext_flow_angle);
          FlowDirMix[1] = sin(ext_flow_angle);
          
          /*--- Computes the total state ---*/
          FluidModel->SetTDState_PT(P_Total, T_Total);
          Enthalpy_e = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
          Entropy_e = FluidModel->GetEntropy();
          
          /*--- Compute the boundary state u_e ---*/
          Velocity2_e = Velocity2_i;
          if (nDim == 2) {
            NormalVelocity= -sqrt(Velocity2_e)*FlowDirMix[0];
            TangVelocity= -sqrt(Velocity2_e)*FlowDirMix[1];
            Velocity_e[0]= UnitNormal[0]*NormalVelocity - UnitNormal[1]*TangVelocity;
            Velocity_e[1]= UnitNormal[1]*NormalVelocity + UnitNormal[0]*TangVelocity;
          }else {
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity_e[iDim] = sqrt(Velocity2_e)*FlowDirMix[iDim];
          }
          StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
          FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
          Density_e = FluidModel->GetDensity();
          StaticEnergy_e = FluidModel->GetStaticEnergy();
          Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
          if (tkeNeeded) Energy_e += GetTke_Inf();
          break;
          
        case DENSITY_VELOCITY:
          
          /*--- Retrieve the specified density and velocity magnitude ---*/
          Density_e  = config->GetRiemann_Var1(Marker_Tag);
          VelMag_e   = config->GetRiemann_Var2(Marker_Tag);
          Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          Density_e /= config->GetDensity_Ref();
          VelMag_e /= config->GetVelocity_Ref();
          
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity_e[iDim] = VelMag_e*Flow_Dir[iDim];
          Energy_e = Energy_i;
          break;
          
        case MIXING_OUT:
          
          /*--- Retrieve the staic pressure for this boundary. ---*/
          Pressure_e = ExtAveragedPressure[val_marker];
          Density_e = Density_i;
          
          /*--- Compute the boundary state u_e ---*/
          FluidModel->SetTDState_Prho(Pressure_e, Density_e);
          Velocity2_e = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_e[iDim] = Velocity_i[iDim];
            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
          }
          Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
          break;
          
        case STATIC_PRESSURE:
          
          /*--- Retrieve the staic pressure for this boundary. ---*/
          Pressure_e = config->GetRiemann_Var1(Marker_Tag);
          Pressure_e /= config->GetPressure_Ref();
          Density_e = Density_i;
          
          /*--- Compute the boundary state u_e ---*/
          FluidModel->SetTDState_Prho(Pressure_e, Density_e);
          Velocity2_e = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_e[iDim] = Velocity_i[iDim];
            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
          }
          Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
          break;
          
        default:
          cout << "Warning! Invalid Riemann input!" << endl;
          exit(EXIT_FAILURE);
          break;
          
      }
      
      /*--- Compute P (matrix of right eigenvectors) ---*/
      conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);
      
      /*--- Compute inverse P (matrix of left eigenvectors)---*/
      conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);
      
      /*--- eigenvalues contribution due to grid motion ---*/
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        
        su2double ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel   += gridVel[iDim]*UnitNormal[iDim];
        ProjVelocity_i -= ProjGridVel;
      }
      
      /*--- Flow eigenvalues ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Lambda_i[iDim] = ProjVelocity_i;
      Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i;
      Lambda_i[nVar-1] = ProjVelocity_i - SoundSpeed_i;
      
      /*--- Compute the boundary state u_e ---*/
      u_e[0] = Density_e;
      for (iDim = 0; iDim < nDim; iDim++)
        u_e[iDim+1] = Velocity_e[iDim]*Density_e;
      u_e[nVar-1] = Energy_e*Density_e;
      
      /*--- Compute the boundary state u_i ---*/
      u_i[0] = Density_i;
      for (iDim = 0; iDim < nDim; iDim++)
        u_i[iDim+1] = Velocity_i[iDim]*Density_i;
      u_i[nVar-1] = Energy_i*Density_i;
      
      /*--- Compute the characteristic jumps ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dw[iVar] = 0;
        for (jVar = 0; jVar < nVar; jVar++)
          dw[iVar] += invP_Tensor[iVar][jVar] * (u_e[jVar] - u_i[jVar]);
        
      }
      
      /*--- Compute the boundary state u_b using characteristics ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        u_b[iVar] = u_i[iVar];
        
        for (jVar = 0; jVar < nVar; jVar++)
        {
          if (Lambda_i[jVar] < 0)
          {
            u_b[iVar] += P_Tensor[iVar][jVar]*dw[jVar];
            
          }
        }
      }
      
      
      /*--- Compute the thermodynamic state in u_b ---*/
      Density_b = u_b[0];
      Velocity2_b = 0;
      for (iDim = 0; iDim < nDim; iDim++)
      {
        Velocity_b[iDim] = u_b[iDim+1]/Density_b;
        Velocity2_b += Velocity_b[iDim]*Velocity_b[iDim];
      }
      Energy_b = u_b[nVar-1]/Density_b;
      StaticEnergy_b = Energy_b - 0.5*Velocity2_b;
      FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);
      Pressure_b = FluidModel->GetPressure();
      Temperature_b = FluidModel->GetTemperature();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;
      Kappa_b = FluidModel->GetdPde_rho() / Density_b;
      Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
      
      /*--- Compute the residuals ---*/
      conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);
      
      /*--- Residual contribution due to grid motion ---*/
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        su2double projVelocity = 0.0;
        
        for (iDim = 0; iDim < nDim; iDim++)
          projVelocity +=  gridVel[iDim]*Normal[iDim];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] -= projVelocity *(u_b[iVar]);
      }
      
      if (implicit) {
        
        Jacobian_b = new su2double*[nVar];
        DubDu = new su2double*[nVar];
        for (iVar = 0; iVar < nVar; iVar++)
        {
          Jacobian_b[iVar] = new su2double[nVar];
          DubDu[iVar] = new su2double[nVar];
        }
        
        /*--- Initialize DubDu to unit matrix---*/
        
        for (iVar = 0; iVar < nVar; iVar++)
        {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0;
          
          DubDu[iVar][iVar]= 1;
        }
        
        /*--- Compute DubDu -= RNL---*/
        for (iVar=0; iVar<nVar; iVar++)
        {
          for (jVar=0; jVar<nVar; jVar++)
          {
            for (kVar=0; kVar<nVar; kVar++)
            {
              if (Lambda_i[kVar]<0)
                DubDu[iVar][jVar] -= P_Tensor[iVar][kVar] * invP_Tensor[kVar][jVar];
            }
          }
        }
        
        /*--- Compute flux Jacobian in state b ---*/
        conv_numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, Normal, 1.0, Jacobian_b);
        
        /*--- Jacobian contribution due to grid motion ---*/
        if (grid_movement)
        {
          gridVel = geometry->node[iPoint]->GetGridVel();
          su2double projVelocity = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            projVelocity +=  gridVel[iDim]*Normal[iDim];
          for (iVar = 0; iVar < nVar; iVar++) {
            Residual[iVar] -= projVelocity *(u_b[iVar]);
            Jacobian_b[iVar][iVar] -= projVelocity;
          }
          
        }
        
        /*--- initiate Jacobian_i to zero matrix ---*/
        for (iVar=0; iVar<nVar; iVar++)
          for (jVar=0; jVar<nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        
        /*--- Compute numerical flux Jacobian at node i ---*/
        for (iVar=0; iVar<nVar; iVar++) {
          for (jVar=0; jVar<nVar; jVar++) {
            for (kVar=0; kVar<nVar; kVar++) {
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
            }
          }
        }
        
        for (iVar = 0; iVar < nVar; iVar++) {
          delete [] Jacobian_b[iVar];
          delete [] DubDu[iVar];
        }
        delete [] Jacobian_b;
        delete [] DubDu;
      }
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Primitive variables, using the derived quantities ---*/
        V_boundary[0] = Temperature_b;
        for (iDim = 0; iDim < nDim; iDim++)
          V_boundary[iDim+1] = Velocity_b[iDim];
        V_boundary[nDim+1] = Pressure_b;
        V_boundary[nDim+2] = Density_b;
        V_boundary[nDim+3] = Enthalpy_b;
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_boundary[nDim+5] = FluidModel->GetLaminarViscosity();
        V_boundary[nDim+6] = node[iPoint]->GetEddyViscosity();
        V_boundary[nDim+7] = FluidModel->GetThermalConductivity();
        V_boundary[nDim+8] = FluidModel->GetCp();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_boundary);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Secondary variables ---*/
        S_domain = node[iPoint]->GetSecondary();
        
        /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
        
        S_boundary[0]= FluidModel->GetdPdrho_e();
        S_boundary[1]= FluidModel->GetdPde_rho();
        
        S_boundary[2]= FluidModel->GetdTdrho_e();
        S_boundary[3]= FluidModel->GetdTde_rho();
        
        /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
        
        S_boundary[4]= FluidModel->Getdmudrho_T();
        S_boundary[5]= FluidModel->GetdmudT_rho();
        
        S_boundary[6]= FluidModel->Getdktdrho_T();
        S_boundary[7]= FluidModel->GetdktdT_rho();
        
        visc_numerics->SetSecondary(S_domain, S_boundary);
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  delete [] Velocity_e;
  delete [] Velocity_b;
  delete [] Velocity_i;
  delete [] FlowDirMix;
  
  delete [] S_boundary;
  delete [] Lambda_i;
  delete [] u_i;
  delete [] u_e;
  delete [] u_b;
  delete [] dw;
  
  
  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}


void CEulerSolver::Mixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker) {
  
  unsigned long iVertex, iPoint, nVert;
  unsigned short iDim, iVar;
  unsigned short mixing_process = config->GetKind_MixingProcess();
  su2double Pressure = 0.0, Density = 0.0, Enthalpy = 0.0,  *Velocity = NULL, *Normal, *gridVel,
  Area, TotalArea, TotalAreaPressure, TotalAreaDensity, *TotalAreaVelocity, UnitNormal[3];
  string Marker_Tag, Monitoring_Tag;
  su2double val_init_pressure;
  bool grid_movement        = config->GetGrid_Movement();
  su2double TotalDensity, TotalPressure, *TotalVelocity, TotalNormal, avgVel2, avgTotalEnthaply;
  
  /*-- Variables declaration and allocation ---*/
  Velocity = new su2double[nDim];
  Normal = new su2double[nDim];
  TotalVelocity = new su2double[nDim];
  TotalAreaVelocity = new su2double[nDim];
  
  for (iDim=0; iDim<nDim; iDim++) {
    TotalVelocity[iDim]=0;
    TotalAreaVelocity[iDim]=0;
  }
  
  TotalDensity = 0.0;
  TotalPressure = 0.0;
  TotalAreaPressure=0.0;
  TotalAreaDensity=0.0;
  TotalArea = 0.0;
  TotalNormal=0.0;
  
  /*--- Forces initialization for Marker vector ---*/
  AveragedPressure[val_Marker] = 0.0;
  AveragedEnthalpy[val_Marker] = 0.0;
  AveragedDensity[val_Marker] = 0.0;
  AveragedSoundSpeed[val_Marker] = 0.0;
  
  for (iDim=0;iDim < nDim;iDim++) {
    AveragedVelocity[val_Marker][iDim] = 0.0;
    AveragedNormal[val_Marker][iDim] = 0.0;
    AveragedGridVel[val_Marker][iDim] = 0.0;
  }
  
  for (iVar=0;iVar<nVar;iVar++)
    TotalFlux[val_Marker][iVar]= 0.0;
  
  /*--- Loop over the vertices to compute the averaged quantities ---*/
  nVert = 0;
  for (iVertex = 0; iVertex < geometry->GetnVertex(val_Marker); iVertex++) {
    
    iPoint = geometry->vertex[val_Marker][iVertex]->GetNode();
    
    /*--- Compute the integral fluxes for the boundaries ---*/
    Pressure = node[iPoint]->GetPressure();
    Density = node[iPoint]->GetDensity();
    Enthalpy = node[iPoint]->GetEnthalpy();
    
    /*--- Note that the fluxes from halo cells are discarded ---*/
    if ( (geometry->node[iPoint]->GetDomain())  ) {
      nVert++;
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_Marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      su2double VelNormal = 0.0, VelSq = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        UnitNormal[iDim] = Normal[iDim]/Area;
        Velocity[iDim] = node[iPoint]->GetPrimitive(iDim+1);
        VelNormal += UnitNormal[iDim]*Velocity[iDim];
        VelSq += Velocity[iDim]*Velocity[iDim];
      }
      
      
      /*--- Compute the integral fluxes for the boundary of interest ---*/
      
      if ((mixing_process == AREA_AVERAGE) || (mixing_process == MIXEDOUT_AVERAGE)) {
        
        TotalFlux[val_Marker][0] += Area*(Density*VelNormal );
        for (iDim = 1; iDim < nDim+1; iDim++)
          TotalFlux[val_Marker][iDim] += Area*(Density*VelNormal*Velocity[iDim -1] + Pressure*UnitNormal[iDim -1] );
        TotalFlux[val_Marker][nDim+1] += Area*(Density*VelNormal*Enthalpy );
        
        TotalArea += Area;
        TotalAreaPressure += Area*Pressure;
        TotalAreaDensity  += Area*Density;
        for (iDim = 0; iDim < nDim; iDim++)
          TotalAreaVelocity[iDim] += Area*Velocity[iDim];
        
      }else {
        
        TotalDensity += Density;
        TotalPressure += Pressure;
        for (iDim = 0; iDim < nDim; iDim++)
          TotalVelocity[iDim] += Velocity[iDim];
        
        
      }
      for (iDim = 0; iDim < nDim; iDim++) AveragedNormal[val_Marker][iDim] +=Normal[iDim];
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          AveragedGridVel[val_Marker][iDim] +=gridVel[iDim];
      }
    }
  }
  
  /*--- Compute the averaged value for the boundary of interest ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    AveragedNormal[val_Marker][iDim] /=nVert;
    TotalNormal+= AveragedNormal[val_Marker][iDim]*AveragedNormal[val_Marker][iDim];
  }
  for (iDim = 0; iDim < nDim; iDim++) AveragedNormal[val_Marker][iDim] /=sqrt(TotalNormal);
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++)
      AveragedGridVel[val_Marker][iDim] /=nVert;
  }
  
  switch(mixing_process) {
      
      
    case ALGEBRAIC_AVERAGE:
      AveragedDensity[val_Marker] = TotalDensity / nVert;
      AveragedPressure[val_Marker] = TotalPressure / nVert;
      for (iDim = 0; iDim < nDim; iDim++)
        AveragedVelocity[val_Marker][iDim] = TotalVelocity[iDim] / nVert;
      break;
      
    case AREA_AVERAGE:
      AveragedDensity[val_Marker] = TotalAreaDensity / TotalArea;
      AveragedPressure[val_Marker] = TotalAreaPressure / TotalArea;
      for (iDim = 0; iDim < nDim; iDim++)
        AveragedVelocity[val_Marker][iDim] = TotalAreaVelocity[iDim] / TotalArea;
      break;
      
    case MIXEDOUT_AVERAGE:
      for (iVar = 0; iVar<nVar; iVar++) {
        AveragedFlux[val_Marker][iVar] = TotalFlux[val_Marker][iVar]/TotalArea;
      }
      val_init_pressure = TotalAreaPressure/TotalArea;
      
      if (abs(AveragedFlux[val_Marker][0])<(10.0e-9)*TotalAreaDensity) {
        cout << "Mass flux is 0.0 so a Area Averaged algorithm is used for the Mixing Procees" << endl;
        AveragedDensity[val_Marker] = TotalAreaDensity / TotalArea;
        AveragedPressure[val_Marker] = TotalAreaPressure / TotalArea;
        for (iDim = 0; iDim < nDim; iDim++)
          AveragedVelocity[val_Marker][iDim] = TotalAreaVelocity[iDim] / TotalArea;
        
      }else {
        MixedOut_Average (val_init_pressure, AveragedFlux[val_Marker], AveragedNormal[val_Marker], &AveragedPressure[val_Marker], &AveragedDensity[val_Marker]);
        for (iDim = 1; iDim < nDim +1;iDim++)
          AveragedVelocity[val_Marker][iDim-1]= ( AveragedFlux[val_Marker][iDim] - AveragedPressure[val_Marker]*AveragedNormal[val_Marker][iDim-1] ) / AveragedFlux[val_Marker][0];
      }
      break;
      
      
    default:
      cout << "Warning! Invalid MIXING_PROCESS input!" << endl;
      exit(EXIT_FAILURE);
      break;
  }
  
  /*--- compute static averaged quantities ---*/
  FluidModel->SetTDState_Prho(AveragedPressure[val_Marker], AveragedDensity[val_Marker]);
  AveragedEnthalpy[val_Marker] = FluidModel->GetStaticEnergy() + AveragedPressure[val_Marker]/AveragedDensity[val_Marker];
  AveragedSoundSpeed[val_Marker] = FluidModel->GetSoundSpeed();
  AveragedEntropy[val_Marker] = FluidModel->GetEntropy();
  AveragedNormalVelocity[val_Marker]= AveragedNormal[val_Marker][0]*AveragedVelocity[val_Marker][0] + AveragedNormal[val_Marker][1]*AveragedVelocity[val_Marker][1];
  AveragedTangVelocity[val_Marker]= AveragedNormal[val_Marker][0]*AveragedVelocity[val_Marker][1] - AveragedNormal[val_Marker][1]*AveragedVelocity[val_Marker][0];
  MassFlow[val_Marker]= AveragedDensity[val_Marker]*AveragedNormalVelocity[val_Marker]*TotalArea;
  FlowAngle[val_Marker]= atan(AveragedTangVelocity[val_Marker]/AveragedNormalVelocity[val_Marker]);
  
  /*--- compute total averaged quantities ---*/
  avgVel2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) avgVel2 += AveragedVelocity[val_Marker][iDim]*AveragedVelocity[val_Marker][iDim];
  
  avgTotalEnthaply = AveragedEnthalpy[val_Marker] + 0.5*avgVel2;
  FluidModel->SetTDState_hs(avgTotalEnthaply,AveragedEntropy[val_Marker]);
  AveragedTotTemperature[val_Marker] = FluidModel->GetTemperature();
  AveragedTotPressure[val_Marker] = FluidModel->GetPressure();
  
  if(grid_movement) {
    AveragedTangGridVelocity[val_Marker] = AveragedNormal[val_Marker][0]*AveragedGridVel[val_Marker][1]-AveragedNormal[val_Marker][1]*AveragedGridVel[val_Marker][0];
    AveragedMach[val_Marker] = sqrt(AveragedNormalVelocity[val_Marker]*AveragedNormalVelocity[val_Marker] + (AveragedTangVelocity[val_Marker] - AveragedTangGridVelocity[val_Marker])*(AveragedTangVelocity[val_Marker] - AveragedTangGridVelocity[val_Marker]));
    AveragedMach[val_Marker] /= AveragedSoundSpeed[val_Marker];
    AveragedTangMach[val_Marker] = (AveragedTangVelocity[val_Marker] - AveragedTangGridVelocity[val_Marker])/AveragedSoundSpeed[val_Marker];
    FlowAngle[val_Marker]= atan((AveragedTangVelocity[val_Marker] - AveragedTangGridVelocity[val_Marker])/AveragedNormalVelocity[val_Marker]);
    
  }else {
    AveragedMach[val_Marker] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      AveragedMach[val_Marker] += AveragedVelocity[val_Marker][iDim]*AveragedVelocity[val_Marker][iDim];
    }
    AveragedMach[val_Marker] = sqrt(AveragedMach[val_Marker])/AveragedSoundSpeed[val_Marker];
    AveragedTangMach[val_Marker] = AveragedTangVelocity[val_Marker]/AveragedSoundSpeed[val_Marker];
    
  }
  
  AveragedNormalMach[val_Marker] = AveragedNormalVelocity[val_Marker]/AveragedSoundSpeed[val_Marker];
  
  
  if ((AveragedDensity[val_Marker] != AveragedDensity[val_Marker]) || (AveragedEnthalpy[val_Marker] !=AveragedEnthalpy[val_Marker]))
    cout<<"nan in mixing process in boundary "<<config->GetMarker_All_TagBound(val_Marker)<< endl;
  
  /*--- Free locally allocated memory ---*/
  delete [] Velocity;
  delete [] Normal;
  delete [] TotalVelocity;
  delete [] TotalAreaVelocity;
}

void CEulerSolver::MixedOut_Average (su2double val_init_pressure, su2double *val_Averaged_Flux, su2double *val_normal,
                                     su2double *pressure_mix, su2double *density_mix) {
  
  unsigned short maxiter = 10;
  unsigned short iter = 0;
  su2double toll = 1.0e-07;
  su2double resdl = 0.0;
  
  su2double *val_func = new su2double, *val_right_func = new su2double, *val_left_func = new su2double;
  su2double deltaP, *p_mix = new su2double, *p_mix_right = new su2double, *p_mix_left = new su2double;
  su2double epsilon = 1.0e-04;
  su2double relax_factor = 1;
  
  *pressure_mix = val_init_pressure;
  
  /*--- Newton-Raphson's method with central difference formula ---*/
  
  while ( iter <= maxiter ) {
    deltaP = 2*epsilon*(*pressure_mix);
    *p_mix_right = *pressure_mix+deltaP/2;
    *p_mix_left = *pressure_mix-deltaP/2;
    *p_mix = *pressure_mix;
    MixedOut_Root_Function(p_mix_right,val_Averaged_Flux,val_normal,val_right_func,density_mix);
    MixedOut_Root_Function(p_mix_left,val_Averaged_Flux,val_normal,val_left_func,density_mix);
    MixedOut_Root_Function(p_mix,val_Averaged_Flux,val_normal,val_func,density_mix);
    su2double der_func = (*val_right_func-*val_left_func) / deltaP;
    deltaP = -*val_func/der_func;
    resdl = deltaP/val_init_pressure;
    *pressure_mix += relax_factor*(deltaP);
    
    iter += 1;
    if ( abs(resdl) <= toll ) {
      break;
    }
    
  }
  
  MixedOut_Root_Function(pressure_mix,val_Averaged_Flux,val_normal,val_func,density_mix);
  
  /*--- Free locally allocated memory ---*/
  delete val_func;
  delete val_right_func;
  delete val_left_func;
  delete p_mix;
  delete p_mix_right;
  delete p_mix_left;
  
}

void CEulerSolver::MixedOut_Root_Function(su2double *pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double *valfunc, su2double *density) {
  
  su2double velnormal, velsq;
  
  su2double *vel;
  vel = new su2double[nDim];
  
  
  *valfunc = 0.0;
  *density = 0.0;
  
  velnormal = 0.0;
  velsq = 0.0;
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    vel[iDim]  = (val_Averaged_Flux[iDim+1] - (*pressure)*val_normal[iDim]) / val_Averaged_Flux[0];
    velnormal += val_normal[iDim]*vel[iDim];
    velsq += vel[iDim]*vel[iDim];
  }
  *density = val_Averaged_Flux[0] / velnormal;
  if (*density <= 0) cout << " desnity in mixedout routine negative : " << endl;
  FluidModel->SetTDState_Prho(*pressure, *density);
  su2double enthalpy = FluidModel->GetStaticEnergy() + (*pressure)/(*density);
  *valfunc = val_Averaged_Flux[nDim+1]/val_Averaged_Flux[0] - enthalpy - velsq/2;
  if (*valfunc!=*valfunc) cout << " mixedout root func gives nan: " << endl;
  
  
  /*--- Free locally allocated memory ---*/
  delete [] vel;
  
}

void CEulerSolver::Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker, vector<std::complex<su2double> >& c4k,signed long& nboundaryvertex) {
  /* Implementation of Fuorier Transformations for non-regfelcting BC will come soon */
}

void CEulerSolver::Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker, vector<std::complex<su2double> >& c2k,vector<std::complex<su2double> >& c3k,signed long& nboundaryvertex) {
  /* Implementation of Fuorier Transformations for non-regfelcting BC will come soon */
}

void CEulerSolver::BC_NonReflecting(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jVar, kVar;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double  Area, UnitNormal[3];
  
  su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b, Temperature_b;
  su2double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
  su2double Pressure_e;
  su2double ProjVelocity_i;
  su2double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **DubDu, *dw, *u_b;
  su2double *gridVel;
  su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  su2double *Normal;
  
  Normal = new su2double[nDim];
  
  Velocity_i = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  
  
  Lambda_i = new su2double[nVar];
  
  u_b = new su2double[nVar];
  dw = new su2double[nVar];
  
  S_boundary = new su2double[8];
  
  P_Tensor = new su2double*[nVar];
  invP_Tensor = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    P_Tensor[iVar] = new su2double[nVar];
    invP_Tensor[iVar] = new su2double[nVar];
  }
  
  /*--- new declarations ---*/
  std::vector<std::complex<su2double> > c4k ;//    std::complex<su2double> c3k[nVertex-OddEven]=0;
  std::vector<std::complex<su2double> > c2k ;//    std::complex<su2double> c3k[nVertex-OddEven]=0;
  std::vector<std::complex<su2double> > c3k ;//    std::complex<su2double> c3k[nVertex-OddEven]=0;
  
  su2double  deltaDensity, deltaPressure, AvgMach, deltaTangVelocity, deltaNormalVelocity, cc,rhoc,c1j,c2j,c3j,c4j,
  avg_c1, avg_c2, avg_c3, avg_c4,TangVelocity, NormalVelocity, GilesBeta, c4js, dc4js, *delta_c, **R_Matrix, *deltaprim;
  
  
  delta_c = new su2double[nVar];
  deltaprim = new su2double[nVar];
  R_Matrix= new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    R_Matrix[iVar] = new su2double[nVar];
  }
  
  
  Mixing_Process(geometry, solver_container,  config, val_marker);
  
  cc = AveragedSoundSpeed[val_marker]*AveragedSoundSpeed[val_marker];
  rhoc = AveragedSoundSpeed[val_marker]*AveragedDensity[val_marker];
  AvgMach = AveragedMach[val_marker];
  
  conv_numerics->GetRMatrix(AveragedSoundSpeed[val_marker], AveragedDensity[val_marker], AveragedNormal[val_marker], R_Matrix);
  
  //  Boundary_Fourier(geometry, solver_container, config, val_marker, c4k, nboundaryvertex);
  //  Boundary_Fourier(geometry, solver_container, config, val_marker, c2k,c3k,nboundaryvertex);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    V_boundary= GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Compute the internal state u_i ---*/
      Velocity2_i = 0;
      for (iDim = 0; iDim < nDim; iDim++)
      {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      
      
      Density_i = node[iPoint]->GetDensity();
      
      Energy_i = node[iPoint]->GetEnergy();
      StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
      
      FluidModel->SetTDState_rhoe(Density_i, StaticEnergy_i);
      
      Pressure_i = FluidModel->GetPressure();
      Enthalpy_i = Energy_i + Pressure_i/Density_i;
      
      SoundSpeed_i = FluidModel->GetSoundSpeed();
      
      Kappa_i = FluidModel->GetdPde_rho() / Density_i;
      Chi_i = FluidModel->GetdPdrho_e() - Kappa_i * StaticEnergy_i;
      
      ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      
      
      switch(config->GetKind_Data_NRBC(Marker_Tag))
      {
          
          //TODO(turbo), generilize for 3D case
          //TODO(turbo), generilize for Inlet and Outlet in for backflow treatment
          //TODO(turbo), implement not uniform inlet and radial equilibrium for the outlet
          
        case MIXING_IN:
          
          /*--- Compute jump of primitive variable  ---*/
          deltaDensity = ExtAveragedDensity[val_marker] - AveragedDensity[val_marker];
          deltaPressure = ExtAveragedPressure[val_marker] - AveragedPressure[val_marker];
          NormalVelocity= UnitNormal[0]*Velocity_i[0] + UnitNormal[1]*Velocity_i[1];
          deltaTangVelocity= ExtAveragedTangVelocity[val_marker]+AveragedTangVelocity[val_marker];
          deltaNormalVelocity= ExtAveragedNormalVelocity[val_marker]+AveragedNormalVelocity[val_marker];
          
          /*--- Compute characteristic jumps  ---*/
          avg_c1= -cc*deltaDensity +deltaPressure;
          avg_c2= (rhoc*deltaTangVelocity);
          avg_c3= (rhoc*deltaNormalVelocity +deltaPressure);
          c4j= -rhoc*(-NormalVelocity +AveragedNormalVelocity[val_marker]) +(Pressure_i - AveragedPressure[val_marker]);
          
          /*--- Impose Inlet BC  ---*/
          delta_c[0] = avg_c1;
          delta_c[1] = avg_c2;
          delta_c[2] = avg_c3;
          delta_c[3] = c4j;
          break;
          
        case MIXING_OUT:
          
          /*--- Compute jump of primitive variable  ---*/
          deltaDensity = Density_i - AveragedDensity[val_marker];
          deltaPressure = Pressure_i - AveragedPressure[val_marker];
          TangVelocity= UnitNormal[0]*Velocity_i[1] - UnitNormal[1]*Velocity_i[0];
          NormalVelocity= UnitNormal[0]*Velocity_i[0] + UnitNormal[1]*Velocity_i[1];
          deltaTangVelocity= TangVelocity - AveragedTangVelocity[val_marker];
          deltaNormalVelocity= NormalVelocity - AveragedNormalVelocity[val_marker];
          
          /*--- Compute characteristic jumps  ---*/
          c1j= -cc*deltaDensity +deltaPressure;
          c2j= rhoc*deltaTangVelocity;
          c3j= rhoc*deltaNormalVelocity + deltaPressure;
          avg_c4 = rhoc*(AveragedNormalVelocity[val_marker]+ExtAveragedNormalVelocity[val_marker]) -(AveragedPressure[val_marker]-ExtAveragedPressure[val_marker]);
          
          /*--- implementation of supersonic NRBC ---*/
          if (AvgMach > 1.001) {
            if (AveragedTangVelocity[val_marker] >= 0.0) {
              GilesBeta= -sqrt(pow(AvgMach,2)-1.0);
            }else {
              GilesBeta= sqrt(pow(AvgMach,2)-1.0);
            }
            c4js= (2.0 * AveragedNormalMach[val_marker])/(GilesBeta - AveragedTangMach[val_marker])*c2j - (GilesBeta+AveragedTangMach[val_marker])/(GilesBeta-AveragedTangMach[val_marker])*c3j;
            dc4js = c4js;
          }else {
            dc4js = 0.0;
          }
          
          /*--- Impose Outlet BC  ---*/
          delta_c[0] = c1j;
          delta_c[1] = c2j;
          delta_c[2] = c3j;
          delta_c[3] = avg_c4 + dc4js;
          break;
          
        case STATIC_PRESSURE:
          
          Pressure_e = config->GetNRBC_Var1(Marker_Tag);
          Pressure_e /= config->GetPressure_Ref();
          
          /*--- Compute jump of primitive variable  ---*/
          deltaDensity = Density_i - AveragedDensity[val_marker];
          deltaPressure = Pressure_i - AveragedPressure[val_marker];
          TangVelocity= UnitNormal[0]*Velocity_i[1] - UnitNormal[1]*Velocity_i[0];
          NormalVelocity= UnitNormal[0]*Velocity_i[0] + UnitNormal[1]*Velocity_i[1];
          deltaTangVelocity= TangVelocity - AveragedTangVelocity[val_marker];
          deltaNormalVelocity= NormalVelocity - AveragedNormalVelocity[val_marker];
          
          /*--- Compute characteristic jumps  ---*/
          c1j= -cc*deltaDensity +deltaPressure;
          c2j= rhoc*deltaTangVelocity;
          c3j=rhoc*deltaNormalVelocity + deltaPressure;
          c4j=-rhoc*deltaNormalVelocity + deltaPressure;
          avg_c4 = -2.0*(AveragedPressure[val_marker]-Pressure_e);
          
          /*--- implementation of supersonic NRBC ---*/
          if (AvgMach > 1.001) {
            if (AveragedTangVelocity[val_marker] >= 0.0) {
              GilesBeta= -sqrt(pow(AvgMach,2)-1.0);
            }else {
              GilesBeta= sqrt(pow(AvgMach,2)-1.0);
            }
            c4js= (2.0 * AveragedNormalMach[val_marker])/(GilesBeta - AveragedTangMach[val_marker])*c2j - (GilesBeta+AveragedTangMach[val_marker])/(GilesBeta-AveragedTangMach[val_marker])*c3j;
            dc4js = c4js;
          }else {
            dc4js = 0.0;
          }
          
          /*--- Impose Outlet BC  ---*/
          delta_c[0] = c1j;
          delta_c[1] = c2j;
          delta_c[2] = c3j;
          delta_c[3] = avg_c4 + dc4js;
          break;
          
        default:
          cout << "Warning! Invalid NRBC input!" << endl;
          exit(EXIT_FAILURE);
          break;
          
      }
      
      /*--- Compute primitive jump from characteristic variables  ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        deltaprim[iVar]=0;
        for (jVar = 0; jVar < nVar; jVar++)
        {
          deltaprim[iVar] +=  R_Matrix[iVar][jVar]*delta_c[jVar];
        }
      }
      
      /*--- Compute P (matrix of right eigenvectors) ---*/
      conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);
      
      /*--- Compute inverse P (matrix of left eigenvectors)---*/
      conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);
      
      /*--- eigenvalues contribution due to grid motion ---*/
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        su2double ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel   += gridVel[iDim]*UnitNormal[iDim];
        ProjVelocity_i -= ProjGridVel;
      }
      
      /*--- Flow eigenvalues ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Lambda_i[iDim] = ProjVelocity_i;
      Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i;
      Lambda_i[nVar-1] = ProjVelocity_i - SoundSpeed_i;
      
      //TODO(turbo), provide the under relaxation factor sigma from cfg file
      su2double sigma;
      sigma = 1.0;
      
      /*--- retrieve boundary variables ---*/
      Density_b = AveragedDensity[val_marker] + sigma*deltaprim[0];
      Pressure_b = AveragedPressure[val_marker] + sigma*deltaprim[3];
      switch(config->GetKind_Data_NRBC(Marker_Tag)) {
        case MIXING_IN:
          NormalVelocity = AveragedNormalVelocity[val_marker] - sigma*deltaprim[1];
          TangVelocity = AveragedTangVelocity[val_marker] - sigma*deltaprim[2];
          break;
        case MIXING_OUT: case STATIC_PRESSURE:
          NormalVelocity = AveragedNormalVelocity[val_marker] + sigma*deltaprim[1];
          TangVelocity = AveragedTangVelocity[val_marker] + sigma*deltaprim[2];
          break;
        default:
          cout << "Warning! Invalid NRBC input!" << endl;
          exit(EXIT_FAILURE);
          break;
      }
      
      Velocity_b[0] = NormalVelocity*UnitNormal[0] - TangVelocity*UnitNormal[1];
      Velocity_b[1] = NormalVelocity*UnitNormal[1] + TangVelocity*UnitNormal[0];
      Velocity2_b = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2_b+= Velocity_b[iDim]*Velocity_b[iDim];
      }
      
      FluidModel->SetTDState_Prho(Pressure_b, Density_b);
      Energy_b = FluidModel->GetStaticEnergy() + 0.5*Velocity2_b;
      StaticEnergy_b = FluidModel->GetStaticEnergy();
      Temperature_b= FluidModel->GetTemperature();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;
      Kappa_b = FluidModel->GetdPde_rho() / Density_b;
      Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
      
      /*--- Compute the thermodynamic state in u_b ---*/
      u_b[0]=Density_b;
      u_b[1]=Density_b*Velocity_b[0];
      u_b[2]=Density_b*Velocity_b[1];
      u_b[3]=Energy_b*Density_b;
      
      /*--- Compute the residuals ---*/
      conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);
      
      /*--- Residual contribution due to grid motion ---*/
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        su2double projVelocity = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          projVelocity +=  gridVel[iDim]*Normal[iDim];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] -= projVelocity *(u_b[iVar]);
      }
      
      if (implicit) {
        /*--- Residual contribution due to grid motion ---*/
        Jacobian_b = new su2double*[nVar];
        DubDu = new su2double*[nVar];
        for (iVar = 0; iVar < nVar; iVar++)
        {
          Jacobian_b[iVar] = new su2double[nVar];
          DubDu[iVar] = new su2double[nVar];
        }
        
        /*--- Initialize DubDu to unit matrix---*/
        for (iVar = 0; iVar < nVar; iVar++)
        {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0;
          
          DubDu[iVar][iVar]= 1;
        }
        
        /*--- Compute DubDu -= RNL---*/
        for (iVar=0; iVar<nVar; iVar++)
        {
          for (jVar=0; jVar<nVar; jVar++)
          {
            for (kVar=0; kVar<nVar; kVar++)
            {
              if (Lambda_i[kVar]<0)
                DubDu[iVar][jVar] -= P_Tensor[iVar][kVar] * invP_Tensor[kVar][jVar];
            }
          }
        }
        
        /*--- Compute flux Jacobian in state b ---*/
        conv_numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, Normal, 1.0, Jacobian_b);
        
        /*--- Jacobian contribution due to grid motion ---*/
        if (grid_movement)
        {
          su2double projVelocity = 0.0;
          gridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            projVelocity +=  gridVel[iDim]*Normal[iDim];
          for (iVar = 0; iVar < nVar; iVar++) {
            Residual[iVar] -= projVelocity *(u_b[iVar]);
            Jacobian_b[iVar][iVar] -= projVelocity;
          }
          
        }
        
        /*--- initiate Jacobian_i to zero matrix ---*/
        for (iVar=0; iVar<nVar; iVar++)
          for (jVar=0; jVar<nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        /*--- Compute numerical flux Jacobian at node i ---*/
        
        for (iVar=0; iVar<nVar; iVar++) {
          for (jVar=0; jVar<nVar; jVar++) {
            for (kVar=0; kVar<nVar; kVar++) {
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
            }
          }
          
        }
        
        for (iVar = 0; iVar < nVar; iVar++) {
          delete [] Jacobian_b[iVar];
          delete [] DubDu[iVar];
        }
        delete [] Jacobian_b;
        delete [] DubDu;
      }
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Primitive variables, using the derived quantities ---*/
        V_boundary[0] = Temperature_b;
        for (iDim = 0; iDim < nDim; iDim++)
          V_boundary[iDim+1] = Velocity_b[iDim];
        V_boundary[nDim+1] = Pressure_b;
        V_boundary[nDim+2] = Density_b;
        V_boundary[nDim+3] = Enthalpy_b;
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_boundary[nDim+5] = FluidModel->GetLaminarViscosity();
        V_boundary[nDim+6] = node[iPoint]->GetEddyViscosity();
        V_boundary[nDim+7] = FluidModel->GetThermalConductivity();
        V_boundary[nDim+8] = FluidModel->GetCp();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_boundary);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Secondary variables ---*/
        S_domain = node[iPoint]->GetSecondary();
        
        /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
        
        S_boundary[0]= FluidModel->GetdPdrho_e();
        S_boundary[1]= FluidModel->GetdPde_rho();
        
        S_boundary[2]= FluidModel->GetdTdrho_e();
        S_boundary[3]= FluidModel->GetdTde_rho();
        
        /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
        
        S_boundary[4]= FluidModel->Getdmudrho_T();
        S_boundary[5]= FluidModel->GetdmudT_rho();
        
        S_boundary[6]= FluidModel->Getdktdrho_T();
        S_boundary[7]= FluidModel->GetdktdT_rho();
        
        visc_numerics->SetSecondary(S_domain, S_boundary);
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
  delete [] Velocity_b;
  delete [] Velocity_i;
  
  delete [] S_boundary;
  delete [] Lambda_i;
  delete [] u_b;
  delete [] dw;
  
  
  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
  
  delete [] delta_c;
  delete [] deltaprim;
  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] R_Matrix[iVar];
  }
  delete [] R_Matrix;
  
  
}

void CEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
  Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
  alpha, aa, bb, cc, dd, Area, UnitNormal[3];
  su2double *V_inlet, *V_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  su2double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
  su2double Gas_Constant       = config->GetGas_ConstantND();
  unsigned short Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the inlet ---*/
    
    V_inlet = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Build the fictitious intlet state based on characteristics ---*/
      

      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/

      switch (Kind_Inlet) {

        /*--- Total properties have been specified at the inlet. ---*/

        case TOTAL_CONDITIONS:

          /*--- Retrieve the specified total conditions for this inlet. ---*/

          if (gravity) P_Total = config->GetInlet_Ptotal(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;
          else P_Total  = config->GetInlet_Ptotal(Marker_Tag);
          T_Total  = config->GetInlet_Ttotal(Marker_Tag);
          Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

          /*--- Non-dim. the inputs if necessary. ---*/

          P_Total /= config->GetPressure_Ref();
          T_Total /= config->GetTemperature_Ref();

          /*--- Store primitives and set some variables for clarity. ---*/

          Density = V_domain[nDim+2];
          Velocity2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = V_domain[iDim+1];
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
          Pressure    = V_domain[nDim+1];
          H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
          SoundSpeed2 = Gamma*Pressure/Density;

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Total speed of sound ---*/

          SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

          /*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/

          alpha = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            alpha += UnitNormal[iDim]*Flow_Dir[iDim];

          /*--- Coefficients in the quadratic equation for the velocity ---*/

          aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
          bb = -1.0*Gamma_Minus_One*alpha*Riemann;
          cc =  0.5*Gamma_Minus_One*Riemann*Riemann
              -2.0*SoundSpeed_Total2/Gamma_Minus_One;

          /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/

          dd = bb*bb - 4.0*aa*cc;
          dd = sqrt(max(0.0, dd));
          Vel_Mag   = (-bb + dd)/(2.0*aa);
          Vel_Mag   = max(0.0, Vel_Mag);
          Velocity2 = Vel_Mag*Vel_Mag;

          /*--- Compute speed of sound from total speed of sound eqn. ---*/

          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/

          Mach2 = Velocity2/SoundSpeed2;
          Mach2 = min(1.0, Mach2);
          Velocity2   = Mach2*SoundSpeed2;
          Vel_Mag     = sqrt(Velocity2);
          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Compute new velocity vector at the inlet ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

          /*--- Static temperature from the speed of sound relation ---*/

          Temperature = SoundSpeed2/(Gamma*Gas_Constant);

          /*--- Static pressure using isentropic relation at a point ---*/

          Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);

          /*--- Density at the inlet from the gas law ---*/

          Density = Pressure/(Gas_Constant*Temperature);

          /*--- Using pressure, density, & velocity, compute the energy ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
          if (tkeNeeded) Energy += GetTke_Inf();

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Temperature;
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Velocity[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;

          /*--- Mass flow has been specified at the inlet. ---*/

        case MASS_FLOW:

          /*--- Retrieve the specified mass flow for the inlet. ---*/

          Density  = config->GetInlet_Ttotal(Marker_Tag);
          Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
          Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

          /*--- Non-dim. the inputs if necessary. ---*/

          Density /= config->GetDensity_Ref();
          Vel_Mag /= config->GetVelocity_Ref();

          /*--- Get primitives from current inlet state. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
          Pressure    = node[iPoint]->GetPressure();
          SoundSpeed2 = Gamma*Pressure/V_domain[nDim+2];

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Speed of sound squared for fictitious inlet state ---*/

          SoundSpeed2 = Riemann;
          for (iDim = 0; iDim < nDim; iDim++)
            SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitNormal[iDim];

          SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
          SoundSpeed2 = SoundSpeed2*SoundSpeed2;

          /*--- Pressure for the fictitious inlet state ---*/

          Pressure = SoundSpeed2*Density/Gamma;

          /*--- Energy for the fictitious inlet state ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Vel_Mag*Vel_Mag;
          if (tkeNeeded) Energy += GetTke_Inf();

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Pressure / ( Gas_Constant * Density);
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;
      }
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Pressure, P_Exit, Velocity[3],
  Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
  Area, UnitNormal[3];
  su2double *V_outlet, *V_domain;
  
  bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  su2double Gas_Constant     = config->GetGas_ConstantND();
  bool grid_movement      = config->GetGrid_Movement();
  string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    V_outlet = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Build the fictitious intlet state based on characteristics ---*/

      /*--- Retrieve the specified back pressure for this outlet. ---*/
      if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;
      else P_Exit = config->GetOutlet_Pressure(Marker_Tag);

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit = P_Exit/config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
         of boundary update depends on this. ---*/
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Mach_Exit  = sqrt(Velocity2)/SoundSpeed;

      if (Mach_Exit >= 1.0) {

        /*--- Supersonic exit flow: there are no incoming characteristics,
           so no boundary condition is necessary. Set outlet state to current
           state so that upwinding handles the direction of propagation. ---*/
        for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = V_domain[iVar];

      } else {

        /*--- Subsonic exit flow: there is one incoming characteristic,
           therefore one variable can be specified (back pressure) and is used
           to update the conservative variables. Compute the entropy and the
           acoustic Riemann variable. These invariants, as well as the
           tangential velocity components, are extrapolated. Adapted from an
           original implementation in the Stanford University multi-block
           (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
           der Weide, last modified 09-10-2007. ---*/

        Entropy = Pressure*pow(1.0/Density, Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Compute the new fictious state at the outlet ---*/
        Density    = pow(P_Exit/Entropy,1.0/Gamma);
        Pressure   = P_Exit;
        SoundSpeed = sqrt(Gamma*P_Exit/Density);
        Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();

        /*--- Conservative variables, using the derived quantities ---*/
        V_outlet[0] = Pressure / ( Gas_Constant * Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_outlet[iDim+1] = Velocity[iDim];
        V_outlet[nDim+1] = Pressure;
        V_outlet[nDim+2] = Density;
        V_outlet[nDim+3] = Energy + Pressure/Density;

      }
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                       CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_inlet, *V_domain;
  
  su2double Density, Pressure, Temperature, Energy, *Velocity, Velocity2;
  su2double Gas_Constant = config->GetGas_ConstantND();
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool viscous              = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
  
  /*--- Supersonic inlet flow: there are no outgoing characteristics,
   so all flow variables can be imposed at the inlet.
   First, retrieve the specified values for the primitive variables. ---*/
  
  Temperature = config->GetInlet_Temperature(Marker_Tag);
  Pressure    = config->GetInlet_Pressure(Marker_Tag);
  Velocity    = config->GetInlet_Velocity(Marker_Tag);
  
  /*--- Density at the inlet from the gas law ---*/
  
  Density = Pressure/(Gas_Constant*Temperature);
  
  /*--- Non-dim. the inputs if necessary. ---*/
  
  Temperature = Temperature/config->GetTemperature_Ref();
  Pressure    = Pressure/config->GetPressure_Ref();
  Density     = Density/config->GetDensity_Ref();
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = Velocity[iDim]/config->GetVelocity_Ref();
  
  /*--- Compute the energy from the specified state ---*/
  
  Velocity2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity2 += Velocity[iDim]*Velocity[iDim];
  Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
  if (tkeNeeded) Energy += GetTke_Inf();
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    
    V_inlet = GetCharacPrimVar(val_marker, iVertex);
    
    /*--- Primitive variables, using the derived quantities ---*/
    
    V_inlet[0] = Temperature;
    for (iDim = 0; iDim < nDim; iDim++)
      V_inlet[iDim+1] = Velocity[iDim];
    V_inlet[nDim+1] = Pressure;
    V_inlet[nDim+2] = Density;
    V_inlet[nDim+3] = Energy + Pressure/Density;
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                                        CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_outlet, *V_domain;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool viscous              = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Supersonic outlet flow: there are no ingoing characteristics,
   so all flow variables can should be interpolated from the domain. ---*/
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Primitive variables, using the derived quantities ---*/
      
      V_outlet[0] = V_domain[0];
      for (iDim = 0; iDim < nDim; iDim++)
        V_outlet[iDim+1] = V_domain[iDim+1];
      V_outlet[nDim+1] = V_domain[nDim+1];
      V_outlet[nDim+2] = V_domain[nDim+2];
      V_outlet[nDim+3] = V_domain[nDim+3];
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Pressure, Inflow_Pressure = 0.0, Velocity[3], Velocity2, Entropy, Target_Inflow_MassFlow = 0.0, Target_Inflow_Mach = 0.0, Density, Energy,
  Riemann, Area, UnitNormal[3], Vn, SoundSpeed, Vn_Exit, Inflow_Pressure_inc, Inflow_Pressure_old, Inflow_Mach_old, Inflow_MassFlow_old;
  su2double *V_inflow, *V_domain;
  
  su2double DampingFactor = config->GetDamp_Engine_Inflow();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  unsigned short Kind_Engine_Inflow = config->GetKind_Engine_Inflow();
  bool viscous              = config->GetViscous();
  su2double Gas_Constant = config->GetGas_ConstantND();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
  bool Engine_HalfModel = config->GetEngine_HalfModel();

  su2double *Normal = new su2double[nDim];
  
  
  if (Kind_Engine_Inflow == FAN_FACE_MACH) {
    
    /*--- Retrieve the specified target fan face mach at the nacelle. ---*/
    
    Target_Inflow_Mach = config->GetEngineInflow_Target(Marker_Tag);
    
    /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
    
    Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
    Inflow_Mach_old = config->GetInflow_Mach(Marker_Tag);
    
    /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/
    
    Inflow_Pressure_inc = - (1.0 - (Inflow_Mach_old/Target_Inflow_Mach)) * Baseline_Press;
    
    /*--- Estimate the new fan face pressure ---*/
    
    Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);
    
  }
  
  if (Kind_Engine_Inflow == FAN_FACE_MDOT) {
    
    /*--- Retrieve the specified target mass flow (non-dimensional) at the nacelle. ---*/
    
    Target_Inflow_MassFlow = config->GetEngineInflow_Target(Marker_Tag) / (config->GetDensity_Ref() * config->GetVelocity_Ref());
    
    if (config->GetSystemMeasurements() == US) Target_Inflow_MassFlow /= 32.174;
    
    if (Engine_HalfModel) Target_Inflow_MassFlow /= 2.0;

    /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
    
    Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
    Inflow_MassFlow_old = config->GetInflow_MassFlow(Marker_Tag);  // same here... it is a non dimensional value
    
    /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/
    
    Inflow_Pressure_inc = - (1.0 - (Inflow_MassFlow_old/Target_Inflow_MassFlow)) * Baseline_Press;
    
    /*--- Estimate the new fan face pressure ---*/
    
    Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);
    
  }
  
  /*--- No iterative scheme if we provide the static pressure ---*/
  
  if (Kind_Engine_Inflow == FAN_FACE_PRESSURE) {
    
    /*--- Retrieve the specified pressure (non-dimensional) at the nacelle. ---*/
    
    Inflow_Pressure = config->GetEngineInflow_Target(Marker_Tag) / config->GetPressure_Ref();
    
  }
  
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    
    V_inflow = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Subsonic nacelle inflow: there is one incoming characteristic,
       therefore one variable can be specified (back pressure) and is used
       to update the conservative variables.
       
       Compute the entropy and the acoustic variable. These
       riemann invariants, as well as the tangential velocity components,
       are extrapolated. ---*/
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Entropy = Pressure*pow(1.0/Density, Gamma);
      Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
      
      /*--- Compute the new fictious state at the outlet ---*/
      
      Density    = pow(Inflow_Pressure/Entropy,1.0/Gamma);
      Pressure   = Inflow_Pressure;
      SoundSpeed = sqrt(Gamma*Inflow_Pressure/Density);
      Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
      Velocity2  = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      
      Energy = Inflow_Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
      if (tkeNeeded) Energy += GetTke_Inf();
      
      /*--- Conservative variables, using the derived quantities ---*/
      
      V_inflow[0] = Pressure / ( Gas_Constant * Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_inflow[iDim+1] = Velocity[iDim];
      V_inflow[nDim+1] = Pressure;
      V_inflow[nDim+2] = Density;
      V_inflow[nDim+3] = Energy + Pressure/Density;
      V_inflow[nDim+4] = SoundSpeed;
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_inflow);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_inflow[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_inflow[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_inflow);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  delete [] Normal;
  
}


void CEulerSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Exhaust_Pressure, Exhaust_Temperature, Velocity[3], Velocity2, H_Exhaust, Temperature, Riemann, Area, UnitNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Exhaust2, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3];
  su2double *V_exhaust, *V_domain, Target_Exhaust_Pressure, Exhaust_Pressure_old, Exhaust_Pressure_inc;
  
  su2double Gas_Constant = config->GetGas_ConstantND();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double DampingFactor = config->GetDamp_Engine_Exhaust();
  su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Retrieve the specified exhaust pressure in the engine (non-dimensional). ---*/
  
  Target_Exhaust_Pressure = config->GetExhaust_Pressure_Target(Marker_Tag) / config->GetPressure_Ref();
  
  /*--- Retrieve the old exhaust pressure in the engine exhaust (this has been computed in a preprocessing). ---*/
  
  Exhaust_Pressure_old = config->GetExhaust_Pressure(Marker_Tag);
  
  /*--- Compute the Pressure increment ---*/
  
  Exhaust_Pressure_inc = (1.0 - (Exhaust_Pressure_old/Target_Exhaust_Pressure)) * Baseline_Press;
  
  /*--- Estimate the new exhaust pressure ---*/
  
  Exhaust_Pressure = (1.0 - DampingFactor) * Exhaust_Pressure_old + DampingFactor * (Exhaust_Pressure_old + Exhaust_Pressure_inc);
  
  /*--- The temperature is given (no iteration is required) ---*/
  
  Exhaust_Temperature  = config->GetExhaust_Temperature_Target(Marker_Tag);
  Exhaust_Temperature /= config->GetTemperature_Ref();
  
  /*--- The pressure is given (no iteration is required) ---*/
  
  Exhaust_Pressure  = config->GetExhaust_Pressure_Target(Marker_Tag);
  Exhaust_Pressure /= config->GetPressure_Ref();
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the exhaust ---*/
    
    V_exhaust = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
       therefore we can specify all but one state variable at the inlet.
       The outgoing Riemann invariant provides the final piece of info. ---*/
      
      /*--- Store primitives and set some variables for clarity. ---*/
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
      Pressure    = V_domain[nDim+1];
      H_Exhaust   = (Gamma*Gas_Constant/Gamma_Minus_One)*Exhaust_Temperature;
      SoundSpeed2 = Gamma*Pressure/Density;
      
      /*--- Compute the acoustic Riemann invariant that is extrapolated
       from the domain interior. ---*/
      
      Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
      for (iDim = 0; iDim < nDim; iDim++)
        Riemann += Velocity[iDim]*UnitNormal[iDim];
      
      /*--- Total speed of sound ---*/
      
      SoundSpeed_Exhaust2 = Gamma_Minus_One*(H_Exhaust - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
      
      /*--- The flow direction is defined by the surface normal ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        Flow_Dir[iDim] = -UnitNormal[iDim];
      
      /*--- Dot product of normal and flow direction. This should
       be negative due to outward facing boundary normal convention. ---*/
      
      alpha = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        alpha += UnitNormal[iDim]*Flow_Dir[iDim];
      
      /*--- Coefficients in the quadratic equation for the velocity ---*/
      
      aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
      bb = -1.0*Gamma_Minus_One*alpha*Riemann;
      cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Exhaust2/Gamma_Minus_One;
      
      /*--- Solve quadratic equation for velocity magnitude. Value must
       be positive, so the choice of root is clear. ---*/
      
      dd      = bb*bb - 4.0*aa*cc;
      dd      = sqrt(max(0.0, dd));
      Vel_Mag = (-bb + dd)/(2.0*aa);
      
      if (Vel_Mag >= 0.0) {
        
        Velocity2 = Vel_Mag*Vel_Mag;
        
        /*--- Compute speed of sound from total speed of sound eqn. ---*/
        
        SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
        Mach2       = Velocity2/SoundSpeed2;
        Velocity2   = Mach2*SoundSpeed2;
        Vel_Mag     = sqrt(Velocity2);
        SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
        
        /*--- Compute new velocity vector at the inlet ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
        
        /*--- Static temperature from the speed of sound relation ---*/
        
        Temperature = SoundSpeed2/(Gamma*Gas_Constant);
        
        /*--- Static pressure using isentropic relation at a point ---*/
        
        Pressure = Exhaust_Pressure*pow((Temperature/Exhaust_Temperature), Gamma/Gamma_Minus_One);
        
        /*--- Density at the exhaust from the gas law ---*/
        
        Density = Pressure/(Gas_Constant*Temperature);
        
        /*--- Using pressure, density, & velocity, compute the energy ---*/
        
        Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();
        
        /*--- Primitive variables, using the derived quantities ---*/
        
        V_exhaust[0] = Temperature;
        for (iDim = 0; iDim < nDim; iDim++)
          V_exhaust[iDim+1] = Velocity[iDim];
        V_exhaust[nDim+1] = Pressure;
        V_exhaust[nDim+2] = Density;
        V_exhaust[nDim+3] = Energy + Pressure/Density;
        V_exhaust[nDim+4] = sqrt(SoundSpeed2);
        
      }
      
      /*--- The flow goes in the wrong direction ---*/
      
      else {
        
        V_exhaust[0] = V_domain[0];
        for (iDim = 0; iDim < nDim; iDim++)
          V_exhaust[iDim+1] = V_domain[iDim+1];
        V_exhaust[nDim+1] = V_domain[nDim+1];
        V_exhaust[nDim+2] = V_domain[nDim+2];
        V_exhaust[nDim+3] = V_domain[nDim+3];
        V_exhaust[nDim+4] = V_domain[nDim+4];
        
      }
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_exhaust);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_exhaust[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_exhaust[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_exhaust);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler residual ---*/
  
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
  
}

void CEulerSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iVar, iMarker;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  su2double P_static, rho_static;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          for (iVar = 0; iVar < nPrimVar; iVar++) {
            PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
            PrimVar_j[iVar] = GetSlidingState(iMarker, iVertex, iVar);
          }

          /*--- Set primitive variables ---*/

          numerics->SetPrimitive( PrimVar_i, PrimVar_j );
          
          if( !( config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS ) ) {
          Secondary_i = node[iPoint]->GetSecondary();

            P_static   = PrimVar_j[nDim+1];
            rho_static = PrimVar_j[nDim+2];           
            FluidModel->SetTDState_Prho(P_static, rho_static);

            Secondary_j[0] = FluidModel->GetdPdrho_e();
            Secondary_j[1] = FluidModel->GetdPde_rho();  

            numerics->SetSecondary(Secondary_i, Secondary_j);
          }

          /*--- Set the normal vector ---*/

          geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
          for (iDim = 0; iDim < nDim; iDim++) 
            Normal[iDim] = -Normal[iDim];

          numerics->SetNormal(Normal);

          if (grid_movement)
            numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

          /*--- Compute the convective residual using an upwind scheme ---*/

          numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

          /*--- Add Residuals and Jacobians ---*/

          LinSysRes.AddBlock(iPoint, Residual);
          if (implicit) 
            Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
}

void CEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, GlobalIndex_iPoint, GlobalIndex_jPoint;
  unsigned short iDim, iVar;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
    GlobalIndex_jPoint = GetDonorGlobalIndex(val_marker, iVertex);
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
      
      /*--- Store the solution for both points ---*/
      
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
        PrimVar_j[iVar] = GetDonorPrimVar(val_marker, iVertex, iVar);
      }
      
      /*--- Set Conservative Variables ---*/
      
      numerics->SetPrimitive(PrimVar_i, PrimVar_j);
      
      /*--- Set Normal ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      numerics->SetNormal(Normal);
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      
      numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add Residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
    
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, GlobalIndex_iPoint, GlobalIndex_jPoint;
  unsigned short iDim, iVar;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
    GlobalIndex_jPoint = GetDonorGlobalIndex(val_marker, iVertex);
    
    if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
      
      /*--- Store the solution for both points ---*/
      
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
        PrimVar_j[iVar] = GetDonorPrimVar(val_marker, iVertex, iVar);
      }
      
      /*--- Set Conservative Variables ---*/
      
      numerics->SetPrimitive(PrimVar_i, PrimVar_j);
      
      /*--- Set Normal ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      numerics->SetNormal(Normal);
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      
      numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add Residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
    
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CEulerSolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
  
  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker, true);
  
}

void CEulerSolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                     CConfig *config, unsigned short val_marker) {
  
  BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker, false);
  
}

void CEulerSolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
        CConfig *config, unsigned short val_marker, bool inlet_surface) {

    unsigned short iDim, iVar, jVar, jDim;
    unsigned long iVertex, iPoint, iPoint_Normal, GlobalIndex_donor, GlobalIndex;
    su2double Pressure, Velocity[3], Target_Press_Jump, Target_Temp_Jump,
    Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Vn_Inlet, Mach_Outlet,
    Area, UnitNormal[3], *V_outlet, *V_domain, *V_inlet, P_Total, T_Total, H_Total, Temperature,
    Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag, alpha, aa, bb, cc, dd;
    su2double Factor, P_static, T_static, SoS_outlet, Rho_outlet, Rho_inlet;
    su2double Vel_normal_inlet[3], Vel_tangent_inlet[3], Vel_inlet[3];
    su2double Vel_normal_outlet[3], Vel_tangent_outlet[3], Vel_outlet[3];
    su2double Vel_normal_inlet_, Vel_tangent_inlet_, Vel_inlet_;
    su2double Vel_normal_outlet_, Vel_outlet_;
    su2double turb_ke, a2, phi;
    bool ReverseFlow;

    bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    su2double Gas_Constant  = config->GetGas_ConstantND();
    bool viscous            = config->GetViscous();
    bool grid_movement      = config->GetGrid_Movement();
    bool tkeNeeded                = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
            (config->GetKind_Turb_Model() == SST));
    bool ratio                = (config->GetActDisk_Jump() == RATIO);
    su2double SecondaryFlow = config->GetSecondaryFlow_ActDisk();

    su2double *Normal = new su2double[nDim];
    su2double *Flow_Dir = new su2double[nDim];

    /*--- Loop over all the vertices on this boundary marker ---*/

    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

        if (inlet_surface) {
            V_inlet = GetCharacPrimVar(val_marker, iVertex);
            V_outlet = GetDonorPrimVar(val_marker, iVertex);
        }
        else {
            V_outlet = GetCharacPrimVar(val_marker, iVertex);
            V_inlet = GetDonorPrimVar(val_marker, iVertex);
        }

        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        iPoint_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
        GlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
        GlobalIndex_donor = GetDonorGlobalIndex(val_marker, iVertex);

        /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

//      if ((geometry->node[iPoint]->GetDomain()) &&
//              (GlobalIndex != GlobalIndex_donor) && (!ActDisk_Perimeter)) {

            if ((geometry->node[iPoint]->GetDomain()) &&
                    (GlobalIndex != GlobalIndex_donor)) {

            /*--- Normal vector for this vertex (negative for outward convention) ---*/

            geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            conv_numerics->SetNormal(Normal);

            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);

            for (iDim = 0; iDim < nDim; iDim++)
                UnitNormal[iDim] = Normal[iDim]/Area;

            /*--- Current solution at this boundary node and jumps values ---*/

            V_domain = node[iPoint]->GetPrimitive();
            Target_Press_Jump = GetActDisk_DeltaP(val_marker, iVertex);
            Target_Temp_Jump = GetActDisk_DeltaT(val_marker, iVertex);

            if (inlet_surface) {
                if (ratio) { P_static = V_outlet[nDim+1]/Target_Press_Jump; T_static = V_outlet[0]/Target_Temp_Jump; }
                else { P_static = V_outlet[nDim+1] - Target_Press_Jump; T_static = V_outlet[0] - Target_Temp_Jump; }
            }
            else {
                if (ratio) { P_static = V_inlet[nDim+1]*Target_Press_Jump; T_static = V_inlet[0]*Target_Temp_Jump; }
                else { P_static = V_inlet[nDim+1] + Target_Press_Jump; T_static = V_inlet[0] + Target_Temp_Jump; }
            }

            /*--- Check the flow direction. Project the flow into the normal to the inlet face ---*/

            Vn = 0.0; ReverseFlow = false;
            for (iDim = 0; iDim < nDim; iDim++) {  Vn += V_domain[iDim+1]*UnitNormal[iDim]; }
            if ((inlet_surface) && (Vn < 0.0)) { ReverseFlow = true; }
            if ((!inlet_surface) && (Vn > 0.0)) { ReverseFlow = true; }

            /*---Depending of inflow outflow and flow reversal it takes a decision about the BC ---*/

            if (!ReverseFlow) {

                /*--- Subsonic inlet ---*/

                if (inlet_surface) {

                    /*--- Build the fictitious intlet state based on characteristics.
         Retrieve the specified back pressure for this inlet ---*/

                    Density = V_domain[nDim+2];
                    Velocity2 = 0.0; Vn = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Velocity[iDim] = V_domain[iDim+1];
                        Velocity2 += Velocity[iDim]*Velocity[iDim];
                        Vn += Velocity[iDim]*UnitNormal[iDim];
                    }
                    Pressure   = V_domain[nDim+1];
                    SoundSpeed = sqrt(Gamma*Pressure/Density);

                    Entropy = Pressure*pow(1.0/Density, Gamma);
                    Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

                    /*--- Compute the new fictious state at the outlet ---*/

                    Pressure   = P_static;
                    Density    = pow(Pressure/Entropy,1.0/Gamma);
                    SoundSpeed = sqrt(Gamma*Pressure/Density);
                    Vn_Inlet    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;

                    Velocity2  = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Velocity[iDim] = Velocity[iDim] + (Vn_Inlet-Vn)*UnitNormal[iDim];
                        Velocity2 += Velocity[iDim]*Velocity[iDim];
                    }
                    Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
                    if (tkeNeeded) Energy += GetTke_Inf();

                    /*--- Conservative variables, using the derived quantities ---*/

                    V_inlet[0] = Pressure / ( Gas_Constant * Density);
                    for (iDim = 0; iDim < nDim; iDim++)
                        V_inlet[iDim+1] = Velocity[iDim];
                    V_inlet[nDim+1] = Pressure;
                    V_inlet[nDim+2] = Density;
                    V_inlet[nDim+3] = Energy + Pressure/Density;
                    V_inlet[nDim+4] = SoundSpeed;
                    conv_numerics->SetPrimitive(V_domain, V_inlet);

                }

                /*--- Subsonic outlet ---*/

                else {

                    FluidModel->SetTDState_PT(P_static, T_static);
                    SoS_outlet = FluidModel->GetSoundSpeed();
                    Rho_outlet = FluidModel->GetDensity();

                    /*--- We use the velocity and the density from the flow inlet
         to evaluate flow direction and mass flow ---*/

                    Rho_inlet = V_inlet[nDim+2];
                    for (iDim = 0; iDim < nDim; iDim++)
                        Vel_inlet[iDim] = V_inlet[iDim+1];

                    Vel_normal_inlet_ = 0.0; Vel_inlet_ = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Vel_normal_inlet[iDim] = -Vel_inlet[iDim]*UnitNormal[iDim];
                        Vel_normal_inlet_ += Vel_normal_inlet[iDim]*Vel_normal_inlet[iDim];
                        Vel_inlet_+= Vel_inlet[iDim]*Vel_inlet[iDim];
                    }
                    Vel_inlet_ = sqrt(Vel_inlet_);
                    Vel_normal_inlet_ = sqrt(Vel_normal_inlet_);

                    Vel_tangent_inlet_ = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Vel_tangent_inlet[iDim] = Vel_inlet[iDim] - Vel_normal_inlet[iDim];
                        Vel_tangent_inlet_ += Vel_tangent_inlet[iDim]*Vel_tangent_inlet[iDim];
                    }
                    Vel_tangent_inlet_ = sqrt(Vel_tangent_inlet_);

                    /*--- Mass flow conservation (normal direction) and
           no jump in the tangential velocity ---*/

                    Vel_normal_outlet_ = (1.0-SecondaryFlow/100.0)*(Rho_inlet*Vel_normal_inlet_)/Rho_outlet;

                    Vel_outlet_ = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Vel_normal_outlet[iDim] = -Vel_normal_outlet_*UnitNormal[iDim];
                        Vel_tangent_outlet[iDim] = Vel_tangent_inlet[iDim];
                        Vel_outlet[iDim] = Vel_normal_outlet[iDim] + Vel_tangent_outlet[iDim];
                        Vel_outlet_ += Vel_outlet[iDim]*Vel_outlet[iDim];
                    }
                    Vel_outlet_ = sqrt(Vel_outlet_);

                    Mach_Outlet = min(Vel_outlet_/SoS_outlet, 1.0);

                    /*--- Reevaluate the Total Pressure and Total Temperature using the
             Fan Face Mach number and the static values from the jum condition ---*/

                    Factor = 1.0 + 0.5*Mach_Outlet*Mach_Outlet*Gamma_Minus_One;
                    P_Total = P_static * pow(Factor, Gamma/Gamma_Minus_One);
                    T_Total = T_static * Factor;

                    /*--- Flow direction using the velocity direction at the outlet  ---*/

                    if (Vel_outlet_ != 0.0) {
                        for (iDim = 0; iDim < nDim; iDim++) Flow_Dir[iDim] = Vel_outlet[iDim]/Vel_outlet_;
                    }
                    else {
                        for (iDim = 0; iDim < nDim; iDim++) Flow_Dir[iDim] = 0.0;
                    }

                    /*--- Store primitives and set some variables for clarity. ---*/

                    Density = V_domain[nDim+2];
                    Velocity2 = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Velocity[iDim] = V_domain[iDim+1];
                        Velocity2 += Velocity[iDim]*Velocity[iDim];
                    }
                    Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
                    Pressure    = V_domain[nDim+1];
                    H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
                    SoundSpeed2 = Gamma*Pressure/Density;

                    /*--- Compute the acoustic Riemann invariant that is extrapolated
           from the domain interior. ---*/

                    Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
                    for (iDim = 0; iDim < nDim; iDim++)
                        Riemann += Velocity[iDim]*UnitNormal[iDim];

                    /*--- Total speed of sound ---*/

                    SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

                    /*--- Dot product of normal and flow direction. This should
           be negative due to outward facing boundary normal convention. ---*/

                    alpha = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++)
                        alpha += UnitNormal[iDim]*Flow_Dir[iDim];

                    /*--- Coefficients in the quadratic equation for the velocity ---*/

                    aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
                    bb = -1.0*Gamma_Minus_One*alpha*Riemann;
                    cc =  0.5*Gamma_Minus_One*Riemann*Riemann
                            -2.0*SoundSpeed_Total2/Gamma_Minus_One;

                    /*--- Solve quadratic equation for velocity magnitude. Value must
           be positive, so the choice of root is clear. ---*/

                    dd = bb*bb - 4.0*aa*cc;
                    dd = sqrt(max(0.0, dd));
                    Vel_Mag   = (-bb + dd)/(2.0*aa);
                    Vel_Mag   = max(0.0, Vel_Mag);
                    Velocity2 = Vel_Mag*Vel_Mag;

                    /*--- Compute speed of sound from total speed of sound eqn. ---*/

                    SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

                    /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/

                    Mach2 = min(1.0, Velocity2/SoundSpeed2);
                    Velocity2   = Mach2*SoundSpeed2;
                    Vel_Mag     = sqrt(Velocity2);
                    SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

                    /*--- Compute new velocity vector at the exit ---*/

                    for (iDim = 0; iDim < nDim; iDim++)
                        Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

                    /*--- Static temperature from the speed of sound relation ---*/

                    Temperature = SoundSpeed2/(Gamma*Gas_Constant);

                    /*--- Static pressure using isentropic relation at a point ---*/

                    Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);

                    /*--- Density at the inlet from the gas law ---*/

                    Density = Pressure/(Gas_Constant*Temperature);

                    /*--- Using pressure, density, & velocity, compute the energy ---*/

                    Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
                    if (tkeNeeded) Energy += GetTke_Inf();

                    /*--- Primitive variables, using the derived quantities ---*/

                    V_outlet[0] = Temperature;
                    for (iDim = 0; iDim < nDim; iDim++)
                        V_outlet[iDim+1] = Velocity[iDim];
                    V_outlet[nDim+1] = Pressure;
                    V_outlet[nDim+2] = Density;
                    V_outlet[nDim+3] = Energy + Pressure/Density;
                    V_outlet[nDim+4] = sqrt(SoundSpeed2);
                    conv_numerics->SetPrimitive(V_domain, V_outlet);

                }

                /*--- Grid Movement ---*/

                if (grid_movement)
                    conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

                /*--- Compute the residual using an upwind scheme ---*/

                conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

                /*--- Update residual value ---*/

                LinSysRes.AddBlock(iPoint, Residual);

                /*--- Jacobian contribution for implicit integration ---*/

                if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

                /*--- Viscous contribution ---*/

                if (viscous) {

                    /*--- Set laminar and eddy viscosity at the infinity ---*/

                    if (inlet_surface) {
                        V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
                        V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
                    }
                    else {
                        V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
                        V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
                    }

                    /*--- Set the normal vector and the coordinates ---*/

                    visc_numerics->SetNormal(Normal);
                    visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint_Normal]->GetCoord());

                    /*--- Primitive variables, and gradient ---*/

                    if (inlet_surface) visc_numerics->SetPrimitive(V_domain, V_inlet);
                    else visc_numerics->SetPrimitive(V_domain, V_outlet);

                    visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

                    /*--- Turbulent kinetic energy ---*/

                    if (config->GetKind_Turb_Model() == SST)
                        visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));

                    /*--- Compute and update residual ---*/

                    visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
                    LinSysRes.SubtractBlock(iPoint, Residual);

                    /*--- Jacobian contribution for implicit integration ---*/

                    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

                }

            }

            /*--- If reverse flow, then we apply Euler standart boundary
       condition (no viscous term) ---*/

            if (ReverseFlow) {

                /*--- Set the value of the charactestisc variables... important for turbulence model ---*/

                for (iVar = 0; iVar < nPrimVar; iVar++) {
                    if (inlet_surface) { V_inlet[iVar] = V_domain[iVar];    }
                    else { V_outlet[iVar] = V_domain[iVar]; }
                }

                /*--- Get the pressure ---*/

                Pressure = node[iPoint]->GetPressure();

                /*--- Add the kinetic energy correction ---*/

                if (tkeNeeded) {
                    turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
                    Pressure += (2.0/3.0)*node[iPoint]->GetDensity()*turb_ke;
                }

                /*--- Compute the residual ---*/

                Residual[0] = 0.0;
                for (iDim = 0; iDim < nDim; iDim++)
                    Residual[iDim+1] = Pressure*UnitNormal[iDim]*Area;

                /*--- Add value to the residual ---*/

                LinSysRes.AddBlock(iPoint, Residual);

                /*--- Form Jacobians for implicit computations ---*/

                if (implicit) {

                    for (iVar = 0; iVar < nVar; iVar++) {
                        for (jVar = 0; jVar < nVar; jVar++)
                            Jacobian_i[iVar][jVar] = 0.0;
                    }

                    a2 = Gamma-1.0;
                    phi = 0.5*a2*node[iPoint]->GetVelocity2();
                    for (iVar = 0; iVar < nVar; iVar++) {
                        Jacobian_i[0][iVar] = 0.0;
                        Jacobian_i[nDim+1][iVar] = 0.0;
                    }
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Jacobian_i[iDim+1][0] = -phi*Normal[iDim];
                        for (jDim = 0; jDim < nDim; jDim++)
                            Jacobian_i[iDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim)*Normal[iDim];
                        Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
                    }

                    Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);

                }

            }

        }

    }

    /*--- Free locally allocated memory ---*/

    delete [] Normal;
    delete [] Flow_Dir;

}

void CEulerSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container,
                                CConfig *config, unsigned short val_marker) { }

void CEulerSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

void CEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!grid_movement) {
    
    /*--- Loop over all nodes (excluding halos) ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/
      
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
    
  }
  
  else {
    
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Get indices for nodes i & j plus the face normal ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      Normal = geometry->edge[iEdge]->GetNormal();
      
      /*--- Grid velocities stored at nodes i & j ---*/
      
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      
      /*--- Compute the GCL term by averaging the grid velocities at the
       edge mid-point and dotting with the face normal. ---*/
      
      Residual_GCL = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_GCL += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      
      /*--- Compute the GCL component of the source term for node i ---*/
      
      U_time_n = node[iPoint]->GetSolution_time_n();
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      U_time_n = node[jPoint]->GetSolution_time_n();
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---   Loop over the boundary edges ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- Get the index for node i plus the boundary face normal ---*/
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Grid velocities stored at boundary node i ---*/
        
        GridVel_i = geometry->node[iPoint]->GetGridVel();
        
        /*--- Compute the GCL term by dotting the grid velocity with the face
         normal. The normal is negated to match the boundary convention. ---*/
        
        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];
        
        /*--- Compute the GCL component of the source term for node i ---*/
        
        U_time_n = node[iPoint]->GetSolution_time_n();
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
        LinSysRes.AddBlock(iPoint, Residual);
        
      }
    }
    
    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/
      
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}

void CEulerSolver::SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement,
                                        CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) {
    unsigned short iDim;
    unsigned long iVertex;
    su2double *Coord, VarCoord[3] = {0,0,0};

  #ifndef HAVE_MPI
    unsigned long iPoint_Donor, iPoint;
    unsigned short iMarker;
    su2double *CoordDonor, *DisplacementDonor;

    for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {

      if (flow_config->GetMarker_All_FSIinterface(iMarker) != 0) {

        for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {

          iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();

          iPoint_Donor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorPoint();

          Coord = flow_geometry[MESH_0]->node[iPoint]->GetCoord();

          CoordDonor = fea_geometry[MESH_0]->node[iPoint_Donor]->GetCoord();

          /*--- The displacements come from the predicted solution ---*/
          DisplacementDonor = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution_Pred();

          for (iDim = 0; iDim < nDim; iDim++)

            VarCoord[iDim] = (CoordDonor[iDim]+DisplacementDonor[iDim])-Coord[iDim];

          flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
        }
      }
    }

  #else

    int rank = MASTER_NODE;
    int size = SINGLE_NODE;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    unsigned long nLocalVertexStruct = 0, nLocalVertexFlow = 0;

    unsigned short nMarkerFSI, nMarkerStruct, nMarkerFlow;      // Number of markers on FSI problem, FEA and Flow side
    unsigned short iMarkerFSI, iMarkerStruct, iMarkerFlow;      // Variables for iteration over markers

    unsigned long MaxLocalVertexStruct = 0, MaxLocalVertexFlow = 0;

    unsigned long nBuffer_StructCoord = 0, nBuffer_FlowNewCoord = 0;
    unsigned long nBuffer_DonorIndices = 0, nBuffer_SetIndex = 0;

    unsigned long Point_Flow, Point_Struct;
    long Point_Flow_Rcv, Processor_Flow_Rcv;
    unsigned long Processor_Flow;

    int Marker_Flow = -1, Marker_Struct = -1;

    int iProcessor, nProcessor = 0;


    su2double *Coord_Struct, *Displacement_Struct;

    /*--- Number of markers on the FSI interface ---*/

    nMarkerFSI     = (flow_config->GetMarker_n_FSIinterface())/2;
    nMarkerStruct  = fea_geometry[MESH_0]->GetnMarker();
    nMarkerFlow    = flow_geometry[MESH_0]->GetnMarker();

    nProcessor = size;

    /*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
    /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/

    for (iMarkerFSI = 1; iMarkerFSI <= nMarkerFSI; iMarkerFSI++) {

        Marker_Struct = -1;
        Marker_Flow = -1;

        /*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
        unsigned long Buffer_Send_nVertexStruct[1], *Buffer_Recv_nVertexStruct = NULL;
        unsigned long Buffer_Send_nVertexFlow[1], *Buffer_Recv_nVertexFlow = NULL;

        /*--- The markers on the fluid and structural side are tagged with the same index.
         *--- This is independent of the MPI domain decomposition.
         *--- We need to loop over all markers on structural side and get the number of nodes
         *--- that belong to each FSI marker for each processor ---*/

        /*--- On the structural side ---*/

        for (iMarkerStruct = 0; iMarkerStruct < nMarkerStruct; iMarkerStruct++) {
            /*--- If the tag GetMarker_All_FSIinterface(iMarkerFEA) equals the index we are looping at ---*/
            if ( fea_config->GetMarker_All_FSIinterface(iMarkerStruct) == iMarkerFSI ) {
                /*--- We have identified the local index of the FEA marker ---*/
                /*--- Store the number of local points that belong to markFEA on each processor ---*/
                /*--- This includes the halo nodes ---*/
                nLocalVertexStruct = fea_geometry[MESH_0]->GetnVertex(iMarkerStruct);
                /*--- Store the identifier for the structural marker ---*/
                Marker_Struct = iMarkerStruct;
                /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
                break;
            }
            else {
                /*--- If the tag hasn't matched any tag within the FEA markers ---*/
                nLocalVertexStruct = 0;
                Marker_Struct = -1;
            }
        }

        /*--- On the fluid side ---*/

        for (iMarkerFlow = 0; iMarkerFlow < nMarkerFlow; iMarkerFlow++) {
            /*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
            if ( flow_config->GetMarker_All_FSIinterface(iMarkerFlow) == iMarkerFSI ) {
                /*--- We have identified the local index of the Flow marker ---*/
                /*--- Store the number of local points that belong to markFlow on each processor ---*/
                /*--- This includes the halo nodes ---*/
                nLocalVertexFlow = flow_geometry[MESH_0]->GetnVertex(iMarkerFlow);
                /*--- Store the identifier for the fluid marker ---*/
                Marker_Flow = iMarkerFlow;
                /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
                break;
            }
            else {
                /*--- If the tag hasn't matched any tag within the Flow markers ---*/
                nLocalVertexFlow = 0;
                Marker_Flow = -1;
            }
        }

        Buffer_Send_nVertexStruct[0] = nLocalVertexStruct;                               // Retrieve total number of vertices on FEA marker
        Buffer_Send_nVertexFlow[0] = nLocalVertexFlow;                               // Retrieve total number of vertices on Flow marker
        if (rank == MASTER_NODE) Buffer_Recv_nVertexStruct = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
        if (rank == MASTER_NODE) Buffer_Recv_nVertexFlow = new unsigned long[size];  // Allocate memory to receive how many vertices are on each rank on the fluid side

        /*--- We receive MaxLocalVertexFEA as the maximum number of vertices in one single processor on the structural side---*/
        SU2_MPI::Allreduce(&nLocalVertexStruct, &MaxLocalVertexStruct, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        /*--- We receive MaxLocalVertexFlow as the maximum number of vertices in one single processor on the fluid side ---*/
        SU2_MPI::Allreduce(&nLocalVertexFlow, &MaxLocalVertexFlow, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

        /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
        SU2_MPI::Gather(&Buffer_Send_nVertexStruct, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexStruct, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
        /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the fluid side ---*/
        SU2_MPI::Gather(&Buffer_Send_nVertexFlow, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexFlow, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

        /*--- We will be gathering the structural coordinates into the master node ---*/
        /*--- Then we will distribute them using a scatter operation into the appropriate fluid processor ---*/
        nBuffer_StructCoord = MaxLocalVertexStruct * nDim;
        nBuffer_FlowNewCoord = MaxLocalVertexFlow * nDim;

        /*--- We will be gathering donor index and donor processor (for structure -> donor = flow) ---*/
        /*--- Then we will pass on to the fluid side the index (flow point) to the appropriate processor ---*/
        nBuffer_DonorIndices = 2 * MaxLocalVertexStruct;
        nBuffer_SetIndex = MaxLocalVertexFlow;

        /*--- Send and Recv buffers ---*/

        /*--- Buffers to send and receive the structural coordinates ---*/
        su2double *Buffer_Send_StructCoord = new su2double[nBuffer_StructCoord];
        su2double *Buffer_Recv_StructCoord = NULL;

        /*--- Buffers to send and receive the donor index and processor ---*/
        long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
        long *Buffer_Recv_DonorIndices = NULL;

        /*--- Buffers to send and receive the new fluid coordinates ---*/
        su2double *Buffer_Send_FlowNewCoord = NULL;
        su2double *Buffer_Recv_FlowNewCoord = new su2double[nBuffer_FlowNewCoord];

        /*--- Buffers to send and receive the fluid index ---*/
        long *Buffer_Send_SetIndex = NULL;
        long *Buffer_Recv_SetIndex = new long[nBuffer_SetIndex];

        /*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/

        if (rank == MASTER_NODE) {
            Buffer_Recv_StructCoord  = new su2double[size*nBuffer_StructCoord];
            Buffer_Recv_DonorIndices = new long[size*nBuffer_DonorIndices];
            Buffer_Send_FlowNewCoord = new su2double[size*nBuffer_FlowNewCoord];
            Buffer_Send_SetIndex     = new long[size*nBuffer_SetIndex];
        }

        /*--- On the structural side ---*/

        /*--- If this processor owns the marker we are looping at on the structural side ---*/

        /*--- First we initialize all of the indices and processors to -1 ---*/
        /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
        for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
            Buffer_Send_DonorIndices[iVertex] = -1;

        if (Marker_Struct >= 0) {

            /*--- We have identified the local index of the FEA marker ---*/
            /*--- We loop over all the vertices in that marker and in that particular processor ---*/

            for (iVertex = 0; iVertex < nLocalVertexStruct; iVertex++) {

                Point_Struct = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetNode();

                Point_Flow = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetDonorPoint();

                Processor_Flow = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetDonorProcessor();

                Coord_Struct = fea_geometry[MESH_0]->node[Point_Struct]->GetCoord();

                /*--- The displacements come from the predicted solution ---*/
                Displacement_Struct = fea_solution[MESH_0][FEA_SOL]->node[Point_Struct]->GetSolution_Pred();

                for (iDim = 0; iDim < nDim; iDim++) {
                    Buffer_Send_StructCoord[iVertex*nDim+iDim] = Coord_Struct[iDim] + Displacement_Struct[iDim];
                }
                /*--- If this processor owns the node ---*/
                if (fea_geometry[MESH_0]->node[Point_Struct]->GetDomain()) {
                    Buffer_Send_DonorIndices[2*iVertex]     = Point_Flow;
                    Buffer_Send_DonorIndices[2*iVertex + 1] = Processor_Flow;
                }
                else {
                    /*--- We set the values to be -1 to be able to identify them later as halo nodes ---*/
                    Buffer_Send_DonorIndices[2*iVertex]     = -1;
                    Buffer_Send_DonorIndices[2*iVertex + 1] = -1;
                }

            }
        }

        /*--- Once all the messages have been sent, we gather them all into the MASTER_NODE ---*/
        SU2_MPI::Gather(Buffer_Send_StructCoord, nBuffer_StructCoord, MPI_DOUBLE, Buffer_Recv_StructCoord, nBuffer_StructCoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);

        /*--- Counter to determine where in the array we have to set the information ---*/
        long *Counter_Processor_Flow = NULL;
        long iProcessor_Struct = 0, iIndex_Struct = 0;
        long iProcessor_Flow = 0, iPoint_Flow = 0, iIndex_Flow = 0;

        /*--- Now we pack the information to send it over to the different processors ---*/

        if (rank == MASTER_NODE) {

            /*--- We set the counter to 0 ---*/
            Counter_Processor_Flow = new long[nProcessor];
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                Counter_Processor_Flow[iProcessor] = 0;
            }

            /*--- First we initialize the index vector to -1 ---*/
            /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
            for (iVertex = 0; iVertex < nProcessor*nBuffer_SetIndex; iVertex++)
                Buffer_Send_SetIndex[iVertex] = -2;

            /*--- As of now we do the loop over the structural points ---*/
            /*--- The number of points for flow and structure does not necessarily have to match ---*/
            /*--- In fact, it's possible that a processor asks for nFlow nodes and there are only ---*/
            /*--- nStruc < nFlow available; this is due to halo nodes ---*/

            /*--- For every processor from which we have received information ---*/
            /*--- (This is, for every processor on the structural side) ---*/
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {

                /*--- This is the initial index on the coordinates buffer for that particular processor on the structural side ---*/
                iProcessor_Struct = iProcessor*nBuffer_StructCoord;
                /*--- This is the initial index on the donor index/processor buffer for that particular processor on the structural side ---*/
                iIndex_Struct = iProcessor*nBuffer_DonorIndices;

                /*--- For every vertex in the information retreived from iProcessor ---*/
                for (iVertex = 0; iVertex < Buffer_Recv_nVertexStruct[iProcessor]; iVertex++) {

                    /*--- The processor and index for the flow are: ---*/
                    Processor_Flow_Rcv = Buffer_Recv_DonorIndices[iIndex_Struct+iVertex*2+1];
                    Point_Flow_Rcv     = Buffer_Recv_DonorIndices[iIndex_Struct+iVertex*2];

                    /*--- Load the buffer at the appropriate position ---*/
                    /*--- This is determined on the fluid side by:
                     *--- Processor_Flow*nBuffer_FlowNewCoord -> Initial position of the processor array (fluid side)
                     *--- +
                     *--- Counter_Processor_Flow*nDim -> Initial position of the nDim array for the particular point on the fluid side
                     *--- +
                     *--- iDim -> Position within the nDim array that corresponds to a point
                     *---
                     *--- While on the structural side is:
                     *--- iProcessor*nBuffer_StructCoord -> Initial position on the processor array (structural side)
                     *--- +
                     *--- iVertex*nDim -> Initial position of the nDim array for the particular point on the structural side
                     */

                    /*--- We check that we are not setting the value for a halo node ---*/
                    if (Point_Flow_Rcv != -1) {
                        iProcessor_Flow = Processor_Flow_Rcv*nBuffer_FlowNewCoord;
                        iIndex_Flow = Processor_Flow_Rcv*nBuffer_SetIndex;
                        iPoint_Flow = Counter_Processor_Flow[Processor_Flow_Rcv]*nDim;

                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_FlowNewCoord[iProcessor_Flow + iPoint_Flow + iDim] = Buffer_Recv_StructCoord[iProcessor_Struct + iVertex*nDim + iDim];

                        /*--- We set the fluid index at an appropriate position matching the coordinates ---*/
                        Buffer_Send_SetIndex[iIndex_Flow + Counter_Processor_Flow[Processor_Flow_Rcv]] = Point_Flow_Rcv;

                        Counter_Processor_Flow[Processor_Flow_Rcv]++;
                    }

                }

            }

        }

        /*--- Once all the messages have been prepared, we scatter them all from the MASTER_NODE ---*/
        SU2_MPI::Scatter(Buffer_Send_FlowNewCoord, nBuffer_FlowNewCoord, MPI_DOUBLE, Buffer_Recv_FlowNewCoord, nBuffer_FlowNewCoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        SU2_MPI::Scatter(Buffer_Send_SetIndex, nBuffer_SetIndex, MPI_LONG, Buffer_Recv_SetIndex, nBuffer_SetIndex, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);

        long indexPoint_iVertex, Point_Flow_Check;

        /*--- For the flow marker we are studying ---*/
        if (Marker_Flow >= 0) {

            /*--- We have identified the local index of the Flow marker ---*/
            /*--- We loop over all the vertices in that marker and in that particular processor ---*/

            for (iVertex = 0; iVertex < nLocalVertexFlow; iVertex++) {

                Point_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNode();

                if (flow_geometry[MESH_0]->node[Point_Flow]->GetDomain()) {
                    /*--- Find the index of the point Point_Flow in the buffer Buffer_Recv_SetIndex ---*/
                    indexPoint_iVertex = std::distance(Buffer_Recv_SetIndex, std::find(Buffer_Recv_SetIndex, Buffer_Recv_SetIndex + MaxLocalVertexFlow, Point_Flow));

                    Point_Flow_Check = Buffer_Recv_SetIndex[indexPoint_iVertex];

                    if (Point_Flow_Check < 0) {
                        cout << "WARNING: A nonphysical point is being considered for mesh deformation." << endl;
                        exit(EXIT_FAILURE);
                    }

                    Coord = flow_geometry[MESH_0]->node[Point_Flow]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                        VarCoord[iDim] = (Buffer_Recv_FlowNewCoord[indexPoint_iVertex*nDim+iDim])-Coord[iDim];

                    flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->SetVarCoord(VarCoord);

                }

            }

        }

        delete [] Buffer_Send_StructCoord;
        delete [] Buffer_Send_DonorIndices;
        delete [] Buffer_Recv_FlowNewCoord;
        delete [] Buffer_Recv_SetIndex;

        if (rank == MASTER_NODE) {
            delete [] Buffer_Recv_nVertexStruct;
            delete [] Buffer_Recv_nVertexFlow;
            delete [] Buffer_Recv_StructCoord;
            delete [] Buffer_Recv_DonorIndices;
            delete [] Buffer_Send_FlowNewCoord;
            delete [] Buffer_Send_SetIndex;
            delete [] Counter_Processor_Flow;
        }


    }

  #endif

    flow_grid_movement->SetVolume_Deformation(flow_geometry[MESH_0], flow_config, true);

}

void CEulerSolver::SetFlow_Displacement_Int(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement,
                                        CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) {
    unsigned short iMarker, iDim, iDonor, nDonor;
    unsigned long iVertex;
    su2double VarCoord[3];

    unsigned long iPoint_Donor;
    su2double *DisplacementDonor, *DisplacementDonor_Prev, coeff;

    for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {

      if (flow_config->GetMarker_All_FSIinterface(iMarker) != 0) {

        for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {

          for (iDim = 0; iDim < nDim; iDim++)
            VarCoord[iDim]=0.0;

          nDonor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetnDonorPoints();

          for (iDonor = 0; iDonor < nDonor; iDonor++) {
            iPoint_Donor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetInterpDonorPoint(iDonor);
            coeff = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorCoeff(iDonor);

            /*--- The displacements come from the predicted solution ---*/
            DisplacementDonor = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution_Pred();

            DisplacementDonor_Prev = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution_Pred_Old();

            for (iDim = 0; iDim < nDim; iDim++)

              VarCoord[iDim] += (DisplacementDonor[iDim] - DisplacementDonor_Prev[iDim])*coeff;
          }

          flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
        }
      }
    }
    flow_grid_movement->SetVolume_Deformation(flow_geometry[MESH_0], flow_config, true);

}

void CEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter) {
  
  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine, dull_val;
  bool grid_movement  = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool steady_restart = config->GetSteadyRestart();
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  string UnstExt, text_line;
  ifstream restart_file;
  
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry[iZone]->GetnZone();

  string restart_filename = config->GetSolution_FlowFileName();

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/
  
  if (dual_time || time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);
  
  /*--- Open the restart file, and throw an error if this fails. ---*/
  
  restart_file.open(restart_filename.data(), ios::in);
  if (restart_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  
  map<unsigned long,unsigned long> Global2Local;
  map<unsigned long,unsigned long>::const_iterator MI;
  
  /*--- Now fill array with the transform values only for local points ---*/
  
  for (iPoint = 0; iPoint < geometry[MESH_0]->GetnPointDomain(); iPoint++) {
    Global2Local[geometry[MESH_0]->node[iPoint]->GetGlobalIndex()] = iPoint;
  }
  
  /*--- Read all lines in the restart file ---*/
  
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  
  /*--- The first line is the header ---*/
  
  getline (restart_file, text_line);
  
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    
    getline (restart_file, text_line);
    
    istringstream point_line(text_line);
    
    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/
    
    MI = Global2Local.find(iPoint_Global);
    if (MI != Global2Local.end()) {
      
      iPoint_Local = Global2Local[iPoint_Global];
      
      if (nDim == 2) point_line >> index >> Coord[0] >> Coord[1] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
      if (nDim == 3) point_line >> index >> Coord[0] >> Coord[1] >> Coord[2] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
      
      node[iPoint_Local]->SetSolution(Solution);
      
      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/
      
      if (grid_movement) {
        
        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/
        
        if (turb_model == SA || turb_model == SA_NEG) {
          point_line >> dull_val;
        } else if (turb_model == SST) {
          point_line >> dull_val >> dull_val;
        }
        
        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/
        
        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {
            if (nDim == 2) point_line >> GridVel[0] >> GridVel[1];
            else point_line >> GridVel[0] >> GridVel[1] >> GridVel[2];
        }

        
        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[MESH_0]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
        
      }
      
    }

  }
  
  /*--- Close the restart file ---*/
  
  restart_file.close();
  
  /*--- MPI solution ---*/
  
  solver[MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  
  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
  }
  
  /*--- Update the geometry for flows on dynamic meshes ---*/
  
  if (grid_movement) {
    
    /*--- Communicate the new coordinates and grid velocities at the halos ---*/
    
    geometry[MESH_0]->Set_MPI_Coord(config);
    geometry[MESH_0]->Set_MPI_GridVel(config);
    
    /*--- Recompute the edges and  dual mesh control volumes in the
     domain and on the boundaries. ---*/
    
    geometry[MESH_0]->SetCoord_CG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);
    
    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/
    
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine], config);
    }
  }
  
  delete [] Coord;
  
}

void CEulerSolver::SetFreeStream_Solution(CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetSolution(0, Density_Inf);
    for (iDim = 0; iDim < nDim; iDim++) {
      node[iPoint]->SetSolution(iDim+1, Density_Inf*Velocity_Inf[iDim]);
    }
    node[iPoint]->SetSolution(nVar-1, Density_Inf*Energy_Inf);
  }
}

CNSSolver::CNSSolver(void) : CEulerSolver() {
  
  /*--- Basic array initialization ---*/
  
  CD_Visc = NULL; CL_Visc = NULL; CSF_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  
  ForceViscous = NULL; MomentViscous = NULL; CSkinFriction = NULL;
  
  /*--- Surface based array initialization ---*/
  
  Surface_CL_Visc = NULL; Surface_CD_Visc = NULL; Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  Surface_HF_Visc = NULL; Surface_MaxHF_Visc = NULL;
  
  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Visc = NULL; CT_Visc = NULL; CQ_Visc = NULL;
  
}

CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CEulerSolver() {
  
  unsigned long iPoint, index, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  su2double Density, Velocity2, Pressure, Temperature, dull_val, StaticEnergy;
  int Unst_RestartIter;
  ifstream restart_file;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  string filename = config->GetSolution_FlowFileName();
  string filename_ = config->GetSolution_FlowFileName();
  su2double AoA_, AoS_, BCThrust_;
  string::size_type position;
  unsigned long ExtIter_;

  unsigned short direct_diff = config->GetDirectDiff();
  unsigned short nMarkerTurboPerf = config->Get_nMarkerTurboPerf();
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Check for a restart file to check if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/
  
  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {
    
    /*--- Multizone problems require the number of the zone to be appended. ---*/
    
    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone);
    
    /*--- Modify file name for a dual-time unsteady restart ---*/
    
    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }
    
    /*--- Modify file name for a simple unsteady restart ---*/
    
    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Open the restart file, throw an error if this fails. ---*/
    
    restart_file.open(filename_.data(), ios::in);
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file!! " << filename_.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }
    
    unsigned long iPoint_Global = 0;
    string text_line;
    
    /*--- The first line is the header (General description) ---*/
    
    getline (restart_file, text_line);
    
    /*--- Space for the solution ---*/
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
    }
    
    /*--- Space for extra info (if any) ---*/
    
    while (getline (restart_file, text_line)) {
      
      /*--- Angle of attack ---*/
      
      position = text_line.find ("AOA=",0);
      if (position != string::npos) {
        text_line.erase (0,4); AoA_ = atof(text_line.c_str());
        if (config->GetDiscard_InFiles() == false) {
          if ((config->GetAoA() != AoA_) &&  (rank == MASTER_NODE)) {
            cout.precision(6);
            cout << fixed <<"WARNING: AoA in the solution file (" << AoA_ << " deg.) +" << endl;
            cout << "         AoA offset in mesh file (" << config->GetAoA_Offset() << " deg.) = " << AoA_ + config->GetAoA_Offset() << " deg." << endl;
          }
          config->SetAoA(AoA_ + config->GetAoA_Offset());
        }
        else {
          if ((config->GetAoA() != AoA_) &&  (rank == MASTER_NODE))
            cout <<"WARNING: Discarding the AoA in the solution file." << endl;
        }
      }
      
      /*--- Sideslip angle ---*/
      
      position = text_line.find ("SIDESLIP_ANGLE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); AoS_ = atof(text_line.c_str());
        if (config->GetDiscard_InFiles() == false) {
          if ((config->GetAoS() != AoS_) &&  (rank == MASTER_NODE)) {
            cout.precision(6);
            cout << fixed <<"WARNING: AoS in the solution file (" << AoS_ << " deg.) +" << endl;
            cout << "         AoS offset in mesh file (" << config->GetAoS_Offset() << " deg.) = " << AoS_ + config->GetAoS_Offset() << " deg." << endl;
          }
          config->SetAoS(AoS_ + config->GetAoS_Offset());
        }
        else {
          if ((config->GetAoS() != AoS_) &&  (rank == MASTER_NODE))
            cout <<"WARNING: Discarding the AoS in the solution file." << endl;
        }
      }
      
      /*--- BCThrust angle ---*/
      
      position = text_line.find ("INITIAL_BCTHRUST=",0);
      if (position != string::npos) {
        text_line.erase (0,17); BCThrust_ = atof(text_line.c_str());
        if (config->GetDiscard_InFiles() == false) {
          if ((config->GetInitial_BCThrust() != BCThrust_) &&  (rank == MASTER_NODE))
            cout <<"WARNING: ACDC will use the initial BC Thrust provided in the solution file: " << BCThrust_ << " lbs." << endl;
          config->SetInitial_BCThrust(BCThrust_);
        }
        else {
          if ((config->GetInitial_BCThrust() != BCThrust_) &&  (rank == MASTER_NODE))
            cout <<"WARNING: Discarding the BC Thrust in the solution file." << endl;
        }
      }
      
      /*--- External iteration ---*/
      
      position = text_line.find ("EXT_ITER=",0);
      if (position != string::npos) {
        text_line.erase (0,9); ExtIter_ = atoi(text_line.c_str());
        if (!config->GetContinuous_Adjoint() && !config->GetDiscrete_Adjoint())
          config->SetExtIter_OffSet(ExtIter_);
      }
      
    }
    
    /*--- Close the restart file... we will open this file again... ---*/
    
    restart_file.close();
    
  }

  /*--- Array initialization ---*/
  
  CD_Visc = NULL; CL_Visc = NULL; CSF_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  
  Surface_CL_Visc = NULL; Surface_CD_Visc = NULL; Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  Surface_HF_Visc = NULL; Surface_MaxHF_Visc = NULL;
  
  CMerit_Visc = NULL;      CT_Visc = NULL;      CQ_Visc = NULL;
  MaxHF_Visc = NULL; ForceViscous = NULL; MomentViscous = NULL;
  CSkinFriction = NULL;    Cauchy_Serie = NULL; HF_Visc = NULL;
  
  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables (T, vx, vy, vz, P, rho, h, c, lamMu, EddyMu, ThCond, Cp).
   ---*/
  
  nDim = geometry->GetnDim();
  
  nVar = nDim+2;
  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  nSecondaryVar = 8; nSecondaryVarGrad = 2;

  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nPrimVarGrad;
  
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
 
  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];
 
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  SetNondimensionalization(geometry, config, iMesh);
  
  /*--- Allocate the node variables ---*/
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  
  Primitive   = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the Secondary solution ---*/
  
  Secondary   = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
  Secondary_i = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
  Secondary_j = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;

  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  
  if (roe_turkel) {
    LowMach_Precontioner = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new su2double[nVar];
  }
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  }
  
  else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Define some auxiliary vectors for computing flow variable
   gradients by least squares, S matrix := inv(R)*traspose(inv(R)),
   c vector := transpose(WA)*(Wb) ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    Cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  CharacPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/
  
  DonorPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      if (rans) {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar+2];
        for (iVar = 0; iVar < nPrimVar + 2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
        for (iVar = 0; iVar < nPrimVar ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
    }
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Delta P at the Actuator Disk ---*/
  
  ActDisk_DeltaP = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaP[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaP[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Delta T at the Actuator Disk ---*/
  
  ActDisk_DeltaT = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ActDisk_DeltaT[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      ActDisk_DeltaT[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Total Pressure at the inlet BC ---*/
  
  Inlet_Ttotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_Ttotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_Ttotal[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Total Temperature at the inlet BC ---*/
  
  Inlet_Ptotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_Ptotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_Ptotal[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Flow direction at the inlet BC ---*/
  
  Inlet_FlowDir = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_FlowDir[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_FlowDir[iMarker][iVertex] = new su2double [nDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        Inlet_FlowDir[iMarker][iVertex][iDim] = 0;
      }
    }
  }

  /*--- Inviscid force definition and coefficient in all the markers ---*/
  
  CPressure = new su2double* [nMarker];
  CPressureTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new su2double [geometry->nVertex[iMarker]];
    CPressureTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Heat flux in all the markers ---*/
  
  HeatFlux = new su2double* [nMarker];
  HeatFluxTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatFlux[iMarker] = new su2double [geometry->nVertex[iMarker]];
    HeatFluxTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      HeatFlux[iMarker][iVertex] = 0.0;
      HeatFluxTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Y plus in all the markers ---*/
  
  YPlus = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    YPlus[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      YPlus[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Skin friction in all the markers ---*/
  
  CSkinFriction = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSkinFriction[iMarker] = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      CSkinFriction[iMarker][iDim] = new su2double[geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        CSkinFriction[iMarker][iDim][iVertex] = 0.0;
      }
    }
  }
  
  /*--- Non dimensional coefficients ---*/
  
  ForceInviscid  = new su2double[3];
  MomentInviscid = new su2double[3];
  CD_Inv         = new su2double[nMarker];
  CL_Inv         = new su2double[nMarker];
  CSF_Inv        = new su2double[nMarker];
  CMx_Inv        = new su2double[nMarker];
  CMy_Inv        = new su2double[nMarker];
  CMz_Inv        = new su2double[nMarker];
  CEff_Inv       = new su2double[nMarker];
  CFx_Inv        = new su2double[nMarker];
  CFy_Inv        = new su2double[nMarker];
  CFz_Inv        = new su2double[nMarker];
  
  ForceMomentum  = new su2double[3];
  MomentMomentum = new su2double[3];
  CD_Mnt         = new su2double[nMarker];
  CL_Mnt         = new su2double[nMarker];
  CSF_Mnt        = new su2double[nMarker];
  CMx_Mnt        = new su2double[nMarker];
  CMy_Mnt        = new su2double[nMarker];
  CMz_Mnt        = new su2double[nMarker];
  CEff_Mnt       = new su2double[nMarker];
  CFx_Mnt        = new su2double[nMarker];
  CFy_Mnt        = new su2double[nMarker];
  CFz_Mnt        = new su2double[nMarker];

  ForceViscous     = new su2double[3];
  MomentViscous    = new su2double[3];
  CD_Visc          = new su2double[nMarker];
  CL_Visc          = new su2double[nMarker];
  CSF_Visc         = new su2double[nMarker];
  CMx_Visc         = new su2double[nMarker];
  CMy_Visc         = new su2double[nMarker];
  CMz_Visc         = new su2double[nMarker];
  CEff_Visc        = new su2double[nMarker];
  CFx_Visc         = new su2double[nMarker];
  CFy_Visc         = new su2double[nMarker];
  CFz_Visc         = new su2double[nMarker];
  
  Surface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Inv     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  
  Surface_CL_Mnt         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Mnt         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF         = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];
  
  Surface_CL_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Visc     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_HF_Visc         = new su2double[config->GetnMarker_Monitoring()];
  Surface_MaxHF_Visc      = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Rotational coefficients ---*/
  
  CMerit_Inv = new su2double[nMarker];
  CT_Inv     = new su2double[nMarker];
  CQ_Inv     = new su2double[nMarker];
  
  CMerit_Mnt = new su2double[nMarker];
  CT_Mnt     = new su2double[nMarker];
  CQ_Mnt     = new su2double[nMarker];

  CMerit_Visc = new su2double[nMarker];
  CT_Visc     = new su2double[nMarker];
  CQ_Visc     = new su2double[nMarker];
  
  /*--- Heat based coefficients ---*/
  
  HF_Visc    = new su2double[nMarker];
  MaxHF_Visc = new su2double[nMarker];
  
  /*--- Supersonic coefficients ---*/
  
  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];
  
  /*--- Engine simulation ---*/
  
  Inflow_MassFlow     = new su2double[nMarker];
  Inflow_Pressure     = new su2double[nMarker];
  Inflow_Mach         = new su2double[nMarker];
  Inflow_Area         = new su2double[nMarker];
  
  Exhaust_MassFlow    = new su2double[nMarker];
  Exhaust_Pressure    = new su2double[nMarker];
  Exhaust_Temperature = new su2double[nMarker];
  Exhaust_Area        = new su2double[nMarker];
  
  /*--- Init total coefficients ---*/
  
  Total_CD         = 0.0;    Total_CL           = 0.0;    Total_CSF          = 0.0;
  Total_CMx        = 0.0;    Total_CMy          = 0.0;    Total_CMz          = 0.0;
  Total_CEff       = 0.0;    Total_CEquivArea   = 0.0;    Total_CNearFieldOF = 0.0;
  Total_CFx        = 0.0;    Total_CFy          = 0.0;    Total_CFz          = 0.0;
  Total_CT         = 0.0;    Total_CQ           = 0.0;    Total_CMerit       = 0.0;
  Total_MaxHeat    = 0.0;   Total_Heat         = 0.0;    Total_ComboObj     = 0.0;
  Total_CpDiff     = 0.0;   Total_HeatFluxDiff = 0.0;    Total_BCThrust_Prev = 0.0;
  Total_NetCThrust = 0.0;   Total_NetCThrust_Prev = 0.0; Total_CL_Prev = 0.0;
  Total_Power      = 0.0;   AoA_Prev           = 0.0;    Total_CD_Prev      = 0.0;
  Total_AeroCD     = 0.0;    Total_RadialDistortion   = 0.0; Total_CircumferentialDistortion           = 0.0;

  /*--- Read farfield conditions from config ---*/
  
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Mach_Inf        = config->GetMach();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  Tke_Inf         = config->GetTke_FreeStreamND();
  
  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  
  switch(direct_diff) {
    case NO_DERIVATIVE:
      break;
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
      break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
      break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
      break;
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Already done in postprocessing of config ---*/
      break;
    default:
      break;
  }
  
  /*--- Initializate fan face pressure, fan face mach number, and mass flow rate ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;
    
    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;

  }
  
  /*--- Initializate quantities for the mixing process*/
  
  AveragedVelocity = new su2double* [nMarker];
  AveragedNormal = new su2double* [nMarker];
  AveragedGridVel = new su2double* [nMarker];
  AveragedFlux = new su2double* [nMarker];
  TotalFlux = new su2double* [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AveragedVelocity[iMarker] = new su2double [nDim];
    AveragedNormal[iMarker] = new su2double [nDim];
    AveragedGridVel[iMarker] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      AveragedVelocity[iMarker][iDim] = 0.0;
      AveragedNormal[iMarker][iDim] = 0.0;
      AveragedGridVel[iMarker][iDim] = 0.0;
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AveragedFlux[iMarker] = new su2double [nVar];
    TotalFlux[iMarker] = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      AveragedFlux[iMarker][iVar] = 0.0;
      TotalFlux[iMarker][iVar] = 0.0;
    }
  }
  
  AveragedNormalVelocity = new su2double[nMarker];
  AveragedTangVelocity = new su2double[nMarker];
  ExtAveragedNormalVelocity = new su2double[nMarker];
  ExtAveragedTangVelocity = new su2double[nMarker];
  MassFlow= new su2double[nMarker];
  FlowAngle= new su2double[nMarker];
  AveragedEnthalpy  = new su2double[nMarker];
  AveragedPressure  = new su2double[nMarker];
  AveragedTotPressure  = new su2double[nMarker];
  AveragedTotTemperature  = new su2double[nMarker];
  ExtAveragedTotPressure  = new su2double[nMarker];
  ExtAveragedTotTemperature  = new su2double[nMarker];
  ExtAveragedPressure  = new su2double[nMarker];
  AveragedDensity   = new su2double[nMarker];
  ExtAveragedDensity   = new su2double[nMarker];
  AveragedSoundSpeed= new su2double[nMarker];
  AveragedEntropy   = new su2double[nMarker];
  AveragedTangGridVelocity = new su2double[nMarker];
  AveragedMach = new su2double[nMarker];
  AveragedNormalMach = new su2double[nMarker];
  AveragedTangMach = new su2double[nMarker];
  
  
  /*--- Initializate quantities for SlidingMesh Interface ---*/
  
  SlidingState = new su2double** [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
  SlidingState[iMarker] = NULL;

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

        SlidingState[iMarker] = new su2double* [geometry->GetnVertex(iMarker)];

        for (iPoint = 0; iPoint < geometry->nVertex[iMarker]; iPoint++) {
          SlidingState[iMarker][iPoint] = new su2double[nPrimVar];
        for (iVar = 0; iVar < nVar; iVar++)
          SlidingState[iMarker][iPoint][iVar] = -1;
      }
    }
    else
      SlidingState[iMarker] = NULL;
  }
  
  
  /*--- Initializate quantities for turboperformace ---*/
  
  TotalStaticEfficiency = new su2double[nMarkerTurboPerf];
  TotalTotalEfficiency = new su2double[nMarkerTurboPerf];
  KineticEnergyLoss= new su2double[nMarkerTurboPerf];
  TotalPressureLoss= new su2double[nMarkerTurboPerf];
  MassFlowIn= new su2double[nMarkerTurboPerf];
  MassFlowOut= new su2double[nMarkerTurboPerf];
  FlowAngleIn= new su2double[nMarkerTurboPerf];
  FlowAngleOut= new su2double[nMarkerTurboPerf];
  EulerianWork= new su2double[nMarkerTurboPerf];
  TotalEnthalpyIn= new su2double[nMarkerTurboPerf];
  PressureRatio= new su2double[nMarkerTurboPerf];
  PressureOut= new su2double[nMarkerTurboPerf];
  EnthalpyOut= new su2double[nMarkerTurboPerf];
  MachIn= new su2double[nMarkerTurboPerf];
  MachOut= new su2double[nMarkerTurboPerf];
  NormalMachIn= new su2double[nMarkerTurboPerf];
  NormalMachOut= new su2double[nMarkerTurboPerf];
  VelocityOutIs= new su2double[nMarkerTurboPerf];
  
  for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++) {
    TotalStaticEfficiency[iMarker]= 0.0;
    TotalTotalEfficiency[iMarker]= 0.0;
    KineticEnergyLoss[iMarker]= 0.0;
    TotalPressureLoss[iMarker]= 0.0;
    MassFlowIn[iMarker]= 0.0;
    MassFlowOut[iMarker]= 0.0;
    FlowAngleIn[iMarker]= 0.0;
    FlowAngleOut[iMarker]= 0.0;
    EulerianWork[iMarker]= 0.0;
    TotalEnthalpyIn[iMarker]= 0.0;
    PressureRatio[iMarker]= 0.0;
    PressureOut[iMarker]= 0.0;
    EnthalpyOut[iMarker]= 0.0;
    MachIn[iMarker]= 0.0;
    MachOut[iMarker]= 0.0;
    NormalMachIn[iMarker]= 0.0;
    NormalMachOut[iMarker]= 0.0;
    VelocityOutIs[iMarker]= 0.0;
  }
  
  
  /*--- Initialize the cauchy critera array for fixed CL mode ---*/
  
  if (config->GetFixed_CL_Mode())
    
    Cauchy_Serie = new su2double [config->GetCauchy_Elems()+1];
  
  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  
  if (!restart || (iMesh != MESH_0)) {
    
    /*--- Restart the solution from the free-stream state ---*/
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
    
  }
  
  else {
    
    /*--- Multizone problems require the number of the zone to be appended. ---*/
    
    if (nZone > 1) filename = config->GetMultizone_FileName(filename, iZone);
    
    /*--- Modify file name for a dual-time unsteady restart ---*/
    
    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }
    
    /*--- Modify file name for a simple unsteady restart ---*/
    
    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }

    /*--- Open the restart file, throw an error if this fails. ---*/
    
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    
    map<unsigned long,unsigned long> Global2Local;
    map<unsigned long,unsigned long>::const_iterator MI;
    
    /*--- Now fill array with the transform values only for local points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    
    /*--- Read all lines in the restart file ---*/
    
    long iPoint_Local;
    unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
    
    /*--- The first line is the header ---*/
    
    getline (restart_file, text_line);
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
      istringstream point_line(text_line);
      
      if (iPoint_Global >= geometry->GetGlobal_nPointDomain()) { sbuf_NotMatching = 1; break; }
      
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/
      
      MI = Global2Local.find(iPoint_Global);
      if (MI != Global2Local.end()) {
        
        iPoint_Local = Global2Local[iPoint_Global];
        
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];

        node[iPoint_Local] = new CNSVariable(Solution, nDim, nVar, config);
        iPoint_Global_Local++;
      }

    }
    
    /*--- Detect a wrong solution file ---*/
    
    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
    
#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
        cout << "It could be empty lines at the end of the file." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Solution, nDim, nVar, config);
    
    /*--- Close the restart file ---*/
    
    restart_file.close();
    
  }
  
  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
    
  counter_local = 0;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    Density = node[iPoint]->GetSolution(0);

    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);

    StaticEnergy= node[iPoint]->GetSolution(nDim+1)/Density - 0.5*Velocity2;

    FluidModel->SetTDState_rhoe(Density, StaticEnergy);
    Pressure= FluidModel->GetPressure();
    Temperature= FluidModel->GetTemperature();

    /*--- Use the values at the infinity ---*/

    if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);
      counter_local++;
    }

  }

  /*--- Warning message about non-physical points ---*/

  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }
  
  /*--- Define solver parameters needed for execution of destructor ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) space_centered = true;
  else space_centered = false;
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
}

CNSSolver::~CNSSolver(void) {
  unsigned short iMarker, iDim;
  unsigned long iVertex;
  
  if (CD_Visc != NULL)          delete [] CD_Visc;
  if (CL_Visc != NULL)          delete [] CL_Visc;
  if (CSF_Visc != NULL)         delete [] CSF_Visc;
  if (CMx_Visc != NULL)         delete [] CMx_Visc;
  if (CMy_Visc != NULL)         delete [] CMy_Visc;
  if (CMz_Visc != NULL)         delete [] CMz_Visc;
  if (CFx_Visc != NULL)         delete [] CFx_Visc;
  if (CFy_Visc != NULL)         delete [] CFy_Visc;
  if (CFz_Visc != NULL)         delete [] CFz_Visc;
  if (CEff_Visc != NULL)        delete [] CEff_Visc;
  if (CMerit_Visc != NULL)      delete [] CMerit_Visc;
  if (CT_Visc != NULL)          delete [] CT_Visc;
  if (CQ_Visc != NULL)          delete [] CQ_Visc;
  if (HF_Visc != NULL)          delete [] HF_Visc;
  if (MaxHF_Visc != NULL)       delete [] MaxHF_Visc;
  if (ForceViscous != NULL)     delete [] ForceViscous;
  if (MomentViscous != NULL)    delete [] MomentViscous;
  
  if (Surface_CL_Visc != NULL)      delete [] Surface_CL_Visc;
  if (Surface_CD_Visc != NULL)      delete [] Surface_CD_Visc;
  if (Surface_CSF_Visc != NULL)     delete [] Surface_CSF_Visc;
  if (Surface_CEff_Visc != NULL)    delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != NULL)     delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != NULL)     delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != NULL)     delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != NULL)     delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)     delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)     delete [] Surface_CMz_Visc;
  if (Surface_HF_Visc != NULL)      delete [] Surface_HF_Visc;
  if (Surface_MaxHF_Visc != NULL)   delete [] Surface_MaxHF_Visc;
  
  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;

  
  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        delete [] CSkinFriction[iMarker][iDim];
      }
      delete [] CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }
  
  if (Inlet_Ttotal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] Inlet_Ttotal[iMarker];
    delete [] Inlet_Ttotal;
  }
  
  if (Inlet_Ptotal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] Inlet_Ptotal[iMarker];
    delete [] Inlet_Ptotal;
  }
  
  if (Inlet_FlowDir != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] Inlet_FlowDir[iMarker][iVertex];
      delete [] Inlet_FlowDir[iMarker];
    }
    delete [] Inlet_FlowDir;
  }

}

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double StrainMag = 0.0, Omega = 0.0, *Vorticity;
    
  unsigned long ExtIter     = config->GetExtIter();
  bool adjoint              = config->GetContinuous_Adjoint();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center               = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst           = center && config->GetKind_Centered_Flow() == JST;
  bool limiter_flow         = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_turb         = ((config->GetSpatialOrder_Turb() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_adjflow      = ((config->GetSpatialOrder_AdjFlow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_visc         = config->GetViscous_Limiter_Flow();
  bool fixed_cl             = config->GetFixed_CL_Mode();
  bool engine               = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk        = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool nearfield            = (config->GetnMarker_NearFieldBound() != 0);
  bool interface            = (config->GetnMarker_InterfaceBound() != 0);
  bool marker_analyze       = (config->GetnMarker_Analyze() != 0);

  /*--- Update the angle of attack at the far-field for fixed CL calculations. ---*/
  
  if (fixed_cl) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Compute the engine properties ---*/

  if (engine) { GetPower_Properties(geometry, config, iMesh, Output); }

  /*--- Compute the control volume properties ---*/

  if (marker_analyze) {
     GetSurface_Properties(geometry, NULL, NULL, config, iMesh, Output);
    GetSurface_Distortion(geometry, config, iMesh, Output);
  }

  /*--- Compute the actuator disk properties and distortion levels ---*/

  if (actuator_disk) {
    Set_MPI_ActDisk(solver_container, geometry, config);
    SetActDisk_BCThrust(geometry, solver_container, config, iMesh, Output);
  }

  /*--- Compute Interface MPI ---*/

  if (interface) { Set_MPI_Interface(geometry, config); }

  /*--- Compute NearField MPI ---*/

  if (nearfield) { Set_MPI_Nearfield(geometry, config); }
 
  /*--- Artificial dissipation ---*/

  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Compute gradient of the primitive variables ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
    //    if (compressible && !ideal_gas) SetSecondary_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
    //    if (compressible && !ideal_gas) SetSecondary_Gradient_LS(geometry, config);
  }

  /*--- Compute the limiter in case we need it in the turbulence model
   or to limit the viscous terms (check this logic with JST and 2nd order turbulence model) ---*/
  
  if ((iMesh == MESH_0) && (limiter_flow || limiter_turb || limiter_adjflow || limiter_visc) && !Output) { SetPrimitive_Limiter(geometry, config);
    //  if (compressible && !ideal_gas) SetSecondary_Limiter(geometry, config);
  }
  
  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  
  StrainMag_Max = 0.0, Omega_Max = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    solver_container[FLOW_SOL]->node[iPoint]->SetVorticity(limiter_visc);
    solver_container[FLOW_SOL]->node[iPoint]->SetStrainMag(limiter_visc);
    
    StrainMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
    Vorticity = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
    Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);
    
    StrainMag_Max = max(StrainMag_Max, StrainMag);
    Omega_Max = max(Omega_Max, Omega);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !config->GetDiscrete_Adjoint()) Jacobian.SetValZero();

  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    su2double MyOmega_Max = Omega_Max; Omega_Max = 0.0;
    su2double MyStrainMag_Max = StrainMag_Max; StrainMag_Max = 0.0;
    
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    if (iMesh == MESH_0) {
      config->SetNonphysical_Points(ErrorCounter);
      solver_container[FLOW_SOL]->SetStrainMag_Max(StrainMag_Max);
      solver_container[FLOW_SOL]->SetOmega_Max(Omega_Max);
    }
    
  }
  
}

unsigned long CNSSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double eddy_visc = 0.0, turb_ke = 0.0;
  unsigned short turb_model = config->GetKind_Turb_Model();
  bool RightSol = true;
  
  bool tkeNeeded            = (turb_model == SST);

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Retrieve the value of the kinetic energy (if need it) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
    }
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
    
    RightSol = node[iPoint]->SetPrimVar(eddy_visc, turb_ke, FluidModel);
    node[iPoint]->SetSecondaryVar(FluidModel);
    
    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
}

void CNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
  
  su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time, Local_Delta_Time_Visc,
  Global_Delta_Time = 1E6, Mean_LaminarVisc = 0.0, Mean_EddyVisc = 0.0, Mean_Density = 0.0, Lambda_1, Lambda_2, K_v = 0.25, Global_Delta_UnstTimeND;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  unsigned short iDim, iMarker;
  su2double ProjVel, ProjVel_i, ProjVel_j;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetMax_Lambda_Inv(0.0);
    node[iPoint]->SetMax_Lambda_Visc(0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j) ;
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
    /*--- Viscous contribution ---*/
    
    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosity() + node[jPoint]->GetEddyViscosity());
    Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));
    
    Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
    //TODO (REAL_GAS) removing Gamma it cannot work with FLUIDPROP
    Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
    Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
    
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      
      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
      
      /*--- Viscous contribution ---*/
      
      Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity();
      Mean_EddyVisc    = node[iPoint]->GetEddyViscosity();
      Mean_Density     = node[iPoint]->GetSolution(0);
      
      Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
      Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
      Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
      
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
      Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
    
  }
  
  
  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;
    
    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }
  
  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }
  
  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
    
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }
  
  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points, coordinates and normal vector in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive and secondary variables ---*/
    
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    numerics->SetSecondary(node[iPoint]->GetSecondary(), node[jPoint]->GetSecondary());
    
    /*--- Gradient and limiters ---*/
    
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());
    numerics->SetPrimVarLimiter(node[iPoint]->GetLimiter_Primitive(), node[jPoint]->GetLimiter_Primitive());
    
    /*--- Turbulent kinetic energy ---*/
    
    if (config->GetKind_Turb_Model() == SST)
      numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                     solver_container[TURB_SOL]->node[jPoint]->GetSolution(0));
    
    /*--- Compute and update residual ---*/
    
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.SubtractBlock(iPoint, Res_Visc);
    LinSysRes.AddBlock(jPoint, Res_Visc);
    
    /*--- Implicit part ---*/
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    }
    
  }
  
}

void CNSSolver::Friction_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint, iPointNormal;
  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
  su2double Viscosity = 0.0, div_vel, *Normal, MomentDist[3] = {0.0, 0.0, 0.0}, WallDist[3] = {0.0, 0.0, 0.0},
  *Coord, *Coord_Normal, Area, WallShearStress, TauNormal, factor, RefTemp, RefVel2,
  RefDensity, GradTemperature, Density = 0.0, WallDistMod, FrictionVel,
  Mach2Vel, Mach_Motion, UnitNormal[3] = {0.0, 0.0, 0.0}, TauElem[3] = {0.0, 0.0, 0.0}, TauTangent[3] = {0.0, 0.0, 0.0},
  Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Force[3] = {0.0, 0.0, 0.0}, Cp, thermal_conductivity, MaxNorm = 8.0,
  Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Grad_Temp[3] = {0.0, 0.0, 0.0},
  delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double AxiFactor;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Visc, MyAllBound_CL_Visc, MyAllBound_CSF_Visc, MyAllBound_CMx_Visc, MyAllBound_CMy_Visc, MyAllBound_CMz_Visc, MyAllBound_CFx_Visc, MyAllBound_CFy_Visc, MyAllBound_CFz_Visc, MyAllBound_CT_Visc, MyAllBound_CQ_Visc, MyAllBound_HF_Visc, MyAllBound_MaxHF_Visc, *MySurface_CL_Visc = NULL, *MySurface_CD_Visc = NULL, *MySurface_CSF_Visc = NULL, *MySurface_CEff_Visc = NULL, *MySurface_CFx_Visc = NULL, *MySurface_CFy_Visc = NULL, *MySurface_CFz_Visc = NULL, *MySurface_CMx_Visc = NULL, *MySurface_CMy_Visc = NULL, *MySurface_CMz_Visc = NULL, *MySurface_HF_Visc = NULL, *MySurface_MaxHF_Visc;
#endif
  
  string Marker_Tag, Monitoring_Tag;
  
  su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefAreaCoeff    = config->GetRefAreaCoeff();
  su2double RefLengthMoment = config->GetRefLengthMoment();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double *Origin         = config->GetRefOriginMoment(0);
  bool grid_movement        = config->GetGrid_Movement();
  su2double Prandtl_Lam     = config->GetPrandtl_Lam();
  bool axisymmetric         = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  RefTemp    = Temperature_Inf;
  RefDensity = Density_Inf;
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  } else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  /*--- Variables initialization ---*/
  
  AllBound_CD_Visc = 0.0;    AllBound_CL_Visc = 0.0;       AllBound_CSF_Visc = 0.0;
  AllBound_CMx_Visc = 0.0;      AllBound_CMy_Visc = 0.0;         AllBound_CMz_Visc = 0.0;
  AllBound_CFx_Visc = 0.0;      AllBound_CFy_Visc = 0.0;         AllBound_CFz_Visc = 0.0;
  AllBound_CT_Visc = 0.0;       AllBound_CQ_Visc = 0.0;          AllBound_CMerit_Visc = 0.0;
  AllBound_HF_Visc = 0.0; AllBound_MaxHF_Visc = 0.0; AllBound_CEff_Visc = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Visc[iMarker_Monitoring]      = 0.0; Surface_CD_Visc[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring] = 0.0; Surface_CEff_Visc[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]        = 0.0; Surface_CFy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]        = 0.0; Surface_CMx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]        = 0.0; Surface_CMz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_HF_Visc[iMarker_Monitoring]              = 0.0; Surface_MaxHF_Visc[iMarker_Monitoring]           = 0.0;
  }
  
  /*--- Loop over the Navier-Stokes markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CD_Visc[iMarker] = 0.0; CL_Visc[iMarker] = 0.0;       CSF_Visc[iMarker] = 0.0;
      CMx_Visc[iMarker] = 0.0;   CMy_Visc[iMarker] = 0.0;         CMz_Visc[iMarker] = 0.0;
      CFx_Visc[iMarker] = 0.0;   CFy_Visc[iMarker] = 0.0;         CFz_Visc[iMarker] = 0.0;
      CT_Visc[iMarker] = 0.0;    CQ_Visc[iMarker] = 0.0;          CMerit_Visc[iMarker] = 0.0;
      HF_Visc[iMarker] = 0.0;  MaxHF_Visc[iMarker] = 0.0; CEff_Visc[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
      MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;
      
      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
        
        Coord = geometry->node[iPoint]->GetCoord();
        Coord_Normal = geometry->node[iPointNormal]->GetCoord();
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
          Grad_Temp[iDim] = node[iPoint]->GetGradient_Primitive(0, iDim);
        }
        
        Viscosity = node[iPoint]->GetLaminarViscosity();
        Density = node[iPoint]->GetDensity();
        
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
        

        for (iDim = 0; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim]/Area;
          MomentDist[iDim] = Coord[iDim] - Origin[iDim];
        }
        
        /*--- Evaluate Tau ---*/
        
        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*delta[iDim][jDim];
          }
        }
        
        /*--- Project Tau in each surface element ---*/
        
        for (iDim = 0; iDim < nDim; iDim++) {
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++) {
            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
          }
        }
        
        /*--- Compute wall shear stress (using the stress tensor). Compute wall skin friction coefficient, and heat flux on the wall ---*/
        
        TauNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitNormal[iDim];
        
        WallShearStress = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
          CSkinFriction[iMarker][iDim][iVertex] = TauTangent[iDim] / (0.5*RefDensity*RefVel2);
          WallShearStress += TauTangent[iDim] * TauTangent[iDim];
        }
        WallShearStress = sqrt(WallShearStress);
        
        for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);
        
        /*--- Compute y+ and non-dimensional velocity ---*/
        
        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);
        
        /*--- Compute total and maximum heat flux on the wall ---*/

        GradTemperature = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          GradTemperature -= Grad_Temp[iDim]*UnitNormal[iDim];

        Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
        thermal_conductivity = Cp * Viscosity/Prandtl_Lam;
        HeatFlux[iMarker][iVertex] = -thermal_conductivity*GradTemperature;
        HF_Visc[iMarker] += HeatFlux[iMarker][iVertex]*Area;
        MaxHF_Visc[iMarker] += pow(HeatFlux[iMarker][iVertex], MaxNorm);

        
        /*--- Note that y+, and heat are computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {
          
          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim] * Area * factor * AxiFactor;
            ForceViscous[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (iDim == 3) {
            MomentViscous[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLengthMoment;
            MomentViscous[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLengthMoment;
          }
          MomentViscous[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLengthMoment;
          
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if (Monitoring == YES) {
        if (nDim == 2) {
          CD_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CL_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker]        = CL_Visc[iMarker] / (CD_Visc[iMarker]+EPS);
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CT_Visc[iMarker]          = -CFx_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker] / (CQ_Visc[iMarker]+EPS);
          MaxHF_Visc[iMarker] = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }
        if (nDim == 3) {
          CD_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          CL_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSF_Visc[iMarker]  = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker]        = CL_Visc[iMarker]/(CD_Visc[iMarker] + EPS);
          CMx_Visc[iMarker]         = MomentViscous[0];
          CMy_Visc[iMarker]         = MomentViscous[1];
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CFz_Visc[iMarker]         = ForceViscous[2];
          CT_Visc[iMarker]          = -CFz_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker] / (CQ_Visc[iMarker] + EPS);
          MaxHF_Visc[iMarker] = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }
        
        AllBound_CD_Visc       += CD_Visc[iMarker];
        AllBound_CL_Visc       += CL_Visc[iMarker];
        AllBound_CSF_Visc  += CSF_Visc[iMarker];
        AllBound_CMx_Visc         += CMx_Visc[iMarker];
        AllBound_CMy_Visc         += CMy_Visc[iMarker];
        AllBound_CMz_Visc         += CMz_Visc[iMarker];
        AllBound_CFx_Visc         += CFx_Visc[iMarker];
        AllBound_CFy_Visc         += CFy_Visc[iMarker];
        AllBound_CFz_Visc         += CFz_Visc[iMarker];
        AllBound_CT_Visc          += CT_Visc[iMarker];
        AllBound_CQ_Visc          += CQ_Visc[iMarker];
        AllBound_HF_Visc          += HF_Visc[iMarker];
        AllBound_MaxHF_Visc       += pow(MaxHF_Visc[iMarker], MaxNorm);
        
        /*--- Compute the coefficients per surface ---*/
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Visc[iMarker_Monitoring]      += CL_Visc[iMarker];
            Surface_CD_Visc[iMarker_Monitoring]      += CD_Visc[iMarker];
            Surface_CSF_Visc[iMarker_Monitoring] += CSF_Visc[iMarker];
            Surface_CEff_Visc[iMarker_Monitoring]       += CEff_Visc[iMarker];
            Surface_CFx_Visc[iMarker_Monitoring]        += CFx_Visc[iMarker];
            Surface_CFy_Visc[iMarker_Monitoring]        += CFy_Visc[iMarker];
            Surface_CFz_Visc[iMarker_Monitoring]        += CFz_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]        += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]        += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]        += CMz_Visc[iMarker];
            Surface_HF_Visc[iMarker_Monitoring]         += HF_Visc[iMarker];
            Surface_MaxHF_Visc[iMarker_Monitoring]      += pow(MaxHF_Visc[iMarker],MaxNorm);
          }
        }
        
      }
      
    }
  }
  
  /*--- Update some global coeffients ---*/
  
  AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);
  
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  MyAllBound_CD_Visc        = AllBound_CD_Visc;                      AllBound_CD_Visc = 0.0;
  MyAllBound_CL_Visc        = AllBound_CL_Visc;                      AllBound_CL_Visc = 0.0;
  MyAllBound_CSF_Visc   = AllBound_CSF_Visc;                 AllBound_CSF_Visc = 0.0;
  AllBound_CEff_Visc = 0.0;
  MyAllBound_CMx_Visc          = AllBound_CMx_Visc;                        AllBound_CMx_Visc = 0.0;
  MyAllBound_CMy_Visc          = AllBound_CMy_Visc;                        AllBound_CMy_Visc = 0.0;
  MyAllBound_CMz_Visc          = AllBound_CMz_Visc;                        AllBound_CMz_Visc = 0.0;
  MyAllBound_CFx_Visc          = AllBound_CFx_Visc;                        AllBound_CFx_Visc = 0.0;
  MyAllBound_CFy_Visc          = AllBound_CFy_Visc;                        AllBound_CFy_Visc = 0.0;
  MyAllBound_CFz_Visc          = AllBound_CFz_Visc;                        AllBound_CFz_Visc = 0.0;
  MyAllBound_CT_Visc           = AllBound_CT_Visc;                         AllBound_CT_Visc = 0.0;
  MyAllBound_CQ_Visc           = AllBound_CQ_Visc;                         AllBound_CQ_Visc = 0.0;
  AllBound_CMerit_Visc = 0.0;
  MyAllBound_HF_Visc     = AllBound_HF_Visc;                   AllBound_HF_Visc = 0.0;
  MyAllBound_MaxHF_Visc  = pow(AllBound_MaxHF_Visc, MaxNorm);  AllBound_MaxHF_Visc = 0.0;
  
  SU2_MPI::Allreduce(&MyAllBound_CD_Visc, &AllBound_CD_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CL_Visc, &AllBound_CL_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CSF_Visc, &AllBound_CSF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CMx_Visc, &AllBound_CMx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMy_Visc, &AllBound_CMy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMz_Visc, &AllBound_CMz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFx_Visc, &AllBound_CFx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFy_Visc, &AllBound_CFy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFz_Visc, &AllBound_CFz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CT_Visc, &AllBound_CT_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CQ_Visc, &AllBound_CQ_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  SU2_MPI::Allreduce(&MyAllBound_HF_Visc, &AllBound_HF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_MaxHF_Visc, &AllBound_MaxHF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);
  
  /*--- Add the forces on the surfaces using all the nodes ---*/
  
  MySurface_CL_Visc         = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Visc         = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_HF_Visc         = new su2double[config->GetnMarker_Monitoring()];
  MySurface_MaxHF_Visc      = new su2double[config->GetnMarker_Monitoring()];
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    
    MySurface_CL_Visc[iMarker_Monitoring]      = Surface_CL_Visc[iMarker_Monitoring];
    MySurface_CD_Visc[iMarker_Monitoring]      = Surface_CD_Visc[iMarker_Monitoring];
    MySurface_CSF_Visc[iMarker_Monitoring] = Surface_CSF_Visc[iMarker_Monitoring];
    MySurface_CEff_Visc[iMarker_Monitoring]       = Surface_CEff_Visc[iMarker_Monitoring];
    MySurface_CFx_Visc[iMarker_Monitoring]        = Surface_CFx_Visc[iMarker_Monitoring];
    MySurface_CFy_Visc[iMarker_Monitoring]        = Surface_CFy_Visc[iMarker_Monitoring];
    MySurface_CFz_Visc[iMarker_Monitoring]        = Surface_CFz_Visc[iMarker_Monitoring];
    MySurface_CMx_Visc[iMarker_Monitoring]        = Surface_CMx_Visc[iMarker_Monitoring];
    MySurface_CMy_Visc[iMarker_Monitoring]        = Surface_CMy_Visc[iMarker_Monitoring];
    MySurface_CMz_Visc[iMarker_Monitoring]        = Surface_CMz_Visc[iMarker_Monitoring];
    MySurface_HF_Visc[iMarker_Monitoring]         = Surface_HF_Visc[iMarker_Monitoring];
    MySurface_MaxHF_Visc[iMarker_Monitoring]      = Surface_MaxHF_Visc[iMarker_Monitoring];
    
    Surface_CL_Visc[iMarker_Monitoring]         = 0.0;
    Surface_CD_Visc[iMarker_Monitoring]         = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CEff_Visc[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_HF_Visc[iMarker_Monitoring]         = 0.0;
    Surface_MaxHF_Visc[iMarker_Monitoring]      = 0.0;
  }
  
  SU2_MPI::Allreduce(MySurface_CL_Visc, Surface_CL_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CD_Visc, Surface_CD_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CSF_Visc, Surface_CSF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
    Surface_CEff_Visc[iMarker_Monitoring] = Surface_CL_Visc[iMarker_Monitoring] / (Surface_CD_Visc[iMarker_Monitoring] + EPS);
  SU2_MPI::Allreduce(MySurface_CFx_Visc, Surface_CFx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFy_Visc, Surface_CFy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFz_Visc, Surface_CFz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMx_Visc, Surface_CMx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMy_Visc, Surface_CMy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMz_Visc, Surface_CMz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_HF_Visc, Surface_HF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_MaxHF_Visc, Surface_MaxHF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  delete [] MySurface_CL_Visc; delete [] MySurface_CD_Visc; delete [] MySurface_CSF_Visc;
  delete [] MySurface_CEff_Visc;  delete [] MySurface_CFx_Visc;   delete [] MySurface_CFy_Visc;
  delete [] MySurface_CFz_Visc;   delete [] MySurface_CMx_Visc;   delete [] MySurface_CMy_Visc;
  delete [] MySurface_CMz_Visc;   delete [] MySurface_HF_Visc; delete [] MySurface_MaxHF_Visc;
  
#endif
  
  /*--- Update the total coefficients (note that all the nodes have the same value)---*/
  
  Total_CD          += AllBound_CD_Visc;
  Total_CL          += AllBound_CL_Visc;
  Total_CSF         += AllBound_CSF_Visc;
  Total_CEff        = Total_CL / (Total_CD + EPS);
  Total_CMx         += AllBound_CMx_Visc;
  Total_CMy         += AllBound_CMy_Visc;
  Total_CMz         += AllBound_CMz_Visc;
  Total_CFx         += AllBound_CFx_Visc;
  Total_CFy         += AllBound_CFy_Visc;
  Total_CFz         += AllBound_CFz_Visc;
  Total_CT          += AllBound_CT_Visc;
  Total_CQ          += AllBound_CQ_Visc;
  Total_CMerit      = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  Total_Heat        = AllBound_HF_Visc;
  Total_MaxHeat     = AllBound_MaxHF_Visc;
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]      += Surface_CL_Visc[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]      += Surface_CD_Visc[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring] += Surface_CSF_Visc[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CL[iMarker_Monitoring] / (Surface_CD[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        += Surface_CFx_Visc[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        += Surface_CFy_Visc[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        += Surface_CFz_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        += Surface_CMz_Visc[iMarker_Monitoring];
  }
  
}

void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, jDim, iVar, jVar;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  su2double Wall_HeatFlux, dist_ij, *Coord_i, *Coord_j, theta2;
  su2double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  su2double ProjGridVel, *GridVel, GridVel2, *Normal, Area, Pressure = 0.0;
  su2double total_viscosity, div_vel, Density, tau_vel[3] = {0.0, 0.0, 0.0}, UnitNormal[3] = {0.0, 0.0, 0.0};
  su2double laminar_viscosity = 0.0, eddy_viscosity = 0.0, Grad_Vel[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}},
  tau[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Identify the boundary by string name ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Get the specified wall heat flux from config ---*/
  
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);
  
  /*--- Loop over all of the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Compute dual-grid area and boundary normal ---*/
      
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }
      
      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there are moving walls (v = u_wall)---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/
      
      node[iPoint]->SetVelocity_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/
      
      Res_Visc[nDim+1] = Wall_HeatFlux * Area;
      
      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      
      if (grid_movement) {
        
        /*--- Get the grid velocity at the current boundary node ---*/
        
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Retrieve other primitive quantities and viscosities ---*/
        
        Density  = node[iPoint]->GetSolution(0);
        Pressure = node[iPoint]->GetPressure();
        laminar_viscosity = node[iPoint]->GetLaminarViscosity();
        eddy_viscosity    = node[iPoint]->GetEddyViscosity();
        total_viscosity   = laminar_viscosity + eddy_viscosity;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
        }
        
        /*--- Divergence of the velocity ---*/
        
        div_vel = 0.0; for (iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim]+Grad_Vel[iDim][jDim] ) - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }
        }
        
        /*--- Dot product of the stress tensor with the grid velocity ---*/
        
        for (iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }
        
        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        
        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Implicit Jacobian contributions due to moving walls ---*/
        
        if (implicit) {
          
          /*--- Jacobian contribution related to the pressure term ---*/
          
          GridVel2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;
          
          /*--- Add the block to the Global Jacobian structure ---*/
          
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
          /*--- Now the Jacobian contribution related to the shear stress ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          /*--- Compute closest normal neighbor ---*/
          
          Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
          
          /*--- Get coordinates of i & nearest normal and compute distance ---*/
          
          Coord_i = geometry->node[iPoint]->GetCoord();
          Coord_j = geometry->node[Point_Normal]->GetCoord();
          
          dist_ij = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
          dist_ij = sqrt(dist_ij);
          
          theta2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            theta2 += UnitNormal[iDim]*UnitNormal[iDim];
          
          factor = total_viscosity*Area/(Density*dist_ij);
          
          if (nDim == 2) {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            
            etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          } else {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;
            
            etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            etax = UnitNormal[1]*UnitNormal[2]/3.0;
            etay = UnitNormal[0]*UnitNormal[2]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }
          
          /*--- Subtract the block from the Global Jacobian structure ---*/
          
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
          
        }
      }
      
      /*--- Convective contribution to the residual at the wall ---*/
      
      LinSysRes.AddBlock(iPoint, Res_Conv);
      
      /*--- Viscous contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      
      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
      
      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
    }
  }
}

void CNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  su2double *Normal, *Coord_i, *Coord_j, Area, dist_ij, theta2;
  su2double Twall, dTdn, dTdrho, thermal_conductivity;
  su2double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  su2double ProjGridVel, *GridVel, GridVel2, Pressure = 0.0, Density, Vel2;
  su2double total_viscosity, div_vel, tau_vel[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};
  su2double laminar_viscosity, eddy_viscosity, Grad_Vel[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}},
  tau[3][3] = {{0.0, 0.0, 0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}, delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Identify the boundary ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Retrieve the specified wall temperature ---*/
  
  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();
  
  /*--- Loop over boundary points ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Compute dual-grid area and boundary normal ---*/
      
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Calculate useful quantities ---*/
      
      theta2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        theta2 += UnitNormal[iDim]*UnitNormal[iDim];
      
      /*--- Compute closest normal neighbor ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Get coordinates of i & nearest normal and compute distance ---*/
      
      Coord_i = geometry->node[iPoint]->GetCoord();
      Coord_j = geometry->node[Point_Normal]->GetCoord();
      dist_ij = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      dist_ij = sqrt(dist_ij);
      
      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there is grid motion (v = u_wall)---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      }
      else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }
      
      /*--- Set the residual, truncation error and velocity value on the boundary ---*/
      
      node[iPoint]->SetVelocity_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
      /*--- Compute the normal gradient in temperature using Twall ---*/
      
      dTdn = -(node[Point_Normal]->GetPrimitive(0) - Twall)/dist_ij;
      
      /*--- Get transport coefficients ---*/
      
      laminar_viscosity    = node[iPoint]->GetLaminarViscosity();
      eddy_viscosity       = node[iPoint]->GetEddyViscosity();
      thermal_conductivity = Cp * ( laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);
      
      // work in progress on real-gases...
      //thermal_conductivity = node[iPoint]->GetThermalConductivity();
      //Cp = node[iPoint]->GetSpecificHeatCp();
      //thermal_conductivity += Cp*eddy_viscosity/Prandtl_Turb;
      
      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/
      
      Res_Visc[nDim+1] = thermal_conductivity * dTdn * Area;
      
      /*--- Calculate Jacobian for implicit time stepping ---*/
      
      if (implicit) {
        
        for (iVar = 0; iVar < nVar; iVar ++)
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_i[iVar][jVar] = 0.0;
        
        /*--- Calculate useful quantities ---*/
        
        Density = node[iPoint]->GetPrimitive(nDim+2);
        Vel2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Vel2 += node[iPoint]->GetPrimitive(iDim+1) * node[iPoint]->GetPrimitive(iDim+1);
        dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );
        
        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
        
        /*--- Add contributions to the Jacobian from the weak enforcement of the energy equations ---*/
        
        Jacobian_i[nDim+1][0]      = -thermal_conductivity*theta2/dist_ij * dTdrho * Area;
        Jacobian_i[nDim+1][nDim+1] = -thermal_conductivity*theta2/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;
        
        /*--- Subtract the block from the Global Jacobian structure ---*/
        
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      
      if (grid_movement) {
        
        /*--- Get the grid velocity at the current boundary node ---*/
        
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Retrieve other primitive quantities and viscosities ---*/
        
        Density  = node[iPoint]->GetSolution(0);
        Pressure = node[iPoint]->GetPressure();
        laminar_viscosity = node[iPoint]->GetLaminarViscosity();
        eddy_viscosity    = node[iPoint]->GetEddyViscosity();
        
        total_viscosity   = laminar_viscosity + eddy_viscosity;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
        }
        
        /*--- Divergence of the velocity ---*/
        
        div_vel = 0.0; for (iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim] ) - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }
        
        /*--- Dot product of the stress tensor with the grid velocity ---*/
        
        for (iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }
        
        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        
        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Implicit Jacobian contributions due to moving walls ---*/
        
        if (implicit) {
          
          /*--- Jacobian contribution related to the pressure term ---*/
          
          GridVel2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;
          
          /*--- Add the block to the Global Jacobian structure ---*/
          
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
          /*--- Now the Jacobian contribution related to the shear stress ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          factor = total_viscosity*Area/(Density*dist_ij);
          
          if (nDim == 2) {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            
            etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          }
          else {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;
            
            etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            etax = UnitNormal[1]*UnitNormal[2]/3.0;
            etay = UnitNormal[0]*UnitNormal[2]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }
          
          /*--- Subtract the block from the Global Jacobian structure ---*/
          
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
        
      }
      
      /*--- Convective contribution to the residual at the wall ---*/
      
      LinSysRes.AddBlock(iPoint, Res_Conv);
      
      /*--- Viscous contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      
      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
      
      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
    }
  }
}
