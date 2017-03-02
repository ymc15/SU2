/*!
 * \file solver_direct_elasticity.cpp
 * \brief Main subroutines for solving direct FEM elasticity problems.
 * \author R. Sanchez
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

CFEM_ElasticitySolver::CFEM_ElasticitySolver(void) : CSolver() {
  
  nElement = 0;
  nDim = 0;
  nMarker = 0;
  
  nPoint = 0;
  nPointDomain = 0;
  
  Total_CFEA = 0.0;
  WAitken_Dyn = 0.0;
  WAitken_Dyn_tn1 = 0.0;
  loadIncrement = 1.0;
  
  element_container = NULL;
  node = NULL;
  
  GradN_X = NULL;
  GradN_x = NULL;
  
  Jacobian_c_ij = NULL;
  Jacobian_s_ij = NULL;
  Jacobian_k_ij = NULL;
  
  MassMatrix_ij = NULL;
  
  mZeros_Aux = NULL;
  mId_Aux = NULL;
  
  Res_Stress_i = NULL;
  Res_Ext_Surf = NULL;
  Res_Time_Cont = NULL;
  Res_FSI_Cont = NULL;
  
  Res_Dead_Load = NULL;
  
  nodeReactions = NULL;
  
  solutionPredictor = NULL;
  
  SolRest = NULL;
  
  normalVertex = NULL;
  stressTensor = NULL;
  
  Solution_Interm = NULL;
  
}

CFEM_ElasticitySolver::CFEM_ElasticitySolver(CGeometry *geometry, CConfig *config) : CSolver() {
  
  unsigned long iPoint;
  unsigned short iVar, jVar, iDim, jDim;
  unsigned short iTerm, iKind;
  
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation
  bool gen_alpha = (config->GetKind_TimeIntScheme_FEA() == GENERALIZED_ALPHA);  // Generalized alpha method requires residual at previous time step.
  
  bool body_forces = config->GetDeadLoad();  // Body forces (dead loads).
  bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  su2double E = config->GetElasticyMod();
  
  nElement      = geometry->GetnElem();
  nDim          = geometry->GetnDim();
  nMarker       = geometry->GetnMarker();
  
  nPoint        = geometry->GetnPoint();
  nPointDomain  = geometry->GetnPointDomain();
  
  /*--- Here is where we assign the kind of each element ---*/
  
  /*--- First level: different possible terms of the equations ---*/
  element_container = new CElement** [MAX_TERMS];
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS];
  
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++) {
    for (iKind = 0; iKind < MAX_FE_KINDS; iKind++) {
      element_container[iTerm][iKind] = NULL;
    }
  }
  
  if (nDim == 2) {
    if (incompressible) {
      element_container[FEA_TERM][EL_TRIA] = new CTRIA1(nDim, config);
      element_container[FEA_TERM][EL_QUAD] = new CQUAD4P1(nDim, config);
    }
    else {
      element_container[FEA_TERM][EL_TRIA] = new CTRIA1(nDim, config);
      element_container[FEA_TERM][EL_QUAD] = new CQUAD4(nDim, config);
    }
  }
  else if (nDim == 3) {
    if (incompressible) {
      element_container[FEA_TERM][EL_TETRA] = new CTETRA1(nDim, config);
      element_container[FEA_TERM][EL_HEXA] = new CHEXA8P1(nDim, config);
    }
    else {
      element_container[FEA_TERM][EL_TETRA] = new CTETRA1(nDim, config);
      element_container[FEA_TERM][EL_HEXA] = new CHEXA8(nDim, config);
    }
  }
  
  node              = new CVariable*[nPoint];
  
  GradN_X = new su2double [nDim];
  GradN_x = new su2double [nDim];
  
  Total_CFEA      = 0.0;
  WAitken_Dyn        = 0.0;
  WAitken_Dyn_tn1    = 0.0;
  loadIncrement     = 0.0;
  
  SetFSI_ConvValue(0,0.0);
  SetFSI_ConvValue(1,0.0);
  
  nVar = nDim;
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual = new su2double[nVar];          for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
  
  Solution_Interm = NULL;
  if (gen_alpha) {
    Solution_Interm = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Solution_Interm[iVar] = 0.0;
  }
  
  nodeReactions = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) nodeReactions[iVar]   = 0.0;
  
  /*--- The length of the solution vector depends on whether the problem is static or dynamic ---*/
  
  unsigned short nSolVar;
  unsigned long index;
  string text_line, filename;
  ifstream restart_file;
  su2double dull_val;
  long Dyn_RestartIter;
  
  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;
  
  SolRest = new su2double[nSolVar];
  
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  
  /*--- Check for a restart, initialize from zero otherwise ---*/
  
  if (!restart) {
    for (iVar = 0; iVar < nSolVar; iVar++) SolRest[iVar] = 0.0;
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CFEM_ElasVariable(SolRest, nDim, nVar, config);
    }
  }
  else {
    
    /*--- Restart the solution from file information ---*/
    
    filename = config->GetSolution_FEMFileName();
    
    /*--- If multizone, append zone name ---*/
    if (nZone > 1)
      filename = config->GetMultizone_FileName(filename, iZone);
    
    if (dynamic) {
      
      Dyn_RestartIter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
      
      filename = config->GetUnsteady_FileName(filename, (int)Dyn_RestartIter);
    }
    
    restart_file.open(filename.data(), ios::in);
    
    /*--- In case there is no file ---*/
    
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no FEM restart file!!" << endl;
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
      
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/
      
      MI = Global2Local.find(iPoint_Global);
      if (MI != Global2Local.end()) {
        
        iPoint_Local = Global2Local[iPoint_Global];
        
        if (dynamic) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> SolRest[0] >> SolRest[1] >> SolRest[2] >> SolRest[3] >> SolRest[4] >> SolRest[5];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> SolRest[0] >> SolRest[1] >> SolRest[2] >> SolRest[3] >> SolRest[4] >> SolRest[5] >> SolRest[6] >> SolRest[7] >> SolRest[8];
        }
        else {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> SolRest[0] >> SolRest[1];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> SolRest[0] >> SolRest[1] >> SolRest[2];
        }
        
        node[iPoint_Local] = new CFEM_ElasVariable(SolRest, nDim, nVar, config);
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
     because a send/recv is performed immediately in the solver (Set_MPI_Solution()). ---*/
    
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CFEM_ElasVariable(SolRest, nDim, nVar, config);
    }
    
    
    /*--- Close the restart file ---*/
    
    restart_file.close();
    
  }
  
  
  bool prestretch_fem = config->GetPrestretch();
  if (prestretch_fem) Set_Prestretch(geometry, config);
  
  
  /*--- Term ij of the Jacobian ---*/
  
  Jacobian_ij = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_ij[iVar] = new su2double [nVar];
    for (jVar = 0; jVar < nVar; jVar++) {
      Jacobian_ij[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Term ij of the Mass Matrix (only if dynamic analysis) ---*/
  MassMatrix_ij = NULL;
  if (dynamic) {
    MassMatrix_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      MassMatrix_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        MassMatrix_ij[iVar][jVar] = 0.0;
      }
    }
  }
  
  Jacobian_c_ij = NULL;
  Jacobian_s_ij = NULL;
  if (nonlinear_analysis) {
    
    /*--- Term ij of the Jacobian (constitutive contribution) ---*/
    
    Jacobian_c_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_c_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_c_ij[iVar][jVar] = 0.0;
      }
    }
    
    /*--- Term ij of the Jacobian (stress contribution) ---*/
    
    Jacobian_s_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_s_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_s_ij[iVar][jVar] = 0.0;
      }
    }
    
  }
  
  /*--- Term ij of the Jacobian (incompressibility term) ---*/
  Jacobian_k_ij = NULL;
  if (incompressible) {
    Jacobian_k_ij = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_k_ij[iVar] = new su2double [nVar];
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_k_ij[iVar][jVar] = 0.0;
      }
    }
  }
  
  /*--- Stress contribution to the node i ---*/
  Res_Stress_i = new su2double[nVar];
  
  /*--- Contribution of the external surface forces to the residual (auxiliary vector) ---*/
  Res_Ext_Surf = new su2double[nVar];
  
  /*--- Contribution of the body forces to the residual (auxiliary vector) ---*/
  Res_Dead_Load = NULL;
  if (body_forces) {
    Res_Dead_Load = new su2double[nVar];
  }
  
  /*--- Contribution of the fluid tractions to the residual (auxiliary vector) ---*/
  Res_FSI_Cont = NULL;
  if (fsi) {
    Res_FSI_Cont = new su2double[nVar];
  }
  
  /*--- Time integration contribution to the residual ---*/
  Res_Time_Cont = NULL;
  if (dynamic) {
    Res_Time_Cont = new su2double [nVar];
  }
  
  /*--- Matrices to impose clamped boundary conditions ---*/
  
  mZeros_Aux = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++)
    mZeros_Aux[iDim] = new su2double[nDim];
  
  mId_Aux = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++)
    mId_Aux[iDim] = new su2double[nDim];
  
  for(iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      mZeros_Aux[iDim][jDim] = 0.0;
      mId_Aux[iDim][jDim] = 0.0;
    }
    mId_Aux[iDim][iDim] = E;
  }
  
  
  /*--- Initialization of matrix structures ---*/
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Non-Linear Elasticity)." << endl;
  
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  
  if (dynamic) {
    MassMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
    TimeRes_Aux.Initialize(nPoint, nPointDomain, nVar, 0.0);
    TimeRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }
  
  /*--- Initialization of linear solver structures ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  LinSysReact.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Initialize the auxiliary vector and matrix for the computation of the nodal Reactions ---*/
  
  normalVertex = new su2double [nDim];
  
  stressTensor = new su2double* [nDim];
  for (iVar = 0; iVar < nVar; iVar++) {
    stressTensor[iVar] = new su2double [nDim];
  }
  
  /*---- Initialize the auxiliary vector for the solution predictor ---*/
  
  solutionPredictor = new su2double [nVar];
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- If dynamic, we also need to communicate the old solution ---*/
  
  if(dynamic) Set_MPI_Solution_Old(geometry, config);
  
}

CFEM_ElasticitySolver::~CFEM_ElasticitySolver(void) {
  
  unsigned short iVar, jVar;
  
  if (element_container != NULL) {
    for (iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (jVar = 0; jVar < MAX_FE_KINDS; jVar++) {
        if (element_container[iVar][jVar] != NULL) delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }
  
  for (iVar = 0; iVar < nVar; iVar++) {
    if (Jacobian_s_ij != NULL) delete [] Jacobian_s_ij[iVar];
    if (Jacobian_c_ij != NULL) delete [] Jacobian_c_ij[iVar];
    if (Jacobian_k_ij != NULL) delete[] Jacobian_k_ij[iVar];
    delete [] mZeros_Aux[iVar];
    delete [] mId_Aux[iVar];
    delete [] stressTensor[iVar];
  }
  
  if (Jacobian_s_ij != NULL) delete [] Jacobian_s_ij;
  if (Jacobian_c_ij != NULL) delete [] Jacobian_c_ij;
  if (Jacobian_k_ij != NULL) delete [] Jacobian_k_ij;
  delete [] Res_Stress_i;
  delete [] Res_Ext_Surf;
  if (Res_Time_Cont != NULL) delete[] Res_Time_Cont;
  if (Res_Dead_Load != NULL) delete[] Res_Dead_Load;
  delete [] SolRest;
  delete [] GradN_X;
  delete [] GradN_x;
  
  delete [] mZeros_Aux;
  delete [] mId_Aux;
  
  delete [] nodeReactions;
  
  delete [] normalVertex;
  delete [] stressTensor;
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  
  unsigned short nSolVar;
  
  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;
  
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
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
        if (dynamic) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel(iVar);
            Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel(iVar);
          }
        }
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
        if (dynamic) {
          for (iVar = nVar; iVar < 3*nVar; iVar++)
            Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        }
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nSolVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, SolRest[iVar]);
        
        if (dynamic) {
          
          for (iVar = 0; iVar < nVar; iVar++) {
            node[iPoint]->SetSolution_Vel(iVar, SolRest[iVar+nVar]);
            node[iPoint]->SetSolution_Accel(iVar, SolRest[iVar+2*nVar]);
          }
          
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  unsigned short nSolVar;
  
  nSolVar = 3 * nVar;
  
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
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++) {
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_time_n(iVar);
          Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel_time_n(iVar);
          Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel_time_n(iVar);
        }
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nSolVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nSolVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          node[iPoint]->SetSolution_time_n(iVar, SolRest[iVar]);
          node[iPoint]->SetSolution_Vel_time_n(iVar, SolRest[iVar+nVar]);
          node[iPoint]->SetSolution_Accel_time_n(iVar, SolRest[iVar+2*nVar]);
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution_DispOnly(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
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
      nBufferS_Vector = nVertexS*nVar;         nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
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
        
        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, SolRest[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution_Pred(CGeometry *geometry, CConfig *config) {
  
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
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
      nBufferS_Vector = nVertexS*nVar;     nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred(iVar);
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
        
        /*--- Copy predicted solution variables back into the variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Pred(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Set_MPI_Solution_Pred_Old(CGeometry *geometry, CConfig *config) {
  
  /*--- We are communicating the solution predicted, current and old, and the old solution ---*/
  /*--- necessary for the Aitken relaxation ---*/
  
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
  /*--- Analogous to the dynamic solution, in this case we need 3 * nVar variables per node ---*/
  unsigned short nSolVar;
  nSolVar = 3 * nVar;
  
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
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++) {
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
          Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred(iVar);
          Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred_Old(iVar);
        }
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nSolVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
          node[iPoint]->SetSolution_Pred(iVar, Buffer_Receive_U[(iVar+nVar)*nVertexR+iVertex]);
          node[iPoint]->SetSolution_Pred_Old(iVar, Buffer_Receive_U[(iVar+2*nVar)*nVertexR+iVertex]);
        }
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}


void CFEM_ElasticitySolver::Set_Prestretch(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint;
  unsigned long index;
  
  unsigned short iVar;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  
  string filename;
  ifstream prestretch_file;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  
  /*--- Restart the solution from file information ---*/
  
  filename = config->GetPrestretch_FEMFileName();
  
  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);
  
  cout << "Filename: " << filename << "." << endl;
  
  prestretch_file.open(filename.data(), ios::in);
  
  /*--- In case there is no file ---*/
  
  if (prestretch_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no FEM prestretch reference file!!" << endl;
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
  
  getline (prestretch_file, text_line);
  
  while (getline (prestretch_file, text_line)) {
    istringstream point_line(text_line);
    
    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/
    
    MI = Global2Local.find(iPoint_Global);
    if (MI != Global2Local.end()) {
      
      iPoint_Local = Global2Local[iPoint_Global];
      
      if (nDim == 2) point_line >> Solution[0] >> Solution[1] >> index;
      if (nDim == 3) point_line >> Solution[0] >> Solution[1] >> Solution[2] >> index;
      
      for (iVar = 0; iVar < nVar; iVar++) node[iPoint_Local]->SetPrestretch(iVar, Solution[iVar]);
      
      iPoint_Global_Local++;
    }
    iPoint_Global++;
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
  
  /*--- TODO: We need to communicate here the prestretched geometry for the halo nodes. ---*/
  
  /*--- Close the restart file ---*/
  
  prestretch_file.close();
  
}


void CFEM_ElasticitySolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {
  
  
  unsigned long iPoint;
  bool initial_calc = (config->GetExtIter() == 0);                  // Checks if it is the first calculation.
  bool first_iter = (config->GetIntIter() == 0);                          // Checks if it is the first iteration
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool restart = config->GetRestart();                        // Restart analysis
  bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter()); // Initial calculation for restart
  
  bool incremental_load = config->GetIncrementalLoad();                // If an incremental load is applied
  
  bool body_forces = config->GetDeadLoad();                      // Body forces (dead loads).
  
  /*--- Set vector entries to zero ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    LinSysAux.SetBlock_Zero(iPoint);
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }
  
  /*--- Set matrix entries to zero ---*/
  
  /*
   * If the problem is linear, we only need one Jacobian matrix in the problem, because
   * it is going to be constant along the calculations. Therefore, we only initialize
   * the Jacobian matrix once, at the beginning of the simulation.
   *
   * We don't need first_iter, because there is only one iteration per time step in linear analysis.
   */
  if ((initial_calc && linear_analysis)||
      (restart && initial_calc_restart && linear_analysis)) {
    Jacobian.SetValZero();
  }
  
  /*
   * If the problem is dynamic, we need a mass matrix, which will be constant along the calculation
   * both for linear and nonlinear analysis. Only initialized once, at the first time step.
   *
   * The same with the integration constants, as for now we consider the time step to be constant.
   *
   * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
   */
  if ((dynamic && initial_calc && first_iter) ||
      (dynamic && restart && initial_calc_restart && first_iter)) {
    MassMatrix.SetValZero();
    Compute_IntegrationConstants(config);
    Compute_MassMatrix(geometry, solver_container, numerics, config);
  }
  
  /*
   * If body forces are taken into account, we need to compute the term that goes into the residual,
   * which will be constant along the calculation both for linear and nonlinear analysis.
   *
   * Only initialized once, at the first iteration or the beginning of the calculation after a restart.
   *
   * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
   */
  
  if ((body_forces && initial_calc && first_iter) ||
      (body_forces && restart && initial_calc_restart && first_iter)) {
    // If the load is incremental, we have to reset the variable to avoid adding up over the increments
    if (incremental_load) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Clear_BodyForces_Res();
    }
    // Compute the dead load term
    Compute_DeadLoad(geometry, solver_container, numerics, config);
  }
  
  /*
   * If the problem is nonlinear, we need to initialize the Jacobian and the stiffness matrix at least at the beginning
   * of each time step. If the solution method is Newton Rapshon, we initialize it also at the beginning of each
   * iteration.
   */
  
  if ((nonlinear_analysis) && ((newton_raphson) || (first_iter)))  {
    Jacobian.SetValZero();
    //    StiffMatrix.SetValZero();
  }
  
  /*
   * Some external forces may be considered constant over the time step.
   */
  if (first_iter)  {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Clear_SurfaceLoad_Res();
  }
  
  /*
   * If we apply pressure forces, we need to clear the residual on each iteration
   */
  unsigned short iMarker;
  unsigned long iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case LOAD_BOUNDARY:
        /*--- Only if the load is nonzero - reduces computational cost ---*/
        if(config->GetLoad_Value(config->GetMarker_All_TagBound(iMarker)) != 0 ) {
          /*--- For all the vertices in the marker iMarker ---*/
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            /*--- Retrieve the point ID ---*/
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            /*--- Clear the residual of the node, to avoid adding on previous values ---*/
            node[iPoint]->Clear_SurfaceLoad_Res();
          }
        }
        break;
    }
  
  
}

void CFEM_ElasticitySolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CFEM_ElasticitySolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, nPoint;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied
  
  nPoint = geometry[MESH_0]->GetnPoint();
  
  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/
  
  if (incremental_load) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_OldSolution();
  }
  
  
}

void CFEM_ElasticitySolver::ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, nPoint;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied
  
  nPoint = geometry[MESH_0]->GetnPoint();
  
  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/
  
  if (incremental_load) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_Solution();
  }
  
}

void CFEM_ElasticitySolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;
  
  su2double *Kab = NULL;
  unsigned short NelNodes, jNode;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }
    
    numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
          }
        }
        
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);
        
      }
      
    }
    
  }
  
  
}

void CFEM_ElasticitySolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  su2double Ks_ab;
  su2double *Kab = NULL;
  su2double *Kk_ab = NULL;
  su2double *Ta = NULL;
  unsigned short NelNodes, jNode;
  
  bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    /*--- If incompressible, we compute the Mean Dilatation term first so the volume is already computed ---*/
    
    if (incompressible) numerics[FEA_TERM]->Compute_MeanDilatation_Term(element_container[FEA_TERM][EL_KIND], config);
    
    numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];
      
      /*--- Check if this is my node or not ---*/
      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);
        Ks_ab = element_container[FEA_TERM][EL_KIND]->Get_Ks_ab(iNode,jNode);
        if (incompressible) Kk_ab = element_container[FEA_TERM][EL_KIND]->Get_Kk_ab(iNode,jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          Jacobian_s_ij[iVar][iVar] = Ks_ab;
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_c_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
            if (incompressible) Jacobian_k_ij[iVar][jVar] = Kk_ab[iVar*nVar+jVar];
          }
        }
        
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
        if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);
        
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;
  
  su2double Mab;
  unsigned short NelNodes, jNode;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    {nNodes = 4; EL_KIND = EL_QUAD;}
    
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }
    
    numerics[FEA_TERM]->Compute_Mass_Matrix(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      for (jNode = 0; jNode < NelNodes; jNode++) {
        
        Mab = element_container[FEA_TERM][EL_KIND]->Get_Mab(iNode, jNode);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          MassMatrix_ij[iVar][iVar] = Mab;
        }
        
        MassMatrix.AddBlock(indexNode[iNode], indexNode[jNode], MassMatrix_ij);
        
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  su2double *Ta = NULL;
  unsigned short NelNodes;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      //      for (iDim = 0; iDim < nDim; iDim++) {
      //        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
      //        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
      //        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      //        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
      //      }
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    numerics[FEA_TERM]->Compute_NodalStress_Term(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];
      
      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
      
    }
    
  }
  
  for (iDim = 0; iDim < nDim; iDim++) {
    val_Coord = geometry->node[0]->GetCoord(iDim);
    val_Sol = node[0]->GetSolution(iDim) + val_Coord;
  }
  
}

void CFEM_ElasticitySolver::Compute_NodalStress(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iPoint, iElem, iVar;
  unsigned short iNode, iDim, iStress;
  unsigned short nNodes = 0, nStress;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol, val_Ref = 0.0;
  int EL_KIND = 0;
  
  bool prestretch_fem = config->GetPrestretch();
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  if (nDim == 2) nStress = 3;
  else nStress = 6;
  
  su2double *Ta = NULL;
  
  unsigned short NelNodes;
  
  /*--- Restart stress to avoid adding results from previous time steps ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iStress = 0; iStress < nStress; iStress++) {
      node[iPoint]->SetStress_FEM(iStress, 0.0);
    }
  }
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      //      for (iDim = 0; iDim < nDim; iDim++) {
      //        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
      //        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
      //        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      //        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
      //      }
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
        if (prestretch_fem) {
          val_Ref = node[indexNode[iNode]]->GetPrestretch(iDim);
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Ref, iNode, iDim);
        }
        else {
          element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        }
      }
    }
    
    numerics[FEA_TERM]->Compute_Averaged_NodalStress(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      /*--- This only works if the problem is nonlinear ---*/
      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];
      
      LinSysReact.AddBlock(indexNode[iNode], Res_Stress_i);
      
      for (iStress = 0; iStress < nStress; iStress++) {
        node[indexNode[iNode]]->AddStress_FEM(iStress,
                                              (element_container[FEA_TERM][EL_KIND]->Get_NodalStress(iNode, iStress) /
                                               geometry->node[indexNode[iNode]]->GetnElem()) );
      }
      
    }
    
  }
  
  su2double *Stress;
  su2double VonMises_Stress, MaxVonMises_Stress = 0.0;
  su2double Sxx,Syy,Szz,Sxy,Sxz,Syz,S1,S2;
  
  /*--- For the number of nodes in the mesh ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Get the stresses, added up from all the elements that connect to the node ---*/
    
    Stress  = node[iPoint]->GetStress_FEM();
    
    /*--- Compute the stress averaged from all the elements connecting to the node and the Von Mises stress ---*/
    
    if (nDim == 2) {
      
      Sxx=Stress[0];
      Syy=Stress[1];
      Sxy=Stress[2];
      
      S1=(Sxx+Syy)/2+sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
      S2=(Sxx+Syy)/2-sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
      
      VonMises_Stress = sqrt(S1*S1+S2*S2-2*S1*S2);
      
    }
    else {
      
      Sxx = Stress[0];
      Syy = Stress[1];
      Szz = Stress[3];
      
      Sxy = Stress[2];
      Sxz = Stress[4];
      Syz = Stress[5];
      
      VonMises_Stress = sqrt(0.5*(   pow(Sxx - Syy, 2.0)
                                  + pow(Syy - Szz, 2.0)
                                  + pow(Szz - Sxx, 2.0)
                                  + 6.0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz)
                                  ));
      
    }
    
    node[iPoint]->SetVonMises_Stress(VonMises_Stress);
    
    /*--- Compute the maximum value of the Von Mises Stress ---*/
    
    MaxVonMises_Stress = max(MaxVonMises_Stress, VonMises_Stress);
    
  }
  
#ifdef HAVE_MPI
  
  /*--- Compute MaxVonMises_Stress using all the nodes ---*/
  
  su2double MyMaxVonMises_Stress = MaxVonMises_Stress; MaxVonMises_Stress = 0.0;
  SU2_MPI::Allreduce(&MyMaxVonMises_Stress, &MaxVonMises_Stress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
#endif
  
  /*--- Set the value of the MaxVonMises_Stress as the CFEA coeffient ---*/
  
  Total_CFEA = MaxVonMises_Stress;
  
  
  bool outputReactions = false;
  
  if (outputReactions) {
    
    ofstream myfile;
    myfile.open ("Reactions.txt");
    
    unsigned short iMarker;
    unsigned long iVertex;
    su2double val_Reaction;
    
    bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
    
    if (!dynamic) {
      /*--- Loop over all the markers  ---*/
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        switch (config->GetMarker_All_KindBC(iMarker)) {
            
            /*--- If it corresponds to a clamped boundary  ---*/
            
          case CLAMPED_BOUNDARY:
            
            myfile << "MARKER " << iMarker << ":" << endl;
            
            /*--- Loop over all the vertices  ---*/
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              
              /*--- Get node index ---*/
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              
              myfile << "Node " << iPoint << "." << " \t ";
              
              for (iDim = 0; iDim < nDim; iDim++) {
                /*--- Retrieve coordinate ---*/
                val_Coord = geometry->node[iPoint]->GetCoord(iDim);
                myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
              }
              
              for (iVar = 0; iVar < nVar; iVar++) {
                /*--- Retrieve reaction ---*/
                val_Reaction = LinSysReact.GetBlock(iPoint, iVar);
                myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
              }
              
              myfile << endl;
            }
            myfile << endl;
            break;
        }
    }
    else if (dynamic) {
      
      switch (config->GetKind_TimeIntScheme_FEA()) {
        case (CD_EXPLICIT):
          cout << "NOT IMPLEMENTED YET" << endl;
          break;
        case (NEWMARK_IMPLICIT):
          
          /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
          if (linear_analysis) {
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
              for (iVar = 0; iVar < nVar; iVar++) {
                Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+    //a0*U(t)
                a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+  //a2*U'(t)
                a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
              }
              TimeRes_Aux.SetBlock(iPoint, Residual);
            }
          }
          else if (nonlinear_analysis) {
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
              for (iVar = 0; iVar < nVar; iVar++) {
                Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
                - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
                + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
                + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
              }
              TimeRes_Aux.SetBlock(iPoint, Residual);
            }
          }
          /*--- Once computed, compute M*TimeRes_Aux ---*/
          MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
          
          /*--- Loop over all the markers  ---*/
          for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            switch (config->GetMarker_All_KindBC(iMarker)) {
                
                /*--- If it corresponds to a clamped boundary  ---*/
                
              case CLAMPED_BOUNDARY:
                
                myfile << "MARKER " << iMarker << ":" << endl;
                
                /*--- Loop over all the vertices  ---*/
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                  
                  /*--- Get node index ---*/
                  iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                  
                  myfile << "Node " << iPoint << "." << " \t ";
                  
                  for (iDim = 0; iDim < nDim; iDim++) {
                    /*--- Retrieve coordinate ---*/
                    val_Coord = geometry->node[iPoint]->GetCoord(iDim);
                    myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
                  }
                  
                  for (iVar = 0; iVar < nVar; iVar++) {
                    /*--- Retrieve the time contribution ---*/
                    Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
                    /*--- Retrieve reaction ---*/
                    val_Reaction = LinSysReact.GetBlock(iPoint, iVar) + Res_Time_Cont[iVar];
                    myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
                  }
                  
                  myfile << endl;
                }
                myfile << endl;
                break;
            }
          
          
          break;
        case (GENERALIZED_ALPHA):
          cout << "NOT IMPLEMENTED YET" << endl;
          break;
      }
      
    }
    
    
    
    myfile.close();
    
  }
  
}

void CFEM_ElasticitySolver::Compute_DeadLoad(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {
  
  unsigned long iElem, iVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;
  
  su2double *Dead_Load = NULL;
  unsigned short NelNodes;
  
  /*--- Loops over all the elements ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
    
    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }
    
    numerics[FEA_TERM]->Compute_Dead_Load(element_container[FEA_TERM][EL_KIND], config);
    
    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();
    
    for (iNode = 0; iNode < NelNodes; iNode++) {
      
      Dead_Load = element_container[FEA_TERM][EL_KIND]->Get_FDL_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Dead_Load[iVar] = Dead_Load[iVar];
      
      node[indexNode[iNode]]->Add_BodyForces_Res(Res_Dead_Load);
      
    }
    
  }
  
  
}

void CFEM_ElasticitySolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
}

void CFEM_ElasticitySolver::Compute_IntegrationConstants(CConfig *config) {
  
  su2double Delta_t= config->GetDelta_DynTime();
  
  su2double delta = config->GetNewmark_delta(), alpha = config->GetNewmark_alpha();
  
  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      cout << "NOT IMPLEMENTED YET" << endl;
      break;
    case (NEWMARK_IMPLICIT):
      
      /*--- Integration constants for Newmark scheme ---*/
      
      a_dt[0]= 1 / (alpha*pow(Delta_t,2.0));
      a_dt[1]= delta / (alpha*Delta_t);
      a_dt[2]= 1 / (alpha*Delta_t);
      a_dt[3]= 1 /(2*alpha) - 1;
      a_dt[4]= delta/alpha - 1;
      a_dt[5]= (Delta_t/2) * (delta/alpha - 2);
      a_dt[6]= Delta_t * (1-delta);
      a_dt[7]= delta * Delta_t;
      a_dt[8]= 0.0;
      
      break;
      
    case (GENERALIZED_ALPHA):
      
      /*--- Integration constants for Generalized Alpha ---*/
      /*--- Needs to be updated if accounting for structural damping ---*/
      
      //      su2double beta = config->Get_Int_Coeffs(0);
      //      //  su2double gamma =  config->Get_Int_Coeffs(1);
      //      su2double alpha_f = config->Get_Int_Coeffs(2), alpha_m =  config->Get_Int_Coeffs(3);
      //
      //      a_dt[0]= (1 / (beta*pow(Delta_t,2.0))) * ((1 - alpha_m) / (1 - alpha_f)) ;
      //      a_dt[1]= 0.0 ;
      //      a_dt[2]= (1 - alpha_m) / (beta*Delta_t);
      //      a_dt[3]= ((1 - 2*beta)*(1-alpha_m) / (2*beta)) - alpha_m;
      //      a_dt[4]= 0.0;
      //      a_dt[5]= 0.0;
      //      a_dt[6]= Delta_t * (1-delta);
      //      a_dt[7]= delta * Delta_t;
      //      a_dt[8]= (1 - alpha_m) / (beta*pow(Delta_t,2.0));
      
      break;
  }
  
  
}


void CFEM_ElasticitySolver::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                       unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, jVar;
  
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Get node index ---*/
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      if (nDim == 2) {
        Solution[0] = 0.0;  Solution[1] = 0.0;
        Residual[0] = 0.0;  Residual[1] = 0.0;
      }
      else {
        Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
        Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
      }
      
      node[iPoint]->SetSolution(Solution);
      
      if (dynamic) {
        node[iPoint]->SetSolution_Vel(Solution);
        node[iPoint]->SetSolution_Accel(Solution);
      }
      
      
      /*--- Initialize the reaction vector ---*/
      LinSysReact.SetBlock(iPoint, Residual);
      
      LinSysRes.SetBlock(iPoint, Residual);
      
      /*--- STRONG ENFORCEMENT OF THE DISPLACEMENT BOUNDARY CONDITION ---*/
      
      /*--- Delete the columns for a particular node ---*/
      
      for (iVar = 0; iVar < nPoint; iVar++) {
        if (iVar==iPoint) {
          Jacobian.SetBlock(iVar,iPoint,mId_Aux);
        }
        else {
          Jacobian.SetBlock(iVar,iPoint,mZeros_Aux);
        }
      }
      
      /*--- Delete the rows for a particular node ---*/
      for (jVar = 0; jVar < nPoint; jVar++) {
        if (iPoint!=jVar) {
          Jacobian.SetBlock(iPoint,jVar,mZeros_Aux);
        }
      }
      
      /*--- If the problem is dynamic ---*/
      /*--- Enforce that in the previous time step all nodes had 0 U, U', U'' ---*/
      
      if(dynamic) {
        
        node[iPoint]->SetSolution_time_n(Solution);
        node[iPoint]->SetSolution_Vel_time_n(Solution);
        node[iPoint]->SetSolution_Accel_time_n(Solution);
        
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                            unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Get node index ---*/
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (nDim == 2) {
      Solution[0] = 0.0;  Solution[1] = 0.0;
    }
    else {
      Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
    }
    
    node[iPoint]->SetSolution(Solution);
    
    if (dynamic) {
      node[iPoint]->SetSolution_Vel(Solution);
      node[iPoint]->SetSolution_Accel(Solution);
    }
    
  }
  
}

void CFEM_ElasticitySolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
                                           unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  bool first_iter = (config->GetIntIter() == 0);
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);    // Nonlinear analysis.
  
  su2double solNorm = 0.0, solNorm_recv = 0.0;
  
#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (nonlinear_analysis) {
    
    /*--- If the problem is nonlinear, we have 3 convergence criteria ---*/
    
    /*--- UTOL = norm(Delta_U(k)) / norm(U(k)) --------------------------*/
    /*--- RTOL = norm(Residual(k)) / norm(Residual(0)) ------------------*/
    /*--- ETOL = Delta_U(k) * Residual(k) / Delta_U(0) * Residual(0) ----*/
    
    if (first_iter) {
      Conv_Ref[0] = 1.0;                      // Position for the norm of the solution
      Conv_Ref[1] = max(LinSysRes.norm(), EPS);          // Position for the norm of the residual
      Conv_Ref[2] = max(dotProd(LinSysSol, LinSysRes), EPS);    // Position for the energy tolerance
      
      /*--- Make sure the computation runs at least 2 iterations ---*/
      Conv_Check[0] = 1.0;
      Conv_Check[1] = 1.0;
      Conv_Check[2] = 1.0;
    }
    else {
      /*--- Compute the norm of the solution vector Uk ---*/
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          solNorm += node[iPoint]->GetSolution(iVar) * node[iPoint]->GetSolution(iVar);
        }
      }
      
      // We need to communicate the norm of the solution and compute the RMS throughout the different processors
      
#ifdef HAVE_MPI
      /*--- We sum the squares of the norms across the different processors ---*/
      SU2_MPI::Allreduce(&solNorm, &solNorm_recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      solNorm_recv         = solNorm;
#endif
      
      Conv_Ref[0] = max(sqrt(solNorm_recv), EPS);            // Norm of the solution vector
      
      Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];          // Norm of the delta-solution vector
      Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];          // Norm of the residual
      Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
      
    }
    
    /*--- MPI solution ---*/
    
    Set_MPI_Solution(geometry, config);
    
  } else {
    
    /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/
    
    /*---  Compute the residual Ax-f ---*/
    
    Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
    
    /*--- Set maximum residual to zero ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      SetRes_RMS(iVar, 0.0);
      SetRes_Max(iVar, 0.0, 0);
    }
    
    /*--- Compute the residual ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
        AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
    
    
    /*--- MPI solution ---*/
    
    Set_MPI_Solution(geometry, config);
    
    /*--- Compute the root mean square residual ---*/
    
    SetResidual_RMS(geometry, config);
  }
  
}

void CFEM_ElasticitySolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                                   unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                           unsigned short val_marker) {
  
  /*--- Retrieve the normal pressure and the application conditions for the considered boundary ---*/
  
  su2double NormalLoad = config->GetLoad_Value(config->GetMarker_All_TagBound(val_marker));
  su2double TotalLoad = 0.0;
  bool Sigmoid_Load = config->GetSigmoid_Load();
  su2double Sigmoid_Time = config->GetSigmoid_Time();
  su2double Sigmoid_K = config->GetSigmoid_K();
  su2double SigAux = 0.0;
  
  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl, NonModAmpl;
  
  bool Ramp_Load = config->GetRamp_Load();
  su2double Ramp_Time = config->GetRamp_Time();
  
  if (Ramp_Load) {
    ModAmpl = NormalLoad*CurrentTime/Ramp_Time;
    NonModAmpl = NormalLoad;
    TotalLoad=min(ModAmpl,NonModAmpl);
  }
  else if (Sigmoid_Load) {
    SigAux = CurrentTime/ Sigmoid_Time;
    ModAmpl = (1 / (1+exp(-1*Sigmoid_K*(SigAux - 0.5)) ) );
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
    TotalLoad=ModAmpl*NormalLoad;
  }
  else {
    TotalLoad = NormalLoad;
  }
  
  /*--- Do only if there is a load applied.
   *--- This reduces the computational cost for cases in which we want boundaries with no load.
   */
  
  if (TotalLoad != 0.0) {
    
    unsigned long iElem;
    unsigned short nNodes = 0;
    su2double Length_Elem_ref = 0.0,  Area_Elem_ref = 0.0;
    su2double Length_Elem_curr = 0.0, Area_Elem_curr = 0.0;
    unsigned long indexNode[4]   = {0,0,0,0};
    unsigned long indexVertex[4] = {0,0,0,0};
    su2double nodeCoord_ref[4][3], nodeCoord_curr[4][3];
    su2double *nodal_normal, nodal_normal_unit[3];
    su2double normal_ref_unit[3], normal_curr_unit[3];
    su2double Norm, dot_Prod;
    su2double val_Coord, val_Sol;
    unsigned short iNode, iDim;
    unsigned long iVertex, iPoint;
    su2double a_ref[3], b_ref[3], AC_ref[3], BD_ref[3];
    su2double a_curr[3], b_curr[3], AC_curr[3], BD_curr[3];
    
    /*--- Determine whether the load conditions are applied in the reference or in the current configuration ---*/
    
    bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS); // Nonlinear analysis.
    
    for (iNode = 0; iNode < 4; iNode++) {
      for (iDim = 0; iDim < 3; iDim++) {
        nodeCoord_ref[iNode][iDim]  = 0.0;
        nodeCoord_curr[iNode][iDim] = 0.0;
      }
    }
    
    for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
      
      /*--- Identify the kind of boundary element ---*/
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == LINE)           nNodes = 2;
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE)       nNodes = 3;
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL)  nNodes = 4;
      
      /*--- Retrieve the boundary reference and current coordinates ---*/
      for (iNode = 0; iNode < nNodes; iNode++) {
        indexNode[iNode] = geometry->bound[val_marker][iElem]->GetNode(iNode);
        for (iDim = 0; iDim < nDim; iDim++) {
          val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
          val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
          /*--- Assign values to the container ---*/
          nodeCoord_ref[iNode][iDim]  = val_Coord;
          nodeCoord_curr[iNode][iDim] = val_Sol;
        }
      }
      
      /*--- We need the indices of the vertices, which are "Dual Grid Info" ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        for (iNode = 0; iNode < nNodes; iNode++) {
          if (iPoint == indexNode[iNode]) indexVertex[iNode] = iVertex;
        }
      }
      
      /*--- Retrieve the reference normal for one of the points. They go INSIDE the structural domain. ---*/
      nodal_normal = geometry->vertex[val_marker][indexVertex[0]]->GetNormal();
      Norm = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Norm += nodal_normal[iDim]*nodal_normal[iDim];
      }
      Norm = sqrt(Norm);
      for (iDim = 0; iDim < nDim; iDim++) {
        nodal_normal_unit[iDim] = nodal_normal[iDim] / Norm;
      }
      
      /*--- Compute area (3D), and length of the surfaces (2D), and the unitary normal vector in current configuration ---*/
      
      if (nDim == 2) {
        
        /*-- Compute the vector a in reference and current configurations ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          a_ref[iDim]  = nodeCoord_ref[0][iDim] -nodeCoord_ref[1][iDim];
          a_curr[iDim] = nodeCoord_curr[0][iDim]-nodeCoord_curr[1][iDim];
        }
        
        /*-- Compute the length of the boundary element in reference and current configurations ---*/
        Length_Elem_curr = sqrt(a_curr[0]*a_curr[0]+a_curr[1]*a_curr[1]);
        Length_Elem_ref  = sqrt(a_ref[0]*a_ref[0]+a_ref[1]*a_ref[1]);
        
        /*-- Compute the length of the boundary element in reference and current configurations ---*/
        normal_ref_unit[0] =   a_ref[1] /Length_Elem_ref;
        normal_ref_unit[1] = -(a_ref[0])/Length_Elem_ref;
        
        normal_curr_unit[0] =   a_curr[1] /Length_Elem_curr;
        normal_curr_unit[1] = -(a_curr[0])/Length_Elem_curr;
        
        /*-- Dot product to check the element orientation in the reference configuration ---*/
        dot_Prod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
        }
        
        /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
        /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
        /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
        if (dot_Prod < 0) {
          for (iDim = 0; iDim < nDim; iDim++) {
            normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
            normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
          }
        }
        
        if (linear_analysis) {
          Residual[0] = (1.0/2.0) * TotalLoad * Length_Elem_ref * normal_ref_unit[0];
          Residual[1] = (1.0/2.0) * TotalLoad * Length_Elem_ref * normal_ref_unit[1];
          
          node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
          node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
        }
        else if (nonlinear_analysis) {
          Residual[0] = (1.0/2.0) * TotalLoad * Length_Elem_curr * normal_curr_unit[0];
          Residual[1] = (1.0/2.0) * TotalLoad * Length_Elem_curr * normal_curr_unit[1];
          
          node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
          node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
        }
        
      }
      
      if (nDim == 3) {
        
        if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
          
          for (iDim = 0; iDim < nDim; iDim++) {
            a_ref[iDim] = nodeCoord_ref[1][iDim]-nodeCoord_ref[0][iDim];
            b_ref[iDim] = nodeCoord_ref[2][iDim]-nodeCoord_ref[0][iDim];
            
            a_curr[iDim] = nodeCoord_curr[1][iDim]-nodeCoord_curr[0][iDim];
            b_curr[iDim] = nodeCoord_curr[2][iDim]-nodeCoord_curr[0][iDim];
          }
          
          su2double Ni=0, Nj=0, Nk=0;
          
          /*--- Reference configuration ---*/
          Ni = a_ref[1]*b_ref[2] - a_ref[2]*b_ref[1];
          Nj = a_ref[2]*b_ref[0] - a_ref[0]*b_ref[2];
          Nk = a_ref[0]*b_ref[1] - a_ref[1]*b_ref[0];
          
          Area_Elem_ref = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_ref_unit[0] = Ni / Area_Elem_ref;
          normal_ref_unit[1] = Nj / Area_Elem_ref;
          normal_ref_unit[2] = Nk / Area_Elem_ref;
          
          /*--- Current configuration ---*/
          Ni = a_curr[1]*b_curr[2] - a_curr[2]*b_curr[1];
          Nj = a_curr[2]*b_curr[0] - a_curr[0]*b_curr[2];
          Nk = a_curr[0]*b_curr[1] - a_curr[1]*b_curr[0];
          
          Area_Elem_curr = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_curr_unit[0] = Ni / Area_Elem_curr;
          normal_curr_unit[1] = Nj / Area_Elem_curr;
          normal_curr_unit[2] = Nk / Area_Elem_curr;
          
          /*-- Dot product to check the element orientation in the reference configuration ---*/
          dot_Prod = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
          }
          
          /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
          /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
          /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
          if (dot_Prod < 0) {
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
              normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
            }
          }
          
          if (linear_analysis) {
            Residual[0] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[0];
            Residual[1] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[1];
            Residual[2] = (1.0/3.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
          }
          else if (nonlinear_analysis) {
            Residual[0] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[0];
            Residual[1] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[1];
            Residual[2] = (1.0/3.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
          }
          
        }
        
        else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
          
          for (iDim = 0; iDim < nDim; iDim++) {
            AC_ref[iDim] = nodeCoord_ref[2][iDim]-nodeCoord_ref[0][iDim];
            BD_ref[iDim] = nodeCoord_ref[3][iDim]-nodeCoord_ref[1][iDim];
            
            AC_curr[iDim] = nodeCoord_curr[2][iDim]-nodeCoord_curr[0][iDim];
            BD_curr[iDim] = nodeCoord_curr[3][iDim]-nodeCoord_curr[1][iDim];
          }
          
          su2double Ni=0, Nj=0, Nk=0;
          
          /*--- Reference configuration ---*/
          Ni=AC_ref[1]*BD_ref[2]-AC_ref[2]*BD_ref[1];
          Nj=-AC_ref[0]*BD_ref[2]+AC_ref[2]*BD_ref[0];
          Nk=AC_ref[0]*BD_ref[1]-AC_ref[1]*BD_ref[0];
          
          Area_Elem_ref = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_ref_unit[0] = Ni / Area_Elem_ref;
          normal_ref_unit[1] = Nj / Area_Elem_ref;
          normal_ref_unit[2] = Nk / Area_Elem_ref;
          
          /*--- Current configuration ---*/
          Ni=AC_curr[1]*BD_curr[2]-AC_curr[2]*BD_curr[1];
          Nj=-AC_curr[0]*BD_curr[2]+AC_curr[2]*BD_curr[0];
          Nk=AC_curr[0]*BD_curr[1]-AC_curr[1]*BD_curr[0];
          
          Area_Elem_curr = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
          
          normal_curr_unit[0] = Ni / Area_Elem_curr;
          normal_curr_unit[1] = Nj / Area_Elem_curr;
          normal_curr_unit[2] = Nk / Area_Elem_curr;
          
          /*-- Dot product to check the element orientation in the reference configuration ---*/
          dot_Prod = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dot_Prod += normal_ref_unit[iDim] * nodal_normal_unit[iDim];
          }
          
          /*--- If dot_Prod > 0, the normal goes inside the structural domain. ---*/
          /*--- If dot_Prod < 0, the normal goes outside the structural domain. ---*/
          /*--- We adopt the criteria of the normal going inside the domain, so if dot_Prod < 1, we change the orientation. ---*/
          if (dot_Prod < 0) {
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_ref_unit[iDim]  = -1.0*normal_ref_unit[iDim];
              normal_curr_unit[iDim] = -1.0*normal_curr_unit[iDim];
            }
          }
          
          if (linear_analysis) {
            Residual[0] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[0];
            Residual[1] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[1];
            Residual[2] = (1.0/4.0) * TotalLoad * Area_Elem_ref * normal_ref_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[3]]->Add_SurfaceLoad_Res(Residual);
          }
          else if (nonlinear_analysis) {
            Residual[0] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[0];
            Residual[1] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[1];
            Residual[2] = (1.0/4.0) * TotalLoad * Area_Elem_curr * normal_curr_unit[2];
            
            node[indexNode[0]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[1]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[2]]->Add_SurfaceLoad_Res(Residual);
            node[indexNode[3]]->Add_SurfaceLoad_Res(Residual);
          }
          
        }
        
      }
      
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {
  
  su2double a[3], b[3], AC[3], BD[3];
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3=0;
  su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
  su2double Length_Elem = 0.0, Area_Elem = 0.0;
  unsigned short iDim;
  
  su2double LoadDirVal = config->GetLoad_Dir_Value(config->GetMarker_All_TagBound(val_marker));
  su2double LoadDirMult = config->GetLoad_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
  su2double *Load_Dir_Local= config->GetLoad_Dir(config->GetMarker_All_TagBound(val_marker));
  
  su2double TotalLoad;
  
  bool Sigmoid_Load = config->GetSigmoid_Load();
  su2double Sigmoid_Time = config->GetSigmoid_Time();
  su2double Sigmoid_K = config->GetSigmoid_K();
  su2double SigAux = 0.0;
  
  su2double CurrentTime=config->GetCurrent_DynTime();
  su2double ModAmpl, NonModAmpl;
  
  bool Ramp_Load = config->GetRamp_Load();
  su2double Ramp_Time = config->GetRamp_Time();
  
  if (Ramp_Load) {
    ModAmpl=LoadDirVal*LoadDirMult*CurrentTime/Ramp_Time;
    NonModAmpl=LoadDirVal*LoadDirMult;
    TotalLoad=min(ModAmpl,NonModAmpl);
  }
  else if (Sigmoid_Load) {
    SigAux = CurrentTime/ Sigmoid_Time;
    ModAmpl = (1 / (1+exp(-1*Sigmoid_K*(SigAux - 0.5)) ) );
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
    TotalLoad=ModAmpl*LoadDirVal*LoadDirMult;
  }
  else {
    TotalLoad=LoadDirVal*LoadDirMult;
  }
  
  /*--- Compute the norm of the vector that was passed in the config file ---*/
  su2double Norm = 1.0;
  if (nDim==2) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]);
  if (nDim==3) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]+Load_Dir_Local[2]*Load_Dir_Local[2]);
  
  for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {
    
    Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);     Coord_0 = geometry->node[Point_0]->GetCoord();
    Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);     Coord_1 = geometry->node[Point_1]->GetCoord();
    if (nDim == 3) {
      
      Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);  Coord_2 = geometry->node[Point_2]->GetCoord();
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);  Coord_3 = geometry->node[Point_3]->GetCoord();
      }
      
    }
    
    /*--- Compute area (3D), and length of the surfaces (2D) ---*/
    
    if (nDim == 2) {
      
      for (iDim = 0; iDim < nDim; iDim++) a[iDim] = Coord_0[iDim]-Coord_1[iDim];
      
      Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
      //      Normal_Elem[0] =   a[1];
      //      Normal_Elem[1] = -(a[0]);
      
    }
    
    if (nDim == 3) {
      
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
        
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_1[iDim]-Coord_0[iDim];
          b[iDim] = Coord_2[iDim]-Coord_0[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=a[1]*b[2]-a[2]*b[1];
        Nj=-a[0]*b[2]+a[2]*b[0];
        Nk=a[0]*b[1]-a[1]*b[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
      }
      
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        
        for (iDim = 0; iDim < nDim; iDim++) {
          AC[iDim] = Coord_2[iDim]-Coord_0[iDim];
          BD[iDim] = Coord_3[iDim]-Coord_1[iDim];
        }
        
        su2double Ni=0 , Nj=0, Nk=0;
        
        Ni=AC[1]*BD[2]-AC[2]*BD[1];
        Nj=-AC[0]*BD[2]+AC[2]*BD[0];
        Nk=AC[0]*BD[1]-AC[1]*BD[0];
        
        Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);
        
      }
    }
    
    if (nDim == 2) {
      
      Residual[0] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
      Residual[1] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
      
      node[Point_0]->Add_SurfaceLoad_Res(Residual);
      node[Point_1]->Add_SurfaceLoad_Res(Residual);
      
    }
    
    else {
      if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE) {
        
        Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);
        node[Point_2]->Add_SurfaceLoad_Res(Residual);
        
      }
      else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == QUADRILATERAL) {
        
        Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
        Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;
        
        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);
        node[Point_2]->Add_SurfaceLoad_Res(Residual);
        node[Point_3]->Add_SurfaceLoad_Res(Residual);
        
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                         unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) { }

void CFEM_ElasticitySolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CFEM_ElasticitySolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint, jPoint;
  unsigned short iVar, jVar;
  
  bool initial_calc = (config->GetExtIter() == 0);                  // Checks if it is the first calculation.
  bool first_iter = (config->GetIntIter() == 0);
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation.
  
  bool body_forces = config->GetDeadLoad();                      // Body forces (dead loads).
  
  bool restart = config->GetRestart();                          // Restart solution
  bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter());  // Restart iteration
  
  bool incremental_load = config->GetIncrementalLoad();
  
  if (!dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /*--- Add the external contribution to the residual    ---*/
      /*--- (the terms that are constant over the time step) ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
    }
    
  }
  
  if (dynamic) {
    
    /*--- Add the mass matrix contribution to the Jacobian ---*/
    
    /*
     * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
     * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
     * iteration, as the Jacobian is recomputed
     *
     * If the problem is linear, we add the Mass Matrix contribution to the Jacobian at the first calculation.
     * From then on, the Jacobian is always the same matrix.
     *
     */
    
    if ((nonlinear_analysis && (newton_raphson || first_iter)) ||
        (linear_analysis && initial_calc) ||
        (linear_analysis && restart && initial_calc_restart)) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (jPoint = 0; jPoint < nPoint; jPoint++) {
          for(iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Jacobian_ij[iVar][jVar] = a_dt[0] * MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
            }
          }
          Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
        }
      }
    }
    
    
    /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
    if (linear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+    //a0*U(t)
          a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+  //a2*U'(t)
          a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
        }
        TimeRes_Aux.SetBlock(iPoint, Residual);
      }
    }
    else if (nonlinear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
          - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
          + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
          + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
        }
        TimeRes_Aux.SetBlock(iPoint, Residual);
      }
      
    }
    
    /*--- Once computed, compute M*TimeRes_Aux ---*/
    MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
    /*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
      /*--- Dynamic contribution ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
      }
      //Res_Time_Cont = TimeRes.GetBlock(iPoint);
      LinSysRes.AddBlock(iPoint, Res_Time_Cont);
      
      /*--- External surface load contribution ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      
      /*--- Body forces contribution (dead load) ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
      /*--- FSI contribution (flow loads) ---*/
      if (fsi) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = loadIncrement * node[iPoint]->Get_FlowTraction(iVar);
          }
        }
        else {
          Res_FSI_Cont = node[iPoint]->Get_FlowTraction();
        }
        LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
      }
    }
  }
  
  
}

void CFEM_ElasticitySolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool linear = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);    // Geometrically linear problems
  bool nonlinear = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Geometrically non-linear problems
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);          // Dynamic simulations.
  
  /*--- Update solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Displacements component of the solution ---*/
      
      /*--- If it's a non-linear problem, the result is the DELTA_U, not U itself ---*/
      
      if (linear) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
      if (nonlinear)  node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
    }
    
  }
  
  if (dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Acceleration component of the solution ---*/
        /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/
        
        Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
                                node[iPoint]->GetSolution_time_n(iVar)) -
        a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
        a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
      }
      
      /*--- Set the acceleration in the node structure ---*/
      
      node[iPoint]->SetSolution_Accel(Solution);
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Velocity component of the solution ---*/
        /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/
        
        Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
        a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
        a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);
        
      }
      
      /*--- Set the velocity in the node structure ---*/
      
      node[iPoint]->SetSolution_Vel(Solution);
      
    }
    
  }
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  
}

void CFEM_ElasticitySolver::ImplicitNewmark_Relaxation(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  su2double *valSolutionPred;
  
  /*--- Update solution and set it to be the solution after applying relaxation---*/
  
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    
    valSolutionPred = node[iPoint]->GetSolution_Pred();
    
    node[iPoint]->SetSolution(valSolutionPred);
  }
  
  /*--- Compute velocities and accelerations ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Acceleration component of the solution ---*/
      /*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/
      
      Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
                              node[iPoint]->GetSolution_time_n(iVar)) -
      a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
      a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
    }
    
    /*--- Set the acceleration in the node structure ---*/
    
    node[iPoint]->SetSolution_Accel(Solution);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Velocity component of the solution ---*/
      /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/
      
      Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
      a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
      a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);
      
    }
    
    /*--- Set the velocity in the node structure ---*/
    
    node[iPoint]->SetSolution_Vel(Solution);
    
  }
  
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- After the solution has been communicated, set the 'old' predicted solution as the solution ---*/
  /*--- Loop over n points (as we have already communicated everything ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->SetSolution_Pred_Old(iVar,node[iPoint]->GetSolution(iVar));
    }
  }
  
  
}


void CFEM_ElasticitySolver::GeneralizedAlpha_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint, jPoint;
  unsigned short iVar, jVar;
  
  bool initial_calc = (config->GetExtIter() == 0);                  // Checks if it is the first calculation.
  bool first_iter = (config->GetIntIter() == 0);
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);    // Newton-Raphson method
  bool fsi = config->GetFSI_Simulation();                        // FSI simulation.
  
  bool body_forces = config->GetDeadLoad();                      // Body forces (dead loads).
  
  bool restart = config->GetRestart();                          // Restart solution
  bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter());  // Restart iteration
  
  su2double alpha_f = config->Get_Int_Coeffs(2);
  
  bool incremental_load = config->GetIncrementalLoad();
  
  if (!dynamic) {
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /*--- Add the external contribution to the residual    ---*/
      /*--- (the terms that are constant over the time step) ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = node[iPoint]->Get_SurfaceLoad_Res(iVar);
        }
        //Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
      }
      
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
    }
    
  }
  
  if (dynamic) {
    
    /*--- Add the mass matrix contribution to the Jacobian ---*/
    
    /*
     * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
     * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
     * iteration, as the Jacobian is recomputed
     *
     * If the problem is linear, we add the Mass Matrix contribution to the Jacobian at the first calculation.
     * From then on, the Jacobian is always the same matrix.
     *
     */
    
    if ((nonlinear_analysis && (newton_raphson || first_iter)) ||
        (linear_analysis && initial_calc) ||
        (linear_analysis && restart && initial_calc_restart)) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (jPoint = 0; jPoint < nPoint; jPoint++) {
          for(iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Jacobian_ij[iVar][jVar] = a_dt[0] * MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
            }
          }
          Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
        }
      }
    }
    
    
    /*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
    if (linear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+    //a0*U(t)
          a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+  //a2*U'(t)
          a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
        }
        TimeRes_Aux.SetBlock(iPoint, Residual);
      }
    }
    else if (nonlinear_analysis) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)        //a0*U(t)
          - a_dt[0]*node[iPoint]->GetSolution(iVar)           //a0*U(t+dt)(k-1)
          + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)    //a2*U'(t)
          + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);  //a3*U''(t)
        }
        TimeRes_Aux.SetBlock(iPoint, Residual);
      }
    }
    /*--- Once computed, compute M*TimeRes_Aux ---*/
    MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
    /*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      /*--- Dynamic contribution ---*/
      //Res_Time_Cont = TimeRes.GetBlock(iPoint);
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Time_Cont[iVar] = TimeRes.GetBlock(iPoint, iVar);
      }
      LinSysRes.AddBlock(iPoint, Res_Time_Cont);
      /*--- External surface load contribution ---*/
      if (incremental_load) {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = loadIncrement * ( (1 - alpha_f) * node[iPoint]->Get_SurfaceLoad_Res(iVar) +
                                                alpha_f  * node[iPoint]->Get_SurfaceLoad_Res_n(iVar) );
        }
      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          Res_Ext_Surf[iVar] = (1 - alpha_f) * node[iPoint]->Get_SurfaceLoad_Res(iVar) +
          alpha_f  * node[iPoint]->Get_SurfaceLoad_Res_n(iVar);
        }
      }
      LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
      
      /*--- Add the contribution to the residual due to body forces.
       *--- It is constant over time, so it's not necessary to distribute it. ---*/
      
      if (body_forces) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = loadIncrement * node[iPoint]->Get_BodyForces_Res(iVar);
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_Dead_Load[iVar] = node[iPoint]->Get_BodyForces_Res(iVar);
          }
          //Res_Dead_Load = node[iPoint]->Get_BodyForces_Res();
        }
        
        LinSysRes.AddBlock(iPoint, Res_Dead_Load);
      }
      
      /*--- Add FSI contribution ---*/
      if (fsi) {
        if (incremental_load) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = loadIncrement * ( (1 - alpha_f) * node[iPoint]->Get_FlowTraction(iVar) +
                                                  alpha_f  * node[iPoint]->Get_FlowTraction_n(iVar) );
          }
        }
        else {
          for (iVar = 0; iVar < nVar; iVar++) {
            Res_FSI_Cont[iVar] = (1 - alpha_f) * node[iPoint]->Get_FlowTraction(iVar) +
            alpha_f  * node[iPoint]->Get_FlowTraction_n(iVar);
          }
        }
        LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
      }
    }
  }
  
}

void CFEM_ElasticitySolver::GeneralizedAlpha_UpdateDisp(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool linear = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);    // Geometrically linear problems
  bool nonlinear = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Geometrically non-linear problems
  
  /*--- Update solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Displacements component of the solution ---*/
      
      /*--- If it's a non-linear problem, the result is the DELTA_U, not U itself ---*/
      
      if (linear) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
      if (nonlinear)  node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      
    }
    
  }
  
  /*--- Perform the MPI communication of the solution, displacements only ---*/
  
  Set_MPI_Solution_DispOnly(geometry, config);
  
}

void CFEM_ElasticitySolver::GeneralizedAlpha_UpdateSolution(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double alpha_f = config->Get_Int_Coeffs(2), alpha_m =  config->Get_Int_Coeffs(3);
  
  /*--- Compute solution at t_n+1, and update velocities and accelerations ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Compute the solution from the previous time step and the solution computed at t+1-alpha_f ---*/
      /*--- U(t+dt) = 1/alpha_f*(U(t+1-alpha_f)-alpha_f*U(t)) ---*/
      
      Solution[iVar]=(1 / (1 - alpha_f))*(node[iPoint]->GetSolution(iVar) -
                                          alpha_f * node[iPoint]->GetSolution_time_n(iVar));
      

    }
    
    /*--- Set the solution in the node structure ---*/
    
    node[iPoint]->SetSolution(Solution);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Acceleration component of the solution ---*/
      /*--- U''(t+dt-alpha_m) = a8*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/
      
      Solution_Interm[iVar]=a_dt[8]*( node[iPoint]->GetSolution(iVar) -
                                     node[iPoint]->GetSolution_time_n(iVar)) -
      a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
      a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
      
      /*--- Compute the solution from the previous time step and the solution computed at t+1-alpha_f ---*/
      /*--- U''(t+dt) = 1/alpha_m*(U''(t+1-alpha_m)-alpha_m*U''(t)) ---*/
      
      Solution[iVar]=(1 / (1 - alpha_m))*(Solution_Interm[iVar] - alpha_m * node[iPoint]->GetSolution_Accel_time_n(iVar));
    }
    
    /*--- Set the acceleration in the node structure ---*/
    
    node[iPoint]->SetSolution_Accel(Solution);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      /*--- Velocity component of the solution ---*/
      /*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/
      
      Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
      a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
      a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);
      
    }
    
    /*--- Set the velocity in the node structure ---*/
    
    node[iPoint]->SetSolution_Vel(Solution);
    
  }
  
  /*--- Perform the MPI communication of the solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
}

void CFEM_ElasticitySolver::GeneralizedAlpha_UpdateLoads(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long iPoint;
  bool fsi = config->GetFSI_Simulation();
  
  /*--- Set the load conditions of the time step n+1 as the load conditions for time step n ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->Set_SurfaceLoad_Res_n();
    if (fsi) node[iPoint]->Set_FlowTraction_n();
  }
  
}

void CFEM_ElasticitySolver::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned long IterLinSol = 0, iPoint, total_index;
  unsigned short iVar;
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
    
  }
  
  CSysSolve femSystem;
  IterLinSol = femSystem.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- The the number of iterations of the linear solver ---*/
  
  SetIterLinSolver(IterLinSol);
  
}



void CFEM_ElasticitySolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry,
                                        CGeometry **flow_geometry, CConfig *fea_config,
                                        CConfig *flow_config, CNumerics *fea_numerics) {
  
  unsigned short nMarkerFSI, nMarkerFlow;    // Number of markers on FSI problem, FEA and Flow side
  unsigned short iMarkerFSI, iMarkerFlow;    // Variables for iteration over markers
  int Marker_Flow = -1;
  
  unsigned long iVertex, iPoint;                // Variables for iteration over vertices and nodes
  
  unsigned short iDim, jDim;
  
  // Check the kind of fluid problem
  bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
  bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
                             (flow_config->GetKind_Solver() == RANS) );
  
  /*--- Redimensionalize the pressure ---*/
  
  su2double *Velocity_ND, *Velocity_Real;
  su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;
  su2double factorForces;
  
  Velocity_Real = flow_config->GetVelocity_FreeStream();
  Density_Real = flow_config->GetDensity_FreeStream();
  
  Velocity_ND = flow_config->GetVelocity_FreeStreamND();
  Density_ND = flow_config->GetDensity_FreeStreamND();
  
  Velocity2_Real = 0.0;
  Velocity2_ND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
    Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
  }
  
  factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);
  
  /*--- Apply a ramp to the transfer of the fluid loads ---*/
  
  su2double ModAmpl;
  su2double CurrentTime = fea_config->GetCurrent_DynTime();
  su2double Static_Time = fea_config->GetStatic_Time();
  
  bool Ramp_Load = fea_config->GetRamp_Load();
  su2double Ramp_Time = fea_config->GetRamp_Time();
  
  if (CurrentTime <= Static_Time) { ModAmpl=0.0; }
  else if((CurrentTime > Static_Time) &&
          (CurrentTime <= (Static_Time + Ramp_Time)) &&
          (Ramp_Load)) {
    ModAmpl = (CurrentTime-Static_Time) / Ramp_Time;
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
  }
  else { ModAmpl = 1.0; }
  
  /*--- Number of markers on the FSI interface ---*/
  
  nMarkerFSI = (fea_config->GetMarker_n_FSIinterface())/2;
  nMarkerFlow = flow_geometry[MESH_0]->GetnMarker();    // Retrieve total number of markers on Fluid side
  
  // Parameters for the calculations
  // Pn: Pressure
  // Pinf: Pressure_infinite
  // div_vel: Velocity divergence
  // Dij: Dirac delta
  su2double Pn = 0.0, Pinf = 0.0, div_vel = 0.0, Dij = 0.0;
  su2double Viscosity = 0.0;
  su2double **Grad_PrimVar = NULL;
  su2double Tau[3][3];
  
  unsigned long Point_Flow, Point_Struct;
  su2double *Normal_Flow;
  
  su2double *tn_f;
  tn_f         = new su2double [nVar];      // Fluid traction
  
#ifndef HAVE_MPI
  
  unsigned long nVertexFlow;            // Number of vertices on FEA and Flow side
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->Clear_FlowTraction();
  }
  
  /*--- Loop over all the markers on the interface ---*/
  
  for (iMarkerFSI = 0; iMarkerFSI < nMarkerFSI; iMarkerFSI++) {
    
    /*--- Identification of the markers ---*/
    
    /*--- Current fluid marker ---*/
    for (iMarkerFlow = 0; iMarkerFlow < nMarkerFlow; iMarkerFlow++) {
      if (flow_config->GetMarker_All_FSIinterface(iMarkerFlow) == (iMarkerFSI+1)) {
        Marker_Flow = iMarkerFlow;
      }
    }
    
    nVertexFlow = flow_geometry[MESH_0]->GetnVertex(Marker_Flow);  // Retrieve total number of vertices on Fluid marker
    
    /*--- Loop over the nodes in the fluid mesh, calculate the tf vector (unitary) ---*/
    /*--- Here, we are looping over the fluid, and we find the pointer to the structure (Point_Struct) ---*/
    for (iVertex = 0; iVertex < nVertexFlow; iVertex++) {
      
      // Node from the flow mesh
      Point_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNode();
      
      // Normals at the vertex: these normals go inside the fluid domain.
      Normal_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNormal();
      
      // Corresponding node on the structural mesh
      Point_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorPoint();
      
      // Retrieve the values of pressure, viscosity and density
      if (incompressible) {
        
        Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
        Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();
        
        if (viscous_flow) {
          
          Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
          Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
        }
      }
      else if (compressible) {
        
        Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
        Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();
        
        if (viscous_flow) {
          
          Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
          Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
        }
      }
      
      // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
      for (iDim = 0; iDim < nDim; iDim++) {
        tn_f[iDim] = -(Pn-Pinf)*Normal_Flow[iDim];
      }
      
      // Calculate tn in the fluid nodes for the viscous term
      
      if (viscous_flow) {
        
        // Divergence of the velocity
        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
        if (incompressible) div_vel = 0.0;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          
          for (jDim = 0 ; jDim < nDim; jDim++) {
            // Dirac delta
            Dij = 0.0; if (iDim == jDim) Dij = 1.0;
            
            // Viscous stress
            Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
            TWO3*Viscosity*div_vel*Dij;
            
            // Viscous component in the tn vector --> Units of force (non-dimensional).
            tn_f[iDim] += Tau[iDim][jDim]*Normal_Flow[jDim];
          }
        }
      }
      
      // Rescale tn to SI units and apply time-dependent coefficient (static structure, ramp load, full load)
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Residual[iDim] = tn_f[iDim]*factorForces*ModAmpl;
      }
      
      /*--- Set the Flow traction ---*/
      //node[Point_Struct]->Set_FlowTraction(Residual);
      /*--- Add to the Flow traction (to add values to corners...) ---*/
      node[Point_Struct]->Add_FlowTraction(Residual);
    }
    
  }
  
#else
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  unsigned long nLocalVertexStruct = 0, nLocalVertexFlow = 0;
  
  unsigned long MaxLocalVertexStruct = 0, MaxLocalVertexFlow = 0;
  
  unsigned long nBuffer_FlowTraction = 0, nBuffer_StructTraction = 0;
  unsigned long nBuffer_DonorIndices = 0, nBuffer_SetIndex = 0;
  
  unsigned long Processor_Struct;
  
  int iProcessor, nProcessor = 0;
  
  unsigned short nMarkerStruct, iMarkerStruct;    // Variables for iteration over markers
  int Marker_Struct = -1;
  
  /*--- Number of markers on the FSI interface ---*/
  
  nMarkerFSI     = (flow_config->GetMarker_n_FSIinterface())/2;
  nMarkerStruct  = fea_geometry[MESH_0]->GetnMarker();
  nMarkerFlow    = flow_geometry[MESH_0]->GetnMarker();
  
  nProcessor = size;
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->Clear_FlowTraction();
  }
  
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
     *--- We need to loop over all markers on both sides and get the number of nodes
     *--- that belong to each FSI marker for each processor ---*/
    
    /*--- On the structural side ---*/
    
    for (iMarkerStruct = 0; iMarkerStruct < nMarkerStruct; iMarkerStruct++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerStruct) equals the index we are looping at ---*/
      if ( fea_config->GetMarker_All_FSIinterface(iMarkerStruct) == iMarkerFSI ) {
        /*--- We have identified the local index of the FEA marker ---*/
        /*--- Store the number of local points that belong to Marker_Struct on each processor ---*/
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
        /*--- Store the number of local points that belong to Marker_Flow on each processor ---*/
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
    
    Buffer_Send_nVertexStruct[0] = nLocalVertexStruct;                  // Retrieve total number of vertices on FEA marker
    Buffer_Send_nVertexFlow[0] = nLocalVertexFlow;                    // Retrieve total number of vertices on Flow marker
    if (rank == MASTER_NODE) Buffer_Recv_nVertexStruct = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
    if (rank == MASTER_NODE) Buffer_Recv_nVertexFlow = new unsigned long[size];     // Allocate memory to receive how many vertices are on each rank on the fluid side
    
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
    nBuffer_FlowTraction = MaxLocalVertexFlow * nDim;
    nBuffer_StructTraction = MaxLocalVertexStruct * nDim;
    
    /*--- We will be gathering donor index and donor processor (for flow -> donor = structure) ---*/
    /*--- Then we will pass on to the structural side the index (fea point) to the appropriate processor ---*/
    nBuffer_DonorIndices = 2 * MaxLocalVertexFlow;
    nBuffer_SetIndex = MaxLocalVertexStruct;
    
    /*--- Send and Recv buffers ---*/
    
    /*--- Buffers to send and receive the structural coordinates ---*/
    su2double *Buffer_Send_FlowTraction = new su2double[nBuffer_FlowTraction];
    su2double *Buffer_Recv_FlowTraction = NULL;
    
    /*--- Buffers to send and receive the donor index and processor ---*/
    long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
    long *Buffer_Recv_DonorIndices = NULL;
    
    /*--- Buffers to send and receive the new fluid coordinates ---*/
    su2double *Buffer_Send_StructTraction = NULL;
    su2double *Buffer_Recv_StructTraction = new su2double[nBuffer_StructTraction];
    
    /*--- Buffers to send and receive the fluid index ---*/
    long *Buffer_Send_SetIndex = NULL;
    long *Buffer_Recv_SetIndex = new long[nBuffer_SetIndex];
    
    /*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/
    
    if (rank == MASTER_NODE) {
      Buffer_Recv_FlowTraction  = new su2double[size*nBuffer_FlowTraction];
      Buffer_Recv_DonorIndices = new long[size*nBuffer_DonorIndices];
      Buffer_Send_StructTraction = new su2double[size*nBuffer_StructTraction];
      Buffer_Send_SetIndex     = new long[size*nBuffer_SetIndex];
    }
    
    /*--- On the fluid side ---*/
    
    /*--- If this processor owns the marker we are looping at on the structural side ---*/
    
    /*--- First we initialize all of the indices and processors to -1 ---*/
    /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
    for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
      Buffer_Send_DonorIndices[iVertex] = -1;
    
    if (Marker_Flow >= 0) {
      
      /*--- We have identified the local index of the FEA marker ---*/
      /*--- We loop over all the vertices in that marker and in that particular processor ---*/
      
      for (iVertex = 0; iVertex < nLocalVertexFlow; iVertex++) {
        
        Point_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNode();
        
        Point_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorPoint();
        
        Processor_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorProcessor();
        
        // Get the normal at the vertex: this normal goes inside the fluid domain.
        Normal_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNormal();
        
        // Retrieve the values of pressure, viscosity and density
        if (incompressible) {
          
          Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
          Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();
          
          if (viscous_flow) {
            
            Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
            Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
          }
        }
        else if (compressible) {
          
          Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
          Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();
          
          if (viscous_flow) {
            
            Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
            Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
          }
        }
        
        // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
        for (iDim = 0; iDim < nDim; iDim++) {
          tn_f[iDim] = -(Pn-Pinf)*Normal_Flow[iDim];
        }
        
        // Calculate tn in the fluid nodes for the viscous term
        
        if ((incompressible || compressible) && viscous_flow) {
          
          // Divergence of the velocity
          div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
          if (incompressible) div_vel = 0.0;
          
          for (iDim = 0; iDim < nDim; iDim++) {
            
            for (jDim = 0 ; jDim < nDim; jDim++) {
              // Dirac delta
              Dij = 0.0; if (iDim == jDim) Dij = 1.0;
              
              // Viscous stress
              Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
              TWO3*Viscosity*div_vel*Dij;
              
              // Viscous component in the tn vector --> Units of force (non-dimensional).
              tn_f[iDim] += Tau[iDim][jDim]*Normal_Flow[jDim];
            }
          }
        }
        
        for (iDim = 0; iDim < nDim; iDim++) {
          Buffer_Send_FlowTraction[iVertex*nDim+iDim] = tn_f[iDim]*factorForces*ModAmpl;
        }
        /*--- If this processor owns the node ---*/
        if (flow_geometry[MESH_0]->node[Point_Flow]->GetDomain()) {
          Buffer_Send_DonorIndices[2*iVertex]     = Point_Struct;
          Buffer_Send_DonorIndices[2*iVertex + 1] = Processor_Struct;
        }
        else {
          /*--- We set the values to be -1 to be able to identify them later as halo nodes ---*/
          Buffer_Send_DonorIndices[2*iVertex]     = -1;
          Buffer_Send_DonorIndices[2*iVertex + 1] = -1;
        }
        
      }
    }
    
    /*--- Once all the messages have been sent, we gather them all into the MASTER_NODE ---*/
    SU2_MPI::Gather(Buffer_Send_FlowTraction, nBuffer_FlowTraction, MPI_DOUBLE, Buffer_Recv_FlowTraction, nBuffer_FlowTraction, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
    //    if (rank == MASTER_NODE) {
    //      cout << endl << "-----------------------------------------------------------" << endl;
    //      cout << "For tag " << iMarkerFSI << ":" << endl;
    //      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    //        cout << "The processor " << iProcessor << " has " << Buffer_Recv_nVertexStruct[iProcessor] << " nodes on the structural side and ";
    //        cout << Buffer_Recv_nVertexFlow[iProcessor] << " nodes on the fluid side " << endl;
    //      }
    //      cout << "The max number of vertices is " << MaxLocalVertexStruct << " on the structural side and ";
    //      cout << MaxLocalVertexFlow << " on the fluid side." << endl;
    //
    //      cout << "---------------- Check received buffers ---------------------" << endl;
    //      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    //        long initialIndex, initialIndex2;
    //        initialIndex = iProcessor*nBuffer_FlowTraction;
    //        initialIndex2 = iProcessor*nBuffer_DonorIndices;
    //        for (long iCheck = 0; iCheck < Buffer_Recv_nVertexStruct[iProcessor]; iCheck++) {
    //          cout << "From processor " << iProcessor << " we get coordinates (";
    //            for (iDim = 0; iDim < nDim; iDim++)
    //              cout << Buffer_Recv_FlowTraction[initialIndex+iCheck*nDim+iDim] << ",";
    //          cout << "), the donor index for the flow " << Buffer_Recv_DonorIndices[initialIndex2+iCheck*2] ;
    //          cout << " and the donor processor " << Buffer_Recv_DonorIndices[initialIndex2+iCheck*2+1] << endl;
    //
    //        }
    //      }
    //
    //    }
    
    /*--- Counter to determine where in the array we have to set the information ---*/
    long *Counter_Processor_Struct = NULL;
    long iProcessor_Flow = 0, iIndex_Flow = 0;
    long iProcessor_Struct = 0, iPoint_Struct = 0, iIndex_Struct = 0;
    long Point_Struct_Send, Processor_Struct_Send;
    
    /*--- Now we pack the information to send it over to the different processors ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- We set the counter to 0 ---*/
      Counter_Processor_Struct = new long[nProcessor];
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        Counter_Processor_Struct[iProcessor] = 0;
      }
      
      /*--- First we initialize the index vector to -1 ---*/
      /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
      for (iVertex = 0; iVertex < nProcessor*nBuffer_SetIndex; iVertex++)
        Buffer_Send_SetIndex[iVertex] = -2;
      
      /*--- As of now we do the loop over the flow points ---*/
      /*--- The number of points for flow and structure does not necessarily have to match ---*/
      /*--- In fact, it's possible that a processor asks for nStruct nodes and there are only ---*/
      /*--- nFlow < nStruct available; this is due to halo nodes ---*/
      
      /*--- For every processor from which we have received information ---*/
      /*--- (This is, for every processor on the structural side) ---*/
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        
        /*--- This is the initial index on the coordinates buffer for that particular processor on the structural side ---*/
        iProcessor_Flow = iProcessor*nBuffer_FlowTraction;
        /*--- This is the initial index on the donor index/processor buffer for that particular processor on the structural side ---*/
        iIndex_Flow = iProcessor*nBuffer_DonorIndices;
        
        /*--- For every vertex in the information retreived from iProcessor ---*/
        for (iVertex = 0; iVertex < Buffer_Recv_nVertexFlow[iProcessor]; iVertex++) {
          
          /*--- The processor and index for the flow are: ---*/
          Processor_Struct_Send = Buffer_Recv_DonorIndices[iIndex_Flow+iVertex*2+1];
          Point_Struct_Send     = Buffer_Recv_DonorIndices[iIndex_Flow+iVertex*2];
          
          /*--- Load the buffer at the appropriate position ---*/
          /*--- This is determined on the fluid side by:
           *--- Processor_Flow*nBuffer_StructTraction -> Initial position of the processor array (fluid side)
           *--- +
           *--- Counter_Processor_Struct*nDim -> Initial position of the nDim array for the particular point on the fluid side
           *--- +
           *--- iDim -> Position within the nDim array that corresponds to a point
           *---
           *--- While on the structural side is:
           *--- iProcessor*nBuffer_FlowTraction -> Initial position on the processor array (structural side)
           *--- +
           *--- iVertex*nDim -> Initial position of the nDim array for the particular point on the structural side
           */
          
          /*--- We check that we are not setting the value for a halo node ---*/
          if (Point_Struct_Send != -1) {
            iProcessor_Struct = Processor_Struct_Send*nBuffer_StructTraction;
            iIndex_Struct = Processor_Struct_Send*nBuffer_SetIndex;
            iPoint_Struct = Counter_Processor_Struct[Processor_Struct_Send]*nDim;
            
            for (iDim = 0; iDim < nDim; iDim++)
              Buffer_Send_StructTraction[iProcessor_Struct + iPoint_Struct + iDim] = Buffer_Recv_FlowTraction[iProcessor_Flow + iVertex*nDim + iDim];
            
            /*--- We set the fluid index at an appropriate position matching the coordinates ---*/
            Buffer_Send_SetIndex[iIndex_Struct + Counter_Processor_Struct[Processor_Struct_Send]] = Point_Struct_Send;
            
            Counter_Processor_Struct[Processor_Struct_Send]++;
          }
          
        }
        
      }
      
      //      cout << "---------------- Check send buffers ---------------------" << endl;
      //
      //      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      //        long initialIndex, initialIndex2;
      //        initialIndex = iProcessor*nBuffer_StructTraction;
      //        initialIndex2 = iProcessor*nBuffer_SetIndex;
      //        for (long iCheck = 0; iCheck < Buffer_Recv_nVertexFlow[iProcessor]; iCheck++) {
      //          cout << "Processor " << iProcessor << " will receive the node " ;
      //          cout << Buffer_Send_SetIndex[initialIndex2+iCheck] << " which corresponds to the coordinates ";
      //          for (iDim = 0; iDim < nDim; iDim++)
      //            cout << "x" << iDim << "=" << Buffer_Send_StructTraction[initialIndex + iCheck*nDim + iDim] << ", ";
      //          cout << endl;
      //        }
      //
      //      }
      
    }
    
    /*--- Once all the messages have been prepared, we scatter them all from the MASTER_NODE ---*/
    SU2_MPI::Scatter(Buffer_Send_StructTraction, nBuffer_StructTraction, MPI_DOUBLE, Buffer_Recv_StructTraction, nBuffer_StructTraction, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Scatter(Buffer_Send_SetIndex, nBuffer_SetIndex, MPI_LONG, Buffer_Recv_SetIndex, nBuffer_SetIndex, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
    long indexPoint_iVertex, Point_Struct_Check;
    long Point_Struct_Recv;
    
    /*--- For the flow marker we are studying ---*/
    if (Marker_Struct >= 0) {
      
      /*--- We have identified the local index of the Structural marker ---*/
      /*--- We loop over all the vertices in that marker and in that particular processor ---*/
      
      for (iVertex = 0; iVertex < nLocalVertexStruct; iVertex++) {
        
        Point_Struct_Recv = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetNode();
        
        if (fea_geometry[MESH_0]->node[Point_Struct_Recv]->GetDomain()) {
          /*--- Find the index of the point Point_Struct in the buffer Buffer_Recv_SetIndex ---*/
          indexPoint_iVertex = std::distance(Buffer_Recv_SetIndex, std::find(Buffer_Recv_SetIndex, Buffer_Recv_SetIndex + MaxLocalVertexStruct, Point_Struct_Recv));
          
          Point_Struct_Check = Buffer_Recv_SetIndex[indexPoint_iVertex];
          
          if (Point_Struct_Check < 0) {
            cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
            exit(EXIT_FAILURE);
          }
          
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[iDim] = Buffer_Recv_StructTraction[indexPoint_iVertex*nDim+iDim];
          
          /*--- Add to the Flow traction ---*/
          node[Point_Struct_Recv]->Add_FlowTraction(Residual);
          
        }
        
      }
      
    }
    
    delete [] Buffer_Send_FlowTraction;
    delete [] Buffer_Send_DonorIndices;
    delete [] Buffer_Recv_StructTraction;
    delete [] Buffer_Recv_SetIndex;
    
    if (rank == MASTER_NODE) {
      delete [] Buffer_Recv_nVertexStruct;
      delete [] Buffer_Recv_nVertexFlow;
      delete [] Buffer_Recv_FlowTraction;
      delete [] Buffer_Recv_DonorIndices;
      delete [] Buffer_Send_StructTraction;
      delete [] Buffer_Send_SetIndex;
      delete [] Counter_Processor_Struct;
    }
    
  }
  
#endif
  
  delete[] tn_f;
  
  
}

void CFEM_ElasticitySolver::SetFEA_Load_Int(CSolver ***flow_solution, CGeometry **fea_geometry,
                                            CGeometry **flow_geometry, CConfig *fea_config,
                                            CConfig *flow_config, CNumerics *fea_numerics) { }

void CFEM_ElasticitySolver::PredictStruct_Displacement(CGeometry **fea_geometry,
                                                       CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned short predOrder = fea_config->GetPredictorOrder();
  su2double Delta_t = fea_config->GetDelta_DynTime();
  unsigned long iPoint, iDim;
  su2double *solDisp, *solVel, *solVel_tn, *valPred;
  
  //To nPointDomain: we need to communicate the predicted solution after setting it
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    if (predOrder==0) fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    else if (predOrder==1) {
      
      solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
      valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      
      for (iDim=0; iDim < nDim; iDim++) {
        valPred[iDim] = solDisp[iDim] + Delta_t*solVel[iDim];
      }
      
    }
    else if (predOrder==2) {
      
      solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
      solVel_tn = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n();
      valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      
      for (iDim=0; iDim < nDim; iDim++) {
        valPred[iDim] = solDisp[iDim] + 0.5*Delta_t*(3*solVel[iDim]-solVel_tn[iDim]);
      }
      
    }
    else {
      cout<< "Higher order predictor not implemented. Solving with order 0." << endl;
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    }
  }
  
}

void CFEM_ElasticitySolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config,
                                                      CSolver ***fea_solution, unsigned long iFSIIter) {
  
  unsigned long iPoint, iDim;
  su2double rbuf_numAitk = 0, sbuf_numAitk = 0;
  su2double rbuf_denAitk = 0, sbuf_denAitk = 0;
  
  su2double *dispPred, *dispCalc, *dispPred_Old, *dispCalc_Old;
  su2double deltaU[3] = {0.0, 0.0, 0.0}, deltaU_p1[3] = {0.0, 0.0, 0.0};
  su2double delta_deltaU[3] = {0.0, 0.0, 0.0};
  su2double CurrentTime=fea_config->GetCurrent_DynTime();
  su2double Static_Time=fea_config->GetStatic_Time();
  su2double WAitkDyn_tn1, WAitkDyn_Max, WAitkDyn_Min, WAitkDyn;
  
  unsigned short RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  ofstream historyFile_FSI;
  bool writeHistFSI = fea_config->GetWrite_Conv_FSI();
  if (writeHistFSI && (rank == MASTER_NODE)) {
    char cstrFSI[200];
    string filenameHistFSI = fea_config->GetConv_FileName_FSI();
    strcpy (cstrFSI, filenameHistFSI.data());
    historyFile_FSI.open (cstrFSI, std::ios_base::app);
  }
  
  
  /*--- Only when there is movement, and a dynamic coefficient is requested, it makes sense to compute the Aitken's coefficient ---*/
  
  if (CurrentTime > Static_Time) {
    
    if (RelaxMethod_FSI == NO_RELAXATION) {
      
      if (writeHistFSI && (rank == MASTER_NODE)) {
        
        SetWAitken_Dyn(1.0);
        
        if (iFSIIter == 0) historyFile_FSI << " " << endl ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
        if (iFSIIter == 0) historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << 1.0 ;
        else historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << 1.0 << "," ;
      }
      
    }
    else if (RelaxMethod_FSI == FIXED_PARAMETER) {
      
      if (writeHistFSI && (rank == MASTER_NODE)) {
        
        SetWAitken_Dyn(fea_config->GetAitkenStatRelax());
        
        if (iFSIIter == 0) historyFile_FSI << " " << endl ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
        historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
        if (iFSIIter == 0) historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << fea_config->GetAitkenStatRelax() ;
        else historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << fea_config->GetAitkenStatRelax() << "," ;
      }
      
    }
    else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {
      
      if (iFSIIter == 0) {
        
        WAitkDyn_tn1 = GetWAitken_Dyn_tn1();
        WAitkDyn_Max = fea_config->GetAitkenDynMaxInit();
        WAitkDyn_Min = fea_config->GetAitkenDynMinInit();
        
        WAitkDyn = min(WAitkDyn_tn1, WAitkDyn_Max);
        WAitkDyn = max(WAitkDyn, WAitkDyn_Min);
        
        SetWAitken_Dyn(WAitkDyn);
        if (writeHistFSI && (rank == MASTER_NODE)) {
          if (iFSIIter == 0) historyFile_FSI << " " << endl ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
          historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn ;
        }
        
      }
      else {
        // To nPointDomain; we need to communicate the values
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          
          dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
          dispPred_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred_Old();
          dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
          dispCalc_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Old();
          
          for (iDim = 0; iDim < nDim; iDim++) {
            
            /*--- Compute the deltaU and deltaU_n+1 ---*/
            deltaU[iDim] = dispCalc_Old[iDim] - dispPred_Old[iDim];
            deltaU_p1[iDim] = dispCalc[iDim] - dispPred[iDim];
            
            /*--- Compute the difference ---*/
            delta_deltaU[iDim] = deltaU_p1[iDim] - deltaU[iDim];
            
            /*--- Add numerator and denominator ---*/
            sbuf_numAitk += deltaU[iDim] * delta_deltaU[iDim];
            sbuf_denAitk += delta_deltaU[iDim] * delta_deltaU[iDim];
            
          }
          
        }
        
#ifdef HAVE_MPI
        SU2_MPI::Allreduce(&sbuf_numAitk, &rbuf_numAitk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(&sbuf_denAitk, &rbuf_denAitk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        rbuf_numAitk = sbuf_numAitk;
        rbuf_denAitk = sbuf_denAitk;
#endif
        
        WAitkDyn = GetWAitken_Dyn();
        
        if (rbuf_denAitk > 1E-15) {
          WAitkDyn = - 1.0 * WAitkDyn * rbuf_numAitk / rbuf_denAitk ;
        }
        
        WAitkDyn = max(WAitkDyn, 0.1);
        WAitkDyn = min(WAitkDyn, 1.0);
        
        SetWAitken_Dyn(WAitkDyn);
        
        if (writeHistFSI && (rank == MASTER_NODE)) {
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
          historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
          historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn << "," ;
        }
        
      }
      
    }
    else {
      if (rank == MASTER_NODE) cout << "No relaxation method used. " << endl;
    }
    
  }
  
  if (writeHistFSI && (rank == MASTER_NODE)) {historyFile_FSI.close();}
  
}

void CFEM_ElasticitySolver::SetAitken_Relaxation(CGeometry **fea_geometry,
                                                 CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned long iPoint, iDim;
  unsigned short RelaxMethod_FSI;
  su2double *dispPred, *dispCalc;
  su2double WAitken;
  su2double CurrentTime=fea_config->GetCurrent_DynTime();
  su2double Static_Time=fea_config->GetStatic_Time();
  
  RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();
  
  /*--- Only when there is movement it makes sense to update the solutions... ---*/
  
  if (CurrentTime > Static_Time) {
    
    if (RelaxMethod_FSI == NO_RELAXATION) {
      WAitken = 1.0;
    }
    else if (RelaxMethod_FSI == FIXED_PARAMETER) {
      WAitken = fea_config->GetAitkenStatRelax();
    }
    else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {
      WAitken = GetWAitken_Dyn();
    }
    else {
      WAitken = 1.0;
    }
    
    // To nPointDomain; we need to communicate the solutions (predicted, old and old predicted) after this routine
    for (iPoint=0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve pointers to the predicted and calculated solutions ---*/
      dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
      dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
      
      /*--- Set predicted solution as the old predicted solution ---*/
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred_Old();
      
      /*--- Set calculated solution as the old solution (needed for dynamic Aitken relaxation) ---*/
      fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Old(dispCalc);
      
      /*--- Apply the Aitken relaxation ---*/
      for (iDim=0; iDim < nDim; iDim++) {
        dispPred[iDim] = (1.0 - WAitken)*dispPred[iDim] + WAitken*dispCalc[iDim];
      }
      
    }
    
  }
  
}

void CFEM_ElasticitySolver::Update_StructSolution(CGeometry **fea_geometry,
                                                  CConfig *fea_config, CSolver ***fea_solution) {
  
  unsigned long iPoint;
  su2double *valSolutionPred;
  
  for (iPoint=0; iPoint < nPointDomain; iPoint++) {
    
    valSolutionPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
    
    fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution(valSolutionPred);
    
  }
  
  /*--- Perform the MPI communication of the solution, displacements only ---*/
  
  Set_MPI_Solution_DispOnly(fea_geometry[MESH_0], fea_config);
  
}

