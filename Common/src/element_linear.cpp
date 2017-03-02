/*!
 * \file element_linear.cpp
 * \brief Definition of the linear element structure for structural applications
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

#include "../include/element_structure.hpp"

CTRIA1::CTRIA1(void) : CElement() {
  
}

CTRIA1::CTRIA1(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).
  
  nNodes = 3;
  nGaussPoints = 1;
  
  nDimSq = nDim*nDim;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }
  
  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[3];
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
  CurrentCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    CurrentCoord [iNode] = new su2double[nDim];
  }
  
  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    RefCoord [iNode] = new su2double[nDim];
  }
  
  GaussWeight = new su2double [nGaussPoints];
  
  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussCoord [iGauss] = new su2double[nDim];
  }
  
  GaussCoord[0][0] = 0.333333333333333;  GaussCoord[0][1] = 0.333333333333333;  GaussWeight[0] = 0.5;
  
  Mab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Mab[iNode] = new su2double [nNodes];
  }
  
  Kab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kab [iNode][jNode] = new su2double[nDimSq];
    }
  }
  
  Ks_ab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Ks_ab[iNode] = new su2double [nNodes];
  }
  
  Kt_a = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kt_a[iNode] = new su2double [nDim];
  }
  
  if (body_forces) {
    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }
  }
  else {
    FDL_a = NULL;
  }
  
  su2double Xi, Eta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    
    val_Ni = Xi;
    GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = Eta;
    GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 1-Xi-Eta;
    GaussPoint[iGauss]->SetNi(val_Ni,2);
  }
  
  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant at a TRIA element ---*/
  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
  
}

CTRIA1::~CTRIA1(void) {
  
}

void CTRIA1::ComputeGrad_Linear(void) {
  
  su2double Jacobian[2][2], dNiXj[3][2];
  su2double detJac, GradNi_Xj;
  su2double ad[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = 1.0; 	dNiXj[0][1] = 0.0;
    dNiXj[1][0] = 0.0; 	dNiXj[1][1] = 1.0;
    dNiXj[2][0] = -1.0; 	dNiXj[2][1] = -1.0;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad[0][0] = Jacobian[1][1];
    ad[0][1] = -Jacobian[0][1];
    ad[1][0] = -Jacobian[1][0];
    ad[1][1] = Jacobian[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
}

void CTRIA1::ComputeGrad_NonLinear(void) {
  
  su2double Jac_Ref[2][2], Jac_Curr[2][2], dNiXj[3][2];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[2][2], ad_Curr[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = 1.0; 	dNiXj[0][1] = 0.0;
    dNiXj[1][0] = 0.0; 	dNiXj[1][1] = 1.0;
    dNiXj[2][0] = -1.0; 	dNiXj[2][1] = -1.0;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1];
    ad_Ref[0][1] = -Jac_Ref[0][1];
    ad_Ref[1][0] = -Jac_Ref[1][0];
    ad_Ref[1][1] = Jac_Ref[0][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1];
    ad_Curr[0][1] = -Jac_Curr[0][1];
    ad_Curr[1][0] = -Jac_Curr[1][0];
    ad_Curr[1][1] = Jac_Curr[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = ad_Ref[0][0]*ad_Ref[1][1]-ad_Ref[0][1]*ad_Ref[1][0];
    detJac_Curr = ad_Curr[0][0]*ad_Curr[1][1]-ad_Curr[0][1]*ad_Curr[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
  
}

CQUAD4::CQUAD4(void) : CElement() {
  
}

CQUAD4::CQUAD4(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).
  
  nNodes = 4;
  nGaussPoints = 4;
  
  nDimSq = nDim*nDim;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }
  
  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[3];
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
  CurrentCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    CurrentCoord [iNode] = new su2double[nDim];
  }
  
  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    RefCoord [iNode] = new su2double[nDim];
  }
  
  GaussWeight = new su2double [nGaussPoints];
  
  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussCoord [iGauss] = new su2double[nDim];
  }
  
  GaussCoord[0][0] = -0.577350269189626;  GaussCoord[0][1] = -0.577350269189626;  GaussWeight[0] = 1.0;
  GaussCoord[1][0] = 0.577350269189626;   GaussCoord[1][1] = -0.577350269189626;  GaussWeight[1] = 1.0;
  GaussCoord[2][0] = 0.577350269189626;   GaussCoord[2][1] = 0.577350269189626;   GaussWeight[2] = 1.0;
  GaussCoord[3][0] = -0.577350269189626;  GaussCoord[3][1] = 0.577350269189626;   GaussWeight[3] = 1.0;
  
  Mab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Mab[iNode] = new su2double [nNodes];
  }
  
  Kab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kab [iNode][jNode] = new su2double[nDimSq];
    }
  }
  
  Ks_ab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Ks_ab[iNode] = new su2double [nNodes];
  }
  
  Kt_a = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kt_a[iNode] = new su2double [nDim];
  }
  
  if (body_forces) {
    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }
  }
  else {
    FDL_a = NULL;
  }
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    
    val_Ni = 0.25*(1.0-Xi)*(1.0-Eta);		GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = 0.25*(1.0+Xi)*(1.0-Eta);		GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 0.25*(1.0+Xi)*(1.0+Eta);		GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = 0.25*(1.0-Xi)*(1.0+Eta);		GaussPoint[iGauss]->SetNi(val_Ni,3);
  }
  
  su2double ExtrapCoord[4][2];
  
  ExtrapCoord[0][0] = -1.732050807568877;  ExtrapCoord[0][1] = -1.732050807568877;
  ExtrapCoord[1][0] = 1.732050807568877;   ExtrapCoord[1][1] = -1.732050807568877;
  ExtrapCoord[2][0] = 1.732050807568877;   ExtrapCoord[2][1] = 1.732050807568877;
  ExtrapCoord[3][0] = -1.732050807568877;  ExtrapCoord[3][1] = 1.732050807568877;
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  for (iNode = 0; iNode < nNodes; iNode++) {
    Xi = ExtrapCoord[iNode][0];
    Eta = ExtrapCoord[iNode][1];
    
    NodalExtrap[iNode][0] = 0.25*(1.0-Xi)*(1.0-Eta);
    NodalExtrap[iNode][1] = 0.25*(1.0+Xi)*(1.0-Eta);
    NodalExtrap[iNode][2] = 0.25*(1.0+Xi)*(1.0+Eta);
    NodalExtrap[iNode][3] = 0.25*(1.0-Xi)*(1.0+Eta);
    
  }
  
}

CQUAD4::~CQUAD4(void) {
  
}

void CQUAD4::ComputeGrad_Linear(void) {
  
  su2double Xi, Eta;
  su2double Jacobian[2][2], dNiXj[4][2];
  su2double detJac, GradNi_Xj;
  su2double ad[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = -0.25*(1.0-Eta); dNiXj[0][1] = -0.25*(1.0-Xi);
    dNiXj[1][0] =  0.25*(1.0-Eta); dNiXj[1][1] = -0.25*(1.0+Xi);
    dNiXj[2][0] =  0.25*(1.0+Eta); dNiXj[2][1] =  0.25*(1.0+Xi);
    dNiXj[3][0] = -0.25*(1.0+Eta); dNiXj[3][1] =  0.25*(1.0-Xi);
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad[0][0] = Jacobian[1][1];
    ad[0][1] = -Jacobian[0][1];
    ad[1][0] = -Jacobian[1][0];
    ad[1][1] = Jacobian[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
}

void CQUAD4::ComputeGrad_NonLinear(void) {
  
  su2double Xi, Eta;
  su2double Jac_Ref[2][2], Jac_Curr[2][2], dNiXj[4][2];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[2][2], ad_Curr[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = -0.25*(1.0-Eta); dNiXj[0][1] = -0.25*(1.0-Xi);
    dNiXj[1][0] =  0.25*(1.0-Eta); dNiXj[1][1] = -0.25*(1.0+Xi);
    dNiXj[2][0] =  0.25*(1.0+Eta); dNiXj[2][1] =  0.25*(1.0+Xi);
    dNiXj[3][0] = -0.25*(1.0+Eta); dNiXj[3][1] =  0.25*(1.0-Xi);
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1];
    ad_Ref[0][1] = -Jac_Ref[0][1];
    ad_Ref[1][0] = -Jac_Ref[1][0];
    ad_Ref[1][1] = Jac_Ref[0][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1];
    ad_Curr[0][1] = -Jac_Curr[0][1];
    ad_Curr[1][0] = -Jac_Curr[1][0];
    ad_Curr[1][1] = Jac_Curr[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = ad_Ref[0][0]*ad_Ref[1][1]-ad_Ref[0][1]*ad_Ref[1][0];
    detJac_Curr = ad_Curr[0][0]*ad_Curr[1][1]-ad_Curr[0][1]*ad_Curr[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
}

CQUAD4P1::CQUAD4P1(void) : CQUAD4() {
  
  GaussPointP = NULL;
  GaussCoordP = NULL;
  GaussWeightP = NULL;
  Kk_ab = NULL;
  nGaussPointsP = 0;
  
}

CQUAD4P1::CQUAD4P1(unsigned short val_nDim, CConfig *config)
: CQUAD4(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  nGaussPointsP = 1;
  
  nDimSq = nDim*nDim;
  
  GaussPointP = new CGaussVariable*[nGaussPointsP];
  for (iGauss = 0; iGauss < nGaussPointsP; iGauss++) {
    GaussPointP[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  GaussWeightP = new su2double [nGaussPointsP];
  
  GaussCoordP = new su2double*[nGaussPointsP];
  for (iGauss = 0; iGauss < nGaussPointsP; iGauss++) {
    GaussCoordP [iGauss] = new su2double[nDim];
  }
  
  GaussCoordP[0][0] = 0.0;  GaussCoordP[0][1] = 0.0;  GaussWeightP[0] = 4.0;
  
  Kk_ab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kk_ab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kk_ab [iNode][jNode] = new su2double[nDimSq];
    }
  }
  
}

CQUAD4P1::~CQUAD4P1(void) {
  
}


void CQUAD4P1::ComputeGrad_Pressure(void) {
  
  su2double Xi, Eta;
  su2double Jac_Ref[2][2], Jac_Curr[2][2], dNiXj[4][2];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[2][2], ad_Curr[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPointsP; iGauss++) {
    
    Xi = GaussCoordP[iGauss][0];
    Eta = GaussCoordP[iGauss][1];
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = -0.25*(1.0-Eta); dNiXj[0][1] = -0.25*(1.0-Xi);
    dNiXj[1][0] =  0.25*(1.0-Eta); dNiXj[1][1] = -0.25*(1.0+Xi);
    dNiXj[2][0] =  0.25*(1.0+Eta); dNiXj[2][1] =  0.25*(1.0+Xi);
    dNiXj[3][0] = -0.25*(1.0+Eta); dNiXj[3][1] =  0.25*(1.0-Xi);
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < 4; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1];
    ad_Ref[0][1] = -Jac_Ref[0][1];
    ad_Ref[1][0] = -Jac_Ref[1][0];
    ad_Ref[1][1] = Jac_Ref[0][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1];
    ad_Curr[0][1] = -Jac_Curr[0][1];
    ad_Curr[1][0] = -Jac_Curr[1][0];
    ad_Curr[1][1] = Jac_Curr[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = ad_Ref[0][0]*ad_Ref[1][1]-ad_Ref[0][1]*ad_Ref[1][0];
    detJac_Curr = ad_Curr[0][0]*ad_Curr[1][1]-ad_Curr[0][1]*ad_Curr[1][0];
    
    GaussPointP[iGauss]->SetJ_X(detJac_Ref);
    GaussPointP[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPointP[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPointP[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
}

CTETRA1::CTETRA1(void) : CElement() {
  
}

CTETRA1::CTETRA1(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).
  
  nNodes = 4;
  nGaussPoints = 1;
  
  nDimSq = nDim*nDim;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }
  
  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[6];
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
  CurrentCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    CurrentCoord [iNode] = new su2double[nDim];
  }
  
  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    RefCoord [iNode] = new su2double[nDim];
  }
  
  GaussWeight = new su2double [nGaussPoints];
  
  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussCoord [iGauss] = new su2double[nDim];
  }
  
  GaussCoord[0][0] = 0.25;  GaussCoord[0][1] = 0.25; GaussCoord[0][2] = 0.25;  GaussWeight[0] = 0.166666666666666;
  
  Mab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Mab[iNode] = new su2double [nNodes];
  }
  
  Kab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kab [iNode][jNode] = new su2double[nDimSq];
    }
  }
  
  Ks_ab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Ks_ab[iNode] = new su2double [nNodes];
  }
  
  Kt_a = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kt_a[iNode] = new su2double [nDim];
  }
  
  if (body_forces) {
    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }
  }
  else {
    FDL_a = NULL;
  }
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, Zeta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];
    
    val_Ni = Xi;						GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = Eta;						GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 1.0 - Xi - Eta - Zeta;	GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = Zeta;					GaussPoint[iGauss]->SetNi(val_Ni,3);
  }
  
  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant at a TETRA element ---*/
  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
  NodalExtrap[3][0] = 1.0;
  
}

CTETRA1::~CTETRA1(void) {
  
}

void CTETRA1::ComputeGrad_Linear(void) {
  
  su2double Jacobian[3][3], dNiXj[4][3];
  su2double detJac, GradNi_Xj;
  su2double ad[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = 1.0;   dNiXj[0][1] = 0.0;   dNiXj[0][2] = 0.0;
    dNiXj[1][0] = 0.0;   dNiXj[1][1] = 1.0;   dNiXj[1][2] = 0.0;
    dNiXj[2][0] = -1.0;  dNiXj[2][1] = -1.0;  dNiXj[2][2] = -1.0;
    dNiXj[3][0] = 0.0;   dNiXj[3][1] = 0.0;   dNiXj[3][2] = 1.0;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad[0][0] = Jacobian[1][1]*Jacobian[2][2]-Jacobian[1][2]*Jacobian[2][1];
    ad[0][1] = Jacobian[0][2]*Jacobian[2][1]-Jacobian[0][1]*Jacobian[2][2];
    ad[0][2] = Jacobian[0][1]*Jacobian[1][2]-Jacobian[0][2]*Jacobian[1][1];
    ad[1][0] = Jacobian[1][2]*Jacobian[2][0]-Jacobian[1][0]*Jacobian[2][2];
    ad[1][1] = Jacobian[0][0]*Jacobian[2][2]-Jacobian[0][2]*Jacobian[2][0];
    ad[1][2] = Jacobian[0][2]*Jacobian[1][0]-Jacobian[0][0]*Jacobian[1][2];
    ad[2][0] = Jacobian[1][0]*Jacobian[2][1]-Jacobian[1][1]*Jacobian[2][0];
    ad[2][1] = Jacobian[0][1]*Jacobian[2][0]-Jacobian[0][0]*Jacobian[2][1];
    ad[2][2] = Jacobian[0][0]*Jacobian[1][1]-Jacobian[0][1]*Jacobian[1][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac = Jacobian[0][0]*ad[0][0]+Jacobian[0][1]*ad[1][0]+Jacobian[0][2]*ad[2][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
}

void CTETRA1::ComputeGrad_NonLinear(void) {
  
  su2double Jac_Ref[3][3], Jac_Curr[3][3], dNiXj[4][3];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[3][3], ad_Curr[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi, dN/d eta ---*/
    
    dNiXj[0][0] = 1.0;   dNiXj[0][1] = 0.0;   dNiXj[0][2] = 0.0;
    dNiXj[1][0] = 0.0;   dNiXj[1][1] = 1.0;   dNiXj[1][2] = 0.0;
    dNiXj[2][0] = -1.0;  dNiXj[2][1] = -1.0;  dNiXj[2][2] = -1.0;
    dNiXj[3][0] = 0.0;   dNiXj[3][1] = 0.0;   dNiXj[3][2] = 1.0;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1]*Jac_Ref[2][2]-Jac_Ref[1][2]*Jac_Ref[2][1];
    ad_Ref[0][1] = Jac_Ref[0][2]*Jac_Ref[2][1]-Jac_Ref[0][1]*Jac_Ref[2][2];
    ad_Ref[0][2] = Jac_Ref[0][1]*Jac_Ref[1][2]-Jac_Ref[0][2]*Jac_Ref[1][1];
    ad_Ref[1][0] = Jac_Ref[1][2]*Jac_Ref[2][0]-Jac_Ref[1][0]*Jac_Ref[2][2];
    ad_Ref[1][1] = Jac_Ref[0][0]*Jac_Ref[2][2]-Jac_Ref[0][2]*Jac_Ref[2][0];
    ad_Ref[1][2] = Jac_Ref[0][2]*Jac_Ref[1][0]-Jac_Ref[0][0]*Jac_Ref[1][2];
    ad_Ref[2][0] = Jac_Ref[1][0]*Jac_Ref[2][1]-Jac_Ref[1][1]*Jac_Ref[2][0];
    ad_Ref[2][1] = Jac_Ref[0][1]*Jac_Ref[2][0]-Jac_Ref[0][0]*Jac_Ref[2][1];
    ad_Ref[2][2] = Jac_Ref[0][0]*Jac_Ref[1][1]-Jac_Ref[0][1]*Jac_Ref[1][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1]*Jac_Curr[2][2]-Jac_Curr[1][2]*Jac_Curr[2][1];
    ad_Curr[0][1] = Jac_Curr[0][2]*Jac_Curr[2][1]-Jac_Curr[0][1]*Jac_Curr[2][2];
    ad_Curr[0][2] = Jac_Curr[0][1]*Jac_Curr[1][2]-Jac_Curr[0][2]*Jac_Curr[1][1];
    ad_Curr[1][0] = Jac_Curr[1][2]*Jac_Curr[2][0]-Jac_Curr[1][0]*Jac_Curr[2][2];
    ad_Curr[1][1] = Jac_Curr[0][0]*Jac_Curr[2][2]-Jac_Curr[0][2]*Jac_Curr[2][0];
    ad_Curr[1][2] = Jac_Curr[0][2]*Jac_Curr[1][0]-Jac_Curr[0][0]*Jac_Curr[1][2];
    ad_Curr[2][0] = Jac_Curr[1][0]*Jac_Curr[2][1]-Jac_Curr[1][1]*Jac_Curr[2][0];
    ad_Curr[2][1] = Jac_Curr[0][1]*Jac_Curr[2][0]-Jac_Curr[0][0]*Jac_Curr[2][1];
    ad_Curr[2][2] = Jac_Curr[0][0]*Jac_Curr[1][1]-Jac_Curr[0][1]*Jac_Curr[1][0];
    
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = Jac_Ref[0][0]*ad_Ref[0][0]+Jac_Ref[0][1]*ad_Ref[1][0]+Jac_Ref[0][2]*ad_Ref[2][0];
    detJac_Curr = Jac_Curr[0][0]*ad_Curr[0][0]+Jac_Curr[0][1]*ad_Curr[1][0]+Jac_Curr[0][2]*ad_Curr[2][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
}

CHEXA8::CHEXA8(void) : CElement() {
  
}

CHEXA8::CHEXA8(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).
  
  nNodes = 8;
  nGaussPoints = 8;
  
  nDimSq = nDim*nDim;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  NodalExtrap = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalExtrap[iNode] = new su2double[nGaussPoints];
  }
  
  NodalStress = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    NodalStress[iNode] = new su2double[6];
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
  CurrentCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    CurrentCoord [iNode] = new su2double[nDim];
  }
  
  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    RefCoord [iNode] = new su2double[nDim];
  }
  
  GaussWeight = new su2double [nGaussPoints];
  
  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussCoord [iGauss] = new su2double[nDim];
  }
  
  GaussCoord[0][0] = -0.577350269189626;  GaussCoord[0][1] = -0.577350269189626;  GaussCoord[0][2] = -0.577350269189626;	GaussWeight[0] = 1.0;
  GaussCoord[1][0] = 0.577350269189626;   GaussCoord[1][1] = -0.577350269189626;  GaussCoord[1][2] = -0.577350269189626;  GaussWeight[1] = 1.0;
  GaussCoord[2][0] = 0.577350269189626;   GaussCoord[2][1] = 0.577350269189626;  	GaussCoord[2][2] = -0.577350269189626;  GaussWeight[2] = 1.0;
  GaussCoord[3][0] = -0.577350269189626;  GaussCoord[3][1] = 0.577350269189626;  	GaussCoord[3][2] = -0.577350269189626;  GaussWeight[3] = 1.0;
  GaussCoord[4][0] = -0.577350269189626;  GaussCoord[4][1] = -0.577350269189626;  GaussCoord[4][2] = 0.577350269189626;  	GaussWeight[4] = 1.0;
  GaussCoord[5][0] = 0.577350269189626;   GaussCoord[5][1] = -0.577350269189626;  GaussCoord[5][2] = 0.577350269189626;  	GaussWeight[5] = 1.0;
  GaussCoord[6][0] = 0.577350269189626;   GaussCoord[6][1] = 0.577350269189626;  	GaussCoord[6][2] = 0.577350269189626;  	GaussWeight[6] = 1.0;
  GaussCoord[7][0] = -0.577350269189626;  GaussCoord[7][1] = 0.577350269189626;  	GaussCoord[7][2] = 0.577350269189626;  	GaussWeight[7] = 1.0;
  
  Mab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Mab[iNode] = new su2double [nNodes];
  }
  
  Kab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kab [iNode][jNode] = new su2double[nDimSq];
    }
  }
  
  Ks_ab = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Ks_ab[iNode] = new su2double [nNodes];
  }
  
  Kt_a = new su2double *[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kt_a[iNode] = new su2double [nDim];
  }
  
  if (body_forces) {
    FDL_a = new su2double *[nNodes];
    for (iNode = 0; iNode < nNodes; iNode++) {
      FDL_a[iNode] = new su2double [nDim];
    }
  }
  else {
    FDL_a = NULL;
  }
  
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, Eta, Zeta, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];
    
    val_Ni = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0-Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0-Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,1);
    val_Ni = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0-Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,2);
    val_Ni = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0-Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,3);
    val_Ni = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0+Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,4);
    val_Ni = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0+Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,5);
    val_Ni = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0+Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,6);
    val_Ni = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0+Zeta);		GaussPoint[iGauss]->SetNi(val_Ni,7);
  }
  
  
  su2double ExtrapCoord[8][3];
  
  ExtrapCoord[0][0] = -1.732050807568877;  ExtrapCoord[0][1] = -1.732050807568877;  	ExtrapCoord[0][2] = -1.732050807568877;
  ExtrapCoord[1][0] = 1.732050807568877;   ExtrapCoord[1][1] = -1.732050807568877;  	ExtrapCoord[1][2] = -1.732050807568877;
  ExtrapCoord[2][0] = 1.732050807568877;   ExtrapCoord[2][1] = 1.732050807568877;  	ExtrapCoord[2][2] = -1.732050807568877;
  ExtrapCoord[3][0] = -1.732050807568877;  ExtrapCoord[3][1] = 1.732050807568877;  	ExtrapCoord[3][2] = -1.732050807568877;
  ExtrapCoord[4][0] = -1.732050807568877;  ExtrapCoord[4][1] = -1.732050807568877;  	ExtrapCoord[4][2] = 1.732050807568877;
  ExtrapCoord[5][0] = 1.732050807568877;   ExtrapCoord[5][1] = -1.732050807568877;  	ExtrapCoord[5][2] = 1.732050807568877;
  ExtrapCoord[6][0] = 1.732050807568877;   ExtrapCoord[6][1] = 1.732050807568877;  	ExtrapCoord[6][2] = 1.732050807568877;
  ExtrapCoord[7][0] = -1.732050807568877;  ExtrapCoord[7][1] = 1.732050807568877;  	ExtrapCoord[7][2] = 1.732050807568877;
  
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  for (iNode = 0; iNode < nNodes; iNode++) {
    Xi = ExtrapCoord[iNode][0];
    Eta = ExtrapCoord[iNode][1];
    Zeta = ExtrapCoord[iNode][2];
    
    NodalExtrap[iNode][0] = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0-Zeta);
    NodalExtrap[iNode][1] = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0-Zeta);
    NodalExtrap[iNode][2] = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0-Zeta);
    NodalExtrap[iNode][3] = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0-Zeta);
    NodalExtrap[iNode][4] = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0+Zeta);
    NodalExtrap[iNode][5] = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0+Zeta);
    NodalExtrap[iNode][6] = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0+Zeta);
    NodalExtrap[iNode][7] = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0+Zeta);
  }
  
}

CHEXA8::~CHEXA8(void) {
  
}

void CHEXA8::ComputeGrad_Linear(void) {
  
  su2double Xi, Eta, Zeta;
  su2double Jacobian[3][3], dNiXj[8][3];
  su2double detJac, GradNi_Xj;
  su2double ad[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];
    
    /*--- dN/d xi ---*/
    
    dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);
    
    /*--- dN/d eta ---*/
    
    dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);
    
    /*--- dN/d mu ---*/
    
    dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);
    
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad[0][0] = Jacobian[1][1]*Jacobian[2][2]-Jacobian[1][2]*Jacobian[2][1];
    ad[0][1] = Jacobian[0][2]*Jacobian[2][1]-Jacobian[0][1]*Jacobian[2][2];
    ad[0][2] = Jacobian[0][1]*Jacobian[1][2]-Jacobian[0][2]*Jacobian[1][1];
    ad[1][0] = Jacobian[1][2]*Jacobian[2][0]-Jacobian[1][0]*Jacobian[2][2];
    ad[1][1] = Jacobian[0][0]*Jacobian[2][2]-Jacobian[0][2]*Jacobian[2][0];
    ad[1][2] = Jacobian[0][2]*Jacobian[1][0]-Jacobian[0][0]*Jacobian[1][2];
    ad[2][0] = Jacobian[1][0]*Jacobian[2][1]-Jacobian[1][1]*Jacobian[2][0];
    ad[2][1] = Jacobian[0][1]*Jacobian[2][0]-Jacobian[0][0]*Jacobian[2][1];
    ad[2][2] = Jacobian[0][0]*Jacobian[1][1]-Jacobian[0][1]*Jacobian[1][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac = Jacobian[0][0]*ad[0][0]+Jacobian[0][1]*ad[1][0]+Jacobian[0][2]*ad[2][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
  
}

void CHEXA8::ComputeGrad_NonLinear(void) {
  
  su2double Xi, Eta, Zeta;
  su2double Jac_Ref[3][3], Jac_Curr[3][3], dNiXj[8][3];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[3][3], ad_Curr[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];
    
    /*--- dN/d xi, dN/d eta ---*/
    
    /*--- dN/d xi ---*/
    
    dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);
    
    /*--- dN/d eta ---*/
    
    dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);
    
    /*--- dN/d mu ---*/
    
    dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1]*Jac_Ref[2][2]-Jac_Ref[1][2]*Jac_Ref[2][1];
    ad_Ref[0][1] = Jac_Ref[0][2]*Jac_Ref[2][1]-Jac_Ref[0][1]*Jac_Ref[2][2];
    ad_Ref[0][2] = Jac_Ref[0][1]*Jac_Ref[1][2]-Jac_Ref[0][2]*Jac_Ref[1][1];
    ad_Ref[1][0] = Jac_Ref[1][2]*Jac_Ref[2][0]-Jac_Ref[1][0]*Jac_Ref[2][2];
    ad_Ref[1][1] = Jac_Ref[0][0]*Jac_Ref[2][2]-Jac_Ref[0][2]*Jac_Ref[2][0];
    ad_Ref[1][2] = Jac_Ref[0][2]*Jac_Ref[1][0]-Jac_Ref[0][0]*Jac_Ref[1][2];
    ad_Ref[2][0] = Jac_Ref[1][0]*Jac_Ref[2][1]-Jac_Ref[1][1]*Jac_Ref[2][0];
    ad_Ref[2][1] = Jac_Ref[0][1]*Jac_Ref[2][0]-Jac_Ref[0][0]*Jac_Ref[2][1];
    ad_Ref[2][2] = Jac_Ref[0][0]*Jac_Ref[1][1]-Jac_Ref[0][1]*Jac_Ref[1][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1]*Jac_Curr[2][2]-Jac_Curr[1][2]*Jac_Curr[2][1];
    ad_Curr[0][1] = Jac_Curr[0][2]*Jac_Curr[2][1]-Jac_Curr[0][1]*Jac_Curr[2][2];
    ad_Curr[0][2] = Jac_Curr[0][1]*Jac_Curr[1][2]-Jac_Curr[0][2]*Jac_Curr[1][1];
    ad_Curr[1][0] = Jac_Curr[1][2]*Jac_Curr[2][0]-Jac_Curr[1][0]*Jac_Curr[2][2];
    ad_Curr[1][1] = Jac_Curr[0][0]*Jac_Curr[2][2]-Jac_Curr[0][2]*Jac_Curr[2][0];
    ad_Curr[1][2] = Jac_Curr[0][2]*Jac_Curr[1][0]-Jac_Curr[0][0]*Jac_Curr[1][2];
    ad_Curr[2][0] = Jac_Curr[1][0]*Jac_Curr[2][1]-Jac_Curr[1][1]*Jac_Curr[2][0];
    ad_Curr[2][1] = Jac_Curr[0][1]*Jac_Curr[2][0]-Jac_Curr[0][0]*Jac_Curr[2][1];
    ad_Curr[2][2] = Jac_Curr[0][0]*Jac_Curr[1][1]-Jac_Curr[0][1]*Jac_Curr[1][0];
    
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = Jac_Ref[0][0]*ad_Ref[0][0]+Jac_Ref[0][1]*ad_Ref[1][0]+Jac_Ref[0][2]*ad_Ref[2][0];
    detJac_Curr = Jac_Curr[0][0]*ad_Curr[0][0]+Jac_Curr[0][1]*ad_Curr[1][0]+Jac_Curr[0][2]*ad_Curr[2][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
  
}


CHEXA8P1::CHEXA8P1(void) : CHEXA8() {
  
}

CHEXA8P1::CHEXA8P1(unsigned short val_nDim, CConfig *config)
: CHEXA8(val_nDim, config) {
  
  unsigned short iNode, iGauss, jNode;
  unsigned short nDimSq;
  
  nGaussPointsP = 1;
  
  nDimSq = nDim*nDim;
  
  GaussPointP = new CGaussVariable*[nGaussPointsP];
  for (iGauss = 0; iGauss < nGaussPointsP; iGauss++) {
    GaussPointP[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  GaussWeightP = new su2double [nGaussPointsP];
  
  GaussCoordP = new su2double*[nGaussPointsP];
  for (iGauss = 0; iGauss < nGaussPointsP; iGauss++) {
    GaussCoordP [iGauss] = new su2double[nDim];
  }
  
  GaussCoordP[0][0] = 0.0;  GaussCoordP[0][1] = 0.0;  GaussCoordP[0][1] = 0.0;  GaussWeightP[0] = 8.0;
  
  Kk_ab = new su2double **[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    Kk_ab [iNode] = new su2double*[nNodes];
    for (jNode = 0; jNode < nNodes; jNode++) {
      Kk_ab [iNode][jNode] = new su2double[nDimSq];
    }
  }
  
}

CHEXA8P1::~CHEXA8P1(void) {
  
}

void CHEXA8P1::ComputeGrad_Pressure(void) {
  
  su2double Xi, Eta, Zeta;
  su2double Jac_Ref[3][3], Jac_Curr[3][3], dNiXj[8][3];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[3][3], ad_Curr[3][3];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPointsP; iGauss++) {
    
    Xi = GaussCoordP[iGauss][0];
    Eta = GaussCoordP[iGauss][1];
    Zeta = GaussCoordP[iGauss][2];
    
    /*--- dN/d xi, dN/d eta ---*/
    
    /*--- dN/d xi ---*/
    
    dNiXj[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
    dNiXj[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
    dNiXj[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
    dNiXj[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
    dNiXj[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);
    
    /*--- dN/d eta ---*/
    
    dNiXj[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
    dNiXj[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
    dNiXj[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
    dNiXj[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
    dNiXj[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);
    
    /*--- dN/d mu ---*/
    
    dNiXj[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
    dNiXj[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
    dNiXj[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
    dNiXj[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
    dNiXj[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1]*Jac_Ref[2][2]-Jac_Ref[1][2]*Jac_Ref[2][1];
    ad_Ref[0][1] = Jac_Ref[0][2]*Jac_Ref[2][1]-Jac_Ref[0][1]*Jac_Ref[2][2];
    ad_Ref[0][2] = Jac_Ref[0][1]*Jac_Ref[1][2]-Jac_Ref[0][2]*Jac_Ref[1][1];
    ad_Ref[1][0] = Jac_Ref[1][2]*Jac_Ref[2][0]-Jac_Ref[1][0]*Jac_Ref[2][2];
    ad_Ref[1][1] = Jac_Ref[0][0]*Jac_Ref[2][2]-Jac_Ref[0][2]*Jac_Ref[2][0];
    ad_Ref[1][2] = Jac_Ref[0][2]*Jac_Ref[1][0]-Jac_Ref[0][0]*Jac_Ref[1][2];
    ad_Ref[2][0] = Jac_Ref[1][0]*Jac_Ref[2][1]-Jac_Ref[1][1]*Jac_Ref[2][0];
    ad_Ref[2][1] = Jac_Ref[0][1]*Jac_Ref[2][0]-Jac_Ref[0][0]*Jac_Ref[2][1];
    ad_Ref[2][2] = Jac_Ref[0][0]*Jac_Ref[1][1]-Jac_Ref[0][1]*Jac_Ref[1][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1]*Jac_Curr[2][2]-Jac_Curr[1][2]*Jac_Curr[2][1];
    ad_Curr[0][1] = Jac_Curr[0][2]*Jac_Curr[2][1]-Jac_Curr[0][1]*Jac_Curr[2][2];
    ad_Curr[0][2] = Jac_Curr[0][1]*Jac_Curr[1][2]-Jac_Curr[0][2]*Jac_Curr[1][1];
    ad_Curr[1][0] = Jac_Curr[1][2]*Jac_Curr[2][0]-Jac_Curr[1][0]*Jac_Curr[2][2];
    ad_Curr[1][1] = Jac_Curr[0][0]*Jac_Curr[2][2]-Jac_Curr[0][2]*Jac_Curr[2][0];
    ad_Curr[1][2] = Jac_Curr[0][2]*Jac_Curr[1][0]-Jac_Curr[0][0]*Jac_Curr[1][2];
    ad_Curr[2][0] = Jac_Curr[1][0]*Jac_Curr[2][1]-Jac_Curr[1][1]*Jac_Curr[2][0];
    ad_Curr[2][1] = Jac_Curr[0][1]*Jac_Curr[2][0]-Jac_Curr[0][0]*Jac_Curr[2][1];
    ad_Curr[2][2] = Jac_Curr[0][0]*Jac_Curr[1][1]-Jac_Curr[0][1]*Jac_Curr[1][0];
    
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = Jac_Ref[0][0]*ad_Ref[0][0]+Jac_Ref[0][1]*ad_Ref[1][0]+Jac_Ref[0][2]*ad_Ref[2][0];
    detJac_Curr = Jac_Curr[0][0]*ad_Curr[0][0]+Jac_Curr[0][1]*ad_Curr[1][0]+Jac_Curr[0][2]*ad_Curr[2][0];
    
    GaussPointP[iGauss]->SetJ_X(detJac_Ref);
    GaussPointP[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPointP[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPointP[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
}


CBOUND2D::CBOUND2D(void) : CElement() {
  
}

CBOUND2D::CBOUND2D(unsigned short val_nDim, CConfig *config)
: CElement(val_nDim, config) {
  
  unsigned short iNode, iGauss;
  
  nNodes = 2;
  nGaussPoints = 2;
  
  GaussPoint = new CGaussVariable*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussPoint[iGauss] = new CGaussVariable(iGauss, nDim, nNodes);
  }
  
  /*--- Initialize structure for current and reference configuration ---*/
  
  CurrentCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    CurrentCoord [iNode] = new su2double[nDim];
  }
  
  RefCoord = new su2double*[nNodes];
  for (iNode = 0; iNode < nNodes; iNode++) {
    RefCoord [iNode] = new su2double[nDim];
  }
  
  GaussWeight = new su2double [nGaussPoints];
  
  GaussCoord = new su2double*[nGaussPoints];
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    GaussCoord [iGauss] = new su2double[1];
  }
  
  GaussCoord[0][0] = -0.577350269189626;  GaussWeight[0] = 1.0;
  GaussCoord[1][0] = 0.577350269189626;   GaussWeight[1] = 1.0;
  
  /*--- Store the shape functions (they only need to be computed once) ---*/
  su2double Xi, val_Ni;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    
    val_Ni = 0.5*(1.0-Xi);   GaussPoint[iGauss]->SetNi(val_Ni,0);
    val_Ni = 0.5*(1.0+Xi);   GaussPoint[iGauss]->SetNi(val_Ni,1);
    
  }
  
}

CBOUND2D::~CBOUND2D(void) {
  
}

void CBOUND2D::ComputeGrad_Linear(void) {
  
  su2double Jacobian[2][2], dNiXj[2][1];
  su2double detJac, GradNi_Xj;
  su2double ad[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi ---*/
    
    dNiXj[0][0] = -0.5;
    dNiXj[1][0] =  0.5;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jacobian[iDim][jDim] = Jacobian[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad[0][0] = Jacobian[1][1];
    ad[0][1] = -Jacobian[0][1];
    ad[1][0] = -Jacobian[1][0];
    ad[1][1] = Jacobian[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jacobian[iDim][jDim] = ad[iDim][jDim]/detJac;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj += Jacobian[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj, iDim, iNode);
      }
    }
  }
  
  
}

void CBOUND2D::ComputeGrad_NonLinear(void) {
  
  su2double Jac_Ref[2][2], Jac_Curr[2][2], dNiXj[4][2];
  su2double detJac_Ref, detJac_Curr, GradNi_Xj_Ref, GradNi_Xj_Curr;
  su2double ad_Ref[2][2], ad_Curr[2][2];
  unsigned short iNode, iDim, jDim, iGauss;
  
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    
    /*--- dN/d xi ---*/
    
    dNiXj[0][0] = -0.5;
    dNiXj[1][0] =  0.5;
    
    /*--- Jacobian transformation ---*/
    /*--- This does dX/dXi transpose ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Jac_Ref[iDim][jDim] = 0.0;
        Jac_Curr[iDim][jDim] = 0.0;
        for (iNode = 0; iNode < nNodes; iNode++) {
          Jac_Ref[iDim][jDim] = Jac_Ref[iDim][jDim]+RefCoord[iNode][jDim]*dNiXj[iNode][iDim];
          Jac_Curr[iDim][jDim] = Jac_Curr[iDim][jDim]+CurrentCoord[iNode][jDim]*dNiXj[iNode][iDim];
        }
      }
    }
    
    /*--- Adjoint to Jacobian ---*/
    
    ad_Ref[0][0] = Jac_Ref[1][1];
    ad_Ref[0][1] = -Jac_Ref[0][1];
    ad_Ref[1][0] = -Jac_Ref[1][0];
    ad_Ref[1][1] = Jac_Ref[0][0];
    
    ad_Curr[0][0] = Jac_Curr[1][1];
    ad_Curr[0][1] = -Jac_Curr[0][1];
    ad_Curr[1][0] = -Jac_Curr[1][0];
    ad_Curr[1][1] = Jac_Curr[0][0];
    
    /*--- Determinant of Jacobian ---*/
    
    detJac_Ref = ad_Ref[0][0]*ad_Ref[1][1]-ad_Ref[0][1]*ad_Ref[1][0];
    detJac_Curr = ad_Curr[0][0]*ad_Curr[1][1]-ad_Curr[0][1]*ad_Curr[1][0];
    
    GaussPoint[iGauss]->SetJ_X(detJac_Ref);
    GaussPoint[iGauss]->SetJ_x(detJac_Curr);
    
    /*--- Jacobian inverse (it was already computed as transpose) ---*/
    
    for (iDim = 0; iDim < 2; iDim++) {
      for (jDim = 0; jDim < 2; jDim++) {
        Jac_Ref[iDim][jDim] = ad_Ref[iDim][jDim]/detJac_Ref;
        Jac_Curr[iDim][jDim] = ad_Curr[iDim][jDim]/detJac_Curr;
      }
    }
    
    /*--- Derivatives with respect to global coordinates ---*/
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradNi_Xj_Ref = 0.0;
        GradNi_Xj_Curr = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          GradNi_Xj_Ref += Jac_Ref[iDim][jDim]*dNiXj[iNode][jDim];
          GradNi_Xj_Curr += Jac_Curr[iDim][jDim]*dNiXj[iNode][jDim];
        }
        GaussPoint[iGauss]->SetGradNi_Xj(GradNi_Xj_Ref, iDim, iNode);
        GaussPoint[iGauss]->SetGradNi_xj(GradNi_Xj_Curr, iDim, iNode);
      }
    }
  }
  
}



