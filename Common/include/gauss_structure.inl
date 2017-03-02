/*!
 * \file gauss_structure.inl
 * \brief In-Line subroutines of the <i>gauss_structure.hpp</i> file.
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
 
#pragma once

inline void CGaussVariable::SetGradNi_Xj(su2double val_GradNi_Xj, unsigned short val_iDim, unsigned short val_Ni) { GradNi_Xj[val_Ni][val_iDim] = val_GradNi_Xj; }

inline void CGaussVariable::SetGradNi_xj(su2double val_GradNi_xj, unsigned short val_iDim, unsigned short val_Ni) { GradNi_xj[val_Ni][val_iDim] = val_GradNi_xj; }

inline void CGaussVariable::SetNi(su2double val_ShapeNi, unsigned short val_Ni) { Ni[val_Ni] = val_ShapeNi; }

inline void CGaussVariable::SetJ_X(su2double valJ_X) { J_X = valJ_X; }

inline void CGaussVariable::SetJ_x(su2double valJ_x) { J_x = valJ_x; }

inline unsigned short CGaussVariable::Get_iGauss(void) { return iGaussPoint; }

inline su2double **CGaussVariable::GetGradNi_Xj(void) { return GradNi_Xj; }

inline su2double CGaussVariable::GetGradNi_Xj(unsigned short val_Ni, unsigned short val_iDim) { return GradNi_Xj[val_Ni][val_iDim]; }

inline su2double CGaussVariable::GetGradNi_xj(unsigned short val_Ni, unsigned short val_iDim) { return GradNi_xj[val_Ni][val_iDim]; }

inline su2double CGaussVariable::GetNi(unsigned short val_Ni) { return Ni[val_Ni]; }

inline su2double CGaussVariable::GetJ_X(void) { return J_X; }

inline su2double CGaussVariable::GetJ_x(void) { return J_x; }
