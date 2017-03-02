/*!
 * \file geometry_structure.cpp
 * \brief Main subroutines for creating the primal grid and multigrid structure.
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

#include "../include/geometry_structure.hpp"
#include "../include/adt_structure.hpp"

/*--- Epsilon definition ---*/

#define EPSILON 0.000001

/*--- Cross product ---*/

#define CROSS(dest,v1,v2) \
(dest)[0] = (v1)[1]*(v2)[2] - (v1)[2]*(v2)[1];	\
(dest)[1] = (v1)[2]*(v2)[0] - (v1)[0]*(v2)[2];	\
(dest)[2] = (v1)[0]*(v2)[1] - (v1)[1]*(v2)[0];

/*--- Cross product ---*/

#define DOT(v1,v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1] + (v1)[2]*(v2)[2]);

/*--- a = b - c ---*/

#define SUB(dest,v1,v2) \
(dest)[0] = (v1)[0] - (v2)[0];	\
(dest)[1] = (v1)[1] - (v2)[1];	\
(dest)[2] = (v1)[2] - (v2)[2];

CGeometry::CGeometry(void) {
  
  nEdge      = 0;
  nPoint     = 0;
  nPointNode = 0;
  nElem      = 0;
  
  nElem_Bound         = NULL;
  Tag_to_Marker       = NULL;
  elem                = NULL;
  face                = NULL;
  bound               = NULL;
  node                = NULL;
  edge                = NULL;
  vertex              = NULL;
  nVertex             = NULL;
  newBound            = NULL;
  nNewElem_Bound      = NULL;
  Marker_All_SendRecv = NULL;
  
  PeriodicPoint[MAX_NUMBER_PERIODIC][2].clear();
  PeriodicElem[MAX_NUMBER_PERIODIC].clear();

  XCoordList.clear();
  Xcoord_plane.clear();
  Ycoord_plane.clear();
  Zcoord_plane.clear();
  FaceArea_plane.clear();
  Plane_points.clear();
  
  /*--- Arrays for defining the linear partitioning ---*/
  
  starting_node = NULL;
  ending_node   = NULL;
  npoint_procs  = NULL;
  
}

CGeometry::~CGeometry(void) {
  
  unsigned long iElem, iElem_Bound, iEdge, iFace, iPoint, iVertex;
  unsigned short iMarker;
  
  if (elem != NULL) {
    for (iElem = 0; iElem < nElem; iElem++)
      if (elem[iElem] != NULL) delete elem[iElem];
    delete[] elem;
  }
  
  if (bound != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (bound[iMarker][iElem_Bound] != NULL) delete bound[iMarker][iElem_Bound];
      }
      if (bound[iMarker] != NULL) delete [] bound[iMarker];
    }
    delete [] bound;
  }
  
  if (face != NULL) {
    for (iFace = 0; iFace < nFace; iFace ++)
      if (face[iFace] != NULL) delete face[iFace];
    delete[] face;
  }
  
  if (node != NULL) {
    for (iPoint = 0; iPoint < nPointNode; iPoint ++)
      if (node[iPoint] != NULL) delete node[iPoint];
    delete[] node;
  }
  
  
  if (edge != NULL) {
    for (iEdge = 0; iEdge < nEdge; iEdge ++)
      if (edge[iEdge] != NULL) delete edge[iEdge];
    delete[] edge;
  }

  if (vertex != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        if (vertex[iMarker][iVertex] != NULL) delete vertex[iMarker][iVertex];
      }
      if (vertex[iMarker] != NULL) delete [] vertex[iMarker];
    }
    delete [] vertex;
  }
  
  if (newBound != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (newBound[iMarker][iElem_Bound] != NULL) delete [] newBound[iMarker][iElem_Bound];
      }
      delete[] newBound[iMarker];
    }
    delete[] newBound;
  }
  
  if (nElem_Bound         != NULL) delete [] nElem_Bound;
  if (nVertex             != NULL) delete [] nVertex;
  if (nNewElem_Bound      != NULL) delete [] nNewElem_Bound;
  if (Marker_All_SendRecv != NULL) delete [] Marker_All_SendRecv;
  if (Tag_to_Marker       != NULL) delete [] Tag_to_Marker;
  
  if (starting_node != NULL) delete [] starting_node;
  if (ending_node   != NULL) delete [] ending_node;
  if (npoint_procs  != NULL) delete [] npoint_procs;
  
}

su2double CGeometry::Point2Plane_Distance(su2double *Coord, su2double *iCoord, su2double *jCoord, su2double *kCoord) {
  su2double CrossProduct[3], iVector[3], jVector[3], distance, modulus;
  unsigned short iDim;
  
  for (iDim = 0; iDim < 3; iDim ++) {
    iVector[iDim] = jCoord[iDim] - iCoord[iDim];
    jVector[iDim] = kCoord[iDim] - iCoord[iDim];
  }
  
  CrossProduct[0] = iVector[1]*jVector[2] - iVector[2]*jVector[1];
  CrossProduct[1] = iVector[2]*jVector[0] - iVector[0]*jVector[2];
  CrossProduct[2] = iVector[0]*jVector[1] - iVector[1]*jVector[0];
  
  modulus = sqrt(CrossProduct[0]*CrossProduct[0]+CrossProduct[1]*CrossProduct[1]+CrossProduct[2]*CrossProduct[2]);
  
  distance = 0.0;
  for (iDim = 0; iDim < 3; iDim ++)
    distance += CrossProduct[iDim]*(Coord[iDim]-iCoord[iDim]);
  distance /= modulus;
  
  return distance;
  
}

long CGeometry::FindEdge(unsigned long first_point, unsigned long second_point) {
  unsigned long iPoint = 0;
  unsigned short iNode;
  for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
    iPoint = node[first_point]->GetPoint(iNode);
    if (iPoint == second_point) break;
  }
  
  if (iPoint == second_point) return node[first_point]->GetEdge(iNode);
  else {
    cout << "\n\n   !!! Error !!!\n" << endl;
    cout <<"Can't find the edge that connects "<< first_point <<" and "<< second_point <<"."<< endl;
    exit(EXIT_FAILURE);
    return -1;
  }
}

bool CGeometry::CheckEdge(unsigned long first_point, unsigned long second_point) {
  unsigned long iPoint = 0;
  unsigned short iNode;
  for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
    iPoint = node[first_point]->GetPoint(iNode);
    if (iPoint == second_point) break;
  }
  
  if (iPoint == second_point) return true;
  else return false;
  
}

void CGeometry::SetEdges(void) {
  unsigned long iPoint, jPoint;
  long iEdge;
  unsigned short jNode, iNode;
  long TestEdge = 0;
  
  nEdge = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
      jPoint = node[iPoint]->GetPoint(iNode);
      for (jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
        if (node[jPoint]->GetPoint(jNode) == iPoint) {
          TestEdge = node[jPoint]->GetEdge(jNode);
          break;
        }
      if (TestEdge == -1) {
        node[iPoint]->SetEdge(nEdge, iNode);
        node[jPoint]->SetEdge(nEdge, jNode);
        nEdge++;
      }
    }
  
  edge = new CEdge*[nEdge];
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
      jPoint = node[iPoint]->GetPoint(iNode);
      iEdge = FindEdge(iPoint, jPoint);
      if (iPoint < jPoint) edge[iEdge] = new CEdge(iPoint, jPoint, nDim);
    }
}

void CGeometry::SetFaces(void) {
  //	unsigned long iPoint, jPoint, iFace;
  //	unsigned short jNode, iNode;
  //	long TestFace = 0;
  //
  //	nFace = 0;
  //	for (iPoint = 0; iPoint < nPoint; iPoint++)
  //		for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
  //			jPoint = node[iPoint]->GetPoint(iNode);
  //			for (jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
  //				if (node[jPoint]->GetPoint(jNode) == iPoint) {
  //					TestFace = node[jPoint]->GetFace(jNode);
  //					break;
  //				}
  //			if (TestFace == -1) {
  //				node[iPoint]->SetFace(nFace, iNode);
  //				node[jPoint]->SetFace(nFace, jNode);
  //				nFace++;
  //			}
  //		}
  //
  //	face = new CFace*[nFace];
  //
  //	for (iPoint = 0; iPoint < nPoint; iPoint++)
  //		for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
  //			jPoint = node[iPoint]->GetPoint(iNode);
  //			iFace = FindFace(iPoint, jPoint);
  //			if (iPoint < jPoint) face[iFace] = new CFace(iPoint, jPoint, nDim);
  //		}
}

void CGeometry::TestGeometry(void) {
  
  ofstream para_file;
  
  para_file.open("test_geometry.dat", ios::out);
  
  su2double *Normal = new su2double[nDim];
  
  for (unsigned long iEdge = 0; iEdge < nEdge; iEdge++) {
    para_file << "Edge index: " << iEdge << endl;
    para_file << "   Point index: " << edge[iEdge]->GetNode(0) << "\t" << edge[iEdge]->GetNode(1) << endl;
    edge[iEdge]->GetNormal(Normal);
    para_file << "      Face normal : ";
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      para_file << Normal[iDim] << "\t";
    para_file << endl;
  }
  
  para_file << endl;
  para_file << endl;
  para_file << endl;
  para_file << endl;
  
  for (unsigned short iMarker =0; iMarker < nMarker; iMarker++) {
    para_file << "Marker index: " << iMarker << endl;
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      para_file << "   Vertex index: " << iVertex << endl;
      para_file << "      Point index: " << vertex[iMarker][iVertex]->GetNode() << endl;
      para_file << "      Point coordinates : ";
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        para_file << node[vertex[iMarker][iVertex]->GetNode()]->GetCoord(iDim) << "\t";}
      para_file << endl;
      vertex[iMarker][iVertex]->GetNormal(Normal);
      para_file << "         Face normal : ";
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        para_file << Normal[iDim] << "\t";
      para_file << endl;
    }
  }
  
  delete [] Normal;
  
}

void CGeometry::SetSpline(vector<su2double> &x, vector<su2double> &y, unsigned long n, su2double yp1, su2double ypn, vector<su2double> &y2) {
  unsigned long i, k;
  su2double p, qn, sig, un, *u;
  
  u = new su2double [n];
  
  if (yp1 > 0.99e30)			// The lower boundary condition is set either to be "nat
    y2[0]=u[0]=0.0;			  // -ural"
  else {									// or else to have a specified first derivative.
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  
  for (i=2; i<=n-1; i++) {									//  This is the decomposition loop of the tridiagonal al-
    sig=(x[i-1]-x[i-2])/(x[i]-x[i-2]);		//	gorithm. y2 and u are used for tem-
    p=sig*y2[i-2]+2.0;										//	porary storage of the decomposed
    y2[i-1]=(sig-1.0)/p;										//	factors.
    
    su2double a1 = (y[i]-y[i-1])/(x[i]-x[i-1]); if (x[i] == x[i-1]) a1 = 1.0;
    su2double a2 = (y[i-1]-y[i-2])/(x[i-1]-x[i-2]); if (x[i-1] == x[i-2]) a2 = 1.0;
    u[i-1]= a1 - a2;
    u[i-1]=(6.0*u[i-1]/(x[i]-x[i-2])-sig*u[i-2])/p;
    
  }
  
  if (ypn > 0.99e30)						// The upper boundary condition is set either to be
    qn=un=0.0;									// "natural"
  else {												// or else to have a specified first derivative.
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-1; k>=1; k--)					// This is the backsubstitution loop of the tridiagonal
    y2[k-1]=y2[k-1]*y2[k]+u[k-1];	  // algorithm.
  
  delete[] u;
  
}

su2double CGeometry::GetSpline(vector<su2double>&xa, vector<su2double>&ya, vector<su2double>&y2a, unsigned long n, su2double x) {
  unsigned long klo, khi, k;
  su2double h, b, a, y;
  
  if (x < xa[0]) x = xa[0];       // Clip max and min values
  if (x > xa[n-1]) x = xa[n-1];
  
  klo = 1;										// We will find the right place in the table by means of
  khi = n;										// bisection. This is optimal if sequential calls to this
  while (khi-klo > 1) {			// routine are at random values of x. If sequential calls
    k = (khi+klo) >> 1;				// are in order, and closely spaced, one would do better
    if (xa[k-1] > x) khi = k;		// to store previous values of klo and khi and test if
    else klo=k;							// they remain appropriate on the next call.
  }								// klo and khi now bracket the input value of x
  h = xa[khi-1] - xa[klo-1];
  if (h == 0.0) h = EPS; // cout << "Bad xa input to routine splint" << endl;	// The xa’s must be distinct.
  a = (xa[khi-1]-x)/h;
  b = (x-xa[klo-1])/h;				// Cubic spline polynomial is now evaluated.
  y = a*ya[klo-1]+b*ya[khi-1]+((a*a*a-a)*y2a[klo-1]+(b*b*b-b)*y2a[khi-1])*(h*h)/6.0;
  
  return y;
}

bool CGeometry::SegmentIntersectsPlane(su2double *Segment_P0, su2double *Segment_P1, su2double Variable_P0, su2double Variable_P1,
                                                           su2double *Plane_P0, su2double *Plane_Normal, su2double *Intersection, su2double &Variable_Interp) {
  su2double u[3], v[3], Denominator, Numerator, Aux, ModU;
  unsigned short iDim;
  
  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P1[iDim] - Segment_P0[iDim];
    v[iDim] = Plane_P0[iDim] - Segment_P0[iDim];
  }
  
  ModU = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
  
  Numerator = Plane_Normal[0]*v[0] + Plane_Normal[1]*v[1] + Plane_Normal[2]*v[2];
  Denominator = Plane_Normal[0]*u[0] + Plane_Normal[1]*u[1] + Plane_Normal[2]*u[2];
  
  if (fabs(Denominator) <= 0.0) return (false); // No intersection.
  
  Aux = Numerator / Denominator;
  
  if (Aux < 0.0 || Aux > 1.0) return (false); // No intersection.
  
  for (iDim = 0; iDim < 3; iDim++)
    Intersection[iDim] = Segment_P0[iDim] + Aux * u[iDim];
  
  
  /*--- Check that the intersection is in the segment ---*/
  
  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P0[iDim] - Intersection[iDim];
    v[iDim] = Segment_P1[iDim] - Intersection[iDim];
  }
  
  Variable_Interp = Variable_P0 + (Variable_P1 - Variable_P0)*sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])/ModU;
  
  Denominator = Plane_Normal[0]*u[0] + Plane_Normal[1]*u[1] + Plane_Normal[2]*u[2];
  Numerator = Plane_Normal[0]*v[0] + Plane_Normal[1]*v[1] + Plane_Normal[2]*v[2];
  
  Aux = Numerator * Denominator;
  
  if (Aux > 0.0) return (false); // Intersection outside the segment.
  
  return (true);
  
}

bool CGeometry::RayIntersectsTriangle(su2double orig[3], su2double dir[3],
                                      su2double vert0[3], su2double vert1[3], su2double vert2[3],
                                      su2double *intersect) {
  
  su2double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  su2double det, inv_det, t, u, v;
  
  /*--- Find vectors for two edges sharing vert0 ---*/
  
  SUB(edge1, vert1, vert0);
  SUB(edge2, vert2, vert0);
  
  /*--- Begin calculating determinant - also used to calculate U parameter ---*/
  
  CROSS(pvec, dir, edge2);
  
  /*--- If determinant is near zero, ray lies in plane of triangle ---*/
  
  det = DOT(edge1, pvec);
  
  
  if (det > -EPSILON && det < EPSILON) return(false);
  
  inv_det = 1.0 / det;
  
  /*--- Calculate distance from vert0 to ray origin ---*/
  
  SUB(tvec, orig, vert0);
  
  /*--- Calculate U parameter and test bounds ---*/
  
  u = inv_det * DOT(tvec, pvec);
  
  if (u < 0.0 || u > 1.0) return(false);
  
  /*--- prepare to test V parameter ---*/
  
  CROSS(qvec, tvec, edge1);
  
  /*--- Calculate V parameter and test bounds ---*/
  
  v = inv_det * DOT(dir, qvec);
  
  if (v < 0.0 || u + v > 1.0) return(false);
  
  /*--- Calculate t, ray intersects triangle ---*/
  
  t = inv_det * DOT(edge2, qvec);
  
  /*--- Compute the intersection point in cartesian coordinates ---*/
  
  intersect[0] = orig[0] + (t * dir[0]);
  intersect[1] = orig[1] + (t * dir[1]);
  intersect[2] = orig[2] + (t * dir[2]);

  return (true);
  
}

bool CGeometry::SegmentIntersectsLine(su2double point0[2], su2double point1[2], su2double vert0[2], su2double vert1[2]) {

  su2double det, diff0_A, diff0_B, diff1_A, diff1_B, intersect[2];

  diff0_A = point0[0] - point1[0];
  diff1_A = point0[1] - point1[1];

  diff0_B = vert0[0] - vert1[0];
  diff1_B = vert0[1] - vert1[1];

  det = (diff0_A)*(diff1_B) - (diff1_A)*(diff0_B);

  if (det == 0) return false;

  /*--- Compute point of intersection ---*/

  intersect[0] = ((point0[0]*point1[1] - point0[1]*point1[0])*diff0_B
                -(vert0[0]* vert1[1]  - vert0[1]* vert1[0])*diff0_A)/det;

  intersect[1] =  ((point0[0]*point1[1] - point0[1]*point1[0])*diff1_B
                  -(vert0[0]* vert1[1]  - vert0[1]* vert1[0])*diff1_A)/det;


  /*--- Check that the point is between the two surface points ---*/

  su2double dist0, dist1, length;

  dist0 = (intersect[0] - point0[0])*(intersect[0] - point0[0])
         +(intersect[1] - point0[1])*(intersect[1] - point0[1]);

  dist1 = (intersect[0] - point1[0])*(intersect[0] - point1[0])
         +(intersect[1] - point1[1])*(intersect[1] - point1[1]);

  length = diff0_A*diff0_A
          +diff1_A*diff1_A;

  if ( (dist0 > length) || (dist1 > length) ) {
    return false;
  }

  return true;
}

bool CGeometry::SegmentIntersectsTriangle(su2double point0[3], su2double point1[3],
                                          su2double vert0[3], su2double vert1[3], su2double vert2[3]) {
  
  su2double dir[3], intersect[3], u[3], v[3], edge1[3], edge2[3], Plane_Normal[3], Denominator, Numerator, Aux;
  
  SUB(dir, point1, point0);
  
  if (RayIntersectsTriangle(point0, dir, vert0, vert1, vert2, intersect)) {
    
    /*--- Check that the intersection is in the segment ---*/
    
    SUB(u, point0, intersect);
    SUB(v, point1, intersect);
    
    SUB(edge1, vert1, vert0);
    SUB(edge2, vert2, vert0);
    CROSS(Plane_Normal, edge1, edge2);
    
    Denominator = DOT(Plane_Normal, u);
    Numerator = DOT(Plane_Normal, v);
    
    Aux = Numerator * Denominator;
    
    /*--- Intersection outside the segment ---*/
    
    if (Aux > 0.0) return (false);
    
  }
  else {
    
    /*--- No intersection with the ray ---*/
    
    return (false);
    
  }
  
  /*--- Intersection inside the segment ---*/

  return (true);
  
}

void CGeometry::ComputeAirfoil_Section(su2double *Plane_P0, su2double *Plane_Normal,
                                       su2double MinXCoord, su2double MaxXCoord, su2double *FlowVariable,
                                       vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil,
                                       vector<su2double> &Zcoord_Airfoil, vector<su2double> &Variable_Airfoil,
                                       bool original_surface, CConfig *config) {
  
  unsigned short iMarker, iNode, jNode, iDim;
  bool intersect;
  long MinDist_Point, MinDistAngle_Point;
  unsigned long iPoint, jPoint, iElem, Trailing_Point, Airfoil_Point, iVertex, jVertex;
  su2double Segment_P0[3] = {0.0, 0.0, 0.0}, Segment_P1[3] = {0.0, 0.0, 0.0}, Variable_P0 = 0.0, Variable_P1 = 0.0, Intersection[3] = {0.0, 0.0, 0.0}, Trailing_Coord, MinDist_Value, MinDistAngle_Value, Dist_Value,
  Airfoil_Tangent[3] = {0.0, 0.0, 0.0}, Segment[3] = {0.0, 0.0, 0.0}, Length, Angle_Value, MaxAngle = 30, *VarCoord = NULL, CosValue, Variable_Interp;
  vector<su2double> Xcoord, Ycoord, Zcoord, Variable;
  vector<unsigned long> Duplicate;
  vector<unsigned long>::iterator it;
  int rank = MASTER_NODE;
  su2double **Coord_Variation = NULL;
  
#ifdef HAVE_MPI
  unsigned long nLocalVertex, nGlobalVertex, MaxLocalVertex, *Buffer_Send_nVertex, *Buffer_Receive_nVertex, nBuffer;
  int nProcessor, iProcessor;
  su2double *Buffer_Send_Coord, *Buffer_Receive_Coord;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  Xcoord_Airfoil.clear();
  Ycoord_Airfoil.clear();
  Zcoord_Airfoil.clear();
  Variable_Airfoil.clear();
  
  /*--- Set the right plane in 2D (note the change in Y-Z plane) ---*/
  
  if (nDim == 2) {
    Plane_P0[0] = 0.0;      Plane_P0[1] = 0.0;      Plane_P0[2] = 0.0;
    Plane_Normal[0] = 0.0;  Plane_Normal[1] = 1.0;  Plane_Normal[2] = 0.0;
  }
  
  /*--- the grid variation is stored using a vertices information,
   we should go from vertex to points ---*/
  
  if (original_surface == false) {
    
    Coord_Variation = new su2double *[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      Coord_Variation[iPoint] = new su2double [nDim];
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_GeoEval(iMarker) == YES) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          VarCoord = vertex[iMarker][iVertex]->GetVarCoord();
          iPoint = vertex[iMarker][iVertex]->GetNode();
          for (iDim = 0; iDim < nDim; iDim++)
            Coord_Variation[iPoint][iDim] = VarCoord[iDim];
        }
      }
    }
    
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_GeoEval(iMarker) == YES) {
      for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
        for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
          iPoint = bound[iMarker][iElem]->GetNode(iNode);
          for (jNode = 0; jNode < bound[iMarker][iElem]->GetnNodes(); jNode++) {
            jPoint = bound[iMarker][iElem]->GetNode(jNode);
            
            if ((jPoint > iPoint) && ((node[iPoint]->GetCoord(0) > MinXCoord) && (node[iPoint]->GetCoord(0) < MaxXCoord))) {
              
              Segment_P0[0] = 0.0;  Segment_P0[1] = 0.0;  Segment_P0[2] = 0.0;  Variable_P0 = 0.0;
              Segment_P1[0] = 0.0;  Segment_P1[1] = 0.0;  Segment_P1[2] = 0.0;  Variable_P1 = 0.0;
              
              for (iDim = 0; iDim < nDim; iDim++) {
                if (original_surface == true) {
                  Segment_P0[iDim] = node[iPoint]->GetCoord(iDim);
                  Segment_P1[iDim] = node[jPoint]->GetCoord(iDim);
                }
                else {
                  Segment_P0[iDim] = node[iPoint]->GetCoord(iDim) + Coord_Variation[iPoint][iDim];
                  Segment_P1[iDim] = node[jPoint]->GetCoord(iDim) + Coord_Variation[jPoint][iDim];
                }
              }
              
              if (FlowVariable != NULL) {
                Variable_P0 = FlowVariable[iPoint];
                Variable_P1 = FlowVariable[jPoint];
              }
              
              /*--- In 2D add the points directly (note the change between Y and Z coordinate) ---*/
              
              if (nDim == 2) {
                Xcoord.push_back(Segment_P0[0]);    Xcoord.push_back(Segment_P1[0]);
                Ycoord.push_back(Segment_P0[2]);    Ycoord.push_back(Segment_P1[2]);
                Zcoord.push_back(Segment_P0[1]);    Zcoord.push_back(Segment_P1[1]);
                Variable.push_back(Variable_P0);    Variable.push_back(Variable_P1);
              }
              /*--- In 3D compute the intersection ---*/
              
              else if (nDim == 3) {
                intersect = SegmentIntersectsPlane(Segment_P0, Segment_P1, Variable_P0, Variable_P1, Plane_P0, Plane_Normal, Intersection, Variable_Interp);
                if (intersect == true) {
                  Xcoord.push_back(Intersection[0]);
                  Ycoord.push_back(Intersection[1]);
                  Zcoord.push_back(Intersection[2]);
                  Variable.push_back(Variable_Interp);
                }
              }
              
            }
            
          }
        }
      }
    }
  }
  
  if (original_surface == false) {
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      delete [] Coord_Variation[iPoint];
    delete [] Coord_Variation;
    
  }
  
  
#ifdef HAVE_MPI
  
  /*--- Copy the coordinates of all the points in the plane to the master node ---*/
  
  nLocalVertex = 0, nGlobalVertex = 0, MaxLocalVertex = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  Buffer_Send_nVertex = new unsigned long [1];
  Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  nLocalVertex = Xcoord.size();
  
  Buffer_Send_nVertex[0] = nLocalVertex;
  
  SU2_MPI::Allreduce(&nLocalVertex, &nGlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  Buffer_Send_Coord = new su2double [MaxLocalVertex*4];
  Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex*4];
  nBuffer = MaxLocalVertex*4;
  
  for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
    Buffer_Send_Coord[iVertex*4 + 0] = Xcoord[iVertex];
    Buffer_Send_Coord[iVertex*4 + 1] = Ycoord[iVertex];
    Buffer_Send_Coord[iVertex*4 + 2] = Zcoord[iVertex];
    Buffer_Send_Coord[iVertex*4 + 3] = Variable[iVertex];
  }
  
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);
  
  /*--- Clean the vectors before adding the new vertices only to the master node ---*/
  
  Xcoord.clear();
  Ycoord.clear();
  Zcoord.clear();
  Variable.clear();
  
  /*--- Copy the boundary to the master node vectors ---*/
  if (rank == MASTER_NODE) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Xcoord.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalVertex*4 + iVertex*4 + 0] );
        Ycoord.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalVertex*4 + iVertex*4 + 1] );
        Zcoord.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalVertex*4 + iVertex*4 + 2] );
        Variable.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalVertex*4 + iVertex*4 + 3] );
      }
    }
  }
  
  delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_nVertex; delete[] Buffer_Receive_nVertex;
  
#endif
  
  if ((rank == MASTER_NODE) && (Xcoord.size() != 0)) {
    
    /*--- Create a list with the duplicated points ---*/
    
    for (iVertex = 0; iVertex < Xcoord.size()-1; iVertex++) {
      for (jVertex = iVertex+1; jVertex < Xcoord.size(); jVertex++) {
        Segment[0] = Xcoord[jVertex] - Xcoord[iVertex];
        Segment[1] = Ycoord[jVertex] - Ycoord[iVertex];
        Segment[2] = Zcoord[jVertex] - Zcoord[iVertex];
        Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
        if (Dist_Value < 1E-6) {
          Duplicate.push_back (jVertex);
        }
      }
    }
    
    sort(Duplicate.begin(), Duplicate.end());
    it = unique(Duplicate.begin(), Duplicate.end());
    Duplicate.resize(it - Duplicate.begin());
    
    /*--- Remove duplicated points (starting from the back) ---*/
    
    for (iVertex = Duplicate.size(); iVertex > 0; iVertex--) {
      Xcoord.erase (Xcoord.begin() + Duplicate[iVertex-1]);
      Ycoord.erase (Ycoord.begin() + Duplicate[iVertex-1]);
      Zcoord.erase (Zcoord.begin() + Duplicate[iVertex-1]);
      Variable.erase (Variable.begin() + Duplicate[iVertex-1]);
    }
    
    if (Xcoord.size() != 1) {
      
      /*--- Find the trailing edge ---*/
      
      Trailing_Point = 0; Trailing_Coord = Xcoord[0];
      for (iVertex = 1; iVertex < Xcoord.size(); iVertex++) {
        if (Xcoord[iVertex] > Trailing_Coord) {
          Trailing_Point = iVertex; Trailing_Coord = Xcoord[iVertex];
        }
      }
      
      /*--- Add the trailing edge to the list, and remove from the original list ---*/
      Xcoord_Airfoil.push_back(Xcoord[Trailing_Point]); Ycoord_Airfoil.push_back(Ycoord[Trailing_Point]); Zcoord_Airfoil.push_back(Zcoord[Trailing_Point]); Variable_Airfoil.push_back(Variable[Trailing_Point]);
      Xcoord.erase (Xcoord.begin() + Trailing_Point); Ycoord.erase (Ycoord.begin() + Trailing_Point); Zcoord.erase (Zcoord.begin() + Trailing_Point); Variable.erase (Variable.begin() + Trailing_Point);
      
      /*--- Find the next point using the right hand side rule ---*/
      MinDist_Value = 1E6; MinDist_Point = 0;
      for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
        Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[0];
        Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[0];
        Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[0];
        Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
        Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;
        
        if ((Dist_Value < MinDist_Value) && (Segment[2] > 0.0)) { MinDist_Point = iVertex; MinDist_Value = Dist_Value; }
      }
      
      Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);  Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);  Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);  Variable_Airfoil.push_back(Variable[MinDist_Point]);
      Xcoord.erase (Xcoord.begin() + MinDist_Point);    Ycoord.erase (Ycoord.begin() + MinDist_Point);    Zcoord.erase (Zcoord.begin() + MinDist_Point);    Variable.erase (Variable.begin() + MinDist_Point);
      
      /*--- Algorithm for the rest of the points ---*/
      do {
        
        /*--- Last added point in the list ---*/
        Airfoil_Point = Xcoord_Airfoil.size() - 1;
        
        /*--- Compute the slope of the curve ---*/
        Airfoil_Tangent[0] = Xcoord_Airfoil[Airfoil_Point] - Xcoord_Airfoil[Airfoil_Point-1];
        Airfoil_Tangent[1] = Ycoord_Airfoil[Airfoil_Point] - Ycoord_Airfoil[Airfoil_Point-1];
        Airfoil_Tangent[2] = Zcoord_Airfoil[Airfoil_Point] - Zcoord_Airfoil[Airfoil_Point-1];
        Length = sqrt(pow(Airfoil_Tangent[0], 2.0) + pow(Airfoil_Tangent[1], 2.0) + pow(Airfoil_Tangent[2], 2.0));
        Airfoil_Tangent[0] /= Length; Airfoil_Tangent[1] /= Length; Airfoil_Tangent[2] /= Length;
        
        /*--- Find the closest point with the right slope ---*/
        MinDist_Value = 1E6; MinDistAngle_Value = 180;
        MinDist_Point = -1; MinDistAngle_Point = -1;
        for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
          
          Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[Airfoil_Point];
          Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[Airfoil_Point];
          Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[Airfoil_Point];
          
          /*--- Compute the distance to each point ---*/
          Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
          
          /*--- Compute the angle of the point ---*/
          Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;
          
          /*--- Clip the value of the cosine, this is important due to the round errors ---*/
          CosValue = Airfoil_Tangent[0]*Segment[0] + Airfoil_Tangent[1]*Segment[1] + Airfoil_Tangent[2]*Segment[2];
          if (CosValue >= 1.0) CosValue = 1.0;
          if (CosValue <= -1.0) CosValue = -1.0;
          
          Angle_Value = acos(CosValue) * 180 / PI_NUMBER;
          
          if (Dist_Value < MinDist_Value) { MinDist_Point = iVertex; MinDist_Value = Dist_Value; }
          if ((Dist_Value < MinDistAngle_Value) && (Angle_Value < MaxAngle)) {MinDistAngle_Point = iVertex; MinDistAngle_Value = Dist_Value;}
          
        }
        
        if ( MinDistAngle_Point != -1) MinDist_Point = MinDistAngle_Point;
        
        /*--- Add and remove the min distance to the list ---*/
        Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);  Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);  Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);    Variable_Airfoil.push_back(Variable[MinDist_Point]);
        Xcoord.erase(Xcoord.begin() + MinDist_Point);     Ycoord.erase(Ycoord.begin() + MinDist_Point);     Zcoord.erase(Zcoord.begin() + MinDist_Point);       Variable.erase(Variable.begin() + MinDist_Point);
        
      } while (Xcoord.size() != 0);
      
      /*--- Clean the vector before using them again for storing the upper or the lower side ---*/
      
      Xcoord.clear(); Ycoord.clear(); Zcoord.clear(); Variable.clear();
      
    }
    
    
  }
  
}

void CGeometry::RegisterCoordinates(CConfig *config) {
  unsigned short iDim;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      AD::RegisterInput(node[iPoint]->GetCoord()[iDim]);
    }
  }
}

void CGeometry::UpdateGeometry(CGeometry **geometry_container, CConfig *config) {

    unsigned short iMesh;
    geometry_container[MESH_0]->Set_MPI_Coord(config);

    geometry_container[MESH_0]->SetCoord_CG();
    geometry_container[MESH_0]->SetControlVolume(config, UPDATE);
    geometry_container[MESH_0]->SetBoundControlVolume(config, UPDATE);

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
        /*--- Update the control volume structures ---*/

        geometry_container[iMesh]->SetControlVolume(config,geometry_container[iMesh-1], UPDATE);
        geometry_container[iMesh]->SetBoundControlVolume(config,geometry_container[iMesh-1], UPDATE);
        geometry_container[iMesh]->SetCoord(geometry_container[iMesh-1]);

    }

    if (config->GetKind_Solver() == DISC_ADJ_RANS)
      geometry_container[MESH_0]->ComputeWall_Distance(config);

}

void CGeometry::ComputeSurf_Curvature(CConfig *config) {
  unsigned short iMarker, iNeigh_Point, iDim, iNode, iNeighbor_Nodes, Neighbor_Node;
  unsigned long Neighbor_Point, iVertex, iPoint, jPoint, iElem_Bound, iEdge, nLocalVertex, MaxLocalVertex , *Buffer_Send_nVertex, *Buffer_Receive_nVertex, TotalnPointDomain;
  int iProcessor, nProcessor;
  vector<unsigned long> Point_NeighborList, Elem_NeighborList, Point_Triangle, Point_Edge, Point_Critical;
  vector<unsigned long>::iterator it;
  su2double U[3] = {0.0,0.0,0.0}, V[3] = {0.0,0.0,0.0}, W[3] = {0.0,0.0,0.0}, Length_U, Length_V, Length_W, CosValue, Angle_Value, *K, *Angle_Defect, *Area_Vertex, *Angle_Alpha, *Angle_Beta, **NormalMeanK, MeanK, GaussK, MaxPrinK, cot_alpha, cot_beta, delta, X1, X2, X3, Y1, Y2, Y3, radius, *Buffer_Send_Coord, *Buffer_Receive_Coord, *Coord, Dist, MinDist, MaxK, MinK, SigmaK;
  bool *Check_Edge;
  int rank;
  
#ifndef HAVE_MPI
  rank = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Allocate surface curvature ---*/
  K = new su2double [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) K[iPoint] = 0.0;
  
  if (nDim == 2) {
    
    /*--- Loop over all the markers ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        
        /*--- Loop through all marker vertices again, this time also
         finding the neighbors of each node.---*/
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint  = vertex[iMarker][iVertex]->GetNode();
          
          if (node[iPoint]->GetDomain()) {
            /*--- Loop through neighbors. In 2-D, there should be 2 nodes on either
             side of this vertex that lie on the same surface. ---*/
            Point_Edge.clear();
            
            for (iNeigh_Point = 0; iNeigh_Point < node[iPoint]->GetnPoint(); iNeigh_Point++) {
              Neighbor_Point = node[iPoint]->GetPoint(iNeigh_Point);
              
              /*--- Check if this neighbor lies on the surface. If so,
               add to the list of neighbors. ---*/
              if (node[Neighbor_Point]->GetPhysicalBoundary()) {
                Point_Edge.push_back(Neighbor_Point);
              }
              
            }
            
            if (Point_Edge.size() == 2) {
              
              /*--- Compute the curvature using three points ---*/
              X1 = node[iPoint]->GetCoord(0);
              X2 = node[Point_Edge[0]]->GetCoord(0);
              X3 = node[Point_Edge[1]]->GetCoord(0);
              Y1 = node[iPoint]->GetCoord(1);
              Y2 = node[Point_Edge[0]]->GetCoord(1);
              Y3 = node[Point_Edge[1]]->GetCoord(1);
              
              radius = sqrt(((X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1))*
                            ((X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3))*
                            ((X3-X1)*(X3-X1) + (Y3-Y1)*(Y3-Y1)))/
              (2.0*fabs(X1*Y2+X2*Y3+X3*Y1-X1*Y3-X2*Y1-X3*Y2)+EPS);
              
              K[iPoint] = 1.0/radius;
              node[iPoint]->SetCurvature(K[iPoint]);
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  else {
    
    Angle_Defect = new su2double [nPoint];
    Area_Vertex = new su2double [nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      Angle_Defect[iPoint] = 2*PI_NUMBER;
      Area_Vertex[iPoint] = 0.0;
    }
    
    Angle_Alpha = new su2double [nEdge];
    Angle_Beta = new su2double [nEdge];
    Check_Edge = new bool [nEdge];
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
      Angle_Alpha[iEdge] = 0.0;
      Angle_Beta[iEdge] = 0.0;
      Check_Edge[iEdge] = true;
    }
    
    NormalMeanK = new su2double *[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      NormalMeanK[iPoint] = new su2double [nDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        NormalMeanK[iPoint][iDim] = 0.0;
      }
    }
    
    /*--- Loop over all the markers ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        
        /*--- Loop over all the boundary elements ---*/
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          
          /*--- Only triangles ---*/
          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE) {
            
            /*--- Loop over all the nodes of the boundary element ---*/
            for (iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++) {
              
              iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);
              
              Point_Triangle.clear();
              
              for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                Neighbor_Point = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                Point_Triangle.push_back(Neighbor_Point);
              }
              
              iEdge = FindEdge(Point_Triangle[0], Point_Triangle[1]);
              
              for (iDim = 0; iDim < nDim; iDim++) {
                U[iDim] = node[Point_Triangle[0]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
                V[iDim] = node[Point_Triangle[1]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
              }
              
              W[0] = 0.5*(U[1]*V[2]-U[2]*V[1]); W[1] = -0.5*(U[0]*V[2]-U[2]*V[0]); W[2] = 0.5*(U[0]*V[1]-U[1]*V[0]);
              
              Length_U = 0.0, Length_V = 0.0, Length_W = 0.0, CosValue = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) { Length_U += U[iDim]*U[iDim]; Length_V += V[iDim]*V[iDim]; Length_W += W[iDim]*W[iDim]; }
              Length_U = sqrt(Length_U); Length_V = sqrt(Length_V); Length_W = sqrt(Length_W);
              for (iDim = 0; iDim < nDim; iDim++) { U[iDim] /= Length_U; V[iDim] /= Length_V; CosValue += U[iDim]*V[iDim]; }
              if (CosValue >= 1.0) CosValue = 1.0;
              if (CosValue <= -1.0) CosValue = -1.0;
              
              Angle_Value = acos(CosValue);
              Area_Vertex[iPoint] += Length_W;
              Angle_Defect[iPoint] -= Angle_Value;
              if (Angle_Alpha[iEdge] == 0.0) Angle_Alpha[iEdge] = Angle_Value;
              else Angle_Beta[iEdge] = Angle_Value;
              
            }
          }
        }
      }
    }
    
    /*--- Compute mean curvature ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE) {
            for (iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++) {
              iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);
              
              for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                jPoint = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                
                iEdge = FindEdge(iPoint, jPoint);
                
                if (Check_Edge[iEdge]) {
                  
                  Check_Edge[iEdge] = false;
                  
                  if (tan(Angle_Alpha[iEdge]) != 0.0) cot_alpha = 1.0/tan(Angle_Alpha[iEdge]); else cot_alpha = 0.0;
                  if (tan(Angle_Beta[iEdge]) != 0.0) cot_beta = 1.0/tan(Angle_Beta[iEdge]); else cot_beta = 0.0;
                  
                  /*--- iPoint, and jPoint ---*/
                  for (iDim = 0; iDim < nDim; iDim++) {
                    if (Area_Vertex[iPoint] != 0.0) NormalMeanK[iPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[iPoint]->GetCoord(iDim) - node[jPoint]->GetCoord(iDim)) / Area_Vertex[iPoint];
                    if (Area_Vertex[jPoint] != 0.0) NormalMeanK[jPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[jPoint]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim)) / Area_Vertex[jPoint];
                  }
                }
                
              }
            }
          }
        }
      }
    }
    
    /*--- Compute Gauss, mean, max and min principal curvature,
     and set the list of critical points ---*/
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint  = vertex[iMarker][iVertex]->GetNode();
          
          if (node[iPoint]->GetDomain()) {
            
            if (Area_Vertex[iPoint] != 0.0) GaussK = 3.0*Angle_Defect[iPoint]/Area_Vertex[iPoint];
            else GaussK = 0.0;
            
            MeanK = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              MeanK += NormalMeanK[iPoint][iDim]*NormalMeanK[iPoint][iDim];
            MeanK = sqrt(MeanK);
            
            delta = max((MeanK*MeanK - GaussK), 0.0);
            
            MaxPrinK = MeanK + sqrt(delta);
            
            /*--- Store the curvature value ---*/
            K[iPoint] = MaxPrinK;
            node[iPoint]->SetCurvature(K[iPoint]);
          }
          
        }
      }
    }
    
    delete [] Angle_Defect;
    delete [] Area_Vertex;
    delete [] Angle_Alpha;
    delete [] Angle_Beta;
    delete [] Check_Edge;
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      delete [] NormalMeanK[iPoint];
    delete [] NormalMeanK;
    
  }
  
  /*--- Sharp edge detection is based in the statistical
   distribution of the curvature ---*/
  
  MaxK = K[0]; MinK = K[0]; MeanK = 0.0; TotalnPointDomain = 0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          MaxK = max(MaxK, fabs(K[iPoint]));
          MinK = min(MinK, fabs(K[iPoint]));
          MeanK += fabs(K[iPoint]);
          TotalnPointDomain++;
        }
      }
    }
  }
  
#ifdef HAVE_MPI
  su2double MyMeanK = MeanK; MeanK = 0.0;
  su2double MyMaxK = MaxK; MaxK = 0.0;
  unsigned long MynPointDomain = TotalnPointDomain; TotalnPointDomain = 0;
  SU2_MPI::Allreduce(&MyMeanK, &MeanK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyMaxK, &MaxK, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MynPointDomain, &TotalnPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Compute the mean ---*/
  MeanK /= su2double(TotalnPointDomain);
  
  /*--- Compute the standard deviation ---*/
  SigmaK = 0.0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          SigmaK += (fabs(K[iPoint]) - MeanK) * (fabs(K[iPoint]) - MeanK);
        }
      }
    }
  }
  
#ifdef HAVE_MPI
  su2double MySigmaK = SigmaK; SigmaK = 0.0;
  SU2_MPI::Allreduce(&MySigmaK, &SigmaK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  SigmaK = sqrt(SigmaK/su2double(TotalnPointDomain));
  
  if (rank == MASTER_NODE)
    cout << "Max K: " << MaxK << ". Mean K: " << MeanK << ". Standard deviation K: " << SigmaK << "." << endl;
  
  Point_Critical.clear();
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          if (fabs(K[iPoint]) > MeanK + config->GetRefSharpEdges()*SigmaK) {
            Point_Critical.push_back(iPoint);
          }
        }
      }
    }
  }
  
  /*--- Variables and buffers needed for MPI ---*/
  
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#else
  nProcessor = 1;
#endif
  
  Buffer_Send_nVertex    = new unsigned long [1];
  Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  /*--- Count the total number of critical edge nodes. ---*/
  
  nLocalVertex = Point_Critical.size();
  Buffer_Send_nVertex[0] = nLocalVertex;
  
  /*--- Communicate to all processors the total number of critical edge nodes. ---*/
  
#ifdef HAVE_MPI
  MaxLocalVertex = 0;
  SU2_MPI::Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  MaxLocalVertex = nLocalVertex;
  Buffer_Receive_nVertex[0] = nLocalVertex;
#endif
  
  
  /*--- Create and initialize to zero some buffers to hold the coordinates
   of the boundary nodes that are communicated from each partition (all-to-all). ---*/
  
  Buffer_Send_Coord     = new su2double [MaxLocalVertex*nDim];
  Buffer_Receive_Coord  = new su2double [nProcessor*MaxLocalVertex*nDim];
  
#ifdef HAVE_MPI
  unsigned long nBuffer               = MaxLocalVertex*nDim;
#endif

  for (iVertex = 0; iVertex < MaxLocalVertex; iVertex++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
    }
  }
  
  /*--- Retrieve and store the coordinates of the sharp edges boundary nodes on
   the local partition and broadcast them to all partitions. ---*/
  
  for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
    iPoint = Point_Critical[iVertex];
    for (iDim = 0; iDim < nDim; iDim++)
      Buffer_Send_Coord[iVertex*nDim+iDim] = node[iPoint]->GetCoord(iDim);
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);
#else
  for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Buffer_Receive_Coord[iVertex*nDim+iDim] = Buffer_Send_Coord[iVertex*nDim+iDim];
    }
  }
#endif
  
  /*--- Loop over all interior mesh nodes on the local partition and compute
   the distances to each of the no-slip boundary nodes in the entire mesh.
   Store the minimum distance to the wall for each interior mesh node. ---*/
  
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    Coord = node[iPoint]->GetCoord();
    
    MinDist = 1E20;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Dist = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Dist += (Coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex+iVertex)*nDim+iDim])*
          (Coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex+iVertex)*nDim+iDim]);
        }
        if (Dist!=0.0) Dist = sqrt(Dist);
        else Dist = 0.0;
        if (Dist < MinDist) MinDist = Dist;
      }
    }
    node[iPoint]->SetSharpEdge_Distance(MinDist);
  }
  
  /*--- Deallocate Max curvature ---*/
  delete[] K;
  
  /*--- Deallocate the buffers needed for the MPI communication. ---*/
  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_nVertex;
  delete[] Buffer_Receive_nVertex;
  
}

CPhysicalGeometry::CPhysicalGeometry() : CGeometry() {

  Local_to_Global_Point  = NULL;
  Local_to_Global_Marker = NULL;
  Global_to_Local_Marker = NULL;

  starting_node = NULL;
  ending_node   = NULL;
  npoint_procs  = NULL;
}

CPhysicalGeometry::CPhysicalGeometry(CConfig *config, unsigned short val_iZone, unsigned short val_nZone) : CGeometry() {
  
  Local_to_Global_Point = NULL;
  Local_to_Global_Marker = NULL;
  Global_to_Local_Marker = NULL;
  
  starting_node = NULL;
  ending_node   = NULL;
  npoint_procs  = NULL;
  
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short iDim, iMarker, iNodes;
  unsigned long iPoint, iElem_Bound;
  su2double *NewCoord;
  nZone = val_nZone;
  ofstream boundary_file;
  string Grid_Marker;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  string val_mesh_filename  = config->GetMesh_FileName();
  unsigned short val_format = config->GetMesh_FileFormat();

  /*--- Initialize counters for local/global points & elements ---*/
  
  if (rank == MASTER_NODE)
    cout << endl <<"---------------------- Read Grid File Information -----------------------" << endl;
  
  switch (val_format) {
    case SU2:
      Read_SU2_Format_Parallel(config, val_mesh_filename, val_iZone, val_nZone);
      break;
    case CGNS:
      Read_CGNS_Format_Parallel(config, val_mesh_filename, val_iZone, val_nZone);
      break;
    default:
      if (rank == MASTER_NODE) cout << "Unrecognized mesh format specified!" << endl;
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
      break;
  }

  /*--- After reading the mesh, assert that the dimension is equal to 2 or 3. ---*/
  
  assert((nDim == 2) || (nDim == 3));
  
  /*--- Loop over the points element to re-scale the mesh, and plot it (only SU2_CFD) ---*/
  
  if (config->GetKind_SU2() == SU2_CFD) {
    
    NewCoord = new su2double [nDim];
    
    /*--- The US system uses feet, but SU2 assumes that the grid is in inches ---*/
    
    if (config->GetSystemMeasurements() == US) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          NewCoord[iDim] = node[iPoint]->GetCoord(iDim)/12.0;
        }
        node[iPoint]->SetCoord(NewCoord);
      }
    }
    
    delete [] NewCoord;
    
  }
  
  /*--- If SU2_DEF then write a file with the boundary information ---*/
  
  if ((config->GetKind_SU2() == SU2_DEF) && (rank == MASTER_NODE)) {

    /*--- Open .su2 grid file ---*/
    
    boundary_file.open("boundary.su2", ios::out);
    
    /*--- Loop through and write the boundary info ---*/
    
    boundary_file << "NMARK= " << nMarker << endl;
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      
      Grid_Marker = config->GetMarker_All_TagBound(iMarker);
      boundary_file << "MARKER_TAG= " << Grid_Marker << endl;
      boundary_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      
      if (nDim == 2) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          boundary_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            boundary_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;

          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == VERTEX) {
            boundary_file << bound[iMarker][iElem_Bound]->GetRotation_Type() << "\t";
          }
          boundary_file	<< iElem_Bound << endl;
        }
      }
      
      if (nDim == 3) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          boundary_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            boundary_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;

          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == VERTEX) {
            boundary_file << bound[iMarker][iElem_Bound]->GetRotation_Type() << "\t";
          }
          boundary_file	<< iElem_Bound << endl;
        }
      }
      
    }
    
    boundary_file.close();

  }
  
}

CPhysicalGeometry::CPhysicalGeometry(CGeometry *geometry, CConfig *config) {
  
  /*--- Initialize several class data members for later. ---*/
  
  Local_to_Global_Point  = NULL;
  Local_to_Global_Marker = NULL;
  Global_to_Local_Marker = NULL;
  
  starting_node = NULL;
  ending_node   = NULL;
  npoint_procs  = NULL;

  /*--- Local variables and counters for the following communications. ---*/
  
  unsigned long iter,  iPoint, jPoint, iElem, iVertex;
  unsigned long iElemTotal, iPointTotal, iPointGhost, iPointDomain, iPointPeriodic, iElemTriangle, iElemQuadrilateral, iElemTetrahedron, iElemHexahedron, iElemPrism, iElemPyramid, iPointCurrent;
  unsigned long nBoundLineTotal = 0, iBoundLineTotal;
  unsigned long nBoundTriangleTotal = 0, iBoundTriangleTotal;
  unsigned long nBoundQuadrilateralTotal = 0, iBoundQuadrilateralTotal;
  unsigned long ReceptorColor = 0, DonorColor = 0, Transformation;
  unsigned long nTotalSendDomain_Periodic = 0, iTotalSendDomain_Periodic, nTotalReceivedDomain_Periodic = 0, iTotalReceivedDomain_Periodic, *nSendDomain_Periodic = NULL, *nReceivedDomain_Periodic = NULL;
  unsigned long Buffer_Send_nPointTotal = 0, Buffer_Send_nPointDomainTotal = 0, Buffer_Send_nPointGhost = 0, Buffer_Send_nPointPeriodic = 0;
  unsigned long Buffer_Send_nElemTotal, Buffer_Send_nElemTriangle = 0, Buffer_Send_nElemQuadrilateral = 0, Buffer_Send_nElemTetrahedron = 0, Buffer_Send_nElemHexahedron = 0, Buffer_Send_nElemPrism = 0, Buffer_Send_nElemPyramid = 0;
  unsigned long Buffer_Send_nTotalSendDomain_Periodic = 0, Buffer_Send_nTotalReceivedDomain_Periodic = 0, *Buffer_Send_nSendDomain_Periodic = NULL, *Buffer_Send_nReceivedDomain_Periodic = NULL;
  unsigned long Buffer_Send_nBoundLineTotal = 0, Buffer_Send_nBoundTriangleTotal = 0, Buffer_Send_nBoundQuadrilateralTotal = 0;
  unsigned long iVertexDomain, iBoundLine, iBoundTriangle, iBoundQuadrilateral;
  
  unsigned long iNode, iDim, iMarker, jMarker, nMarkerDomain = 0, iMarkerDomain;
  unsigned long nDomain = 0, iDomain, jDomain, nPeriodic = 0, iPeriodic, Buffer_Send_nMarkerDomain = 0, Buffer_Send_nDim = 0, Buffer_Send_nZone = 0, Buffer_Send_nPeriodic = 0;
  
  bool *MarkerIn = NULL, **VertexIn = NULL, *ElemIn = NULL;
  long vnodes_local[8];
  
  vector<long> DomainList;
  short *Marker_All_SendRecv_Copy = NULL;
  string *Marker_All_TagBound_Copy = NULL;
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  unsigned short nMarker_Max = config->GetnMarker_Max();
  
   
  /*--- Some dynamic arrays so we're not allocating too much on the stack ---*/
  
  unsigned long *nVertexDomain       = new unsigned long[nMarker_Max];
  unsigned long *nBoundLine          = new unsigned long[nMarker_Max];
  unsigned long *nBoundTriangle      = new unsigned long[nMarker_Max];
  unsigned long *nBoundQuadrilateral = new unsigned long[nMarker_Max];
  
  unsigned long *Buffer_Send_nVertexDomain       = new unsigned long[nMarker_Max];
  unsigned long *Buffer_Send_nBoundLine          = new unsigned long[nMarker_Max];
  unsigned long *Buffer_Send_nBoundTriangle      = new unsigned long[nMarker_Max];
  unsigned long *Buffer_Send_nBoundQuadrilateral = new unsigned long[nMarker_Max];
  
  short *Buffer_Send_Marker_All_SendRecv = new short[nMarker_Max];
  
  char *Marker_All_TagBound             = new char[nMarker_Max*MAX_STRING_SIZE];
  char *Buffer_Send_Marker_All_TagBound = new char[nMarker_Max*MAX_STRING_SIZE];
  
#ifdef HAVE_MPI
  
  /*--- MPI initialization ---*/
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  MPI_Status status, status2;
  unsigned long source;
  int recv_count=0;
  
  int offset = 17;
  MPI_Status *send_stat = new MPI_Status[offset+size];
  MPI_Status *recv_stat = new MPI_Status[offset+size];
  
  MPI_Request *send_req = new MPI_Request[offset+size];
  MPI_Request *recv_req = new MPI_Request[offset+size];
  
#endif
  
  if (rank == MASTER_NODE && size > SINGLE_NODE)
    cout << "Communicating partition data and creating halo layers." << endl;
  
  /*--- Define buffer vectors for the interior points / elements ---*/
  
  su2double *Buffer_Send_Coord = NULL;
  
  unsigned long *Buffer_Send_Color            = NULL;
  unsigned long *Buffer_Send_GlobalPointIndex = NULL;
  unsigned long *Buffer_Send_Triangle         = NULL;
  unsigned long *Buffer_Send_Quadrilateral    = NULL;
  unsigned long *Buffer_Send_Tetrahedron      = NULL;
  unsigned long *Buffer_Send_Hexahedron       = NULL;
  unsigned long *Buffer_Send_Prism            = NULL;
  unsigned long *Buffer_Send_Pyramid          = NULL;
  unsigned long *Buffer_Send_GlobElem         = NULL;
  
  /*--- Define buffer vectors for boundary information ---*/
  
  unsigned long *Buffer_Send_BoundLine = NULL,           *Buffer_Receive_BoundLine = NULL;
  unsigned long *Buffer_Send_BoundTriangle = NULL,       *Buffer_Receive_BoundTriangle = NULL;
  unsigned long *Buffer_Send_BoundQuadrilateral = NULL,  *Buffer_Receive_BoundQuadrilateral = NULL;
  unsigned long *Buffer_Send_Local2Global_Marker = NULL, *Buffer_Receive_Local2Global_Marker = NULL;
  
  /*--- Define buffer vectors for periodic boundary conditions ---*/
  
  su2double *Buffer_Send_Center = NULL,    *Buffer_Receive_Center = NULL;
  su2double *Buffer_Send_Rotation = NULL,  *Buffer_Receive_Rotation = NULL;
  su2double *Buffer_Send_Translate = NULL, *Buffer_Receive_Translate = NULL;
  
  /*--- Define buffer vector periodic boundary conditions ---*/
  
  unsigned long *Buffer_Send_SendDomain_Periodic = NULL,          *Buffer_Receive_SendDomain_Periodic = NULL;
  unsigned long *Buffer_Send_SendDomain_PeriodicTrans = NULL,     *Buffer_Receive_SendDomain_PeriodicTrans = NULL;
  unsigned long *Buffer_Send_SendDomain_PeriodicReceptor = NULL,  *Buffer_Receive_SendDomain_PeriodicReceptor = NULL;
  unsigned long *Buffer_Send_ReceivedDomain_Periodic = NULL,      *Buffer_Receive_ReceivedDomain_Periodic = NULL;
  unsigned long *Buffer_Send_ReceivedDomain_PeriodicTrans = NULL, *Buffer_Receive_ReceivedDomain_PeriodicTrans = NULL;
  unsigned long *Buffer_Send_ReceivedDomain_PeriodicDonor = NULL, *Buffer_Receive_ReceivedDomain_PeriodicDonor = NULL;
  
  /*--- Variables below are needed specifically for the ParMETIS version ---*/
  
  unsigned short *nDim_s  = new unsigned short[size];
  unsigned short *nDim_r  = new unsigned short[size];
  unsigned short *nZone_s = new unsigned short[size];
  unsigned short *nZone_r = new unsigned short[size];
  
  unsigned long *nPointTotal_s        = new unsigned long[size];
  unsigned long *nPointDomainTotal_s  = new unsigned long[size];
  unsigned long *nPointGhost_s        = new unsigned long[size];
  unsigned long *nPointPeriodic_s     = new unsigned long[size];
  unsigned long *nElemTotal_s         = new unsigned long[size];
  unsigned long *nElemTriangle_s      = new unsigned long[size];
  unsigned long *nElemQuadrilateral_s = new unsigned long[size];
  unsigned long *nElemTetrahedron_s   = new unsigned long[size];
  unsigned long *nElemHexahedron_s    = new unsigned long[size];
  unsigned long *nElemPrism_s         = new unsigned long[size];
  unsigned long *nElemPyramid_s       = new unsigned long[size];
  
  unsigned long *nPointTotal_r        = new unsigned long[size];
  unsigned long *nPointDomainTotal_r  = new unsigned long[size];
  unsigned long *nPointGhost_r        = new unsigned long[size];
  unsigned long *nPointPeriodic_r     = new unsigned long[size];
  unsigned long *nElemTotal_r         = new unsigned long[size];
  unsigned long *nElemTriangle_r      = new unsigned long[size];
  unsigned long *nElemQuadrilateral_r = new unsigned long[size];
  unsigned long *nElemTetrahedron_r   = new unsigned long[size];
  unsigned long *nElemHexahedron_r    = new unsigned long[size];
  unsigned long *nElemPrism_r         = new unsigned long[size];
  unsigned long *nElemPyramid_r       = new unsigned long[size];
  
  /*--- Counters needed to track numbers of points, elements, etc. ---*/
  
  unsigned long nPointTotal_r_tot=0;
  unsigned long nPointDomainTotal_r_tot=0;
  unsigned long nPointGhost_r_tot=0;
  unsigned long nPointPeriodic_r_tot=0;
  unsigned long nElemTotal_r_tot=0;
  unsigned long nElemTriangle_r_tot=0;
  unsigned long nElemQuadrilateral_r_tot=0;
  unsigned long nElemTetrahedron_r_tot=0;
  unsigned long nElemHexahedron_r_tot=0;
  unsigned long nElemPrism_r_tot=0;
  unsigned long nElemPyramid_r_tot=0;
  
  unsigned long Buffer_Size_Coord = 0;
  unsigned long Buffer_Size_Color = 0;
  unsigned long Buffer_Size_GlobalPointIndex = 0;
  unsigned long Buffer_Size_Triangle = 0;
  unsigned long Buffer_Size_Quadrilateral = 0;
  unsigned long Buffer_Size_Tetrahedron = 0;
  unsigned long Buffer_Size_Hexahedron = 0;
  unsigned long Buffer_Size_Prism = 0;
  unsigned long Buffer_Size_Pyramid = 0;
  unsigned long Buffer_Size_GlobElem = 0;
  
  unsigned long ElemTotal_Counter = 0;
  unsigned long PointTotal_Counter = 0;
  unsigned long PointDomain_Counter = 0;
  
  /*--- WARNING: check the next two counters ---*/
  unsigned long PointPeriodic_Counter = 0;
  unsigned long PointGhost_Counter = 0;
  unsigned long ElemTriangle_Counter = 0;
  unsigned long ElemQuadrilateral_Counter = 0;
  unsigned long ElemTetrahedron_Counter = 0;
  unsigned long ElemHexahedron_Counter = 0;
  unsigned long ElemPrism_Counter = 0;
  unsigned long ElemPyramid_Counter = 0;
  
  unsigned long *Local_to_global_Triangle      = NULL;
  unsigned long *Local_to_global_Quadrilateral = NULL;
  unsigned long *Local_to_global_Tetrahedron   = NULL;
  unsigned long *Local_to_global_Hexahedron    = NULL;
  unsigned long *Local_to_global_Prism         = NULL;
  unsigned long *Local_to_global_Pyramid       = NULL;
  
  map<unsigned long,bool> Triangle_presence;
  map<unsigned long,bool> Quadrilateral_presence;
  map<unsigned long,bool> Tetrahedron_presence;
  map<unsigned long,bool> Hexahedron_presence;
  map<unsigned long,bool> Prism_presence;
  map<unsigned long,bool> Pyramid_presence;
  
  su2double *Buffer_Receive_Coord_loc = NULL;
  
  unsigned long *Buffer_Receive_Color_loc            = NULL;
  unsigned long *Buffer_Receive_GlobalPointIndex_loc = NULL;
  unsigned long *Buffer_Receive_Triangle_loc         = NULL;
  unsigned long *Buffer_Receive_Quadrilateral_loc    = NULL;
  unsigned long *Buffer_Receive_Tetrahedron_loc      = NULL;
  unsigned long *Buffer_Receive_Hexahedron_loc       = NULL;
  unsigned long *Buffer_Receive_Prism_loc            = NULL;
  unsigned long *Buffer_Receive_Pyramid_loc          = NULL;
  
  unsigned long *Buffer_Receive_GlobElem_loc               = NULL;
  unsigned long *Buffer_Receive_Triangle_presence_loc      = NULL;
  unsigned long *Buffer_Receive_Quadrilateral_presence_loc = NULL;
  unsigned long *Buffer_Receive_Tetrahedron_presence_loc   = NULL;
  unsigned long *Buffer_Receive_Hexahedron_presence_loc    = NULL;
  unsigned long *Buffer_Receive_Prism_presence_loc         = NULL;
  unsigned long *Buffer_Receive_Pyramid_presence_loc       = NULL;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
#ifdef HAVE_MPI
  
  su2double *Buffer_Receive_Coord = NULL;
  
  unsigned long *Buffer_Receive_Color            = NULL;
  unsigned long *Buffer_Receive_GlobalPointIndex = NULL;
  unsigned long *Buffer_Receive_Triangle         = NULL;
  unsigned long *Buffer_Receive_Quadrilateral    = NULL;
  unsigned long *Buffer_Receive_Tetrahedron      = NULL;
  unsigned long *Buffer_Receive_Hexahedron       = NULL;
  unsigned long *Buffer_Receive_Prism            = NULL;
  unsigned long *Buffer_Receive_Pyramid          = NULL;
  unsigned long *Buffer_Receive_GlobElem         = NULL;
  
  unsigned long **Buffer_Receive_Triangle_presence      = new unsigned long*[size];
  unsigned long **Buffer_Receive_Quadrilateral_presence = new unsigned long*[size];
  unsigned long **Buffer_Receive_Tetrahedron_presence   = new unsigned long*[size];
  unsigned long **Buffer_Receive_Hexahedron_presence    = new unsigned long*[size];
  unsigned long **Buffer_Receive_Prism_presence         = new unsigned long*[size];
  unsigned long **Buffer_Receive_Pyramid_presence       = new unsigned long*[size];
  
#endif
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  Marker_All_SendRecv      = new short[nMarker_Max];
  nSendDomain_Periodic     = new unsigned long [nDomain];
  nReceivedDomain_Periodic = new unsigned long [nDomain];
  
  /*--- Auxiliar vector based on the original geometry ---*/
  
  ElemIn = new bool[geometry->GetnElem()];
  
  /*--- Define some mapping variables ---*/
  map<unsigned long,bool> PointIn;
  map<unsigned long,unsigned long> Global_to_local_Point_recv;
  
  Buffer_Send_nDim  = geometry->GetnDim();
  Buffer_Send_nZone = geometry->GetnZone();
  
  /*--- Divide the elements in color list to speed up the grid partitioning ---*/
  
  map<unsigned long,unsigned long> Local_to_global_elem;
  for (unsigned long i=0; i<geometry->GetGlobal_nElem(); i++) {
    map<unsigned long, unsigned long>::const_iterator MI = geometry->Global_to_Local_Elem.find(i);
    if (MI != geometry->Global_to_Local_Elem.end()) {
      Local_to_global_elem[geometry->Global_to_Local_Elem[i]] = i;
    }
  }
  
  /*--- MEMORY WARNING: Bad usage of memory here for local_colour_values. Not scalable. 
        In the future, we should avoid sharing a single array of all colors in all nodes. ---*/
  unsigned long *local_colour_values = new unsigned long[geometry->GetGlobal_nPoint()];
  unsigned long *local_colour_temp   = new unsigned long[geometry->ending_node[rank]-geometry->starting_node[rank]];
  
  for (unsigned long i=0; i<geometry->ending_node[rank]-geometry->starting_node[rank]; i++) {
    local_colour_temp[i]=geometry->node[i]->GetColor();
    local_colour_values[geometry->starting_node[rank]+i]=local_colour_temp[i];
  }
  
  /*--- Communicate the grid coloring to all partitions. This information
   will be repeatedly used throughout the organization of the partitions
   and sorting out their ghost points/elements. ---*/
  
#ifdef HAVE_MPI
  
  int comm_counter=0;
  for (iDomain=0; iDomain < (unsigned long)size; iDomain++) {
    if (iDomain != (unsigned long)rank) {
      SU2_MPI::Isend(local_colour_temp, geometry->ending_node[rank]-geometry->starting_node[rank],
                     MPI_UNSIGNED_LONG, iDomain, iDomain,  MPI_COMM_WORLD, &send_req[comm_counter]);
      comm_counter++;
    }
  }
  
  for (iDomain=0; iDomain < (unsigned long)size-1; iDomain++) {
    MPI_Probe(MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status2);
    source = status2.MPI_SOURCE;
    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
    SU2_MPI::Recv(&local_colour_values[geometry->starting_node[source]], recv_count,
                  MPI_UNSIGNED_LONG, source, rank, MPI_COMM_WORLD, &status2);
  }
  
  /*--- Wait for the sends to complete (will be true since we're using
   blocking recv's above. ---*/
  
  SU2_MPI::Waitall(size-1, send_req, send_stat);
  
#endif
  
  /*--- Free temporary buffer for communicating colors. ---*/
  
  delete [] local_colour_temp;
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- This loop gets the array sizes of points, elements, etc. for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Interior dimensionalization. Loop over the original grid to
     perform the dimensionalizaton of the domain variables ---*/
    
    Buffer_Send_nElemTotal         = 0;
    Buffer_Send_nPointTotal        = 0;
    Buffer_Send_nPointGhost        = 0;
    Buffer_Send_nPointDomainTotal  = 0;
    Buffer_Send_nPointPeriodic     = 0;
    Buffer_Send_nElemTriangle      = 0;
    Buffer_Send_nElemQuadrilateral = 0;
    Buffer_Send_nElemTetrahedron   = 0;
    Buffer_Send_nElemHexahedron    = 0;
    Buffer_Send_nElemPrism         = 0;
    Buffer_Send_nElemPyramid       = 0;
    
    /*--- Initialize the global to local mapping ---*/
    
    PointIn.clear();
    
    /*--- Loop over all of the local elements and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      /*--- Check if the element belongs to the domain ---*/
      
      ElemIn[iElem] = false;
      for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        if (local_colour_values[iPoint] == iDomain) {
          ElemIn[iElem] = true; break;
        }
      }
      
      /*--- If this element is needed by iDomain, get information
       about the number of points and element type. ---*/
      
      if (ElemIn[iElem]) {
        
        for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
          iPoint = geometry->elem[iElem]->GetNode(iNode);
          
          /*--- If we haven't already found this point... ---*/
          
          map<unsigned long, bool>::const_iterator MI = PointIn.find(iPoint);
          if (MI == PointIn.end()) {
            
            /*--- Mark point as found and collect information ---*/
            
            PointIn[iPoint] = true;
            
            if ((iPoint >= geometry->starting_node[rank]) &&
                (iPoint < geometry->ending_node[rank])) {
              
              Buffer_Send_nPointTotal++;
              
              /*--- Increment our counters ---*/
              if ( local_colour_values[iPoint] == iDomain ) {
                if ( iPoint > geometry->GetGlobal_nPointDomain() - 1)
                  Buffer_Send_nPointPeriodic++;
                else
                  Buffer_Send_nPointDomainTotal++;
              }
              else Buffer_Send_nPointGhost++;
              
              
            }
          }
        }
        
        /*--- Increment the counter for the current type of element ---*/
        
        switch(geometry->elem[iElem]->GetVTK_Type()) {
          case TRIANGLE:      Buffer_Send_nElemTriangle++;      break;
          case QUADRILATERAL: Buffer_Send_nElemQuadrilateral++; break;
          case TETRAHEDRON:   Buffer_Send_nElemTetrahedron++;   break;
          case HEXAHEDRON:    Buffer_Send_nElemHexahedron++;    break;
          case PRISM:         Buffer_Send_nElemPrism++;         break;
          case PYRAMID:       Buffer_Send_nElemPyramid++;       break;
        }
        
        /*--- Increment the total number of elements for iDomain ---*/
        
        Buffer_Send_nElemTotal++;
        
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nDim_s[iDomain]               = geometry->GetnDim();
    nZone_s[iDomain]              = Buffer_Send_nZone;
    nPointTotal_s[iDomain]        = Buffer_Send_nPointTotal;
    nPointDomainTotal_s[iDomain]  = Buffer_Send_nPointDomainTotal;
    nPointGhost_s[iDomain]        = Buffer_Send_nPointGhost;
    nPointPeriodic_s[iDomain]     = Buffer_Send_nPointPeriodic;
    nElemTotal_s[iDomain]         = Buffer_Send_nElemTotal;
    nElemTriangle_s[iDomain]      = Buffer_Send_nElemTriangle;
    nElemQuadrilateral_s[iDomain] = Buffer_Send_nElemQuadrilateral;
    nElemTetrahedron_s[iDomain]   = Buffer_Send_nElemTetrahedron;
    nElemHexahedron_s[iDomain]    = Buffer_Send_nElemHexahedron;
    nElemPrism_s[iDomain]         = Buffer_Send_nElemPrism;
    nElemPyramid_s[iDomain]       = Buffer_Send_nElemPyramid;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_Coord            += nPointTotal_s[iDomain]*nDim_s[iDomain];
    Buffer_Size_Color            += nPointTotal_s[iDomain];
    Buffer_Size_GlobalPointIndex += nPointTotal_s[iDomain];
    Buffer_Size_Triangle         += nElemTriangle_s[iDomain];
    Buffer_Size_Quadrilateral    += nElemQuadrilateral_s[iDomain];
    Buffer_Size_Tetrahedron      += nElemTetrahedron_s[iDomain];
    Buffer_Size_Hexahedron       += nElemHexahedron_s[iDomain];
    Buffer_Size_Prism            += nElemPrism_s[iDomain];
    Buffer_Size_Pyramid          += nElemPyramid_s[iDomain];
    Buffer_Size_GlobElem         += nElemTotal_s[iDomain];
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_Coord = new su2double[Buffer_Size_Coord];
  
  Buffer_Send_Color             = new unsigned long[Buffer_Size_Color];
  Buffer_Send_GlobalPointIndex  = new unsigned long[Buffer_Size_GlobalPointIndex];
  Buffer_Send_Triangle          = new unsigned long[Buffer_Size_Triangle*N_POINTS_TRIANGLE];
  Buffer_Send_Quadrilateral     = new unsigned long[Buffer_Size_Quadrilateral*N_POINTS_QUADRILATERAL];
  Buffer_Send_Tetrahedron       = new unsigned long[Buffer_Size_Tetrahedron*N_POINTS_TETRAHEDRON];
  Buffer_Send_Hexahedron        = new unsigned long[Buffer_Size_Hexahedron*N_POINTS_HEXAHEDRON];
  Buffer_Send_Prism             = new unsigned long[Buffer_Size_Prism*N_POINTS_PRISM];
  Buffer_Send_Pyramid           = new unsigned long[Buffer_Size_Pyramid*N_POINTS_PYRAMID];
  Buffer_Send_GlobElem          = new unsigned long[Buffer_Size_GlobElem];
  
  Local_to_global_Triangle      = new unsigned long[Buffer_Size_Triangle];
  Local_to_global_Quadrilateral = new unsigned long[Buffer_Size_Quadrilateral];
  Local_to_global_Tetrahedron   = new unsigned long[Buffer_Size_Tetrahedron];
  Local_to_global_Hexahedron    = new unsigned long[Buffer_Size_Hexahedron];
  Local_to_global_Prism         = new unsigned long[Buffer_Size_Prism];
  Local_to_global_Pyramid       = new unsigned long[Buffer_Size_Pyramid];
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  ElemTotal_Counter         = 0;
  PointTotal_Counter        = 0;
  PointDomain_Counter       = 0;
  /*--- WARNING: check the next two counters ---*/
  PointPeriodic_Counter     = 0;
  PointGhost_Counter        = 0;
  ElemTriangle_Counter      = 0;
  ElemQuadrilateral_Counter = 0;
  ElemTetrahedron_Counter   = 0;
  ElemHexahedron_Counter    = 0;
  ElemPrism_Counter         = 0;
  ElemPyramid_Counter       = 0;
  
  /*--- Now that we know the sizes of the point, elem, etc. arrays, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if ((unsigned long)rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Isend(&nDim_s[iDomain], 1, MPI_UNSIGNED_SHORT, iDomain,
                     iDomain*13+0, MPI_COMM_WORLD, &send_req[0]);
      
      SU2_MPI::Isend(&nZone_s[iDomain], 1, MPI_UNSIGNED_SHORT, iDomain,
                     iDomain*13+1, MPI_COMM_WORLD, &send_req[1]);
      
      SU2_MPI::Isend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+2, MPI_COMM_WORLD, &send_req[2]);
      
      SU2_MPI::Isend(&nPointDomainTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+3, MPI_COMM_WORLD, &send_req[3]);
      
      SU2_MPI::Isend(&nPointGhost_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+4, MPI_COMM_WORLD, &send_req[4]);
      
      SU2_MPI::Isend(&nPointPeriodic_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+5, MPI_COMM_WORLD, &send_req[5]);
      
      SU2_MPI::Isend(&nElemTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+6, MPI_COMM_WORLD, &send_req[6]);
      
      SU2_MPI::Isend(&nElemTriangle_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+7, MPI_COMM_WORLD, &send_req[7]);
      
      SU2_MPI::Isend(&nElemQuadrilateral_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+8, MPI_COMM_WORLD, &send_req[8]);
      
      SU2_MPI::Isend(&nElemTetrahedron_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+9, MPI_COMM_WORLD, &send_req[9]);
      
      SU2_MPI::Isend(&nElemHexahedron_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+10, MPI_COMM_WORLD, &send_req[10]);
      
      SU2_MPI::Isend(&nElemPrism_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+11, MPI_COMM_WORLD, &send_req[11]);
      
      SU2_MPI::Isend(&nElemPyramid_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*13+12, MPI_COMM_WORLD, &send_req[12]);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nDim              = nDim_s[iDomain];
      nZone             = nZone_s[iDomain];
      
      //      nPointTotal        = nPointTotal_s[iDomain];
      //      nPointDomainTotal  = nPointDomainTotal_s[iDomain];
      //      nPointGhost        = nPointGhost_s[iDomain];
      //      nPointPeriodic     = nPointPeriodic_s[iDomain];
      //      nElemTotal         = nElemTotal_s[iDomain];
      //      nElemTriangle      = nElemTriangle_s[iDomain];
      //      nElemQuadrilateral = nElemQuadrilateral_s[iDomain];
      //      nElemTetrahedron   = nElemTetrahedron_s[iDomain];
      //      nElemHexahedron    = nElemHexahedron_s[iDomain];
      //      nElemPrism         = nElemPrism_s[iDomain];
      //      nElemPyramid       = nElemPyramid_s[iDomain];
      
      nDim_r[iDomain]               = nDim_s[iDomain];
      nZone_r[iDomain]              = nZone_s[iDomain];
      nPointTotal_r[iDomain]        = nPointTotal_s[iDomain];
      nPointDomainTotal_r[iDomain]  = nPointDomainTotal_s[iDomain];
      nPointPeriodic_r[iDomain]     = nPointPeriodic_s[iDomain];
      nElemTotal_r[iDomain]         = nElemTotal_s[iDomain];
      nElemTriangle_r[iDomain]      = nElemTriangle_s[iDomain];
      nElemQuadrilateral_r[iDomain] = nElemQuadrilateral_s[iDomain];
      nElemTetrahedron_r[iDomain]   = nElemTetrahedron_s[iDomain];
      nElemHexahedron_r[iDomain]    = nElemHexahedron_s[iDomain];
      nElemPrism_r[iDomain]         = nElemPrism_s[iDomain];
      nElemPyramid_r[iDomain]       = nElemPyramid_s[iDomain];
      
      nPointTotal_r_tot        += nPointTotal_r[iDomain];
      nPointDomainTotal_r_tot  += nPointDomainTotal_r[iDomain];
      nPointGhost_r_tot        += nPointGhost_r[iDomain];
      nPointPeriodic_r_tot     += nPointPeriodic_r[iDomain];
      nElemTotal_r_tot         += nElemTotal_r[iDomain];
      nElemTriangle_r_tot      += nElemTriangle_r[iDomain];
      nElemQuadrilateral_r_tot += nElemQuadrilateral_r[iDomain];
      nElemTetrahedron_r_tot   += nElemTetrahedron_r[iDomain];
      nElemHexahedron_r_tot    += nElemHexahedron_r[iDomain];
      nElemPrism_r_tot         += nElemPrism_r[iDomain];
      nElemPyramid_r_tot       += nElemPyramid_r[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if ((unsigned long)rank == iDomain) {
      
      for (jDomain = 0; jDomain < (unsigned long)size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if ((unsigned long)rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          MPI_Probe(jDomain, 13*rank+0, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nDim_r[jDomain], 1, MPI_UNSIGNED_SHORT, jDomain,
                        rank*13+0, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+1, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nZone_r[jDomain], 1, MPI_UNSIGNED_SHORT, jDomain,
                        rank*13+1, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+2, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+2, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+3, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nPointDomainTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+3, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+4, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nPointGhost_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+4, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+5, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nPointPeriodic_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+5, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+6, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nElemTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+6, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+7, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nElemTriangle_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+7, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+8, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nElemQuadrilateral_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+8, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+9, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nElemTetrahedron_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+9, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+10, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nElemHexahedron_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+10, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+11, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nElemPrism_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+11, MPI_COMM_WORLD, &status2);
          
          MPI_Probe(jDomain, 13*rank+12, MPI_COMM_WORLD, &status2);
          SU2_MPI::Recv(&nElemPyramid_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                        rank*13+12, MPI_COMM_WORLD, &status2);
          
#endif
          
          /*--- These are the cumulative totals that we will recv below. ----*/
          
          nPointTotal_r_tot        += nPointTotal_r[jDomain];
          nPointDomainTotal_r_tot  += nPointDomainTotal_r[jDomain];
          nPointGhost_r_tot        += nPointGhost_r[jDomain];
          nPointPeriodic_r_tot     += nPointPeriodic_r[jDomain];
          nElemTotal_r_tot         += nElemTotal_r[jDomain];
          nElemTriangle_r_tot      += nElemTriangle_r[jDomain];
          nElemQuadrilateral_r_tot += nElemQuadrilateral_r[jDomain];
          nElemTetrahedron_r_tot   += nElemTetrahedron_r[jDomain];
          nElemHexahedron_r_tot    += nElemHexahedron_r[jDomain];
          nElemPrism_r_tot         += nElemPrism_r[jDomain];
          nElemPyramid_r_tot       += nElemPyramid_r[jDomain];
          
        }
      }
      
    }
  }
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Wait for the non-blocking sends to complete. ---*/
    
#ifdef HAVE_MPI
    if ((unsigned long)rank != iDomain) SU2_MPI::Waitall(13, send_req, send_stat);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  }
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Above was number of elements to send and receive, and here is where
     we send/recv the actual elements. Here you're sending global index values,
     which are later changed to local. ---*/
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iElemTotal       = 0;
    iPointTotal      = 0;
    iPointDomain     = 0;
    iPointPeriodic   = nPointDomainTotal_s[iDomain];
    iPointGhost      = nPointDomainTotal_s[iDomain] + nPointPeriodic_s[iDomain];
    iElemTriangle    = 0;
    iElemQuadrilateral   = 0;
    iElemTetrahedron = 0;
    iElemHexahedron  = 0;
    iElemPrism       = 0;
    iElemPyramid     = 0;
    
    /*--- Initialize the global to local mapping ---*/
    
    PointIn.clear();
    
    /*--- Load up the actual elements into the buffers for sending. ---*/
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      /*--- Check if the element belongs to the domain ---*/
      
      ElemIn[iElem] = false;
      for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        if (local_colour_values[iPoint] == iDomain) {
          ElemIn[iElem] = true; break;
        }
      }
      
      /*--- If this element should be sent ---*/
      
      if (ElemIn[iElem]) {
        
        /*--- We need to send this element, so add it to the send buffer. The
         local to global mapping has already been done as a class data member. ---*/
        
        Buffer_Send_GlobElem[ElemTotal_Counter+iElemTotal] = Local_to_global_elem[iElem];
        
        /*--- Loop through the nodes of the current element ---*/
        
        for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
          
          /*--- Get the global index for this node in the element ---*/
          iPoint = geometry->elem[iElem]->GetNode(iNode);
          
          /*--- Store the connectivity for this element for each node ---*/
          vnodes_local[iNode] = iPoint;
          
          /*--- Check if this point has been found previously ---*/
          
          map<unsigned long, bool>::const_iterator MI = PointIn.find(iPoint);
          if (MI == PointIn.end()) {
            
            /*--- Check if this node lives on the current rank based on the
             initial linear partitioning. We are only ever sending nodes that
             we own in the linear partitioning (no duplicate nodes are sent) ---*/
            
            if ((iPoint >= geometry->starting_node[rank]) &&
                (iPoint < geometry->ending_node[rank])) {
              
              /*--- Decide whether this is an interior, periodic, or ghost node ---*/
              
              if (local_colour_values[iPoint] == iDomain) {
                
                /*--- If iDomain owns the point, it must be either an interior
                 node (iPoint < nPointDomain) or a periodic node. ---*/
                
                if (iPoint > geometry->GetGlobal_nPointDomain() - 1)
                  iPointCurrent = iPointPeriodic;
                else
                  iPointCurrent = iPointDomain;
                
              } else {
                
                /*--- Otherwise, it must be a ghost point for iDomain ---*/
                iPointCurrent = iPointGhost;
                
              }
              
              /*--- Setting global to local, the color, and index. ---*/
              
              PointIn[iPoint] = true;
              
              Buffer_Send_Color[PointTotal_Counter+iPointCurrent] = local_colour_values[iPoint];
              Buffer_Send_GlobalPointIndex[PointTotal_Counter+iPointCurrent] = iPoint;
              
              /*--- Get the coordinates for this point ---*/
              
              for (iDim = 0; iDim < nDim_s[iDomain]; iDim++) {
                
                /*--- iPoint is the global index, but we store everything local
                 to this rank. So we need to subtract the starting index. All
                 ranks re-index their points from zero. ---*/
                Buffer_Send_Coord[nDim_s[iDomain]*(PointTotal_Counter+iPointCurrent)+iDim] = geometry->node[iPoint-geometry->starting_node[rank]]->GetCoord(iDim);
              }
              
              /*--- Increment our counters ---*/
              if ( local_colour_values[iPoint] == iDomain ) {
                if ( iPoint > geometry->GetGlobal_nPointDomain() - 1)
                  iPointPeriodic++;
                else
                  iPointDomain++;
              }
              else iPointGhost++;
              
              /*--- Increment the total number of points we're sending ---*/
              iPointTotal++;
              
            }
          }
        }
        
        /*--- Load the connectivity for the current element into the send buffer.
         Also store the local to global mapping for the elements.
         Note that we are using the vnode_local array we filled above to store
         the connectivity. Loop through each element type. ---*/
        
        switch(geometry->elem[iElem]->GetVTK_Type()) {
          case TRIANGLE:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Triangle[3*(ElemTriangle_Counter+iElemTriangle)+iNode] = vnodes_local[iNode];
            Local_to_global_Triangle[ElemTriangle_Counter+iElemTriangle] = Buffer_Send_GlobElem[ElemTotal_Counter+iElemTotal];
            iElemTriangle++; break;
          case QUADRILATERAL:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Quadrilateral[4*(ElemQuadrilateral_Counter+iElemQuadrilateral)+iNode] = vnodes_local[iNode];
            Local_to_global_Quadrilateral[ElemQuadrilateral_Counter+iElemQuadrilateral] =Buffer_Send_GlobElem[ElemTotal_Counter+iElemTotal];
            iElemQuadrilateral++; break;
          case TETRAHEDRON:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Tetrahedron[4*(ElemTetrahedron_Counter+iElemTetrahedron)+iNode] = vnodes_local[iNode];
            Local_to_global_Tetrahedron[ElemTetrahedron_Counter+iElemTetrahedron] =Buffer_Send_GlobElem[ElemTotal_Counter+iElemTotal];
            iElemTetrahedron++; break;
          case HEXAHEDRON:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Hexahedron[8*(ElemHexahedron_Counter+iElemHexahedron)+iNode] = vnodes_local[iNode];
            Local_to_global_Hexahedron[ElemHexahedron_Counter+iElemHexahedron] =Buffer_Send_GlobElem[ElemTotal_Counter+iElemTotal];
            iElemHexahedron++; break;
          case PRISM:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Prism[6*(ElemPrism_Counter+iElemPrism)+iNode] = vnodes_local[iNode];
            Local_to_global_Prism[ElemPrism_Counter+iElemPrism] =Buffer_Send_GlobElem[ElemTotal_Counter+iElemTotal];
            iElemPrism++; break;
          case PYRAMID:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Pyramid[5*(ElemPyramid_Counter+iElemPyramid)+iNode] = vnodes_local[iNode];
            Local_to_global_Pyramid[ElemPyramid_Counter+iElemPyramid] = Buffer_Send_GlobElem[ElemTotal_Counter+iElemTotal];
            iElemPyramid++; break;
        }
        
        /*--- Regardless of the type, increment the total count ---*/
        iElemTotal++;
        
      }
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != (unsigned long)rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Isend(&Buffer_Send_Coord[PointTotal_Counter*nDim_s[iDomain]],
                     nPointTotal_s[iDomain]*nDim_s[iDomain], MPI_DOUBLE, iDomain,
                     iDomain*16+0,  MPI_COMM_WORLD, &send_req[0]);
      
      SU2_MPI::Isend(&Buffer_Send_GlobalPointIndex[PointTotal_Counter],
                     nPointTotal_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+1,  MPI_COMM_WORLD, &send_req[1]);
      
      SU2_MPI::Isend(&Buffer_Send_Color[PointTotal_Counter],
                     nPointTotal_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+2,  MPI_COMM_WORLD, &send_req[2]);
      
      SU2_MPI::Isend(&Buffer_Send_Triangle[ElemTriangle_Counter*3],
                     nElemTriangle_s[iDomain]*3, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+3,  MPI_COMM_WORLD, &send_req[3]);
      
      SU2_MPI::Isend(&Buffer_Send_Quadrilateral[ElemQuadrilateral_Counter*4],
                     nElemQuadrilateral_s[iDomain]*4, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+4,  MPI_COMM_WORLD, &send_req[4]);
      
      SU2_MPI::Isend(&Buffer_Send_Tetrahedron[ElemTetrahedron_Counter*4],
                     nElemTetrahedron_s[iDomain]*4, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+5,  MPI_COMM_WORLD, &send_req[5]);
      
      SU2_MPI::Isend(&Buffer_Send_Hexahedron[ElemHexahedron_Counter*8],
                     nElemHexahedron_s[iDomain]*8, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+6,  MPI_COMM_WORLD, &send_req[6]);
      
      SU2_MPI::Isend(&Buffer_Send_Prism[ElemPrism_Counter*6],
                     nElemPrism_s[iDomain]*6, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+7,  MPI_COMM_WORLD, &send_req[7]);
      
      SU2_MPI::Isend(&Buffer_Send_Pyramid[ElemPyramid_Counter*5],
                     nElemPyramid_s[iDomain]*5, MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+8,  MPI_COMM_WORLD, &send_req[8]);
      
      SU2_MPI::Isend(&Buffer_Send_GlobElem[ElemTotal_Counter],
                     nElemTotal_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+9,  MPI_COMM_WORLD, &send_req[9]);
      
      SU2_MPI::Isend(&Local_to_global_Triangle[ElemTriangle_Counter],
                     nElemTriangle_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+10,  MPI_COMM_WORLD, &send_req[10]);
      
      SU2_MPI::Isend(&Local_to_global_Quadrilateral[ElemQuadrilateral_Counter],
                     nElemQuadrilateral_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+11,  MPI_COMM_WORLD, &send_req[11]);
      
      SU2_MPI::Isend(&Local_to_global_Tetrahedron[ElemTetrahedron_Counter],
                     nElemTetrahedron_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+12,  MPI_COMM_WORLD, &send_req[12]);
      
      SU2_MPI::Isend(&Local_to_global_Hexahedron[ElemHexahedron_Counter],
                     nElemHexahedron_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+13,  MPI_COMM_WORLD, &send_req[13]);
      
      SU2_MPI::Isend(&Local_to_global_Prism[ElemPrism_Counter],
                     nElemPrism_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+14,  MPI_COMM_WORLD, &send_req[14]);
      
      SU2_MPI::Isend(&Local_to_global_Pyramid[ElemPyramid_Counter],
                     nElemPyramid_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                     iDomain*16+15,  MPI_COMM_WORLD, &send_req[15]);
      
#endif
      
    } else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_Coord_loc  = new su2double[nPointTotal_s[iDomain]*nDim_s[iDomain]];
      
      Buffer_Receive_GlobalPointIndex_loc = new unsigned long[nPointTotal_s[iDomain]];
      Buffer_Receive_Color_loc            = new unsigned long[nPointTotal_s[iDomain]];
      Buffer_Receive_Triangle_loc         = new unsigned long[nElemTriangle_s[iDomain]*N_POINTS_TRIANGLE];
      Buffer_Receive_Quadrilateral_loc    = new unsigned long[nElemQuadrilateral_s[iDomain]*N_POINTS_QUADRILATERAL];
      Buffer_Receive_Tetrahedron_loc      = new unsigned long[nElemTetrahedron_s[iDomain]*N_POINTS_TETRAHEDRON];
      Buffer_Receive_Hexahedron_loc       = new unsigned long[nElemHexahedron_s[iDomain]*N_POINTS_HEXAHEDRON];
      Buffer_Receive_Prism_loc            = new unsigned long[nElemPrism_s[iDomain]*N_POINTS_PRISM];
      Buffer_Receive_Pyramid_loc          = new unsigned long[nElemPyramid_s[iDomain]*N_POINTS_PYRAMID];
      Buffer_Receive_GlobElem_loc         = new unsigned long[nElemTotal_s[iDomain]];
      
      Buffer_Receive_Triangle_presence_loc      = new unsigned long[nElemTriangle_s[iDomain]];
      Buffer_Receive_Quadrilateral_presence_loc = new unsigned long[nElemQuadrilateral_s[iDomain]];
      Buffer_Receive_Tetrahedron_presence_loc   = new unsigned long[nElemTetrahedron_s[iDomain]];
      Buffer_Receive_Hexahedron_presence_loc    = new unsigned long[nElemHexahedron_s[iDomain]];
      Buffer_Receive_Prism_presence_loc         = new unsigned long[nElemPrism_s[iDomain]];
      Buffer_Receive_Pyramid_presence_loc       = new unsigned long[nElemPyramid_s[iDomain]];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*nDim_s[iDomain]; iter++)
        Buffer_Receive_Coord_loc[iter] = Buffer_Send_Coord[PointTotal_Counter*nDim_s[iDomain]+iter];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]; iter++) {
        Buffer_Receive_GlobalPointIndex_loc[iter] = Buffer_Send_GlobalPointIndex[PointTotal_Counter+iter];
        Buffer_Receive_Color_loc[iter] = Buffer_Send_Color[PointTotal_Counter+iter];
      }
      
      for (iter = 0; iter < nElemTriangle_s[iDomain]*N_POINTS_TRIANGLE; iter++)
        Buffer_Receive_Triangle_loc[iter] =  Buffer_Send_Triangle[ElemTriangle_Counter*N_POINTS_TRIANGLE+iter];
      
      for (iter = 0; iter < nElemQuadrilateral_s[iDomain]*N_POINTS_QUADRILATERAL; iter++)
        Buffer_Receive_Quadrilateral_loc[iter] =  Buffer_Send_Quadrilateral[ElemQuadrilateral_Counter*N_POINTS_QUADRILATERAL+iter];
      
      for (iter = 0; iter < nElemTetrahedron_s[iDomain]*N_POINTS_TETRAHEDRON; iter++)
        Buffer_Receive_Tetrahedron_loc[iter] =  Buffer_Send_Tetrahedron[ElemTetrahedron_Counter*N_POINTS_TETRAHEDRON+iter];
      
      for (iter = 0; iter < nElemHexahedron_s[iDomain]*N_POINTS_HEXAHEDRON; iter++)
        Buffer_Receive_Hexahedron_loc[iter] =  Buffer_Send_Hexahedron[ElemHexahedron_Counter*N_POINTS_HEXAHEDRON+iter];
      
      for (iter = 0; iter < nElemPrism_s[iDomain]*N_POINTS_PRISM; iter++)
        Buffer_Receive_Prism_loc[iter] =  Buffer_Send_Prism[ElemPrism_Counter*N_POINTS_PRISM+iter];
      
      for (iter = 0; iter < nElemPyramid_s[iDomain]*N_POINTS_PYRAMID; iter++)
        Buffer_Receive_Pyramid_loc[iter] =  Buffer_Send_Pyramid[ElemPyramid_Counter*N_POINTS_PYRAMID+iter];
      
      for (unsigned long i=0; i<nElemTotal_s[iDomain]; i++) {
        Buffer_Receive_GlobElem_loc[i]=Buffer_Send_GlobElem[ElemTotal_Counter+i];
      }
      
      for (unsigned long i=0; i<nElemTriangle_s[iDomain]; i++) {
        Buffer_Receive_Triangle_presence_loc[i]=Local_to_global_Triangle[ElemTriangle_Counter+i];
      }
      
      for (unsigned long i=0; i<nElemQuadrilateral_s[iDomain]; i++) {
        Buffer_Receive_Quadrilateral_presence_loc[i]=Local_to_global_Quadrilateral[ElemQuadrilateral_Counter+i];
      }
      
      for (unsigned long i=0; i<nElemTetrahedron_s[iDomain]; i++) {
        Buffer_Receive_Tetrahedron_presence_loc[i]=Local_to_global_Tetrahedron[ElemTetrahedron_Counter+i];
      }
      
      for (unsigned long i=0; i<nElemHexahedron_s[iDomain]; i++) {
        Buffer_Receive_Hexahedron_presence_loc[i]=Local_to_global_Hexahedron[ElemHexahedron_Counter+i];
      }
      
      for (unsigned long i=0; i<nElemPrism_s[iDomain]; i++) {
        Buffer_Receive_Prism_presence_loc[i]=Local_to_global_Prism[ElemPrism_Counter+i];
      }
      
      for (unsigned long i=0; i<nElemPyramid_s[iDomain]; i++) {
        Buffer_Receive_Pyramid_presence_loc[i]=Local_to_global_Pyramid[ElemPyramid_Counter+i];
      }
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    ElemTotal_Counter         += iElemTotal;
    PointTotal_Counter        += iPointTotal;
    PointDomain_Counter       += iPointDomain;
    /*--- WARNING: check the next two counters ---*/
    PointPeriodic_Counter     += iPointPeriodic;
    PointGhost_Counter        += iPointGhost;
    ElemTriangle_Counter      += iElemTriangle;
    ElemQuadrilateral_Counter += iElemQuadrilateral;
    ElemTetrahedron_Counter   += iElemTetrahedron;
    ElemHexahedron_Counter    += iElemHexahedron;
    ElemPrism_Counter         += iElemPrism;
    ElemPyramid_Counter       += iElemPyramid;
    
  }
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank ---*/
  
  nPoint = nPointTotal_r_tot;
  nPointDomain = nPointDomainTotal_r_tot;
  nPointNode = nPoint;
  node = new CPoint*[nPoint];
  Local_to_Global_Point = new long[nPoint];
  
  /*--- Array initialization ---*/
  
  for (iPoint = 0; iPoint < nPointTotal_r_tot; iPoint++) {
    Local_to_Global_Point[iPoint] = -1;
  }
  
  /*--- Initialize some counters ---*/
  
  unsigned long temp_node_count = 0;
  unsigned long temp_node_count_periodic = nPointDomainTotal_r_tot;
  unsigned long temp_node_count_ghost = nPointDomainTotal_r_tot+nPointPeriodic_r_tot;
  
  
  /*--- First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < (unsigned long)size; iDomain++) {
    
    if ((unsigned long)rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_Coord =  new su2double [nPointTotal_r[iDomain]*nDim_r[iDomain]];
      Buffer_Receive_Color =  new unsigned long [nPointTotal_r[iDomain]];
      Buffer_Receive_GlobalPointIndex = new unsigned long [nPointTotal_r[iDomain]];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      MPI_Probe(iDomain, rank*16+0, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_DOUBLE, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_Coord, recv_count , MPI_DOUBLE,
                    source, rank*16+0, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+1, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_GlobalPointIndex, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+1, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+2, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_Color, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+2, MPI_COMM_WORLD, &status2);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index, and colors ---*/
      
      unsigned long index=0;
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        /*--- If this rank owns the current point ---*/
        
        if (Buffer_Receive_Color[iPoint] == (unsigned long)rank) {
          
          /*--- If iDomain owns the point, it must be either an interior
           node (iPoint < nPointDomain) or a periodic node. ---*/
          
          if (Buffer_Receive_GlobalPointIndex[iPoint] > geometry->GetGlobal_nPointDomain() - 1) {
            
            /*--- Set the starting point for the local index of the recv points.
             The temp_node_count increments for the interior nodes, between 0 up
             to nPointDomain-1. ---*/
            index = temp_node_count_periodic;
            
            /*--- Get the global index ---*/
            Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex[iPoint];
            
            /*--- Allocating the Point object ---*/
            if ( nDim == 2 ) node[index] = new CPoint(Buffer_Receive_Coord[iPoint*nDim+0],
                                                      Buffer_Receive_Coord[iPoint*nDim+1],
                                                      Local_to_Global_Point[index], config);
            if ( nDim == 3 ) node[index] = new CPoint(Buffer_Receive_Coord[iPoint*nDim+0],
                                                      Buffer_Receive_Coord[iPoint*nDim+1],
                                                      Buffer_Receive_Coord[iPoint*nDim+2],
                                                      Local_to_Global_Point[index], config);
            
            /*--- Set the color ---*/
            node[index]->SetColor(Buffer_Receive_Color[iPoint]);
            
            /*--- Increment the interior node counter ---*/
            temp_node_count_periodic++;
            
            
          }
          
          else {
            
            
            /*--- Set the starting point for the local index of the recv points.
             The temp_node_count increments for the interior nodes, between 0 up
             to nPointDomain-1. ---*/
            index = temp_node_count;
            
            /*--- Get the global index ---*/
            Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex[iPoint];
            
            /*--- Allocating the Point object ---*/
            if ( nDim == 2 ) node[index] = new CPoint(Buffer_Receive_Coord[iPoint*nDim+0],
                                                      Buffer_Receive_Coord[iPoint*nDim+1],
                                                      Local_to_Global_Point[index], config);
            if ( nDim == 3 ) node[index] = new CPoint(Buffer_Receive_Coord[iPoint*nDim+0],
                                                      Buffer_Receive_Coord[iPoint*nDim+1],
                                                      Buffer_Receive_Coord[iPoint*nDim+2],
                                                      Local_to_Global_Point[index], config);
            
            /*--- Set the color ---*/
            node[index]->SetColor(Buffer_Receive_Color[iPoint]);
            
            /*--- Increment the interior node counter ---*/
            temp_node_count++;
            
            
            
            
          }
          
          
        } else {
          
          /*--- Set the starting point for the local index of the recv points.
           The temp_node_count_domain increments for the ghost nodes, between
           nPointDomain up to nPoint. ---*/
          
          index=temp_node_count_ghost;
          
          /*--- Get the global index ---*/
          Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex[iPoint];
          
          /*--- Allocating the Point object ---*/
          if ( nDim == 2 ) node[index] = new CPoint(Buffer_Receive_Coord[iPoint*nDim+0],
                                                    Buffer_Receive_Coord[iPoint*nDim+1],
                                                    Local_to_Global_Point[index], config);
          if ( nDim == 3 ) node[index] = new CPoint(Buffer_Receive_Coord[iPoint*nDim+0],
                                                    Buffer_Receive_Coord[iPoint*nDim+1],
                                                    Buffer_Receive_Coord[iPoint*nDim+2],
                                                    Local_to_Global_Point[index], config);
          
          /*--- Set the color ---*/
          node[index]->SetColor(Buffer_Receive_Color[iPoint]);
          
          /*--- Increment the ghost node counter ---*/
          temp_node_count_ghost++;
          
        }
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      delete [] Buffer_Receive_Coord;
      delete [] Buffer_Receive_Color;
      delete [] Buffer_Receive_GlobalPointIndex;
      
#endif
      
    } else {
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      unsigned long index = 0;
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        if (Buffer_Receive_Color_loc[iPoint] == (unsigned long)rank) {
          
          /*--- If iDomain owns the point, it must be either an interior
           node (iPoint < nPointDomain) or a periodic node. ---*/
          
          if (Buffer_Receive_GlobalPointIndex_loc[iPoint] > geometry->GetGlobal_nPointDomain() - 1) {
            
            index = temp_node_count_periodic;
            
            Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex_loc[iPoint];
            if ( nDim == 2 ) node[index] = new CPoint(Buffer_Receive_Coord_loc[iPoint*nDim+0],
                                                      Buffer_Receive_Coord_loc[iPoint*nDim+1],
                                                      Local_to_Global_Point[index], config);
            if ( nDim == 3 ) node[index] = new CPoint(Buffer_Receive_Coord_loc[iPoint*nDim+0],
                                                      Buffer_Receive_Coord_loc[iPoint*nDim+1],
                                                      Buffer_Receive_Coord_loc[iPoint*nDim+2],
                                                      Local_to_Global_Point[index], config);
            node[index]->SetColor(Buffer_Receive_Color_loc[iPoint]);
            temp_node_count_periodic++;
            
            
            
            
          }
          else {
            
            index = temp_node_count;
            Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex_loc[iPoint];
            if ( nDim == 2 ) node[index] = new CPoint(Buffer_Receive_Coord_loc[iPoint*nDim+0],
                                                      Buffer_Receive_Coord_loc[iPoint*nDim+1],
                                                      Local_to_Global_Point[index], config);
            if ( nDim == 3 ) node[index] = new CPoint(Buffer_Receive_Coord_loc[iPoint*nDim+0],
                                                      Buffer_Receive_Coord_loc[iPoint*nDim+1],
                                                      Buffer_Receive_Coord_loc[iPoint*nDim+2],
                                                      Local_to_Global_Point[index], config);
            node[index]->SetColor(Buffer_Receive_Color_loc[iPoint]);
            temp_node_count++;
            
            
            
          }
          
          
        } else {
          
          index=temp_node_count_ghost;
          Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex_loc[iPoint];
          if ( nDim == 2 ) node[index] = new CPoint(Buffer_Receive_Coord_loc[iPoint*nDim+0],
                                                    Buffer_Receive_Coord_loc[iPoint*nDim+1],
                                                    Local_to_Global_Point[index], config);
          if ( nDim == 3 ) node[index] = new CPoint(Buffer_Receive_Coord_loc[iPoint*nDim+0],
                                                    Buffer_Receive_Coord_loc[iPoint*nDim+1],
                                                    Buffer_Receive_Coord_loc[iPoint*nDim+2],
                                                    Local_to_Global_Point[index], config);
          node[index]->SetColor(Buffer_Receive_Color_loc[iPoint]);
          temp_node_count_ghost++;
          
        }
      }
      
      delete [] Buffer_Receive_Coord_loc;
      delete [] Buffer_Receive_Color_loc;
      delete [] Buffer_Receive_GlobalPointIndex_loc;
      
    }
  }
  
  /*--- Get the global to local mapping ---*/
  
  for (iPoint = 0; iPoint < nPointTotal_r_tot; iPoint++) {
    Global_to_local_Point_recv[Local_to_Global_Point[iPoint]] = iPoint;
  }
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  /*--- Recv all of the element data. First decide which elements we need to own on each proc ---*/
  
  iElem = 0;
  for (iDomain = 0; iDomain < (unsigned long)size; iDomain++) {
    
    if ((unsigned long)rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate memory for the element recv ---*/
      
      Buffer_Receive_Triangle_presence[iDomain]      = new unsigned long[nElemTriangle_r[iDomain]];
      Buffer_Receive_Quadrilateral_presence[iDomain] = new unsigned long[nElemQuadrilateral_r[iDomain]];
      Buffer_Receive_Tetrahedron_presence[iDomain]   = new unsigned long[nElemTetrahedron_r[iDomain]];
      Buffer_Receive_Hexahedron_presence[iDomain]    = new unsigned long[nElemHexahedron_r[iDomain]];
      Buffer_Receive_Prism_presence[iDomain]         = new unsigned long[nElemPrism_r[iDomain]];
      Buffer_Receive_Pyramid_presence[iDomain]       = new unsigned long[nElemPyramid_r[iDomain]];
      
      /*--- Recv the element data ---*/
      
      MPI_Probe(iDomain, rank*16+10, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(&Buffer_Receive_Triangle_presence[iDomain][0],
                    recv_count, MPI_UNSIGNED_LONG, source,
                    rank*16+10, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+11, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(&Buffer_Receive_Quadrilateral_presence[iDomain][0],
                    recv_count, MPI_UNSIGNED_LONG, source,
                    rank*16+11, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+12, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(&Buffer_Receive_Tetrahedron_presence[iDomain][0],
                    recv_count, MPI_UNSIGNED_LONG, source,
                    rank*16+12, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+13, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(&Buffer_Receive_Hexahedron_presence[iDomain][0],
                    recv_count, MPI_UNSIGNED_LONG, source,
                    rank*16+13, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+14, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(&Buffer_Receive_Prism_presence[iDomain][0],
                    recv_count, MPI_UNSIGNED_LONG, source,
                    rank*16+14, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+15, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(&Buffer_Receive_Pyramid_presence[iDomain][0],
                    recv_count, MPI_UNSIGNED_LONG, source,
                    rank*16+15, MPI_COMM_WORLD, &status2);
      
      /*--- Allocating the elements after the recv ---*/
      
      for (iElemTriangle = 0; iElemTriangle < nElemTriangle_r[iDomain]; iElemTriangle++) {
        map<unsigned long, bool>::const_iterator MI = Triangle_presence.find(Buffer_Receive_Triangle_presence[iDomain][iElemTriangle]);
        if (MI == Triangle_presence.end()) {
          Triangle_presence[Buffer_Receive_Triangle_presence[iDomain][iElemTriangle]] = true;
          iElem++;
        }
      }
      
      for (iElemQuadrilateral = 0; iElemQuadrilateral < nElemQuadrilateral_r[iDomain]; iElemQuadrilateral++) {
        map<unsigned long, bool>::const_iterator MI = Quadrilateral_presence.find(Buffer_Receive_Quadrilateral_presence[iDomain][iElemQuadrilateral]);
        if (MI == Quadrilateral_presence.end()) {
          Quadrilateral_presence[Buffer_Receive_Quadrilateral_presence[iDomain][iElemQuadrilateral]] = true;
          iElem++;
        }
      }
      
      for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron_r[iDomain]; iElemTetrahedron++) {
        map<unsigned long, bool>::const_iterator MI = Tetrahedron_presence.find(Buffer_Receive_Tetrahedron_presence[iDomain][iElemTetrahedron]);
        if (MI == Tetrahedron_presence.end()) {
          Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence[iDomain][iElemTetrahedron]] = true;
          iElem++;
        }
      }
      
      for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron_r[iDomain]; iElemHexahedron++) {
        map<unsigned long, bool>::const_iterator MI = Hexahedron_presence.find(Buffer_Receive_Hexahedron_presence[iDomain][iElemHexahedron]);
        if (MI == Hexahedron_presence.end()) {
          Hexahedron_presence[Buffer_Receive_Hexahedron_presence[iDomain][iElemHexahedron]] = true;
          iElem++;
        }
      }
      
      for (iElemPrism = 0; iElemPrism < nElemPrism_r[iDomain]; iElemPrism++) {
        map<unsigned long, bool>::const_iterator MI = Prism_presence.find(Buffer_Receive_Prism_presence[iDomain][iElemPrism]);
        if (MI == Prism_presence.end()) {
          Prism_presence[Buffer_Receive_Prism_presence[iDomain][iElemPrism]] = true;
          iElem++;
        }
      }
      
      for (iElemPyramid = 0; iElemPyramid < nElemPyramid_r[iDomain]; iElemPyramid++) {
        map<unsigned long, bool>::const_iterator MI = Pyramid_presence.find(Buffer_Receive_Pyramid_presence[iDomain][iElemPyramid]);
        if (MI == Pyramid_presence.end()) {
          Pyramid_presence[Buffer_Receive_Pyramid_presence[iDomain][iElemPyramid]] = true;
          iElem++;
        }
      }
      
#endif
      
    } else {
      
      /*--- Store the element data from our own local rank info ---*/
      
      for (iElemTriangle = 0; iElemTriangle < nElemTriangle_r[iDomain]; iElemTriangle++) {
        map<unsigned long, bool>::const_iterator MI = Triangle_presence.find(Buffer_Receive_Triangle_presence_loc[iElemTriangle]);
        if (MI == Triangle_presence.end()) {
          Triangle_presence[Buffer_Receive_Triangle_presence_loc[iElemTriangle]] = true;
          iElem++;
        }
      }
      
      for (iElemQuadrilateral = 0; iElemQuadrilateral < nElemQuadrilateral_r[iDomain]; iElemQuadrilateral++) {
        map<unsigned long, bool>::const_iterator MI = Quadrilateral_presence.find(Buffer_Receive_Quadrilateral_presence_loc[iElemQuadrilateral]);
        if (MI == Quadrilateral_presence.end()) {
          Quadrilateral_presence[Buffer_Receive_Quadrilateral_presence_loc[iElemQuadrilateral]] = true;
          iElem++;
        }
      }
      
      for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron_r[iDomain]; iElemTetrahedron++) {
        map<unsigned long, bool>::const_iterator MI = Tetrahedron_presence.find(Buffer_Receive_Tetrahedron_presence_loc[iElemTetrahedron]);
        if (MI == Tetrahedron_presence.end()) {
          Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence_loc[iElemTetrahedron]] = true;
          iElem++;
        }
      }
      
      for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron_r[iDomain]; iElemHexahedron++) {
        map<unsigned long, bool>::const_iterator MI = Hexahedron_presence.find(Buffer_Receive_Hexahedron_presence_loc[iElemHexahedron]);
        if (MI == Hexahedron_presence.end()) {
          Hexahedron_presence[Buffer_Receive_Hexahedron_presence_loc[iElemHexahedron]] = true;
          iElem++;
        }
      }
      
      for (iElemPrism = 0; iElemPrism < nElemPrism_r[iDomain]; iElemPrism++) {
        map<unsigned long, bool>::const_iterator MI = Prism_presence.find(Buffer_Receive_Prism_presence_loc[iElemPrism]);
        if (MI == Prism_presence.end()) {
          Prism_presence[Buffer_Receive_Prism_presence_loc[iElemPrism]] = true;
          iElem++;
        }
      }
      
      for (iElemPyramid = 0; iElemPyramid < nElemPyramid_r[iDomain]; iElemPyramid++) {
        map<unsigned long, bool>::const_iterator MI = Pyramid_presence.find(Buffer_Receive_Pyramid_presence_loc[iElemPyramid]);
        if (MI == Pyramid_presence.end()) {
          Pyramid_presence[Buffer_Receive_Pyramid_presence_loc[iElemPyramid]] = true;
          iElem++;
        }
      }
      
    }
  }
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- iElem now contains the number of elements that this processor needs in
   total. Now we can complete the recv of the element connectivity and only
   store the elements that we need on this particular rank. Initialize space
   for the elements on this rank. ---*/
  
  nElem = iElem; iElem = 0;
  elem = new CPrimalGrid*[nElem];
  unsigned long iElemTria = 0;
  unsigned long iElemRect = 0;
  unsigned long iElemTetr = 0;
  unsigned long iElemHexa = 0;
  unsigned long iElemPris = 0;
  unsigned long iElemPyra = 0;
  
  unsigned long iElemRecv = 0;
  
  /*--- Reset presence before storing elems now that we know nElem ---*/
  
  Triangle_presence.clear();
  Quadrilateral_presence.clear();
  Tetrahedron_presence.clear();
  Hexahedron_presence.clear();
  Prism_presence.clear();
  Pyramid_presence.clear();
  
  /*--- Now recv all of the element connectivity data ---*/
  
  for (iDomain = 0; iDomain < (unsigned long)size; iDomain++) {
    
    if ((unsigned long)rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate memory for the element recv ---*/
      
      Buffer_Receive_Triangle      = new unsigned long[nElemTriangle_r[iDomain]*N_POINTS_TRIANGLE];
      Buffer_Receive_Quadrilateral = new unsigned long[nElemQuadrilateral_r[iDomain]*N_POINTS_QUADRILATERAL];
      Buffer_Receive_Tetrahedron   = new unsigned long[nElemTetrahedron_r[iDomain]*N_POINTS_TETRAHEDRON];
      Buffer_Receive_Hexahedron    = new unsigned long[nElemHexahedron_r[iDomain]*N_POINTS_HEXAHEDRON];
      Buffer_Receive_Prism         = new unsigned long[nElemPrism_r[iDomain]*N_POINTS_PRISM];
      Buffer_Receive_Pyramid       = new unsigned long[nElemPyramid_r[iDomain]*N_POINTS_PYRAMID];
      Buffer_Receive_GlobElem      = new unsigned long[nElemTotal_r[iDomain]];
      
      /*--- Recv the element data ---*/
      
      MPI_Probe(iDomain, rank*16+3, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_Triangle, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+3, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+4, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_Quadrilateral, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+4, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+5, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_Tetrahedron, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+5, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+6, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_Hexahedron, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+6, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+7, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_Prism, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+7, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+8, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_Pyramid, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+8, MPI_COMM_WORLD, &status2);
      
      MPI_Probe(iDomain, rank*16+9, MPI_COMM_WORLD, &status2);
      source = status2.MPI_SOURCE;
      MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
      SU2_MPI::Recv(Buffer_Receive_GlobElem, recv_count, MPI_UNSIGNED_LONG,
                    source, rank*16+9, MPI_COMM_WORLD, &status2);
      
      /*--- Allocating the elements after the recv. Note that here we are
       reusing the presence arrays to make sure that we find the exact same
       set of elements that were counted above to get nElem. ---*/
      
      iElemRecv = 0;
      
      for (iElemTriangle = 0; iElemTriangle < nElemTriangle_r[iDomain]; iElemTriangle++) {
        map<unsigned long, bool>::const_iterator MI = Triangle_presence.find(Buffer_Receive_Triangle_presence[iDomain][iElemTriangle]);
        if (MI == Triangle_presence.end()) {
          Triangle_presence[Buffer_Receive_Triangle_presence[iDomain][iElemTriangle]] = true;
          elem[iElem] = new CTriangle(Global_to_local_Point_recv[Buffer_Receive_Triangle[iElemTriangle*3+0]],
                                      Global_to_local_Point_recv[Buffer_Receive_Triangle[iElemTriangle*3+1]],
                                      Global_to_local_Point_recv[Buffer_Receive_Triangle[iElemTriangle*3+2]], 2);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem[iElemRecv]);
          iElem++; iElemTria++;
        }
        iElemRecv++;
      }
      
      for (iElemQuadrilateral = 0; iElemQuadrilateral < nElemQuadrilateral_r[iDomain]; iElemQuadrilateral++) {
        map<unsigned long, bool>::const_iterator MI = Quadrilateral_presence.find(Buffer_Receive_Quadrilateral_presence[iDomain][iElemQuadrilateral]);
        if (MI == Quadrilateral_presence.end()) {
          Quadrilateral_presence[Buffer_Receive_Quadrilateral_presence[iDomain][iElemQuadrilateral]] = true;
          elem[iElem] = new CQuadrilateral(Global_to_local_Point_recv[Buffer_Receive_Quadrilateral[iElemQuadrilateral*4+0]],
                                           Global_to_local_Point_recv[Buffer_Receive_Quadrilateral[iElemQuadrilateral*4+1]],
                                           Global_to_local_Point_recv[Buffer_Receive_Quadrilateral[iElemQuadrilateral*4+2]],
                                           Global_to_local_Point_recv[Buffer_Receive_Quadrilateral[iElemQuadrilateral*4+3]], 2);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem[iElemRecv]);
          iElem++; iElemRect++;
        }
        iElemRecv++;
      }
      
      for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron_r[iDomain]; iElemTetrahedron++) {
        map<unsigned long, bool>::const_iterator MI = Tetrahedron_presence.find(Buffer_Receive_Tetrahedron_presence[iDomain][iElemTetrahedron]);
        if (MI == Tetrahedron_presence.end()) {
          Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence[iDomain][iElemTetrahedron]] = true;
          elem[iElem] = new CTetrahedron(Global_to_local_Point_recv[Buffer_Receive_Tetrahedron[iElemTetrahedron*4+0]],
                                         Global_to_local_Point_recv[Buffer_Receive_Tetrahedron[iElemTetrahedron*4+1]],
                                         Global_to_local_Point_recv[Buffer_Receive_Tetrahedron[iElemTetrahedron*4+2]],
                                         Global_to_local_Point_recv[Buffer_Receive_Tetrahedron[iElemTetrahedron*4+3]]);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem[iElemRecv]);
          iElem++; iElemTetr++;
        }
        iElemRecv++;
      }
      
      for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron_r[iDomain]; iElemHexahedron++) {
        map<unsigned long, bool>::const_iterator MI = Hexahedron_presence.find(Buffer_Receive_Hexahedron_presence[iDomain][iElemHexahedron]);
        if (MI == Hexahedron_presence.end()) {
          Hexahedron_presence[Buffer_Receive_Hexahedron_presence[iDomain][iElemHexahedron]] = true;
          elem[iElem] = new CHexahedron(Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron*8+0]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron*8+1]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron*8+2]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron*8+3]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron*8+4]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron*8+5]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron*8+6]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron*8+7]]);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem[iElemRecv]);
          iElem++; iElemHexa++;
        }
        iElemRecv++;
      }
      
      for (iElemPrism = 0; iElemPrism < nElemPrism_r[iDomain]; iElemPrism++) {
        map<unsigned long, bool>::const_iterator MI = Prism_presence.find(Buffer_Receive_Prism_presence[iDomain][iElemPrism]);
        if (MI == Prism_presence.end()) {
          Prism_presence[Buffer_Receive_Prism_presence[iDomain][iElemPrism]] = true;
          elem[iElem] = new CPrism(Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism*6+0]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism*6+1]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism*6+2]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism*6+3]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism*6+4]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism*6+5]]);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem[iElemRecv]);
          iElem++; iElemPris++;
        }
        iElemRecv++;
      }
      
      for (iElemPyramid = 0; iElemPyramid < nElemPyramid_r[iDomain]; iElemPyramid++) {
        map<unsigned long, bool>::const_iterator MI = Pyramid_presence.find(Buffer_Receive_Pyramid_presence[iDomain][iElemPyramid]);
        if (MI == Pyramid_presence.end()) {
          Pyramid_presence[Buffer_Receive_Pyramid_presence[iDomain][iElemPyramid]] = true;
          elem[iElem] = new CPyramid(Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid*5+0]],
                                     Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid*5+1]],
                                     Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid*5+2]],
                                     Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid*5+3]],
                                     Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid*5+4]]);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem[iElemRecv]);
          iElem++; iElemPyra++;
        }
        iElemRecv++;
      }
      
      /*--- Free memory for the element data --*/
      
      delete[] Buffer_Receive_Triangle;
      delete[] Buffer_Receive_Quadrilateral;
      delete[] Buffer_Receive_Tetrahedron;
      delete[] Buffer_Receive_Hexahedron;
      delete[] Buffer_Receive_Prism;
      delete[] Buffer_Receive_Pyramid;
      delete[] Buffer_Receive_GlobElem;
      
      delete[] Buffer_Receive_Triangle_presence[iDomain];
      delete[] Buffer_Receive_Quadrilateral_presence[iDomain];
      delete[] Buffer_Receive_Tetrahedron_presence[iDomain];
      delete[] Buffer_Receive_Hexahedron_presence[iDomain];
      delete[] Buffer_Receive_Prism_presence[iDomain];
      delete[] Buffer_Receive_Pyramid_presence[iDomain];
      
#endif
      
    } else {
      
      /*--- Store the element data from our local rank ---*/
      
      iElemRecv = 0;
      
      for (iElemTriangle = 0; iElemTriangle < nElemTriangle_r[iDomain]; iElemTriangle++) {
        map<unsigned long, bool>::const_iterator MI = Triangle_presence.find(Buffer_Receive_Triangle_presence_loc[iElemTriangle]);
        if (MI == Triangle_presence.end()) {
          Triangle_presence[Buffer_Receive_Triangle_presence_loc[iElemTriangle]] = true;
          elem[iElem] = new CTriangle(Global_to_local_Point_recv[Buffer_Receive_Triangle_loc[iElemTriangle*3+0]],
                                      Global_to_local_Point_recv[Buffer_Receive_Triangle_loc[iElemTriangle*3+1]],
                                      Global_to_local_Point_recv[Buffer_Receive_Triangle_loc[iElemTriangle*3+2]], 2);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem_loc[iElemRecv]);
          iElem++; iElemTria++;
        }
        iElemRecv++;
      }
      
      for (iElemQuadrilateral = 0; iElemQuadrilateral < nElemQuadrilateral_r[iDomain]; iElemQuadrilateral++) {
        map<unsigned long, bool>::const_iterator MI = Quadrilateral_presence.find(Buffer_Receive_Quadrilateral_presence_loc[iElemQuadrilateral]);
        if (MI == Quadrilateral_presence.end()) {
          Quadrilateral_presence[Buffer_Receive_Quadrilateral_presence_loc[iElemQuadrilateral]] = true;
          elem[iElem] = new CQuadrilateral(Global_to_local_Point_recv[Buffer_Receive_Quadrilateral_loc[iElemQuadrilateral*4+0]],
                                           Global_to_local_Point_recv[Buffer_Receive_Quadrilateral_loc[iElemQuadrilateral*4+1]],
                                           Global_to_local_Point_recv[Buffer_Receive_Quadrilateral_loc[iElemQuadrilateral*4+2]],
                                           Global_to_local_Point_recv[Buffer_Receive_Quadrilateral_loc[iElemQuadrilateral*4+3]], 2);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem_loc[iElemRecv]);
          iElem++; iElemRect++;
        }
        iElemRecv++;
      }
      
      for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron_r[iDomain]; iElemTetrahedron++) {
        map<unsigned long, bool>::const_iterator MI = Tetrahedron_presence.find(Buffer_Receive_Tetrahedron_presence_loc[iElemTetrahedron]);
        if (MI == Tetrahedron_presence.end()) {
          Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence_loc[iElemTetrahedron]] = true;
          elem[iElem] = new CTetrahedron(Global_to_local_Point_recv[Buffer_Receive_Tetrahedron_loc[iElemTetrahedron*4+0]],
                                         Global_to_local_Point_recv[Buffer_Receive_Tetrahedron_loc[iElemTetrahedron*4+1]],
                                         Global_to_local_Point_recv[Buffer_Receive_Tetrahedron_loc[iElemTetrahedron*4+2]],
                                         Global_to_local_Point_recv[Buffer_Receive_Tetrahedron_loc[iElemTetrahedron*4+3]]);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem_loc[iElemRecv]);
          iElem++; iElemTetr++;
        }
        iElemRecv++;
      }
      
      for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron_r[iDomain]; iElemHexahedron++) {
        map<unsigned long, bool>::const_iterator MI = Hexahedron_presence.find(Buffer_Receive_Hexahedron_presence_loc[iElemHexahedron]);
        if (MI == Hexahedron_presence.end()) {
          Hexahedron_presence[Buffer_Receive_Hexahedron_presence_loc[iElemHexahedron]] = true;
          elem[iElem] = new CHexahedron(Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron*8+0]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron*8+1]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron*8+2]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron*8+3]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron*8+4]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron*8+5]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron*8+6]],
                                        Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron*8+7]]);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem_loc[iElemRecv]);
          iElem++; iElemHexa++;
        }
        iElemRecv++;
      }
      
      for (iElemPrism = 0; iElemPrism < nElemPrism_r[iDomain]; iElemPrism++) {
        map<unsigned long, bool>::const_iterator MI = Prism_presence.find(Buffer_Receive_Prism_presence_loc[iElemPrism]);
        if (MI == Prism_presence.end()) {
          Prism_presence[Buffer_Receive_Prism_presence_loc[iElemPrism]] = true;
          elem[iElem] = new CPrism(Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism*6+0]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism*6+1]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism*6+2]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism*6+3]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism*6+4]],
                                   Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism*6+5]]);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem_loc[iElemRecv]);
          iElem++; iElemPris++;
        }
        iElemRecv++;
      }
      
      for (iElemPyramid = 0; iElemPyramid < nElemPyramid_r[iDomain]; iElemPyramid++) {
        map<unsigned long, bool>::const_iterator MI = Pyramid_presence.find(Buffer_Receive_Pyramid_presence_loc[iElemPyramid]);
        if (MI == Pyramid_presence.end()) {
          Pyramid_presence[Buffer_Receive_Pyramid_presence_loc[iElemPyramid]] = true;
          elem[iElem] = new CPyramid(Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid*5+0]],
                                     Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid*5+1]],
                                     Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid*5+2]],
                                     Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid*5+3]],
                                     Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid*5+4]]);
          elem[iElem]->SetGlobalIndex(Buffer_Receive_GlobElem_loc[iElemRecv]);
          iElem++; iElemPyra++;
        }
        iElemRecv++;
      }
      
      /*--- Free memory for element data ---*/
      
      delete[] Buffer_Receive_Triangle_loc;
      delete[] Buffer_Receive_Quadrilateral_loc;
      delete[] Buffer_Receive_Tetrahedron_loc;
      delete[] Buffer_Receive_Hexahedron_loc;
      delete[] Buffer_Receive_Prism_loc;
      delete[] Buffer_Receive_Pyramid_loc;
      delete[] Buffer_Receive_GlobElem_loc;
      
      delete[] Buffer_Receive_Triangle_presence_loc;
      delete[] Buffer_Receive_Quadrilateral_presence_loc;
      delete[] Buffer_Receive_Tetrahedron_presence_loc;
      delete[] Buffer_Receive_Hexahedron_presence_loc;
      delete[] Buffer_Receive_Prism_presence_loc;
      delete[] Buffer_Receive_Pyramid_presence_loc;
      
    }
  }
  
#ifdef HAVE_MPI
  for (iDomain = 0; iDomain < (unsigned long)size; iDomain++) {
    if ((unsigned long)rank != iDomain) SU2_MPI::Waitall(16, send_req, send_stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete [] Buffer_Send_Coord;
  delete [] Buffer_Send_GlobalPointIndex;
  delete [] Buffer_Send_Color;
  delete [] Buffer_Send_Triangle;
  delete [] Buffer_Send_Quadrilateral;
  delete [] Buffer_Send_Tetrahedron;
  delete [] Buffer_Send_Hexahedron;
  delete [] Buffer_Send_Prism;
  delete [] Buffer_Send_Pyramid;
  delete [] Buffer_Send_GlobElem;
  delete [] Buffer_Send_BoundLine;
  delete [] Buffer_Send_BoundTriangle;
  delete [] Buffer_Send_BoundQuadrilateral;
  delete [] Buffer_Send_Local2Global_Marker;
  
  delete [] Buffer_Send_SendDomain_Periodic;
  delete [] Buffer_Send_SendDomain_PeriodicTrans;
  delete [] Buffer_Send_SendDomain_PeriodicReceptor;
  delete [] Buffer_Send_ReceivedDomain_Periodic;
  delete [] Buffer_Send_ReceivedDomain_PeriodicTrans;
  delete [] Buffer_Send_ReceivedDomain_PeriodicDonor;
  
#ifdef HAVE_MPI
  delete [] Buffer_Receive_Triangle_presence;
  delete [] Buffer_Receive_Quadrilateral_presence;
  delete [] Buffer_Receive_Tetrahedron_presence;
  delete [] Buffer_Receive_Hexahedron_presence;
  delete [] Buffer_Receive_Prism_presence;
  delete [] Buffer_Receive_Pyramid_presence;
#endif
  
  delete [] Local_to_global_Triangle;
  delete [] Local_to_global_Quadrilateral;
  delete [] Local_to_global_Tetrahedron;
  delete [] Local_to_global_Hexahedron;
  delete [] Local_to_global_Prism;
  delete [] Local_to_global_Pyramid;
  
  
  /*--- Communicate the number of each element type to all processors. These
   values are important for merging and writing output later. ---*/
  
#ifdef HAVE_MPI
  unsigned long Local_nElem = nElem;
  SU2_MPI::Allreduce(&Local_nElem, &Global_nElem, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nElem = nElem;
#endif
  
  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << Global_nElem << " interior elements including halo cells. " << endl;
  
  /*--- Set the value of Global_nElemDomain (stored in the geometry container that is passed in) ---*/

  Global_nElemDomain = geometry->GetGlobal_nElemDomain();

  /*--- Store total number of each element type after incrementing the
   counters in the recv loop above (to make sure there aren't repeats). ---*/
  
  nelem_triangle = iElemTria;
  nelem_quad     = iElemRect;
  nelem_tetra    = iElemTetr;
  nelem_hexa     = iElemHexa;
  nelem_prism    = iElemPris;
  nelem_pyramid  = iElemPyra;
  
#ifdef HAVE_MPI
  unsigned long Local_nElemTri     = nelem_triangle;
  unsigned long Local_nElemQuad    = nelem_quad;
  unsigned long Local_nElemTet     = nelem_tetra;
  unsigned long Local_nElemHex     = nelem_hexa;
  unsigned long Local_nElemPrism   = nelem_prism;
  unsigned long Local_nElemPyramid = nelem_pyramid;
  SU2_MPI::Allreduce(&Local_nElemTri, &Global_nelem_triangle, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemQuad, &Global_nelem_quad, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemTet, &Global_nelem_tetra, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemHex, &Global_nelem_hexa, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemPrism, &Global_nelem_prism, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemPyramid, &Global_nelem_pyramid, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nelem_triangle = nelem_triangle;
  Global_nelem_quad     = nelem_quad;
  Global_nelem_tetra    = nelem_tetra;
  Global_nelem_hexa     = nelem_hexa;
  Global_nelem_prism    = nelem_prism;
  Global_nelem_pyramid  = nelem_pyramid;
#endif
  
  /*--- Print information about the elements to the console ---*/
  
  if (rank == MASTER_NODE) {
    if (Global_nelem_triangle > 0)  cout << Global_nelem_triangle << " triangles."      << endl;
    if (Global_nelem_quad > 0)      cout << Global_nelem_quad     << " quadrilaterals." << endl;
    if (Global_nelem_tetra > 0)     cout << Global_nelem_tetra    << " tetrahedra."     << endl;
    if (Global_nelem_hexa > 0)      cout << Global_nelem_hexa     << " hexahedra."      << endl;
    if (Global_nelem_prism > 0)     cout << Global_nelem_prism    << " prisms."         << endl;
    if (Global_nelem_pyramid > 0)   cout << Global_nelem_pyramid  << " pyramids."       << endl;
  }
  
  /*--- Now partition the boundary elements on the markers. Note that, for
   now, we are still performing the boundary partitioning using the master
   node alone. The boundaries should make up a much smaller portion of the
   mesh, so this is ok for now, but we will transition to a parallel version
   of this soon that follows the same procedure above for the interior. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Create auxiliary vectors based on the original geometry ---*/
    
    MarkerIn = new bool[geometry->GetnMarker()];
    VertexIn = new bool*[geometry->GetnMarker()];
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
      VertexIn[iMarker] = new bool[geometry->GetnElem_Bound(iMarker)];
    
    Buffer_Send_nDim      = geometry->GetnDim();
    Buffer_Send_nZone     = geometry->GetnZone();
    Buffer_Send_nPeriodic = config->GetnPeriodicIndex();
    Buffer_Send_Center    = new su2double[Buffer_Send_nPeriodic*3];
    Buffer_Send_Rotation  = new su2double[Buffer_Send_nPeriodic*3];
    Buffer_Send_Translate = new su2double[Buffer_Send_nPeriodic*3];
    
    Buffer_Send_nSendDomain_Periodic     = new unsigned long[nDomain];
    Buffer_Send_nReceivedDomain_Periodic = new unsigned long[nDomain];
    
    /*--- Create a local copy of config->GetMarker_All_SendRecv and
     config->GetMarker_All_TagBound in the master node ---*/
    
    Marker_All_SendRecv_Copy = new short[geometry->GetnMarker()];
    Marker_All_TagBound_Copy = new string[geometry->GetnMarker()];
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      Marker_All_SendRecv_Copy[iMarker] = config->GetMarker_All_SendRecv(iMarker);
      Marker_All_TagBound_Copy[iMarker] = config->GetMarker_All_TagBound(iMarker);
    }
    
  }
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    if (rank == MASTER_NODE) {
      
      /*--- Interior dimensionalization. Loop over the original grid
       to perform the dimensionalizaton of the domain variables ---*/
      
      //      Buffer_Send_nElemTotal         = 0;
      //      Buffer_Send_nPointTotal        = 0;
      //      Buffer_Send_nPointGhost        = 0;
      //      Buffer_Send_nPointDomainTotal  = 0;
      //      Buffer_Send_nPointPeriodic     = 0;
      //      Buffer_Send_nElemTriangle      = 0;
      //      Buffer_Send_nElemQuadrilateral = 0;
      //      Buffer_Send_nElemTetrahedron   = 0;
      //      Buffer_Send_nElemHexahedron    = 0;
      //      Buffer_Send_nElemPrism         = 0;
      //      Buffer_Send_nElemPyramid       = 0;
      
      /*--- Boundary dimensionalization. Dimensionalization with physical
       boundaries, compute Buffer_Send_nMarkerDomain,
       Buffer_Send_nVertexDomain[nMarkerDomain] ---*/
      
      Buffer_Send_nMarkerDomain        = 0;
      Buffer_Send_nBoundLineTotal      = 0;
      Buffer_Send_nBoundTriangleTotal  = 0;
      Buffer_Send_nBoundQuadrilateralTotal = 0;
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        Buffer_Send_nVertexDomain[iMarker]   = 0;
        Buffer_Send_nBoundLine[iMarker]      = 0;
        Buffer_Send_nBoundTriangle[iMarker]  = 0;
        Buffer_Send_nBoundQuadrilateral[iMarker] = 0;
        Buffer_Send_Marker_All_SendRecv[iMarker] = Marker_All_SendRecv_Copy[iMarker];
        SPRINTF(&Buffer_Send_Marker_All_TagBound[iMarker*MAX_STRING_SIZE], "%s",
                Marker_All_TagBound_Copy[iMarker].c_str());
      }
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
          
          MarkerIn[iMarker] = false;
          Buffer_Send_nVertexDomain[Buffer_Send_nMarkerDomain] = 0;
          
          for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
            VertexIn[iMarker][iVertex] = false;
            for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++) {
              iPoint = geometry->bound[iMarker][iVertex]->GetNode(iNode);
              if (local_colour_values[iPoint] == iDomain) VertexIn[iMarker][iVertex] = true;
            }
            
            /*--- If this vertex should be sent, increment the element type ---*/
            if (VertexIn[iMarker][iVertex]) {
              switch(geometry->bound[iMarker][iVertex]->GetVTK_Type()) {
                case LINE:
                  Buffer_Send_nBoundLine[Buffer_Send_nMarkerDomain]++;
                  Buffer_Send_nBoundLineTotal++;
                  break;
                case TRIANGLE:
                  Buffer_Send_nBoundTriangle[Buffer_Send_nMarkerDomain]++;
                  Buffer_Send_nBoundTriangleTotal++;
                  break;
                case QUADRILATERAL:
                  Buffer_Send_nBoundQuadrilateral[Buffer_Send_nMarkerDomain]++;
                  Buffer_Send_nBoundQuadrilateralTotal++;
                  break;
              }
              
              /*--- Increment the total number of vertices to be sent ---*/
              Buffer_Send_nVertexDomain[Buffer_Send_nMarkerDomain]++;
              MarkerIn[iMarker] = true;
              
            }
          }
          
          /*--- Increment the number of markers to be sent ---*/
          if (MarkerIn[iMarker]) { Buffer_Send_nMarkerDomain++; }
          
        }
      }
      
      /*--- Copy periodic information from the config file ---*/
      
      for (iPeriodic = 0; iPeriodic < Buffer_Send_nPeriodic; iPeriodic++) {
        for (iDim = 0; iDim < 3; iDim++) {
          Buffer_Send_Center[iDim+iPeriodic*3]    = config->GetPeriodicCenter(iPeriodic)[iDim];
          Buffer_Send_Rotation[iDim+iPeriodic*3]  = config->GetPeriodicRotation(iPeriodic)[iDim];
          Buffer_Send_Translate[iDim+iPeriodic*3] = config->GetPeriodicTranslate(iPeriodic)[iDim];
        }
      }
      
      /*--- Dimensionalization of the periodic auxiliary vectors ---*/
      
      for (jDomain = 0; jDomain < nDomain; jDomain++) {
        Buffer_Send_nSendDomain_Periodic[jDomain]     = 0;
        Buffer_Send_nReceivedDomain_Periodic[jDomain] = 0;
      }
      Buffer_Send_nTotalSendDomain_Periodic     = 0;
      Buffer_Send_nTotalReceivedDomain_Periodic = 0;
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
          for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
            iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
            if (iDomain == local_colour_values[iPoint]) {
              
              if (config->GetMarker_All_SendRecv(iMarker) > 0) {
                
                /*--- Identify the color of the receptor ---*/
                
                for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                  if ((config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) &&
                      (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker))) {
                    jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                    ReceptorColor = local_colour_values[jPoint];
                  }
                }
                
                Buffer_Send_nSendDomain_Periodic[ReceptorColor]++;
                Buffer_Send_nTotalSendDomain_Periodic++;
                
              }
              if (config->GetMarker_All_SendRecv(iMarker) < 0) {
                
                /*--- Identify the color of the donor ---*/
                
                for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                  if ((config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) &&
                      (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker))) {
                    jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                    DonorColor = local_colour_values[jPoint];
                  }
                }
                
                Buffer_Send_nReceivedDomain_Periodic[DonorColor]++;
                Buffer_Send_nTotalReceivedDomain_Periodic++;
                
              }
            }
          }
        }
      }
      
      /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
      
      Buffer_Send_BoundLine           = new unsigned long[Buffer_Send_nBoundLineTotal*N_POINTS_LINE];
      Buffer_Send_BoundTriangle       = new unsigned long[Buffer_Send_nBoundTriangleTotal*N_POINTS_TRIANGLE];
      Buffer_Send_BoundQuadrilateral  = new unsigned long[Buffer_Send_nBoundQuadrilateralTotal*N_POINTS_QUADRILATERAL];
      Buffer_Send_Local2Global_Marker = new unsigned long[Buffer_Send_nMarkerDomain];
      
      Buffer_Send_SendDomain_Periodic           = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
      Buffer_Send_SendDomain_PeriodicTrans      = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
      Buffer_Send_SendDomain_PeriodicReceptor   = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
      Buffer_Send_ReceivedDomain_Periodic       = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];
      Buffer_Send_ReceivedDomain_PeriodicTrans  = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];
      Buffer_Send_ReceivedDomain_PeriodicDonor  = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];
      
      if (iDomain != (unsigned long)MASTER_NODE) {
        
#ifdef HAVE_MPI
        
        SU2_MPI::Isend(&Buffer_Send_nBoundLineTotal, 1,
                       MPI_UNSIGNED_LONG, iDomain,
                       0, MPI_COMM_WORLD, &send_req[0]);
        
        SU2_MPI::Isend(&Buffer_Send_nBoundTriangleTotal, 1,
                       MPI_UNSIGNED_LONG, iDomain,
                       1, MPI_COMM_WORLD, &send_req[1]);
        
        SU2_MPI::Isend(&Buffer_Send_nBoundQuadrilateralTotal, 1,
                       MPI_UNSIGNED_LONG,  iDomain,
                       2, MPI_COMM_WORLD, &send_req[2]);
        
        SU2_MPI::Isend(&Buffer_Send_nMarkerDomain, 1,
                       MPI_UNSIGNED_SHORT, iDomain,
                       3, MPI_COMM_WORLD, &send_req[3]);
        
        SU2_MPI::Isend(Buffer_Send_nVertexDomain,
                       nMarker_Max, MPI_UNSIGNED_LONG, iDomain,
                       4, MPI_COMM_WORLD, &send_req[4]);
        
        SU2_MPI::Isend(Buffer_Send_nBoundLine,
                       nMarker_Max, MPI_UNSIGNED_LONG, iDomain,
                       5, MPI_COMM_WORLD, &send_req[5]);
        
        SU2_MPI::Isend(Buffer_Send_nBoundTriangle,
                       nMarker_Max, MPI_UNSIGNED_LONG, iDomain,
                       6, MPI_COMM_WORLD, &send_req[6]);
        
        SU2_MPI::Isend(Buffer_Send_nBoundQuadrilateral,
                       nMarker_Max, MPI_UNSIGNED_LONG, iDomain,
                       7, MPI_COMM_WORLD, &send_req[7]);
        
        SU2_MPI::Isend(Buffer_Send_Marker_All_SendRecv,
                       nMarker_Max, MPI_SHORT, iDomain,
                       8, MPI_COMM_WORLD, &send_req[8]);
        
        SU2_MPI::Isend(Buffer_Send_Marker_All_TagBound,
                       nMarker_Max*MAX_STRING_SIZE, MPI_CHAR, iDomain,
                       9, MPI_COMM_WORLD, &send_req[9]);
        
        SU2_MPI::Isend(&Buffer_Send_nPeriodic,
                       1, MPI_UNSIGNED_SHORT, iDomain,
                       10, MPI_COMM_WORLD, &send_req[10]);
        
        SU2_MPI::Isend(Buffer_Send_Center,
                       nPeriodic*3, MPI_DOUBLE, iDomain,
                       11, MPI_COMM_WORLD, &send_req[11]);
        
        SU2_MPI::Isend(Buffer_Send_Rotation,
                       nPeriodic*3, MPI_DOUBLE, iDomain,
                       12, MPI_COMM_WORLD, &send_req[12]);
        
        SU2_MPI::Isend(Buffer_Send_Translate,
                       nPeriodic*3, MPI_DOUBLE, iDomain,
                       13, MPI_COMM_WORLD, &send_req[13]);
        
        SU2_MPI::Isend(&Buffer_Send_nTotalSendDomain_Periodic,
                       1, MPI_UNSIGNED_LONG, iDomain,
                       14, MPI_COMM_WORLD, &send_req[14]);
        
        SU2_MPI::Isend(&Buffer_Send_nTotalReceivedDomain_Periodic,
                       1, MPI_UNSIGNED_LONG, iDomain,
                       15, MPI_COMM_WORLD, &send_req[15]);
        
        SU2_MPI::Isend(Buffer_Send_nSendDomain_Periodic,
                       nDomain, MPI_UNSIGNED_LONG, iDomain,
                       16, MPI_COMM_WORLD, &send_req[16]);
        
        SU2_MPI::Isend(Buffer_Send_nReceivedDomain_Periodic,
                       nDomain, MPI_UNSIGNED_LONG, iDomain,
                       17, MPI_COMM_WORLD, &send_req[17]);
        
        /*--- Wait for this set of non-blocking comm. to complete ---*/
        
        SU2_MPI::Waitall(18, send_req, send_stat);
        
#endif
        
      } else {
        
        /*--- We are the master node, so simply copy values into place ---*/
        
        nDim  = Buffer_Send_nDim;
        nZone = Buffer_Send_nZone;
        
        nPeriodic      = Buffer_Send_nPeriodic;
        //        nPointGhost    = Buffer_Send_nPointGhost;
        //        nPointPeriodic = Buffer_Send_nPointPeriodic;
        
        nBoundLineTotal      = Buffer_Send_nBoundLineTotal;
        nBoundTriangleTotal  = Buffer_Send_nBoundTriangleTotal;
        nBoundQuadrilateralTotal = Buffer_Send_nBoundQuadrilateralTotal;
        nMarkerDomain        = Buffer_Send_nMarkerDomain;
        
        for (iMarker = 0; iMarker < nMarker_Max; iMarker++) {
          nVertexDomain[iMarker] = Buffer_Send_nVertexDomain[iMarker];
          nBoundLine[iMarker] = Buffer_Send_nBoundLine[iMarker];
          nBoundTriangle[iMarker] = Buffer_Send_nBoundTriangle[iMarker];
          nBoundQuadrilateral[iMarker] = Buffer_Send_nBoundQuadrilateral[iMarker];
          Marker_All_SendRecv[iMarker] = Buffer_Send_Marker_All_SendRecv[iMarker];
          for (iter = 0; iter < MAX_STRING_SIZE; iter++)
            Marker_All_TagBound[iMarker*MAX_STRING_SIZE+iter] = Buffer_Send_Marker_All_TagBound[iMarker*MAX_STRING_SIZE+iter];
        }
        
        Buffer_Receive_Center    = new su2double[nPeriodic*3];
        Buffer_Receive_Rotation  = new su2double[nPeriodic*3];
        Buffer_Receive_Translate = new su2double[nPeriodic*3];
        
        for (iter = 0; iter < nPeriodic*3; iter++) {
          Buffer_Receive_Center[iter]    =  Buffer_Send_Center[iter];
          Buffer_Receive_Rotation[iter]  =  Buffer_Send_Rotation[iter];
          Buffer_Receive_Translate[iter] =  Buffer_Send_Translate[iter];
        }
        
        nTotalSendDomain_Periodic     = Buffer_Send_nTotalSendDomain_Periodic;
        nTotalReceivedDomain_Periodic = Buffer_Send_nTotalReceivedDomain_Periodic;
        
        for (iter = 0; iter < nDomain; iter++) {
          nSendDomain_Periodic[iter] = Buffer_Send_nSendDomain_Periodic[iter];
          nReceivedDomain_Periodic[iter] = Buffer_Send_nReceivedDomain_Periodic[iter];
        }
        
      }
    }
    
    /*--- Each rank now begins to receive information from the master ---*/
    
    if ((unsigned long)rank == iDomain) {
      
      /*--- First, receive the size of buffers before receiving the data ---*/
      
      if (rank != MASTER_NODE) {
        
#ifdef HAVE_MPI
        
        MPI_Probe(MASTER_NODE, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(&nBoundLineTotal, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 0, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(&nBoundTriangleTotal, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 1, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 2, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(&nBoundQuadrilateralTotal, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 2, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 3, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_SHORT, &recv_count);
        SU2_MPI::Recv(&nMarkerDomain, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 3, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 4, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(nVertexDomain, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 4, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 5, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(nBoundLine, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 5, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 6, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(nBoundTriangle, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 6, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 7, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(nBoundQuadrilateral, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 7, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 8, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_SHORT, &recv_count);
        SU2_MPI::Recv(Marker_All_SendRecv, recv_count, MPI_SHORT,
                      MASTER_NODE, 8, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 9, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &recv_count);
        SU2_MPI::Recv(Marker_All_TagBound, recv_count, MPI_CHAR,
                      MASTER_NODE, 9, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 10, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_SHORT, &recv_count);
        SU2_MPI::Recv(&nPeriodic, recv_count, MPI_UNSIGNED_SHORT,
                      MASTER_NODE, 10, MPI_COMM_WORLD, &status);
        
#endif
        
        /*--- Marker_All_TagBound and Marker_All_SendRecv, set the same
         values in the config files of all the files ---*/
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          config->SetMarker_All_SendRecv(iMarker,
                                         Marker_All_SendRecv[iMarker]);
          config->SetMarker_All_TagBound(iMarker,
                                         string(&Marker_All_TagBound[iMarker*MAX_STRING_SIZE]));
        }
        
        
        /*--- Periodic boundary conditions ---*/
        
        Buffer_Receive_Center    = new su2double[nPeriodic*3];
        Buffer_Receive_Rotation  = new su2double[nPeriodic*3];
        Buffer_Receive_Translate = new su2double[nPeriodic*3];
        
#ifdef HAVE_MPI
        
        MPI_Probe(MASTER_NODE, 11, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_Center, recv_count, MPI_DOUBLE,
                      MASTER_NODE, 11, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 12, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_Rotation, recv_count, MPI_DOUBLE,
                      MASTER_NODE, 12, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 13, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_Translate, recv_count, MPI_DOUBLE,
                      MASTER_NODE, 13, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 14, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(&nTotalSendDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 14, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 15, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(&nTotalReceivedDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 15, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 16, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(nSendDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 16, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 17, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(nReceivedDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 17, MPI_COMM_WORLD, &status);
        
#endif
        
        config->SetnPeriodicIndex(nPeriodic);
        
        for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
          
          su2double* center    = new su2double[3];
          su2double* rotation  = new su2double[3];
          su2double* translate = new su2double[3];
          
          for (iDim = 0; iDim < 3; iDim++) {
            center[iDim] = Buffer_Receive_Center[iDim+iPeriodic*3];
            rotation[iDim] = Buffer_Receive_Rotation[iDim+iPeriodic*3];
            translate[iDim] = Buffer_Receive_Translate[iDim+iPeriodic*3];
          }
          config->SetPeriodicCenter(iPeriodic, center);
          config->SetPeriodicRotation(iPeriodic, rotation);
          config->SetPeriodicTranslate(iPeriodic, translate);
        
          delete [] center; delete [] rotation; delete [] translate;
          
        }
        
      }
      
      delete [] Buffer_Receive_Center;
      delete [] Buffer_Receive_Rotation;
      delete [] Buffer_Receive_Translate;
      
      /*--- Allocate the receive buffer vector ---*/
      
      Buffer_Receive_BoundLine           = new unsigned long[nBoundLineTotal*2];
      Buffer_Receive_BoundTriangle       = new unsigned long[nBoundTriangleTotal*3];
      Buffer_Receive_BoundQuadrilateral      = new unsigned long[nBoundQuadrilateralTotal*4];
      Buffer_Receive_Local2Global_Marker = new unsigned long[nMarkerDomain];
      
      Buffer_Receive_SendDomain_Periodic          = new unsigned long[nTotalSendDomain_Periodic];
      Buffer_Receive_SendDomain_PeriodicTrans     = new unsigned long[nTotalSendDomain_Periodic];
      Buffer_Receive_SendDomain_PeriodicReceptor  = new unsigned long[nTotalSendDomain_Periodic];
      Buffer_Receive_ReceivedDomain_Periodic      = new unsigned long[nTotalReceivedDomain_Periodic];
      Buffer_Receive_ReceivedDomain_PeriodicTrans = new unsigned long[nTotalReceivedDomain_Periodic];
      Buffer_Receive_ReceivedDomain_PeriodicDonor = new unsigned long[nTotalReceivedDomain_Periodic];
      
    }
    
    /*--- Set the value of the Send buffers ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- Set the value of the boundary geometry ---*/
      
      iMarkerDomain = 0;
      iBoundLineTotal = 0; iBoundTriangleTotal = 0; iBoundQuadrilateralTotal = 0;
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) && (MarkerIn[iMarker])) {
          for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
            
            if (VertexIn[iMarker][iVertex]) {
              
              /*--- Send global index here and then convert to local on the recv ---*/
              
              for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++) {
                vnodes_local[iNode] = geometry->bound[iMarker][iVertex]->GetNode(iNode);
              }
              
              switch(geometry->bound[iMarker][iVertex]->GetVTK_Type()) {
                case LINE:
                  Buffer_Send_BoundLine[N_POINTS_LINE*iBoundLineTotal+0] = vnodes_local[0];
                  Buffer_Send_BoundLine[N_POINTS_LINE*iBoundLineTotal+1] = vnodes_local[1];
                  iBoundLineTotal++;
                  break;
                case TRIANGLE:
                  Buffer_Send_BoundTriangle[N_POINTS_TRIANGLE*iBoundTriangleTotal+0] = vnodes_local[0];
                  Buffer_Send_BoundTriangle[N_POINTS_TRIANGLE*iBoundTriangleTotal+1] = vnodes_local[1];
                  Buffer_Send_BoundTriangle[N_POINTS_TRIANGLE*iBoundTriangleTotal+2] = vnodes_local[2];
                  iBoundTriangleTotal++;
                  break;
                case QUADRILATERAL:
                  Buffer_Send_BoundQuadrilateral[N_POINTS_QUADRILATERAL*iBoundQuadrilateralTotal+0] = vnodes_local[0];
                  Buffer_Send_BoundQuadrilateral[N_POINTS_QUADRILATERAL*iBoundQuadrilateralTotal+1] = vnodes_local[1];
                  Buffer_Send_BoundQuadrilateral[N_POINTS_QUADRILATERAL*iBoundQuadrilateralTotal+2] = vnodes_local[2];
                  Buffer_Send_BoundQuadrilateral[N_POINTS_QUADRILATERAL*iBoundQuadrilateralTotal+3] = vnodes_local[3];
                  iBoundQuadrilateralTotal++;
                  break;
              }
            }
          }
          
          Buffer_Send_Local2Global_Marker[iMarkerDomain] = iMarker;
          iMarkerDomain++;
          
        }
      }
      
      /*--- Evaluate the number of already existing periodic boundary conditions ---*/
      
      iTotalSendDomain_Periodic = 0;
      iTotalReceivedDomain_Periodic = 0;
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        
        if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
          
          for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
            
            iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
            Transformation = geometry->bound[iMarker][iVertex]->GetRotation_Type();
            
            if (iDomain == local_colour_values[iPoint]) {
              
              /*--- If the information is going to be sended, find the
               domain of the receptor ---*/
              
              if (config->GetMarker_All_SendRecv(iMarker) > 0) {
                
                /*--- Identify the color of the receptor ---*/
                
                for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                  if ((config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) &&
                      (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker))) {
                    jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                    ReceptorColor = local_colour_values[jPoint];
                  }
                }
                
                /*--- For each color of the receptor we will han an extra marker (+) ---*/
                
                Buffer_Send_SendDomain_Periodic[iTotalSendDomain_Periodic] = iPoint;
                Buffer_Send_SendDomain_PeriodicTrans[iTotalSendDomain_Periodic] = Transformation;
                Buffer_Send_SendDomain_PeriodicReceptor[iTotalSendDomain_Periodic] = ReceptorColor;
                
                iTotalSendDomain_Periodic++;
                
              }
              
              /*--- If the information is goint to be received, find the domain if the donor ---*/
              
              if (config->GetMarker_All_SendRecv(iMarker) < 0) {
                
                /*--- Identify the color of the donor ---*/
                
                for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                  if ((config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) &&
                      (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker) )) {
                    jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                    DonorColor = local_colour_values[jPoint];
                  }
                }
                
                /*--- For each color of the donor we will han an extra marker (-) ---*/
                
                Buffer_Send_ReceivedDomain_Periodic[iTotalReceivedDomain_Periodic] = iPoint;
                Buffer_Send_ReceivedDomain_PeriodicTrans[iTotalReceivedDomain_Periodic] = Transformation;
                Buffer_Send_ReceivedDomain_PeriodicDonor[iTotalReceivedDomain_Periodic] = DonorColor;
                
                iTotalReceivedDomain_Periodic++;
                
              }
            }
          }
        }
      }
      
      /*--- Send the buffers with the geometrical information ---*/
      
      if (iDomain != (unsigned long)MASTER_NODE) {
        
#ifdef HAVE_MPI
        
        SU2_MPI::Isend(Buffer_Send_BoundLine,
                       Buffer_Send_nBoundLineTotal*N_POINTS_LINE, MPI_UNSIGNED_LONG, iDomain,
                       0, MPI_COMM_WORLD, &send_req[0]);
        
        SU2_MPI::Isend(Buffer_Send_BoundTriangle,
                       Buffer_Send_nBoundTriangleTotal*N_POINTS_TRIANGLE, MPI_UNSIGNED_LONG, iDomain,
                       1, MPI_COMM_WORLD, &send_req[1]);
        
        SU2_MPI::Isend(Buffer_Send_BoundQuadrilateral,
                       Buffer_Send_nBoundQuadrilateralTotal*N_POINTS_QUADRILATERAL, MPI_UNSIGNED_LONG, iDomain,
                       2, MPI_COMM_WORLD, &send_req[2]);
        
        SU2_MPI::Isend(Buffer_Send_Local2Global_Marker,
                       Buffer_Send_nMarkerDomain, MPI_UNSIGNED_LONG, iDomain,
                       3, MPI_COMM_WORLD, &send_req[3]);
        
        SU2_MPI::Isend(Buffer_Send_SendDomain_Periodic,
                       Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                       4, MPI_COMM_WORLD, &send_req[4]);
        
        SU2_MPI::Isend(Buffer_Send_SendDomain_PeriodicTrans,
                       Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                       5, MPI_COMM_WORLD, &send_req[5]);
        
        SU2_MPI::Isend(Buffer_Send_SendDomain_PeriodicReceptor,
                       Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                       6, MPI_COMM_WORLD, &send_req[6]);
        
        SU2_MPI::Isend(Buffer_Send_ReceivedDomain_Periodic,
                       Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                       7, MPI_COMM_WORLD, &send_req[7]);
        
        SU2_MPI::Isend(Buffer_Send_ReceivedDomain_PeriodicTrans,
                       Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                       8, MPI_COMM_WORLD, &send_req[8]);
        
        SU2_MPI::Isend(Buffer_Send_ReceivedDomain_PeriodicDonor,
                       Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                       9, MPI_COMM_WORLD, &send_req[9]);
        
        /*--- Wait for this set of non-blocking comm. to complete ---*/
        
        SU2_MPI::Waitall(10, send_req, send_stat);
        
#endif
        
      } else {
        
        /*--- Copy the data directly from our own rank ---*/
        
        for (iter = 0; iter < Buffer_Send_nBoundLineTotal*N_POINTS_LINE; iter++)
          Buffer_Receive_BoundLine[iter] =  Buffer_Send_BoundLine[iter];
        
        for (iter = 0; iter < Buffer_Send_nBoundTriangleTotal*N_POINTS_TRIANGLE; iter++)
          Buffer_Receive_BoundTriangle[iter] =  Buffer_Send_BoundTriangle[iter];
        
        for (iter = 0; iter < Buffer_Send_nBoundQuadrilateralTotal*N_POINTS_QUADRILATERAL; iter++)
          Buffer_Receive_BoundQuadrilateral[iter] =  Buffer_Send_BoundQuadrilateral[iter];
        
        for (iter = 0; iter < Buffer_Send_nMarkerDomain; iter++)
          Buffer_Receive_Local2Global_Marker[iter] =  Buffer_Send_Local2Global_Marker[iter];
        
        for (iter = 0; iter < Buffer_Send_nTotalSendDomain_Periodic; iter++) {
          Buffer_Receive_SendDomain_Periodic[iter] = Buffer_Send_SendDomain_Periodic[iter];
          Buffer_Receive_SendDomain_PeriodicTrans[iter] = Buffer_Send_SendDomain_PeriodicTrans[iter];
          Buffer_Receive_SendDomain_PeriodicReceptor[iter] = Buffer_Send_SendDomain_PeriodicReceptor[iter];
        }
        
        for (iter = 0; iter < Buffer_Send_nTotalReceivedDomain_Periodic; iter++) {
          Buffer_Receive_ReceivedDomain_Periodic[iter] = Buffer_Send_ReceivedDomain_Periodic[iter];
          Buffer_Receive_ReceivedDomain_PeriodicTrans[iter] = Buffer_Send_ReceivedDomain_PeriodicTrans[iter];
          Buffer_Receive_ReceivedDomain_PeriodicDonor[iter] = Buffer_Send_ReceivedDomain_PeriodicDonor[iter];
        }
        
      }
      
      delete[] Buffer_Send_BoundLine;
      delete[] Buffer_Send_BoundTriangle;
      delete[] Buffer_Send_BoundQuadrilateral;
      delete[] Buffer_Send_Local2Global_Marker;
      
      delete[] Buffer_Send_SendDomain_Periodic;
      delete[] Buffer_Send_SendDomain_PeriodicTrans;
      delete[] Buffer_Send_SendDomain_PeriodicReceptor;
      delete[] Buffer_Send_ReceivedDomain_Periodic;
      delete[] Buffer_Send_ReceivedDomain_PeriodicTrans;
      delete[] Buffer_Send_ReceivedDomain_PeriodicDonor;
      
    }
    
    if ((unsigned long)rank == iDomain) {
      
      if (rank != MASTER_NODE) {
        
        /*--- Receive the buffers with the geometrical information ---*/
        
#ifdef HAVE_MPI
        
        MPI_Probe(MASTER_NODE, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_BoundLine, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 0, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_BoundTriangle, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 1, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 2, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_BoundQuadrilateral, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 2, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 3, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_Local2Global_Marker, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 3, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 4, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_SendDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 4, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 5, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_SendDomain_PeriodicTrans, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 5, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 6, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_SendDomain_PeriodicReceptor, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 6, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 7, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_ReceivedDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 7, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 8, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_ReceivedDomain_PeriodicTrans, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 8, MPI_COMM_WORLD, &status);
        
        MPI_Probe(MASTER_NODE, 9, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
        SU2_MPI::Recv(Buffer_Receive_ReceivedDomain_PeriodicDonor, recv_count, MPI_UNSIGNED_LONG,
                      MASTER_NODE, 9, MPI_COMM_WORLD, &status);
        
#endif
        
      }
      
      /*--- Create the domain structures for the boundaries ---*/
      
      nMarker                = nMarkerDomain;
      nElem_Bound            = new unsigned long[nMarker_Max];
      Local_to_Global_Marker = new unsigned short[nMarker_Max];
      Tag_to_Marker          = new string[nMarker_Max];
      string *TagBound_Copy  = new string[nMarker_Max];
      short *SendRecv_Copy   = new short[nMarker_Max];
      
      for (iMarker = 0; iMarker < nMarker; iMarker++)
        nElem_Bound[iMarker] = nVertexDomain[iMarker];
      
      bound = new CPrimalGrid**[nMarker+(OVERHEAD*size)];
      for (iMarker = 0; iMarker < nMarker+(OVERHEAD*size); iMarker++)
        bound[iMarker] = NULL;
    
      for (iMarker = 0; iMarker < nMarker; iMarker++)
        bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];
      
      /*--- Initialize boundary element counters ---*/
      iBoundLineTotal      = 0;
      iBoundTriangleTotal  = 0;
      iBoundQuadrilateralTotal = 0;
      
      /*--- Store the boundary element connectivity. Note here that we have
       communicated the global index values for the elements, so we need to
       convert this to the local index when instantiating the element. ---*/
      
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        
        iVertexDomain = 0;
        
        for (iBoundLine = 0; iBoundLine < nBoundLine[iMarker]; iBoundLine++) {
          bound[iMarker][iVertexDomain] = new CLine(Global_to_local_Point_recv[Buffer_Receive_BoundLine[iBoundLineTotal*2+0]],
                                                    Global_to_local_Point_recv[Buffer_Receive_BoundLine[iBoundLineTotal*2+1]], 2);
          iVertexDomain++; iBoundLineTotal++;
        }
        for (iBoundTriangle = 0; iBoundTriangle < nBoundTriangle[iMarker]; iBoundTriangle++) {
          bound[iMarker][iVertexDomain] = new CTriangle(Global_to_local_Point_recv[Buffer_Receive_BoundTriangle[iBoundTriangleTotal*3+0]],
                                                        Global_to_local_Point_recv[Buffer_Receive_BoundTriangle[iBoundTriangleTotal*3+1]],
                                                        Global_to_local_Point_recv[Buffer_Receive_BoundTriangle[iBoundTriangleTotal*3+2]], 3);
          iVertexDomain++; iBoundTriangleTotal++;
        }
        for (iBoundQuadrilateral = 0; iBoundQuadrilateral < nBoundQuadrilateral[iMarker]; iBoundQuadrilateral++) {
          bound[iMarker][iVertexDomain] = new CQuadrilateral(Global_to_local_Point_recv[Buffer_Receive_BoundQuadrilateral[iBoundQuadrilateralTotal*4+0]],
                                                             Global_to_local_Point_recv[Buffer_Receive_BoundQuadrilateral[iBoundQuadrilateralTotal*4+1]],
                                                             Global_to_local_Point_recv[Buffer_Receive_BoundQuadrilateral[iBoundQuadrilateralTotal*4+2]],
                                                             Global_to_local_Point_recv[Buffer_Receive_BoundQuadrilateral[iBoundQuadrilateralTotal*4+3]], 3);
          iVertexDomain++; iBoundQuadrilateralTotal++;
        }
        
        Local_to_Global_Marker[iMarker] = Buffer_Receive_Local2Global_Marker[iMarker];
        
        /*--- Now each domain has the right information ---*/
        
        string Grid_Marker = config->GetMarker_All_TagBound(Local_to_Global_Marker[iMarker]);
        short SendRecv = config->GetMarker_All_SendRecv(Local_to_Global_Marker[iMarker]);
        TagBound_Copy[iMarker] = Grid_Marker;
        SendRecv_Copy[iMarker] = SendRecv;
        
      }
      
      /*--- Store total number of each boundary element type ---*/
      
      nelem_edge_bound     = iBoundLineTotal;
      nelem_triangle_bound = iBoundTriangleTotal;
      nelem_quad_bound     = iBoundQuadrilateralTotal;
      
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        config->SetMarker_All_TagBound(iMarker, TagBound_Copy[iMarker]);
        config->SetMarker_All_SendRecv(iMarker, SendRecv_Copy[iMarker]);
      }
      
      /*--- Add the new periodic markers to the domain ---*/
      
      //      iTotalSendDomain_Periodic = 0;
      //      iTotalReceivedDomain_Periodic = 0;
      
      for (jDomain = 0; jDomain < nDomain; jDomain++) {
        
        if (nSendDomain_Periodic[jDomain] != 0) {
          nVertexDomain[nMarker] = 0;
          bound[nMarker] = new CPrimalGrid* [nSendDomain_Periodic[jDomain]];
          
          iVertex = 0;
          for (iTotalSendDomain_Periodic = 0; iTotalSendDomain_Periodic < nTotalSendDomain_Periodic; iTotalSendDomain_Periodic++) {
            if (Buffer_Receive_SendDomain_PeriodicReceptor[iTotalSendDomain_Periodic] == jDomain) {
              bound[nMarker][iVertex] = new CVertexMPI(Global_to_local_Point_recv[Buffer_Receive_SendDomain_Periodic[iTotalSendDomain_Periodic]], nDim);
              bound[nMarker][iVertex]->SetRotation_Type(Buffer_Receive_SendDomain_PeriodicTrans[iTotalSendDomain_Periodic]);
              nVertexDomain[nMarker]++; iVertex++;
            }
          }
          
          Marker_All_SendRecv[nMarker] = jDomain+1;
          nElem_Bound[nMarker] = nVertexDomain[nMarker];
          nMarker++;
        }
        
        if (nReceivedDomain_Periodic[jDomain] != 0) {
          nVertexDomain[nMarker] = 0;
          bound[nMarker] = new CPrimalGrid* [nReceivedDomain_Periodic[jDomain]];
          
          iVertex = 0;
          for (iTotalReceivedDomain_Periodic = 0; iTotalReceivedDomain_Periodic < nTotalReceivedDomain_Periodic; iTotalReceivedDomain_Periodic++) {
            if (Buffer_Receive_ReceivedDomain_PeriodicDonor[iTotalReceivedDomain_Periodic] == jDomain) {
              bound[nMarker][iVertex] = new CVertexMPI(Global_to_local_Point_recv[Buffer_Receive_ReceivedDomain_Periodic[iTotalReceivedDomain_Periodic]], nDim);
              bound[nMarker][iVertex]->SetRotation_Type(Buffer_Receive_ReceivedDomain_PeriodicTrans[iTotalReceivedDomain_Periodic]);
              nVertexDomain[nMarker]++; iVertex++;
            }
          }
          
          Marker_All_SendRecv[nMarker] = -(jDomain+1);
          nElem_Bound[nMarker] = nVertexDomain[nMarker];
          nMarker++;
        }
        
      }
      
      delete[] TagBound_Copy;
      delete[] SendRecv_Copy;
      
      delete[] Buffer_Receive_BoundLine;
      delete[] Buffer_Receive_BoundTriangle;
      delete[] Buffer_Receive_BoundQuadrilateral;
      delete[] Buffer_Receive_Local2Global_Marker;
      
      delete[] Buffer_Receive_SendDomain_Periodic;
      delete[] Buffer_Receive_SendDomain_PeriodicTrans;
      delete[] Buffer_Receive_SendDomain_PeriodicReceptor;
      delete[] Buffer_Receive_ReceivedDomain_Periodic;
      delete[] Buffer_Receive_ReceivedDomain_PeriodicTrans;
      delete[] Buffer_Receive_ReceivedDomain_PeriodicDonor;
      
    }
    
  }
  
  /*--- The MASTER should wait for the sends above to complete ---*/
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- Set the value of Marker_All_SendRecv and Marker_All_TagBound in the config structure ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
  }
  
  /*--- Set the value of Global_nPoint and Global_nPointDomain ---*/
  
  unsigned long Local_nPoint = nPoint;
  unsigned long Local_nPointDomain = nPointDomain;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&Local_nPoint, &Global_nPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nPoint = Local_nPoint;
  Global_nPointDomain = Local_nPointDomain;
#endif
  
  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << Global_nPoint << " vertices including ghost points. " << endl;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nDim_s;
  delete [] nDim_r;
  
  delete [] nPointTotal_s;
  delete [] nPointDomainTotal_s;
  delete [] nPointGhost_s;
  delete [] nPointPeriodic_s;
  delete [] nElemTotal_s;
  delete [] nElemTriangle_s;
  delete [] nElemQuadrilateral_s;
  delete [] nElemTetrahedron_s;
  delete [] nElemHexahedron_s;
  delete [] nElemPrism_s;
  delete [] nElemPyramid_s;
  delete [] nZone_s;
  
  delete [] nPointTotal_r;
  delete [] nPointDomainTotal_r;
  delete [] nPointGhost_r;
  delete [] nPointPeriodic_r;
  delete [] nElemTotal_r;
  delete [] nElemTriangle_r;
  delete [] nElemQuadrilateral_r;
  delete [] nElemTetrahedron_r;
  delete [] nElemHexahedron_r;
  delete [] nElemPrism_r;
  delete [] nElemPyramid_r;
  delete [] nZone_r;
  
  if (rank == MASTER_NODE) {
    delete [] MarkerIn;
    delete [] Buffer_Send_Center;
    delete [] Buffer_Send_Rotation;
    delete [] Buffer_Send_Translate;
    delete [] Buffer_Send_nSendDomain_Periodic;
    delete [] Buffer_Send_nReceivedDomain_Periodic;
    delete [] Marker_All_SendRecv_Copy;
    delete [] Marker_All_TagBound_Copy;
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
      delete [] VertexIn[iMarker];
    delete[] VertexIn;
  }
  
  delete [] Marker_All_TagBound;
  delete [] Buffer_Send_Marker_All_TagBound;
  
  delete [] nSendDomain_Periodic;
  delete [] nReceivedDomain_Periodic;
  delete [] nVertexDomain;
  delete [] nBoundLine;
  delete [] nBoundTriangle;
  delete [] nBoundQuadrilateral;
  delete [] Buffer_Send_nVertexDomain;
  delete [] Buffer_Send_nBoundLine;
  delete [] Buffer_Send_nBoundTriangle;
  delete [] Buffer_Send_nBoundQuadrilateral;
  delete [] Buffer_Send_Marker_All_SendRecv;
  
#ifdef HAVE_MPI
  delete [] send_stat;
  delete [] recv_stat;
  delete [] send_req;
  delete [] recv_req;
#endif
  
  delete [] local_colour_values;
  delete [] ElemIn;
  
}


CPhysicalGeometry::~CPhysicalGeometry(void) {
  
  if (Local_to_Global_Point  != NULL) delete [] Local_to_Global_Point;
  if (Global_to_Local_Marker != NULL) delete [] Global_to_Local_Marker;
  if (Local_to_Global_Marker != NULL) delete [] Local_to_Global_Marker;
  
}



void CPhysicalGeometry::SetSendReceive(CConfig *config) {
  
  unsigned short Counter_Send, Counter_Receive, iMarkerSend, iMarkerReceive;
  unsigned long iVertex, LocalNode;
  unsigned short nMarker_Max = config->GetnMarker_Max();
  unsigned long  iPoint, jPoint, iElem;
  unsigned long *nVertexDomain = new unsigned long[nMarker_Max];
  unsigned short nDomain, iNode, iDomain, jDomain, jNode;
  vector<unsigned long>::iterator it;
  
  vector<vector<unsigned long> > SendTransfLocal;	/*!< \brief Vector to store the type of transformation for this send point. */
  vector<vector<unsigned long> > ReceivedTransfLocal;	/*!< \brief Vector to store the type of transformation for this received point. */
	vector<vector<unsigned long> > SendDomainLocal; /*!< \brief SendDomain[from domain][to domain] and return the point index of the node that must me sended. */
	vector<vector<unsigned long> > ReceivedDomainLocal; /*!< \brief SendDomain[from domain][to domain] and return the point index of the node that must me sended. */
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  /*--- MPI initialization ---*/
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (rank == MASTER_NODE && size > SINGLE_NODE)
    cout << "Establishing MPI communication patterns." << endl;

  nDomain = size;
  
  SendTransfLocal.resize(nDomain);
  ReceivedTransfLocal.resize(nDomain);
  SendDomainLocal.resize(nDomain);
  ReceivedDomainLocal.resize(nDomain);
  
  /*--- Loop over the all the points of the element
   to find the points with different colours, and create the send/received list ---*/
  for (iElem = 0; iElem < nElem; iElem++) {
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      iDomain = node[iPoint]->GetColor();
      
      if (iDomain == rank) {
        for (jNode = 0; jNode < elem[iElem]->GetnNodes(); jNode++) {
          jPoint = elem[iElem]->GetNode(jNode);
          jDomain = node[jPoint]->GetColor();
          
          /*--- If different color and connected by an edge, then we add them to the list ---*/
          if (iDomain != jDomain) {
            
            /*--- We send from iDomain to jDomain the value of iPoint, we save the
             global value becuase we need to sort the lists ---*/
            SendDomainLocal[jDomain].push_back(Local_to_Global_Point[iPoint]);
            /*--- We send from jDomain to iDomain the value of jPoint, we save the
             global value becuase we need to sort the lists ---*/
            ReceivedDomainLocal[jDomain].push_back(Local_to_Global_Point[jPoint]);
            
          }
        }
      }
    }
  }
  
  /*--- Sort the points that must be sended and delete repeated points, note
   that the sorting should be done with the global point (not the local) ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    sort( SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
    it = unique( SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
    SendDomainLocal[iDomain].resize( it - SendDomainLocal[iDomain].begin() );
  }
  
  /*--- Sort the points that must be received and delete repeated points, note
   that the sorting should be done with the global point (not the local) ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    sort( ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
    it = unique( ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
    ReceivedDomainLocal[iDomain].resize( it - ReceivedDomainLocal[iDomain].begin() );
  }
  
  /*--- Create Global to Local Point array, note that the array is smaller (Max_GlobalPoint) than the total
   number of points in the simulation  ---*/
  Max_GlobalPoint = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (Local_to_Global_Point[iPoint] > (long)Max_GlobalPoint)
      Max_GlobalPoint = Local_to_Global_Point[iPoint];
  }
  
  /*--- Set the value of some of the points ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    Global_to_Local_Point[Local_to_Global_Point[iPoint]] = iPoint;
  
  map<long, long>::const_iterator MI;

  /*--- Add the new MPI send receive boundaries, reset the transformation, and save the local value ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (SendDomainLocal[iDomain].size() != 0) {
      nVertexDomain[nMarker] = SendDomainLocal[iDomain].size();
      for (iVertex = 0; iVertex < nVertexDomain[nMarker]; iVertex++) {
        
        MI = Global_to_Local_Point.find(SendDomainLocal[iDomain][iVertex]);
        if (MI != Global_to_Local_Point.end()) iPoint = Global_to_Local_Point[SendDomainLocal[iDomain][iVertex]];
        else iPoint = -1;
          
        SendDomainLocal[iDomain][iVertex] = iPoint;
        SendTransfLocal[iDomain].push_back(0);
      }
      nElem_Bound[nMarker] = nVertexDomain[nMarker];
      bound[nMarker] = new CPrimalGrid*[nElem_Bound[nMarker]];
      nMarker++;
    }
  }
  
  /*--- Add the new MPI receive boundaries, reset the transformation, and save the local value ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (ReceivedDomainLocal[iDomain].size() != 0) {
      nVertexDomain[nMarker] = ReceivedDomainLocal[iDomain].size();
      for (iVertex = 0; iVertex < nVertexDomain[nMarker]; iVertex++) {
        
        MI = Global_to_Local_Point.find(ReceivedDomainLocal[iDomain][iVertex]);
        if (MI != Global_to_Local_Point.end()) iPoint = Global_to_Local_Point[ReceivedDomainLocal[iDomain][iVertex]];
        else iPoint = -1;
        
        ReceivedDomainLocal[iDomain][iVertex] = iPoint;
        ReceivedTransfLocal[iDomain].push_back(0);
      }
      nElem_Bound[nMarker] = nVertexDomain[nMarker];
      bound[nMarker] = new CPrimalGrid*[nElem_Bound[nMarker]];
      nMarker++;
    }
  }
  
  /*--- First compute the Send/Receive boundaries ---*/
  Counter_Send = 0; 	Counter_Receive = 0;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    if (SendDomainLocal[iDomain].size() != 0) Counter_Send++;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    if (ReceivedDomainLocal[iDomain].size() != 0) Counter_Receive++;
  
  iMarkerSend = nMarker - Counter_Send - Counter_Receive;
  iMarkerReceive = nMarker - Counter_Receive;
  
  /*--- First we do the send ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (SendDomainLocal[iDomain].size() != 0) {
      for (iVertex = 0; iVertex < GetnElem_Bound(iMarkerSend); iVertex++) {
        LocalNode = SendDomainLocal[iDomain][iVertex];
        bound[iMarkerSend][iVertex] = new CVertexMPI(LocalNode, nDim);
        bound[iMarkerSend][iVertex]->SetRotation_Type(SendTransfLocal[iDomain][iVertex]);
      }
      Marker_All_SendRecv[iMarkerSend] = iDomain+1;
      iMarkerSend++;
    }
  }
  
  /*--- Second we do the receive ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (ReceivedDomainLocal[iDomain].size() != 0) {
      for (iVertex = 0; iVertex < GetnElem_Bound(iMarkerReceive); iVertex++) {
        LocalNode = ReceivedDomainLocal[iDomain][iVertex];
        bound[iMarkerReceive][iVertex] = new CVertexMPI(LocalNode, nDim);
        bound[iMarkerReceive][iVertex]->SetRotation_Type(ReceivedTransfLocal[iDomain][iVertex]);
      }
      Marker_All_SendRecv[iMarkerReceive] = -(iDomain+1);
      iMarkerReceive++;
    }
  }
 
  /*--- Free memory ---*/
  delete [] nVertexDomain;
}


void CPhysicalGeometry::SetBoundaries(CConfig *config) {
  
  unsigned long iElem_Bound, TotalElem, *nElem_Bound_Copy, iVertex_;
  string Grid_Marker;
  unsigned short iDomain, nDomain, iMarkersDomain, iLoop, *DomainCount, nMarker_Physical, Duplicate_SendReceive, *DomainSendCount, **DomainSendMarkers, *DomainReceiveCount, **DomainReceiveMarkers, nMarker_SendRecv, iMarker, iMarker_;
  CPrimalGrid*** bound_Copy;
  short *Marker_All_SendRecv_Copy;
  bool CheckStart;
  
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  /*--- MPI initialization ---*/
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  nDomain = size+1;
  
  /*--- Count the number of physical markers
   in the boundaries ---*/
  
  nMarker_Physical = 0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      nMarker_Physical++;
    }
  }
  
  /*--- Identify if there are markers that send/received with the same domain,
   they should be together---*/
  
  Duplicate_SendReceive = 0;
  for (iLoop = 0; iLoop < 2; iLoop++) {
    
    DomainCount = new unsigned short [nDomain];
    
    for (iDomain = 0; iDomain < nDomain; iDomain++)
      DomainCount[iDomain] = 0;
    
    if (iLoop == 0) {
      for (iDomain = 0; iDomain < nDomain; iDomain++)
        for (iMarker = 0; iMarker < nMarker; iMarker++)
          if (bound[iMarker][0]->GetVTK_Type() == VERTEX)
            if (Marker_All_SendRecv[iMarker] == iDomain) DomainCount[iDomain]++;
    }
    else {
      for (iDomain = 0; iDomain < nDomain; iDomain++)
        for (iMarker = 0; iMarker < nMarker; iMarker++)
          if (bound[iMarker][0]->GetVTK_Type() == VERTEX)
            if (Marker_All_SendRecv[iMarker] == -iDomain) DomainCount[iDomain]++;
    }
    
    for (iDomain = 0; iDomain < nDomain; iDomain++)
      if (DomainCount[iDomain] > 1) Duplicate_SendReceive++;
    
    delete [] DomainCount;
    
  }
  
  DomainSendCount = new unsigned short [nDomain];
  DomainSendMarkers = new unsigned short *[nDomain];
  DomainReceiveCount = new unsigned short [nDomain];
  DomainReceiveMarkers = new unsigned short *[nDomain];
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    DomainSendCount[iDomain] = 0;
    DomainSendMarkers[iDomain] = new unsigned short [nMarker];
    
    DomainReceiveCount[iDomain] = 0;
    DomainReceiveMarkers[iDomain] = new unsigned short [nMarker];
  }
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
        if (Marker_All_SendRecv[iMarker] == iDomain) {
          DomainSendMarkers[iDomain][DomainSendCount[iDomain]] = iMarker;
          DomainSendCount[iDomain]++;
        }
        if (Marker_All_SendRecv[iMarker] == -iDomain) {
          DomainReceiveMarkers[iDomain][DomainReceiveCount[iDomain]] = iMarker;
          DomainReceiveCount[iDomain]++;
        }
      }
    }
  }
  
  /*--- Create an structure to store the Send/Receive
   boundaries, because they require some reorganization ---*/
  
  nMarker_SendRecv = nMarker - nMarker_Physical - Duplicate_SendReceive;
  bound_Copy = new CPrimalGrid**[nMarker_Physical + nMarker_SendRecv];
  nElem_Bound_Copy = new unsigned long [nMarker_Physical + nMarker_SendRecv];
  Marker_All_SendRecv_Copy = new short [nMarker_Physical + nMarker_SendRecv];
  iMarker_ = nMarker_Physical;
  iVertex_ = 0;
  CheckStart = false;
  
  /*--- Copy and allocate the physical markers in the data structure ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      
      nElem_Bound_Copy[iMarker] = nElem_Bound[iMarker];
      bound_Copy[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == LINE)
          bound_Copy[iMarker][iElem_Bound] = new CLine(bound[iMarker][iElem_Bound]->GetNode(0),
                                                       bound[iMarker][iElem_Bound]->GetNode(1), 2);
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE)
          
          bound_Copy[iMarker][iElem_Bound] = new CTriangle(bound[iMarker][iElem_Bound]->GetNode(0),
                                                           bound[iMarker][iElem_Bound]->GetNode(1),
                                                           bound[iMarker][iElem_Bound]->GetNode(2), 3);
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == QUADRILATERAL)
          bound_Copy[iMarker][iElem_Bound] = new CQuadrilateral(bound[iMarker][iElem_Bound]->GetNode(0),
                                                            bound[iMarker][iElem_Bound]->GetNode(1),
                                                            bound[iMarker][iElem_Bound]->GetNode(2),
                                                            bound[iMarker][iElem_Bound]->GetNode(3), 3);
      }
    }
  }
  
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Compute the total number of elements (adding all the
     boundaries with the same Send/Receive ---*/
    
    if (DomainSendCount[iDomain] != 0) {
      TotalElem = 0;
      for (iMarkersDomain = 0; iMarkersDomain < DomainSendCount[iDomain]; iMarkersDomain++) {
        iMarker = DomainSendMarkers[iDomain][iMarkersDomain];
        TotalElem += nElem_Bound[iMarker];
      }
      if (CheckStart) iMarker_++;
      CheckStart = true;
      iVertex_ = 0;
      nElem_Bound_Copy[iMarker_] = TotalElem;
      bound_Copy[iMarker_] = new CPrimalGrid*[TotalElem];
    }
    
    for (iMarkersDomain = 0; iMarkersDomain < DomainSendCount[iDomain]; iMarkersDomain++) {
      iMarker = DomainSendMarkers[iDomain][iMarkersDomain];
      Marker_All_SendRecv_Copy[iMarker_] = Marker_All_SendRecv[iMarker];
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        bound_Copy[iMarker_][iVertex_] = new CVertexMPI(bound[iMarker][iElem_Bound]->GetNode(0), nDim);
        bound_Copy[iMarker_][iVertex_]->SetRotation_Type(bound[iMarker][iElem_Bound]->GetRotation_Type());
        iVertex_++;
      }
      
    }
    
    /*--- Compute the total number of elements (adding all the
     boundaries with the same Send/Receive ---*/
    
    if (DomainReceiveCount[iDomain] != 0) {
      TotalElem = 0;
      for (iMarkersDomain = 0; iMarkersDomain < DomainReceiveCount[iDomain]; iMarkersDomain++) {
        iMarker = DomainReceiveMarkers[iDomain][iMarkersDomain];
        TotalElem += nElem_Bound[iMarker];
      }
      if (CheckStart) iMarker_++;
      CheckStart = true;
      iVertex_ = 0;
      nElem_Bound_Copy[iMarker_] = TotalElem;
      bound_Copy[iMarker_] = new CPrimalGrid*[TotalElem];
      
    }
    
    for (iMarkersDomain = 0; iMarkersDomain < DomainReceiveCount[iDomain]; iMarkersDomain++) {
      iMarker = DomainReceiveMarkers[iDomain][iMarkersDomain];
      Marker_All_SendRecv_Copy[iMarker_] = Marker_All_SendRecv[iMarker];
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        bound_Copy[iMarker_][iVertex_] = new CVertexMPI(bound[iMarker][iElem_Bound]->GetNode(0), nDim);
        bound_Copy[iMarker_][iVertex_]->SetRotation_Type(bound[iMarker][iElem_Bound]->GetRotation_Type());
        iVertex_++;
      }
      
    }
    
  }
  
  delete [] DomainSendCount;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    delete [] DomainSendMarkers[iDomain];
  delete[] DomainSendMarkers;
  
  delete [] DomainReceiveCount;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    delete [] DomainReceiveMarkers[iDomain];
  delete[] DomainReceiveMarkers;
  
   /*--- Deallocate the bound variables ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
   for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
     if (bound[iMarker][iElem_Bound] != NULL) delete bound[iMarker][iElem_Bound];
    if (bound[iMarker] != NULL) delete [] bound[iMarker];
  }
  if (bound != NULL) delete [] bound;
 
  /*--- Allocate the new bound variables, and set the number of markers ---*/
  
  bound = bound_Copy;
  nMarker = nMarker_Physical + nMarker_SendRecv;
  
  config->SetnMarker_All(nMarker);

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    nElem_Bound[iMarker] = nElem_Bound_Copy[iMarker];
  }
  for (iMarker = nMarker_Physical; iMarker < nMarker; iMarker++) {
    Marker_All_SendRecv[iMarker] = Marker_All_SendRecv_Copy[iMarker];
    config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
    config->SetMarker_All_TagBound(iMarker, "SEND_RECEIVE");
  }
  
  /*--- Update config information storing the boundary information in the right place ---*/
  
  for (iMarker = 0 ; iMarker < nMarker; iMarker++) {
    
    string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    
    if (Marker_Tag != "SEND_RECEIVE") {
      
      /*--- Update config information storing the boundary information in the right place ---*/
      
      Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
      config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
      config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
      config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
      config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
      config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
      config->SetMarker_All_Analyze(iMarker, config->GetMarker_CfgFile_Analyze(Marker_Tag));
      config->SetMarker_All_FSIinterface(iMarker, config->GetMarker_CfgFile_FSIinterface(Marker_Tag));
      config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
      config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
      config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
      config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));

    }
    
    /*--- Send-Receive boundaries definition ---*/
    
    else {
      
      config->SetMarker_All_KindBC(iMarker, SEND_RECEIVE);
      config->SetMarker_All_Monitoring(iMarker, NO);
      config->SetMarker_All_GeoEval(iMarker, NO);
      config->SetMarker_All_Designing(iMarker, NO);
      config->SetMarker_All_Plotting(iMarker, NO);
      config->SetMarker_All_Analyze(iMarker, NO);
  	  config->SetMarker_All_FSIinterface(iMarker, NO);
      config->SetMarker_All_DV(iMarker, NO);
      config->SetMarker_All_Moving(iMarker, NO);
      config->SetMarker_All_PerBound(iMarker, NO);
      config->SetMarker_All_Out_1D(iMarker, NO);
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (config->GetMarker_All_SendRecv(iMarker) < 0)
          node[bound[iMarker][iElem_Bound]->GetNode(0)]->SetDomain(false);
      }
      
    }
    
    /*--- Loop over the surface element to set the boundaries ---*/
    
    unsigned long Point_Surface, iElem_Surface;
    unsigned short iNode_Surface;
    
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
        Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
        node[Point_Surface]->SetBoundary(nMarker);
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[Point_Surface]->SetPhysicalBoundary(true);
        
        if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL &&
            config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX &&
            config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)
          node[Point_Surface]->SetSolidBoundary(true);
      }
    }
    
  }
  
  delete [] Marker_All_SendRecv_Copy;
  delete [] nElem_Bound_Copy;
}

void CPhysicalGeometry::Read_SU2_Format_Parallel(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone) {
  
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short nMarker_Max = config->GetnMarker_Max();
  unsigned long VTK_Type, iMarker, iChar;
  unsigned long iCount = 0;
  unsigned long iElem_Bound = 0, iPoint = 0, ielem = 0;
  unsigned long vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4];
  unsigned long vnodes_tetra[4], vnodes_hexa[8], vnodes_prism[6],
  vnodes_pyramid[5], dummyLong, GlobalIndex;
  unsigned long i;
  long local_index;
  vector<unsigned long>::iterator it;
  char cstr[200];
  su2double Coord_2D[2], Coord_3D[3], AoA_Offset, AoS_Offset, AoA_Current, AoS_Current;
  string::size_type position;
  int rank = MASTER_NODE, size = SINGLE_NODE;
  bool domain_flag = false;
  bool found_transform = false;
  bool harmonic_balance = config->GetUnsteady_Simulation() == HARMONIC_BALANCE;
  bool actuator_disk  = (((config->GetnMarker_ActDiskInlet() != 0) ||
                          (config->GetnMarker_ActDiskOutlet() != 0)) &&
                         ((config->GetKind_SU2() == SU2_CFD) ||
                          ((config->GetKind_SU2() == SU2_DEF) && (config->GetActDisk_SU2_DEF()))));
  if (config->GetActDisk_DoubleSurface()) actuator_disk = false;

  nZone = val_nZone;
  
  /*--- Initialize some additional counters for the parallel partitioning ---*/
  
  unsigned long total_pt_accounted = 0;
  unsigned long rem_points = 0;
  unsigned long element_count = 0;
  unsigned long boundary_marker_count = 0;
  unsigned long node_count = 0;
  unsigned long local_element_count = 0;
  
  /*--- Initialize counters for local/global points & elements ---*/
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  unsigned long LocalIndex, j;
#endif
  
  /*--- Actuator disk preprocesing ---*/
  
  string Marker_Tag_Duplicate;
  bool *ActDisk_Bool 			= NULL, *MapVolumePointBool = NULL, InElem, Perimeter;
  unsigned long *ActDiskPoint_Back = NULL, *VolumePoint_Inv = NULL, *ActDiskPoint_Front_Inv = NULL, ActDiskNewPoints = 0, Counter = 0;
  su2double *CoordXActDisk = NULL, *CoordYActDisk = NULL, *CoordZActDisk = NULL;
  su2double *CoordXVolumePoint = NULL, *CoordYVolumePoint = NULL, *CoordZVolumePoint = NULL;
  su2double Xloc = 0.0, Yloc = 0.0, Zloc = 0.0, Xcg = 0.0;
  unsigned long nElem_Bound_, kPoint;

  vector<unsigned long long> EdgeBegin, EdgeEnd;
  
  unsigned long AuxEdge, iEdge, jEdge, nEdges, nPointVolume, iElem;
  unsigned long long FirstEdgeIndex, SecondEdgeIndex;
  
  vector<unsigned long> ActDiskPoint_Front, VolumePoint, PerimeterPoint;
  
  /*--- If actuator disk, we should split the surface, the first step is to identify the existing boundary ---*/
  
  if (actuator_disk) {
    
    /*--- Open grid file ---*/
    
    strcpy (cstr, val_mesh_filename.c_str());
    mesh_file.open(cstr, ios::in);
    
    /*--- Check the grid ---*/
    
    if (mesh_file.fail()) {
      cout << "There is no mesh file (CPhysicalGeometry)!! " << cstr << endl;
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
    
    /*--- Read grid file with format SU2 ---*/
    
    while (getline (mesh_file, text_line)) {
      
      position = text_line.find ("NDIME=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nDim = atoi(text_line.c_str());
      }
      
      position = text_line.find ("NPOIN=",0);
      if (position != string::npos) {
        text_line.erase (0,6); stringstream test_line(text_line);
        iCount = 0; while (test_line >> dummyLong) iCount++;
        stringstream  stream_line(text_line);
        if (iCount == 2) {  stream_line >> nPoint;  stream_line >> nPointDomain; }
        else if (iCount == 1) { stream_line >> nPoint; }
        for (iPoint = 0; iPoint < nPoint; iPoint++) getline (mesh_file, text_line);
      }
      
      position = text_line.find ("NELEM=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nElem = atoi(text_line.c_str());
        for (iElem = 0; iElem < nElem; iElem++) getline (mesh_file, text_line);
      }
      
      position = text_line.find ("NMARK=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nMarker = atoi(text_line.c_str());
        
        for (iMarker = 0 ; iMarker < nMarker; iMarker++) {
          
          getline (mesh_file, text_line);
          text_line.erase (0,11); string::size_type position;
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find( " ", 0 );
            if (position != string::npos) text_line.erase (position,1); position = text_line.find( "\r", 0 );
            if (position != string::npos) text_line.erase (position,1); position = text_line.find( "\n", 0 );
            if (position != string::npos) text_line.erase (position,1);
          }
          Marker_Tag = text_line.c_str();
          
          getline (mesh_file, text_line);
          text_line.erase (0,13); nElem_Bound_ = atoi(text_line.c_str());
          
          if (( Marker_Tag  == config->GetMarker_ActDiskInlet_TagBound(0)) && (rank == MASTER_NODE))
            cout << "Splitting the surface " << Marker_Tag << "( " << nElem_Bound_  << " boundary elements )." << endl;
          
          if (Marker_Tag  != config->GetMarker_ActDiskInlet_TagBound(0)) {
            for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) { getline (mesh_file, text_line); }
          }
          else {
            
            /*--- Create a list of edges ---*/
            
            for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
              
              getline(mesh_file, text_line);
              
              istringstream bound_line(text_line); bound_line >> VTK_Type;
              
              switch(VTK_Type) {
                case LINE:
                  bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                  EdgeBegin.push_back(vnodes_edge[0]); EdgeEnd.push_back(vnodes_edge[1]);
                  break;
                case TRIANGLE:
                  bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                  EdgeBegin.push_back(vnodes_triangle[0]); EdgeEnd.push_back(vnodes_triangle[1]);
                  EdgeBegin.push_back(vnodes_triangle[1]); EdgeEnd.push_back(vnodes_triangle[2]);
                  EdgeBegin.push_back(vnodes_triangle[2]); EdgeEnd.push_back(vnodes_triangle[0]);
                  break;
                case QUADRILATERAL:
                  bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                  EdgeBegin.push_back(vnodes_quad[0]); EdgeEnd.push_back(vnodes_quad[1]);
                  EdgeBegin.push_back(vnodes_quad[1]); EdgeEnd.push_back(vnodes_quad[2]);
                  EdgeBegin.push_back(vnodes_quad[2]); EdgeEnd.push_back(vnodes_quad[3]);
                  EdgeBegin.push_back(vnodes_quad[3]); EdgeEnd.push_back(vnodes_quad[0]);
                  break;
              }
              
            }
            
            /*--- Set the total number of edges ---*/
            
            nEdges = EdgeBegin.size();
            
            /*--- Sort edges based on local point index, first index is always the largest ---*/
            
            for (iEdge = 0; iEdge <  nEdges; iEdge++) {
              if (EdgeEnd[iEdge] < EdgeBegin[iEdge]) {
                AuxEdge = EdgeEnd[iEdge]; EdgeEnd[iEdge] = EdgeBegin[iEdge]; EdgeBegin[iEdge] = AuxEdge;
              }
            }
            
            /*--- Bubble sort of the points based on the first index   ---*/
            
            for (iEdge = 0; iEdge < nEdges; iEdge++) {
              for (jEdge = iEdge+1; jEdge < nEdges; jEdge++) {
                
                FirstEdgeIndex = EdgeBegin[jEdge] << 31;
                FirstEdgeIndex += EdgeEnd[jEdge];
                
                SecondEdgeIndex = EdgeBegin[iEdge] << 31;
                SecondEdgeIndex += EdgeEnd[iEdge];
                
                if (FirstEdgeIndex <= SecondEdgeIndex) {
                  AuxEdge = EdgeBegin[iEdge]; EdgeBegin[iEdge] = EdgeBegin[jEdge]; EdgeBegin[jEdge] = AuxEdge;
                  AuxEdge = EdgeEnd[iEdge];  EdgeEnd[iEdge] = EdgeEnd[jEdge]; EdgeEnd[jEdge] = AuxEdge;
                }
              }
            }
            
            if (nDim == 3) {
              
              /*--- Check the begning of the list ---*/
              
              if (!((EdgeBegin[0] == EdgeBegin[1]) && (EdgeEnd[0] == EdgeEnd[1]))) {
                PerimeterPoint.push_back(EdgeBegin[0]);
                PerimeterPoint.push_back(EdgeEnd[0]);
              }
              
              for (iEdge = 1; iEdge < nEdges-1; iEdge++) {
                bool Check_1 = !((EdgeBegin[iEdge] == EdgeBegin[iEdge-1]) && (EdgeEnd[iEdge] == EdgeEnd[iEdge-1]));
                bool Check_2 = !((EdgeBegin[iEdge] == EdgeBegin[iEdge+1]) && (EdgeEnd[iEdge] == EdgeEnd[iEdge+1]));
                if ((Check_1 && Check_2)) {
                  PerimeterPoint.push_back(EdgeBegin[iEdge]);
                  PerimeterPoint.push_back(EdgeEnd[iEdge]);
                }
              }
              
              /*--- Check the  end of the list ---*/
              
              if (!((EdgeBegin[nEdges-1] == EdgeBegin[nEdges-2]) && (EdgeEnd[nEdges-1] == EdgeEnd[nEdges-2]))) {
                PerimeterPoint.push_back(EdgeBegin[nEdges-1]);
                PerimeterPoint.push_back(EdgeEnd[nEdges-1]);
              }
              
            }
            
            else {
              
              
              /*--- Create a list with all the points ---*/
              
              for (iEdge = 0; iEdge < nEdges; iEdge++) {
                ActDiskPoint_Front.push_back(EdgeBegin[iEdge]);
                ActDiskPoint_Front.push_back(EdgeEnd[iEdge]);
              }
              
              sort(ActDiskPoint_Front.begin(), ActDiskPoint_Front.end());
              it = unique(ActDiskPoint_Front.begin(), ActDiskPoint_Front.end());
              ActDiskPoint_Front.resize(it - ActDiskPoint_Front.begin());
              
              /*--- Check the begning of the list ---*/
              
              if (!(ActDiskPoint_Front[0] == ActDiskPoint_Front[1]) ) { PerimeterPoint.push_back(ActDiskPoint_Front[0]); }
              
              for (iPoint = 1; iPoint < ActDiskPoint_Front.size()-1; iPoint++) {
                bool Check_1 = !((ActDiskPoint_Front[iPoint] == ActDiskPoint_Front[iPoint-1]) );
                bool Check_2 = !((ActDiskPoint_Front[iPoint] == ActDiskPoint_Front[iPoint+1]) );
                if ((Check_1 && Check_2)) { PerimeterPoint.push_back(ActDiskPoint_Front[iEdge]); }
              }
              
              /*--- Check the  end of the list ---*/
              
              if (!((EdgeBegin[ActDiskPoint_Front.size()-1] == EdgeBegin[ActDiskPoint_Front.size()-2]) )) {
                PerimeterPoint.push_back(ActDiskPoint_Front[ActDiskPoint_Front.size()-1]);
              }
              
              ActDiskPoint_Front.clear();
              
            }
            
            sort(PerimeterPoint.begin(), PerimeterPoint.end());
            it = unique(PerimeterPoint.begin(), PerimeterPoint.end());
            PerimeterPoint.resize(it - PerimeterPoint.begin());
            
            for (iEdge = 0; iEdge < nEdges; iEdge++) {
              
              Perimeter = false;
              for (iPoint = 0; iPoint < PerimeterPoint.size(); iPoint++) {
                if (EdgeBegin[iEdge] == PerimeterPoint[iPoint]) {
                  Perimeter = true; break;
                }
              }
              
              if (!Perimeter) ActDiskPoint_Front.push_back(EdgeBegin[iEdge]);
              
              Perimeter = false;
              for (iPoint = 0; iPoint < PerimeterPoint.size(); iPoint++) {
                if (EdgeEnd[iEdge] == PerimeterPoint[iPoint]) {
                  Perimeter = true; break;
                }
              }
              
              if (!Perimeter) ActDiskPoint_Front.push_back(EdgeEnd[iEdge]);
              
            }
            
            /*--- Sort, and remove repeated points from the disk list of points ---*/
            
            sort(ActDiskPoint_Front.begin(), ActDiskPoint_Front.end());
            it = unique(ActDiskPoint_Front.begin(), ActDiskPoint_Front.end());
            ActDiskPoint_Front.resize(it - ActDiskPoint_Front.begin());
            ActDiskNewPoints = ActDiskPoint_Front.size();
            
            if (rank == MASTER_NODE)
              cout << "Splitting the surface " << Marker_Tag << "( " << ActDiskPoint_Front.size()  << " internal points )." << endl;
            
            /*--- Create a map from original point to the new ones (back plane) ---*/
            
            ActDiskPoint_Back = new unsigned long [nPoint];
            ActDisk_Bool = new bool [nPoint];
            ActDiskPoint_Front_Inv= new unsigned long [nPoint];
            
            for (iPoint = 0; iPoint < nPoint; iPoint++) {
              ActDisk_Bool[iPoint] = false;
              ActDiskPoint_Back[iPoint] = 0;
            }
            
            kPoint = nPoint;
            for (iPoint = 0; iPoint < ActDiskPoint_Front.size(); iPoint++) {
              ActDiskPoint_Front_Inv[ActDiskPoint_Front[iPoint]] = iPoint;
              ActDisk_Bool[ActDiskPoint_Front[iPoint]] = true;
              ActDiskPoint_Back[ActDiskPoint_Front[iPoint]] = kPoint;
             	kPoint++;
            }
            
          }
        }
      }
    }
    
    mesh_file.close();
    
    /*--- Store the coordinates of the new points ---*/
    
    CoordXActDisk = new su2double[ActDiskNewPoints];
    CoordYActDisk = new su2double[ActDiskNewPoints];
    CoordZActDisk = new su2double[ActDiskNewPoints];
    
    strcpy (cstr, val_mesh_filename.c_str());
    mesh_file.open(cstr, ios::in);
    
    
    /*--- Read the coordinates of the points ---*/
    
    while (getline (mesh_file, text_line)) {
      
      position = text_line.find ("NPOIN=",0);
      if (position != string::npos) {
        text_line.erase (0,6); stringstream test_line(text_line);
        iCount = 0; while (test_line >> dummyLong) iCount++;
        stringstream  stream_line(text_line);
        if (iCount == 2) {  stream_line >> nPoint;  stream_line >> nPointDomain; }
        else if (iCount == 1) { stream_line >> nPoint; }
        
        Counter = 0;
        for (iPoint = 0; iPoint < nPoint; iPoint++) {
          getline (mesh_file, text_line);
          istringstream point_line(text_line);
          if (nDim == 2) {point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; }
          else { point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; }
          
          /*--- Compute the CG of the actuator disk surface ---*/
          
          if (ActDisk_Bool[iPoint]) {
            CoordXActDisk[ActDiskPoint_Front_Inv[iPoint]] = Coord_3D[0];
            CoordYActDisk[ActDiskPoint_Front_Inv[iPoint]] = Coord_3D[1];
            Xloc += Coord_3D[0]; Yloc += Coord_3D[1];
            if (nDim == 3) {
              CoordZActDisk[ActDiskPoint_Front_Inv[iPoint]] = Coord_3D[2];
              Zloc += Coord_3D[2];
            }
            Counter++;
          }
          
        }
      }
      
      /*--- Find points that touch the actuator disk surface ---*/
      
      position = text_line.find ("NELEM=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nElem = atoi(text_line.c_str());
        for (iElem = 0; iElem < nElem; iElem++) {
          
          getline(mesh_file, text_line);
          istringstream elem_line(text_line);
          
          elem_line >> VTK_Type;
          
          switch(VTK_Type) {
            case TRIANGLE:
              elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
              InElem = false;
              for (i = 0; i < (unsigned long)N_POINTS_TRIANGLE; i++) {
                if (ActDisk_Bool[vnodes_triangle[i]]) { InElem = true; break; } }
              if (InElem) {
                for (i = 0; i < (unsigned long)N_POINTS_TRIANGLE; i++) {
                  VolumePoint.push_back(vnodes_triangle[i]); } }
              break;
            case QUADRILATERAL:
              elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
              InElem = false;
              for (i = 0; i < (unsigned long)N_POINTS_QUADRILATERAL; i++) {
                if (ActDisk_Bool[vnodes_quad[i]]) { InElem = true; break; } }
              if (InElem) {
                for (i = 0; i < (unsigned long)N_POINTS_QUADRILATERAL; i++) {
                  VolumePoint.push_back(vnodes_quad[i]); } }
              break;
            case TETRAHEDRON:
              elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
              InElem = false;
              for (i = 0; i < (unsigned long)N_POINTS_TETRAHEDRON; i++) {
                if (ActDisk_Bool[vnodes_tetra[i]]) { InElem = true; break; }              }
              if (InElem) {
                for (i = 0; i < (unsigned long)N_POINTS_TETRAHEDRON; i++) {
                  VolumePoint.push_back(vnodes_tetra[i]); } }
              break;
            case HEXAHEDRON:
              elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2];
              elem_line >> vnodes_hexa[3]; elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5];
              elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
              InElem = false;
              for (i = 0; i < (unsigned long)N_POINTS_HEXAHEDRON; i++) {
                if (ActDisk_Bool[vnodes_hexa[i]]) { InElem = true; break; } }
              if (InElem) {
                for (i = 0; i < (unsigned long)N_POINTS_HEXAHEDRON; i++) {
                  VolumePoint.push_back(vnodes_hexa[i]); } }
              break;
            case PRISM:
              elem_line >> vnodes_prism[0]; elem_line >> vnodes_prism[1]; elem_line >> vnodes_prism[2];
              elem_line >> vnodes_prism[3]; elem_line >> vnodes_prism[4]; elem_line >> vnodes_prism[5];
              InElem = false;
              for (i = 0; i < (unsigned long)N_POINTS_PRISM; i++) {
                if (ActDisk_Bool[vnodes_prism[i]]) { InElem = true; break; } }
              if (InElem) {
                for (i = 0; i < (unsigned long)N_POINTS_PRISM; i++) {
                  VolumePoint.push_back(vnodes_prism[i]); } }
              break;
            case PYRAMID:
              elem_line >> vnodes_pyramid[0]; elem_line >> vnodes_pyramid[1]; elem_line >> vnodes_pyramid[2];
              elem_line >> vnodes_pyramid[3]; elem_line >> vnodes_pyramid[4];
              InElem = false;
              for (i = 0; i < (unsigned long)N_POINTS_PYRAMID; i++) {
                if (ActDisk_Bool[vnodes_pyramid[i]]) { InElem = true; break; } }
              if (InElem) {
                for (i = 0; i < (unsigned long)N_POINTS_PYRAMID; i++) {
                  VolumePoint.push_back(vnodes_pyramid[i]); } }
              break;
          }
        }
      }
      
    }
    
    mesh_file.close();
    
    /*--- Compute the CG of the surface ---*/
    
    Xloc /= su2double(Counter);  Yloc /= su2double(Counter);  Zloc /= su2double(Counter);
    
    /*--- Sort,and remove repeated points from the disk list of points ---*/
    
    sort(VolumePoint.begin(), VolumePoint.end());
    it = unique(VolumePoint.begin(), VolumePoint.end());
    VolumePoint.resize(it - VolumePoint.begin());
    nPointVolume = VolumePoint.size();
    
    
    CoordXVolumePoint = new su2double[nPointVolume];
    CoordYVolumePoint = new su2double[nPointVolume];
    CoordZVolumePoint = new su2double[nPointVolume];
    MapVolumePointBool = new bool[nPoint];
    VolumePoint_Inv = new unsigned long[nPoint];
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      MapVolumePointBool[iPoint] = false;
    }
    
    for (iPoint = 0; iPoint < nPointVolume; iPoint++) {
      MapVolumePointBool[VolumePoint[iPoint]] = true;
      VolumePoint_Inv[VolumePoint[iPoint]] = iPoint;
    }
    
    strcpy (cstr, val_mesh_filename.c_str());
    mesh_file.open(cstr, ios::in);
    
    /*--- Store the coordinates of all the surface and volume
     points that touch the actuator disk ---*/
    
    while (getline (mesh_file, text_line)) {
      
      position = text_line.find ("NPOIN=",0);
      if (position != string::npos) {
        text_line.erase (0,6); stringstream test_line(text_line);
        iCount = 0; while (test_line >> dummyLong) iCount++;
        stringstream  stream_line(text_line);
        if (iCount == 2) {  stream_line >> nPoint;  stream_line >> nPointDomain; }
        else if (iCount == 1) { stream_line >> nPoint; }
        
        Counter =0;
        for (iPoint = 0; iPoint < nPoint; iPoint++) {
          getline (mesh_file, text_line);
          istringstream point_line(text_line);
          if (nDim == 2) {point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; }
          else { point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; }
          
          if (MapVolumePointBool[iPoint]) {
            CoordXVolumePoint[VolumePoint_Inv[iPoint]] = Coord_3D[0];
            CoordYVolumePoint[VolumePoint_Inv[iPoint]] = Coord_3D[1];
            if (nDim == 3) { CoordZVolumePoint[VolumePoint_Inv[iPoint]] = Coord_3D[2]; }
          }
        }
      }
      
    }
    
    mesh_file.close();
    
    /*--- Deallocate memory ---*/
    
    delete [] MapVolumePointBool;
    delete [] ActDiskPoint_Front_Inv;
    
    //    }
    
    //  	/*--- Allocate and Send-Receive some of the vectors that we have computed on the MASTER_NODE ---*/
    
    //    SU2_MPI::Bcast(&ActDiskNewPoints, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(&nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(&nPointVolume, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(&Xloc, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(&Yloc, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(&Zloc, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
    //    if (rank != MASTER_NODE) {
    //      MapActDisk 				= new unsigned long [nPoint];
    //      ActDisk_Bool 		= new bool [nPoint];
    //      VolumePoint_Inv 	= new unsigned long [nPoint];
    //      CoordXVolumePoint = new su2double [nPointVolume];
    //      CoordYVolumePoint = new su2double [nPointVolume];
    //      CoordZVolumePoint = new su2double [nPointVolume];
    //      CoordXActDisk 		= new su2double[ActDiskNewPoints];
    //      CoordYActDisk 		= new su2double[ActDiskNewPoints];
    //      CoordZActDisk 		= new su2double[ActDiskNewPoints];
    //    }
    
    //    SU2_MPI::Bcast(MapActDisk, nPoint, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(ActDisk_Bool, nPoint, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(VolumePoint_Inv, nPoint, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(CoordXVolumePoint, nPointVolume, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(CoordYVolumePoint, nPointVolume, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(CoordZVolumePoint, nPointVolume, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(CoordXActDisk, ActDiskNewPoints, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(CoordYActDisk, ActDiskNewPoints, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    //    SU2_MPI::Bcast(CoordZActDisk, ActDiskNewPoints, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    
  }
  
  Global_nPoint  = 0; Global_nPointDomain   = 0; Global_nElem = 0; Global_nElemDomain = 0;
  nelem_edge     = 0; Global_nelem_edge     = 0;
  nelem_triangle = 0; Global_nelem_triangle = 0;
  nelem_quad     = 0; Global_nelem_quad     = 0;
  nelem_tetra    = 0; Global_nelem_tetra    = 0;
  nelem_hexa     = 0; Global_nelem_hexa     = 0;
  nelem_prism    = 0; Global_nelem_prism    = 0;
  nelem_pyramid  = 0; Global_nelem_pyramid  = 0;
  
  /*--- Allocate memory for the linear partition of the mesh. These
   arrays are the size of the number of ranks. ---*/
  
  starting_node = new unsigned long[size];
  ending_node   = new unsigned long[size];
  npoint_procs  = new unsigned long[size];
  
  /*--- Open grid file ---*/
  
  strcpy (cstr, val_mesh_filename.c_str());
  mesh_file.open(cstr, ios::in);
  
  /*--- Check the grid ---*/
  
  if (mesh_file.fail()) {
    cout << "There is no mesh file (CPhysicalGeometry)!! " << cstr << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  /*--- If more than one, find the zone in the mesh file ---*/
  
  if (val_nZone > 1 || harmonic_balance) {
    if (harmonic_balance) {
      if (rank == MASTER_NODE) cout << "Reading time instance " << val_iZone+1 << ":" << endl;
    } else {
      while (getline (mesh_file,text_line)) {
        /*--- Search for the current domain ---*/
        position = text_line.find ("IZONE=",0);
        if (position != string::npos) {
          text_line.erase (0,6);
          unsigned short jDomain = atoi(text_line.c_str());
          if (jDomain == val_iZone+1) {
            if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone+1 << " points:" << endl;
            break;
          }
        }
      }
    }
  }
  
  /*--- Read grid file with format SU2 ---*/
  
  while (getline (mesh_file, text_line)) {
    
    /*--- Read the dimension of the problem ---*/
    
    position = text_line.find ("NDIME=",0);
    if (position != string::npos) {
      if (domain_flag == false) {
        text_line.erase (0,6); nDim = atoi(text_line.c_str());
        if (rank == MASTER_NODE) {
          if (nDim == 2) cout << "Two dimensional problem." << endl;
          if (nDim == 3) cout << "Three dimensional problem." << endl;
        }
        domain_flag = true;
      } else { break; }
    }
    
    /*--- Read if there is any offset in the mesh parameters ---*/
    
    position = text_line.find ("AOA_OFFSET=",0);
    if (position != string::npos) {
      AoA_Offset = 0.0;
      text_line.erase (0,11); AoA_Offset = atof(text_line.c_str());
      
      /*--- The offset is in deg ---*/
      
      AoA_Current = config->GetAoA() + AoA_Offset;
      
      if (config->GetDiscard_InFiles() == false) {
        if ((rank == MASTER_NODE) && (AoA_Offset != 0.0))  {
          cout.precision(6);
          cout << fixed <<"WARNING: AoA in the config file (" << config->GetAoA() << " deg.) +" << endl;
          cout << "         AoA offset in mesh file (" << AoA_Offset << " deg.) = " << AoA_Current << " deg." << endl;
        }
        config->SetAoA_Offset(AoA_Offset);
        config->SetAoA(AoA_Current);
      }
      else {
        if ((rank == MASTER_NODE) && (AoA_Offset != 0.0))
          cout <<"WARNING: Discarding the AoA offset in the geometry file." << endl;
      }
      
    }
    
    position = text_line.find ("AOS_OFFSET=",0);
    if (position != string::npos) {
      AoS_Offset = 0.0;
      text_line.erase (0,11); AoS_Offset = atof(text_line.c_str());
      
      /*--- The offset is in deg ---*/
      
      AoS_Current = config->GetAoS() + AoS_Offset;
      
      if (config->GetDiscard_InFiles() == false) {
        if ((rank == MASTER_NODE) && (AoS_Offset != 0.0))  {
          cout.precision(6);
          cout << fixed <<"WARNING: AoS in the config file (" << config->GetAoS() << " deg.) +" << endl;
          cout << "         AoS offset in mesh file (" << AoS_Offset << " deg.) = " << AoS_Current << " deg." << endl;
        }
        config->SetAoS_Offset(AoS_Offset);
        config->SetAoS(AoS_Current);
      }
      else {
        if ((rank == MASTER_NODE) && (AoS_Offset != 0.0))
          cout <<"WARNING: Discarding the AoS offset in the geometry file." << endl;
      }
      
    }
    
    /*--- Read number of points ---*/
    
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      text_line.erase (0,6);
      
      /*--- Check for ghost points. ---*/
      stringstream test_line(text_line);
      iCount = 0; while (test_line >> dummyLong) iCount++;
      
      /*--- Now read and store the number of points and possible ghost points. ---*/
      
      stringstream  stream_line(text_line);
      if (iCount == 2) {
        
        stream_line >> nPoint;
        stream_line >> nPointDomain;
        
        if (actuator_disk) { nPoint += ActDiskNewPoints;  nPointDomain += ActDiskNewPoints; }
        
        /*--- Set some important point information for parallel simulations. ---*/
        
        Global_nPoint = nPoint;
        Global_nPointDomain = nPointDomain;
        if (rank == MASTER_NODE && size > SINGLE_NODE) {
          cout << Global_nPointDomain << " points and " << Global_nPoint-Global_nPointDomain;
          cout << " ghost points before parallel partitioning." << endl;
        } else if (rank == MASTER_NODE) {
          cout << Global_nPointDomain << " points and " << Global_nPoint-Global_nPointDomain;
          cout << " ghost points." << endl;
        }
        
      } else if (iCount == 1) {
        stream_line >> nPoint;
        
        if (actuator_disk) { nPoint += ActDiskNewPoints; }
        
        nPointDomain = nPoint;
        Global_nPointDomain = nPoint;
        Global_nPoint = nPoint;
        if (rank == MASTER_NODE && size > SINGLE_NODE) {
          cout << nPoint << " points before parallel partitioning." << endl;
        } else if (rank == MASTER_NODE) {
          cout << nPoint << " points." << endl;
        }
      }
      else {
        cout << "NPOIN improperly specified!!" << endl;
#ifndef HAVE_MPI
        exit(EXIT_FAILURE);
#else
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#endif
      }
      
      if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
        cout << "Performing linear partitioning of the grid nodes." << endl;
      
      /*--- Compute the number of points that will be on each processor.
       This is a linear partitioning with the addition of a simple load
       balancing for any remainder points. ---*/
      
      total_pt_accounted = 0;
      for (i = 0; i < (unsigned long)size; i++) {
        npoint_procs[i] = nPoint/size;
        total_pt_accounted = total_pt_accounted + npoint_procs[i];
      }
      
      /*--- Get the number of remainder points after the even division ---*/
      
      rem_points = nPoint-total_pt_accounted;
      for (i = 0; i<rem_points; i++) {
        npoint_procs[i]++;
      }
      
      /*--- Store the local number of nodes and the beginning/end index.
       nPoint is always used to store the local number of points. ---*/
      
      nPoint = npoint_procs[rank];
      starting_node[0] = 0;
      ending_node[0]   = starting_node[0] + npoint_procs[0];
      for (unsigned long i = 1; i < (unsigned long)size; i++) {
        starting_node[i] = ending_node[i-1];
        ending_node[i]   = starting_node[i] + npoint_procs[i] ;
      }
      
      /*--- Here we check if a point in the mesh file lies in the domain
       and if so, then store it on the local processor. We only create enough
       space in the node container for the local nodes at this point. ---*/
      
      nPointNode = nPoint; 
      node = new CPoint*[nPoint];
      iPoint = 0; node_count = 0;
      while (node_count < Global_nPoint) {
        
        if (!actuator_disk) { getline(mesh_file, text_line); }
        else {
          if (node_count < Global_nPoint-ActDiskNewPoints) {
            getline(mesh_file, text_line);
          }
          else {
            ostringstream strsX, strsY, strsZ;
            unsigned long BackActDisk_Index = node_count;
            unsigned long LocalIndex = BackActDisk_Index - (Global_nPoint-ActDiskNewPoints);
            strsX.precision(20); strsY.precision(20); strsZ.precision(20);
            su2double CoordX = CoordXActDisk[LocalIndex]; strsX << scientific << CoordX;
            su2double CoordY = CoordYActDisk[LocalIndex]; strsY << scientific << CoordY;
            su2double CoordZ = CoordZActDisk[LocalIndex]; strsZ << scientific << CoordZ;
            text_line = strsX.str() + "\t" + strsY.str() + "\t" + strsZ.str();
          }
        }

        istringstream point_line(text_line);
        
        /*--- We only read information for this node if it is owned by this
         rank based upon our initial linear partitioning. ---*/
        
        if ((node_count >= starting_node[rank]) && (node_count < ending_node[rank])) {
          switch(nDim) {
            case 2:
              GlobalIndex = node_count;
#ifndef HAVE_MPI
              point_line >> Coord_2D[0]; point_line >> Coord_2D[1];
#else
              if (size > SINGLE_NODE) { point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; point_line >> LocalIndex; point_line >> GlobalIndex; }
              else { point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; LocalIndex = iPoint; GlobalIndex = node_count; }
#endif
              node[iPoint] = new CPoint(Coord_2D[0], Coord_2D[1], GlobalIndex, config);
              iPoint++; break;
            case 3:
              GlobalIndex = node_count;
#ifndef HAVE_MPI
              point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2];
#else
              if (size > SINGLE_NODE) { point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; point_line >> LocalIndex; point_line >> GlobalIndex; }
              else { point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; LocalIndex = iPoint; GlobalIndex = node_count; }
#endif
              node[iPoint] = new CPoint(Coord_3D[0], Coord_3D[1], Coord_3D[2], GlobalIndex, config);
              iPoint++; break;
          }
        }
        node_count++;
      }
    }
  }
  
  mesh_file.close();
  strcpy (cstr, val_mesh_filename.c_str());
  
  /*--- Read the elements in the file with two loops. The first finds
   elements that live in the local partition and builds the adjacency
   for ParMETIS (if parallel). Once we know how many elements we have
   on the local partition, we allocate memory and store those elements. ---*/
  
  map<unsigned long,bool> ElemIn;
  map<unsigned long, bool>::const_iterator MI;
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  /*--- Initialize a vector for the adjacency information (ParMETIS). ---*/
  vector< vector<unsigned long> > adj_nodes(nPoint, vector<unsigned long>(0));
#endif
#endif
  
  /*--- Open the mesh file and find the section with the elements. ---*/
  
  mesh_file.open(cstr, ios::in);
  
  /*--- If more than one, find the zone in the mesh file  ---*/

  if (val_nZone > 1 && !harmonic_balance) {
    while (getline (mesh_file,text_line)) {
      /*--- Search for the current domain ---*/
      position = text_line.find ("IZONE=",0);
      if (position != string::npos) {
        text_line.erase (0,6);
        unsigned short jDomain = atoi(text_line.c_str());
        if (jDomain == val_iZone+1) {
          if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone+1 << " elements:" << endl;
          break;
        }
      }
    }
  }
  
  while (getline (mesh_file, text_line)) {
    
    /*--- Read the information about inner elements ---*/
    
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nElem = atoi(text_line.c_str());
      
      /*--- Store total number of elements in the original mesh ---*/
      
      Global_nElem = nElem;
      if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
        cout << Global_nElem << " interior elements before parallel partitioning." << endl;
      
      /*--- Loop over all the volumetric elements and store any element that
       contains at least one of an owned node for this rank (i.e., there will
       be element redundancy, since multiple ranks will store the same elems
       on the boundaries of the initial linear partitioning. ---*/
      
      element_count = 0; local_element_count = 0;
      while (element_count < Global_nElem) {
        getline(mesh_file, text_line);
        istringstream elem_line(text_line);
        
        /*--- Decide whether this rank needs each element. If so, build the
         adjacency arrays needed by ParMETIS and store the element connectivity.
         Note that every proc starts it's node indexing from zero. ---*/
        
        elem_line >> VTK_Type;
        switch(VTK_Type) {
            
          case TRIANGLE:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> vnodes_triangle[0];
            elem_line >> vnodes_triangle[1];
            elem_line >> vnodes_triangle[2];
            
            if (actuator_disk) {
              for (unsigned short i = 0; i<N_POINTS_TRIANGLE; i++) {
                if (ActDisk_Bool[vnodes_triangle[i]]) {
                  
                  Xcg = 0.0; Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_TRIANGLE; j++) {
                    if (vnodes_triangle[j] < Global_nPoint-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_triangle[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) {
                      vnodes_triangle[i] = ActDiskPoint_Back[vnodes_triangle[i]];
                    }
                    else { vnodes_triangle[i] = vnodes_triangle[i]; }
                  }
                  
                }
              }
            }
            
            /*--- Decide whether we need to store this element, i.e., check if
             any of the nodes making up this element have a global index value
             that falls within the range of our linear partitioning. ---*/
            
            for (unsigned short i = 0; i < N_POINTS_TRIANGLE; i++) {
              
              local_index = vnodes_triangle[i]-starting_node[rank];
              
              if ((local_index >= 0) && (local_index < (long)nPoint)) {
                
                /*--- This node is within our linear partition. Mark this
                 entire element to be added to our list for this rank, and
                 add the neighboring nodes to this nodes' adjacency list. ---*/
                
                ElemIn[element_count] = true;
                
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
                /*--- Build adjacency assuming the VTK connectivity ---*/
                for (unsigned short j=0; j<N_POINTS_TRIANGLE; j++) {
                  if (i != j) adj_nodes[local_index].push_back(vnodes_triangle[j]);
                }
#endif
#endif
              }
            }
            
            MI = ElemIn.find(element_count);
            if (MI != ElemIn.end()) local_element_count++;
            
            break;
            
          case QUADRILATERAL:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> vnodes_quad[0];
            elem_line >> vnodes_quad[1];
            elem_line >> vnodes_quad[2];
            elem_line >> vnodes_quad[3];
            
            if (actuator_disk) {
              for (unsigned short  i = 0; i<N_POINTS_QUADRILATERAL; i++) {
                if (ActDisk_Bool[vnodes_quad[i]]) {
                  
                  Xcg = 0.0; Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_QUADRILATERAL; j++) {
                    if (vnodes_quad[j] < Global_nPoint-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_quad[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) {
                      vnodes_quad[i] = ActDiskPoint_Back[vnodes_quad[i]];
                    }
                    else { vnodes_quad[i] = vnodes_quad[i]; }
                  }
                  
                }
              }
            }

            /*--- Decide whether we need to store this element, i.e., check if
             any of the nodes making up this element have a global index value
             that falls within the range of our linear partitioning. ---*/
            
            for (unsigned short i = 0; i < N_POINTS_QUADRILATERAL; i++) {
              
              local_index = vnodes_quad[i]-starting_node[rank];
              
              if ((local_index >= 0) && (local_index < (long)nPoint)) {
                
                /*--- This node is within our linear partition. Mark this
                 entire element to be added to our list for this rank, and
                 add the neighboring nodes to this nodes' adjacency list. ---*/
                
                ElemIn[element_count] = true;
                
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
                /*--- Build adjacency assuming the VTK connectivity ---*/
                adj_nodes[local_index].push_back(vnodes_quad[(i+1)%4]);
                adj_nodes[local_index].push_back(vnodes_quad[(i+3)%4]);
#endif
#endif
              }
            }
            
            MI = ElemIn.find(element_count);
            if (MI != ElemIn.end()) local_element_count++;
            
            break;
            
          case TETRAHEDRON:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> vnodes_tetra[0];
            elem_line >> vnodes_tetra[1];
            elem_line >> vnodes_tetra[2];
            elem_line >> vnodes_tetra[3];
            
            if (actuator_disk) {
              for (unsigned short  i = 0; i<N_POINTS_TETRAHEDRON; i++) {
                if (ActDisk_Bool[vnodes_tetra[i]]) {
                  
                  Xcg = 0.0; Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_TETRAHEDRON; j++) {
                    if (vnodes_tetra[j] < Global_nPoint-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_tetra[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) {
                      vnodes_tetra[i] = ActDiskPoint_Back[vnodes_tetra[i]];
                    }
                    else { vnodes_tetra[i] = vnodes_tetra[i]; }
                  }
                  
                }
              }
            }
            
            /*--- Decide whether we need to store this element, i.e., check if
             any of the nodes making up this element have a global index value
             that falls within the range of our linear partitioning. ---*/
            
            for (unsigned short i = 0; i < N_POINTS_TETRAHEDRON; i++) {
              
              local_index = vnodes_tetra[i]-starting_node[rank];
              
              if ((local_index >= 0) && (local_index < (long)nPoint)) {
                
                /*--- This node is within our linear partition. Mark this
                 entire element to be added to our list for this rank, and
                 add the neighboring nodes to this nodes' adjacency list. ---*/
                
                ElemIn[element_count] = true;
                
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
                /*--- Build adjacency assuming the VTK connectivity ---*/
                for (unsigned short j=0; j<N_POINTS_TETRAHEDRON; j++) {
                  if (i != j) adj_nodes[local_index].push_back(vnodes_tetra[j]);
                }
#endif
#endif
              }
            }
            
            MI = ElemIn.find(element_count);
            if (MI != ElemIn.end()) local_element_count++;
            
            break;
            
          case HEXAHEDRON:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> vnodes_hexa[0];
            elem_line >> vnodes_hexa[1];
            elem_line >> vnodes_hexa[2];
            elem_line >> vnodes_hexa[3];
            elem_line >> vnodes_hexa[4];
            elem_line >> vnodes_hexa[5];
            elem_line >> vnodes_hexa[6];
            elem_line >> vnodes_hexa[7];
            
            if (actuator_disk) {
              for (unsigned short  i = 0; i<N_POINTS_HEXAHEDRON; i++) {
                if (ActDisk_Bool[vnodes_hexa[i]]) {
                  
                  Xcg = 0.0; Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_HEXAHEDRON; j++) {
                    if (vnodes_hexa[j] < Global_nPoint-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_hexa[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) { vnodes_hexa[i] = ActDiskPoint_Back[vnodes_hexa[i]]; }
                    else { vnodes_hexa[i] = vnodes_hexa[i]; }
                  }
                }
              }
            }

            /*--- Decide whether we need to store this element, i.e., check if
             any of the nodes making up this element have a global index value
             that falls within the range of our linear partitioning. ---*/
            
            for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
              
              local_index = vnodes_hexa[i]-starting_node[rank];
              
              if ((local_index >= 0) && (local_index < (long)nPoint)) {
                
                /*--- This node is within our linear partition. Mark this
                 entire element to be added to our list for this rank, and
                 add the neighboring nodes to this nodes' adjacency list. ---*/
                
                ElemIn[element_count] = true;
                
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
                /*--- Build adjacency assuming the VTK connectivity ---*/
                if (i < 4) {
                  adj_nodes[local_index].push_back(vnodes_hexa[(i+1)%4]);
                  adj_nodes[local_index].push_back(vnodes_hexa[(i+3)%4]);
                } else {
                  adj_nodes[local_index].push_back(vnodes_hexa[(i-3)%4+4]);
                  adj_nodes[local_index].push_back(vnodes_hexa[(i-1)%4+4]);
                }
                adj_nodes[local_index].push_back(vnodes_hexa[(i+4)%8]);
#endif
#endif
              }
            }
            
            MI = ElemIn.find(element_count);
            if (MI != ElemIn.end()) local_element_count++;
            
            break;
            
          case PRISM:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> vnodes_prism[0];
            elem_line >> vnodes_prism[1];
            elem_line >> vnodes_prism[2];
            elem_line >> vnodes_prism[3];
            elem_line >> vnodes_prism[4];
            elem_line >> vnodes_prism[5];
            
            if (actuator_disk) {
              for (unsigned short i = 0; i<N_POINTS_PRISM; i++) {
                if (ActDisk_Bool[vnodes_prism[i]]) {
                  
                  Xcg = 0.0; Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_PRISM; j++) {
                    if (vnodes_prism[j] < Global_nPoint-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_prism[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) { vnodes_prism[i] = ActDiskPoint_Back[vnodes_prism[i]]; }
                    else { vnodes_prism[i] = vnodes_prism[i]; }
                  }
                }
              }
            }

            /*--- Decide whether we need to store this element, i.e., check if
             any of the nodes making up this element have a global index value
             that falls within the range of our linear partitioning. ---*/
            
            for (unsigned short i = 0; i < N_POINTS_PRISM; i++) {
              
              local_index = vnodes_prism[i]-starting_node[rank];
              
              if ((local_index >= 0) && (local_index < (long)nPoint)) {
                
                /*--- This node is within our linear partition. Mark this
                 entire element to be added to our list for this rank, and
                 add the neighboring nodes to this nodes' adjacency list. ---*/
                
                ElemIn[element_count] = true;
                
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
                /*--- Build adjacency assuming the VTK connectivity ---*/
                if (i < 3) {
                  adj_nodes[local_index].push_back(vnodes_prism[(i+1)%3]);
                  adj_nodes[local_index].push_back(vnodes_prism[(i+2)%3]);
                } else {
                  adj_nodes[local_index].push_back(vnodes_prism[(i-2)%3+3]);
                  adj_nodes[local_index].push_back(vnodes_prism[(i-1)%3+3]);
                }
                adj_nodes[local_index].push_back(vnodes_prism[(i+3)%6]);
#endif
#endif
              }
            }
            
            MI = ElemIn.find(element_count);
            if (MI != ElemIn.end()) local_element_count++;
            
            break;
            
          case PYRAMID:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> vnodes_pyramid[0];
            elem_line >> vnodes_pyramid[1];
            elem_line >> vnodes_pyramid[2];
            elem_line >> vnodes_pyramid[3];
            elem_line >> vnodes_pyramid[4];
            
            if (actuator_disk) {
              for (unsigned short i = 0; i<N_POINTS_PYRAMID; i++) {
                if (ActDisk_Bool[vnodes_pyramid[i]]) {
                  
                  Xcg = 0.0; Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_PYRAMID; j++) {
                    if (vnodes_pyramid[j] < Global_nPoint-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_pyramid[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) { vnodes_pyramid[i] = ActDiskPoint_Back[vnodes_pyramid[i]]; }
                    else { vnodes_pyramid[i] = vnodes_pyramid[i]; }
                  }
                }
              }
            }

            /*--- Decide whether we need to store this element, i.e., check if
             any of the nodes making up this element have a global index value
             that falls within the range of our linear partitioning. ---*/
            
            for (unsigned short i = 0; i < N_POINTS_PYRAMID; i++) {
              
              local_index = vnodes_pyramid[i]-starting_node[rank];
              
              if ((local_index >= 0) && (local_index < (long)nPoint)) {
                
                /*--- This node is within our linear partition. Mark this
                 entire element to be added to our list for this rank, and
                 add the neighboring nodes to this nodes' adjacency list. ---*/
                
                ElemIn[element_count] = true;
                
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
                /*--- Build adjacency assuming the VTK connectivity ---*/
                if (i < 4) {
                  adj_nodes[local_index].push_back(vnodes_pyramid[(i+1)%4]);
                  adj_nodes[local_index].push_back(vnodes_pyramid[(i+3)%4]);
                  adj_nodes[local_index].push_back(vnodes_pyramid[4]);
                } else {
                  adj_nodes[local_index].push_back(vnodes_pyramid[0]);
                  adj_nodes[local_index].push_back(vnodes_pyramid[1]);
                  adj_nodes[local_index].push_back(vnodes_pyramid[2]);
                  adj_nodes[local_index].push_back(vnodes_pyramid[3]);
                }
#endif
#endif
              }
            }
            
            MI = ElemIn.find(element_count);
            if (MI != ElemIn.end()) local_element_count++;
            
            break;
        }
        element_count++;
      }
      if (element_count == Global_nElem) break;
    }
  }
  
  mesh_file.close();
  
  /*--- Store the number of elements on the whole domain, excluding halos. ---*/

  Global_nElemDomain = element_count;

  /*--- Store the number of local elements on each rank after determining
   which elements must be kept in the loop above. ---*/
  
  nElem = local_element_count;
  
  /*--- Begin dealing with the partitioning by adjusting the adjacency
   information and clear out memory where possible. ---*/
  
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  
  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << "Calling the partitioning functions." << endl;
  
  /*--- Post process the adjacency information in order to get it into the
   proper format before sending the data to ParMETIS. We need to remove
   repeats and adjust the size of the array for each local node. ---*/
  
  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << "Building the graph adjacency structure." << endl;
  
  unsigned long loc_adjc_size=0;
  vector<unsigned long> adjac_vec;
  unsigned long adj_elem_size;
  
  xadj = new idx_t [npoint_procs[rank]+1];
  xadj[0]=0;
  vector<unsigned long> temp_adjacency;
  unsigned long local_count=0;
  
  /*--- Here, we transfer the adjacency information from a multi-dim vector
   on a node-by-node basis into a single vector container. First, we sort
   the entries and remove the duplicates we find for each node, then we
   copy it into the single vect and clear memory from the multi-dim vec. ---*/
  
  for (unsigned long i = 0; i < nPoint; i++) {
    
    for (j = 0; j<adj_nodes[i].size(); j++) {
      temp_adjacency.push_back(adj_nodes[i][j]);
    }
    
    sort(temp_adjacency.begin(), temp_adjacency.end());
    it = unique(temp_adjacency.begin(), temp_adjacency.end());
    loc_adjc_size = it - temp_adjacency.begin();
    
    temp_adjacency.resize(loc_adjc_size);
    xadj[local_count+1]=xadj[local_count]+loc_adjc_size;
    local_count++;
    
    for (j = 0; j<loc_adjc_size; j++) {
      adjac_vec.push_back(temp_adjacency[j]);
    }
    
    temp_adjacency.clear();
    adj_nodes[i].clear();
    
  }
  
  /*--- Now that we know the size, create the final adjacency array. This
   is the array that we will feed to ParMETIS for partitioning. ---*/
  
  adj_elem_size = xadj[npoint_procs[rank]];
  adjacency = new idx_t [adj_elem_size];
  copy(adjac_vec.begin(), adjac_vec.end(), adjacency);
  
  xadj_size = npoint_procs[rank]+1;
  adjacency_size = adj_elem_size;
  
  /*--- Free temporary memory used to build the adjacency. ---*/
  
  adjac_vec.clear();
  adj_nodes.clear();
  
#endif
#endif
  
  /*--- Open the mesh file again and now that we know the number of
   elements needed on each partition, allocate memory for them. ---*/
  
  mesh_file.open(cstr, ios::in);
  
  /*--- If more than one, find the zone in the mesh file  ---*/
  
  if (val_nZone > 1 && !harmonic_balance) {
    while (getline (mesh_file,text_line)) {
      /*--- Search for the current domain ---*/
      position = text_line.find ("IZONE=",0);
      if (position != string::npos) {
        text_line.erase (0,6);
        unsigned short jDomain = atoi(text_line.c_str());
        if (jDomain == val_iZone+1) {
          if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone+1 << " elements:" << endl;
          break;
        }
      }
    }
  }
  
  while (getline (mesh_file, text_line)) {
    
    /*--- Read the information about inner elements ---*/
    
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      
      /*--- Allocate space for elements ---*/
      elem = new CPrimalGrid*[nElem];
      
      /*--- Set up the global to local element mapping. ---*/
      Global_to_Local_Elem.clear();
      
      if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
        cout << "Distributing elements across all ranks." << endl;
      
      /*--- Loop over all the volumetric elements and store any element that
       contains at least one of an owned node for this rank (i.e., there will
       be element redundancy, since multiple ranks will store the same elems
       on the boundaries of the initial linear partitioning. ---*/
      
      element_count = 0; local_element_count = 0;
      while (element_count < Global_nElem) {
        getline(mesh_file, text_line);
        istringstream elem_line(text_line);
        
        /*--- If this element was marked as required, check type and store. ---*/
        
        map<unsigned long, bool>::const_iterator MI = ElemIn.find(element_count);
        if (MI != ElemIn.end()) {
          
          elem_line >> VTK_Type;
          switch(VTK_Type) {
              
            case TRIANGLE:
              
              /*--- Load the connectivity for this element. ---*/
              
              elem_line >> vnodes_triangle[0];
              elem_line >> vnodes_triangle[1];
              elem_line >> vnodes_triangle[2];
              
              if (actuator_disk) {
                for (unsigned short i = 0; i<N_POINTS_TRIANGLE; i++) {
                  if (ActDisk_Bool[vnodes_triangle[i]]) {
                    
                    Xcg = 0.0; Counter = 0;
                    for (unsigned short j = 0; j<N_POINTS_TRIANGLE; j++) {
                      if (vnodes_triangle[j] < Global_nPoint-ActDiskNewPoints) {
                        Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_triangle[j]]];
                        Counter++;
                      }
                    }
                    Xcg = Xcg / su2double(Counter);
                    
                    if (Counter != 0)  {
                      if (Xcg > Xloc) {
                        vnodes_triangle[i] = ActDiskPoint_Back[vnodes_triangle[i]];
                      }
                      else { vnodes_triangle[i] = vnodes_triangle[i]; }
                    }
                    
                  }
                }
              }

              /*--- If any of the nodes were within the linear partition, the
               element is added to our element data structure. ---*/
              
              Global_to_Local_Elem[element_count] = local_element_count;
              elem[local_element_count] = new CTriangle(vnodes_triangle[0],
                                                        vnodes_triangle[1],
                                                        vnodes_triangle[2], 2);
              local_element_count++;
              nelem_triangle++;
              break;
              
            case QUADRILATERAL:
              
              /*--- Load the connectivity for this element. ---*/
              
              elem_line >> vnodes_quad[0];
              elem_line >> vnodes_quad[1];
              elem_line >> vnodes_quad[2];
              elem_line >> vnodes_quad[3];
              
              if (actuator_disk) {
                for (unsigned short i = 0; i<N_POINTS_QUADRILATERAL; i++) {
                  if (ActDisk_Bool[vnodes_quad[i]]) {
                    
                    Xcg = 0.0; Counter = 0;
                    for (unsigned short j = 0; j<N_POINTS_QUADRILATERAL; j++) {
                      if (vnodes_quad[j] < Global_nPoint-ActDiskNewPoints) {
                        Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_quad[j]]];
                        Counter++;
                      }
                    }
                    Xcg = Xcg / su2double(Counter);
                    
                    
                    if (Counter != 0)  {
                      if (Xcg > Xloc) {
                        vnodes_quad[i] = ActDiskPoint_Back[vnodes_quad[i]];
                      }
                      else { vnodes_quad[i] = vnodes_quad[i]; }
                    }
                    
                  }
                }
              }
              
              /*--- If any of the nodes were within the linear partition, the
               element is added to our element data structure. ---*/
              
              Global_to_Local_Elem[element_count] = local_element_count;
              elem[local_element_count] = new CQuadrilateral(vnodes_quad[0],
                                                             vnodes_quad[1],
                                                             vnodes_quad[2],
                                                             vnodes_quad[3], 2);
              local_element_count++;
              nelem_quad++;
              break;
              
            case TETRAHEDRON:
              
              /*--- Load the connectivity for this element. ---*/
              
              elem_line >> vnodes_tetra[0];
              elem_line >> vnodes_tetra[1];
              elem_line >> vnodes_tetra[2];
              elem_line >> vnodes_tetra[3];
              
              if (actuator_disk) {
                for (unsigned short i = 0; i<N_POINTS_TETRAHEDRON; i++) {
                  if (ActDisk_Bool[vnodes_tetra[i]]) {
                    
                    Xcg = 0.0; Counter = 0;
                    for (unsigned short j = 0; j<N_POINTS_TETRAHEDRON; j++) {
                      if (vnodes_tetra[j] < Global_nPoint-ActDiskNewPoints) {
                        Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_tetra[j]]];
                        Counter++;
                      }
                    }
                    Xcg = Xcg / su2double(Counter);
                    
                    
                    if (Counter != 0)  {
                      if (Xcg > Xloc) {
                        vnodes_tetra[i] = ActDiskPoint_Back[vnodes_tetra[i]];
                      }
                      else { vnodes_tetra[i] = vnodes_tetra[i]; }
                    }
                    
                  }
                }
              }
              
              /*--- If any of the nodes were within the linear partition, the
               element is added to our element data structure. ---*/
              
              Global_to_Local_Elem[element_count] = local_element_count;
              elem[local_element_count] = new CTetrahedron(vnodes_tetra[0],
                                                           vnodes_tetra[1],
                                                           vnodes_tetra[2],
                                                           vnodes_tetra[3]);
              local_element_count++;
              nelem_tetra++;
              break;
              
            case HEXAHEDRON:
              
              /*--- Load the connectivity for this element. ---*/
              
              elem_line >> vnodes_hexa[0];
              elem_line >> vnodes_hexa[1];
              elem_line >> vnodes_hexa[2];
              elem_line >> vnodes_hexa[3];
              elem_line >> vnodes_hexa[4];
              elem_line >> vnodes_hexa[5];
              elem_line >> vnodes_hexa[6];
              elem_line >> vnodes_hexa[7];
              
              if (actuator_disk) {
                for (unsigned short i = 0; i<N_POINTS_HEXAHEDRON; i++) {
                  if (ActDisk_Bool[vnodes_hexa[i]]) {
                    
                    Xcg = 0.0; Counter = 0;
                    for (unsigned short j = 0; j<N_POINTS_HEXAHEDRON; j++) {
                      if (vnodes_hexa[j] < Global_nPoint-ActDiskNewPoints) {
                        Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_hexa[j]]];
                        Counter++;
                      }
                    }
                    Xcg = Xcg / su2double(Counter);
                    
                    if (Counter != 0)  {
                      if (Xcg > Xloc) { vnodes_hexa[i] = ActDiskPoint_Back[vnodes_hexa[i]]; }
                      else { vnodes_hexa[i] = vnodes_hexa[i]; }
                    }
                  }
                }
              }

              /*--- If any of the nodes were within the linear partition, the
               element is added to our element data structure. ---*/
              
              Global_to_Local_Elem[element_count] = local_element_count;
              elem[local_element_count] = new CHexahedron(vnodes_hexa[0],
                                                          vnodes_hexa[1],
                                                          vnodes_hexa[2],
                                                          vnodes_hexa[3],
                                                          vnodes_hexa[4],
                                                          vnodes_hexa[5],
                                                          vnodes_hexa[6],
                                                          vnodes_hexa[7]);
              local_element_count++;
              nelem_hexa++;
              break;
              
            case PRISM:
              
              /*--- Load the connectivity for this element. ---*/
              
              elem_line >> vnodes_prism[0];
              elem_line >> vnodes_prism[1];
              elem_line >> vnodes_prism[2];
              elem_line >> vnodes_prism[3];
              elem_line >> vnodes_prism[4];
              elem_line >> vnodes_prism[5];
              
              if (actuator_disk) {
                for (unsigned short i = 0; i<N_POINTS_PRISM; i++) {
                  if (ActDisk_Bool[vnodes_prism[i]]) {
                    
                    Xcg = 0.0; Counter = 0;
                    for (unsigned short j = 0; j<N_POINTS_PRISM; j++) {
                      if (vnodes_prism[j] < Global_nPoint-ActDiskNewPoints) {
                        Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_prism[j]]];
                        Counter++;
                      }
                    }
                    Xcg = Xcg / su2double(Counter);
                    
                    if (Counter != 0)  {
                      if (Xcg > Xloc) { vnodes_prism[i] = ActDiskPoint_Back[vnodes_prism[i]]; }
                      else { vnodes_prism[i] = vnodes_prism[i]; }
                    }
                  }
                }
              }

              /*--- If any of the nodes were within the linear partition, the
               element is added to our element data structure. ---*/
              
              Global_to_Local_Elem[element_count] = local_element_count;
              elem[local_element_count] = new CPrism(vnodes_prism[0],
                                                     vnodes_prism[1],
                                                     vnodes_prism[2],
                                                     vnodes_prism[3],
                                                     vnodes_prism[4],
                                                     vnodes_prism[5]);
              local_element_count++;
              nelem_prism++;
              break;
              
            case PYRAMID:
              
              /*--- Load the connectivity for this element. ---*/
              
              elem_line >> vnodes_pyramid[0];
              elem_line >> vnodes_pyramid[1];
              elem_line >> vnodes_pyramid[2];
              elem_line >> vnodes_pyramid[3];
              elem_line >> vnodes_pyramid[4];
              
              if (actuator_disk) {
                for (unsigned short i = 0; i<N_POINTS_PYRAMID; i++) {
                  if (ActDisk_Bool[vnodes_pyramid[i]]) {
                    
                    Xcg = 0.0; Counter = 0;
                    for (unsigned short j = 0; j<N_POINTS_PYRAMID; j++) {
                      if (vnodes_pyramid[j] < Global_nPoint-ActDiskNewPoints) {
                        Xcg += CoordXVolumePoint[VolumePoint_Inv[vnodes_pyramid[j]]];
                        Counter++;
                      }
                    }
                    Xcg = Xcg / su2double(Counter);
                    
                    if (Counter != 0)  {
                      if (Xcg > Xloc) { vnodes_pyramid[i] = ActDiskPoint_Back[vnodes_pyramid[i]]; }
                      else { vnodes_pyramid[i] = vnodes_pyramid[i]; }
                    }
                  }
                }
              }

              /*--- If any of the nodes were within the linear partition, the
               element is added to our element data structure. ---*/
              
              Global_to_Local_Elem[element_count]=local_element_count;
              elem[local_element_count] = new CPyramid(vnodes_pyramid[0],
                                                       vnodes_pyramid[1],
                                                       vnodes_pyramid[2],
                                                       vnodes_pyramid[3],
                                                       vnodes_pyramid[4]);
              local_element_count++;
              nelem_pyramid++;
              break;
              
          }
        }
        element_count++;
      }
      if (element_count == Global_nElem) break;
    }
  }
  
  mesh_file.close();
  
  /*--- For now, the boundary marker information is still read by the
   master node alone (and eventually distributed by the master as well).
   In the future, this component will also be performed in parallel. ---*/
  
  mesh_file.open(cstr, ios::in);
  
  /*--- If more than one, find the zone in the mesh file ---*/
  

  if (val_nZone > 1 && !harmonic_balance) {
    while (getline (mesh_file,text_line)) {
      /*--- Search for the current domain ---*/
      position = text_line.find ("IZONE=",0);
      if (position != string::npos) {
        text_line.erase (0,6);
        unsigned short jDomain = atoi(text_line.c_str());
        if (jDomain == val_iZone+1) {
          if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone+1 << " markers:" << endl;
          break;
        }
      }
    }
  }
    
    while (getline (mesh_file, text_line)) {
      
      /*--- Read number of markers ---*/
      
      position = text_line.find ("NMARK=",0);
      boundary_marker_count = 0;
      
      if (position != string::npos) {
        text_line.erase (0,6); nMarker = atoi(text_line.c_str());
        
        if (actuator_disk) { nMarker++;  }
        
        if (rank == MASTER_NODE) cout << nMarker << " surface markers." << endl;
        config->SetnMarker_All(nMarker);
        bound = new CPrimalGrid**[nMarker];
        nElem_Bound = new unsigned long [nMarker];
        Tag_to_Marker = new string [nMarker_Max];
        
        bool duplicate = false;
        iMarker=0;
        do {
          
          getline (mesh_file, text_line);
          text_line.erase (0,11);
          string::size_type position;
          
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find( " ", 0 );
            if (position != string::npos) text_line.erase (position,1);
            position = text_line.find( "\r", 0 );
            if (position != string::npos) text_line.erase (position,1);
            position = text_line.find( "\n", 0 );
            if (position != string::npos) text_line.erase (position,1);
          }
          Marker_Tag = text_line.c_str();
          
          duplicate = false;
          if ((actuator_disk) && ( Marker_Tag  == config->GetMarker_ActDiskInlet_TagBound(0))) {
            duplicate = true;
            Marker_Tag_Duplicate  = config->GetMarker_ActDiskOutlet_TagBound(0);
          }
          
          /*--- Physical boundaries definition ---*/
          
          if (Marker_Tag != "SEND_RECEIVE") {
            getline (mesh_file, text_line);
            text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
            if (duplicate)  nElem_Bound[iMarker+1]  = nElem_Bound[iMarker];

            if (rank == MASTER_NODE) {
              cout << nElem_Bound[iMarker]  << " boundary elements in index "<< iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;
              if (duplicate)  cout << nElem_Bound[iMarker+1]  << " boundary elements in index "<< iMarker+1 <<" (Marker = " <<Marker_Tag_Duplicate<< ")." << endl;
            }
            
            /*--- Allocate space for elements ---*/
            
            bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
            
            if (duplicate) bound[iMarker+1] = new CPrimalGrid* [nElem_Bound[iMarker+1]];

            nelem_edge_bound = 0; nelem_triangle_bound = 0; nelem_quad_bound = 0; ielem = 0;
            for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
              getline(mesh_file, text_line);
              istringstream bound_line(text_line);
              bound_line >> VTK_Type;
              switch(VTK_Type) {
                case LINE:
                  
                  if (nDim == 3) {
                    cout << "Please remove line boundary conditions from the mesh file!" << endl;
#ifndef HAVE_MPI
                    exit(EXIT_FAILURE);
#else
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Abort(MPI_COMM_WORLD,1);
                    MPI_Finalize();
#endif
                  }
                  
                  bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                  bound[iMarker][ielem] = new CLine(vnodes_edge[0], vnodes_edge[1],2);
                  
                  if (duplicate) {
                    if (ActDisk_Bool[vnodes_edge[0]]) { vnodes_edge[0] = ActDiskPoint_Back[vnodes_edge[0]]; }
                    if (ActDisk_Bool[vnodes_edge[1]]) { vnodes_edge[1] = ActDiskPoint_Back[vnodes_edge[1]]; }
                    bound[iMarker+1][ielem] = new CLine(vnodes_edge[0], vnodes_edge[1],2);
                  }
                  
                  ielem++; nelem_edge_bound++; break;
                  
                case TRIANGLE:
                  bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                  bound[iMarker][ielem] = new CTriangle(vnodes_triangle[0], vnodes_triangle[1], vnodes_triangle[2],3);
                  
                  if (duplicate) {
                    if (ActDisk_Bool[vnodes_triangle[0]]) { vnodes_triangle[0] = ActDiskPoint_Back[vnodes_triangle[0]]; }
                    if (ActDisk_Bool[vnodes_triangle[1]]) { vnodes_triangle[1] = ActDiskPoint_Back[vnodes_triangle[1]]; }
                    if (ActDisk_Bool[vnodes_triangle[2]]) { vnodes_triangle[2] = ActDiskPoint_Back[vnodes_triangle[2]]; }
                    bound[iMarker+1][ielem] = new CTriangle(vnodes_triangle[0], vnodes_triangle[1], vnodes_triangle[2],3);
                    
                  }
                  
                  ielem++; nelem_triangle_bound++; break;
                  
                case QUADRILATERAL:
                  
                  bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                  
                  bound[iMarker][ielem] = new CQuadrilateral(vnodes_quad[0], vnodes_quad[1], vnodes_quad[2], vnodes_quad[3],3);
                  
                  if (duplicate) {
                    if (ActDisk_Bool[vnodes_quad[0]]) { vnodes_quad[0] = ActDiskPoint_Back[vnodes_quad[0]]; }
                    if (ActDisk_Bool[vnodes_quad[1]]) { vnodes_quad[1] = ActDiskPoint_Back[vnodes_quad[1]]; }
                    if (ActDisk_Bool[vnodes_quad[2]]) { vnodes_quad[2] = ActDiskPoint_Back[vnodes_quad[2]]; }
                    if (ActDisk_Bool[vnodes_quad[3]]) { vnodes_quad[3] = ActDiskPoint_Back[vnodes_quad[3]]; }
                    bound[iMarker+1][ielem] = new CQuadrilateral(vnodes_quad[0], vnodes_quad[1], vnodes_quad[2], vnodes_quad[3],3);
                  }
                  
                  ielem++; nelem_quad_bound++;
                  
                  break;
                  
                  
              }
            }
            
            /*--- Update config information storing the boundary information in the right place ---*/
            
            Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
            config->SetMarker_All_TagBound(iMarker, Marker_Tag);
            config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
            config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
            config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
            config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
            config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
            config->SetMarker_All_Analyze(iMarker, config->GetMarker_CfgFile_Analyze(Marker_Tag));
            config->SetMarker_All_FSIinterface(iMarker, config->GetMarker_CfgFile_FSIinterface(Marker_Tag));
            config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
            config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
            config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
            config->SetMarker_All_SendRecv(iMarker, NONE);
            config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
            
            if (duplicate) {
              Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag_Duplicate)] = Marker_Tag_Duplicate;
              config->SetMarker_All_TagBound(iMarker+1, Marker_Tag_Duplicate);
              config->SetMarker_All_KindBC(iMarker+1, config->GetMarker_CfgFile_KindBC(Marker_Tag_Duplicate));
              config->SetMarker_All_Monitoring(iMarker+1, config->GetMarker_CfgFile_Monitoring(Marker_Tag_Duplicate));
              config->SetMarker_All_GeoEval(iMarker+1, config->GetMarker_CfgFile_GeoEval(Marker_Tag_Duplicate));
              config->SetMarker_All_Designing(iMarker+1, config->GetMarker_CfgFile_Designing(Marker_Tag_Duplicate));
              config->SetMarker_All_Plotting(iMarker+1, config->GetMarker_CfgFile_Plotting(Marker_Tag_Duplicate));
              config->SetMarker_All_Analyze(iMarker+1, config->GetMarker_CfgFile_Analyze(Marker_Tag_Duplicate));
              config->SetMarker_All_FSIinterface(iMarker+1, config->GetMarker_CfgFile_FSIinterface(Marker_Tag_Duplicate));
              config->SetMarker_All_DV(iMarker+1, config->GetMarker_CfgFile_DV(Marker_Tag_Duplicate));
              config->SetMarker_All_Moving(iMarker+1, config->GetMarker_CfgFile_Moving(Marker_Tag_Duplicate));
              config->SetMarker_All_PerBound(iMarker+1, config->GetMarker_CfgFile_PerBound(Marker_Tag_Duplicate));
              config->SetMarker_All_SendRecv(iMarker+1, NONE);
              config->SetMarker_All_Out_1D(iMarker+1, config->GetMarker_CfgFile_Out_1D(Marker_Tag_Duplicate));

              boundary_marker_count++;
              iMarker++;
              
            }
            
          }
          
          /*--- Send-Receive boundaries definition ---*/
          
          else {
            
            unsigned long nelem_vertex = 0, vnodes_vertex;
            unsigned short transform;
            getline (mesh_file, text_line);
            text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
            bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
            
            nelem_vertex = 0; ielem = 0;
            getline (mesh_file, text_line); text_line.erase (0,8);
            config->SetMarker_All_KindBC(iMarker, SEND_RECEIVE);
            config->SetMarker_All_SendRecv(iMarker, atoi(text_line.c_str()));
            
            for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
              getline(mesh_file, text_line);
              istringstream bound_line(text_line);
              bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;
              
              bound[iMarker][ielem] = new CVertexMPI(vnodes_vertex, nDim);
              bound[iMarker][ielem]->SetRotation_Type(transform);
              ielem++; nelem_vertex++;
            }
            
          }
          
          boundary_marker_count++;
          iMarker++;
          
        } while (iMarker < nMarker);
        
        if (boundary_marker_count == nMarker) break;
        
      }
    }
    
    while (getline (mesh_file, text_line)) {
      
      /*--- Read periodic transformation info (center, rotation, translation) ---*/
      
      position = text_line.find ("NPERIODIC=",0);
      if (position != string::npos) {
        unsigned short nPeriodic, iPeriodic, iIndex;
        
        /*--- Set bool signifying that periodic transormations were found ---*/
        found_transform = true;
        
        /*--- Read and store the number of transformations. ---*/
        text_line.erase (0,10); nPeriodic = atoi(text_line.c_str());
        if (rank == MASTER_NODE) {
          if (nPeriodic - 1 != 0)
            cout << nPeriodic - 1 << " periodic transformations." << endl;
        }
        config->SetnPeriodicIndex(nPeriodic);
        
        /*--- Store center, rotation, & translation in that order for each. ---*/
        for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
          getline (mesh_file, text_line);
          position = text_line.find ("PERIODIC_INDEX=",0);
          if (position != string::npos) {
            text_line.erase (0,15); iIndex = atoi(text_line.c_str());
            if (iIndex != iPeriodic) {
              cout << "PERIODIC_INDEX out of order in SU2 file!!" << endl;
#ifndef HAVE_MPI
              exit(EXIT_FAILURE);
#else
              MPI_Barrier(MPI_COMM_WORLD);
              MPI_Abort(MPI_COMM_WORLD,1);
              MPI_Finalize();
#endif
            }
          }
          su2double* center    = new su2double[3];
          su2double* rotation  = new su2double[3];
          su2double* translate = new su2double[3];
          getline (mesh_file, text_line);
          istringstream cent(text_line);
          cent >> center[0]; cent >> center[1]; cent >> center[2];
          config->SetPeriodicCenter(iPeriodic, center);
          getline (mesh_file, text_line);
          istringstream rot(text_line);
          rot >> rotation[0]; rot >> rotation[1]; rot >> rotation[2];
          config->SetPeriodicRotation(iPeriodic, rotation);
          getline (mesh_file, text_line);
          istringstream tran(text_line);
          tran >> translate[0]; tran >> translate[1]; tran >> translate[2];
          config->SetPeriodicTranslate(iPeriodic, translate);
          
          delete [] center; delete [] rotation; delete [] translate;
        }
      }
    }
    
    /*--- If no periodic transormations were found, store default zeros ---*/
    
    if (!found_transform) {
      unsigned short nPeriodic = 1, iPeriodic = 0;
      config->SetnPeriodicIndex(nPeriodic);
      su2double* center    = new su2double[3];
      su2double* rotation  = new su2double[3];
      su2double* translate = new su2double[3];
      for (unsigned short iDim = 0; iDim < 3; iDim++) {
        center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
      }
      config->SetPeriodicCenter(iPeriodic,    center);
      config->SetPeriodicRotation(iPeriodic,  rotation);
      config->SetPeriodicTranslate(iPeriodic, translate);
      delete [] center; delete [] rotation; delete [] translate;
    }
  
  /*--- Close the input file ---*/
  
  mesh_file.close();
  
  /*--- Release actuator disk memory ---*/

  if (actuator_disk) {
    delete [] ActDisk_Bool;
    delete [] ActDiskPoint_Back;
    delete [] VolumePoint_Inv;
    delete [] CoordXVolumePoint;
    delete [] CoordYVolumePoint;
    delete [] CoordZVolumePoint;
    delete [] CoordXActDisk;
    delete [] CoordYActDisk;
    delete [] CoordZActDisk;
  }
  
}

void CPhysicalGeometry::Read_CGNS_Format_Parallel(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone) {
  
  /*--- Original CGNS reader implementation by Thomas D. Economon,
   Francisco Palacios. Improvements for mixed-element meshes generated
   by ICEM added by Martin Spel (3D) & Shlomy Shitrit (2D), April 2014.
   Parallel version by Thomas D. Economon, February 2015. ---*/
  
#ifdef HAVE_CGNS
    
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short VTK_Type = 0, iMarker = 0;
  unsigned short nMarker_Max = config->GetnMarker_Max();
  unsigned long iPoint = 0, iProcessor = 0, ielem = 0, GlobalIndex = 0;
  unsigned long globalOffset = 0;
  nZone = val_nZone;
  
  /*--- Local variables needed when calling the CGNS mid-level API. ---*/
  
  unsigned long vnodes_cgns[8] = {0,0,0,0,0,0,0,0};
  su2double Coord_cgns[3] = {0.0,0.0,0.0};
  int fn, nbases = 0, nzones = 0, ngrids = 0, ncoords = 0, nsections = 0;
  int *vertices = NULL, *cells = NULL, nMarkers = 0, *boundVerts = NULL, npe;
  int interiorElems = 0, totalVerts = 0;
  int cell_dim = 0, phys_dim = 0, nbndry, parent_flag, file_type;
  char basename[CGNS_STRING_SIZE], zonename[CGNS_STRING_SIZE];
  char coordname[CGNS_STRING_SIZE];
  cgsize_t* cgsize; cgsize = new cgsize_t[3];
  ZoneType_t zonetype;
  DataType_t datatype;
  passivedouble** coordArray = NULL;
  passivedouble*** gridCoords = NULL;
  ElementType_t elemType;
  cgsize_t range_min, range_max, startE, endE;
  range_min = 1;
  string currentElem;
  int** elemTypeVTK = NULL;
  int** elemIndex = NULL;
  int** elemBegin = NULL;
  int** elemEnd = NULL;
  int** nElems = NULL;
  cgsize_t**** connElems = NULL;
  cgsize_t* connElemCGNS = NULL;
  cgsize_t* connElemTemp = NULL;
  cgsize_t ElementDataSize = 0;
  cgsize_t* parentData = NULL;
  int** dataSize = NULL;
  bool** isInternal = NULL;
  char*** sectionNames = NULL;
  
  /*--- Initialize counters for local/global points & elements ---*/
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  unsigned long Local_nElem;
  unsigned long Local_nElemTri, Local_nElemQuad, Local_nElemTet;
  unsigned long Local_nElemHex, Local_nElemPrism, Local_nElemPyramid;
  MPI_Request *send_req, *recv_req;
  MPI_Status  status;
  int ind;
#endif
  
  /*--- Initialize counters for local/global points & elements ---*/
  
  Global_nPoint  = 0; Global_nPointDomain = 0; Global_nElem = 0;
  nelem_edge     = 0; Global_nelem_edge     = 0;
  nelem_triangle = 0; Global_nelem_triangle = 0;
  nelem_quad     = 0; Global_nelem_quad     = 0;
  nelem_tetra    = 0; Global_nelem_tetra    = 0;
  nelem_hexa     = 0; Global_nelem_hexa     = 0;
  nelem_prism    = 0; Global_nelem_prism    = 0;
  nelem_pyramid  = 0; Global_nelem_pyramid  = 0;
  
  /*--- Initialize some additional counters for the parallel partitioning ---*/
  
  unsigned long total_pt_accounted = 0;
  unsigned long rem_points         = 0;
  unsigned long element_count      = 0;
  unsigned long element_remainder  = 0;
  unsigned long total_elems        = 0;
  
  /*--- Allocate memory for the linear partitioning of the mesh. These
   arrays are the size of the number of ranks. ---*/
  
  starting_node = new unsigned long[size];
  ending_node   = new unsigned long[size];
  npoint_procs  = new unsigned long[size];
  
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  unsigned long *nElem_Linear  = new unsigned long[size];
  
  unsigned long *elemB = new unsigned long[size];
  unsigned long *elemE = new unsigned long[size];
  
  unsigned long *elemGlobalID = NULL;
  
  unsigned short *nPoinPerElem = NULL;
  unsigned short *elemTypes = NULL;
  
  bool *isMixed = NULL;
  
  unsigned short connSize = 10;
  
  /*--- Check whether the supplied file is truly a CGNS file. ---*/
  if (cg_is_cgns(val_mesh_filename.c_str(), &file_type) != CG_OK) {
    if (rank == MASTER_NODE) {
    printf( "\n\n   !!! Error !!!\n" );
    printf( " %s is not a CGNS file.\n", val_mesh_filename.c_str());
    printf( " Now exiting...\n\n");
    }
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  /*--- Open the CGNS file for reading. The value of fn returned
   is the specific index number for this file and will be
   repeatedly used in the function calls. ---*/
  
  if (cg_open(val_mesh_filename.c_str(), CG_MODE_READ, &fn)) cg_error_exit();
  if (rank == MASTER_NODE) {
    cout << "Reading the CGNS file: ";
    cout << val_mesh_filename.c_str() << "." << endl;
  }
  
  /*--- Get the number of databases. This is the highest node
   in the CGNS heirarchy. ---*/
  
  if ( cg_nbases(fn, &nbases) ) cg_error_exit();
  if (rank == MASTER_NODE)
    cout << "CGNS file contains " << nbases << " database(s)." << endl;
  
  /*--- Check if there is more than one database. Throw an
   error if there is because this reader can currently
   only handle one database. ---*/
  
  if ( nbases > 1 ) {
    if (rank == MASTER_NODE) {
    printf("\n\n   !!! Error !!!\n" );
    printf("CGNS reader currently incapable of handling more than 1 database.");
    printf("Now exiting...\n\n");
    }
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  /*--- Read the databases. Note that the CGNS indexing starts at 1. ---*/
  
  for (int i = 1; i <= nbases; i++) {
    
    if (cg_base_read(fn, i, basename, &cell_dim, &phys_dim)) cg_error_exit();
    
    /*--- Get the number of zones for this base. ---*/
    
    if ( cg_nzones(fn, i, &nzones) ) cg_error_exit();
    if (rank == MASTER_NODE) {
      cout << "Database " << i << ", " << basename << ": " << nzones;
      cout << " zone(s), cell dimension of " << cell_dim << ", physical ";
      cout << "dimension of " << phys_dim << "." << endl;
    }
    
    /*--- Check if there is more than one zone. Throw an
     error if there is, because this reader can currently
     only handle one zone. This could be extended in the future. ---*/
    
    if ( nzones > 1 ) {
      if (rank == MASTER_NODE) {
      printf("\n\n   !!! Error !!!\n" );
      printf("CGNS reader currently incapable of handling more than 1 zone.");
      printf("Now exiting...\n\n");
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
    
    /*--- Initialize some data structures for  all zones. ---*/
    
    vertices     = new int[nzones];
    cells        = new int[nzones];
    boundVerts   = new int[nzones];
    coordArray   = new passivedouble*[nzones];
    gridCoords   = new passivedouble**[nzones];
    elemTypeVTK  = new int*[nzones];
    elemIndex    = new int*[nzones];
    elemBegin    = new int*[nzones];
    elemEnd      = new int*[nzones];
    nElems       = new int*[nzones];
    dataSize     = new int*[nzones];
    isInternal   = new bool*[nzones];
    nMarkers     = 0;
    sectionNames = new char**[nzones];
    connElems    = new cgsize_t***[nzones];
    
    /*--- Loop over all zones in this base. Again, indexing starts at 1. ---*/
    
    for (int j = 1; j <= nzones; j++) {
     
      connElems[j-1] = NULL;
 
      /*--- Read the basic information for this zone, including
       the name and the number of vertices, cells, and
       boundary cells which are stored in the cgsize variable. ---*/
      
      if (cg_zone_read(fn, i, j, zonename, cgsize)) cg_error_exit();
      
      /*--- Rename the zone size information for clarity.
       NOTE: The number of cells here may be only the number of
       interior elements or it may be the total. This needs to
       be counted explicitly later. ---*/
      
      vertices[j-1]   = cgsize[0];
      cells[j-1]      = cgsize[1];
      boundVerts[j-1] = cgsize[2];
      
      /*--- Increment the total number of vertices from all zones. ---*/
      
      nPoint       = vertices[j-1];
      nPointDomain = vertices[j-1];
      
      Global_nPoint       = vertices[j-1];
      Global_nPointDomain = vertices[j-1];
      
      totalVerts += vertices[j-1];
      
      /*--- Print some information about the current zone. ---*/
      
      if (cg_zone_type(fn, i, j, &zonetype)) cg_error_exit();
      if (rank == MASTER_NODE) {
        cout << "Zone " << j << ", " << zonename << ": " << vertices[j-1];
        cout << " vertices, " << cells[j-1] << " cells, " << boundVerts[j-1];
        cout << " boundary vertices." << endl;
      }
      
      /*--- Retrieve the number of grids in this zone. For now, we know
       this is one, but to be more general, this will need to check and
       allow for a loop over all grids. ---*/
      
      if (cg_ngrids(fn, i, j, &ngrids)) cg_error_exit();
      if (ngrids > 1) {
        if (rank == MASTER_NODE) {
          printf("\n\n   !!! Error !!!\n" );
          printf("CGNS reader currently handles only 1 grid per zone.");
          printf("Now exiting...\n\n");
        }
#ifndef HAVE_MPI
        exit(EXIT_FAILURE);
#else
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#endif
      }
      
      /*--- Check the number of coordinate arrays stored in this zone.
       Should be 2 for 2-D grids and 3 for 3-D grids. ---*/
      
      if (cg_ncoords( fn, i, j, &ncoords)) cg_error_exit();
      if (rank == MASTER_NODE) {
        cout << "Reading grid coordinates." << endl;
        cout << "Number of coordinate dimensions is " << ncoords << "." << endl;
      }
      
      /*--- Compute the number of points that will be on each processor.
       This is a linear partitioning with the addition of a simple load
       balancing for any remainder points. ---*/
      
      total_pt_accounted = 0;
      for (int ii = 0; ii < size; ii++) {
        npoint_procs[ii] = vertices[j-1]/size;
        total_pt_accounted = total_pt_accounted + npoint_procs[ii];
      }
      
      /*--- Get the number of remainder points after the even division ---*/
      
      rem_points = vertices[j-1]-total_pt_accounted;
      for (unsigned long ii = 0; ii < rem_points; ii++) {
        npoint_procs[ii]++;
      }
      
      /*--- Store the local number of nodes and the beginning/end index ---*/
      
      nPoint = npoint_procs[rank];
      starting_node[0] = 0;
      ending_node[0]   = starting_node[0] + npoint_procs[0];
      nPoint_Linear[0] = 0;
      for (int ii = 1; ii < size; ii++) {
        starting_node[ii] = ending_node[ii-1];
        ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
        nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
      }
      nPoint_Linear[size] = vertices[j-1];
      
      /*--- Set the value of range_max to the total number of nodes in
       the unstructured mesh. Also allocate memory for the temporary array
       that will hold the grid coordinates as they are extracted. Note the
       +1 for CGNS convention. ---*/
      
      range_min = (cgsize_t)starting_node[rank]+1;
      range_max = (cgsize_t)ending_node[rank];
      coordArray[j-1] = new passivedouble[nPoint];
      
      /*--- Allocate memory for the 2-D array that will store the x, y,
       & z (if required) coordinates for writing into the SU2 mesh. ---*/
      
      gridCoords[j-1] = new passivedouble*[ncoords];
      for (int ii = 0; ii < ncoords; ii++) {
        *(gridCoords[j-1]+ii) = new passivedouble[nPoint];
      }
      
      /*--- Loop over each set of coordinates. Note again
       that the indexing starts at 1. ---*/
      
      for (int k = 1; k <= ncoords; k++) {
        
        /*--- Read the coordinate info. This will retrieve the
         data type (either RealSingle or RealDouble) as
         well as the coordname which will specifiy the
         type of data that it is based in the SIDS convention.
         This might be "CoordinateX," for instance. ---*/
        
        if (cg_coord_info(fn, i, j, k, &datatype, coordname))
          cg_error_exit();
        if (rank == MASTER_NODE) {
          cout << "Loading " << coordname;
          if (size > SINGLE_NODE) {
            cout << " values into linear partitions." << endl;
          } else {
            cout << " values." << endl;
          }
        }
        
        /*--- Always retrieve the grid coords in su2double precision. ---*/
        
        if (datatype != RealDouble) {
          printf("\n\n   !!! Error !!!\n" );
          printf(" CGNS coordinates are not su2double precision.\n");
          printf(" Now exiting...\n\n");
#ifndef HAVE_MPI
          exit(EXIT_FAILURE);
#else
          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Abort(MPI_COMM_WORLD,1);
          MPI_Finalize();
#endif
        }
        if ( cg_coord_read(fn, i, j, coordname, datatype, &range_min,
                           &range_max, coordArray[j-1]) ) cg_error_exit();
        
        /*--- Copy these coords into the array for storage until
         writing the SU2 mesh. ---*/
        
        for (unsigned long m = 0; m < nPoint; m++ ) {
          gridCoords[j-1][k-1][m] = coordArray[j-1][m];
        }
        
      }
      
      /*--- Begin section for retrieving the connectivity info. ---*/
      
      if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
        cout << "Distributing connectivity across all ranks." << endl;
      
      /*--- First check the number of sections. ---*/
      
      if ( cg_nsections(fn, i, j, &nsections) ) cg_error_exit();
      if (rank == MASTER_NODE) {
        cout << "Number of connectivity sections is ";
        cout << nsections << "." << endl;
      }
      
      /*--- Allocate several data structures to hold the various
       pieces of information describing each section. It is
       stored in this manner so that it can be written to
       SU2 memory later. ---*/
      
      elemTypeVTK[j-1] = new int[nsections];
      elemIndex[j-1]   = new int[nsections];
      elemBegin[j-1]   = new int[nsections];
      elemEnd[j-1]     = new int[nsections];
      nElems[j-1]      = new int[nsections];
      dataSize[j-1]    = new int[nsections];
      isInternal[j-1]  = new bool[nsections];
      
      sectionNames[j-1] = new char*[nsections];
      for (int ii = 0; ii < nsections; ii++) {
        sectionNames[j-1][ii]= new char[CGNS_STRING_SIZE];
      }
      
      connElems[j-1] = new cgsize_t**[nsections];
      
      /*--- Loop over each section. This will include the main
       connectivity information for the grid cells, as well
       as any boundaries which were labeled before export. ---*/
      
      for (int s = 1; s <= nsections; s++) {
      
        connElems[j-1][s-1] = NULL; 
        /*--- Read the connectivity details for this section.
         Store the total number of elements in this section
         to be used later for memory allocation. ---*/
        
        if (cg_section_read(fn, i, j, s, sectionNames[j-1][s-1],
                            &elemType, &startE, &endE, &nbndry,
                            &parent_flag)) cg_error_exit();
        
        /*--- Store the beginning and ending index for this section. ---*/
        
        elemBegin[j-1][s-1] = (int)startE;
        elemEnd[j-1][s-1]   = (int)endE;
        
        /*--- Compute element linear partitioning ---*/
        
        element_count = (int) (endE-startE+1);
        total_elems = 0;
        for (int ii = 0; ii < size; ii++) {
          nElem_Linear[ii] = element_count/size;
          total_elems += nElem_Linear[ii];
        }
        
        /*--- Get the number of remainder elements after even division ---*/
        
        element_remainder = element_count-total_elems;
        for (unsigned long ii = 0; ii < element_remainder; ii++) {
          nElem_Linear[ii]++;
        }
        
        /*--- Store the number of elements that this rank is responsible for
         in the current section. ---*/
        
        nElems[j-1][s-1] = (int)nElem_Linear[rank];
        
        /*--- Get starting and end element index for my rank. ---*/
        
        elemB[0] = startE;
        elemE[0] = startE + nElem_Linear[0] - 1;
        for (unsigned long ii = 1; ii < (unsigned long)size; ii++) {
          elemB[ii] = elemE[ii-1]+1;
          elemE[ii] = elemB[ii] + nElem_Linear[ii] - 1;
        }
        
        /*--- Allocate some memory for the handling the connectivity
         and auxiliary data that we are need to communicate. ---*/
        
        connElemCGNS = new cgsize_t[nElems[j-1][s-1]*connSize];
        nPoinPerElem = new unsigned short[nElems[j-1][s-1]];
        elemGlobalID = new unsigned long[nElems[j-1][s-1]];
        elemTypes    = new unsigned short[nElems[j-1][s-1]];
        
        isMixed = new bool[nElems[j-1][s-1]];
        for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) isMixed[ii] = false;

        /*--- Protect against the situation where there are fewer elements
        in a section than number of ranks, or the linear partitioning will
        fail. For now, assume that these must be surfaces, and we will 
        avoid a parallel read and have the master read this section (the
        master processes all of the markers anyway). ---*/

        if (nElems[j-1][s-1] < rank+1) {

          isInternal[j-1][s-1] = false;

        } else {        

        /*--- Retrieve the connectivity information and store. Note that
         we are only accessing our rank's piece of the data here in the
         partial read function in the CGNS API. ---*/

        if (cg_elements_partial_read(fn, i, j, s, (cgsize_t)elemB[rank],
                                    (cgsize_t)elemE[rank], connElemCGNS,
                                    parentData) != CG_OK) cg_error_exit();
        
        /*--- Find the number of nodes required to represent
         this type of element. ---*/
        
        ElementType_t elmt_type;
        if (cg_npe(elemType, &npe)) cg_error_exit();
        
        /*--- Loop through all of the elements in this section to get more
         information and to decide whether it has internal elements. ---*/
        
        int counter = 0;
        for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) {
          
          /*--- If we have a mixed element section, we need to check the elem
           type one by one. Set the flag to true if mixed. ---*/
          
          if (elemType == MIXED) {
            elmt_type = ElementType_t(connElemCGNS[counter]);
            cg_npe(elmt_type, &npe);
            counter++; for ( int jj = 0; jj < npe; jj++ ) counter++;
            isMixed[ii] = true;
          } else {
            elmt_type = elemType;
          }
          
          /*--- Store the number of verts per elem for the current elem. ---*/
          
          nPoinPerElem[ii] = npe;
          
          /*--- Store the global ID for this element. Note the -1 to move
           from CGNS convention to SU2 convention. We also subtract off
           an additional offset in case we have found boundary sections
           prior to this one, in order to keep the internal element global
           IDs indexed starting from zero. ---*/
          
          elemGlobalID[ii] = elemB[rank] + ii - 1 - globalOffset;
          
          /*--- Need to check the element type and correctly specify the
           VTK identifier for that element. SU2 recognizes elements by
           their VTK number. ---*/
          
          switch (elmt_type) {
            case NODE:
              currentElem   = "Vertex";
              elemTypes[ii] = 1;
              break;
            case BAR_2:
              currentElem   = "Line";
              elemTypes[ii] = 3;
              break;
            case BAR_3:
              currentElem   = "Line";
              elemTypes[ii] = 3;
              break;
            case TRI_3:
              currentElem   = "Triangle";
              elemTypes[ii] = 5;
              break;
            case QUAD_4:
              currentElem   = "Quadrilateral";
              elemTypes[ii] = 9;
              break;
            case TETRA_4:
              currentElem   = "Tetrahedron";
              elemTypes[ii] = 10;
              break;
            case HEXA_8:
              currentElem   = "Hexahedron";
              elemTypes[ii] = 12;
              break;
            case PENTA_6:
              currentElem   = "Prism";
              elemTypes[ii] = 13;
              break;
            case PYRA_5:
              currentElem   = "Pyramid";
              elemTypes[ii] = 14;
              break;
            case HEXA_20:
              if (rank == MASTER_NODE) {
                printf("\n\n   !!! Error !!!\n" );
                printf(" HEXA-20 element type not supported\n");
                printf(" Section %d, npe=%d\n", s, npe);
                printf(" startE %d, endE %d\n", (int)startE, (int)endE);
                printf(" Now exiting...\n\n");
              }
#ifndef HAVE_MPI
              exit(EXIT_FAILURE);
#else
              MPI_Barrier(MPI_COMM_WORLD);
              MPI_Abort(MPI_COMM_WORLD,1);
              MPI_Finalize();
#endif
              break;
            default:
              if (rank == MASTER_NODE) {
                printf("\n\n   !!! Error !!!\n" );
                printf(" Unknown elem: (type %d, npe=%d)\n", elemType, npe);
                printf(" Section %d\n", s);
                printf(" startE %d, endE %d\n", (int)startE, (int)endE);
                printf(" Now exiting...\n\n");
              }
#ifndef HAVE_MPI
              exit(EXIT_FAILURE);
#else
              MPI_Abort(MPI_COMM_WORLD,1);
              MPI_Finalize();
#endif
              break;
          }
          
          /*--- Check if the elements in this section are part
           of the internal domain or are part of the boundary
           surfaces. This will be used to separate the
           internal connectivity from the boundary connectivity.
           We will check for quad and tri elements for 3-D meshes
           because these will be the boundaries. Similarly, line
           elements will be boundaries to 2-D problems. ---*/
          
          if ( cell_dim == 2 ) {
            
            /*--- In 2-D check for line elements, VTK type 3. ---*/
            
            if (elemTypes[ii] == 3) {
              isInternal[j-1][s-1] = false;
            } else {
              isInternal[j-1][s-1] = true;
              interiorElems++;
            }
            
          } else if (cell_dim == 3) {
            
            /*--- In 3-D check for tri/quad elements, VTK types 5 or 9. ---*/
            
            switch (elemTypes[ii]) {
              case 5:
              case 9:
                isInternal[j-1][s-1] = false;
                break;
              default:
                isInternal[j-1][s-1] = true;
                interiorElems++;
                break;
            }
            
          }
        }
        
        /*--- Print some information to the console. ---*/
        
        if (rank == MASTER_NODE) {
          for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ )
            if (isMixed[ii]) {currentElem = "Mixed"; break;}
          cout << "Loading section " << sectionNames[j-1][s-1];
          cout << " of element type " << currentElem << "." << endl;
        }
       
        } 
        
         /*--- If we have found that this is a boundary section (we assume
         that internal cells and boundary cells do not exist in the same
         section together), the master node reads the boundary section.
         Otherwise, we have all ranks read and communicate the internals. ---*/
        
        if (!isInternal[j-1][s-1]) {
          
          /*--- Master node should read this entire marker section. Free
           the memory for the conn. from the CGNS file since we are going
           to read the section again with the master. ---*/
          
          delete [] connElemCGNS;
          delete [] nPoinPerElem;
          delete [] elemTypes;
          delete [] elemGlobalID;
          delete [] isMixed;
          
          /*--- Since we found an internal section, we should adjust the
           element global ID offset by the total size of the section. ---*/
          
          globalOffset += element_count;
          
          if (rank == MASTER_NODE) {
            
            /*--- First increment the markers ---*/
            
            nMarkers++;
            
            /*--- Read the section info again ---*/
            
            if ( cg_section_read(fn, i, j, s, sectionNames[j-1][s-1],
                                 &elemType, &startE, &endE, &nbndry,
                                 &parent_flag) ) cg_error_exit();
            
            /*--- Store the number of elems (all on the master). ---*/
            
            nElems[j-1][s-1] = (int) (endE-startE+1);
            
            /*--- Read and store the total amount of data that will be
             listed when reading this section. ---*/
            
            if (cg_ElementDataSize(fn, i, j, s, &ElementDataSize))
              cg_error_exit();
            dataSize[j-1][s-1] = ElementDataSize;
            
            /*--- Find the number of nodes required to represent
             this type of element. ---*/
            
            if (cg_npe(elemType, &npe)) cg_error_exit();
            elemIndex[j-1][s-1] = npe;
            
            /*--- Need to check the element type and correctly
             specify the VTK identifier for that element.
             SU2 recognizes elements by their VTK number. ---*/
            
            switch (elemType) {
              case NODE:
                elemTypeVTK[j-1][s-1] = 1;
                break;
              case BAR_2:
                elemTypeVTK[j-1][s-1] = 3;
                break;
              case BAR_3:
                elemTypeVTK[j-1][s-1] = 3;
                break;
              case TRI_3:
                elemTypeVTK[j-1][s-1] = 5;
                break;
              case QUAD_4:
                elemTypeVTK[j-1][s-1] = 9;
                break;
              case TETRA_4:
                elemTypeVTK[j-1][s-1] = 10;
                break;
              case HEXA_8:
                elemTypeVTK[j-1][s-1] = 12;
                break;
              case PENTA_6:
                elemTypeVTK[j-1][s-1] = 13;
                break;
              case PYRA_5:
                elemTypeVTK[j-1][s-1] = 14;
                break;
              case HEXA_20:
                printf( "\n\n   !!! Error !!!\n" );
                printf( " HEXA-20 element type not supported\n");
                printf(" Section %d, npe=%d\n", s, npe);
                printf(" startE %d, endE %d\n", (int)startE, (int)endE);
                printf( " Now exiting...\n\n");
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Abort(MPI_COMM_WORLD,1);
                MPI_Finalize();
#endif
                break;
              case MIXED:
                currentElem = "Mixed";
                elemTypeVTK[j-1][s-1] = -1;
                break;
              default:
                printf( "\n\n   !!! Error !!!\n" );
                printf( " Unknown elem: (type %d, npe=%d)\n", elemType, npe);
                printf(" Section %d\n", s);
                printf(" startE %d, endE %d\n", (int)startE, (int)endE);
                printf( " Now exiting...\n\n");
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD,1);
                MPI_Finalize();
#endif
                break;
            }
            
            /*--- In case of mixed data type, allocate place for 8 nodes
             maximum (hex), plus element type. ---*/
            
            if (elemTypeVTK[j-1][s-1] == -1) elemIndex[j-1][s-1] = 9;
            
            /*--- Allocate memory for accessing the connectivity and to
             store it in the proper data structure for post-processing. ---*/
            
            connElemTemp = new cgsize_t[dataSize[j-1][s-1]];
            connElems[j-1][s-1] = new cgsize_t*[elemIndex[j-1][s-1]];
            for (int jj = 0; jj < elemIndex[j-1][s-1]; jj++) {
              connElems[j-1][s-1][jj] = new cgsize_t[nElems[j-1][s-1]];
            }
            
            /*--- Retrieve the connectivity information and store. ---*/
            
            if (cg_elements_read(fn, i, j, s, connElemTemp, parentData))
              cg_error_exit();
            
            /*--- Copy these values into the larger array for
             storage until writing the SU2 file. ---*/
            
            if (elemTypeVTK[j-1][s-1] == -1) {
              int counter = 0;
              for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) {
                ElementType_t elmt_type = ElementType_t(connElemTemp[counter]);
                cg_npe( elmt_type, &npe);
                counter++;
                connElems[j-1][s-1][0][ii] = elmt_type;
                for ( int jj = 0; jj < npe; jj++ ) {
                  connElems[j-1][s-1][jj+1][ii] = connElemTemp[counter] - 1;
                  counter++;
                }
              }
            } else {
              int counter = 0;
              for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) {
                for ( int jj = 0; jj < elemIndex[j-1][s-1]; jj++ ) {
                  connElems[j-1][s-1][jj][ii] = connElemTemp[counter] - 1;
                  counter++;
                }
              }
            }
            delete[] connElemTemp;
          
          } // end master
          
        } else {
          
          /*--- These are internal elems. Allocate memory on each proc. ---*/
          
          connElemTemp = new cgsize_t[nElems[j-1][s-1]*connSize];
          
           /*--- Copy these values into the larger array for
           storage until writing the SU2 file. ---*/

          int counterTemp = 0, counterCGNS = 0;
          for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) {
            
            /*--- Store the conn in chunks of connSize for simplicity. ---*/
            
            counterTemp = ii*connSize;
            
            /*--- Store the connectivity values. Note we subtract one from
             the CGNS 1-based convention. We may also need to remove the first
             entry is this is a mixed element section. ---*/
            
            if (isMixed[ii]) counterCGNS++;
            for ( int jj = 0; jj < nPoinPerElem[ii]; jj++) {
              connElemTemp[counterTemp] = connElemCGNS[counterCGNS + jj] - 1;
              counterTemp++;
            }
            counterCGNS += nPoinPerElem[ii];
            
          }

          /*--- Free the memory for the conn. from the CGNS file. ---*/
          
          delete [] connElemCGNS;
          delete [] isMixed;
          
          /*--- We now have the connectivity stored in linearly partitioned
           chunks. We need to loop through and decide how many elements we
           must send to each rank in order to have all elements that
           surround a particular "owned" node on each rank (i.e., elements
           will appear on multiple ranks). First, initialize a counter
           and flag. ---*/
          
          int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
          int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
          int *nElem_Flag = new int[size];
          
          for (int ii=0; ii < size; ii++) {
            nElem_Send[ii] = 0;
            nElem_Recv[ii] = 0;
            nElem_Flag[ii]= -1;
          }
          nElem_Send[size] = 0; nElem_Recv[size] = 0;
          
          for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) {
            for ( int jj = 0; jj < nPoinPerElem[ii]; jj++ ) {
              
              /*--- Get the index of the current point. ---*/
              
              iPoint = connElemTemp[ii*connSize + jj];
              
              /*--- Search for the processor that owns this point ---*/
              
              iProcessor = iPoint/npoint_procs[0];
              if (iProcessor >= (unsigned long)size) iProcessor = (unsigned long)size-1;
              if (iPoint >= nPoint_Linear[iProcessor])
                while(iPoint >= nPoint_Linear[iProcessor+1]) iProcessor++;
              else
                while(iPoint <  nPoint_Linear[iProcessor])   iProcessor--;
              
              /*--- If we have not visted this element yet, increment our
               number of elements that must be sent to a particular proc. ---*/
              
              if (nElem_Flag[iProcessor] != ii) {
                nElem_Flag[iProcessor] = ii;
                nElem_Send[iProcessor+1]++;
              }
              
            }
          }
          
          /*--- Communicate the number of cells to be sent/recv'd amongst
           all processors. After this communication, each proc knows how
           many cells it will receive from each other processor. ---*/
          
#ifdef HAVE_MPI
          MPI_Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                       &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
          nElem_Recv[1] = nElem_Send[1];
#endif
          
          /*--- Prepare to send connectivities. First check how many
           messages we will be sending and receiving. Here we also put
           the counters into cumulative storage format to make the
           communications simpler. ---*/
          
          int nSends = 0, nRecvs = 0;
          for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
          
          for (int ii = 0; ii < size; ii++) {
            
            if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
            if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
            
            nElem_Send[ii+1] += nElem_Send[ii];
            nElem_Recv[ii+1] += nElem_Recv[ii];
          }

          /*--- Allocate memory to hold the connectivity that we are
           sending. Note that we are also sending the VTK element type
           in the first position and also the global ID. We have assumed
           a constant message size of a hex element + 2 extra vals. ---*/
          
          unsigned long *connSend = NULL;
          connSend = new unsigned long[connSize*nElem_Send[size]];
          for (int ii = 0; ii < connSize*nElem_Send[size]; ii++)
            connSend[ii] = 0;
          
          /*--- Create an index variable to keep track of our index
           position as we load up the send buffer. ---*/
          
          unsigned long *index = new unsigned long[size];
          for (int ii=0; ii < size; ii++) index[ii] = connSize*nElem_Send[ii];
          
          /*--- Loop through our elements and load the elems and their
           additional data that we will send to the other procs. ---*/
          
          for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) {
            for ( int jj = 0; jj < nPoinPerElem[ii]; jj++ ) {
              
              /*--- Get the index of the current point. ---*/
              
              iPoint = connElemTemp[ii*connSize + jj];
              
              /*--- Search for the processor that owns this point ---*/
              
              iProcessor = iPoint/npoint_procs[0];
              if (iProcessor >= (unsigned long)size) iProcessor = (unsigned long)size-1;
              if (iPoint >= nPoint_Linear[iProcessor])
                while(iPoint >= nPoint_Linear[iProcessor+1]) iProcessor++;
              else
                while(iPoint <  nPoint_Linear[iProcessor])   iProcessor--;
              
              /*--- Load connectivity into the buffer for sending ---*/
              
              if (nElem_Flag[iProcessor] != ii) {
                
                nElem_Flag[iProcessor] = ii;
                unsigned long nn = index[iProcessor];
                
                /*--- Load the VTK type first into the conn array,
                 then the connectivity vals, and last, the global ID. ---*/
                
                connSend[nn] = elemTypes[ii]; nn++;
                for ( int kk = 0; kk < nPoinPerElem[ii]; kk++ ) {
                  connSend[nn] = connElemTemp[ii*connSize + kk]; nn++;
                }
                connSend[nn] = (cgsize_t)elemGlobalID[ii];
                
                /*--- Increment the index by the message length ---*/
                
                index[iProcessor] += connSize;
                
              }
            }
          }

          /*--- Free memory after loading up the send buffer. ---*/
          
          delete [] connElemTemp;
          delete [] elemTypes;
          delete [] nPoinPerElem;
          delete [] elemGlobalID;
          delete [] index;
          
          /*--- Allocate the memory that we need for receiving the conn
           values and then cue up the non-blocking receives. Note that
           we do not include our own rank in the communications. We will
           directly copy our own data later. ---*/
          
          unsigned long *connRecv = NULL;
          connRecv = new unsigned long[connSize*nElem_Recv[size]];
          for (int ii = 0; ii < connSize*nElem_Recv[size]; ii++)
            connRecv[ii] = 0;
            
#ifdef HAVE_MPI
          send_req = new MPI_Request[nSends];
          recv_req = new MPI_Request[nRecvs];
          unsigned long iMessage = 0;
          for (int ii=0; ii<size; ii++) {
            if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
              int ll     = connSize*nElem_Recv[ii];
              int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
              int count  = connSize*kk;
              int source = ii;
              int tag    = ii + 1;
              SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                        MPI_COMM_WORLD, &(recv_req[iMessage]));
              iMessage++;
            }
          }
          
          /*--- Launch the non-blocking sends of the connectivity. ---*/
          
          iMessage = 0;
          for (int ii=0; ii<size; ii++) {
            if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
              int ll = connSize*nElem_Send[ii];
              int kk = nElem_Send[ii+1] - nElem_Send[ii];
              int count  = connSize*kk;
              int dest = ii;
              int tag    = rank + 1;
              SU2_MPI::Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                        MPI_COMM_WORLD, &(send_req[iMessage]));
              iMessage++;
            }
          }
#endif
          
          /*--- Copy my own rank's data into the recv buffer directly. ---*/
          
          int mm = connSize*nElem_Recv[rank];
          int ll = connSize*nElem_Send[rank];
          int kk = connSize*nElem_Send[rank+1];

          for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];
          
          /*--- Wait for the non-blocking sends and recvs to complete ---*/
          
#ifdef HAVE_MPI
          int number = nSends;
          for (int ii = 0; ii < nSends; ii++)
            SU2_MPI::Waitany(number, send_req, &ind, &status);
          
          number = nRecvs;
          for (int ii = 0; ii < nRecvs; ii++)
            SU2_MPI::Waitany(number, recv_req, &ind, &status);
          
          delete [] send_req;
          delete [] recv_req;
#endif
          
          /*--- Store the connectivity for this rank in the proper data
           structure before post-processing below. First, allocate the
           appropriate amount of memory for this section. ---*/

          connElems[j-1][s-1] = new cgsize_t*[connSize];
          for (int jj = 0; jj < connSize; jj++) {
            connElems[j-1][s-1][jj] = new cgsize_t[nElem_Recv[size]];
          }
          for (int ii = 0; ii < nElem_Recv[size]; ii++) {
            for (int jj = 0; jj < connSize; jj++) {
              connElems[j-1][s-1][jj][ii] = (cgsize_t)connRecv[ii*connSize+jj];
            }
          }
          
          /*--- Store the total number of elements I now have for
           the current section after completing the communications. ---*/
          
          nElems[j-1][s-1] = nElem_Recv[size];

          /*--- Free temporary memory from communications ---*/
          
          delete [] connSend;
          delete [] connRecv;
          delete [] nElem_Recv;
          delete [] nElem_Send;
          delete [] nElem_Flag;
          
        }
        
      } // end section
      
    } // end zone
    
  } // end database
  
  /*--- Close the CGNS file. ---*/
  
  if ( cg_close(fn) ) cg_error_exit();
  if (rank == MASTER_NODE)
    cout << "Successfully closed the CGNS file." << endl;
  
  /*--- Load the data from the CGNS file into SU2 memory. ---*/
  
  if (rank == MASTER_NODE)
    cout << endl << "Loading CGNS data into SU2 data structures." << endl;
  
  /*--- Read the dimension of the problem ---*/
  
  nDim = cell_dim;
  if (rank == MASTER_NODE) {
    if (nDim == 2) cout << "Two dimensional problem." << endl;
    if (nDim == 3) cout << "Three dimensional problem." << endl;
  }
  
  /*--- Initialize an array for the adjacency information (ParMETIS). ---*/
  
  vector< vector<unsigned long> > adj_nodes(nPoint, vector<unsigned long>(0));
  
  /*--- Loop to check total number of elements we have locally. ---*/
  
  ielem = 0;
  for (int k = 0; k < nzones; k++) {
    for (int s = 0; s < nsections; s++) {
      if (isInternal[k][s]) {
        for ( int i = 0; i < nElems[k][s]; i++) {
          ielem++;
        }
      }
    }
  }
  nElem = ielem;
  
  /*--- Store the total number of interior elements (global). ---*/
  
#ifdef HAVE_MPI
  Local_nElem = interiorElems;
  SU2_MPI::Allreduce(&Local_nElem, &Global_nElem, 1, MPI_UNSIGNED_LONG,
                MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nElem = interiorElems;
  nElem        = Global_nElem;
#endif
  
  if ((rank == MASTER_NODE) && (size > SINGLE_NODE)) {
    cout << Global_nElem << " interior elements before linear partitioning." << endl;
  } else if (rank == MASTER_NODE) {
    cout << Global_nElem << " interior elements." << endl;
  }
  
  /*--- Set up the global to local element mapping. ---*/
  Global_to_Local_Elem.clear();
  
  /*--- Allocate space for elements. We allocate enough for all interior
   elements globally, but we will only instantiate our local set. ---*/
  
  elem = new CPrimalGrid*[nElem];
  ielem = 0;
  unsigned long global_id = 0;
  
  /*--- Loop over all the internal, local volumetric elements. ---*/
  
  for (int k = 0; k < nzones; k++) {
    for (int s = 0; s < nsections; s++) {
      if (isInternal[k][s]) {
        for ( int i = 0; i < nElems[k][s]; i++) {
          
          /*--- Get the VTK type for this element. This is stored in the
           first entry of the connectivity structure. ---*/
          
          VTK_Type = connElems[k][s][0][i];
          
          /*--- Instantiate this element and build adjacency structure. ---*/
          
          switch(VTK_Type) {
              
            case TRIANGLE:
              
              for ( unsigned short j = 0; j < N_POINTS_TRIANGLE; j++ ) {
                vnodes_cgns[j] = connElems[k][s][j+1][i];
              }
              global_id = connElems[k][s][N_POINTS_TRIANGLE+1][i];
              for (unsigned short ii=0; ii<N_POINTS_TRIANGLE; ii++) {
                if ((vnodes_cgns[ii]>=starting_node[rank])&&(vnodes_cgns[ii]<ending_node[rank])) {
                  for (unsigned short j=0; j<N_POINTS_TRIANGLE; j++) {
                    if (ii!=j) {
                      adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[j]);
                    }
                  }
                }
              }
              Global_to_Local_Elem[global_id]=ielem;
              elem[ielem] = new CTriangle(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], nDim);
              ielem++; nelem_triangle++;
              break;
              
            case QUADRILATERAL:
              
              for ( unsigned short j = 0; j < N_POINTS_QUADRILATERAL; j++ ) {
                vnodes_cgns[j] = connElems[k][s][j+1][i];
              }
              global_id = connElems[k][s][N_POINTS_QUADRILATERAL+1][i];
              
              for (unsigned short ii=0; ii<N_POINTS_QUADRILATERAL; ii++) {
                if ((vnodes_cgns[ii]>=starting_node[rank])&&(vnodes_cgns[ii]<ending_node[rank])) {
                  
                  /*--- Build adjacency assuming the VTK connectivity ---*/
                  
                  adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+1)%4]);
                  adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+3)%4]);
                  
                }
              }
              
              Global_to_Local_Elem[global_id]=ielem;
              elem[ielem] = new CQuadrilateral(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], nDim);
              ielem++; nelem_quad++;
              break;
              
            case TETRAHEDRON:
              
              for ( unsigned short j = 0; j < N_POINTS_TETRAHEDRON; j++ ) {
                vnodes_cgns[j] = connElems[k][s][j+1][i];
              }
              global_id = connElems[k][s][N_POINTS_TETRAHEDRON+1][i];
              for (unsigned short ii=0; ii<N_POINTS_TETRAHEDRON; ii++) {
                if ((vnodes_cgns[ii]>=starting_node[rank])&&(vnodes_cgns[ii]<ending_node[rank])) {
                  for (unsigned short j=0; j<N_POINTS_TETRAHEDRON; j++) {
                    if (ii!=j) {
                      adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[j]);
                    }
                  }
                }
              }
              Global_to_Local_Elem[global_id]=ielem;
              elem[ielem] = new CTetrahedron(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3]);
              ielem++; nelem_tetra++;
              break;
              
            case HEXAHEDRON:
              
              for ( unsigned short j = 0; j < N_POINTS_HEXAHEDRON; j++ ) {
                vnodes_cgns[j] = connElems[k][s][j+1][i];
              }
              global_id = connElems[k][s][N_POINTS_HEXAHEDRON+1][i];
              
              for (unsigned short ii=0; ii<N_POINTS_HEXAHEDRON; ii++) {
                if ((vnodes_cgns[ii]>=starting_node[rank])&&(vnodes_cgns[ii]<ending_node[rank])) {
                  
                  /*--- Build adjacency assuming the VTK connectivity ---*/
                  
                  if (ii < 4) {
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+1)%4]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+3)%4]);
                  } else {
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii-3)%4+4]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii-1)%4+4]);
                  }
                  adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+4)%8]);

                }
              }
              
              Global_to_Local_Elem[global_id]=ielem;
              elem[ielem] = new CHexahedron(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], vnodes_cgns[4], vnodes_cgns[5], vnodes_cgns[6], vnodes_cgns[7]);
              ielem++; nelem_hexa++;
              break;
              
            case PRISM:
              
              for ( unsigned short j = 0; j < N_POINTS_PRISM; j++ ) {
                vnodes_cgns[j] = connElems[k][s][j+1][i];
              }
              global_id = connElems[k][s][N_POINTS_PRISM+1][i];
              
              for (unsigned short ii=0; ii<N_POINTS_PRISM; ii++) {
                if ((vnodes_cgns[ii]>=starting_node[rank])&&(vnodes_cgns[ii]<ending_node[rank])) {
                  
                  /*--- Build adjacency assuming the VTK connectivity ---*/
                  
                  if (ii < 3) {
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+1)%3]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+2)%3]);
                  } else {
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii-2)%3+3]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii-1)%3+3]);
                  }
                  adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+3)%6]);
                  
                }
              }
              
              Global_to_Local_Elem[global_id]=ielem;
              elem[ielem] = new CPrism(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], vnodes_cgns[4], vnodes_cgns[5]);
              ielem++; nelem_prism++;
              break;
              
            case PYRAMID:
              
              for ( unsigned short j = 0; j < N_POINTS_PYRAMID; j++ ) {
                vnodes_cgns[j] = connElems[k][s][j+1][i];
              }
              global_id = connElems[k][s][N_POINTS_PYRAMID+1][i];
              
              for (unsigned short ii=0; ii<N_POINTS_PYRAMID; ii++) {
                if ((vnodes_cgns[ii]>=starting_node[rank])&&(vnodes_cgns[ii]<ending_node[rank])) {
                  
                  /*--- Build adjacency assuming the VTK connectivity ---*/
                  
                  if (ii < 4) {
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+1)%4]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[(ii+3)%4]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[4]);
                  } else {
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[0]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[1]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[2]);
                    adj_nodes[vnodes_cgns[ii]-starting_node[rank]].push_back(vnodes_cgns[3]);
                  }
                  
                }
              }
              
              Global_to_Local_Elem[global_id]=ielem;
              elem[ielem] = new CPyramid(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], vnodes_cgns[4]);
              ielem++; nelem_pyramid++;
              break;
              
            default:
              if (rank == MASTER_NODE)
                cout << "Element type not suppported!" << endl;
#ifndef HAVE_MPI
              exit(EXIT_FAILURE);
#else
              MPI_Barrier(MPI_COMM_WORLD);
              MPI_Abort(MPI_COMM_WORLD,1);
              MPI_Finalize();
#endif
              break;
          }
        }
      }
    }
  }
  
#ifdef HAVE_MPI
  Local_nElemTri     = nelem_triangle;
  Local_nElemQuad    = nelem_quad;
  Local_nElemTet     = nelem_tetra;
  Local_nElemHex     = nelem_hexa;
  Local_nElemPrism   = nelem_prism;
  Local_nElemPyramid = nelem_pyramid;
  SU2_MPI::Allreduce(&Local_nElemTri,     &Global_nelem_triangle,  1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemQuad,    &Global_nelem_quad,      1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemTet,     &Global_nelem_tetra,     1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemHex,     &Global_nelem_hexa,      1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemPrism,   &Global_nelem_prism,     1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemPyramid, &Global_nelem_pyramid,   1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nelem_triangle = nelem_triangle;
  Global_nelem_quad     = nelem_quad;
  Global_nelem_tetra    = nelem_tetra;
  Global_nelem_hexa     = nelem_hexa;
  Global_nelem_prism    = nelem_prism;
  Global_nelem_pyramid  = nelem_pyramid;
#endif

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  
  /*--- Post process the adjacency information in order to get it into the
   proper format before sending the data to ParMETIS. We need to remove
   repeats and adjust the size of the array for each local node. ---*/
  
  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << "Building the graph adjacency structure." << endl;
  
  unsigned long loc_adjc_size=0;
  vector<unsigned long> adjac_vec;
  unsigned long adj_elem_size;
  vector<unsigned long>::iterator it;
  
  xadj = new idx_t[npoint_procs[rank]+1];
  xadj[0]=0;
  vector<unsigned long> temp_adjacency;
  unsigned long local_count=0;
  
  for (unsigned long i = 0; i < nPoint; i++) {
    
    for (unsigned long j=0; j<adj_nodes[i].size(); j++) {
      temp_adjacency.push_back(adj_nodes[i][j]);
    }
    
    sort(temp_adjacency.begin(), temp_adjacency.end());
    it = unique( temp_adjacency.begin(), temp_adjacency.end());
    loc_adjc_size=it - temp_adjacency.begin();
    
    temp_adjacency.resize( loc_adjc_size);
    xadj[local_count+1]=xadj[local_count]+loc_adjc_size;
    local_count++;
    
    for (unsigned long j=0; j<loc_adjc_size; j++) {
      adjac_vec.push_back(temp_adjacency[j]);
    }
    temp_adjacency.clear();
    adj_nodes[i].clear();
  }
  
  /*--- Now that we know the size, create the final adjacency array ---*/
  
  adj_elem_size = xadj[npoint_procs[rank]];
  adjacency = new idx_t [adj_elem_size];
  copy(adjac_vec.begin(), adjac_vec.end(), adjacency);
  
  xadj_size = npoint_procs[rank]+1;
  adjacency_size = adj_elem_size;
  
  /*--- Free temporary memory used to build the adjacency. ---*/
  
  adjac_vec.clear();
  
#endif
#endif
  
  adj_nodes.clear();
  
  /*--- Store the nodal coordinates from the linear partitioning. ---*/

  if ((rank == MASTER_NODE) && (size > SINGLE_NODE)) {
    cout << Global_nPoint << " grid points before linear partitioning." << endl;
  } else if (rank == MASTER_NODE) {
    cout << Global_nPoint << " grid points." << endl;
  }
  
  iPoint = 0;
  nPointNode = nPoint;
  node = new CPoint*[nPoint];
  GlobalIndex = starting_node[rank];
  for (int k = 0; k < nzones; k++ ) {
    for (unsigned long i = 0; i < nPoint; i++ ) {
      for (int j = 0; j < cell_dim; j++ ) Coord_cgns[j] = gridCoords[k][j][i];
      switch(nDim) {
        case 2:
          node[iPoint] = new CPoint(Coord_cgns[0], Coord_cgns[1], GlobalIndex, config);
          iPoint++; break;
        case 3:
          node[iPoint] = new CPoint(Coord_cgns[0], Coord_cgns[1], Coord_cgns[2], GlobalIndex, config);
          iPoint++; break;
      }
      GlobalIndex++;
    }
  }
  
  /*--- For now, the master node takes care of all markers. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Read number of markers ---*/

    nMarker = nMarkers;
    cout << nMarker << " surface markers." << endl;
    config->SetnMarker_All(nMarker);
    bound = new CPrimalGrid**[nMarker];
    nElem_Bound = new unsigned long [nMarker];
    Tag_to_Marker = new string [nMarker_Max];
    
    iMarker = 0;
    for ( int k = 0; k < nzones; k ++ ) {
      for ( int s = 0; s < nsections; s++ ) {
        if ( !isInternal[k][s] ) {
          
          /*--- Initialize some counter variables ---*/
          
          nelem_edge_bound = 0; nelem_triangle_bound = 0;
          nelem_quad_bound = 0; ielem = 0;
          
          Marker_Tag = sectionNames[k][s];
          
          /*--- Remove whitespaces from the marker names ---*/
          Marker_Tag.erase(remove(Marker_Tag.begin(), Marker_Tag.end(),' '), Marker_Tag.end());
          
          if (Marker_Tag != "SEND_RECEIVE") {
            nElem_Bound[iMarker] = nElems[k][s];
            if (rank == MASTER_NODE) {
              cout << nElem_Bound[iMarker]  << " boundary elements in index ";
              cout << iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;
            }
            bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];
            
            for ( int i = 0; i < nElems[k][s]; i++ ) {
              
              /*--- Get the VTK type for this element. Check for mixed
               elements. ---*/
              
              if (elemTypeVTK[k][s] == -1 ) {
                
                /*--- Mixed-element support. Check the elem type. ---*/
                
                ElementType_t elmt_type = ElementType_t(connElems[k][s][0][i]);
                cg_npe( elmt_type, &npe);
                
                switch (elmt_type) {
                  case NODE:    VTK_Type = 1;  break;
                  case BAR_2:   VTK_Type = 3;  break;
                  case BAR_3:   VTK_Type = 3;  break;
                  case TRI_3:   VTK_Type = 5;  break;
                  case QUAD_4:  VTK_Type = 9;  break;
                  case TETRA_4: VTK_Type = 10; break;
                  case HEXA_8:  VTK_Type = 12; break;
                  case PENTA_6: VTK_Type = 13; break;
                  case PYRA_5:  VTK_Type = 14; break;
                  default:
                    cout << "Kind of element not suppported!" << endl;
#ifndef HAVE_MPI
                    exit(EXIT_FAILURE);
#else
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Abort(MPI_COMM_WORLD,1);
                    MPI_Finalize();
#endif
                    break;
                }
                
                /*--- Transfer the nodes for this element. ---*/
                
                for ( int j = 1; j < npe+1; j++ ) {
                  vnodes_cgns[j-1] = connElems[k][s][j][i];
                }
                
              } else {
                
                /*--- Not a mixed section. We know the element type. ---*/
                
                VTK_Type = elemTypeVTK[k][s];
                
                /*--- Transfer the nodes for this element. ---*/
                
                for ( int j = 0; j < elemIndex[k][s]; j++ ) {
                  vnodes_cgns[j] = connElems[k][s][j][i];
                }
                
              }
              
              /*--- Instantiate the boundary elements. ---*/
              
              switch(VTK_Type) {
                case LINE:
                  if (nDim == 3) {
                    if (rank == MASTER_NODE) {
                      printf( "\n\n   !!! Error !!!\n" );
                      printf( "Remove line boundary elems from the mesh.");
                    }
#ifndef HAVE_MPI
                    exit(EXIT_FAILURE);
#else
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Abort(MPI_COMM_WORLD,1);
                    MPI_Finalize();
#endif
                  }
                  bound[iMarker][ielem] = new CLine(vnodes_cgns[0], vnodes_cgns[1],2);
                  ielem++; nelem_edge_bound++; break;
                case TRIANGLE:
                  bound[iMarker][ielem] = new CTriangle(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2],3);
                  ielem++; nelem_triangle_bound++; break;
                case QUADRILATERAL:
                  bound[iMarker][ielem] = new CQuadrilateral(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3],3);
                  ielem++; nelem_quad_bound++; break;
              }
            }
            
            /*--- Update config information storing the boundary information in the right place ---*/
            
            Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
            config->SetMarker_All_TagBound(iMarker, Marker_Tag);
            config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
            config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
            config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
            config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
            config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
            config->SetMarker_All_Analyze(iMarker, config->GetMarker_CfgFile_Analyze(Marker_Tag));
            config->SetMarker_All_FSIinterface(iMarker, config->GetMarker_CfgFile_FSIinterface(Marker_Tag));
            config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
            config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
            config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
            config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
            config->SetMarker_All_SendRecv(iMarker, NONE);
            
          }
          iMarker++;
        }
      }
    }
    
    /*--- Periodic transormations is not implement, store default zeros ---*/
    unsigned short nPeriodic = 1, iPeriodic = 0;
    config->SetnPeriodicIndex(nPeriodic);
    su2double* center    = new su2double[3];
    su2double* rotation  = new su2double[3];
    su2double* translate = new su2double[3];
    for (unsigned short iDim = 0; iDim < 3; iDim++) {
      center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
    }
    config->SetPeriodicCenter(iPeriodic, center);
    config->SetPeriodicRotation(iPeriodic, rotation);
    config->SetPeriodicTranslate(iPeriodic, translate);
    delete [] center; delete [] rotation; delete [] translate;
    
  }
  
  /*--- Deallocate temporary memory. ---*/
  
  delete[] vertices;
  delete[] cells;
  delete[] boundVerts;
  
  for ( int kk = 0; kk < nzones; kk++) {
    for (int ii = 0; ii < nsections; ii++) {
      if (isInternal[kk][ii]) {
        for (int jj = 0; jj < connSize; jj++) {
          if (connElems[kk][ii][jj] != NULL) delete [] connElems[kk][ii][jj];
        }
        if (connElems[kk][ii] != NULL) delete []  connElems[kk][ii];
      } else if (!isInternal[kk][ii] && rank == MASTER_NODE) {
        for (int jj = 0; jj < elemIndex[kk][ii]; jj++) {
          if (connElems[kk][ii][jj] != NULL) delete [] connElems[kk][ii][jj];
        }
        if (connElems[kk][ii] != NULL) delete [] connElems[kk][ii];
      }
    }
    if (connElems[kk] != NULL) delete [] connElems[kk];
  }
  if (connElems != NULL) delete[] connElems;
  
  for ( int j = 0; j < nzones; j++) {
    delete [] coordArray[j];
    delete [] elemTypeVTK[j];
    delete [] elemIndex[j];
    delete [] nElems[j];
    delete [] dataSize[j];
    delete [] isInternal[j];
    delete [] elemBegin[j];
    delete [] elemEnd[j];
    for (int ii = 0; ii < nsections; ii++) {
      delete[] sectionNames[j][ii];
    }
    delete[] sectionNames[j];
  }
  
  delete [] coordArray;
  delete [] elemTypeVTK;
  delete [] elemIndex;
  delete [] nElems;
  delete [] dataSize;
  delete [] isInternal;
  delete [] sectionNames;
  delete [] elemBegin;
  delete [] elemEnd;
  
  for ( int j = 0; j < nzones; j++) {
    for ( int i = 0; i < ncoords; i++ ) {
      delete [] gridCoords[j][i];
    }
    delete [] gridCoords[j];
  }
  delete [] gridCoords;
  
  delete [] nPoint_Linear;
  delete [] nElem_Linear;
  
  delete [] elemB;
  delete [] elemE;
  
  delete [] cgsize;

#else
  cout << "SU2 built without CGNS support!!" << endl;
  cout << "To use CGNS, remove the -DNO_CGNS directive ";
  cout << "from the makefile and supply the correct path";
  cout << " to the CGNS library." << endl;
  exit(EXIT_FAILURE);
#endif
  
}

void CPhysicalGeometry::Check_IntElem_Orientation(CConfig *config) {
  
  unsigned long Point_1, Point_2, Point_3, Point_4, Point_5, Point_6,
  iElem, triangle_flip = 0, quad_flip = 0, tet_flip = 0, prism_flip = 0,
  hexa_flip = 0, pyram_flip = 0;
  su2double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
  *Coord_5, *Coord_6, a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0}, c[3] = {0.0,0.0,0.0}, n[3] = {0.0,0.0,0.0}, test;
  unsigned short iDim;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Loop over all the elements ---*/
  
  for (iElem = 0; iElem < nElem; iElem++) {
    
    /*--- 2D grid, triangle case ---*/
    
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
      test = a[0]*b[1]-b[0]*a[1];
      
      if (test < 0.0) {
    	  elem[iElem]->Change_Orientation();
    	  triangle_flip++;
      }
    }
    
    /*--- 2D grid, quadrilateral case ---*/
    
    if (elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
      test_1 = a[0]*b[1]-b[0]*a[1];
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
        b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]); }
      test_2 = a[0]*b[1]-b[0]*a[1];
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
        b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]); }
      test_3 = a[0]*b[1]-b[0]*a[1];
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]); }
      test_4 = a[0]*b[1]-b[0]*a[1];
      
      if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0)) {
        elem[iElem]->Change_Orientation();
        quad_flip++;
      }
    }
    
    /*--- 3D grid, tetrahedron case ---*/
    
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      if (test < 0.0) {
    	  elem[iElem]->Change_Orientation();
    	  tet_flip++;
      }
      
    }
    
    /*--- 3D grid, prism case ---*/
    
    if (elem[iElem]->GetVTK_Type() == PRISM) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      Point_5 = elem[iElem]->GetNode(4); Coord_5 = node[Point_5]->GetCoord();
      Point_6 = elem[iElem]->GetNode(5); Coord_6 = node[Point_6]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        c[iDim] = (Coord_4[iDim]-Coord_1[iDim])+
        (Coord_5[iDim]-Coord_2[iDim])+
        (Coord_6[iDim]-Coord_3[iDim]); }
      
      /*--- The normal vector should point to the interior of the element ---*/
      
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_5[iDim]-Coord_4[iDim]);
        b[iDim] = 0.5*(Coord_6[iDim]-Coord_4[iDim]);
        c[iDim] = (Coord_1[iDim]-Coord_4[iDim])+
        (Coord_2[iDim]-Coord_5[iDim])+
        (Coord_3[iDim]-Coord_6[iDim]); }
      
      /*--- The normal vector should point to the interior of the element ---*/
      
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0)) {
          elem[iElem]->Change_Orientation();
          prism_flip++;
      }
      
    }
    
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(5); Coord_4 = node[Point_4]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(7); Coord_4 = node[Point_4]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(1); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(2); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(3); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(6); Coord_4 = node[Point_4]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(3); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(0); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(1); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0) || (test_3 < 0.0)
          || (test_4 < 0.0)) {
    	  elem[iElem]->Change_Orientation();
      	  hexa_flip++;
      }
      
    }
    
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0)) {
          elem[iElem]->Change_Orientation();
      	  pyram_flip++;
      }
      
    }
    
  }
  
#ifdef HAVE_MPI
  unsigned long Mytriangle_flip  = triangle_flip;
  unsigned long Myquad_flip      = quad_flip;
  unsigned long Mytet_flip       = tet_flip;
  unsigned long Myprism_flip     = prism_flip;
  unsigned long Myhexa_flip      = hexa_flip;
  unsigned long Mypyram_flip     = pyram_flip;

  SU2_MPI::Allreduce(&Mytriangle_flip, &triangle_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Myquad_flip, &quad_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Mytet_flip, &tet_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Myprism_flip, &prism_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Myhexa_flip, &hexa_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Mypyram_flip, &pyram_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (rank == MASTER_NODE) {
  	if (triangle_flip > 0) cout << "There has been a re-orientation of the TRIANGLE volume elements." << endl;
  	if (quad_flip > 0) cout << "There has been a re-orientation of the QUADRILATERAL volume elements." << endl;
  	if (tet_flip > 0) cout << "There has been a re-orientation of the TETRAHEDRON volume elements." << endl;
  	if (prism_flip > 0) cout << "There has been a re-orientation of the PRISM volume elements." << endl;
  	if (hexa_flip > 0) cout << "There has been a re-orientation of the HEXAHEDRON volume elements." << endl;
  	if (pyram_flip > 0) cout << "There has been a re-orientation of the PYRAMID volume elements." << endl;
  }

}

void CPhysicalGeometry::Check_BoundElem_Orientation(CConfig *config) {
  
  unsigned long Point_1_Surface, Point_2_Surface, Point_3_Surface, Point_4_Surface,
  iElem_Domain, Point_Domain = 0, Point_Surface, iElem_Surface,
  line_flip = 0, triangle_flip = 0, quad_flip = 0;
  su2double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
  *Coord_5, a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0}, c[3] = {0.0,0.0,0.0}, n[3] = {0.0,0.0,0.0}, test;
  unsigned short iDim, iMarker, iNode_Domain, iNode_Surface;
  bool find;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {
      
      for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
        
        iElem_Domain = bound[iMarker][iElem_Surface]->GetDomainElement();
        for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
          Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
          find = false;
          for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
            Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
            if (Point_Surface == Point_Domain) {find = true; break;}
          }
          if (!find) break;
        }
        
        /*--- 2D grid, line case ---*/
        
        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == LINE) {

          Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
          Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
          Coord_3 = node[Point_Domain]->GetCoord();

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
            b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
          }
          test = a[0]*b[1]-b[0]*a[1];

          if (test < 0.0) {
            bound[iMarker][iElem_Surface]->Change_Orientation();
            node[Point_1_Surface]->SetFlip_Orientation();
            node[Point_2_Surface]->SetFlip_Orientation();
            line_flip++;
          }

        }
        
        /*--- 3D grid, triangle case ---*/
        
        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TRIANGLE) {

          Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
          Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
          Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
          Coord_4 = node[Point_Domain]->GetCoord();

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
            b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
            c[iDim] = Coord_4[iDim]-Coord_1[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];

          test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
          if (test < 0.0) {
            bound[iMarker][iElem_Surface]->Change_Orientation();
            node[Point_1_Surface]->SetFlip_Orientation();
            node[Point_2_Surface]->SetFlip_Orientation();
            node[Point_3_Surface]->SetFlip_Orientation();
            triangle_flip++;
          }

        }
        
        /*--- 3D grid, rectangle case ---*/
        
        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == QUADRILATERAL) {

          Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
          Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
          Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
          Point_4_Surface = bound[iMarker][iElem_Surface]->GetNode(3); Coord_4 = node[Point_4_Surface]->GetCoord();
          Coord_5 = node[Point_Domain]->GetCoord();

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
            b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
            c[iDim] = Coord_5[iDim]-Coord_1[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];
          test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
            b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]);
            c[iDim] = Coord_5[iDim]-Coord_2[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];
          test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
            b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]);
            c[iDim] = Coord_5[iDim]-Coord_3[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];
          test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
            b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]);
            c[iDim] = Coord_5[iDim]-Coord_4[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];
          test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

          if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0)) {
            bound[iMarker][iElem_Surface]->Change_Orientation();
            node[Point_1_Surface]->SetFlip_Orientation();
            node[Point_2_Surface]->SetFlip_Orientation();
            node[Point_3_Surface]->SetFlip_Orientation();
            node[Point_4_Surface]->SetFlip_Orientation();
            quad_flip++;
          }

        }
      }
    }
  }

#ifdef HAVE_MPI
  unsigned long Myline_flip   = line_flip;
  unsigned long Mytriangle_flip  = triangle_flip;
  unsigned long Myquad_flip   = quad_flip;
  SU2_MPI::Allreduce(&Myline_flip, &line_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Mytriangle_flip, &triangle_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Myquad_flip, &quad_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (rank == MASTER_NODE) {
	if (line_flip > 0) cout << "There has been a re-orientation of the LINE surface elements." << endl;
  	if (triangle_flip > 0) cout << "There has been a re-orientation of the TRIANGLE surface elements." << endl;
  	if (quad_flip > 0) cout << "There has been a re-orientation of the QUADRILATERAL surface elements." << endl;
  }

}

void CPhysicalGeometry::ComputeWall_Distance(CConfig *config) {

  /*--- Compute the total number of nodes on no-slip boundaries ---*/

  unsigned long nVertex_SolidWall = 0;
  for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if( (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)  ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {
      nVertex_SolidWall += GetnVertex(iMarker);
    }
  }

  /*--- Allocate the vectors to hold boundary node coordinates
   and its local ID. ---*/

  vector<su2double>     Coord_bound(nDim*nVertex_SolidWall);
  vector<unsigned long> PointIDs(nVertex_SolidWall);

  /*--- Retrieve and store the coordinates of the no-slip boundary nodes
   and their local point IDs. ---*/

  unsigned long ii = 0, jj = 0;
  for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if( (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)  ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {
      for(unsigned long iVertex=0; iVertex<GetnVertex(iMarker); ++iVertex) {
        unsigned long iPoint = vertex[iMarker][iVertex]->GetNode();
        PointIDs[jj++] = iPoint;
        for(unsigned short iDim=0; iDim<nDim; ++iDim)
          Coord_bound[ii++] = node[iPoint]->GetCoord(iDim);
      }
    }
  }

  /*--- Build the ADT of the boundary nodes. ---*/

  su2_adtPointsOnlyClass WallADT(nDim, nVertex_SolidWall, Coord_bound.data(), PointIDs.data());

  /*--- Loop over all interior mesh nodes and compute the distances to each
   of the no-slip boundary nodes. Store the minimum distance to the wall
   for each interior mesh node. ---*/

  if( WallADT.IsEmpty() ) {
  
    /*--- No solid wall boundary nodes in the entire mesh.
     Set the wall distance to zero for all nodes. ---*/
    
    for(unsigned long iPoint=0; iPoint<GetnPoint(); ++iPoint)
      node[iPoint]->SetWall_Distance(0.0);
  }
  else {

    /*--- Solid wall boundary nodes are present. Compute the wall
     distance for all nodes. ---*/
    
    for(unsigned long iPoint=0; iPoint<GetnPoint(); ++iPoint) {
      
      su2double dist;
      unsigned long pointID;
      int rankID;
      
      WallADT.DetermineNearestNode(node[iPoint]->GetCoord(), dist,
                                   pointID, rankID);
      node[iPoint]->SetWall_Distance(dist);
    }
  }
  
}

void CPhysicalGeometry::SetPositive_ZArea(CConfig *config) {
  unsigned short iMarker, Boundary, Monitoring;
  unsigned long iVertex, iPoint;
  su2double *Normal, PositiveXArea, PositiveYArea, PositiveZArea, WettedArea;
  su2double TotalPositiveXArea = 0.0, TotalPositiveYArea = 0.0, TotalPositiveZArea = 0.0, TotalWettedArea = 0.0, AxiFactor;

  bool axisymmetric = config->GetAxisymmetric();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  PositiveXArea = 0.0;
  PositiveYArea = 0.0;
  PositiveZArea = 0.0;
  WettedArea = 0.0;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if (((Boundary == EULER_WALL)              ||
         (Boundary == HEAT_FLUX)               ||
         (Boundary == ISOTHERMAL)              ||
         (Boundary == LOAD_BOUNDARY)           ||
         (Boundary == DISPLACEMENT_BOUNDARY)) && (Monitoring == YES))

      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();

        if (node[iPoint]->GetDomain()) {
          Normal = vertex[iMarker][iVertex]->GetNormal();

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          if (nDim == 2) WettedArea += AxiFactor * sqrt (Normal[0]*Normal[0] + Normal[1]*Normal[1]);
          if (nDim == 3) WettedArea += sqrt (Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);

          if (Normal[0] < 0) PositiveXArea -= Normal[0];
          if (Normal[1] < 0) PositiveYArea -= Normal[1];
          if ((nDim == 3) && (Normal[2] < 0)) PositiveZArea -= Normal[2];

        }
      }
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&PositiveXArea, &TotalPositiveXArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&PositiveYArea, &TotalPositiveYArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&PositiveZArea, &TotalPositiveZArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&WettedArea, &TotalWettedArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  TotalPositiveXArea = PositiveXArea;
  TotalPositiveYArea = PositiveYArea;
  TotalPositiveZArea = PositiveZArea;
  TotalWettedArea    = WettedArea;
#endif
    
  if (config->GetRefAreaCoeff() == 0.0) {
  	if (nDim == 3) config->SetRefAreaCoeff(TotalPositiveZArea);
  	else config->SetRefAreaCoeff(TotalPositiveYArea);
  }
  
  if (rank == MASTER_NODE) {

  	cout << "Wetted area = "<< TotalWettedArea;
    if ((nDim == 3) || (axisymmetric)) { if (config->GetSystemMeasurements() == SI) cout <<" m^2." << endl; else cout <<" ft^2." << endl; }
    else { if (config->GetSystemMeasurements() == SI) cout <<" m." << endl; else cout <<" ft." << endl; }

    cout << "Area projection in the x-plane = "<< TotalPositiveXArea;
    if (nDim == 3) { if (config->GetSystemMeasurements() == SI) cout <<" m^2,"; else cout <<" ft^2,"; }
    else { if (config->GetSystemMeasurements() == SI) cout <<" m,"; else cout <<" ft,"; }

    cout << " y-plane = "<< TotalPositiveYArea;
    if (nDim == 3) { if (config->GetSystemMeasurements() == SI) cout <<" m^2,"; else cout <<" ft^2,"; }
    else { if (config->GetSystemMeasurements() == SI) cout <<" m." << endl; else cout <<" ft." << endl; }

    if (nDim == 3) { cout << " z-plane = "<< TotalPositiveZArea;
    if (config->GetSystemMeasurements() == SI) cout <<" m^2." << endl; else cout <<" ft^2."<< endl; }

  }
  
}

void CPhysicalGeometry::SetPoint_Connectivity(void) {
  
  unsigned short Node_Neighbor, iNode, iNeighbor;
  unsigned long jElem, Point_Neighbor, iPoint, iElem;
  
  /*--- Loop over all the elements ---*/
  
  for (iElem = 0; iElem < nElem; iElem++)
    
  /*--- Loop over all the nodes of an element ---*/
    
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      
      /*--- Store the element into the point ---*/
      
      node[iPoint]->SetElem(iElem);
    }

  /*--- Loop over all the points ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    
  /*--- Loop over all elements shared by the point ---*/
    
    for (iElem = 0; iElem < node[iPoint]->GetnElem(); iElem++) {
      
      jElem = node[iPoint]->GetElem(iElem);
      
      /*--- If we find the point iPoint in the surronding element ---*/
      
      for (iNode = 0; iNode < elem[jElem]->GetnNodes(); iNode++)
        
        if (elem[jElem]->GetNode(iNode) == iPoint)
          
        /*--- Localize the local index of the neighbor of iPoint in the element ---*/
          
          for (iNeighbor = 0; iNeighbor < elem[jElem]->GetnNeighbor_Nodes(iNode); iNeighbor++) {
            Node_Neighbor = elem[jElem]->GetNeighbor_Nodes(iNode, iNeighbor);
            Point_Neighbor = elem[jElem]->GetNode(Node_Neighbor);
            
            /*--- Store the point into the point ---*/
            
            node[iPoint]->SetPoint(Point_Neighbor);
          }
    }
  
  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetnNeighbor(node[iPoint]->GetnPoint());
  
}

void CPhysicalGeometry::SetRCM_Ordering(CConfig *config) {
  unsigned long iPoint, AdjPoint, AuxPoint, AddPoint, iElem, iNode, jNode;
  vector<unsigned long> Queue, AuxQueue, Result;
  unsigned short Degree, MinDegree, iDim, iMarker;
  bool *inQueue;
  
  inQueue = new bool [nPoint];
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    inQueue[iPoint] = false;
  
  /*--- Select the node with the lowest degree in the grid. ---*/
  
  MinDegree = node[0]->GetnNeighbor(); AddPoint = 0;
  for (iPoint = 1; iPoint < nPointDomain; iPoint++) {
    Degree = node[iPoint]->GetnPoint();
    if (Degree < MinDegree) { MinDegree = Degree; AddPoint = iPoint; }
  }
  
  /*--- Add the node in the first free position. ---*/
  
  Result.push_back(AddPoint); inQueue[AddPoint] = true;
  
  /*--- Loop until reorganize all the nodes ---*/
  
  do {
    
    /*--- Add to the queue all the nodes adjacent in the increasing
     order of their degree, checking if the element is already
     in the Queue. ---*/
    
    AuxQueue.clear();
    for (iNode = 0; iNode < node[AddPoint]->GetnPoint(); iNode++) {
      AdjPoint = node[AddPoint]->GetPoint(iNode);
      if ((!inQueue[AdjPoint]) && (AdjPoint < nPointDomain)) {
        AuxQueue.push_back(AdjPoint);
      }
    }
    
    if (AuxQueue.size() != 0) {
      
      /*--- Sort the auxiliar queue based on the number of neighbors ---*/
      
      for (iNode = 0; iNode < AuxQueue.size(); iNode++) {
        for (jNode = 0; jNode < AuxQueue.size() - 1 - iNode; jNode++) {
          if (node[AuxQueue[jNode]]->GetnPoint() > node[AuxQueue[jNode+1]]->GetnPoint()) {
            AuxPoint = AuxQueue[jNode];
            AuxQueue[jNode] = AuxQueue[jNode+1];
            AuxQueue[jNode+1] = AuxPoint;
          }
        }
      }
      
      Queue.insert(Queue.end(), AuxQueue.begin(), AuxQueue.end());
      for (iNode = 0; iNode < AuxQueue.size(); iNode++) {
        inQueue[AuxQueue[iNode]] = true;
      }
      
    }
    
    /*--- Extract the first node from the queue and add it in the first free
     position. ---*/
    
    if (Queue.size() != 0) {
      AddPoint = Queue[0];
      Result.push_back(Queue[0]);
      Queue.erase (Queue.begin(), Queue.begin()+1);
    }
    
    /*--- Add to the queue all the nodes adjacent in the increasing
     order of their degree, checking if the element is already
     in the Queue. ---*/
    
  } while (Queue.size() != 0);
  
  /*--- Check that all the points have been added ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (inQueue[iPoint] == false) Result.push_back(iPoint);
  }
  
  delete[] inQueue;
  
  reverse(Result.begin(), Result.end());
  
  /*--- Add the MPI points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    Result.push_back(iPoint);
  }
  
  /*--- Reset old data structures ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->ResetElem();
    node[iPoint]->ResetPoint();
    node[iPoint]->ResetBoundary();
    node[iPoint]->SetPhysicalBoundary(false);
    node[iPoint]->SetSolidBoundary(false);
    node[iPoint]->SetDomain(true);
  }
  
  /*--- Set the new coordinates ---*/
  
  su2double **AuxCoord;
  unsigned long *AuxGlobalIndex;
  
  AuxGlobalIndex = new unsigned long [nPoint];
  AuxCoord = new su2double* [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    AuxCoord[iPoint] = new su2double [nDim];
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    AuxGlobalIndex[iPoint] = node[iPoint]->GetGlobalIndex();
    for (iDim = 0; iDim < nDim; iDim++) {
      AuxCoord[iPoint][iDim] = node[iPoint]->GetCoord(iDim);
    }
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetGlobalIndex(AuxGlobalIndex[Result[iPoint]]);
    for (iDim = 0; iDim < nDim; iDim++)
      node[iPoint]->SetCoord(iDim, AuxCoord[Result[iPoint]][iDim]);
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    delete[] AuxCoord[iPoint];
  delete[] AuxCoord;
  delete[] AuxGlobalIndex;
  
  /*--- Set the new conectivities ---*/
  
  unsigned long *InvResult;
  InvResult = new unsigned long [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    InvResult[Result[iPoint]] = iPoint;
  
  for (iElem = 0; iElem < nElem; iElem++) {
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      elem[iElem]->SetNode(iNode, InvResult[iPoint]);
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
      
      string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
      if (Marker_Tag == "SEND_RECEIVE") {
        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          if (config->GetMarker_All_SendRecv(iMarker) < 0)
            node[bound[iMarker][iElem_Bound]->GetNode(0)]->SetDomain(false);
        }
      }
      
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        bound[iMarker][iElem]->SetNode(iNode, InvResult[iPoint]);
        node[InvResult[iPoint]]->SetBoundary(nMarker);
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[InvResult[iPoint]]->SetPhysicalBoundary(true);
        
        if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL ||
            config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX ||
            config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)
          node[InvResult[iPoint]]->SetSolidBoundary(true);
      }
    }
  }
  
  
  delete[] InvResult;
  
}

void CPhysicalGeometry::SetElement_Connectivity(void) {
  unsigned short first_elem_face, second_elem_face, iFace, iNode, jElem;
  unsigned long face_point, Test_Elem, iElem;
  
  /*--- Loop over all the elements, faces and nodes ---*/
  
  for (iElem = 0; iElem < nElem; iElem++)
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++)
      for (iNode = 0; iNode < elem[iElem]->GetnNodesFace(iFace); iNode++) {
        face_point = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iNode));
        
        /*--- Loop over all elements sharing the face point ---*/
        
        for (jElem = 0; jElem < node[face_point]->GetnElem(); jElem++) {
          Test_Elem = node[face_point]->GetElem(jElem);
          
          /*--- If it is a new element in this face ---*/
          
          if ((elem[iElem]->GetNeighbor_Elements(iFace) == -1) && (iElem < Test_Elem) &&
              (FindFace(iElem, Test_Elem, first_elem_face, second_elem_face))) {
            
            /*--- Localice which faces are sharing both elements ---*/
            
            elem[iElem]->SetNeighbor_Elements(Test_Elem, first_elem_face);
            
            /*--- Store the element for both elements ---*/
            
            elem[Test_Elem]->SetNeighbor_Elements(iElem, second_elem_face);
            
          }
        }
      }
}

void CPhysicalGeometry::SetBoundVolume(void) {
  unsigned short cont, iMarker, iElem, iNode_Domain, iNode_Surface;
  unsigned long Point_Domain, Point_Surface, Point, iElem_Surface, iElem_Domain;
  bool CheckVol;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      
      /*--- Choose and arbitrary point from the surface --*/
      Point = bound[iMarker][iElem_Surface]->GetNode(0);
      CheckVol = false;
      
      for (iElem = 0; iElem < node[Point]->GetnElem(); iElem++) {
        /*--- Look for elements surronding that point --*/
        cont = 0; iElem_Domain = node[Point]->GetElem(iElem);
        for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
          Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
          for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
            Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
            if (Point_Surface == Point_Domain) cont++;
            if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
          }
          if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
        }
        
        if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) {
          bound[iMarker][iElem_Surface]->SetDomainElement(iElem_Domain);
          CheckVol = true;
          break;
        }
      }
      if (!CheckVol) {
        cout << "The surface element ("<< iMarker <<", "<< iElem_Surface << ") doesn't have an associated volume element." << endl;
        exit(EXIT_FAILURE);
      }
    }
}

void CPhysicalGeometry::SetVertex(CConfig *config) {
  unsigned long  iPoint, iVertex, iElem;
  unsigned short iMarker, iNode;
  
  /*--- Initialize the Vertex vector for each node of the grid ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1, iMarker);
  
  /*--- Create and compute the vector with the number of vertex per marker ---*/
  
  nVertex = new unsigned long [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    /*--- Initialize the number of Bound Vertex for each Marker ---*/
    
    nVertex[iMarker] = 0;
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        
        /*--- Set the vertex in the node information ---*/
        
        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE)) {
          node[iPoint]->SetVertex(nVertex[iMarker], iMarker);
          nVertex[iMarker]++;
        }
      }
  }
  
  /*--- Initialize the Vertex vector for each node, the previous result is deleted ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1, iMarker);
  
  /*--- Create the bound vertex structure, note that the order
   is the same as in the input file, this is important for Send/Receive part ---*/
  
  vertex = new CVertex**[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [nVertex[iMarker]];
    nVertex[iMarker] = 0;
    
    /*--- Initialize the number of Bound Vertex for each Marker ---*/
    
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);

        /*--- Set the vertex in the node information ---*/
        
        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE)) {
          iVertex = nVertex[iMarker];
          vertex[iMarker][iVertex] = new CVertex(iPoint, nDim);
          
          if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
            vertex[iMarker][iVertex]->SetRotation_Type(bound[iMarker][iElem]->GetRotation_Type());
          }
          node[iPoint]->SetVertex(nVertex[iMarker], iMarker);
          nVertex[iMarker]++;
        }
      }
  }
}

void CPhysicalGeometry::SetCoord_CG(void) {
  unsigned short nNode, iDim, iMarker, iNode;
  unsigned long elem_poin, edge_poin, iElem, iEdge;
  su2double **Coord;
  
  /*--- Compute the center of gravity for elements ---*/
  
  for (iElem = 0; iElem<nElem; iElem++) {
    nNode = elem[iElem]->GetnNodes();
    Coord = new su2double* [nNode];
    
    /*--- Store the coordinates for all the element nodes ---*/
    
    for (iNode = 0; iNode < nNode; iNode++) {
      elem_poin = elem[iElem]->GetNode(iNode);
      Coord[iNode] = new su2double [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Coord[iNode][iDim]=node[elem_poin]->GetCoord(iDim);
    }
    
    /*--- Compute the element CG coordinates ---*/
    
    elem[iElem]->SetCoord_CG(Coord);
    
    for (iNode = 0; iNode < nNode; iNode++)
      if (Coord[iNode] != NULL) delete[] Coord[iNode];
    if (Coord != NULL) delete[] Coord;
  }
  
  /*--- Center of gravity for face elements ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
      nNode = bound[iMarker][iElem]->GetnNodes();
      Coord = new su2double* [nNode];
      
      /*--- Store the coordinates for all the element nodes ---*/
      
      for (iNode = 0; iNode < nNode; iNode++) {
        elem_poin = bound[iMarker][iElem]->GetNode(iNode);
        Coord[iNode] = new su2double [nDim];
        for (iDim = 0; iDim < nDim; iDim++)
          Coord[iNode][iDim]=node[elem_poin]->GetCoord(iDim);
      }
      /*--- Compute the element CG coordinates ---*/
      
      bound[iMarker][iElem]->SetCoord_CG(Coord);
      for (iNode = 0; iNode < nNode; iNode++)
        if (Coord[iNode] != NULL) delete[] Coord[iNode];
      if (Coord != NULL) delete[] Coord;
    }
  
  /*--- Center of gravity for edges ---*/
  
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    nNode = edge[iEdge]->GetnNodes();
    Coord = new su2double* [nNode];
    
    /*--- Store the coordinates for all the element nodes ---*/
    
    for (iNode = 0; iNode < nNode; iNode++) {
      edge_poin=edge[iEdge]->GetNode(iNode);
      Coord[iNode] = new su2double [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Coord[iNode][iDim]=node[edge_poin]->GetCoord(iDim);
    }
    
    /*--- Compute the edge CG coordinates ---*/
    
    edge[iEdge]->SetCoord_CG(Coord);
    
    for (iNode = 0; iNode < nNode; iNode++)
      if (Coord[iNode] != NULL) delete[] Coord[iNode];
    if (Coord != NULL) delete[] Coord;
  }
}

void CPhysicalGeometry::SetBoundControlVolume(CConfig *config, unsigned short action) {
  unsigned short Neighbor_Node, iMarker, iNode, iNeighbor_Nodes, iDim;
  unsigned long Neighbor_Point, iVertex, iPoint, iElem;
  long iEdge;
  su2double Area, *NormalFace = NULL;
  
  /*--- Update values of faces of the edge ---*/
  
  if (action != ALLOCATE)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        vertex[iMarker][iVertex]->SetZeroValues();
  
  su2double *Coord_Edge_CG = new su2double [nDim];
  su2double *Coord_Elem_CG = new su2double [nDim];
  su2double *Coord_Vertex = new su2double [nDim];
  
  /*--- Loop over all the markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    
  /*--- Loop over all the boundary elements ---*/
    
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      
    /*--- Loop over all the nodes of the boundary ---*/
      
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        iVertex = node[iPoint]->GetVertex(iMarker);
        
        /*--- Loop over the neighbor nodes, there is a face for each one ---*/
        
        for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
          Neighbor_Node = bound[iMarker][iElem]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
          Neighbor_Point = bound[iMarker][iElem]->GetNode(Neighbor_Node);
          
          /*--- Shared edge by the Neighbor Point and the point ---*/
          
          iEdge = FindEdge(iPoint, Neighbor_Point);
          for (iDim = 0; iDim < nDim; iDim++) {
            Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
            Coord_Elem_CG[iDim] = bound[iMarker][iElem]->GetCG(iDim);
            Coord_Vertex[iDim] = node[iPoint]->GetCoord(iDim);
          }
          switch (nDim) {
            case 2:
              
              /*--- Store the 2D face ---*/
              
              if (iNode == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Vertex);
              if (iNode == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Vertex, Coord_Elem_CG);
              break;
            case 3:
              
              /*--- Store the 3D face ---*/
              
              if (iNeighbor_Nodes == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, Coord_Vertex);
              if (iNeighbor_Nodes == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG, Coord_Vertex);
              break;
          }
        }
      }
  
  delete[] Coord_Edge_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_Vertex;
  
  /*--- Check if there is a normal with null area ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
      Area = sqrt(Area);
      if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }
  
}

void CPhysicalGeometry::MatchInterface(CConfig *config) {
  
  su2double epsilon = 1e-1;
  
  unsigned short nMarker_InterfaceBound = config->GetnMarker_InterfaceBound();
  
  if (nMarker_InterfaceBound != 0) {
    
    unsigned short iMarker, iDim, jMarker, pMarker = 0;
    unsigned long iVertex, iPoint, pVertex = 0, pPoint = 0, jVertex, jPoint, iPointGlobal, jPointGlobal, jVertex_, pPointGlobal = 0;
    su2double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local, maxdist_global;
    int iProcessor, pProcessor = 0;
    unsigned long nLocalVertex_Interface = 0, MaxLocalVertex_Interface = 0;
    int rank, nProcessor;
    
#ifndef HAVE_MPI
    rank = MASTER_NODE;
    nProcessor = SINGLE_NODE;
#else
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
    
    unsigned long *Buffer_Send_nVertex = new unsigned long [1];
    unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
    
    if (rank == MASTER_NODE) cout << "Set Interface boundary conditions (if any)." << endl;
    
    /*--- Compute the number of vertex that have interfase boundary condition
     without including the ghost nodes ---*/
    
    nLocalVertex_Interface = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) nLocalVertex_Interface ++;
        }
    
    Buffer_Send_nVertex[0] = nLocalVertex_Interface;
    
    /*--- Send Interface vertex information --*/
    
#ifndef HAVE_MPI
    MaxLocalVertex_Interface = nLocalVertex_Interface;
    Buffer_Receive_nVertex[0] = Buffer_Send_nVertex[0];
#else
    SU2_MPI::Allreduce(&nLocalVertex_Interface, &MaxLocalVertex_Interface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
    
    su2double *Buffer_Send_Coord = new su2double [MaxLocalVertex_Interface*nDim];
    unsigned long *Buffer_Send_Point = new unsigned long [MaxLocalVertex_Interface];
    unsigned long *Buffer_Send_GlobalIndex  = new unsigned long [MaxLocalVertex_Interface];
    unsigned long *Buffer_Send_Vertex  = new unsigned long [MaxLocalVertex_Interface];
    unsigned long *Buffer_Send_Marker  = new unsigned long [MaxLocalVertex_Interface];
    
    su2double *Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_Interface*nDim];
    unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_Interface];
    unsigned long *Buffer_Receive_GlobalIndex = new unsigned long [nProcessor*MaxLocalVertex_Interface];
    unsigned long *Buffer_Receive_Vertex = new unsigned long [nProcessor*MaxLocalVertex_Interface];
    unsigned long *Buffer_Receive_Marker = new unsigned long [nProcessor*MaxLocalVertex_Interface];
    
    unsigned long nBuffer_Coord = MaxLocalVertex_Interface*nDim;
    unsigned long nBuffer_Point = MaxLocalVertex_Interface;
    unsigned long nBuffer_GlobalIndex = MaxLocalVertex_Interface;
    unsigned long nBuffer_Vertex = MaxLocalVertex_Interface;
    unsigned long nBuffer_Marker = MaxLocalVertex_Interface;
    
    for (iVertex = 0; iVertex < MaxLocalVertex_Interface; iVertex++) {
      Buffer_Send_Point[iVertex] = 0;
      Buffer_Send_GlobalIndex[iVertex] = 0;
      Buffer_Send_Vertex[iVertex] = 0;
      Buffer_Send_Marker[iVertex] = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
    }
    
    /*--- Copy coordinates and point to the auxiliar vector --*/
    
    nLocalVertex_Interface = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          iPointGlobal = node[iPoint]->GetGlobalIndex();
          if (node[iPoint]->GetDomain()) {
            Buffer_Send_Point[nLocalVertex_Interface] = iPoint;
            Buffer_Send_GlobalIndex[nLocalVertex_Interface] = iPointGlobal;
            Buffer_Send_Vertex[nLocalVertex_Interface] = iVertex;
            Buffer_Send_Marker[nLocalVertex_Interface] = iMarker;
            for (iDim = 0; iDim < nDim; iDim++)
              Buffer_Send_Coord[nLocalVertex_Interface*nDim+iDim] = node[iPoint]->GetCoord(iDim);
            nLocalVertex_Interface++;
          }
        }
    
#ifndef HAVE_MPI
    for (unsigned long iBuffer_Coord = 0; iBuffer_Coord < nBuffer_Coord; iBuffer_Coord++)
      Buffer_Receive_Coord[iBuffer_Coord] = Buffer_Send_Coord[iBuffer_Coord];
    for (unsigned long iBuffer_Point = 0; iBuffer_Point < nBuffer_Point; iBuffer_Point++)
      Buffer_Receive_Point[iBuffer_Point] = Buffer_Send_Point[iBuffer_Point];
    for (unsigned long iBuffer_GlobalIndex = 0; iBuffer_GlobalIndex < nBuffer_GlobalIndex; iBuffer_GlobalIndex++)
      Buffer_Receive_GlobalIndex[iBuffer_GlobalIndex] = Buffer_Send_GlobalIndex[iBuffer_GlobalIndex];
    for (unsigned long iBuffer_Vertex = 0; iBuffer_Vertex < nBuffer_Vertex; iBuffer_Vertex++)
      Buffer_Receive_Vertex[iBuffer_Vertex] = Buffer_Send_Vertex[iBuffer_Vertex];
    for (unsigned long iBuffer_Marker = 0; iBuffer_Marker < nBuffer_Marker; iBuffer_Marker++)
      Buffer_Receive_Marker[iBuffer_Marker] = Buffer_Send_Marker[iBuffer_Marker];
#else
    SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, Buffer_Receive_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, Buffer_Receive_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
    
    
    /*--- Compute the closest point to a Near-Field boundary point ---*/
    
    maxdist_local = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
        
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          iPointGlobal = node[iPoint]->GetGlobalIndex();
          
          if (node[iPoint]->GetDomain()) {
            
            /*--- Coordinates of the boundary point ---*/
            
            Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;
            
            /*--- Loop over all the boundaries to find the pair ---*/
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
              for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_Interface+jVertex];
                jPointGlobal = Buffer_Receive_GlobalIndex[iProcessor*MaxLocalVertex_Interface+jVertex];
                jVertex_ = Buffer_Receive_Vertex[iProcessor*MaxLocalVertex_Interface+jVertex];
                jMarker = Buffer_Receive_Marker[iProcessor*MaxLocalVertex_Interface+jVertex];
                
                if (jPointGlobal != iPointGlobal) {
                  
                  /*--- Compute the distance ---*/
                  
                  dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
                    Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Interface+jVertex)*nDim+iDim];
                    dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
                  } dist = sqrt(dist);
                  
                  if (((dist < mindist) && (iProcessor != rank)) ||
                      ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
                    mindist = dist; pProcessor = iProcessor; pPoint = jPoint; pPointGlobal = jPointGlobal;
                    pVertex = jVertex_; pMarker = jMarker;
                    if (dist == 0.0) break;
                  }
                }
              }
            
            /*--- Store the value of the pair ---*/
            
            maxdist_local = max(maxdist_local, mindist);
            vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pPointGlobal, pVertex, pMarker, pProcessor);
            
            if (mindist > epsilon) {
              cout.precision(10);
              cout << endl;
              cout << "   Bad match for point " << iPoint << ".\tNearest";
              cout << " donor distance: " << scientific << mindist << ".";
              vertex[iMarker][iVertex]->SetDonorPoint(iPoint, iPointGlobal, pVertex, pMarker, pProcessor);
              maxdist_local = min(maxdist_local, 0.0);
            }
            
          }
        }
      }
    }
    
#ifndef HAVE_MPI
    maxdist_global = maxdist_local;
#else
    SU2_MPI::Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
#endif
    
    if (rank == MASTER_NODE) cout <<"The max distance between points is: " << maxdist_global <<"."<< endl;
    
    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Point;
    
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Point;
    
    delete[] Buffer_Send_nVertex;
    delete[] Buffer_Receive_nVertex;
    
  }
  
}

void CPhysicalGeometry::MatchNearField(CConfig *config) {
  
  su2double epsilon = 1e-1;
  
  unsigned short nMarker_NearFieldBound = config->GetnMarker_NearFieldBound();
  
  if (nMarker_NearFieldBound != 0) {
    
    unsigned short iMarker, iDim, jMarker, pMarker = 0;
    unsigned long iVertex, iPoint, pVertex = 0, pPoint = 0, jVertex, jPoint, iPointGlobal, jPointGlobal, jVertex_, pPointGlobal = 0;
    su2double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local, maxdist_global;
    int iProcessor, pProcessor = 0;
    unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
    int rank, nProcessor;
    
#ifndef HAVE_MPI
    rank = MASTER_NODE;
    nProcessor = SINGLE_NODE;
#else
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
    
    unsigned long *Buffer_Send_nVertex = new unsigned long [1];
    unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
    
    if (rank == MASTER_NODE) cout << "Set NearField boundary conditions (if any)." << endl;
    
    /*--- Compute the number of vertex that have interfase boundary condition
     without including the ghost nodes ---*/
    
    nLocalVertex_NearField = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) nLocalVertex_NearField ++;
        }
    
    Buffer_Send_nVertex[0] = nLocalVertex_NearField;
    
    /*--- Send NearField vertex information --*/
    
#ifndef HAVE_MPI
    MaxLocalVertex_NearField = nLocalVertex_NearField;
    Buffer_Receive_nVertex[0] = Buffer_Send_nVertex[0];
#else
    SU2_MPI::Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
    
    su2double *Buffer_Send_Coord = new su2double [MaxLocalVertex_NearField*nDim];
    unsigned long *Buffer_Send_Point = new unsigned long [MaxLocalVertex_NearField];
    unsigned long *Buffer_Send_GlobalIndex  = new unsigned long [MaxLocalVertex_NearField];
    unsigned long *Buffer_Send_Vertex  = new unsigned long [MaxLocalVertex_NearField];
    unsigned long *Buffer_Send_Marker  = new unsigned long [MaxLocalVertex_NearField];
    
    su2double *Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_NearField*nDim];
    unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_NearField];
    unsigned long *Buffer_Receive_GlobalIndex = new unsigned long [nProcessor*MaxLocalVertex_NearField];
    unsigned long *Buffer_Receive_Vertex = new unsigned long [nProcessor*MaxLocalVertex_NearField];
    unsigned long *Buffer_Receive_Marker = new unsigned long [nProcessor*MaxLocalVertex_NearField];
    
    unsigned long nBuffer_Coord = MaxLocalVertex_NearField*nDim;
    unsigned long nBuffer_Point = MaxLocalVertex_NearField;
    unsigned long nBuffer_GlobalIndex = MaxLocalVertex_NearField;
    unsigned long nBuffer_Vertex = MaxLocalVertex_NearField;
    unsigned long nBuffer_Marker = MaxLocalVertex_NearField;
    
    for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
      Buffer_Send_Point[iVertex] = 0;
      Buffer_Send_GlobalIndex[iVertex] = 0;
      Buffer_Send_Vertex[iVertex] = 0;
      Buffer_Send_Marker[iVertex] = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
    }
    
    /*--- Copy coordinates and point to the auxiliar vector --*/
    
    nLocalVertex_NearField = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          iPointGlobal = node[iPoint]->GetGlobalIndex();
          if (node[iPoint]->GetDomain()) {
            Buffer_Send_Point[nLocalVertex_NearField] = iPoint;
            Buffer_Send_GlobalIndex[nLocalVertex_NearField] = iPointGlobal;
            Buffer_Send_Vertex[nLocalVertex_NearField] = iVertex;
            Buffer_Send_Marker[nLocalVertex_NearField] = iMarker;
            for (iDim = 0; iDim < nDim; iDim++)
              Buffer_Send_Coord[nLocalVertex_NearField*nDim+iDim] = node[iPoint]->GetCoord(iDim);
            nLocalVertex_NearField++;
          }
        }
    
#ifndef HAVE_MPI
    for (unsigned long iBuffer_Coord = 0; iBuffer_Coord < nBuffer_Coord; iBuffer_Coord++)
      Buffer_Receive_Coord[iBuffer_Coord] = Buffer_Send_Coord[iBuffer_Coord];
    for (unsigned long iBuffer_Point = 0; iBuffer_Point < nBuffer_Point; iBuffer_Point++)
      Buffer_Receive_Point[iBuffer_Point] = Buffer_Send_Point[iBuffer_Point];
    for (unsigned long iBuffer_GlobalIndex = 0; iBuffer_GlobalIndex < nBuffer_GlobalIndex; iBuffer_GlobalIndex++)
      Buffer_Receive_GlobalIndex[iBuffer_GlobalIndex] = Buffer_Send_GlobalIndex[iBuffer_GlobalIndex];
    for (unsigned long iBuffer_Vertex = 0; iBuffer_Vertex < nBuffer_Vertex; iBuffer_Vertex++)
      Buffer_Receive_Vertex[iBuffer_Vertex] = Buffer_Send_Vertex[iBuffer_Vertex];
    for (unsigned long iBuffer_Marker = 0; iBuffer_Marker < nBuffer_Marker; iBuffer_Marker++)
      Buffer_Receive_Marker[iBuffer_Marker] = Buffer_Send_Marker[iBuffer_Marker];
#else
    SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, Buffer_Receive_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, Buffer_Receive_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
    
    
    /*--- Compute the closest point to a Near-Field boundary point ---*/
    
    maxdist_local = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          iPointGlobal = node[iPoint]->GetGlobalIndex();
          
          if (node[iPoint]->GetDomain()) {
            
            /*--- Coordinates of the boundary point ---*/
            
            Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;
            
            /*--- Loop over all the boundaries to find the pair ---*/
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
              for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_NearField+jVertex];
                jPointGlobal = Buffer_Receive_GlobalIndex[iProcessor*MaxLocalVertex_NearField+jVertex];
                jVertex_ = Buffer_Receive_Vertex[iProcessor*MaxLocalVertex_NearField+jVertex];
                jMarker = Buffer_Receive_Marker[iProcessor*MaxLocalVertex_NearField+jVertex];
                
                if (jPointGlobal != iPointGlobal) {
                  
                  /*--- Compute the distance ---*/
                  
                  dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
                    Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NearField+jVertex)*nDim+iDim];
                    dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
                  } dist = sqrt(dist);
                  
                  if (((dist < mindist) && (iProcessor != rank)) ||
                      ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
                    mindist = dist; pProcessor = iProcessor; pPoint = jPoint; pPointGlobal = jPointGlobal;
                    pVertex = jVertex_; pMarker = jMarker;
                    if (dist == 0.0) break;
                  }
                }
              }
            
            /*--- Store the value of the pair ---*/
            
            maxdist_local = max(maxdist_local, mindist);
            vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pPointGlobal, pVertex, pMarker, pProcessor);
            
            if (mindist > epsilon) {
              cout.precision(10);
              cout << endl;
              cout << "   Bad match for point " << iPoint << ".\tNearest";
              cout << " donor distance: " << scientific << mindist << ".";
              vertex[iMarker][iVertex]->SetDonorPoint(iPoint, iPointGlobal, pVertex, pMarker, pProcessor);
              maxdist_local = min(maxdist_local, 0.0);
            }
            
          }
        }
      }
    }
    
#ifndef HAVE_MPI
    maxdist_global = maxdist_local;
#else
    SU2_MPI::Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
#endif
    
    if (rank == MASTER_NODE) cout <<"The max distance between points is: " << maxdist_global <<"."<< endl;
    
    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Point;
    
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Point;
    
    delete[] Buffer_Send_nVertex;
    delete[] Buffer_Receive_nVertex;
    
  }
  
}

void CPhysicalGeometry::MatchActuator_Disk(CConfig *config) {
  
  su2double epsilon = 1e-1;
  
  unsigned short nMarker_ActDiskInlet = config->GetnMarker_ActDiskInlet();
  
  if (nMarker_ActDiskInlet != 0) {
    
    unsigned short iMarker, iDim;
    unsigned long iVertex, iPoint, iPointGlobal, pPoint = 0, pPointGlobal = 0, pVertex = 0, pMarker = 0, jVertex, jVertex_, jPoint, jPointGlobal, jMarker;
    su2double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local = 0.0, maxdist_global = 0.0;
    int iProcessor, pProcessor = 0;
    unsigned long nLocalVertex_ActDisk = 0, MaxLocalVertex_ActDisk = 0;
    int rank, nProcessor;
    unsigned short Beneficiary = 0, Donor = 0, iBC;
    bool Perimeter;
    
    for (iBC = 0; iBC < 2; iBC++) {
      
      if (iBC == 0) { Beneficiary = ACTDISK_INLET; Donor = ACTDISK_OUTLET; }
      if (iBC == 1) { Beneficiary = ACTDISK_OUTLET; Donor = ACTDISK_INLET; }
      
#ifndef HAVE_MPI
      rank = MASTER_NODE;
      nProcessor = SINGLE_NODE;
#else
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
      
      unsigned long *Buffer_Send_nVertex = new unsigned long [1];
      unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
      
      if ((iBC == 0) && (rank == MASTER_NODE)) cout << "Set Actuator Disk inlet boundary conditions." << endl;
      if ((iBC == 1) && (rank == MASTER_NODE)) cout << "Set Actuator Disk outlet boundary conditions." << endl;
      
      /*--- Compute the number of vertex that have an actuator disk outlet boundary condition
       without including the ghost nodes ---*/
      
      nLocalVertex_ActDisk = 0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Donor) {
          for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            if (node[iPoint]->GetDomain()) nLocalVertex_ActDisk ++;
          }
        }
      }
      
      Buffer_Send_nVertex[0] = nLocalVertex_ActDisk;
      
      /*--- Send actuator disk vertex information --*/
      
#ifndef HAVE_MPI
      MaxLocalVertex_ActDisk = nLocalVertex_ActDisk;
      Buffer_Receive_nVertex[0] = Buffer_Send_nVertex[0];
#else
      SU2_MPI::Allreduce(&nLocalVertex_ActDisk, &MaxLocalVertex_ActDisk, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
      
      /*--- Array dimensionalization --*/
      
      su2double *Buffer_Send_Coord = new su2double [MaxLocalVertex_ActDisk*nDim];
      unsigned long *Buffer_Send_Point  = new unsigned long [MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Send_GlobalIndex  = new unsigned long [MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Send_Vertex  = new unsigned long [MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Send_Marker  = new unsigned long [MaxLocalVertex_ActDisk];
      
      su2double *Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_ActDisk*nDim];
      unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Receive_GlobalIndex = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Receive_Vertex = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Receive_Marker = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];
      
      unsigned long nBuffer_Coord = MaxLocalVertex_ActDisk*nDim;
      unsigned long nBuffer_Point = MaxLocalVertex_ActDisk;
      unsigned long nBuffer_GlobalIndex = MaxLocalVertex_ActDisk;
      unsigned long nBuffer_Vertex = MaxLocalVertex_ActDisk;
      unsigned long nBuffer_Marker = MaxLocalVertex_ActDisk;
      
      for (iVertex = 0; iVertex < MaxLocalVertex_ActDisk; iVertex++) {
        Buffer_Send_Point[iVertex] = 0;
        Buffer_Send_GlobalIndex[iVertex] = 0;
        Buffer_Send_Vertex[iVertex] = 0;
        Buffer_Send_Marker[iVertex] = 0;
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
      }
      
      /*--- Copy coordinates and point to the auxiliar vector --*/
      
      nLocalVertex_ActDisk = 0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Donor) {
          for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            iPointGlobal = node[iPoint]->GetGlobalIndex();
            if (node[iPoint]->GetDomain()) {
              Buffer_Send_Point[nLocalVertex_ActDisk] = iPoint;
              Buffer_Send_GlobalIndex[nLocalVertex_ActDisk] = iPointGlobal;
              Buffer_Send_Vertex[nLocalVertex_ActDisk] = iVertex;
              Buffer_Send_Marker[nLocalVertex_ActDisk] = iMarker;
              for (iDim = 0; iDim < nDim; iDim++)
                Buffer_Send_Coord[nLocalVertex_ActDisk*nDim+iDim] = node[iPoint]->GetCoord(iDim);
              nLocalVertex_ActDisk++;
            }
          }
        }
      }
      
#ifndef HAVE_MPI
      for (unsigned long iBuffer_Coord = 0; iBuffer_Coord < nBuffer_Coord; iBuffer_Coord++)
        Buffer_Receive_Coord[iBuffer_Coord] = Buffer_Send_Coord[iBuffer_Coord];
      for (unsigned long iBuffer_Point = 0; iBuffer_Point < nBuffer_Point; iBuffer_Point++)
        Buffer_Receive_Point[iBuffer_Point] = Buffer_Send_Point[iBuffer_Point];
      for (unsigned long iBuffer_GlobalIndex = 0; iBuffer_GlobalIndex < nBuffer_GlobalIndex; iBuffer_GlobalIndex++)
        Buffer_Receive_GlobalIndex[iBuffer_GlobalIndex] = Buffer_Send_GlobalIndex[iBuffer_GlobalIndex];
      for (unsigned long iBuffer_Vertex = 0; iBuffer_Vertex < nBuffer_Vertex; iBuffer_Vertex++)
        Buffer_Receive_Vertex[iBuffer_Vertex] = Buffer_Send_Vertex[iBuffer_Vertex];
      for (unsigned long iBuffer_Marker = 0; iBuffer_Marker < nBuffer_Marker; iBuffer_Marker++)
        Buffer_Receive_Marker[iBuffer_Marker] = Buffer_Send_Marker[iBuffer_Marker];
      
#else
      SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, Buffer_Receive_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, Buffer_Receive_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
      
      /*--- Compute the closest point to an actuator disk inlet point ---*/
      
      maxdist_local = 0.0;
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Beneficiary) {
          
          for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            iPointGlobal = node[iPoint]->GetGlobalIndex();
            
            
            if (node[iPoint]->GetDomain()) {
              
              /*--- Coordinates of the boundary point ---*/
              
              Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;
              
              /*--- Loop over all the boundaries to find the pair ---*/
              
              Perimeter = false;
              
              for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                  jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_ActDisk+jVertex];
                  jPointGlobal = Buffer_Receive_GlobalIndex[iProcessor*MaxLocalVertex_ActDisk+jVertex];
                  jVertex_ = Buffer_Receive_Vertex[iProcessor*MaxLocalVertex_ActDisk+jVertex];
                  jMarker = Buffer_Receive_Marker[iProcessor*MaxLocalVertex_ActDisk+jVertex];
                  
 //                 if (jPointGlobal != iPointGlobal) {
 //                 ActDisk_Perimeter

                    /*--- Compute the distance ---*/
                    
                    dist = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                      Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_ActDisk+jVertex)*nDim+iDim];
                      dist += pow(Coord_j[iDim]-Coord_i[iDim], 2.0);
                    }
                    dist = sqrt(dist);
                    
                    if (dist < mindist) {
                      mindist = dist; pProcessor = iProcessor; pPoint = jPoint; pPointGlobal = jPointGlobal;
                      pVertex = jVertex_; pMarker = jMarker;
                      if (dist == 0.0) break;
                    }
                    
//                  }
//                  else { Perimeter = true; mindist = 0.0; dist = 0.0; break; }
               }
              }
              
              /*--- Store the value of the pair ---*/
              
              maxdist_local = max(maxdist_local, mindist);
              vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pPointGlobal, pVertex, pMarker, pProcessor);
              vertex[iMarker][iVertex]->SetActDisk_Perimeter(Perimeter);
              
              if (mindist > epsilon) {
                cout.precision(10);
                cout << endl;
                cout << "   Bad match for point " << iPoint << ".\tNearest";
                cout << " donor distance: " << scientific << mindist << ".";
                vertex[iMarker][iVertex]->SetDonorPoint(iPoint, iPointGlobal, pVertex, pMarker, pProcessor);
                maxdist_local = min(maxdist_local, 0.0);
              }
              
            }
          }
          
        }
      }
      
#ifndef HAVE_MPI
      maxdist_global = maxdist_local;
#else
      SU2_MPI::Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
#endif
      
      if (rank == MASTER_NODE) cout <<"The max distance between points is: " << maxdist_global <<"."<< endl;
      
      delete[] Buffer_Send_Coord;
      delete[] Buffer_Send_Point;
      
      delete[] Buffer_Receive_Coord;
      delete[] Buffer_Receive_Point;
      
      delete[] Buffer_Send_nVertex;
      delete[] Buffer_Receive_nVertex;
      
    }
  }
  
}

void CPhysicalGeometry::MatchZone(CConfig *config, CGeometry *geometry_donor, CConfig *config_donor,
                                  unsigned short val_iZone, unsigned short val_nZone) {
  
#ifndef HAVE_MPI
  
  unsigned short iMarker, jMarker;
  unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0, pGlobalPoint = 0;
  su2double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
  
//  if (val_iZone == ZONE_0) cout << "Set zone boundary conditions (if any)." << endl;
  
  maxdist = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iPoint = vertex[iMarker][iVertex]->GetNode();
      Coord_i = node[iPoint]->GetCoord();
      
      mindist = 1E6;
      for (jMarker = 0; jMarker < config_donor->GetnMarker_All(); jMarker++)
        for (jVertex = 0; jVertex < geometry_donor->GetnVertex(jMarker); jVertex++) {
          jPoint = geometry_donor->vertex[jMarker][jVertex]->GetNode();
          Coord_j = geometry_donor->node[jPoint]->GetCoord();
          if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
          if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
//          if (dist < mindist) { mindist = dist; pPoint = jPoint; pGlobalPoint = node[jPoint]->GetGlobalIndex();}
          if (dist < mindist) { mindist = dist; pPoint = jPoint; pGlobalPoint = geometry_donor->node[jPoint]->GetGlobalIndex();}
        }
      
      maxdist = max(maxdist, mindist);
      vertex[iMarker][iVertex]->SetDonorPoint(pPoint, MASTER_NODE, pGlobalPoint);
      
    }
  }
  
#else
  
  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint, jGlobalPoint = 0, pGlobalPoint = 0;
  su2double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist;
  int iProcessor, pProcessor = 0;
  unsigned long nLocalVertex_Zone = 0, nGlobalVertex_Zone = 0, MaxLocalVertex_Zone = 0;
  int rank, nProcessor;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  unsigned long *Buffer_Send_nVertex = new unsigned long [1];
  unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
//  if (val_iZone == ZONE_0 && rank == MASTER_NODE) cout << "Set zone boundary conditions (if any)." << endl;
  
  nLocalVertex_Zone = 0;
  for (iMarker = 0; iMarker < config_donor->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < geometry_donor->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry_donor->vertex[iMarker][iVertex]->GetNode();
      if (geometry_donor->node[iPoint]->GetDomain()) nLocalVertex_Zone ++;
    }
  
  Buffer_Send_nVertex[0] = nLocalVertex_Zone;
  
  /*--- Send Interface vertex information --*/
  
  SU2_MPI::Allreduce(&nLocalVertex_Zone, &nGlobalVertex_Zone, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalVertex_Zone, &MaxLocalVertex_Zone, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  su2double *Buffer_Send_Coord = new su2double [MaxLocalVertex_Zone*nDim];
  unsigned long *Buffer_Send_Point = new unsigned long [MaxLocalVertex_Zone];
  unsigned long *Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Zone];
  
  su2double *Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_Zone*nDim];
  unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_Zone];
  unsigned long *Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Zone];
  
  unsigned long nBuffer_Coord = MaxLocalVertex_Zone*nDim;
  unsigned long nBuffer_Point = MaxLocalVertex_Zone;
  

  for (iVertex = 0; iVertex < MaxLocalVertex_Zone; iVertex++) {
    Buffer_Send_Point[iVertex] = 0;
    Buffer_Send_GlobalPoint[iVertex] = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
  }
  
  /*--- Copy coordinates and point to the auxiliar vector --*/
  nLocalVertex_Zone = 0;
  for (iMarker = 0; iMarker < config_donor->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < geometry_donor->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry_donor->vertex[iMarker][iVertex]->GetNode();
      if (geometry_donor->node[iPoint]->GetDomain()) {
        Buffer_Send_Point[nLocalVertex_Zone] = iPoint;
        Buffer_Send_GlobalPoint[nLocalVertex_Zone] = geometry_donor->node[iPoint]->GetGlobalIndex();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord[nLocalVertex_Zone*nDim+iDim] = geometry_donor->node[iPoint]->GetCoord(iDim);
        nLocalVertex_Zone++;
      }
    }
  
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /*--- Compute the closest point to a Near-Field boundary point ---*/
  maxdist = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iPoint = vertex[iMarker][iVertex]->GetNode();
      
      if (node[iPoint]->GetDomain()) {
        
        /*--- Coordinates of the boundary point ---*/
        Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;
        
        /*--- Loop over all the boundaries to find the pair ---*/
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
          for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
            jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_Zone+jVertex];
            jGlobalPoint = Buffer_Receive_GlobalPoint[iProcessor*MaxLocalVertex_Zone+jVertex];

            /*--- Compute the distance ---*/
            dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
              Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Zone+jVertex)*nDim+iDim];
              dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
            } dist = sqrt(dist);
            
            if (((dist < mindist) && (iProcessor != rank)) ||
                ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
              mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
              pGlobalPoint = jGlobalPoint;
            }
          }
        
        /*--- Store the value of the pair ---*/
        maxdist = max(maxdist, mindist);
        vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor, pGlobalPoint);
        
        
      }
    }
  }
  
  delete[] Buffer_Send_Coord;
  delete[] Buffer_Send_Point;
  delete[] Buffer_Send_GlobalPoint;
  
  delete[] Buffer_Receive_Coord;
  delete[] Buffer_Receive_Point;
  delete[] Buffer_Receive_GlobalPoint;
  
  delete[] Buffer_Send_nVertex;
  delete[] Buffer_Receive_nVertex;
  
#endif
  
}


void CPhysicalGeometry::SetControlVolume(CConfig *config, unsigned short action) {
  unsigned long face_iPoint = 0, face_jPoint = 0, iPoint, iElem;
  long iEdge;
  unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
  su2double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint, *Coord_FacejPoint, Area,
  Volume, DomainVolume, my_DomainVolume, *NormalFace = NULL;
  bool change_face_orientation;
  int rank;
  
#ifndef HAVE_MPI
  rank = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Update values of faces of the edge ---*/
  if (action != ALLOCATE) {
    for (iEdge = 0; iEdge < (long)nEdge; iEdge++)
      edge[iEdge]->SetZeroValues();
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetVolume (0.0);
  }
  
  Coord_Edge_CG = new su2double [nDim];
  Coord_FaceElem_CG = new su2double [nDim];
  Coord_Elem_CG = new su2double [nDim];
  Coord_FaceiPoint = new su2double [nDim];
  Coord_FacejPoint = new su2double [nDim];
  
  my_DomainVolume = 0.0;
  for (iElem = 0; iElem < nElem; iElem++)
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {
      
      /*--- In 2D all the faces have only one edge ---*/
      if (nDim == 2) nEdgesFace = 1;
      /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
      if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);
      
      /*-- Loop over the edges of a face ---*/
      for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {
        
        /*--- In 2D only one edge (two points) per edge ---*/
        if (nDim == 2) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
          face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
        }
        
        /*--- In 3D there are several edges in each face ---*/
        if (nDim == 3) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace));
          if (iEdgesFace != nEdgesFace-1)
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace+1));
          else
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
        }
        
        /*--- We define a direction (from the smalest index to the greatest) --*/
        change_face_orientation = false;
        if (face_iPoint > face_jPoint) change_face_orientation = true;
        iEdge = FindEdge(face_iPoint, face_jPoint);
        
        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
          Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
          Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace, iDim);
          Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
          Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
        }
        
        switch (nDim) {
          case 2:
            /*--- Two dimensional problem ---*/
            if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG);
            else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG);
            Area = edge[iEdge]->GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_Elem_CG);
            node[face_iPoint]->AddVolume(Area); my_DomainVolume +=Area;
            Area = edge[iEdge]->GetVolume(Coord_FacejPoint, Coord_Edge_CG, Coord_Elem_CG);
            node[face_jPoint]->AddVolume(Area); my_DomainVolume +=Area;
            break;
          case 3:
            /*--- Three dimensional problem ---*/
            if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_FaceElem_CG, Coord_Edge_CG, Coord_Elem_CG);
            else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
            Volume = edge[iEdge]->GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
            node[face_iPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
            Volume = edge[iEdge]->GetVolume(Coord_FacejPoint, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
            node[face_jPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
            break;
        }
      }
    }
  
  /*--- Check if there is a normal with null area ---*/
  for (iEdge = 0; iEdge < (long)nEdge; iEdge++) {
    NormalFace = edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
    Area = sqrt(Area);
    if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
  }
  
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_DomainVolume, &DomainVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  DomainVolume = my_DomainVolume;
#endif
  
  if ((rank == MASTER_NODE) && (action == ALLOCATE)) {
    if (nDim == 2) cout <<"Area of the computational grid: "<< DomainVolume <<"."<< endl;
    if (nDim == 3) cout <<"Volume of the computational grid: "<< DomainVolume <<"."<< endl;
  }
  
  config->SetDomainVolume(DomainVolume);
  
  delete[] Coord_Edge_CG;
  delete[] Coord_FaceElem_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_FaceiPoint;
  delete[] Coord_FacejPoint;
}

void CPhysicalGeometry::VisualizeControlVolume(CConfig *config, unsigned short action) {
  
  /*--- This routine is only meant for visualization in serial currently ---*/
#ifndef HAVE_MPI
  
  unsigned long face_iPoint = 0, face_jPoint = 0, iElem, iPoint_Viz;
  long iEdge;
  unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
  su2double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint,
  *Coord_FacejPoint;
  int counter = 0;
  char cstr[MAX_STRING_SIZE], buffer[50];
  ofstream Tecplot_File;
  string mesh_filename;
  vector<su2double> X, Y, Z, X_n, Y_n, Z_n;
  su2double r1[3], r2[3], CrossProduct[3];
  
  /*--- Access the point number for control volume we want to vizualize ---*/
  
  iPoint_Viz = config->GetVisualize_CV();
  
  /*--- Allocate some structures for building the dual CVs ---*/
  
  Coord_Edge_CG     = new su2double [nDim];
  Coord_FaceElem_CG = new su2double [nDim];
  Coord_Elem_CG     = new su2double [nDim];
  Coord_FaceiPoint  = new su2double [nDim];
  Coord_FacejPoint  = new su2double [nDim];
  
  /*--- Loop over each face of each element ---*/
  
  CrossProduct[0] = 0.0; CrossProduct[1] = 0.0; CrossProduct[2] = 0.0;
  
  for (iElem = 0; iElem < nElem; iElem++) {
    
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {
      
      /*--- In 2D all the faces have only one edge ---*/
      if (nDim == 2) nEdgesFace = 1;
      /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
      if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);
      
      /*-- Loop over the edges of a face ---*/
      for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {
        
        /*--- In 2D only one edge (two points) per edge ---*/
        if (nDim == 2) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
          face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
        }
        
        /*--- In 3D there are several edges in each face ---*/
        if (nDim == 3) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace));
          if (iEdgesFace != nEdgesFace-1)
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace+1));
          else
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
        }
        
        /*--- We define a direction (from the smallest index to the greatest) --*/
        iEdge = FindEdge(face_iPoint, face_jPoint);
        
        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
          Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
          Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace, iDim);
          Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
          Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
        }
        
        /*--- Print out the coordinates for a set of triangles making
         up a single dual control volume for visualization. ---*/

        if (face_iPoint == iPoint_Viz || face_jPoint == iPoint_Viz) {

          if (nDim == 2) {
            X.push_back(Coord_Elem_CG[0]); X.push_back(Coord_Edge_CG[0]);
            Y.push_back(Coord_Elem_CG[1]); Y.push_back(Coord_Edge_CG[1]);
          } else if (nDim == 3) {
            X.push_back(Coord_FaceElem_CG[0]); X.push_back(Coord_Edge_CG[0]); X.push_back(Coord_Elem_CG[0]);
            Y.push_back(Coord_FaceElem_CG[1]); Y.push_back(Coord_Edge_CG[1]); Y.push_back(Coord_Elem_CG[1]);
            Z.push_back(Coord_FaceElem_CG[2]); Z.push_back(Coord_Edge_CG[2]); Z.push_back(Coord_Elem_CG[2]);

            for (iDim = 0; iDim < nDim; iDim++) {
              r1[iDim] = Coord_FaceElem_CG[iDim]-Coord_Elem_CG[iDim];
              r2[iDim] = Coord_Edge_CG[iDim]-Coord_Elem_CG[iDim];
            }
            CrossProduct[0] += 0.5*(r1[1]*r2[2] - r1[2]*r2[1]);
            CrossProduct[1] += 0.5*(r1[2]*r2[0] - r1[0]*r2[2]);
            CrossProduct[2] += 0.5*(r1[0]*r2[1] - r1[1]*r2[0]);
          }
          counter++;
        }
      }
    }
  }
  
  /*--- Write a Tecplot file to visualize the CV ---*/
  
  strcpy(cstr,"dual_cv");
  SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iPoint_Viz));
  strcat(cstr, buffer);
  
  Tecplot_File.open(cstr, ios::out);
  Tecplot_File << "TITLE= \"Visualization of the control volume\"" << endl;
  
  if (nDim == 2) {
    Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    Tecplot_File << "ZONE NODES= "<< counter*2 <<", ELEMENTS= ";
    Tecplot_File << counter <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
  } if (nDim == 3) {
    Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
    Tecplot_File << "ZONE NODES= "<< counter*3 <<", ELEMENTS= ";
    Tecplot_File << counter <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
  }
  
  /*--- Write coordinates for the nodes in the order that they were found
   for each of the edges/triangles making up a dual control volume. ---*/
  
  for (vector<su2double>::size_type i = 0; i != X.size(); i++) {
    Tecplot_File << X[i] << "\t" << Y[i];
    if (nDim == 3) Tecplot_File << "\t" << Z[i];
    Tecplot_File << "\n";
  }
  
  /*--- Create a new connectivity table in the order the faces were found ---*/
  
  int j;
  for (int i= 0; i < counter; i++) {
    if (nDim == 2) {
      j = i*2;
      Tecplot_File << j+1 <<"\t"<<j+2 <<"\t"<<j+2 <<"\t"<<j+2 << endl;
    } if (nDim == 3) {
      j = i*3;
      Tecplot_File << j+1 <<"\t"<<j+2 <<"\t"<<j+3 <<"\t"<<j+3 <<"\t";
      Tecplot_File << j+3<<"\t" <<j+3 <<"\t"<<j+3 <<"\t"<<j+3 << endl;
    }
  }
  
  Tecplot_File.close();
  X.clear();
  Y.clear();
  Z.clear();
  
  delete[] Coord_Edge_CG;
  delete[] Coord_FaceElem_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_FaceiPoint;
  delete[] Coord_FacejPoint;
  
#endif

}

void CPhysicalGeometry::SetMeshFile (CConfig *config, string val_mesh_out_filename) {
  unsigned long iElem, iPoint, iElem_Bound;
  unsigned short iMarker, iNodes, iDim;
  unsigned short iPeriodic, nPeriodic = 0;
  ofstream output_file;
  string Grid_Marker;
  char *cstr;
  su2double *center, *angles, *transl;
  
  cstr = new char [val_mesh_out_filename.size()+1];
  strcpy (cstr, val_mesh_out_filename.c_str());
  
  /*--- Open .su2 grid file ---*/
  
  output_file.precision(15);
  output_file.open(cstr, ios::out);
  
  /*--- Write dimension, number of elements and number of points ---*/
  
  output_file << "NDIME= " << nDim << endl;
  output_file << "NELEM= " << nElem << endl;
  for (iElem = 0; iElem < nElem; iElem++) {
    output_file << elem[iElem]->GetVTK_Type();
    for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
      output_file << "\t" << elem[iElem]->GetNode(iNodes);
    output_file << "\t"<<iElem<< endl;
  }
  
  /*--- Write the node coordinates ---*/
  
  output_file << "NPOIN= " << nPoint << "\t" << nPointDomain << endl;
  output_file.precision(15);
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      output_file << scientific << "\t" << node[iPoint]->GetCoord(iDim) ;
#ifndef HAVE_MPI
    output_file << "\t" << iPoint << endl;
#else
    output_file << "\t" << iPoint << "\t" << node[iPoint]->GetGlobalIndex() << endl;
#endif
    
  }
  
  /*--- Loop through and write the boundary info ---*/
  
  output_file << "NMARK= " << nMarker << endl;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    /*--- Ignore SEND_RECEIVE for the moment ---*/
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      
      Grid_Marker = config->GetMarker_All_TagBound(iMarker);
      output_file << "MARKER_TAG= " << Grid_Marker << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      
      if (nDim == 2) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
          output_file	<< iElem_Bound << endl;
        }
      }
      
      if (nDim == 3) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
          output_file	<< iElem_Bound << endl;
        }
      }
      
    } else if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
      output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
        bound[iMarker][iElem_Bound]->GetNode(0) << "\t" <<
        bound[iMarker][iElem_Bound]->GetRotation_Type() << endl;
      }
      
    }
  }
  
  /*--- Get the total number of periodic transformations ---*/
  
  nPeriodic = config->GetnPeriodicIndex();
  output_file << "NPERIODIC= " << nPeriodic << endl;
  
  /*--- From iPeriodic obtain the iMarker ---*/
  
  for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
    
    /*--- Retrieve the supplied periodic information. ---*/
    
    center = config->GetPeriodicCenter(iPeriodic);
    angles = config->GetPeriodicRotation(iPeriodic);
    transl = config->GetPeriodicTranslate(iPeriodic);
    
    output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
    output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
    output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
    output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
    
  }
  
  
  output_file.close();
}

void CPhysicalGeometry::SetCoord_Smoothing (unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) {
  unsigned short iSmooth, nneigh, iMarker;
  su2double *Coord_Old, *Coord_Sum, *Coord, *Coord_i, *Coord_j, Position_Plane = 0.0;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  su2double eps = 1E-6;
  bool NearField = false;
  
  Coord = new su2double [nDim];
  
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    su2double *Coord = node[iPoint]->GetCoord();
    node[iPoint]->SetCoord_Old(Coord);
  }
  
  /*--- Jacobi iterations ---*/
  for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetCoord_SumZero();
    
    
    /*--- Loop over Interior edges ---*/
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
      iPoint = edge[iEdge]->GetNode(0);
      Coord_i = node[iPoint]->GetCoord();
      
      jPoint = edge[iEdge]->GetNode(1);
      Coord_j = node[jPoint]->GetCoord();
      
      /*--- Accumulate nearest neighbor Coord to Res_sum for each variable ---*/
      node[iPoint]->AddCoord_Sum(Coord_j);
      node[jPoint]->AddCoord_Sum(Coord_i);
      
    }
    
    /*--- Loop over all mesh points (Update Coords with averaged sum) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      nneigh = node[iPoint]->GetnPoint();
      Coord_Sum = node[iPoint]->GetCoord_Sum();
      Coord_Old = node[iPoint]->GetCoord_Old();
      
      if (nDim == 2) {
        Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        if ((NearField) && ((Coord_Old[1] > Position_Plane-eps) && (Coord_Old[1] < Position_Plane+eps)))
          Coord[1] = Coord_Old[1];
      }
      
      if (nDim == 3) {
        Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        Coord[2] =(Coord_Old[2] + val_smooth_coeff*Coord_Sum[2]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        if ((NearField) && ((Coord_Old[2] > Position_Plane-eps) && (Coord_Old[2] < Position_Plane+eps)))
          Coord[2] = Coord_Old[2];
      }
      
      node[iPoint]->SetCoord(Coord);
    }
    
    /*--- Copy boundary values ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Coord_Old = node[iPoint]->GetCoord_Old();
        node[iPoint]->SetCoord(Coord_Old);
      }
  }
  
  delete[] Coord;
}

bool CPhysicalGeometry::FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem,
                                 unsigned short &face_second_elem) {
  
  /*--- Find repeated nodes between two elements to identify the common face ---*/
  unsigned long iPoint = 0, jPoint = 0;
  unsigned short face_node, iFace, iNode, jNode, nNodesFace;
  vector<unsigned long> CommonPoints, PointFaceFirst, PointFaceSecond;
  vector<unsigned long>::iterator IterPoint;
  pair<vector <unsigned long>::iterator, vector <unsigned long>::iterator> mypair;
  bool face_first_found = false, face_second_found =false;
  
  if (first_elem == second_elem) return false;
  
  for (iNode = 0; iNode < elem[first_elem]->GetnNodes(); iNode++) {
    iPoint = elem[first_elem]->GetNode(iNode);
    for (jNode = 0; jNode < elem[second_elem]->GetnNodes(); jNode++) {
      jPoint = elem[second_elem]->GetNode(jNode);
      if (iPoint == jPoint) {
        CommonPoints.push_back(iPoint);
        break;
      }
    }
  }
  
  /*--- Sort point in face and check that the list is unique ---*/
  sort( CommonPoints.begin(), CommonPoints.end());
  IterPoint = unique( CommonPoints.begin(), CommonPoints.end());
  CommonPoints.resize( distance(CommonPoints.begin(), IterPoint) );
  
  /*--- In 2D, the two elements must share two points that make up
   an edge, as all "faces" are edges in 2D. In 3D, we need to find
   exactly 3 (tri) or 4 (quad) common points. Return immediately to
   avoid a memory issue due to vectors of different lengths below. ---*/
  
  if ((nDim == 2) && (CommonPoints.size() != 2)) return false;
  if ((nDim == 3) && ((CommonPoints.size() != 3) &&
                      (CommonPoints.size() != 4))) return false;
  
  /*--- Search the sequence in the first element ---*/
  for (iFace = 0; iFace < elem[first_elem]->GetnFaces(); iFace++) {
    nNodesFace = elem[first_elem]->GetnNodesFace(iFace);
    
    if (nNodesFace == CommonPoints.size()) {
    for (iNode = 0; iNode < nNodesFace; iNode++) {
      face_node = elem[first_elem]->GetFaces(iFace, iNode);
      PointFaceFirst.push_back(elem[first_elem]->GetNode(face_node));
    }
    
    /*--- Sort face_poin to perform comparison ---*/
    sort( PointFaceFirst.begin(), PointFaceFirst.end());
    
    /*--- List comparison ---*/
    mypair = mismatch (PointFaceFirst.begin(), PointFaceFirst.end(), CommonPoints.begin());
    if (mypair.first == PointFaceFirst.end()) {
      face_first_elem = iFace;
      face_first_found = true;
      break;
    }
    
    PointFaceFirst.erase (PointFaceFirst.begin(), PointFaceFirst.end());
  }
  }
  
  /*--- Search the secuence in the second element ---*/
  for (iFace = 0; iFace < elem[second_elem]->GetnFaces(); iFace++) {
    nNodesFace = elem[second_elem]->GetnNodesFace(iFace);
    
    if (nNodesFace == CommonPoints.size()) {
    for (iNode = 0; iNode < nNodesFace; iNode++) {
      face_node = elem[second_elem]->GetFaces(iFace, iNode);
      PointFaceSecond.push_back(elem[second_elem]->GetNode(face_node));
    }
    
    /*--- Sort face_poin to perform comparison ---*/
    sort( PointFaceSecond.begin(), PointFaceSecond.end());
    
    /*--- List comparison ---*/
    mypair = mismatch (PointFaceSecond.begin(), PointFaceSecond.end(), CommonPoints.begin());
    if (mypair.first == PointFaceSecond.end()) {
      face_second_elem = iFace;
      face_second_found = true;
      break;
    }
    
    PointFaceSecond.erase (PointFaceSecond.begin(), PointFaceSecond.end());
  }
  }
  
  if (face_first_found && face_second_found) return true;
  else return false;
  
}

void CPhysicalGeometry::SetTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file) {
  
  unsigned long iElem, iPoint;
  unsigned short iDim;
  ofstream Tecplot_File;
  
  /*--- Open the tecplot file and write the header ---*/
  
  if (new_file) {
    Tecplot_File.open(mesh_filename, ios::out);
    Tecplot_File << "TITLE= \"Visualization of the volumetric grid\"" << endl;
    if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
  }
  else Tecplot_File.open(mesh_filename, ios::out | ios::app);
  
  Tecplot_File << "ZONE T= ";
  if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
  else Tecplot_File << "\"Deformed grid\", C=RED, ";
  Tecplot_File << "NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING= POINT";
  if (nDim == 2) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL"<< endl;
  if (nDim == 3) Tecplot_File << ", ZONETYPE= FEBRICK"<< endl;
  
  /*--- Adding coordinates ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
    Tecplot_File << "\n";
  }
  
  /*--- Adding conectivity ---*/
  
  for (iElem = 0; iElem < nElem; iElem++) {
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
      elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PRISM) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
    }
  }
  
  Tecplot_File.close();
}

void CPhysicalGeometry::SetBoundTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config) {
  
  ofstream Tecplot_File;
  unsigned long iPoint, Total_nElem_Bound, iElem, *PointSurface = NULL, nPointSurface = 0;
  unsigned short Coord_i, iMarker;
  
  /*--- It is important to do a renumbering to don't add points
   that do not belong to the surfaces ---*/
  
  PointSurface = new unsigned long[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    if (node[iPoint]->GetBoundary()) {
      PointSurface[iPoint] = nPointSurface;
      nPointSurface++;
    }
  
  /*--- Compute the total number of elements ---*/
  
  Total_nElem_Bound = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      Total_nElem_Bound += nElem_Bound[iMarker];
    }
  }
  
  /*--- Open the tecplot file and write the header ---*/
  
  if (new_file) {
    Tecplot_File.open(mesh_filename, ios::out);
    Tecplot_File << "TITLE= \"Visualization of the surface grid\"" << endl;
    if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
  }
  else Tecplot_File.open(mesh_filename, ios::out | ios::app);
  
  if (Total_nElem_Bound != 0) {
    
    /*--- Write the header of the file ---*/
    
    Tecplot_File << "ZONE T= ";
    if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
    else Tecplot_File << "\"Deformed grid\", C=RED, ";
    Tecplot_File << "NODES= "<< nPointSurface <<", ELEMENTS= "<< Total_nElem_Bound <<", DATAPACKING= POINT";
    if (nDim == 2) Tecplot_File << ", ZONETYPE= FELINESEG"<< endl;
    if (nDim == 3) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL"<< endl;
    
    /*--- Only write the coordiantes of the points that are on the surfaces ---*/
    
    if (nDim == 3) {
      for (iPoint = 0; iPoint < nPoint; iPoint++)
        if (node[iPoint]->GetBoundary()) {
          for (Coord_i = 0; Coord_i < nDim-1; Coord_i++)
            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
          Tecplot_File << node[iPoint]->GetCoord(nDim-1) << "\n";
        }
    }
    else {
      for (iPoint = 0; iPoint < nPoint; iPoint++)
        if (node[iPoint]->GetBoundary()) {
          for (Coord_i = 0; Coord_i < nDim; Coord_i++)
            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
          Tecplot_File << "\n";
        }
    }
    
    /*--- Write the cells using the new numbering ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
          if (nDim == 2) {
            Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
            << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << endl;
          }
          if (nDim == 3) {
            if (bound[iMarker][iElem]->GetnNodes() == 3) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << endl;
            }
            if (bound[iMarker][iElem]->GetnNodes() == 4) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(3)]+1 << endl;
            }
          }
        }
  }
  else {
    
    /*--- No elements in the surface ---*/
    
    if (nDim == 2) {
      Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
      Tecplot_File << "0.0 0.0"<< endl;
      Tecplot_File << "1 1"<< endl;
    }
    if (nDim == 3) {
      Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
      Tecplot_File << "0.0 0.0 0.0"<< endl;
      Tecplot_File << "1 1 1 1"<< endl;
    }
  }
  
  /*--- Dealocate memory and close the file ---*/
  
  delete[] PointSurface;
  Tecplot_File.close();
  
}

void CPhysicalGeometry::SetBoundSTL(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config) {
  
  ofstream STL_File;
  unsigned long this_node, iNode, nNode, iElem;
  unsigned short iDim, iMarker;
  su2double p[3] = {0.0,0.0,0.0}, u[3] = {0.0,0.0,0.0}, v[3] = {0.0,0.0,0.0}, n[3] = {0.0,0.0,0.0}, a;
  
  /*---	STL format:
   solid NAME
   ...
   facet normal 0.00 0.00 1.00
   outer loop
   vertex  2.00  2.00  0.00
   vertex -1.00  1.00  0.00
   vertex  0.00 -1.00  0.00
   endloop
   endfacet
   ...
   end solid
   ---*/
  
  /*--- Open the STL file ---*/
  
  if (new_file) STL_File.open(mesh_filename, ios::out);
  else STL_File.open(mesh_filename, ios::out | ios::app);
  
  /*--- Write the header of the file ---*/
  
  STL_File << "solid surface_mesh" << endl;
  
  /*--- Write facets of surface markers ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
        
        /*--- number of nodes for this elemnt ---*/
        
        nNode = bound[iMarker][iElem]->GetnNodes();
        
        /*--- Calculate Normal Vector ---*/
        
        for (iDim = 0; iDim < nDim; iDim++) {
          p[0] = node[bound[iMarker][iElem]->GetNode(0)]      ->GetCoord(iDim);
          p[1] = node[bound[iMarker][iElem]->GetNode(1)]      ->GetCoord(iDim);
          p[2] = node[bound[iMarker][iElem]->GetNode(nNode-1)]->GetCoord(iDim);
          u[iDim] = p[1]-p[0];
          v[iDim] = p[2]-p[0];
        }
        
        n[0] = u[1]*v[2]-u[2]*v[1];
        n[1] = u[2]*v[0]-u[0]*v[2];
        n[2] = u[0]*v[1]-u[1]*v[0];
        a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        
        /*--- Print normal vector ---*/
        
        STL_File << "  facet normal ";
        for (iDim = 0; iDim < nDim; iDim++) {
          STL_File << n[iDim]/a << " ";
        }
        STL_File << endl;
        
        /*--- STL Facet Loop --*/
        
        STL_File << "    outer loop" << endl;
        
        /*--- Print Nodes for Facet ---*/
        
        for (iNode = 0; iNode < nNode; iNode++) {
          this_node = bound[iMarker][iElem]->GetNode(iNode);
          STL_File << "      vertex ";
          for (iDim = 0; iDim < nDim; iDim++)
            STL_File << node[this_node]->GetCoord(iDim) << " ";
          if (nDim==2)
            STL_File << 0.0 << " ";
          STL_File <<  endl;
        }
        STL_File << "    endloop" << endl;
        STL_File << "  endfacet" << endl;
      }
  
  /*--- Done with Surface Mesh ---*/
  
  STL_File << "endsolid" << endl;
  
  /*--- Close the file ---*/
  
  STL_File.close();
  
}

void CPhysicalGeometry::SetColorGrid(CConfig *config) {
  
#ifdef HAVE_MPI
#ifdef HAVE_METIS
  
  unsigned long iPoint, iElem, iElem_Triangle, iElem_Tetrahedron, nElem_Triangle,
  nElem_Tetrahedron;
  idx_t ne = 0, nn, *elmnts = NULL, *epart = NULL, *npart = NULL, nparts, edgecut, *eptr;
  int rank, size;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if (size != SINGLE_ZONE)
    cout << endl <<"---------------------------- Grid partitioning --------------------------" << endl;
  
  unsigned short nDomain = size;
  
  nElem_Triangle = 0;
  nElem_Tetrahedron = 0;
  for (iElem = 0; iElem < GetnElem(); iElem++) {
    if (elem[iElem]->GetVTK_Type() == TRIANGLE)    nElem_Triangle = nElem_Triangle + 1;
    if (elem[iElem]->GetVTK_Type() == QUADRILATERAL)   nElem_Triangle = nElem_Triangle + 2;
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) nElem_Tetrahedron = nElem_Tetrahedron + 1;
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON)  nElem_Tetrahedron = nElem_Tetrahedron + 5;
    if (elem[iElem]->GetVTK_Type() == PYRAMID)     nElem_Tetrahedron = nElem_Tetrahedron + 2;
    if (elem[iElem]->GetVTK_Type() == PRISM)       nElem_Tetrahedron = nElem_Tetrahedron + 3;
  }
  
  if (GetnDim() == 2) {
    ne = nElem_Triangle;
    elmnts = new idx_t [ne*3];
  }
  if (GetnDim() == 3) {
    ne = nElem_Tetrahedron;
    elmnts = new idx_t [ne*4];
  }
  
  nn = nPoint;
  nparts = nDomain;
  epart = new idx_t [ne];
  npart = new idx_t [nn];
  eptr  = new idx_t[ne+1];
  
  /*--- Initialize the color vector ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetColor(0);
  
  if (nparts > 1) {
    
    iElem_Triangle = 0; iElem_Tetrahedron = 0;
    for (iElem = 0; iElem < GetnElem(); iElem++) {
      if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
        elmnts[3*iElem_Triangle+0]= elem[iElem]->GetNode(0);
        elmnts[3*iElem_Triangle+1]= elem[iElem]->GetNode(1);
        elmnts[3*iElem_Triangle+2]= elem[iElem]->GetNode(2);
        eptr[iElem_Triangle] = 3*iElem_Triangle;
        iElem_Triangle++;
      }
      if (elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
        elmnts[3*iElem_Triangle+0]= elem[iElem]->GetNode(0);
        elmnts[3*iElem_Triangle+1]= elem[iElem]->GetNode(1);
        elmnts[3*iElem_Triangle+2]= elem[iElem]->GetNode(2);
        eptr[iElem_Triangle] = 3*iElem_Triangle;
        iElem_Triangle++;
        elmnts[3*iElem_Triangle+0]= elem[iElem]->GetNode(0);
        elmnts[3*iElem_Triangle+1]= elem[iElem]->GetNode(2);
        elmnts[3*iElem_Triangle+2]= elem[iElem]->GetNode(3);
        eptr[iElem_Triangle] = 3*iElem_Triangle;
        iElem_Triangle++;
      }
      if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(1);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(3);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
      }
      if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(1);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(5);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(3);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(7);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(5);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(7);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(4);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(7);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(5);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(6);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(7);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(5);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
      }
      if (elem[iElem]->GetVTK_Type() == PYRAMID) {
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(1);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(4);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(3);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(4);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
      }
      if (elem[iElem]->GetVTK_Type() == PRISM) {
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(1);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(4);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(2);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(3);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(4);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(3);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(4);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(5);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(2);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
      }
    }
    
    /*--- Add final value to element pointer array ---*/
    
    if (GetnDim() == 2) eptr[ne] = 3*ne;
    else eptr[ne] = 4*ne;
    
    METIS_PartMeshNodal(&ne, &nn, eptr, elmnts, NULL, NULL, &nparts, NULL, NULL, &edgecut, epart, npart);
    
    cout << "Finished partitioning using METIS. ("  << edgecut << " edge cuts)." << endl;
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetColor(npart[iPoint]);
  }
  
  delete[] epart;
  delete[] npart;
  delete[] elmnts;
  delete[] eptr;
  
#endif
  
#endif
  
}

void CPhysicalGeometry::SetColorGrid_Parallel(CConfig *config) {
  
  /*--- Initialize the color vector ---*/
  
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetColor(0);
  
  /*--- This routine should only ever be called if we have parallel support
   with MPI and have the ParMETIS library compiled and linked. ---*/
  
#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  
  unsigned long iPoint;
  int rank, size;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  /*--- Only call ParMETIS if we have more than one rank to avoid errors ---*/
  
  if (size > SINGLE_NODE) {
    
    /*--- Create some structures that ParMETIS needs for partitioning. ---*/
    
    idx_t numflag, nparts, edgecut, wgtflag, ncon;
    
    idx_t *vtxdist = new idx_t[size+1];
    idx_t *part    = new idx_t[nPoint];
    
    real_t ubvec;
    real_t *tpwgts = new real_t[size];
    
    /*--- Some recommended defaults for the various ParMETIS options. ---*/
    
    wgtflag = 0;
    numflag = 0;
    ncon    = 1;
    ubvec   = 1.05;
    nparts  = (idx_t)size;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[1] = 0;
    
    /*--- Fill the necessary ParMETIS data arrays. Note that xadj_size and
     adjacency_size are class data members that have been defined and set
     earlier in the partitioning process. ---*/
    
    for (int i = 0; i < size; i++) {
      tpwgts[i] = 1.0/((real_t)size);
    }
    
    vtxdist[0] = 0;
    for (int i = 0; i < size; i++) {
      vtxdist[i+1] = (idx_t)ending_node[i];
    }
    
    /*--- Calling ParMETIS ---*/
    if (rank == MASTER_NODE) cout << "Calling ParMETIS..." << endl;
    ParMETIS_V3_PartKway(vtxdist,xadj, adjacency, NULL, NULL, &wgtflag,
                         &numflag, &ncon, &nparts, tpwgts, &ubvec, options,
                         &edgecut, part, &comm);
    if (rank == MASTER_NODE) {
      cout << "Finished partitioning using ParMETIS (";
      cout << edgecut << " edge cuts)." << endl;
    }
    
    /*--- Store the results of the partitioning (note that this is local
     since each processor is calling ParMETIS in parallel and storing the
     results for its initial piece of the grid. ---*/
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint]->SetColor(part[iPoint]);
    }
    
    /*--- Free all memory needed for the ParMETIS structures ---*/
    
    delete [] vtxdist;
    delete [] part;
    delete [] tpwgts;
    
  }
  
  /*--- Delete the memory from the geometry class that carried the
   adjacency structure. ---*/
  
  delete [] xadj;
  delete [] adjacency;
  
#endif
#endif
  
}

void CPhysicalGeometry::GetQualityStatistics(su2double *statistics) {
  unsigned long jPoint, Point_2, Point_3, iElem;
  su2double *Coord_j, *Coord_2, *Coord_3;
  unsigned short iDim;
  
  statistics[0] = 1e06;
  statistics[1] = 0;
  
  /*--- Loop interior edges ---*/
  for (iElem = 0; iElem < this->GetnElem(); iElem++) {
    
    if ((this->GetnDim() == 2) && (elem[iElem]->GetVTK_Type() == TRIANGLE)) {
      
      jPoint = elem[iElem]->GetNode(0); Coord_j = node[jPoint]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      
      /*--- Compute sides of the triangle ---*/
      su2double a = 0, b = 0, c = 0;
      for (iDim = 0; iDim < nDim; iDim++) {
        a += (Coord_2[iDim]-Coord_j[iDim])*(Coord_2[iDim]-Coord_j[iDim]);
        b += (Coord_3[iDim]-Coord_j[iDim])*(Coord_3[iDim]-Coord_j[iDim]);
        c += (Coord_3[iDim]-Coord_2[iDim])*(Coord_3[iDim]-Coord_2[iDim]);
      }
      a = sqrt(a); b = sqrt(b); c = sqrt(c);
      
      /*--- Compute semiperimeter (s) and area ---*/
      su2double s = 0.5*(a + b + c);
      su2double Area = sqrt(s*(s-a)*(s-b)*(s-c));
      
      /*--- Compute radius of the circumcircle (R) and of the incircle (r) ---*/
      su2double R = (a*b*c) / (4.0*Area);
      su2double r = Area / s;
      su2double roR = r / R;
      
      /*--- Update statistics ---*/
      if (roR < statistics[0])
        statistics[0] = roR;
      statistics[1] += roR;
      
    }
  }
  statistics[1] /= this->GetnElem();
  
}

void CPhysicalGeometry::SetRotationalVelocity(CConfig *config, unsigned short val_iZone) {
  
  unsigned long iPoint;
  su2double RotVel[3], Distance[3], *Coord, Center[3], Omega[3], L_Ref;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Center of rotation & angular velocity vector from config ---*/
  
  Center[0] = config->GetMotion_Origin_X(val_iZone);
  Center[1] = config->GetMotion_Origin_Y(val_iZone);
  Center[2] = config->GetMotion_Origin_Z(val_iZone);
  Omega[0]  = config->GetRotation_Rate_X(val_iZone)/config->GetOmega_Ref();
  Omega[1]  = config->GetRotation_Rate_Y(val_iZone)/config->GetOmega_Ref();
  Omega[2]  = config->GetRotation_Rate_Z(val_iZone)/config->GetOmega_Ref();
  L_Ref     = config->GetLength_Ref();
  
  /*--- Print some information to the console ---*/
  
  if (rank == MASTER_NODE) {
    cout << " Rotational origin (x, y, z): ( " << Center[0] << ", " << Center[1];
    cout << ", " << Center[2] << " )" << endl;
    cout << " Angular velocity about x, y, z axes: ( " << Omega[0] << ", ";
    cout << Omega[1] << ", " << Omega[2] << " ) rad/s" << endl;
  }
  
  /*--- Loop over all nodes and set the rotational velocity ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Get the coordinates of the current node ---*/
    
    Coord = node[iPoint]->GetCoord();
    
    /*--- Calculate the non-dim. distance from the rotation center ---*/
    
    Distance[0] = (Coord[0]-Center[0])/L_Ref;
    Distance[1] = (Coord[1]-Center[1])/L_Ref;
    Distance[2] = (Coord[2]-Center[2])/L_Ref;
    
    /*--- Calculate the angular velocity as omega X r ---*/
    
    RotVel[0] = Omega[1]*(Distance[2]) - Omega[2]*(Distance[1]);
    RotVel[1] = Omega[2]*(Distance[0]) - Omega[0]*(Distance[2]);
    RotVel[2] = Omega[0]*(Distance[1]) - Omega[1]*(Distance[0]);
    
    /*--- Store the grid velocity at this node ---*/
    
    node[iPoint]->SetGridVel(RotVel);
    
  }
  
}

void CPhysicalGeometry::SetTranslationalVelocity(CConfig *config) {
  
  unsigned short iDim;
  unsigned long iPoint;
  su2double xDot[3] = {0.0,0.0,0.0};
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Get the translational velocity vector from config ---*/
  
  xDot[0] = config->GetTranslation_Rate_X(ZONE_0)/config->GetVelocity_Ref();
  xDot[1] = config->GetTranslation_Rate_Y(ZONE_0)/config->GetVelocity_Ref();
  xDot[2] = config->GetTranslation_Rate_Z(ZONE_0)/config->GetVelocity_Ref();
  
  /*--- Print some information to the console ---*/
  
  if (rank == MASTER_NODE) {
    cout << " Non-dim. translational velocity: (" << xDot[0] << ", " << xDot[1];
    cout << ", " << xDot[2] << ")." << endl;
  }
  
  /*--- Loop over all nodes and set the translational velocity ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Store the grid velocity at this node ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      node[iPoint]->SetGridVel(iDim,xDot[iDim]);
    }
  
  }
  
}

void CPhysicalGeometry::SetGridVelocity(CConfig *config, unsigned long iter) {
  
  /*--- Local variables ---*/
  
  su2double *Coord_nP1 = NULL, *Coord_n = NULL, *Coord_nM1 = NULL;
  su2double TimeStep, GridVel = 0.0;
  unsigned long iPoint;
  unsigned short iDim;
  
  /*--- Compute the velocity of each node in the volume mesh ---*/
  
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/
    
    Coord_nM1 = node[iPoint]->GetCoord_n1();
    Coord_n   = node[iPoint]->GetCoord_n();
    Coord_nP1 = node[iPoint]->GetCoord();

    /*--- Unsteady time step ---*/
    
    TimeStep = config->GetDelta_UnstTimeND();
    
    /*--- Compute mesh velocity with 1st or 2nd-order approximation ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        GridVel = ( Coord_nP1[iDim] - Coord_n[iDim] ) / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        GridVel = ( 3.0*Coord_nP1[iDim] - 4.0*Coord_n[iDim]
                   + 1.0*Coord_nM1[iDim] ) / (2.0*TimeStep);
      
      /*--- Store grid velocity for this point ---*/
      
      node[iPoint]->SetGridVel(iDim, GridVel);
    }
  }
  
}

void CPhysicalGeometry::Set_MPI_Coord(CConfig *config) {
  
  unsigned short iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_Coord = NULL, *Buffer_Send_Coord = NULL, *Coord = NULL, *newCoord = NULL;
  su2double *translation;
  newCoord = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = nVertex[MarkerS];  nVertexR = nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nDim;        nBufferR_Vector = nVertexR*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive_Coord = new su2double [nBufferR_Vector];
      Buffer_Send_Coord = new su2double[nBufferS_Vector];
      
      /*--- Copy the coordinates that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = vertex[MarkerS][iVertex]->GetNode();
        Coord = node[iPoint]->GetCoord();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord[iDim*nVertexS+iVertex] = Coord[iDim];
      }
      
#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Coord, nBufferS_Vector, MPI_DOUBLE, send_to,0,
                   Buffer_Receive_Coord, nBufferR_Vector, MPI_DOUBLE, receive_from,0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Receive_Coord[iDim*nVertexR+iVertex] = Buffer_Send_Coord[iDim*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send_Coord;
      
      /*--- Do the coordinate transformation ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        
        iPoint = vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        translation = config->GetPeriodicTranslate(iPeriodic_Index);
        
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
        
        /*--- Copy coordinates before performing transformation. ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          newCoord[iDim] = Buffer_Receive_Coord[iDim*nVertexR+iVertex];
        
        /*--- Rotate the coordinates. ---*/
        
        if (nDim == 2) {
          newCoord[0] = (rotMatrix[0][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_Coord[1*nVertexR+iVertex]) - translation[0];
          newCoord[1] = (rotMatrix[1][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_Coord[1*nVertexR+iVertex]) - translation[1];
        }
        else {
          newCoord[0] = (rotMatrix[0][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_Coord[1*nVertexR+iVertex] +
                         rotMatrix[0][2]*Buffer_Receive_Coord[2*nVertexR+iVertex]);
          newCoord[1] = (rotMatrix[1][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_Coord[1*nVertexR+iVertex] +
                         rotMatrix[1][2]*Buffer_Receive_Coord[2*nVertexR+iVertex]);
          newCoord[2] = (rotMatrix[2][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[2][1]*Buffer_Receive_Coord[1*nVertexR+iVertex] +
                         rotMatrix[2][2]*Buffer_Receive_Coord[2*nVertexR+iVertex]);
        }
        
        /*--- Copy transformed coordinates back into buffer. ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          node[iPoint]->SetCoord(iDim, newCoord[iDim]);
        
      }
      
      /*--- Deallocate receive buffer. ---*/
      
      delete [] Buffer_Receive_Coord;
      
    }
    
  }
  
  delete [] newCoord;
  
}

void CPhysicalGeometry::Set_MPI_GridVel(CConfig *config) {
  
  unsigned short iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_GridVel = NULL, *Buffer_Send_GridVel = NULL, *GridVel = NULL, *newGridVel = NULL;
  
  newGridVel = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
     
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
   
      nVertexS = nVertex[MarkerS];  nVertexR = nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nDim;        nBufferR_Vector = nVertexR*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive_GridVel = new su2double [nBufferR_Vector];
      Buffer_Send_GridVel = new su2double[nBufferS_Vector];
      
      /*--- Copy the grid velocity that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = vertex[MarkerS][iVertex]->GetNode();
        GridVel = node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_GridVel[iDim*nVertexS+iVertex] = GridVel[iDim];
      }
      
#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_GridVel, nBufferS_Vector, MPI_DOUBLE, send_to,0,
                   Buffer_Receive_GridVel, nBufferR_Vector, MPI_DOUBLE, receive_from,0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Receive_GridVel[iDim*nVertexR+iVertex] = Buffer_Send_GridVel[iDim*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send_GridVel;
      
      /*--- Do the coordinate transformation ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        
        iPoint = vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = vertex[MarkerR][iVertex]->GetRotation_Type();
        
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
        
        /*--- Copy grid velocity before performing transformation. ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          newGridVel[iDim] = Buffer_Receive_GridVel[iDim*nVertexR+iVertex];
        
        if (nDim == 2) {
          newGridVel[0] = (rotMatrix[0][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[0][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex]);
          newGridVel[1] = (rotMatrix[1][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[1][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex]);
        }
        else {
          newGridVel[0] = (rotMatrix[0][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[0][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex] +
                           rotMatrix[0][2]*Buffer_Receive_GridVel[2*nVertexR+iVertex]);
          newGridVel[1] = (rotMatrix[1][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[1][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex] +
                           rotMatrix[1][2]*Buffer_Receive_GridVel[2*nVertexR+iVertex]);
          newGridVel[2] = (rotMatrix[2][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[2][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex] +
                           rotMatrix[2][2]*Buffer_Receive_GridVel[2*nVertexR+iVertex]);
        }
        
        /*--- Copy transformed grid velocity back into buffer. ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      
      delete [] Buffer_Receive_GridVel;
      
    }
    
  }
  
  delete [] newGridVel;
  
}

void CPhysicalGeometry::Set_MPI_OldCoord(CConfig *config) {

  unsigned short iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;

  su2double *Buffer_Receive_Coord_n = NULL, *Buffer_Send_Coord_n = NULL, *Coord_n = NULL, *newCoord_n = NULL;

  newCoord_n = new su2double[nDim];

#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = nVertex[MarkerS];  nVertexR = nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nDim;        nBufferR_Vector = nVertexR*nDim;

      /*--- Allocate Receive and send buffers  ---*/

      Buffer_Receive_Coord_n = new su2double [nBufferR_Vector];
      Buffer_Send_Coord_n = new su2double[nBufferS_Vector];

      /*--- Copy the coordinates that should be sended ---*/

      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = vertex[MarkerS][iVertex]->GetNode();
        Coord_n = node[iPoint]->GetCoord_n();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord_n[iDim*nVertexS+iVertex] = Coord_n[iDim];
      }

#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Coord_n, nBufferS_Vector, MPI_DOUBLE, send_to,0,
                   Buffer_Receive_Coord_n, nBufferR_Vector, MPI_DOUBLE, receive_from,0, MPI_COMM_WORLD, &status);
#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Receive_Coord_n[iDim*nVertexR+iVertex] = Buffer_Send_Coord_n[iDim*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/

      delete [] Buffer_Send_Coord_n;

      /*--- Do the coordinate transformation ---*/

      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/

        iPoint = vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = vertex[MarkerR][iVertex]->GetRotation_Type();

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

        /*--- Copy coordinates before performing transformation. ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          newCoord_n[iDim] = Buffer_Receive_Coord_n[iDim*nVertexR+iVertex];

        /*--- Rotate the coordinates. ---*/

        if (nDim == 2) {
          newCoord_n[0] = (rotMatrix[0][0]*Buffer_Receive_Coord_n[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_Coord_n[1*nVertexR+iVertex]);
          newCoord_n[1] = (rotMatrix[1][0]*Buffer_Receive_Coord_n[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_Coord_n[1*nVertexR+iVertex]);
        }
        else {
          newCoord_n[0] = (rotMatrix[0][0]*Buffer_Receive_Coord_n[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_Coord_n[1*nVertexR+iVertex] +
                         rotMatrix[0][2]*Buffer_Receive_Coord_n[2*nVertexR+iVertex]);
          newCoord_n[1] = (rotMatrix[1][0]*Buffer_Receive_Coord_n[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_Coord_n[1*nVertexR+iVertex] +
                         rotMatrix[1][2]*Buffer_Receive_Coord_n[2*nVertexR+iVertex]);
          newCoord_n[2] = (rotMatrix[2][0]*Buffer_Receive_Coord_n[0*nVertexR+iVertex] +
                         rotMatrix[2][1]*Buffer_Receive_Coord_n[1*nVertexR+iVertex] +
                         rotMatrix[2][2]*Buffer_Receive_Coord_n[2*nVertexR+iVertex]);
        }

        /*--- Copy transformed coordinates back into buffer. ---*/

        node[iPoint]->SetCoord_n(newCoord_n);

      }

      /*--- Deallocate receive buffer. ---*/

      delete [] Buffer_Receive_Coord_n;

    }

  }

  delete [] newCoord_n;

  /*--------------------------------------------------------------------------------------------------*/
  /*--- We repeat the process for the coordinate n-1, in the case that the simulation is 2nd order ---*/
  /*--------------------------------------------------------------------------------------------------*/

  if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) {

	  su2double *Buffer_Receive_Coord_n1 = NULL, *Buffer_Send_Coord_n1 = NULL, *Coord_n1 = NULL, *newCoord_n1 = NULL;
	  newCoord_n1 = new su2double[nDim];

	  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		  if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
				  (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			  MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
			  send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			  receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

			  nVertexS = nVertex[MarkerS];  nVertexR = nVertex[MarkerR];
			  nBufferS_Vector = nVertexS*nDim;        nBufferR_Vector = nVertexR*nDim;

			  /*--- Allocate Receive and send buffers  ---*/

			  Buffer_Receive_Coord_n1 = new su2double [nBufferR_Vector];
			  Buffer_Send_Coord_n1 = new su2double[nBufferS_Vector];

			  /*--- Copy the coordinates that should be sended ---*/

			  for (iVertex = 0; iVertex < nVertexS; iVertex++) {
				  iPoint = vertex[MarkerS][iVertex]->GetNode();
				  Coord_n1 = node[iPoint]->GetCoord_n1();
				  for (iDim = 0; iDim < nDim; iDim++)
					  Buffer_Send_Coord_n1[iDim*nVertexS+iVertex] = Coord_n1[iDim];
			  }

#ifdef HAVE_MPI
			  /*--- Send/Receive information using Sendrecv ---*/
			  SU2_MPI::Sendrecv(Buffer_Send_Coord_n1, nBufferS_Vector, MPI_DOUBLE, send_to,0,
					  Buffer_Receive_Coord_n1, nBufferR_Vector, MPI_DOUBLE, receive_from,0, MPI_COMM_WORLD, &status);
#else

			  /*--- Receive information without MPI ---*/
			  for (iVertex = 0; iVertex < nVertexR; iVertex++) {
				  for (iDim = 0; iDim < nDim; iDim++)
					  Buffer_Receive_Coord_n1[iDim*nVertexR+iVertex] = Buffer_Send_Coord_n1[iDim*nVertexR+iVertex];
			  }

#endif

			  /*--- Deallocate send buffer ---*/

			  delete [] Buffer_Send_Coord_n1;

			  /*--- Do the coordinate transformation ---*/

			  for (iVertex = 0; iVertex < nVertexR; iVertex++) {

				  /*--- Find point and its type of transformation ---*/

				  iPoint = vertex[MarkerR][iVertex]->GetNode();
				  iPeriodic_Index = vertex[MarkerR][iVertex]->GetRotation_Type();

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

				  /*--- Copy coordinates before performing transformation. ---*/

				  for (iDim = 0; iDim < nDim; iDim++)
					  newCoord_n1[iDim] = Buffer_Receive_Coord_n1[iDim*nVertexR+iVertex];

				  /*--- Rotate the coordinates. ---*/

				  if (nDim == 2) {
					  newCoord_n1[0] = (rotMatrix[0][0]*Buffer_Receive_Coord_n1[0*nVertexR+iVertex] +
							  rotMatrix[0][1]*Buffer_Receive_Coord_n1[1*nVertexR+iVertex]);
					  newCoord_n1[1] = (rotMatrix[1][0]*Buffer_Receive_Coord_n1[0*nVertexR+iVertex] +
							  rotMatrix[1][1]*Buffer_Receive_Coord_n1[1*nVertexR+iVertex]);
				  }
				  else {
					  newCoord_n1[0] = (rotMatrix[0][0]*Buffer_Receive_Coord_n1[0*nVertexR+iVertex] +
							  rotMatrix[0][1]*Buffer_Receive_Coord_n1[1*nVertexR+iVertex] +
							  rotMatrix[0][2]*Buffer_Receive_Coord_n1[2*nVertexR+iVertex]);
					  newCoord_n1[1] = (rotMatrix[1][0]*Buffer_Receive_Coord_n1[0*nVertexR+iVertex] +
							  rotMatrix[1][1]*Buffer_Receive_Coord_n1[1*nVertexR+iVertex] +
							  rotMatrix[1][2]*Buffer_Receive_Coord_n1[2*nVertexR+iVertex]);
					  newCoord_n1[2] = (rotMatrix[2][0]*Buffer_Receive_Coord_n1[0*nVertexR+iVertex] +
							  rotMatrix[2][1]*Buffer_Receive_Coord_n1[1*nVertexR+iVertex] +
							  rotMatrix[2][2]*Buffer_Receive_Coord_n1[2*nVertexR+iVertex]);
				  }

				  /*--- Copy transformed coordinates back into buffer. ---*/

				  node[iPoint]->SetCoord_n1(newCoord_n1);

			  }

			  /*--- Deallocate receive buffer. ---*/

			  delete [] Buffer_Receive_Coord_n1;

		  }

	  }

	  delete [] newCoord_n1;

  }

  /*--------------------------------------------------------------------------------------------------*/

}

void CPhysicalGeometry::SetPeriodicBoundary(CConfig *config) {
  
  unsigned short iMarker, jMarker, kMarker = 0, iPeriodic, iDim, nPeriodic = 0, VTK_Type;
  unsigned long iNode, iIndex, iVertex, iPoint, iElem, kElem;
  unsigned long jElem, kPoint = 0, jVertex = 0, jPoint = 0, pPoint = 0, nPointPeriodic, newNodes[4] = {0,0,0,0};
  vector<unsigned long>::iterator IterElem, IterPoint[MAX_NUMBER_PERIODIC][2];
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
  dx, dy, dz, rotCoord[3], epsilon = 1e-10, mindist = 1e6, *Coord_i, *Coord_j, dist = 0.0;
  bool isBadMatch = false;

  /*--- Check this dimensionalization ---*/

  vector<unsigned long> OldBoundaryElems[100];
  vector<unsigned long>::iterator IterNewElem[100];

  /*--- We only create the mirror structure for the second boundary ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      /*--- Evaluate the number of periodic boundary conditions ---*/
      nPeriodic++;
    }
  }
  bool *CreateMirror = new bool[nPeriodic+1];
  CreateMirror[0] = false;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (iPeriodic <= nPeriodic/2) CreateMirror[iPeriodic] = false;
    else CreateMirror[iPeriodic] = true;
  }
  
  /*--- Send an initial message to the console. ---*/
  cout << "Setting the periodic boundary conditions." << endl;
  
  /*--- Loop through each marker to find any periodic boundaries. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      /*--- Get marker index of the periodic donor boundary. ---*/
      jMarker = config->GetMarker_Periodic_Donor(config->GetMarker_All_TagBound(iMarker));
      
      /*--- Write some info to the console. ---*/
      cout << "Checking " << config->GetMarker_All_TagBound(iMarker);
      cout << " boundary against periodic donor, " << config->GetMarker_All_TagBound(jMarker) << ". ";
      
      /*--- Retrieve the supplied periodic information. ---*/
      center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
      angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
      trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
      
      /*--- Store (center+trans) as it is constant and will be added on. ---*/
      translation[0] = center[0] + trans[0];
      translation[1] = center[1] + trans[1];
      translation[2] = center[2] + trans[2];
      
      /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
      theta = angles[0];
      phi   = angles[1];
      psi   = angles[2];
      
      cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
      sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
      
      /*--- Compute the rotation matrix. Note that the implicit
       ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
      rotMatrix[0][0] = cosPhi*cosPsi;
      rotMatrix[1][0] = cosPhi*sinPsi;
      rotMatrix[2][0] = -sinPhi;
      
      rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
      rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
      rotMatrix[2][1] = sinTheta*cosPhi;
      
      rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
      rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
      rotMatrix[2][2] = cosTheta*cosPhi;
      
      /*--- Loop through all vertices and find/set the periodic point. ---*/
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        /*--- Retrieve node information for this boundary point. ---*/
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        Coord_i = node[iPoint]->GetCoord();
        
        /*--- Get the position vector from rot center to point. ---*/
        dx = Coord_i[0] - center[0];
        dy = Coord_i[1] - center[1];
        if (nDim == 3) {
          dz = Coord_i[2] - center[2];
        } else {
          dz = 0.0;
        }
        
        /*--- Compute transformed point coordinates. ---*/
        rotCoord[0] = rotMatrix[0][0]*dx
        + rotMatrix[0][1]*dy
        + rotMatrix[0][2]*dz + translation[0];
        
        rotCoord[1] = rotMatrix[1][0]*dx
        + rotMatrix[1][1]*dy
        + rotMatrix[1][2]*dz + translation[1];
        
        rotCoord[2] = rotMatrix[2][0]*dx
        + rotMatrix[2][1]*dy
        + rotMatrix[2][2]*dz + translation[2];
        
        /*--- Perform a search to find the closest donor point. ---*/
        mindist = 1e10;
        for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
          
          /*--- Retrieve information for this jPoint. ---*/
          jPoint = vertex[jMarker][jVertex]->GetNode();
          Coord_j = node[jPoint]->GetCoord();
          
          /*--- Check the distance between the computed periodic
           location and this jPoint. ---*/
          dist = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            dist += (Coord_j[iDim]-rotCoord[iDim])*(Coord_j[iDim]-rotCoord[iDim]);
          }
          dist = sqrt(dist);
          
          /*---  Store vertex information if this is the closest
           point found thus far. ---*/
          if (dist < mindist) { mindist = dist; pPoint = jPoint; }
        }
        
        /*--- Set the periodic point for this iPoint. ---*/
        vertex[iMarker][iVertex]->SetDonorPoint(pPoint, MASTER_NODE);
        
        /*--- Print warning if the nearest point was not within
         the specified tolerance. Computation will continue. ---*/
        if (mindist > epsilon) {
          isBadMatch = true;
          cout.precision(10);
          cout << endl;
          cout << "   Bad match for point " << iPoint << ".\tNearest";
          cout << " donor distance: " << scientific << mindist << ".";
        }
      }
      
      /*--- Print final warning when finding bad matches. ---*/
      if (isBadMatch) {
        cout << endl;
        cout << "\n !!! Warning !!!" << endl;
        cout << "Bad matches found. Computation will continue, but be cautious.\n";
      }
      cout << endl;
      isBadMatch = false;
      
    }
  
  /*--- Create a vector to identify the points that belong to each periodic boundary condition ---*/
  bool *PeriodicBC = new bool [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) PeriodicBC[iPoint] = false;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        PeriodicBC[iPoint] = true;
      }
  
  /*--- Determine the new points that must be added to each periodic boundary,
   note that only one of the boundaries require the extra data ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      
      /*--- An integer identify the periodic boundary condition --*/
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        /*--- iPoint is the original point on the surface and jPoint is the
         equivalent point in the other periodic surface ---*/
        iPoint = vertex[iMarker][iVertex]->GetNode();
        jPoint = vertex[iMarker][iVertex]->GetDonorPoint();
        
        /*--- First the case in which it is necessary to create a mirror set of elements ---*/
        if (CreateMirror[iPeriodic]) {
          /*--- Now we must determine the neighbor points (including indirect ones) to the periodic points
           and store all the information (in this list we do not include the points
           that already belong to the periodic boundary), we also add the elements that
           share a point with the periodic boundary condition ---*/
          for (iIndex = 0; iIndex < node[jPoint]->GetnElem(); iIndex++) {
            iElem = node[jPoint]->GetElem(iIndex);
            PeriodicElem[iPeriodic].push_back(iElem);
            for (unsigned short iNode = 0; iNode <	elem[iElem]->GetnNodes(); iNode ++) {
              kPoint = elem[iElem]->GetNode(iNode);
              if (!PeriodicBC[kPoint]) PeriodicPoint[iPeriodic][0].push_back(kPoint);
            }
          }
        }
        /*--- Second the case where no new element is added, neither points ---*/
        else {
          PeriodicPoint[iPeriodic][0].push_back(jPoint);
          PeriodicPoint[iPeriodic][1].push_back(iPoint);
        }
      }
    }
  }
  
  /*--- Sort the points that must be sended and delete repeated points ---*/
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      sort( PeriodicPoint[iPeriodic][0].begin(), PeriodicPoint[iPeriodic][0].end());
      IterPoint[iPeriodic][0] = unique( PeriodicPoint[iPeriodic][0].begin(), PeriodicPoint[iPeriodic][0].end());
      PeriodicPoint[iPeriodic][0].resize( IterPoint[iPeriodic][0] - PeriodicPoint[iPeriodic][0].begin() );
    }
  }
  
  /*--- Create a list of the points that receive the values (only the new points) ---*/
  nPointPeriodic = nPoint;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      for (iPoint = 0; iPoint < PeriodicPoint[iPeriodic][0].size(); iPoint++) {
        PeriodicPoint[iPeriodic][1].push_back(nPointPeriodic);
        nPointPeriodic++;
      }
    }
  }
  
  /*--- Sort the elements that must be replicated in the periodic boundary
   and delete the repeated elements ---*/
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      sort( PeriodicElem[iPeriodic].begin(), PeriodicElem[iPeriodic].end());
      IterElem = unique( PeriodicElem[iPeriodic].begin(), PeriodicElem[iPeriodic].end());
      PeriodicElem[iPeriodic].resize( IterElem - PeriodicElem[iPeriodic].begin() );
    }
  }
  
  /*--- Check all SEND points to see if they also lie on another boundary. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      
      /*--- iPoint is a node that lies on the current marker. ---*/
      iPoint = vertex[iMarker][iVertex]->GetNode();
      
      /*--- Search through SEND points to check for iPoint. ---*/
      for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
          
          /*--- jPoint is the SEND point. ---*/
          for (iElem = 0; iElem < PeriodicPoint[iPeriodic][0].size(); iElem++) {
            jPoint = PeriodicPoint[iPeriodic][0][iElem];
            
            /*--- If the two match, then jPoint lies on this boundary.
             However, we are concerned with the new points, so we
             will store kPoint instead. ---*/
            if (iPoint == jPoint) {
//              kPoint = PeriodicPoint[iPeriodic][1][iElem];
              
              /*--- We also want the type of boundary element that this point
               was within, so that we know what type of element to add
               built from the new points. ---*/
              bool isJPoint, isPeriodic;
              for (jElem = 0; jElem < nElem_Bound[iMarker]; jElem++) {
                isJPoint = false; isPeriodic = false;
                for (iNode = 0; iNode < bound[iMarker][jElem]->GetnNodes(); iNode++) {
                  if (bound[iMarker][jElem]->GetNode(iNode) == jPoint) isJPoint = true;
                  if (PeriodicBC[bound[iMarker][jElem]->GetNode(iNode)]) isPeriodic = true;
                }
                
                /*--- If both points were found, store this element. ---*/
                if (isJPoint && isPeriodic) {
                  OldBoundaryElems[iMarker].push_back(jElem);
                }
                
              }
              
            }
          }
        }
      }
    }
  }
  
  /*--- Sort the elements that must be added and remove duplicates. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    sort( OldBoundaryElems[iMarker].begin(), OldBoundaryElems[iMarker].end());
    IterNewElem[iMarker] = unique( OldBoundaryElems[iMarker].begin(), OldBoundaryElems[iMarker].end());
    OldBoundaryElems[iMarker].resize( IterNewElem[iMarker] - OldBoundaryElems[iMarker].begin() );
  }
  
  /*--- Create the new boundary elements. Points making up these new
   elements must either be SEND/RECEIVE or periodic points. ---*/
  nNewElem_Bound = new unsigned long[nMarker];
  newBound = new CPrimalGrid**[nMarker];
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    nNewElem_Bound[iMarker] = OldBoundaryElems[iMarker].size();
    newBound[iMarker] = new CPrimalGrid*[nNewElem_Bound[iMarker]];
    
    /*--- Loop through all new elements to be added. ---*/
    for (iElem = 0; iElem < nNewElem_Bound[iMarker]; iElem++) {
      jElem = OldBoundaryElems[iMarker][iElem];
      
      /*--- Loop through all nodes of this element. ---*/
      for (iNode = 0; iNode < bound[iMarker][jElem]->GetnNodes(); iNode++) {
        pPoint = bound[iMarker][jElem]->GetNode(iNode);
        
        /*--- Check if this node is a send point. If so, the corresponding
         receive point will be used in making the new boundary element. ---*/
        for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
          for (kElem = 0; kElem < PeriodicPoint[iPeriodic][0].size(); kElem++) {
            if (pPoint == PeriodicPoint[iPeriodic][0][kElem]) newNodes[iNode] = PeriodicPoint[iPeriodic][1][kElem];
          }
        }
        
        /*--- Check if this node is a periodic point. If so, the corresponding
         periodic point will be used in making the new boundary element. ---*/
        if (PeriodicBC[pPoint]) {
          
          /*--- Find the corresponding periodic point. ---*/
          for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++) {
            if (config->GetMarker_All_KindBC(jMarker) == PERIODIC_BOUNDARY) {
              for (iVertex = 0; iVertex < nVertex[jMarker]; iVertex++) {
                if (pPoint == vertex[jMarker][iVertex]->GetNode()) {kMarker = jMarker; jVertex = iVertex;}
              }
            }
          }
          newNodes[iNode] = vertex[kMarker][jVertex]->GetDonorPoint();
        }
      }
      
      /*--- Now instantiate the new element. ---*/
      VTK_Type = bound[iMarker][jElem]->GetVTK_Type();
      switch(VTK_Type) {
        case LINE:
          newBound[iMarker][iElem] = new CLine(newNodes[0], newNodes[1],2);
          break;
        case TRIANGLE:
          newBound[iMarker][iElem] = new CTriangle(newNodes[0], newNodes[1], newNodes[2],3);
          break;
        case QUADRILATERAL:
          newBound[iMarker][iElem] = new CQuadrilateral(newNodes[0], newNodes[1], newNodes[2], newNodes[3],3);
          break;
      }
      
    }
  }
  
  delete [] PeriodicBC;
  delete [] CreateMirror;
  
}

void CPhysicalGeometry::FindNormal_Neighbor(CConfig *config) {
  su2double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord, *Normal;
  unsigned long Point_Normal, jPoint;
  unsigned short iNeigh, iMarker, iDim;
  unsigned long iPoint, iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY ) {
      
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Normal = vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Compute closest normal neighbor, note that the normal are oriented inwards ---*/
        Point_Normal = 0; cos_max = -1.0;
        for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
          jPoint = node[iPoint]->GetPoint(iNeigh);
          scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
            scalar_prod += diff_coord*Normal[iDim];
            norm_vect += diff_coord*diff_coord;
            norm_Normal += Normal[iDim]*Normal[iDim];
          }
          norm_vect = sqrt(norm_vect);
          norm_Normal = sqrt(norm_Normal);
          cos_alpha = scalar_prod/(norm_vect*norm_Normal);
          
          /*--- Get maximum cosine ---*/
          if (cos_alpha >= cos_max) {
            Point_Normal = jPoint;
            cos_max = cos_alpha;
          }
        }
        vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
      }
    }
  }
}

void CPhysicalGeometry::SetGeometryPlanes(CConfig *config) {
  
  bool loop_on;
  unsigned short iMarker = 0;
  su2double auxXCoord, auxYCoord, auxZCoord,	*Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
  unsigned long jVertex, iVertex, ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;
  
  /*--- Compute the total number of points on the near-field ---*/
  nVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      nVertex_Wall += nVertex[iMarker];
  
  
  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  Xcoord = new su2double[nVertex_Wall];
  Ycoord = new su2double[nVertex_Wall];
  if (nDim == 3)	Zcoord = new su2double[nVertex_Wall];
  FaceArea = new su2double[nVertex_Wall];
  
  /*--- Copy the boundary information to an array ---*/
  iVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
        Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
        if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
        Face_Normal = vertex[iMarker][iVertex]->GetNormal();
        FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
        iVertex_Wall ++;
      }
  
  
  //vector<su2double> XCoordList;
  vector<su2double>::iterator IterXCoordList;
  
  for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
    XCoordList.push_back(Xcoord[iVertex]);
  
  sort( XCoordList.begin(), XCoordList.end());
  IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
  XCoordList.resize( IterXCoordList - XCoordList.begin() );
  
  /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
  Xcoord_plane.resize(XCoordList.size());
  Ycoord_plane.resize(XCoordList.size());
  if (nDim==3) Zcoord_plane.resize(XCoordList.size());
  FaceArea_plane.resize(XCoordList.size());
  Plane_points.resize(XCoordList.size());
  
  
  su2double dist_ratio;
  unsigned long iCoord;
  
  /*--- Distribute the values among the different PhiAngles ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (node[iPoint]->GetDomain()) {
      loop_on = true;
      for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
        dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
        if (dist_ratio >= 0 && dist_ratio <= 1.0) {
          if (dist_ratio <= 0.5) iCoord = ixCoord;
          else iCoord = ixCoord+1;
          Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
          Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
          if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
          FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
          Plane_points[iCoord].push_back(iPoint );
          loop_on = false;
        }
      }
    }
  }
  
  unsigned long auxPoint;
  /*--- Order the arrays in ascending values of y ---*/
  for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
    for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
      for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
        if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
          auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex+1]; Xcoord_plane[ixCoord][jVertex+1] = auxXCoord;
          auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex+1]; Ycoord_plane[ixCoord][jVertex+1] = auxYCoord;
          auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex+1]; Plane_points[ixCoord][jVertex+1] = auxPoint;
          if (nDim==3) {
            auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex+1]; Zcoord_plane[ixCoord][jVertex+1] = auxZCoord;
          }
          auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex+1]; FaceArea_plane[ixCoord][jVertex+1] = auxArea;
        }
  
  /*--- Delete structures ---*/
  delete[] Xcoord; delete[] Ycoord;
  if (Zcoord != NULL) delete[] Zcoord;
  delete[] FaceArea;
}

void CPhysicalGeometry::SetBoundSensitivity(CConfig *config) {
  unsigned short iMarker, icommas;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double Sensitivity;
  bool *PointInDomain;
  
#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  nPointLocal = nPoint;
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nPointGlobal = nPointLocal;
#endif
  
  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];
  
  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    if (config->GetMarker_All_DV(iMarker) == YES)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        /*--- The sensitivity file uses the global numbering ---*/
        iPoint = node[vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();

        if (vertex[iMarker][iVertex]->GetNode() < GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          vertex[iMarker][iVertex]->SetAuxVar(0.0);
        }
      }
  
  /*--- Time-average any unsteady surface sensitivities ---*/
  
  unsigned long iExtIter, nExtIter;
  su2double delta_T, total_T;
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    nExtIter = config->GetUnst_AdjointIter();
    delta_T  = config->GetDelta_UnstTime();
    total_T  = (su2double)nExtIter*delta_T;
  } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    
    /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
    
    su2double period = config->GetHarmonicBalance_Period();
    nExtIter  = config->GetnTimeInstances();
    delta_T   = period/(su2double)nExtIter;
    total_T   = period;
    
  } else {
    nExtIter = 1;
    delta_T  = 1.0;
    total_T  = 1.0;
  }
  
  for (iExtIter = 0; iExtIter < nExtIter; iExtIter++) {
    
    /*--- Prepare to read surface sensitivity files (CSV) ---*/
    
    string text_line;
    ifstream Surface_file;
    char buffer[50];
    char cstr[MAX_STRING_SIZE];
    string surfadj_filename = config->GetSurfAdjCoeff_FileName();
    strcpy (cstr, surfadj_filename.c_str());
    
    /*--- Write file name with extension if unsteady or steady ---*/
    if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE)
    	SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));

    if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
        (config->GetUnsteady_Simulation() == HARMONIC_BALANCE)) {
      if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
      if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
    }
    else
      SPRINTF (buffer, ".csv");
    
    strcat (cstr, buffer);
    
    /*--- Read the sensitivity file ---*/
    
    string::size_type position;
    
    Surface_file.open(cstr, ios::in);
    
    /*--- Read extra inofmration ---*/
    
    getline(Surface_file, text_line);
    text_line.erase (0,9);
    su2double AoASens = atof(text_line.c_str());
    config->SetAoA_Sens(AoASens);
    
    /*--- File header ---*/
    
    getline(Surface_file, text_line);
    
    while (getline(Surface_file, text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if (position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);
      point_line >> iPoint >> Sensitivity;
      
      if (PointInDomain[iPoint]) {
        
        /*--- Find the vertex for the Point and Marker ---*/
        
        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];
        
        /*--- Increment the auxiliary variable with the contribution of
         this unsteady timestep. For steady problems, this reduces to
         a single sensitivity value multiplied by 1.0. ---*/
        
        vertex[iMarker][iVertex]->AddAuxVar(Sensitivity*(delta_T/total_T));
      }
      
    }
    Surface_file.close();
  }
  
  delete[] Point2Vertex;
  delete[] PointInDomain;
  
}

void CPhysicalGeometry::SetSensitivity(CConfig *config) {
  
  ifstream restart_file;
  string filename = config->GetSolution_AdjFileName();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool sst = config->GetKind_Turb_Model() == SST;
  bool sa = config->GetKind_Turb_Model() == SA;
  bool grid_movement = config->GetGrid_Movement();
  bool wrt_residuals = config->GetWrt_Residuals();
  su2double Sens, dull_val, AoASens;
  unsigned short nExtIter, iDim;
  unsigned long iPoint, index;
  string::size_type position;

  Sensitivity = new su2double[nPoint*nDim];

  if (config->GetUnsteady_Simulation()) {
    nExtIter = config->GetnExtIter();
  }else {
    nExtIter = 1;
  }
    int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short skipVar = nDim, skipMult = 1;

  if (wrt_residuals) { skipMult = 2; }
  if (incompressible) { skipVar += skipMult*(nDim+1); }
  if (compressible)   { skipVar += skipMult*(nDim+2); }
  if (sst)            { skipVar += skipMult*2;}
  if (sa)             { skipVar += skipMult*1;}
  if (grid_movement)  { skipVar += nDim;}
  
  /*--- Sensitivity in normal direction ---*/
  
  skipVar += 1;
  
  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  
  map<unsigned long,unsigned long> Global2Local;
  map<unsigned long,unsigned long>::const_iterator MI;
  
  /*--- Now fill array with the transform values only for local points ---*/
  for(iPoint = 0; iPoint < nPointDomain; iPoint++)
    Global2Local[node[iPoint]->GetGlobalIndex()] = iPoint;
  
  /*--- Read all lines in the restart file ---*/
  long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;
  
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Sensitivity[iPoint*nDim+iDim] = 0.0;
    }
  }
  

  iPoint_Global = 0;

  filename = config->GetSolution_AdjFileName();

  filename = config->GetObjFunc_Extension(filename);

  if (config->GetUnsteady_Simulation()) {
    filename = config->GetUnsteady_FileName(filename, nExtIter-1);
  }

  restart_file.open(filename.data(), ios::in);
  if (restart_file.fail()) {
    cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
    exit(EXIT_FAILURE);
  }
  
  if (rank == MASTER_NODE)
    cout << "Reading in sensitivity at iteration " << nExtIter-1 << "."<< endl;
  /*--- The first line is the header ---*/
  getline (restart_file, text_line);
  
  for (iPoint_Global = 0; iPoint_Global < GetGlobal_nPointDomain(); iPoint_Global++ ) {

    getline (restart_file, text_line);

  	istringstream point_line(text_line);

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/
    
    MI = Global2Local.find(iPoint_Global);
    if (MI != Global2Local.end()) {
      
      iPoint_Local = Global2Local[iPoint_Global];
      
      point_line >> index;
      for (iDim = 0; iDim < skipVar; iDim++) { point_line >> dull_val;}
      for (iDim = 0; iDim < nDim; iDim++) {
        point_line >> Sens;
        Sensitivity[iPoint_Local*nDim+iDim] = Sens;
      }
    }

  }
  
  /*--- Read AoA sensitivity ---*/
  
  while (getline (restart_file, text_line)) {
    position = text_line.find ("SENS_AOA=",0);
    if (position != string::npos) {
      text_line.erase (0,9); AoASens = atof(text_line.c_str());
      config->SetAoA_Sens(AoASens);
    }
  }
  
  restart_file.close();

}

void CPhysicalGeometry::Check_Periodicity(CConfig *config) {
  
  bool isPeriodic = false;
  unsigned long iVertex;
  unsigned short iMarker, RotationKind, nPeriodicR = 0, nPeriodicS = 0;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Check for the presence of any periodic BCs ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        RotationKind = vertex[iMarker][iVertex]->GetRotation_Type();
        if (RotationKind > 0) nPeriodicS++;
      }
    }
  }
#ifndef HAVE_MPI
  nPeriodicR = nPeriodicS;
#else
  SU2_MPI::Allreduce(&nPeriodicS, &nPeriodicR, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (nPeriodicR != 0) isPeriodic = true;
  
  if (isPeriodic && (config->GetnMGLevels() > 0)) {
    if (rank == MASTER_NODE)
      cout << "WARNING: Periodicity has been detected. Disabling multigrid. "<< endl;
    config->SetMGLevels(0);
  }
  
}

su2double CPhysicalGeometry::Compute_MaxThickness(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, jVertex, n, Trailing_Point, Leading_Point;
  su2double Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, zp1, zpn, MaxThickness_Value = 0, Thickness, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, MaxDistance, Distance, AoA;
  vector<su2double> Xcoord, Ycoord, Zcoord, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  
  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }
  
  AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;
  
  /*--- Translate to the origin ---*/
  
  Xcoord_Trailing = Xcoord_Airfoil[0];
  Ycoord_Trailing = Ycoord_Airfoil[0];
  Zcoord_Trailing = Zcoord_Airfoil[0];
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
    Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
    Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
  }
  
  /*--- Rotate the airfoil ---*/
  
  ValCos = cos(AoA*PI_NUMBER/180.0);
  ValSin = sin(AoA*PI_NUMBER/180.0);
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    XValue = Xcoord_Airfoil_[iVertex];
    ZValue = Zcoord_Airfoil_[iVertex];
    
    Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
    Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
    
  }
  
  /*--- Identify upper and lower side, and store the value of the normal --*/
  
  for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
    Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
    Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
    Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
    Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
    
    BiNormal[0] = Plane_Normal[0];
    BiNormal[1] = Plane_Normal[1];
    BiNormal[2] = Plane_Normal[2];
    Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
    BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
    
    Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
    Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
    Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
    
    Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
    
    unsigned short index = 2;
    if ((config->GetAxis_Stations() == Z_AXIS) && (nDim == 3)) index = 0;
    
    if (Normal[index] >= 0.0) {
      Xcoord.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord.push_back(Zcoord_Airfoil_[iVertex]);
    }
    
  }
  
  /*--- Order the arrays using the X component ---*/
  
  for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord.size() - 1 - iVertex; jVertex++) {
      if (Xcoord[jVertex] > Xcoord[jVertex+1]) {
        auxXCoord = Xcoord[jVertex]; Xcoord[jVertex] = Xcoord[jVertex+1]; Xcoord[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord[jVertex]; Ycoord[jVertex] = Ycoord[jVertex+1]; Ycoord[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord[jVertex]; Zcoord[jVertex] = Zcoord[jVertex+1]; Zcoord[jVertex+1] = auxZCoord;
      }
    }
  }
  
  n = Xcoord.size();
  zp1 = (Zcoord[1]-Zcoord[0])/(Xcoord[1]-Xcoord[0]);
  zpn = (Zcoord[n-1]-Zcoord[n-2])/(Xcoord[n-1]-Xcoord[n-2]);
  Z2coord.resize(n+1);
  SetSpline(Xcoord, Zcoord, n, zp1, zpn, Z2coord);
  
  /*--- Compute the thickness (we add a fabs because we can not guarantee the
   right sorting of the points and the upper and/or lower part of the airfoil is not well defined) ---*/
  
  MaxThickness_Value = 0.0;
  for (iVertex = 0; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    if (Zcoord_Normal[iVertex] < 0.0) {
      Thickness = fabs(Zcoord_Airfoil_[iVertex] - GetSpline(Xcoord, Zcoord, Z2coord, n, Xcoord_Airfoil_[iVertex]));
      if (Thickness > MaxThickness_Value) { MaxThickness_Value = Thickness; }
    }
  }
  
  return MaxThickness_Value;
  
}

su2double CPhysicalGeometry::Compute_Dihedral(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                                              su2double *LeadingEdge_i, su2double *TrailingEdge_i) {

  su2double Dihedral = 0.0, Dihedral_Leading = 0.0, Dihedral_Trailing = 0.0;

  Dihedral_Leading = atan((LeadingEdge_i[2] - LeadingEdge_im1[2]) / (LeadingEdge_i[1] - LeadingEdge_im1[1]))*180/PI_NUMBER;
  Dihedral_Trailing = atan((TrailingEdge_i[2] - TrailingEdge_im1[2]) / (TrailingEdge_i[1] - TrailingEdge_im1[1]))*180/PI_NUMBER;

  Dihedral = 0.5*(Dihedral_Leading + Dihedral_Trailing);

  return Dihedral;

}

su2double CPhysicalGeometry::Compute_Curvature(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                                               su2double *LeadingEdge_i, su2double *TrailingEdge_i,
                                               su2double *LeadingEdge_ip1, su2double *TrailingEdge_ip1) {

  su2double Curvature = 0.0, Curvature_Leading = 0.0, Curvature_Trailing = 0.0;

  su2double A[2], B[2], C[2], BC[2], AB[2], AC[2], BC_MOD, AB_MOD,  AC_MOD, AB_CROSS_AC;

  A[0] = LeadingEdge_im1[1];     A[1] = LeadingEdge_im1[2];
  B[0] = LeadingEdge_i[1];           B[1] = LeadingEdge_i[2];
  C[0] = LeadingEdge_ip1[1];      C[1] = LeadingEdge_ip1[2];

  BC[0] = C[0] - B[0]; BC[1] = C[1] - B[1];
  AB[0] = B[0] - A[0]; AB[1] = B[1] - A[1];
  AC[0] = C[0] - A[0]; AC[1] = C[1] - A[1];
  BC_MOD = sqrt(BC[0]*BC[0] + BC[1]*BC[1] );
  AB_MOD = sqrt(AB[0]*AB[0] + AB[1]*AB[1] );
  AC_MOD = sqrt(AC[0]*AC[0] + AC[1]*AC[1] );
  AB_CROSS_AC = AB[0]* AC[1] - AB[1]* AC[0];

  Curvature_Leading = fabs(1.0/(0.5*BC_MOD*AB_MOD*AC_MOD/AB_CROSS_AC));

  A[0] = TrailingEdge_im1[1];      A[1] = TrailingEdge_im1[2];
  B[0] = TrailingEdge_i[1];                B[1] = TrailingEdge_i[2];
  C[0] = TrailingEdge_ip1[1];      C[1] = TrailingEdge_ip1[2];

  BC[0] = C[0] - B[0]; BC[1] = C[1] - B[1];
  AB[0] = B[0] - A[0]; AB[1] = B[1] - A[1];
  AC[0] = C[0] - A[0]; AC[1] = C[1] - A[1];
  BC_MOD = sqrt(BC[0]*BC[0] + BC[1]*BC[1] );
  AB_MOD = sqrt(AB[0]*AB[0] + AB[1]*AB[1] );
  AC_MOD = sqrt(AC[0]*AC[0] + AC[1]*AC[1] );
  AB_CROSS_AC = AB[0]* AC[1] - AB[1]* AC[0];

  Curvature_Trailing = fabs(1.0/(0.5*BC_MOD*AB_MOD*AC_MOD/AB_CROSS_AC));

  Curvature = 0.5*(Curvature_Leading + Curvature_Trailing);

  return Curvature;

}

su2double CPhysicalGeometry::Compute_Twist(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {
  unsigned long iVertex, Trailing_Point, Leading_Point;
  su2double MaxDistance, Distance, Twist = 0.0;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/

  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }
  
  Twist = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;

  if (fabs(Twist) < 0.01) Twist = 0.0;

  return Twist;

}

void CPhysicalGeometry::Compute_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection,
                                                vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, Trailing_Point, Leading_Point;
  su2double MaxDistance, Distance;

  /*--- Find the leading and trailing edges and compute the angle of attack ---*/

  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {

    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }

  LeadingEdge[0] = Xcoord_Airfoil[Leading_Point];
  LeadingEdge[1] = Ycoord_Airfoil[Leading_Point];
  LeadingEdge[2] = Zcoord_Airfoil[Leading_Point];
  
  TrailingEdge[0] = Xcoord_Airfoil[Trailing_Point];
  TrailingEdge[1] = Ycoord_Airfoil[Trailing_Point];
  TrailingEdge[2] = Zcoord_Airfoil[Trailing_Point];
  
}

su2double CPhysicalGeometry::Compute_Chord(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {
  unsigned long iVertex, Trailing_Point;
  su2double MaxDistance, Distance, Chord = 0.0;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  MaxDistance = 0.0; Trailing_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; }
  }
  
  Chord = MaxDistance;
  
  return Chord;
  
}

su2double CPhysicalGeometry::Compute_Thickness(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, su2double Location, CConfig *config, vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, jVertex, n_Upper, n_Lower, Trailing_Point, Leading_Point;
  su2double Thickness_Location, Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, Thickness_Value = 0.0, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, zp1, zpn, Chord, MaxDistance, Distance, AoA;
  vector<su2double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Z2coord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord_Lower, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  
  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }
  
  AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;
  Chord = MaxDistance;
  
  /*--- Translate to the origin ---*/
  
  Xcoord_Trailing = Xcoord_Airfoil[0];
  Ycoord_Trailing = Ycoord_Airfoil[0];
  Zcoord_Trailing = Zcoord_Airfoil[0];
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
    Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
    Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
  }
  
  /*--- Rotate the airfoil ---*/
  
  ValCos = cos(AoA*PI_NUMBER/180.0);
  ValSin = sin(AoA*PI_NUMBER/180.0);
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    XValue = Xcoord_Airfoil_[iVertex];
    ZValue = Zcoord_Airfoil_[iVertex];
    
    Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
    Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
  }
  
  /*--- Identify upper and lower side, and store the value of the normal --*/
  
  for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
    Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
    Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
    Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
    Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
    
    BiNormal[0] = Plane_Normal[0];
    BiNormal[1] = Plane_Normal[1];
    BiNormal[2] = Plane_Normal[2];
    Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
    BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
    
    Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
    Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
    Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
    
    Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
    
    unsigned short index = 2;
    if ((config->GetAxis_Stations() == Z_AXIS) && (nDim == 3)) index = 0;
    
    if (Normal[index] >= 0.0) {
      Xcoord_Upper.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord_Upper.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord_Upper.push_back(Zcoord_Airfoil_[iVertex]);
    }
    else {
      Xcoord_Lower.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord_Lower.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord_Lower.push_back(Zcoord_Airfoil_[iVertex]);
    }
    
  }
  
  /*--- Order the arrays using the X component ---*/
  
  for (iVertex = 0; iVertex < Xcoord_Upper.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord_Upper.size() - 1 - iVertex; jVertex++) {
      if (Xcoord_Upper[jVertex] > Xcoord_Upper[jVertex+1]) {
        auxXCoord = Xcoord_Upper[jVertex]; Xcoord_Upper[jVertex] = Xcoord_Upper[jVertex+1]; Xcoord_Upper[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord_Upper[jVertex]; Ycoord_Upper[jVertex] = Ycoord_Upper[jVertex+1]; Ycoord_Upper[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord_Upper[jVertex]; Zcoord_Upper[jVertex] = Zcoord_Upper[jVertex+1]; Zcoord_Upper[jVertex+1] = auxZCoord;
      }
    }
  }
  
  /*--- Order the arrays using the X component ---*/
  
  for (iVertex = 0; iVertex < Xcoord_Lower.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord_Lower.size() - 1 - iVertex; jVertex++) {
      if (Xcoord_Lower[jVertex] > Xcoord_Lower[jVertex+1]) {
        auxXCoord = Xcoord_Lower[jVertex]; Xcoord_Lower[jVertex] = Xcoord_Lower[jVertex+1]; Xcoord_Lower[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord_Lower[jVertex]; Ycoord_Lower[jVertex] = Ycoord_Lower[jVertex+1]; Ycoord_Lower[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord_Lower[jVertex]; Zcoord_Lower[jVertex] = Zcoord_Lower[jVertex+1]; Zcoord_Lower[jVertex+1] = auxZCoord;
      }
    }
  }
  
  n_Upper = Xcoord_Upper.size();
  zp1 = (Zcoord_Upper[1]-Zcoord_Upper[0])/(Xcoord_Upper[1]-Xcoord_Upper[0]);
  zpn = (Zcoord_Upper[n_Upper-1]-Zcoord_Upper[n_Upper-2])/(Xcoord_Upper[n_Upper-1]-Xcoord_Upper[n_Upper-2]);
  Z2coord_Upper.resize(n_Upper+1);
  SetSpline(Xcoord_Upper, Zcoord_Upper, n_Upper, zp1, zpn, Z2coord_Upper);
  
  n_Lower = Xcoord_Lower.size();
  zp1 = (Zcoord_Lower[1]-Zcoord_Lower[0])/(Xcoord_Lower[1]-Xcoord_Lower[0]);
  zpn = (Zcoord_Lower[n_Lower-1]-Zcoord_Lower[n_Lower-2])/(Xcoord_Lower[n_Lower-1]-Xcoord_Lower[n_Lower-2]);
  Z2coord_Lower.resize(n_Lower+1);
  SetSpline(Xcoord_Lower, Zcoord_Lower, n_Lower, zp1, zpn, Z2coord_Lower);
  
  /*--- Compute the thickness (we add a fabs because we can not guarantee the
   right sorting of the points and the upper and/or lower part of the airfoil is not well defined) ---*/
  
  Thickness_Location = - Chord*(1.0-Location);
  
  Thickness_Value = fabs(GetSpline(Xcoord_Upper, Zcoord_Upper, Z2coord_Upper, n_Upper, Thickness_Location) - GetSpline(Xcoord_Lower, Zcoord_Lower, Z2coord_Lower, n_Lower, Thickness_Location));
  
  return Thickness_Value;
  
}

su2double CPhysicalGeometry::Compute_Area(su2double *Plane_P0, su2double *Plane_Normal, unsigned short iSection, CConfig *config,
                                          vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {
  unsigned long iVertex;
  su2double Area_Value = 0.0;
  vector<su2double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
  su2double DeltaZ, DeltaX, X, Z;
  
  /*--- Use the Green theorem to evaluate the area (the points have been sortered),
   we assume that the airfoil is in the X-Z plane  ---*/
  
  Area_Value = 0.0;
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size()-1; iVertex++) {
    X = 0.5*(Xcoord_Airfoil[iVertex]+Xcoord_Airfoil[iVertex+1]);
    Z = 0.5*(Zcoord_Airfoil[iVertex]+Zcoord_Airfoil[iVertex+1]);
    DeltaX = Xcoord_Airfoil[iVertex+1] - Xcoord_Airfoil[iVertex];
    DeltaZ = Zcoord_Airfoil[iVertex+1] - Zcoord_Airfoil[iVertex];
    Area_Value += 0.5*( X*DeltaZ-Z*DeltaX);
  }
  
  X = 0.5*(Xcoord_Airfoil[Xcoord_Airfoil.size()-1]+Xcoord_Airfoil[0]);
  Z = 0.5*(Zcoord_Airfoil[Xcoord_Airfoil.size()-1]+Zcoord_Airfoil[0]);
  DeltaX = Xcoord_Airfoil[0] - Xcoord_Airfoil[Xcoord_Airfoil.size()-1];
  DeltaZ = Zcoord_Airfoil[0] - Zcoord_Airfoil[Xcoord_Airfoil.size()-1];
  Area_Value += 0.5 * (X*DeltaZ-Z*DeltaX);
  
  Area_Value = fabs(Area_Value);
  
  return Area_Value;
  
}

void CPhysicalGeometry::Compute_Wing(CConfig *config, bool original_surface,
                                     su2double &Wing_Volume, su2double &Wing_MinMaxThickness, su2double &Wing_MaxChord,
                                     su2double &Wing_MinToC, su2double &Wing_MaxTwist, su2double &Wing_MaxCurvature,
                                     su2double &Wing_MaxDihedral) {

  unsigned short iPlane, iDim, nPlane = 0;
  unsigned long iVertex;
  su2double MinPlane, MaxPlane, MinXCoord, MaxXCoord, dPlane, *Area, *MaxThickness, *ToC, *Chord, *Twist, *Curvature, *Dihedral, SemiSpan;
  vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  ofstream Wing_File, Section_File;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  
  /*--- Make a large number of section cuts for approximating volume ---*/
  
  nPlane = config->GetnWingStations();
  SemiSpan = config->GetSemiSpan();
  
  /*--- Allocate memory for the section cutting ---*/
  
  Area = new su2double [nPlane];
  MaxThickness = new su2double [nPlane];
  Chord = new su2double [nPlane];
  ToC = new su2double [nPlane];
  Twist = new su2double [nPlane];
  Curvature = new su2double [nPlane];
  Dihedral = new su2double [nPlane];

  su2double **LeadingEdge = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    LeadingEdge[iPlane] = new su2double[nDim];
  
  su2double **TrailingEdge = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    TrailingEdge[iPlane] = new su2double[nDim];

  su2double **Plane_P0 = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    Plane_P0[iPlane] = new su2double[nDim];
  
  su2double **Plane_Normal = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    Plane_Normal[iPlane] = new su2double[nDim];
  
  MinPlane = config->GetSection_WingBounds(0); MaxPlane = config->GetSection_WingBounds(1);
  MinXCoord = -1E6; MaxXCoord = 1E6;
  dPlane = fabs((MaxPlane - MinPlane)/su2double(nPlane-1));

  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
    Plane_Normal[iPlane][1] = 0.0;    Plane_P0[iPlane][1] = 0.0;
    Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;
    Plane_Normal[iPlane][config->GetAxis_Stations()] = 1.0;
    Plane_P0[iPlane][config->GetAxis_Stations()] = MinPlane + iPlane*dPlane;
  }
  
  /*--- Allocate some vectors for storing airfoil coordinates ---*/
  
  Xcoord_Airfoil   = new vector<su2double>[nPlane];
  Ycoord_Airfoil   = new vector<su2double>[nPlane];
  Zcoord_Airfoil   = new vector<su2double>[nPlane];
  Variable_Airfoil = new vector<su2double>[nPlane];
  
  /*--- Create the section slices through the geometry ---*/

  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane],
                           MinXCoord, MaxXCoord, NULL, Xcoord_Airfoil[iPlane],
                           Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                           Variable_Airfoil[iPlane], original_surface, config);
  }

  /*--- Compute the area at each section ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Write an output file---*/

    if (config->GetOutput_FileFormat() == PARAVIEW) {
      Wing_File.open("Wing_Distribution.csv", ios::out);
      if (config->GetSystemMeasurements() == US)
        Wing_File << "\"yCoord/SemiSpan\",\"Area (in^2)\",\"Max. Thickness (in)\",\"Chord (in)\",\"t_max/c\",\"Twist (deg)\",\"Curvature (1/in)\",\"Dihedral (deg)\",\"Leading Edge X (in)\",\"Leading Edge Y (in)\",\"Leading Edge Z (in)\",\"Trailing Edge X (in)\",\"Trailing Edge Y (in)\",\"Trailing Edge Z (in)\"" << endl;
      else
        Wing_File << "\"yCoord/SemiSpan\",\"Area (m^2)\",\"Max. Thickness (m)\",\"Chord (m)\",\"t_max/c\",\"Twist (deg)\",\"Curvature (1/in)\",\"Dihedral (deg)\",\"Leading Edge X (m)\",\"Leading Edge Y (m)\",\"Leading Edge Z (m)\",\"Trailing Edge X (m)\",\"Trailing Edge Y (m)\",\"Trailing Edge Z (m)\"" << endl;
    }
    else {
      Wing_File.open("wing_description.dat", ios::out);
      Wing_File << "TITLE = \"Wing description\"" << endl;
      if (config->GetSystemMeasurements() == US)
        Wing_File << "VARIABLES = \"<greek>h</greek>\",\"Area (in<sup>2</sup>)\",\"Max. Thickness (in)\",\"Chord (in)\",\"t<sub>max</sub>/c\",\"Twist (deg)\",\"Curvature (1/in)\",\"Dihedral (deg)\",\"Leading Edge X (in)\",\"Leading Edge Y (in)\",\"Leading Edge Z (in)\",\"Trailing Edge X (in)\",\"Trailing Edge Y (in)\",\"Trailing Edge Z (in)\"" << endl;
      else
        Wing_File << "VARIABLES = \"<greek>h</greek>\",\"Area (m<sup>2</sup>)\",\"Max. Thickness (m)\",\"Chord (m)\",\"t<sub>max</sub>/c\",\"Twist (deg)\",\"Curvature (1/m)\",\"Dihedral (deg)\",\"Leading Edge X (m)\",\"Leading Edge Y (m)\",\"Leading Edge Z (m)\",\"Trailing Edge X (m)\",\"Trailing Edge Y (m)\",\"Trailing Edge Z (m)\"" << endl;
      Wing_File << "ZONE T= \"Baseline wing\"" << endl;
    }


    /*--- Evaluate  geometrical quatities that do not require any kind of filter, local to each point ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        LeadingEdge[iPlane][iDim] = 0.0;
        TrailingEdge[iPlane][iDim] = 0.0;
      }

      Area[iPlane] = 0.0;
      MaxThickness[iPlane] = 0.0;
      Chord[iPlane] = 0.0;
      ToC[iPlane] = 0.0;
      Twist[iPlane] = 0.0;

      if (Xcoord_Airfoil[iPlane].size() != 0) {

        Compute_LeadingTrailing(LeadingEdge[iPlane], TrailingEdge[iPlane], Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Area[iPlane] = Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        MaxThickness[iPlane]= Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Chord[iPlane] = Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        ToC[iPlane] = MaxThickness[iPlane] / Chord[iPlane];

        Twist[iPlane] = Compute_Twist(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

      }

    }

    /*--- Evaluate  geometrical quatities that have been computed using a filtered value (they depend on more than one point) ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      Curvature[iPlane] = 0.0;
      Dihedral[iPlane] = 0.0;

      if (Xcoord_Airfoil[iPlane].size() != 0) {

        if ((iPlane == 0) || (iPlane == nPlane-1)) Curvature[iPlane] = 0.0;
        else Curvature[iPlane] = Compute_Curvature(LeadingEdge[iPlane-1], TrailingEdge[iPlane-1],
                                                   LeadingEdge[iPlane], TrailingEdge[iPlane],
                                                   LeadingEdge[iPlane+1], TrailingEdge[iPlane+1]);

        if (iPlane == 0) Dihedral[iPlane] = 0.0;
        else Dihedral[iPlane] = Compute_Dihedral(LeadingEdge[iPlane-1], TrailingEdge[iPlane-1],
                                                 LeadingEdge[iPlane], TrailingEdge[iPlane]);

      }

    }

    /*--- Set the curvature and dihedral angles at the extremes ---*/

    if (nPlane > 1) {
      if ((Xcoord_Airfoil[0].size() != 0) && (Xcoord_Airfoil[1].size() != 0)) {
        Curvature[0] = Curvature[1]; Dihedral[0] = Dihedral[1];
      }
      if ((Xcoord_Airfoil[nPlane-1].size() != 0) && (Xcoord_Airfoil[nPlane-2].size() != 0)) {
        Curvature[nPlane-1] = Curvature[nPlane-2];
      }
    }

    /*--- Plot the geometrical quatities ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
    	if (Xcoord_Airfoil[iPlane].size() != 0) {
    		if (config->GetOutput_FileFormat() == PARAVIEW) {
    			Wing_File  << Ycoord_Airfoil[iPlane][0]/SemiSpan <<", "<< Area[iPlane]  <<", "<< MaxThickness[iPlane]  <<", "<< Chord[iPlane]  <<", "<< ToC[iPlane]
    			                                                                                                                                            <<", "<< Twist[iPlane] <<", "<< Curvature[iPlane] <<", "<< Dihedral[iPlane]
    			                                                                                                                                                                                                                <<", "<< LeadingEdge[iPlane][0] <<", "<< LeadingEdge[iPlane][1]  <<", "<< LeadingEdge[iPlane][2]
    			                                                                                                                                                                                                                                                                                                              <<", "<< TrailingEdge[iPlane][0] <<", "<< TrailingEdge[iPlane][1]  <<", "<< TrailingEdge[iPlane][2]  << endl;
    		}
    		else  {
    			Wing_File  << Ycoord_Airfoil[iPlane][0]/SemiSpan <<" "<< Area[iPlane]  <<" "<< MaxThickness[iPlane]  <<" "<< Chord[iPlane]  <<" "<< ToC[iPlane]
    			                                                                                                                                        <<" "<< Twist[iPlane] <<" "<< Curvature[iPlane]  <<" "<< Dihedral[iPlane]
    			                                                                                                                                                                                                          <<" "<< LeadingEdge[iPlane][0] <<" "<< LeadingEdge[iPlane][1]  <<" "<< LeadingEdge[iPlane][2]
    			                                                                                                                                                                                                                                                                                                     <<" "<< TrailingEdge[iPlane][0] <<" "<< TrailingEdge[iPlane][1]  <<" "<< TrailingEdge[iPlane][2] << endl;
    		}
    	}
    }

    Wing_File.close();


    Section_File.open("wing_slices.dat", ios::out);

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      if (iPlane == 0) {
        Section_File << "TITLE = \"Aircraft Slices\"" << endl;
        if (config->GetSystemMeasurements() == US)
          Section_File << "VARIABLES = \"x (in)\", \"y (in)\", \"z (in)\", \"x<sub>2D</sub>/c\", \"y<sub>2D</sub>/c\"" << endl;
        else Section_File << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\", \"x<sub>2D</sub>/c\", \"y<sub>2D</sub>/c\"" << endl;
      }

      if (Xcoord_Airfoil[iPlane].size() !=0) {

        Section_File << "ZONE T=\"<greek>h</greek> = " << Ycoord_Airfoil[iPlane][0]/SemiSpan << " \", I= " << Xcoord_Airfoil[iPlane].size() << ", F=POINT" << endl;

        for (iVertex = 0; iVertex < Xcoord_Airfoil[iPlane].size(); iVertex++) {

          /*--- Move to the origin  ---*/

          su2double XValue_ = Xcoord_Airfoil[iPlane][iVertex] - LeadingEdge[iPlane][0];
          su2double ZValue_ = Zcoord_Airfoil[iPlane][iVertex] - LeadingEdge[iPlane][2];

          /*--- Rotate the airfoil and divide by the chord ---*/

          su2double ValCos = cos(Twist[iPlane]*PI_NUMBER/180.0);
          su2double ValSin = sin(Twist[iPlane]*PI_NUMBER/180.0);

          su2double XValue = (XValue_*ValCos - ZValue_*ValSin) / Chord[iPlane];
          su2double ZValue = (ZValue_*ValCos + XValue_*ValSin) / Chord[iPlane];

          /*--- Write the file ---*/

          Section_File  << Xcoord_Airfoil[iPlane][iVertex] << " " << Ycoord_Airfoil[iPlane][iVertex] << " " << Zcoord_Airfoil[iPlane][iVertex] << " " << XValue  << " " << ZValue << endl;
        }
      }

    }

    Section_File.close();


    /*--- Compute the wing volume using a composite Simpson's rule ---*/

    Wing_Volume = 0.0;
    for (iPlane = 0; iPlane < nPlane-2; iPlane+=2) {
      if (Xcoord_Airfoil[iPlane].size() != 0) {
        Wing_Volume += (1.0/3.0)*dPlane*(Area[iPlane] + 4.0*Area[iPlane+1] + Area[iPlane+2]);
      }
    }
    if (Xcoord_Airfoil[nPlane-1].size() != 0) {
      Wing_Volume /= (Ycoord_Airfoil[nPlane-1][0] - Ycoord_Airfoil[0][0]);
    }

    /*--- Evaluate Max and Min quantities ---*/

    Wing_MinMaxThickness = 1E6; Wing_MaxChord = -1E6; Wing_MinToC = 1E6;
    Wing_MaxTwist = -1E6; Wing_MaxCurvature = -1E6; Wing_MaxDihedral = -1E6;

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      Wing_MinMaxThickness = min(Wing_MinMaxThickness, MaxThickness[iPlane]);
      Wing_MaxChord = max(Wing_MaxChord, Chord[iPlane]);
      Wing_MinToC = min(Wing_MinToC, ToC[iPlane]);
      Wing_MaxTwist = max(Wing_MaxTwist, fabs(Twist[iPlane]));
      Wing_MaxCurvature = max(Wing_MaxCurvature, Curvature[iPlane]);
      Wing_MaxDihedral = max(Wing_MaxDihedral, fabs(Dihedral[iPlane]));
    }

  }

  /*--- Free memory for the section cuts ---*/

  delete [] Xcoord_Airfoil;
  delete [] Ycoord_Airfoil;
  delete [] Zcoord_Airfoil;
  delete [] Variable_Airfoil;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] LeadingEdge[iPlane];
  delete [] LeadingEdge;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] TrailingEdge[iPlane];
  delete [] TrailingEdge;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] Plane_P0[iPlane];
  delete [] Plane_P0;
  
  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] Plane_Normal[iPlane];
  delete [] Plane_Normal;
  
  delete [] Area;
  delete [] MaxThickness;
  delete [] Chord;
  delete [] ToC;
  delete [] Twist;
  delete [] Curvature;
  delete [] Dihedral;

}

CMultiGridGeometry::CMultiGridGeometry(CGeometry ***geometry, CConfig **config_container, unsigned short iMesh, unsigned short iZone) : CGeometry() {
  
  /*--- CGeometry & CConfig pointers to the fine grid level for clarity. We may
   need access to the other zones in the mesh for zone boundaries. ---*/
  
  CGeometry *fine_grid = geometry[iZone][iMesh-1];
  CConfig *config = config_container[iZone];
  
  /*--- Local variables ---*/
  
  unsigned long iPoint, Index_CoarseCV, CVPoint, iElem, iVertex, jPoint, iteration, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, iParent, jVertex, *Buffer_Receive_Parent = NULL, *Buffer_Send_Parent = NULL, *Buffer_Receive_Children = NULL, *Buffer_Send_Children = NULL, *Parent_Remote = NULL, *Children_Remote = NULL, *Parent_Local = NULL, *Children_Local = NULL, Local_nPointCoarse, Local_nPointFine, Global_nPointCoarse, Global_nPointFine;;
  short marker_seed;
  bool agglomerate_seed = true;
  unsigned short nChildren, iNode, counter, iMarker, jMarker, priority, MarkerS, MarkerR, *nChildren_MPI;
  vector<unsigned long> Suitable_Indirect_Neighbors, Aux_Parent;
  vector<unsigned long>::iterator it;
  int rank;

  unsigned short nMarker_Max = config->GetnMarker_Max();

  unsigned short *copy_marker = new unsigned short [nMarker_Max];
  
#ifndef HAVE_MPI
  rank = MASTER_NODE;
#else
  int send_to, receive_from;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;
#endif
  
  nDim = fine_grid->GetnDim(); // Write the number of dimensions of the coarse grid.
  
  /*--- Create a queue system to deo the agglomeration
   1st) More than two markers ---> Vertices (never agglomerate)
   2nd) Two markers ---> Edges (agglomerate if same BC, never agglomerate if different BC)
   3rd) One marker ---> Surface (always agglomarate)
   4th) No marker ---> Internal Volume (always agglomarate) ---*/
  
  /*--- Set a marker to indicate indirect agglomeration ---*/
  
  if (iMesh == MESH_1) {
    
    for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++)
      fine_grid->node[iPoint]->SetAgglomerate_Indirect(false);
    
    for (iElem = 0; iElem < fine_grid->GetnElem(); iElem++) {
      if ((fine_grid->elem[iElem]->GetVTK_Type() == HEXAHEDRON) ||
          (fine_grid->elem[iElem]->GetVTK_Type() == QUADRILATERAL)) {
        for (iNode = 0; iNode < fine_grid->elem[iElem]->GetnNodes(); iNode++) {
          iPoint = fine_grid->elem[iElem]->GetNode(iNode);
          fine_grid->node[iPoint]->SetAgglomerate_Indirect(true);
        }
      }
    }
    
  }
  
  /*--- Create the coarse grid structure using as baseline the fine grid ---*/
  
  CMultiGridQueue MGQueue_InnerCV(fine_grid->GetnPoint());
 
  nPointNode = fine_grid->GetnPoint(); 
  node = new CPoint*[fine_grid->GetnPoint()];
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    
    /*--- Create node structure ---*/
    
    node[iPoint] = new CPoint(nDim, iPoint, config);
    
    /*--- Set the indirect agglomeration to false ---*/
    
    node[iPoint]->SetAgglomerate_Indirect(false);
  }
  
  Index_CoarseCV = 0;
  
  /*--- The first step is the boundary agglomeration. ---*/
  
  for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++) {
    
    for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      
      /*--- If the element has not being previously agglomerated and it belongs
       to the physical domain, then the agglomeration is studied ---*/
      
      if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
          (fine_grid->node[iPoint]->GetDomain()) &&
          (GeometricalCheck(iPoint, fine_grid, config))) {
        
        nChildren = 1;
        
        /*--- We set an index for the parent control volume ---*/
        
        fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
        
        /*--- We add the seed point (child) to the parent control volume ---*/
        
        node[Index_CoarseCV]->SetChildren_CV(0, iPoint);
        agglomerate_seed = true; counter = 0; marker_seed = iMarker;
        
        /*--- For a particular point in the fine grid we save all the markers
         that are in that point ---*/
        
        for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
          if (fine_grid->node[iPoint]->GetVertex(jMarker) != -1) {
            copy_marker[counter] = jMarker;
            counter++;
          }
        
        /*--- To aglomerate a vertex it must have only one physical bc!!
         This can be improved. If there is only a marker, it is a good
         candidate for agglomeration ---*/
        
        if (counter == 1) agglomerate_seed = true;
        
        /*--- If there are two markers, we will aglomerate if one of the
         marker is SEND_RECEIVE ---*/
        
        if (counter == 2) {
          if ((config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE) ||
              (config->GetMarker_All_KindBC(copy_marker[1]) == SEND_RECEIVE)) agglomerate_seed = true;
          else agglomerate_seed = false;
        }
        
        /*--- If there are more than 2 markers, the aglomeration will be discarted ---*/
        
        if (counter > 2) agglomerate_seed = false;
        
        /*--- If the seed can be agglomerated, we try to agglomerate more points ---*/
        
        if (agglomerate_seed) {
          
          /*--- Now we do a sweep over all the nodes that surround the seed point ---*/
          
          for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
            
            CVPoint = fine_grid->node[iPoint]->GetPoint(iNode);
            
            /*--- The new point can be agglomerated ---*/
            
            if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config)) {
              
              /*--- We set the value of the parent ---*/
              
              fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
              
              /*--- We set the value of the child ---*/
              
              node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
              nChildren++;
            }
            
          }
          
          Suitable_Indirect_Neighbors.clear();
          
          if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
            SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);
          
          /*--- Now we do a sweep over all the indirect nodes that can be added ---*/
          
          for (iNode = 0; iNode <	Suitable_Indirect_Neighbors.size(); iNode ++) {
            
            CVPoint = Suitable_Indirect_Neighbors[iNode];
            
            /*--- The new point can be agglomerated ---*/
            
            if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config)) {
              
              /*--- We set the value of the parent ---*/
              
              fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
              
              /*--- We set the indirect agglomeration information ---*/
              
              if (fine_grid->node[CVPoint]->GetAgglomerate_Indirect())
                node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
              
              /*--- We set the value of the child ---*/
              
              node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
              nChildren++;
            }
          }
          
          
        }
        
        /*--- Update the number of child of the control volume ---*/
        
        node[Index_CoarseCV]->SetnChildren_CV(nChildren);
        Index_CoarseCV++;
      }
    }
  }
  
  /*--- Agglomerate all the nodes that have more than one physical boundary condition,
   Maybe here we can add the posibility of merging the vertex that have the same number,
   and kind  of markers---*/
  
  for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++)
    for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
          (fine_grid->node[iPoint]->GetDomain())) {
        fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
        node[Index_CoarseCV]->SetChildren_CV(0, iPoint);
        node[Index_CoarseCV]->SetnChildren_CV(1);
        Index_CoarseCV++;
      }
    }
  
  /*--- Update the queue with the results from the boundary agglomeration ---*/
  
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    
    /*--- The CV has been agglomerated, remove form the list ---*/
    
    if (fine_grid->node[iPoint]->GetAgglomerate() == true) {
      
      MGQueue_InnerCV.RemoveCV(iPoint);
      
    }
    
    else {
      
      /*--- Count the number of agglomerated neighbors, and modify the queue ---*/
      
      priority = 0;
      for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
        jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
        if (fine_grid->node[jPoint]->GetAgglomerate() == true) priority++;
      }
      MGQueue_InnerCV.MoveCV(iPoint, priority);
    }
  }
  
  /*--- Agglomerate the domain nodes ---*/
  
  iteration = 0;
  while (!MGQueue_InnerCV.EmptyQueue() && (iteration < fine_grid->GetnPoint())) {
    
    iPoint = MGQueue_InnerCV.NextCV();
    iteration ++;
    
    /*--- If the element has not being previously agglomerated, belongs to the physical domain,
     and satisfies several geometrical criteria then the seed CV is acepted for agglomeration ---*/
    
    if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
        (fine_grid->node[iPoint]->GetDomain()) &&
        (GeometricalCheck(iPoint, fine_grid, config))) {
      
      nChildren = 1;
      
      /*--- We set an index for the parent control volume ---*/
      
      fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
      
      /*--- We add the seed point (child) to the parent control volume ---*/
      
      node[Index_CoarseCV]->SetChildren_CV(0, iPoint);
      
      /*--- Update the queue with the seed point (remove the seed and
       increase the priority of the neighbors) ---*/
      
      MGQueue_InnerCV.Update(iPoint, fine_grid);
      
      /*--- Now we do a sweep over all the nodes that surround the seed point ---*/
      
      for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
        
        CVPoint = fine_grid->node[iPoint]->GetPoint(iNode);
        
        /*--- Determine if the CVPoint can be agglomerated ---*/
        
        if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
            (fine_grid->node[CVPoint]->GetDomain()) &&
            (GeometricalCheck(CVPoint, fine_grid, config))) {
          
          /*--- We set the value of the parent ---*/
          
          fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
          
          /*--- We set the value of the child ---*/
          
          node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
          nChildren++;
          
          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of the neighbors) ---*/
          
          MGQueue_InnerCV.Update(CVPoint, fine_grid);
          
        }
        
      }
      
      /*--- Subrotuine to identify the indirect neighbors ---*/
      
      Suitable_Indirect_Neighbors.clear();
      if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
        SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);
      
      /*--- Now we do a sweep over all the indirect nodes that can be added ---*/
      
      for (iNode = 0; iNode <	Suitable_Indirect_Neighbors.size(); iNode ++) {
        
        CVPoint = Suitable_Indirect_Neighbors[iNode];
        
        /*--- The new point can be agglomerated ---*/
        
        if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
            (fine_grid->node[CVPoint]->GetDomain())) {
          
          /*--- We set the value of the parent ---*/
          
          fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
          
          /*--- We set the indirect agglomeration information ---*/
          
          if (fine_grid->node[CVPoint]->GetAgglomerate_Indirect())
            node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
          
          /*--- We set the value of the child ---*/
          
          node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
          nChildren++;
          
          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of the neighbors) ---*/
          
          MGQueue_InnerCV.Update(CVPoint, fine_grid);
          
        }
      }
      
      /*--- Update the number of control of childrens ---*/
      
      node[Index_CoarseCV]->SetnChildren_CV(nChildren);
      Index_CoarseCV++;
    }
    else {
      
      /*--- The seed point can not be agglomerated because of size, domain, streching, etc.
       move the point to the lowest priority ---*/
      
      MGQueue_InnerCV.MoveCV(iPoint, -1);
    }
    
  }
  
  /*--- Add all the elements that have not being agglomerated, in the previous stage ---*/
  
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    if ((fine_grid->node[iPoint]->GetAgglomerate() == false) && (fine_grid->node[iPoint]->GetDomain())) {
      
      nChildren = 1;
      fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
      if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
        node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
      node[Index_CoarseCV]->SetChildren_CV(0, iPoint);
      node[Index_CoarseCV]->SetnChildren_CV(nChildren);
      Index_CoarseCV++;
      
    }
  }
  
  nPointDomain = Index_CoarseCV;
  
  /*--- Check that there are no hanging nodes ---*/
  
  unsigned long iFinePoint, iFinePoint_Neighbor, iCoarsePoint, iCoarsePoint_Complete;
  unsigned short iChildren;
  
  /*--- Find the point surrounding a point ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPointDomain; iCoarsePoint ++) {
    for (iChildren = 0; iChildren <  node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
        iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
        iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
        if (iParent != iCoarsePoint) node[iCoarsePoint]->SetPoint(iParent);
      }
    }
  }
  
  /*--- Detect isolated points and merge them with its correct neighbor ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPointDomain; iCoarsePoint ++) {
    
    if (node[iCoarsePoint]->GetnPoint() == 1) {
      
      /*--- Find the neighbor of the isolated point. This neighbor is the right control volume ---*/
      
      iCoarsePoint_Complete = node[iCoarsePoint]->GetPoint(0);
      
      /*--- Add the children to the connected control volume (and modify it parent indexing).
       Identify the child CV from the finest grid and added to the correct control volume.
       Set the parent CV of iFinePoint. Instead of using the original
       (iCoarsePoint) one use the new one (iCoarsePoint_Complete) ---*/
      
      nChildren = node[iCoarsePoint_Complete]->GetnChildren_CV();
      
      for (iChildren = 0; iChildren <  node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        node[iCoarsePoint_Complete]->SetChildren_CV(nChildren, iFinePoint);
        nChildren++;
        fine_grid->node[iFinePoint]->SetParent_CV(iCoarsePoint_Complete);
      }
      
      /*--- Update the number of children control volumes ---*/
      
      node[iCoarsePoint_Complete]->SetnChildren_CV(nChildren);
      node[iCoarsePoint]->SetnChildren_CV(0);
      
    }
  }
  
  //  unsigned long iPointFree = nPointDomain-1;
  //  iCoarsePoint = 0;
  //
  //  do {
  //
  //    if (node[iCoarsePoint]->GetnChildren_CV() == 0) {
  //
  //      while (node[iPointFree]->GetnChildren_CV() == 0) {
  //        Index_CoarseCV--;
  //        iPointFree--;
  //      }
  //
  //      nChildren = node[iPointFree]->GetnChildren_CV();
  //      for (iChildren = 0; iChildren <  nChildren; iChildren ++) {
  //        iFinePoint = node[iPointFree]->GetChildren_CV(iChildren);
  //        node[iCoarsePoint]->SetChildren_CV(iChildren, iFinePoint);
  //        fine_grid->node[iFinePoint]->SetParent_CV(iCoarsePoint);
  //      }
  //      node[iCoarsePoint]->SetnChildren_CV(nChildren);
  //      node[iPointFree]->SetnChildren_CV(0);
  //
  //      Index_CoarseCV--;
  //      iPointFree--;
  //
  //    }
  //
  //    iCoarsePoint++;
  //
  //  } while ((iCoarsePoint-1) < Index_CoarseCV);
  //
  //  nPointDomain = Index_CoarseCV;
  
  /*--- Reset the point surrounding a point ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPointDomain; iCoarsePoint ++) {
    node[iCoarsePoint]->ResetPoint();
  }
  
  /*--- Dealing with MPI parallelization, the objective is that the received nodes must be agglomerated
   in the same way as the donor nodes. Send the node agglomeration information of the donor
   (parent and children), Sending only occurs with MPI ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = fine_grid->nVertex[MarkerS];   nVertexR = fine_grid->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;               nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive_Children = new unsigned long [nBufferR_Vector];
      Buffer_Send_Children = new unsigned long [nBufferS_Vector];
      
      Buffer_Receive_Parent = new unsigned long [nBufferR_Vector];
      Buffer_Send_Parent = new unsigned long [nBufferS_Vector];
      
      /*--- Copy the information that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = fine_grid->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Children[iVertex] = iPoint;
        Buffer_Send_Parent[iVertex] = fine_grid->node[iPoint]->GetParent_CV();
      }
      
#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Children, nBufferS_Vector, MPI_UNSIGNED_LONG, send_to,0,
                   Buffer_Receive_Children, nBufferR_Vector, MPI_UNSIGNED_LONG, receive_from,0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_Parent, nBufferS_Vector, MPI_UNSIGNED_LONG, send_to,1,
                   Buffer_Receive_Parent, nBufferR_Vector, MPI_UNSIGNED_LONG, receive_from,1, MPI_COMM_WORLD, &status);
#else
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Children[iVertex] = Buffer_Send_Children[iVertex];
        Buffer_Receive_Parent[iVertex] = Buffer_Send_Parent[iVertex];
      }
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send_Children;
      delete [] Buffer_Send_Parent;
      
      /*--- Create a list of the parent nodes without repeated parents ---*/
      
      Aux_Parent.clear();
      for (iVertex = 0; iVertex < nVertexR; iVertex++)
        Aux_Parent.push_back (Buffer_Receive_Parent[iVertex]);
      
      sort(Aux_Parent.begin(), Aux_Parent.end());
      it = unique(Aux_Parent.begin(), Aux_Parent.end());
      Aux_Parent.resize(it - Aux_Parent.begin());
      
      /*--- Allocate some structures ---*/
      
      Parent_Remote = new unsigned long[nVertexR];
      Children_Remote = new unsigned long[nVertexR];
      Parent_Local = new unsigned long[nVertexR];
      Children_Local = new unsigned long[nVertexR];
      
      /*--- Create the local vector and remote for the parents and the children ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        Parent_Remote[iVertex] = Buffer_Receive_Parent[iVertex];
        
        /*--- We use the same sorting as in the donor domain ---*/
        
        for (jVertex = 0; jVertex < Aux_Parent.size(); jVertex++) {
          if (Parent_Remote[iVertex] == Aux_Parent[jVertex]) {
            Parent_Local[iVertex] = jVertex + Index_CoarseCV;
            break;
          }
        }
        
        Children_Remote[iVertex] = Buffer_Receive_Children[iVertex];
        Children_Local[iVertex] = fine_grid->vertex[MarkerR][iVertex]->GetNode();
        
      }
      
      Index_CoarseCV += Aux_Parent.size();
      
      nChildren_MPI = new unsigned short [Index_CoarseCV];
      for (iParent = 0; iParent < Index_CoarseCV; iParent++)
        nChildren_MPI[iParent] = 0;
      
      /*--- Create the final structure ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Be careful, it is possible that a node change the agglomeration configuration, the priority
         is always, when receive the information ---*/
        
        fine_grid->node[Children_Local[iVertex]]->SetParent_CV(Parent_Local[iVertex]);
        node[Parent_Local[iVertex]]->SetChildren_CV(nChildren_MPI[Parent_Local[iVertex]], Children_Local[iVertex]);
        nChildren_MPI[Parent_Local[iVertex]]++;
        node[Parent_Local[iVertex]]->SetnChildren_CV(nChildren_MPI[Parent_Local[iVertex]]);
        node[Parent_Local[iVertex]]->SetDomain(false);
        
      }
      
      /*--- Deallocate auxiliar structures ---*/
      
      delete[] nChildren_MPI;
      delete[] Parent_Remote;
      delete[] Children_Remote;
      delete[] Parent_Local;
      delete[] Children_Local;
      
      /*--- Deallocate receive buffer ---*/
      
      delete [] Buffer_Receive_Children;
      delete [] Buffer_Receive_Parent;
      
    }
    
  }
  
  /*--- Update the number of points after the MPI agglomeration ---*/
  
  nPoint = Index_CoarseCV;
  
  /*--- Console output with the summary of the agglomeration ---*/
  
  Local_nPointCoarse = nPoint;
  Local_nPointFine = fine_grid->GetnPoint();
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&Local_nPointCoarse, &Global_nPointCoarse, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nPointFine, &Global_nPointFine, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nPointCoarse = Local_nPointCoarse;
  Global_nPointFine = Local_nPointFine;
#endif
  
  su2double Coeff = 1.0, CFL = 0.0, factor = 1.5;
  
  if (iMesh != MESH_0) {
    if (nDim == 2) Coeff = pow(su2double(Global_nPointFine)/su2double(Global_nPointCoarse), 1./2.);
    if (nDim == 3) Coeff = pow(su2double(Global_nPointFine)/su2double(Global_nPointCoarse), 1./3.);
    CFL = factor*config->GetCFL(iMesh-1)/Coeff;
    config->SetCFL(iMesh, CFL);
  }
  
  su2double ratio = su2double(Global_nPointFine)/su2double(Global_nPointCoarse);
  
  if (((nDim == 2) && (ratio < 2.5)) ||
      ((nDim == 3) && (ratio < 2.5))) {
    config->SetMGLevels(iMesh-1);
  }
  else {
    if (rank == MASTER_NODE) {
      if (iMesh == 1) cout <<"MG level: "<< iMesh-1 <<" -> CVs: " << Global_nPointFine << ". Agglomeration rate 1/1.00. CFL "<< config->GetCFL(iMesh-1) <<"." << endl;
      cout <<"MG level: "<< iMesh <<" -> CVs: " << Global_nPointCoarse << ". Agglomeration rate 1/" << ratio <<". CFL "<< CFL <<"." << endl;
    }
  }
 
  delete [] copy_marker;
  
}


CMultiGridGeometry::~CMultiGridGeometry(void) {
  
}

bool CMultiGridGeometry::SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, CGeometry *fine_grid, CConfig *config) {
  
  bool agglomerate_CV = false;
  unsigned short counter, jMarker;
  
  unsigned short nMarker_Max = config->GetnMarker_Max();
  
  unsigned short *copy_marker = new unsigned short [nMarker_Max];

  /*--- Basic condition, the element has not being previously agglomerated, it belongs to the domain,
   and has passed some basic geometrical check ---*/
  
  if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
      (fine_grid->node[CVPoint]->GetDomain()) &&
      (GeometricalCheck(CVPoint, fine_grid, config))) {
    
    /*--- If the element belong to the boundary, we must be careful ---*/
    
    if (fine_grid->node[CVPoint]->GetBoundary()) {
      
      /*--- Identify the markers of the vertex that we want to agglomerate ---*/
      
      counter = 0;
      for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
        if (fine_grid->node[CVPoint]->GetVertex(jMarker) != -1) {
          copy_marker[counter] = jMarker;
          counter++;
        }
      
      /*--- The basic condition is that the aglomerated vertex must have the same physical marker,
       but eventually a send-receive condition ---*/
      
      /*--- Only one marker in the vertex that is going to be aglomerated ---*/
      
      if (counter == 1) {
        
        /*--- We agglomerate if there is only a marker and is the same marker as the seed marker ---*/
        
        if (copy_marker[0] == marker_seed)
          agglomerate_CV = true;
        
        /*--- If there is only a marker, but the marker is the SEND_RECEIVE ---*/
        
        if (config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE)
          agglomerate_CV = true;
        
      }
      
      /*--- If there are two markers in the vertex that is going to be aglomerated ---*/
      
      if (counter == 2) {
        
        /*--- First we verify that the seed is a physical boundary ---*/
        
        if (config->GetMarker_All_KindBC(marker_seed) != SEND_RECEIVE) {
          
          /*--- Then we check that one of the marker is equal to the seed marker, and the other is send/receive ---*/
          
          if (((copy_marker[0] == marker_seed) && (config->GetMarker_All_KindBC(copy_marker[1]) == SEND_RECEIVE)) ||
              ((config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE) && (copy_marker[1] == marker_seed)))
            agglomerate_CV = true;
        }
        
      }
      
    }
    
    /*--- If the element belong to the domain, it is allways aglomerated ---*/
    
    else { agglomerate_CV = true; }
    
  }
  
  delete [] copy_marker;

  return agglomerate_CV;

}


bool CMultiGridGeometry::GeometricalCheck(unsigned long iPoint, CGeometry *fine_grid, CConfig *config) {
  
  su2double max_dimension = 1.2;
  
  /*--- Evaluate the total size of the element ---*/
  
  bool Volume = true;
  su2double ratio = pow(fine_grid->node[iPoint]->GetVolume(), 1.0/su2double(nDim))*max_dimension;
  su2double limit = pow(config->GetDomainVolume(), 1.0/su2double(nDim));
  if ( ratio > limit ) Volume = false;
  
  /*--- Evaluate the stretching of the element ---*/
  
  bool Stretching = true;
  
  /*	unsigned short iNode, iDim;
   unsigned long jPoint;
   su2double *Coord_i = fine_grid->node[iPoint]->GetCoord();
   su2double max_dist = 0.0 ; su2double min_dist = 1E20;
   for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
   jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
   su2double *Coord_j = fine_grid->node[jPoint]->GetCoord();
   su2double distance = 0.0;
   for (iDim = 0; iDim < nDim; iDim++)
   distance += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
   distance = sqrt(distance);
   max_dist = max(distance, max_dist);
   min_dist = min(distance, min_dist);
   }
   if ( max_dist/min_dist > 100.0 ) Stretching = false;*/
  
  return (Stretching && Volume);
  
}

void CMultiGridGeometry::SetSuitableNeighbors(vector<unsigned long> *Suitable_Indirect_Neighbors, unsigned long iPoint,
                                              unsigned long Index_CoarseCV, CGeometry *fine_grid) {
  
  unsigned long jPoint, kPoint, lPoint;
  unsigned short iNode, jNode, iNeighbor, jNeighbor, kNode;
  bool SecondNeighborSeed, ThirdNeighborSeed;
  vector<unsigned long>::iterator it;
  
  /*--- Create a list with the first neighbors, including the seed ---*/
  
  vector<unsigned long> First_Neighbor_Points;
  First_Neighbor_Points.push_back(iPoint);
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    First_Neighbor_Points.push_back(jPoint);
  }
  
  /*--- Create a list with the second neighbors, without first, and seed neighbors ---*/
  
  vector<unsigned long> Second_Neighbor_Points, Second_Origin_Points, Suitable_Second_Neighbors;
  
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    
    for (jNode = 0; jNode <	fine_grid->node[jPoint]->GetnPoint(); jNode ++) {
      kPoint = fine_grid->node[jPoint]->GetPoint(jNode);
      
      /*--- Check that the second neighbor do not belong to the first neighbor or the seed ---*/
      
      SecondNeighborSeed = true;
      for (iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
        if (kPoint == First_Neighbor_Points[iNeighbor]) {
          SecondNeighborSeed = false; break;
        }
      
      if (SecondNeighborSeed) {
        Second_Neighbor_Points.push_back(kPoint);
        Second_Origin_Points.push_back(jPoint);
      }
      
    }
  }
  
  /*---  Identify those second neighbors that are repeated (candidate to be added) ---*/
  
  for (iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
    
    for (jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++)
      
    /*--- Repeated second neighbor with different origin ---*/
      
      if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) &&
          (Second_Origin_Points[iNeighbor] != Second_Origin_Points[jNeighbor]) &&
          (iNeighbor < jNeighbor)) {
        
        Suitable_Indirect_Neighbors->push_back(Second_Neighbor_Points[iNeighbor]);
        
        /*--- Create alist with the suitable second neighbor, that we will use
         to compute the third neighbors --*/
        
        Suitable_Second_Neighbors.push_back(Second_Neighbor_Points[iNeighbor]);
        
      }
  
  
  /*--- Remove repeated from the suitable second neighbors ---*/
  
  sort(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  it = unique(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  Suitable_Second_Neighbors.resize(it - Suitable_Second_Neighbors.begin());
  
  /*--- Remove repeated from first neighbors ---*/
  
  sort(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
  it = unique(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
  First_Neighbor_Points.resize(it - First_Neighbor_Points.begin());
  
  /*--- Create a list with the third neighbors, without first, second, and seed neighbors ---*/
  
  vector<unsigned long> Third_Neighbor_Points, Third_Origin_Points;
  
  for (jNode = 0; jNode <	Suitable_Second_Neighbors.size(); jNode ++) {
    kPoint = Suitable_Second_Neighbors[jNode];
    
    for (kNode = 0; kNode <	fine_grid->node[kPoint]->GetnPoint(); kNode ++) {
      lPoint = fine_grid->node[kPoint]->GetPoint(kNode);
      
      /*--- Check that the third neighbor do not belong to the first neighbors or the seed ---*/
      
      ThirdNeighborSeed = true;
      
      for (iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
        if (lPoint == First_Neighbor_Points[iNeighbor]) {
          ThirdNeighborSeed = false;
          break;
        }
      
      /*--- Check that the third neighbor do not belong to the second neighbors ---*/
      
      for (iNeighbor = 0; iNeighbor <	Suitable_Second_Neighbors.size(); iNeighbor ++)
        if (lPoint == Suitable_Second_Neighbors[iNeighbor]) {
          ThirdNeighborSeed = false;
          break;
        }
      
      if (ThirdNeighborSeed) {
        Third_Neighbor_Points.push_back(lPoint);
        Third_Origin_Points.push_back(kPoint);
      }
      
    }
  }
  
  /*---  Identify those third neighbors that are repeated (candidate to be added) ---*/
  
  for (iNeighbor = 0; iNeighbor <	Third_Neighbor_Points.size(); iNeighbor ++)
    for (jNeighbor = 0; jNeighbor <	Third_Neighbor_Points.size(); jNeighbor ++)
      
    /*--- Repeated second neighbor with different origin ---*/
      
      if ((Third_Neighbor_Points[iNeighbor] == Third_Neighbor_Points[jNeighbor]) &&
          (Third_Origin_Points[iNeighbor] != Third_Origin_Points[jNeighbor]) &&
          (iNeighbor < jNeighbor)) {
        
        Suitable_Indirect_Neighbors->push_back(Third_Neighbor_Points[iNeighbor]);
        
      }
  
  /*--- Remove repeated from Suitable Indirect Neighbors List ---*/
  
  sort(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
  it = unique(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
  Suitable_Indirect_Neighbors->resize(it - Suitable_Indirect_Neighbors->begin());
  
}

void CMultiGridGeometry::SetPoint_Connectivity(CGeometry *fine_grid) {
  
  unsigned long iFinePoint, iFinePoint_Neighbor, iParent, iCoarsePoint;
  unsigned short iChildren, iNode;
  
  /*--- Set the point surrounding a point ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    for (iChildren = 0; iChildren <  node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
        iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
        iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
        if (iParent != iCoarsePoint) node[iCoarsePoint]->SetPoint(iParent);
      }
    }
  }
  
  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    node[iCoarsePoint]->SetnNeighbor(node[iCoarsePoint]->GetnPoint());
  
}

void CMultiGridGeometry::SetVertex(CGeometry *fine_grid, CConfig *config) {
  unsigned long  iVertex, iFinePoint, iCoarsePoint;
  unsigned short iMarker, iMarker_Tag, iChildren;
  
  nMarker = fine_grid->GetnMarker();
  unsigned short nMarker_Max = config->GetnMarker_Max();

  /*--- If any children node belong to the boundary then the entire control
   volume will belong to the boundary ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      if (fine_grid->node[iFinePoint]->GetBoundary()) {
        node[iCoarsePoint]->SetBoundary(nMarker);
        break;
      }
    }
  
  vertex = new CVertex**[nMarker];
  nVertex = new unsigned long [nMarker];
  
  Tag_to_Marker = new string [nMarker_Max];
  for (iMarker_Tag = 0; iMarker_Tag < nMarker_Max; iMarker_Tag++)
    Tag_to_Marker[iMarker_Tag] = fine_grid->GetMarker_Tag(iMarker_Tag);
  
  /*--- Compute the number of vertices to do the dimensionalization ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;
  
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    if (node[iCoarsePoint]->GetBoundary()) {
      for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        for (iMarker = 0; iMarker < nMarker; iMarker ++) {
          if ((fine_grid->node[iFinePoint]->GetVertex(iMarker) != -1) && (node[iCoarsePoint]->GetVertex(iMarker) == -1)) {
            iVertex = nVertex[iMarker];
            node[iCoarsePoint]->SetVertex(iVertex, iMarker);
            nVertex[iMarker]++;
          }
        }
      }
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [fine_grid->GetnVertex(iMarker)+1];
    nVertex[iMarker] = 0;
  }
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    if (node[iCoarsePoint]->GetBoundary())
      for (iMarker = 0; iMarker < nMarker; iMarker ++)
        node[iCoarsePoint]->SetVertex(-1, iMarker);
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    if (node[iCoarsePoint]->GetBoundary())
      for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker ++) {
          if ((fine_grid->node[iFinePoint]->GetVertex(iMarker) != -1) && (node[iCoarsePoint]->GetVertex(iMarker) == -1)) {
            iVertex = nVertex[iMarker];
            vertex[iMarker][iVertex] = new CVertex(iCoarsePoint, nDim);
            node[iCoarsePoint]->SetVertex(iVertex, iMarker);
            
            /*--- Set the transformation to apply ---*/
            unsigned long ChildVertex = fine_grid->node[iFinePoint]->GetVertex(iMarker);
            unsigned short RotationKind = fine_grid->vertex[iMarker][ChildVertex]->GetRotation_Type();
            vertex[iMarker][iVertex]->SetRotation_Type(RotationKind);
            nVertex[iMarker]++;
          }
        }
      }
}

void CMultiGridGeometry::MatchNearField(CConfig *config) {
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  int iProcessor;
  
#ifndef HAVE_MPI
  iProcessor = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &iProcessor);
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, node[iPoint]->GetGlobalIndex(), iVertex, iMarker, iProcessor);
        }
      }
    }
  }
  
}

void CMultiGridGeometry::MatchActuator_Disk(CConfig *config) {
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  int iProcessor;
  
#ifndef HAVE_MPI
  iProcessor = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &iProcessor);
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
        (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, node[iPoint]->GetGlobalIndex(), iVertex, iMarker, iProcessor);
        }
      }
    }
  }
  
}

void CMultiGridGeometry::MatchInterface(CConfig *config) {
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  int iProcessor;
  
#ifndef HAVE_MPI
  iProcessor = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &iProcessor);
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, node[iPoint]->GetGlobalIndex(), iVertex, iMarker, iProcessor);
        }
      }
    }
  }
  
}

void CMultiGridGeometry::SetControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {
  
  unsigned long iFinePoint, iFinePoint_Neighbor, iCoarsePoint, iEdge, iParent;
  long FineEdge, CoarseEdge;
  unsigned short iChildren, iNode, iDim;
  bool change_face_orientation;
  su2double *Normal, Coarse_Volume, Area, *NormalFace = NULL;
  Normal = new su2double [nDim];
  
  /*--- Compute the area of the coarse volume ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    node[iCoarsePoint]->SetVolume(0.0);
    Coarse_Volume = 0.0;
    for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      Coarse_Volume += fine_grid->node[iFinePoint]->GetVolume();
    }
    node[iCoarsePoint]->SetVolume(Coarse_Volume);
  }
  
  /*--- Update or not the values of faces at the edge ---*/
  if (action != ALLOCATE) {
    for (iEdge=0; iEdge < nEdge; iEdge++)
      edge[iEdge]->SetZeroValues();
  }
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      
      for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
        iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
        iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
        if ((iParent != iCoarsePoint) && (iParent < iCoarsePoint)) {
          
          FineEdge = fine_grid->FindEdge(iFinePoint, iFinePoint_Neighbor);
          
          change_face_orientation = false;
          if (iFinePoint < iFinePoint_Neighbor) change_face_orientation = true;
          
          CoarseEdge = FindEdge(iParent, iCoarsePoint);
          
          fine_grid->edge[FineEdge]->GetNormal(Normal);
          
          if (change_face_orientation) {
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            edge[CoarseEdge]->AddNormal(Normal);
          }
          else {
            edge[CoarseEdge]->AddNormal(Normal);
          }
        }
      }
    }
  delete[] Normal;
  
  /*--- Check if there is a normal with null area ---*/
  
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    NormalFace = edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
    Area = sqrt(Area);
    if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
  }
  
}

void CMultiGridGeometry::SetBoundControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {
  unsigned long iCoarsePoint, iFinePoint, FineVertex, iVertex;
  unsigned short iMarker, iChildren, iDim;
  su2double *Normal, Area, *NormalFace = NULL;
  
  Normal = new su2double [nDim];
  
  if (action != ALLOCATE) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        vertex[iMarker][iVertex]->SetZeroValues();
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iCoarsePoint = vertex[iMarker][iVertex]->GetNode();
      for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        if (fine_grid->node[iFinePoint]->GetVertex(iMarker)!=-1) {
          FineVertex = fine_grid->node[iFinePoint]->GetVertex(iMarker);
          fine_grid->vertex[iMarker][FineVertex]->GetNormal(Normal);
          vertex[iMarker][iVertex]->AddNormal(Normal);
        }
      }
    }
  
  delete[] Normal;
  
  /*--- Check if there is a normal with null area ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
      Area = sqrt(Area);
      if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }
  
}

void CMultiGridGeometry::SetCoord(CGeometry *geometry) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iChildren, iDim;
  su2double Area_Parent, Area_Children;
  su2double *Coordinates_Fine, *Coordinates;
  Coordinates = new su2double[nDim];
  
  for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
    Area_Parent = node[Point_Coarse]->GetVolume();
    for (iDim = 0; iDim < nDim; iDim++) Coordinates[iDim] = 0.0;
    for (iChildren = 0; iChildren < node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geometry->node[Point_Fine]->GetVolume();
      Coordinates_Fine = geometry->node[Point_Fine]->GetCoord();
      for (iDim = 0; iDim < nDim; iDim++)
        Coordinates[iDim] += Coordinates_Fine[iDim]*Area_Children/Area_Parent;
    }
    for (iDim = 0; iDim < nDim; iDim++)
      node[Point_Coarse]->SetCoord(iDim, Coordinates[iDim]);
  }
  delete[] Coordinates;
}

void CMultiGridGeometry::SetRotationalVelocity(CConfig *config, unsigned short val_iZone) {
  
  unsigned long iPoint_Coarse;
  su2double *RotVel, Distance[3] = {0.0,0.0,0.0}, *Coord;
  su2double Center[3] = {0.0,0.0,0.0}, Omega[3] = {0.0,0.0,0.0}, L_Ref;
  RotVel = new su2double [3];
  
  /*--- Center of rotation & angular velocity vector from config. ---*/
  
  Center[0] = config->GetMotion_Origin_X(val_iZone);
  Center[1] = config->GetMotion_Origin_Y(val_iZone);
  Center[2] = config->GetMotion_Origin_Z(val_iZone);
  Omega[0]  = config->GetRotation_Rate_X(val_iZone)/config->GetOmega_Ref();
  Omega[1]  = config->GetRotation_Rate_Y(val_iZone)/config->GetOmega_Ref();
  Omega[2]  = config->GetRotation_Rate_Z(val_iZone)/config->GetOmega_Ref();
  L_Ref     = config->GetLength_Ref();
  
  /*--- Loop over all nodes and set the rotational velocity. ---*/
  
  for (iPoint_Coarse = 0; iPoint_Coarse < GetnPoint(); iPoint_Coarse++) {
    
    /*--- Get the coordinates of the current node ---*/
    
    Coord = node[iPoint_Coarse]->GetCoord();
    
    /*--- Calculate the non-dim. distance from the rotation center ---*/
    
    Distance[0] = (Coord[0]-Center[0])/L_Ref;
    Distance[1] = (Coord[1]-Center[1])/L_Ref;
    Distance[2] = (Coord[2]-Center[2])/L_Ref;
    
    /*--- Calculate the angular velocity as omega X r ---*/
    
    RotVel[0] = Omega[1]*(Distance[2]) - Omega[2]*(Distance[1]);
    RotVel[1] = Omega[2]*(Distance[0]) - Omega[0]*(Distance[2]);
    RotVel[2] = Omega[0]*(Distance[1]) - Omega[1]*(Distance[0]);
    
    /*--- Store the grid velocity at this node ---*/
    
    node[iPoint_Coarse]->SetGridVel(RotVel);
    
  }
  
  delete [] RotVel;
  
}

void CMultiGridGeometry::SetTranslationalVelocity(CConfig *config) {
  
  unsigned iDim;
  unsigned long iPoint_Coarse;
  su2double xDot[3] = {0.0,0.0,0.0};
  
  /*--- Get the translational velocity vector from config ---*/
  
  xDot[0]   = config->GetTranslation_Rate_X(ZONE_0)/config->GetVelocity_Ref();
  xDot[1]   = config->GetTranslation_Rate_Y(ZONE_0)/config->GetVelocity_Ref();
  xDot[2]   = config->GetTranslation_Rate_Z(ZONE_0)/config->GetVelocity_Ref();
  
  /*--- Loop over all nodes and set the translational velocity ---*/
  
  for (iPoint_Coarse = 0; iPoint_Coarse < nPoint; iPoint_Coarse++) {
    
    /*--- Store the grid velocity at this node ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      node[iPoint_Coarse]->SetGridVel(iDim,xDot[iDim]);
    
  }
  
}

void CMultiGridGeometry::SetGridVelocity(CConfig *config, unsigned long iter) {
  
  /*--- Local variables ---*/
  
  su2double *Coord_nP1 = NULL, *Coord_n = NULL, *Coord_nM1 = NULL;
  su2double TimeStep, GridVel = 0.0;
  unsigned long Point_Coarse;
  unsigned short iDim;
  
  /*--- Compute the velocity of each node in the volume mesh ---*/
  
  for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
    
    /*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/
    
    Coord_nM1 = node[Point_Coarse]->GetCoord_n1();
    Coord_n   = node[Point_Coarse]->GetCoord_n();
    Coord_nP1 = node[Point_Coarse]->GetCoord();
    
    /*--- Unsteady time step ---*/
    
    TimeStep = config->GetDelta_UnstTimeND();
    
    /*--- Compute mesh velocity with 1st or 2nd-order approximation ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        GridVel = ( Coord_nP1[iDim] - Coord_n[iDim] ) / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        GridVel = ( 3.0*Coord_nP1[iDim] - 4.0*Coord_n[iDim]
                   +  1.0*Coord_nM1[iDim] ) / (2.0*TimeStep);
      
      /*--- Store grid velocity for this point ---*/
      
      node[Point_Coarse]->SetGridVel(iDim, GridVel);
      
    }
  }
}

void CMultiGridGeometry::SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config) {
  
  /*--- Local variables ---*/
  unsigned short iDim, iChild;
  unsigned long Point_Coarse, Point_Fine;
  su2double Area_Parent, Area_Child, Grid_Vel[3], *Grid_Vel_Fine;
  
  /*--- Loop over all coarse mesh points ---*/
  for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
    Area_Parent = node[Point_Coarse]->GetVolume();
    
    /*--- Zero out the grid velocity ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Grid_Vel[iDim] = 0.0;
    
    /*--- Loop over all of the children for this coarse CV and compute
     a grid velocity based on the values in the child CVs (fine mesh). ---*/
    for (iChild = 0; iChild < node[Point_Coarse]->GetnChildren_CV(); iChild++) {
      Point_Fine    = node[Point_Coarse]->GetChildren_CV(iChild);
      Area_Child    = fine_mesh->node[Point_Fine]->GetVolume();
      Grid_Vel_Fine = fine_mesh->node[Point_Fine]->GetGridVel();
      for (iDim = 0; iDim < nDim; iDim++)
        Grid_Vel[iDim] += Grid_Vel_Fine[iDim]*Area_Child/Area_Parent;
    }
    
    /*--- Set the grid velocity for this coarse node. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      node[Point_Coarse]->SetGridVel(iDim, Grid_Vel[iDim]);
  }
}


void CMultiGridGeometry::FindNormal_Neighbor(CConfig *config) {
  
  unsigned short iMarker, iDim;
  unsigned long iPoint, iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY ) {
      
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        iPoint = vertex[iMarker][iVertex]->GetNode();
        
        /*--- If the node belong to the domain ---*/
        if (node[iPoint]->GetDomain()) {
          
          /*--- Compute closest normal neighbor ---*/
          su2double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord;
          unsigned long Point_Normal = 0, jPoint;
          unsigned short iNeigh;
          su2double *Normal = vertex[iMarker][iVertex]->GetNormal();
          cos_max = -1.0;
          for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
            jPoint = node[iPoint]->GetPoint(iNeigh);
            scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
              scalar_prod += diff_coord*Normal[iDim];
              norm_vect += diff_coord*diff_coord;
              norm_Normal += Normal[iDim]*Normal[iDim];
            }
            norm_vect = sqrt(norm_vect);
            norm_Normal = sqrt(norm_Normal);
            cos_alpha = scalar_prod/(norm_vect*norm_Normal);
            
            /*--- Get maximum cosine (not minimum because normals are oriented inwards) ---*/
            if (cos_alpha >= cos_max) {
              Point_Normal = jPoint;
              cos_max = cos_alpha;
            }
          }
          vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
        }
      }
    }
  }
}


void CMultiGridGeometry::SetGeometryPlanes(CConfig *config) {
  bool loop_on;
  unsigned short iMarker = 0;
  su2double auxXCoord, auxYCoord, auxZCoord,	*Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
  unsigned long jVertex, iVertex, ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;
  
  /*--- Compute the total number of points on the near-field ---*/
  nVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      nVertex_Wall += nVertex[iMarker];
  
  
  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  Xcoord = new su2double[nVertex_Wall];
  Ycoord = new su2double[nVertex_Wall];
  if (nDim == 3)	Zcoord = new su2double[nVertex_Wall];
  FaceArea = new su2double[nVertex_Wall];
  
  /*--- Copy the boundary information to an array ---*/
  iVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
        Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
        if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
        Face_Normal = vertex[iMarker][iVertex]->GetNormal();
        FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
        iVertex_Wall ++;
      }
  
  
  //vector<su2double> XCoordList;
  vector<su2double>::iterator IterXCoordList;
  
  for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
    XCoordList.push_back(Xcoord[iVertex]);
  
  sort( XCoordList.begin(), XCoordList.end());
  IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
  XCoordList.resize( IterXCoordList - XCoordList.begin() );
  
  /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
  Xcoord_plane.resize(XCoordList.size());
  Ycoord_plane.resize(XCoordList.size());
  if (nDim==3) Zcoord_plane.resize(XCoordList.size());
  FaceArea_plane.resize(XCoordList.size());
  Plane_points.resize(XCoordList.size());
  
  
  su2double dist_ratio;
  unsigned long iCoord;
  
  /*--- Distribute the values among the different PhiAngles ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (node[iPoint]->GetDomain()) {
      loop_on = true;
      for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
        dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
        if (dist_ratio >= 0 && dist_ratio <= 1.0) {
          if (dist_ratio <= 0.5) iCoord = ixCoord;
          else iCoord = ixCoord+1;
          Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
          Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
          if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
          FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
          Plane_points[iCoord].push_back(iPoint );
          loop_on = false;
        }
      }
    }
  }
  
  unsigned long auxPoint;
  /*--- Order the arrays in ascending values of y ---*/
  for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
    for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
      for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
        if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
          auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex+1]; Xcoord_plane[ixCoord][jVertex+1] = auxXCoord;
          auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex+1]; Ycoord_plane[ixCoord][jVertex+1] = auxYCoord;
          auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex+1]; Plane_points[ixCoord][jVertex+1] = auxPoint;
          if (nDim==3) {
            auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex+1]; Zcoord_plane[ixCoord][jVertex+1] = auxZCoord;
          }
          auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex+1]; FaceArea_plane[ixCoord][jVertex+1] = auxArea;
        }
  
  /*--- Delete structures ---*/
  delete[] Xcoord; delete[] Ycoord;
  if (Zcoord != NULL) delete[] Zcoord;
  delete[] FaceArea;
}

CPeriodicGeometry::CPeriodicGeometry(CGeometry *geometry, CConfig *config) {
  unsigned long nElem_new, nPoint_new, jPoint, iPoint, iElem, jElem, iVertex,
  nelem_triangle = 0, nelem_quad = 0, nelem_tetra = 0, nelem_hexa = 0, nelem_prism = 0,
  nelem_pyramid = 0, iIndex, newElementsBound = 0;
  unsigned short  iMarker, nPeriodic = 0, iPeriodic;
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
  dx, dy, dz, rotCoord[3], *Coord_i;
  unsigned short nMarker_Max = config->GetnMarker_Max();

  /*--- We only create the mirror structure for the second boundary ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      /*--- Evaluate the number of periodic boundary conditions ---*/
      nPeriodic++;
    }
  }
  bool *CreateMirror = new bool[nPeriodic+1];
  CreateMirror[0] = false;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (iPeriodic <= nPeriodic/2) CreateMirror[iPeriodic] = false;
    else CreateMirror[iPeriodic] = true;
  }
  
  /*--- Write the number of dimensions of the problem ---*/
  nDim = geometry->GetnDim();
  
  /*--- Copy the new boundary element information from the geometry class.
   Be careful, as these are pointers to vectors/objects. ---*/
  nNewElem_BoundPer = geometry->nNewElem_Bound;
  newBoundPer       = geometry->newBound;
  
  /*--- Count the number of new boundary elements. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    newElementsBound += nNewElem_BoundPer[iMarker];
  
  /*--- Loop over the original grid to perform the dimensionalizaton of the new vectors ---*/
  nElem_new = 0; nPoint_new = 0;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      nElem_new += geometry->PeriodicElem[iPeriodic].size();
      nPoint_new += geometry->PeriodicPoint[iPeriodic][0].size();
    }
  }
  
  cout << "Number of new points: " << nPoint_new << "." << endl;
  cout << "Number of new interior elements: " << nElem_new << "." << endl;
  cout << "Number of new boundary elements added to preexisting markers: " << newElementsBound << "." << endl;
  
  /*--- Create a copy of the original grid ---*/
  elem = new CPrimalGrid*[geometry->GetnElem() + nElem_new];
  for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
    switch(geometry->elem[iElem]->GetVTK_Type()) {
      case TRIANGLE:
        elem[iElem] = new CTriangle(geometry->elem[iElem]->GetNode(0),
                                    geometry->elem[iElem]->GetNode(1),
                                    geometry->elem[iElem]->GetNode(2), 2);
        nelem_triangle++;
        break;
        
      case QUADRILATERAL:
        elem[iElem] = new CQuadrilateral(geometry->elem[iElem]->GetNode(0),
                                     geometry->elem[iElem]->GetNode(1),
                                     geometry->elem[iElem]->GetNode(2),
                                     geometry->elem[iElem]->GetNode(3), 2);
        nelem_quad++;
        break;
        
      case TETRAHEDRON:
        elem[iElem] = new CTetrahedron(geometry->elem[iElem]->GetNode(0),
                                       geometry->elem[iElem]->GetNode(1),
                                       geometry->elem[iElem]->GetNode(2),
                                       geometry->elem[iElem]->GetNode(3));
        nelem_tetra++;
        break;
        
      case HEXAHEDRON:
        elem[iElem] = new CHexahedron(geometry->elem[iElem]->GetNode(0),
                                      geometry->elem[iElem]->GetNode(1),
                                      geometry->elem[iElem]->GetNode(2),
                                      geometry->elem[iElem]->GetNode(3),
                                      geometry->elem[iElem]->GetNode(4),
                                      geometry->elem[iElem]->GetNode(5),
                                      geometry->elem[iElem]->GetNode(6),
                                      geometry->elem[iElem]->GetNode(7));
        nelem_hexa++;
        break;
        
      case PRISM:
        elem[iElem] = new CPrism(geometry->elem[iElem]->GetNode(0),
                                 geometry->elem[iElem]->GetNode(1),
                                 geometry->elem[iElem]->GetNode(2),
                                 geometry->elem[iElem]->GetNode(3),
                                 geometry->elem[iElem]->GetNode(4),
                                 geometry->elem[iElem]->GetNode(5));
        nelem_prism++;
        break;
        
      case PYRAMID:
        elem[iElem] = new CPyramid(geometry->elem[iElem]->GetNode(0),
                                   geometry->elem[iElem]->GetNode(1),
                                   geometry->elem[iElem]->GetNode(2),
                                   geometry->elem[iElem]->GetNode(3),
                                   geometry->elem[iElem]->GetNode(4));
        nelem_pyramid++;
        break;
        
    }
  }
  
  /*--- Create a list with all the points and the new index ---*/
  unsigned long *Index = new unsigned long [geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) Index[iPoint] = 0;
  
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
        iPoint =  geometry->PeriodicPoint[iPeriodic][0][iIndex];
        Index[iPoint] = geometry->PeriodicPoint[iPeriodic][1][iIndex];
      }
    }
  }
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
        Index[iPoint] = jPoint;
      }
  
  /*--- Add the new elements due to the periodic boundary condtion ---*/
  iElem = geometry->GetnElem();
  
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      for (iIndex = 0; iIndex < geometry->PeriodicElem[iPeriodic].size(); iIndex++) {
        jElem = geometry->PeriodicElem[iPeriodic][iIndex];
        
        switch(geometry->elem[jElem]->GetVTK_Type()) {
          case TRIANGLE:
            elem[iElem] = new CTriangle(Index[geometry->elem[jElem]->GetNode(0)],
                                        Index[geometry->elem[jElem]->GetNode(1)],
                                        Index[geometry->elem[jElem]->GetNode(2)], 2);
            iElem++; nelem_triangle++;
            break;
            
          case QUADRILATERAL:
            elem[iElem] = new CQuadrilateral(Index[geometry->elem[jElem]->GetNode(0)],
                                         Index[geometry->elem[jElem]->GetNode(1)],
                                         Index[geometry->elem[jElem]->GetNode(2)],
                                         Index[geometry->elem[jElem]->GetNode(3)], 2);
            iElem++; nelem_quad++;
            break;
            
          case TETRAHEDRON:
            elem[iElem] = new CTetrahedron(Index[geometry->elem[jElem]->GetNode(0)],
                                           Index[geometry->elem[jElem]->GetNode(1)],
                                           Index[geometry->elem[jElem]->GetNode(2)],
                                           Index[geometry->elem[jElem]->GetNode(3)]);
            iElem++; nelem_tetra++;
            break;
            
          case HEXAHEDRON:
            elem[iElem] = new CHexahedron(Index[geometry->elem[jElem]->GetNode(0)],
                                          Index[geometry->elem[jElem]->GetNode(1)],
                                          Index[geometry->elem[jElem]->GetNode(2)],
                                          Index[geometry->elem[jElem]->GetNode(3)],
                                          Index[geometry->elem[jElem]->GetNode(4)],
                                          Index[geometry->elem[jElem]->GetNode(5)],
                                          Index[geometry->elem[jElem]->GetNode(6)],
                                          Index[geometry->elem[jElem]->GetNode(7)]);
            iElem++; nelem_hexa++;
            break;
            
          case PRISM:
            elem[iElem] = new CPrism(Index[geometry->elem[jElem]->GetNode(0)],
                                     Index[geometry->elem[jElem]->GetNode(1)],
                                     Index[geometry->elem[jElem]->GetNode(2)],
                                     Index[geometry->elem[jElem]->GetNode(3)],
                                     Index[geometry->elem[jElem]->GetNode(4)],
                                     Index[geometry->elem[jElem]->GetNode(5)]);
            iElem++; nelem_prism++;
            break;
            
          case PYRAMID:
            elem[iElem] = new CPyramid(Index[geometry->elem[jElem]->GetNode(0)],
                                       Index[geometry->elem[jElem]->GetNode(1)],
                                       Index[geometry->elem[jElem]->GetNode(2)],
                                       Index[geometry->elem[jElem]->GetNode(3)],
                                       Index[geometry->elem[jElem]->GetNode(4)]);
            iElem++; nelem_pyramid++;
            break;
            
        }
      }
    }
  }
  
  nElem = geometry->GetnElem() + nElem_new;
  
  /*--- Add the old points ---*/

  nPointNode = geometry->GetnPoint() + nPoint_new;
  node = new CPoint*[geometry->GetnPoint() + nPoint_new];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    if (geometry->GetnDim() == 2)
      node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
                                geometry->node[iPoint]->GetCoord(1), iPoint, config);
    if (geometry->GetnDim() == 3)
      node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
                                geometry->node[iPoint]->GetCoord(1),
                                geometry->node[iPoint]->GetCoord(2), iPoint, config);
  }
  
  /*--- Add the new points due to the periodic boundary condtion (only in the mirror part) ---*/
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
        
        /*--- From iPeriodic obtain the iMarker ---*/
        for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
          if (iPeriodic == config->GetMarker_All_PerBound(iMarker)) break;
        
        /*--- Retrieve the supplied periodic information. ---*/
        center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
        angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
        trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
        
        /*--- Store center - trans as it is constant and will be added on.
         Note the subtraction, as this is the inverse translation. ---*/
        translation[0] = center[0] - trans[0];
        translation[1] = center[1] - trans[1];
        translation[2] = center[2] - trans[2];
        
        /*--- Store angles separately for clarity. Compute sines/cosines.
         Note the negative sign, as this is the inverse rotation. ---*/
        theta = -angles[0];
        phi   = -angles[1];
        psi   = -angles[2];
        
        cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = cosPhi*sinPsi;
        rotMatrix[2][0] = -sinPhi;
        
        rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = sinTheta*cosPhi;
        
        rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Retrieve node information for this boundary point. ---*/
        iPoint = geometry->PeriodicPoint[iPeriodic][0][iIndex];
        jPoint = geometry->PeriodicPoint[iPeriodic][1][iIndex];
        Coord_i = geometry->node[iPoint]->GetCoord();
        
        /*--- Get the position vector from rot center to point. ---*/
        dx = Coord_i[0] - center[0];
        dy = Coord_i[1] - center[1];
        if (nDim == 3) {
          dz = Coord_i[2] - center[2];
        } else {
          dz = 0.0;
        }
        
        /*--- Compute transformed point coordinates. ---*/
        rotCoord[0] = rotMatrix[0][0]*dx + rotMatrix[0][1]*dy + rotMatrix[0][2]*dz + translation[0];
        rotCoord[1] = rotMatrix[1][0]*dx + rotMatrix[1][1]*dy + rotMatrix[1][2]*dz + translation[1];
        rotCoord[2] = rotMatrix[2][0]*dx + rotMatrix[2][1]*dy + rotMatrix[2][2]*dz + translation[2];
        
        /*--- Save the new points with the new coordinates. ---*/
        if (geometry->GetnDim() == 2)
          node[jPoint] = new CPoint(rotCoord[0], rotCoord[1], jPoint, config);
        if (geometry->GetnDim() == 3)
          node[jPoint] = new CPoint(rotCoord[0], rotCoord[1], rotCoord[2], jPoint, config);
        
      }
    }
  }
  
  nPoint = geometry->GetnPoint() + nPoint_new;
  
  /*--- Add the old boundary, reserving space for two new bc (send/recive periodic bc) ---*/
  nMarker = geometry->GetnMarker() + 2;
  nElem_Bound = new unsigned long [nMarker];
  bound = new CPrimalGrid**[nMarker];
  Tag_to_Marker = new string [nMarker_Max];
  config->SetnMarker_All(nMarker);
  
  /*--- Copy the olf boundary ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    
    bound[iMarker] = new CPrimalGrid* [geometry->GetnElem_Bound(iMarker)];
    
    for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
      if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == LINE)
        bound[iMarker][iVertex] = new CLine(geometry->bound[iMarker][iVertex]->GetNode(0),
                                            geometry->bound[iMarker][iVertex]->GetNode(1), 2);
      if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TRIANGLE)
        bound[iMarker][iVertex] = new CTriangle(geometry->bound[iMarker][iVertex]->GetNode(0),
                                                geometry->bound[iMarker][iVertex]->GetNode(1),
                                                geometry->bound[iMarker][iVertex]->GetNode(2), 3);
      if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == QUADRILATERAL)
        bound[iMarker][iVertex] = new CQuadrilateral(geometry->bound[iMarker][iVertex]->GetNode(0),
                                                 geometry->bound[iMarker][iVertex]->GetNode(1),
                                                 geometry->bound[iMarker][iVertex]->GetNode(2),
                                                 geometry->bound[iMarker][iVertex]->GetNode(3), 3);
    }
    
    nElem_Bound[iMarker] = geometry->GetnElem_Bound(iMarker);
    Tag_to_Marker[iMarker] = geometry->GetMarker_Tag(iMarker);
    
  }
  
  delete [] Index;
  delete [] CreateMirror;
  
}

CPeriodicGeometry::~CPeriodicGeometry(void) {
  unsigned long iElem_Bound;
  unsigned short iMarker;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
      if (newBoundPer[iMarker][iElem_Bound] != NULL) delete [] newBoundPer[iMarker][iElem_Bound];
    }
  }
  if (newBoundPer != NULL) delete[] newBoundPer;
  
  if (nNewElem_BoundPer != NULL) delete[] nNewElem_BoundPer;
  
}

void CPeriodicGeometry::SetPeriodicBoundary(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, iPeriodic, nPeriodic = 0, iMarkerSend, iMarkerReceive;
  unsigned long iVertex, Counter_Send = 0, Counter_Receive = 0, iIndex;
  
  /*--- Compute the number of periodic bc on the geometry ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      nPeriodic++;
  
  /*--- First compute the Send/Receive boundaries, count the number of points ---*/
  Counter_Send = 0; 	Counter_Receive = 0;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (geometry->PeriodicPoint[iPeriodic][0].size() != 0)
      Counter_Send += geometry->PeriodicPoint[iPeriodic][0].size();
    if (geometry->PeriodicPoint[iPeriodic][1].size() != 0)
      Counter_Receive += geometry->PeriodicPoint[iPeriodic][1].size();
  }
  
  /*--- Adimensionalization of the new boundaries ---*/
  iMarkerSend = nMarker - 2; iMarkerReceive = nMarker - 1;
  config->SetMarker_All_SendRecv(iMarkerSend,1);
  config->SetMarker_All_SendRecv(iMarkerReceive,-1);
  nElem_Bound[iMarkerSend] = Counter_Send;
  nElem_Bound[iMarkerReceive] = Counter_Receive;
  bound[iMarkerSend] = new CPrimalGrid* [Counter_Send];
  bound[iMarkerReceive] = new CPrimalGrid* [Counter_Receive];
  
  /*--- First we do the send ---*/
  iVertex = 0;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
    if (geometry->PeriodicPoint[iPeriodic][0].size() != 0)
      for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
        bound[iMarkerSend][iVertex] = new CVertexMPI(geometry->PeriodicPoint[iPeriodic][0][iIndex], nDim);
        bound[iMarkerSend][iVertex]->SetRotation_Type(iPeriodic);
        iVertex++;
      }
  
  /*--- Second we do the receive ---*/
  iVertex = 0;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
    if (geometry->PeriodicPoint[iPeriodic][1].size() != 0)
      for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][1].size(); iIndex++) {
        bound[iMarkerReceive][iVertex] = new CVertexMPI(geometry->PeriodicPoint[iPeriodic][1][iIndex], nDim);
        bound[iMarkerReceive][iVertex]->SetRotation_Type(iPeriodic);
        iVertex++;
      }
  
}

void CPeriodicGeometry::SetMeshFile(CGeometry *geometry, CConfig *config, string val_mesh_out_filename) {
  unsigned long iElem, iPoint, iElem_Bound, GhostPoints;
  unsigned short iMarker, iNodes, iDim;
  unsigned short iMarkerReceive, iPeriodic, nPeriodic = 0;
  ofstream output_file;
  string Grid_Marker;
  char *cstr;
  su2double *center, *angles, *transl;
  
  cstr = new char [val_mesh_out_filename.size()+1];
  strcpy (cstr, val_mesh_out_filename.c_str());
  
  /*--- Open .su2 grid file ---*/
  output_file.precision(15);
  output_file.open(cstr, ios::out);
  
  /*--- Ghost points, look at the nodes in the send receive ---*/
  iMarkerReceive = nMarker - 1;
  GhostPoints = nElem_Bound[iMarkerReceive];
  
  /*--- Change the numbering to guarantee that the all the receive
   points are at the end of the file ---*/
  unsigned long OldnPoint = geometry->GetnPoint();
  unsigned long *NewSort = new unsigned long[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    NewSort[iPoint] = iPoint;
  }
  
  unsigned long Index = OldnPoint-1;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
      if (config->GetMarker_All_SendRecv(iMarker) < 0) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          if (bound[iMarker][iElem_Bound]->GetNode(0) < geometry->GetnPoint()) {
            NewSort[bound[iMarker][iElem_Bound]->GetNode(0)] = Index;
            NewSort[Index] = bound[iMarker][iElem_Bound]->GetNode(0);
            Index--;
          }
        }
      }
    }
  }
  
  
  /*--- Write dimension, number of elements and number of points ---*/
  output_file << "NDIME= " << nDim << endl;
  output_file << "NELEM= " << nElem << endl;
  for (iElem = 0; iElem < nElem; iElem++) {
    output_file << elem[iElem]->GetVTK_Type();
    for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
      output_file << "\t" << NewSort[elem[iElem]->GetNode(iNodes)];
    output_file << "\t"<<iElem<< endl;
  }
  
  output_file << "NPOIN= " << nPoint << "\t" << nPoint - GhostPoints << endl;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      output_file << scientific << "\t" << node[NewSort[iPoint]]->GetCoord(iDim) ;
    output_file << "\t" << iPoint << endl;
  }
  
  output_file << "NMARK= " << nMarker << endl;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      
      Grid_Marker = config->GetMarker_All_TagBound(iMarker);
      output_file << "MARKER_TAG= " << Grid_Marker << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker] + nNewElem_BoundPer[iMarker] << endl;
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
        for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes()-1; iNodes++)
          output_file << NewSort[bound[iMarker][iElem_Bound]->GetNode(iNodes)] << "\t" ;
        iNodes = bound[iMarker][iElem_Bound]->GetnNodes()-1;
        output_file << NewSort[bound[iMarker][iElem_Bound]->GetNode(iNodes)] << endl;
      }
      
      /*--- Write any new elements at the end of the list. ---*/
      if (nNewElem_BoundPer[iMarker] > 0) {
        for (iElem_Bound = 0; iElem_Bound < nNewElem_BoundPer[iMarker]; iElem_Bound++) {
          output_file << newBoundPer[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < newBoundPer[iMarker][iElem_Bound]->GetnNodes()-1; iNodes++)
            output_file << NewSort[newBoundPer[iMarker][iElem_Bound]->GetNode(iNodes)] << "\t" ;
          iNodes = newBoundPer[iMarker][iElem_Bound]->GetnNodes()-1;
          output_file << NewSort[newBoundPer[iMarker][iElem_Bound]->GetNode(iNodes)] << endl;
        }
      }
      
    }
    
    if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
      output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
        NewSort[bound[iMarker][iElem_Bound]->GetNode(0)] << "\t" <<
        bound[iMarker][iElem_Bound]->GetRotation_Type()  << endl;
      }
    }
  }
  
  /*--- Compute the number of periodic bc on the geometry ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      nPeriodic++;
  
  output_file << "NPERIODIC= " << nPeriodic + 1 << endl;
  
  /*--- Periodic 0 correspond with no movement of the surface ---*/
  output_file << "PERIODIC_INDEX= 0" << endl;
  output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
  output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
  output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
  
  /*--- From iPeriodic obtain the iMarker ---*/
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (iPeriodic == config->GetMarker_All_PerBound(iMarker)) break;
    
    /*--- Retrieve the supplied periodic information. ---*/
    center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
    angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
    transl = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
    
    output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
    output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
    output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
    output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
    
  }
  
  
  output_file.close();
  
  /*--- Free memory ---*/
  delete [] NewSort;
  
}

void CPeriodicGeometry::SetTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file) {
  
  unsigned long iElem, iPoint;
  unsigned short iDim;
  ofstream Tecplot_File;
  
  Tecplot_File.open(mesh_filename, ios::out);
  Tecplot_File << "TITLE= \"Visualization of the volumetric grid\"" << endl;
  
  if (nDim == 2) {
    Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
  }
  if (nDim == 3) {
    Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
    Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
    Tecplot_File << "\n";
  }
  
  for (iElem = 0; iElem < nElem; iElem++) {
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
      elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PRISM) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
    }
  }
  
  Tecplot_File.close();
}

CMultiGridQueue::CMultiGridQueue(unsigned long val_npoint) {
  unsigned long iPoint;
  
  nPoint = val_npoint;
  Priority = new short[nPoint];
  RightCV = new bool[nPoint];
  
  QueueCV.resize(1);
  
  /*--- Queue initialization with all the points in the finer grid ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    QueueCV[0].push_back(iPoint);
    Priority[iPoint] = 0;
    RightCV[iPoint] = true;
  }
  
}

CMultiGridQueue::~CMultiGridQueue(void) {
  
  delete[] Priority;
  delete[] RightCV;
  
}

void CMultiGridQueue::AddCV(unsigned long val_new_point, unsigned short val_number_neighbors) {
  
  unsigned short Max_Neighbors = QueueCV.size()-1;
  
  /*--- Basic check ---*/
  if (val_new_point > nPoint) {
    cout << "The index of the CV is greater than the size of the priority list." << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Resize the list ---*/
  if (val_number_neighbors > Max_Neighbors)
    QueueCV.resize(val_number_neighbors+1);
  
  /*--- Find the point in the queue ---*/
  bool InQueue = false;
  if (Priority[val_new_point] == val_number_neighbors) InQueue = true;
  
  if (!InQueue) {
    /*--- Add the control volume, and update the priority list ---*/
    QueueCV[val_number_neighbors].push_back(val_new_point);
    Priority[val_new_point] = val_number_neighbors;
  }
  
}

void CMultiGridQueue::RemoveCV(unsigned long val_remove_point) {
  unsigned short iPoint;
  bool check;
  
  /*--- Basic check ---*/
  if (val_remove_point > nPoint) {
    cout << "The index of the CV is greater than the size of the priority list." << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Find priority of the Control Volume ---*/
  short Number_Neighbors = Priority[val_remove_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_remove_point <<" is not in the priority list. (RemoveCV)" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Find the point in the queue ---*/
  vector<unsigned long>::iterator ItQueue = find(QueueCV[Number_Neighbors].begin(),
                                                 QueueCV[Number_Neighbors].end(),
                                                 val_remove_point);
  if ( ItQueue != QueueCV[Number_Neighbors].end() ) QueueCV[Number_Neighbors].erase(ItQueue);
  
  Priority[val_remove_point] = -1;
  
  /*--- Check that the size of the queue is the right one ---*/
  unsigned short Size_QueueCV = 0;
  check = false;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++)
    if (QueueCV[iPoint].size() != 0) { Size_QueueCV = iPoint; check = true;}
  
  /*--- Resize the queue, if check = false, the queue is empty, at least
   we need one element in the queue ---*/
  if (check) QueueCV.resize(Size_QueueCV+1);
  else QueueCV.resize(1);
  
}

void CMultiGridQueue::MoveCV(unsigned long val_move_point, short val_number_neighbors) {
  
  if (val_number_neighbors < 0) {
    val_number_neighbors = 0;
    RightCV[val_move_point] = false;
  }
  else {
    RightCV[val_move_point] = true;
  }
  
  /*--- Remove the control volume ---*/
  RemoveCV(val_move_point);
  
  /*--- Add a new control volume ---*/
  AddCV(val_move_point, val_number_neighbors);
  
}

void CMultiGridQueue::IncrPriorityCV(unsigned long val_incr_point) {
  
  /*--- Find the priority list ---*/
  short Number_Neighbors = Priority[val_incr_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_incr_point <<" is not in the priority list. (IncrPriorityCV)" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Remove the control volume ---*/
  RemoveCV(val_incr_point);
  
  /*--- Increase the priority ---*/
  AddCV(val_incr_point, Number_Neighbors+1);
  
}

void CMultiGridQueue::RedPriorityCV(unsigned long val_red_point) {
  
  /*--- Find the priority list ---*/
  short Number_Neighbors = Priority[val_red_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_red_point <<" is not in the priority list. (RedPriorityCV)" << endl;
    exit(EXIT_FAILURE);
  }
  
  if (Number_Neighbors != 0) {
    
    /*--- Remove the control volume ---*/
    RemoveCV(val_red_point);
    
    /*--- Increase the priority ---*/
    AddCV(val_red_point, Number_Neighbors-1);
    
  }
  
}

void CMultiGridQueue::VisualizeQueue(void) {
  unsigned short iPoint;
  unsigned long jPoint;
  
  cout << endl;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++) {
    cout << "Number of neighbors " << iPoint <<": ";
    for (jPoint = 0; jPoint < QueueCV[iPoint].size(); jPoint ++) {
      cout << QueueCV[iPoint][jPoint] << " ";
    }
    cout << endl;
  }
  
}

void CMultiGridQueue::VisualizePriority(void) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++)
    cout << "Control Volume: " << iPoint <<" Priority: " << Priority[iPoint] << endl;
  
}

long CMultiGridQueue::NextCV(void) {
  if (QueueCV.size() != 0) return QueueCV[QueueCV.size()-1][0];
  else return -1;
}

bool CMultiGridQueue::EmptyQueue(void) {
  unsigned short iPoint;
  
  /*--- In case there is only the no agglomerated elements,
   check if they can be agglomerated or we have already finished ---*/
  bool check = true;
  
  if ( QueueCV.size() == 1 ) {
    for (iPoint = 0; iPoint < QueueCV[0].size(); iPoint ++) {
      if (RightCV[QueueCV[0][iPoint]]) { check = false; break; }
    }
  }
  else {
    for (iPoint = 1; iPoint < QueueCV.size(); iPoint ++)
      if (QueueCV[iPoint].size() != 0) { check = false; break;}
  }
  
  return check;
}

unsigned long CMultiGridQueue::TotalCV(void) {
  unsigned short iPoint;
  unsigned long TotalCV;
  
  TotalCV = 0;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++)
    if (QueueCV[iPoint].size() != 0) { TotalCV += QueueCV[iPoint].size(); }
  
  return TotalCV;
}

void CMultiGridQueue::Update(unsigned long iPoint, CGeometry *fine_grid) {
  unsigned short iNode;
  unsigned long jPoint;
  
  RemoveCV(iPoint);
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    if (fine_grid->node[jPoint]->GetAgglomerate() == false)
      IncrPriorityCV(jPoint);
  }
  
}
