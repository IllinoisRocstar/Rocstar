/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   RI_DFT.cpp
  \brief  

  \author Thomas Leurent

  \date   2004-04-12
*/

#include <cstdio>
#include "RI_DFT.hpp"

using namespace Mesquite;

/*****
  Functions for:
	               ||T*T' - I||_F^(2b)
	a * ---------------------------------------
	    (det(T) + sqrt(det(T)^2 + 4*delta^2))^c

  Note that for this metric, the function is equivalent to:
	               ||T'*T - I||_F^(2b)
	a * ---------------------------------------
	    (det(T) + sqrt(det(T)^2 + 4*delta^2))^c

  Note: c can be any value since the denominator is always restricted 
  to be positive.  However, b can only take on values of one, two, or
  any number greater than 2 because the gradient and/or Hessian may
  not exist in these other cases.  In particular, you cannot guarantee
  that the numerator is bounded away from zero.
*****/

inline bool m_fcn_ridft2(double &obj, 
			 const Vector3D x[3], const Vector3D &n,
			 const Matrix3D &invW,
			 const double a, const Exponent b, const Exponent c,
			 const double delta)
{
  static double matr[9], f, t1, t2;
  static double fmat[6], g;

  /* Calculate M = A*inv(W). */
  f       = x[1][0] - x[0][0];
  g       = x[2][0] - x[0][0];
  t1      = n[0];
  matr[0] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[1] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[2] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][1] - x[0][1];
  g       = x[2][1] - x[0][1];
  t1      = n[1];
  matr[3] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[4] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[5] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][2] - x[0][2];
  g       = x[2][2] - x[0][2];
  t1      = n[2];
  matr[6] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[7] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[8] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  /* Calculate det(M). */
  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  t2 = sqrt(t1*t1 + 4.0*delta*delta);
  g = t1 + t2;
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] - 1.0;
  fmat[1] = matr[0]*matr[3] + matr[1]*matr[4] + matr[2]*matr[5];
  fmat[2] = matr[0]*matr[6] + matr[1]*matr[7] + matr[2]*matr[8];

  fmat[3] = matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] - 1.0;
  fmat[4] = matr[3]*matr[6] + matr[4]*matr[7] + matr[5]*matr[8];

  fmat[5] = matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  /* Calculate objective function. */
  obj = a * pow(f, b) * pow(g, c);
  return true;
}

inline bool g_fcn_ridft2(double &obj, Vector3D g_obj[3], 
			 const Vector3D x[3], const Vector3D &n,
			 const Matrix3D &invW,
			 const double a, const Exponent b, const Exponent c,
			 const double delta)
{
  static double matr[9], f, t1, t2;
  static double fmat[6], g;
  static double adj_m[9], df[9], loc1, loc2, loc3;

  /* Calculate M = A*inv(W). */
  f       = x[1][0] - x[0][0];
  g       = x[2][0] - x[0][0];
  t1      = n[0];
  matr[0] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[1] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[2] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][1] - x[0][1];
  g       = x[2][1] - x[0][1];
  t1      = n[1];
  matr[3] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[4] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[5] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][2] - x[0][2];
  g       = x[2][2] - x[0][2];
  t1      = n[2];
  matr[6] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[7] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[8] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  /* Calculate det(M). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[5]*matr[6] - matr[3]*matr[8];
  loc3 = matr[3]*matr[7] - matr[4]*matr[6];
  t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;
  t2 = sqrt(t1*t1 + 4.0*delta*delta);
  g = t1 + t2;
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] - 1.0;
  fmat[1] = matr[0]*matr[3] + matr[1]*matr[4] + matr[2]*matr[5];
  fmat[2] = matr[0]*matr[6] + matr[1]*matr[7] + matr[2]*matr[8];

  fmat[3] = matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] - 1.0;
  fmat[4] = matr[3]*matr[6] + matr[4]*matr[7] + matr[5]*matr[8];

  fmat[5] = matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  /* Calculate objective function. */
  obj = a * pow(f, b) * pow(g, c);

  /* Calculate the derivative of the objective function. */
  f = b * obj / f * 4.0;              /* Constant on nabla f */
  g = c * obj / g * (1 + t1 / t2);    /* Constant on nabla g */

  df[0] = fmat[0]*matr[0] + fmat[1]*matr[3] + fmat[2]*matr[6];
  df[1] = fmat[0]*matr[1] + fmat[1]*matr[4] + fmat[2]*matr[7];
  df[2] = fmat[0]*matr[2] + fmat[1]*matr[5] + fmat[2]*matr[8];

  df[3] = fmat[1]*matr[0] + fmat[3]*matr[3] + fmat[4]*matr[6];
  df[4] = fmat[1]*matr[1] + fmat[3]*matr[4] + fmat[4]*matr[7];
  df[5] = fmat[1]*matr[2] + fmat[3]*matr[5] + fmat[4]*matr[8];

  df[6] = fmat[2]*matr[0] + fmat[4]*matr[3] + fmat[5]*matr[6];
  df[7] = fmat[2]*matr[1] + fmat[4]*matr[4] + fmat[5]*matr[7];
  df[8] = fmat[2]*matr[2] + fmat[4]*matr[5] + fmat[5]*matr[8];

  adj_m[0] = df[0]*f + loc1*g;
  adj_m[1] = df[1]*f + loc2*g;
  adj_m[2] = df[2]*f + loc3*g;

  loc1 = matr[0]*g;
  loc2 = matr[1]*g;
  loc3 = matr[2]*g;

  adj_m[3] = df[3]*f + loc3*matr[7] - loc2*matr[8];
  adj_m[4] = df[4]*f + loc1*matr[8] - loc3*matr[6];
  adj_m[5] = df[5]*f + loc2*matr[6] - loc1*matr[7];

  adj_m[6] = df[6]*f + loc2*matr[5] - loc3*matr[4];
  adj_m[7] = df[7]*f + loc3*matr[3] - loc1*matr[5];
  adj_m[8] = df[8]*f + loc1*matr[4] - loc2*matr[3];

  loc1 = invW[0][0]*adj_m[0];
  loc2 = invW[1][0]*adj_m[0];
  g_obj[0][0] = -loc1 - loc2;
  g_obj[1][0] =  loc1;
  g_obj[2][0] =  loc2;

  loc1 = invW[0][1]*adj_m[1];
  loc2 = invW[1][1]*adj_m[1];
  g_obj[0][0] -= loc1 + loc2;
  g_obj[1][0] += loc1;
  g_obj[2][0] += loc2;

  loc1 = invW[0][2]*adj_m[2];
  loc2 = invW[1][2]*adj_m[2];
  g_obj[0][0] -= loc1 + loc2;
  g_obj[1][0] += loc1;
  g_obj[2][0] += loc2;

  loc1 = invW[0][0]*adj_m[3];
  loc2 = invW[1][0]*adj_m[3];
  g_obj[0][1] = -loc1 - loc2;
  g_obj[1][1] =  loc1;
  g_obj[2][1] =  loc2;

  loc1 = invW[0][1]*adj_m[4];
  loc2 = invW[1][1]*adj_m[4];
  g_obj[0][1] -= loc1 + loc2;
  g_obj[1][1] += loc1;
  g_obj[2][1] += loc2;

  loc1 = invW[0][2]*adj_m[5];
  loc2 = invW[1][2]*adj_m[5];
  g_obj[0][1] -= loc1 + loc2;
  g_obj[1][1] += loc1;
  g_obj[2][1] += loc2;

  loc1 = invW[0][0]*adj_m[6];
  loc2 = invW[1][0]*adj_m[6];
  g_obj[0][2] = -loc1 - loc2;
  g_obj[1][2] =  loc1;
  g_obj[2][2] =  loc2;

  loc1 = invW[0][1]*adj_m[7];
  loc2 = invW[1][1]*adj_m[7];
  g_obj[0][2] -= loc1 + loc2;
  g_obj[1][2] += loc1;
  g_obj[2][2] += loc2;

  loc1 = invW[0][2]*adj_m[8];
  loc2 = invW[1][2]*adj_m[8];
  g_obj[0][2] -= loc1 + loc2;
  g_obj[1][2] += loc1;
  g_obj[2][2] += loc2;
  return true;
}

inline bool h_fcn_ridft2(double &obj, Vector3D g_obj[3], Matrix3D h_obj[6],
			 const Vector3D x[3], const Vector3D &n,
			 const Matrix3D &invW,
			 const double a, const Exponent b, const Exponent c,
			 const double delta)
{
  static double matr[9], f, t1, t2;
  static double fmat[6], ftmat[6], g, t3, t4;
  static double adj_m[9], df[9], dg[9], loc0, loc1, loc2, loc3, loc4;
  static double A[12], J_A[6], J_B[9], J_C[9], cross;
  static double aux[45];

  /* Calculate M = A*inv(W). */
  f       = x[1][0] - x[0][0];
  g       = x[2][0] - x[0][0];
  t1      = n[0];
  matr[0] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[1] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[2] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][1] - x[0][1];
  g       = x[2][1] - x[0][1];
  t1      = n[1];
  matr[3] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[4] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[5] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][2] - x[0][2];
  g       = x[2][2] - x[0][2];
  t1      = n[2];
  matr[6] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[7] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[8] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  /* Calculate products for M*M' */
  aux[ 0] = matr[0]*matr[0];
  aux[ 1] = matr[0]*matr[1];
  aux[ 2] = matr[0]*matr[2];
  aux[ 3] = matr[0]*matr[3];
  aux[ 4] = matr[0]*matr[4];
  aux[ 5] = matr[0]*matr[5];
  aux[ 6] = matr[0]*matr[6];
  aux[ 7] = matr[0]*matr[7];
  aux[ 8] = matr[0]*matr[8];
  aux[ 9] = matr[1]*matr[1];
  aux[10] = matr[1]*matr[2];
  aux[11] = matr[1]*matr[3];
  aux[12] = matr[1]*matr[4];
  aux[13] = matr[1]*matr[5];
  aux[14] = matr[1]*matr[6];
  aux[15] = matr[1]*matr[7];
  aux[16] = matr[1]*matr[8];
  aux[17] = matr[2]*matr[2];
  aux[18] = matr[2]*matr[3];
  aux[19] = matr[2]*matr[4];
  aux[20] = matr[2]*matr[5];
  aux[21] = matr[2]*matr[6];
  aux[22] = matr[2]*matr[7];
  aux[23] = matr[2]*matr[8];
  aux[24] = matr[3]*matr[3];
  aux[25] = matr[3]*matr[4];
  aux[26] = matr[3]*matr[5];
  aux[27] = matr[3]*matr[6];
  aux[28] = matr[3]*matr[7];
  aux[29] = matr[3]*matr[8];
  aux[30] = matr[4]*matr[4];
  aux[31] = matr[4]*matr[5];
  aux[32] = matr[4]*matr[6];
  aux[33] = matr[4]*matr[7];
  aux[34] = matr[4]*matr[8];
  aux[35] = matr[5]*matr[5];
  aux[36] = matr[5]*matr[6];
  aux[37] = matr[5]*matr[7];
  aux[38] = matr[5]*matr[8];
  aux[39] = matr[6]*matr[6];
  aux[40] = matr[6]*matr[7];
  aux[41] = matr[6]*matr[8];
  aux[42] = matr[7]*matr[7];
  aux[43] = matr[7]*matr[8];
  aux[44] = matr[8]*matr[8];

  /* Calculate det(M). */
  dg[0] = aux[34] - aux[37];
  dg[1] = aux[36] - aux[29];
  dg[2] = aux[28] - aux[32];
  t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  t2 = t1*t1 + 4.0*delta*delta;
  t3 = sqrt(t2);
  g = t1 + t3;
  if (g < MSQ_MIN) { obj = g; return false; }

  fmat[0] = aux[ 0] + aux[ 9] + aux[17] - 1.0;
  fmat[1] = aux[ 3] + aux[12] + aux[20];
  fmat[2] = aux[ 6] + aux[15] + aux[23];

  fmat[3] = aux[24] + aux[30] + aux[35] - 1.0;
  fmat[4] = aux[27] + aux[33] + aux[38];

  fmat[5] = aux[39] + aux[42] + aux[44] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  loc3 = f;
  loc4 = g;

  /* Calculate objective function. */
  obj  = a * pow(f, b) * pow(g, c);

  /* Calculate the derivative of the objective function. */
  t4 = 1 + t1 / t3;

  f = b * obj / f * 4.0;                /* Constant on nabla f */
  g = c * obj / g * t4;                 /* Constant on nabla g */

  df[0] = fmat[0]*matr[0] + fmat[1]*matr[3] + fmat[2]*matr[6];
  df[1] = fmat[0]*matr[1] + fmat[1]*matr[4] + fmat[2]*matr[7];
  df[2] = fmat[0]*matr[2] + fmat[1]*matr[5] + fmat[2]*matr[8];

  df[3] = fmat[1]*matr[0] + fmat[3]*matr[3] + fmat[4]*matr[6];
  df[4] = fmat[1]*matr[1] + fmat[3]*matr[4] + fmat[4]*matr[7];
  df[5] = fmat[1]*matr[2] + fmat[3]*matr[5] + fmat[4]*matr[8];

  df[6] = fmat[2]*matr[0] + fmat[4]*matr[3] + fmat[5]*matr[6];
  df[7] = fmat[2]*matr[1] + fmat[4]*matr[4] + fmat[5]*matr[7];
  df[8] = fmat[2]*matr[2] + fmat[4]*matr[5] + fmat[5]*matr[8];

  dg[3] = aux[22] - aux[16];
  dg[4] = aux[ 8] - aux[21];
  dg[5] = aux[14] - aux[ 7];

  dg[6] = aux[13] - aux[19];
  dg[7] = aux[18] - aux[ 5];
  dg[8] = aux[ 4] - aux[11];

  adj_m[0] = df[0]*f + dg[0]*g;
  adj_m[1] = df[1]*f + dg[1]*g;
  adj_m[2] = df[2]*f + dg[2]*g;
  adj_m[3] = df[3]*f + dg[3]*g;
  adj_m[4] = df[4]*f + dg[4]*g;
  adj_m[5] = df[5]*f + dg[5]*g;
  adj_m[6] = df[6]*f + dg[6]*g;
  adj_m[7] = df[7]*f + dg[7]*g;
  adj_m[8] = df[8]*f + dg[8]*g;

  loc0 = invW[0][0]*adj_m[0];
  loc1 = invW[1][0]*adj_m[0];
  g_obj[0][0] = -loc0 - loc1;
  g_obj[1][0] =  loc0;
  g_obj[2][0] =  loc1;

  loc0 = invW[0][1]*adj_m[1];
  loc1 = invW[1][1]*adj_m[1];
  g_obj[0][0] -= loc0 + loc1;
  g_obj[1][0] += loc0;
  g_obj[2][0] += loc1;

  loc0 = invW[0][2]*adj_m[2];
  loc1 = invW[1][2]*adj_m[2];
  g_obj[0][0] -= loc0 + loc1;
  g_obj[1][0] += loc0;
  g_obj[2][0] += loc1;

  loc0 = invW[0][0]*adj_m[3];
  loc1 = invW[1][0]*adj_m[3];
  g_obj[0][1] = -loc0 - loc1;
  g_obj[1][1] =  loc0;
  g_obj[2][1] =  loc1;

  loc0 = invW[0][1]*adj_m[4];
  loc1 = invW[1][1]*adj_m[4];
  g_obj[0][1] -= loc0 + loc1;
  g_obj[1][1] += loc0;
  g_obj[2][1] += loc1;

  loc0 = invW[0][2]*adj_m[5];
  loc1 = invW[1][2]*adj_m[5];
  g_obj[0][1] -= loc0 + loc1;
  g_obj[1][1] += loc0;
  g_obj[2][1] += loc1;

  loc0 = invW[0][0]*adj_m[6];
  loc1 = invW[1][0]*adj_m[6];
  g_obj[0][2] = -loc0 - loc1;
  g_obj[1][2] =  loc0;
  g_obj[2][2] =  loc1;

  loc0 = invW[0][1]*adj_m[7];
  loc1 = invW[1][1]*adj_m[7];
  g_obj[0][2] -= loc0 + loc1;
  g_obj[1][2] += loc0;
  g_obj[2][2] += loc1;

  loc0 = invW[0][2]*adj_m[8];
  loc1 = invW[1][2]*adj_m[8];
  g_obj[0][2] -= loc0 + loc1;
  g_obj[1][2] += loc0;
  g_obj[2][2] += loc1;

  /* Start of the Hessian evaluation */

  ftmat[0] = aux[ 0] + aux[24] + aux[39];
  ftmat[1] = aux[ 1] + aux[25] + aux[40];
  ftmat[2] = aux[ 2] + aux[26] + aux[41];

  ftmat[3] = aux[ 9] + aux[30] + aux[42];
  ftmat[4] = aux[10] + aux[31] + aux[43];

  ftmat[5] = aux[17] + aux[35] + aux[44];

  loc0 = f;                           /* Constant on nabla^2 f */
  loc1 = g;                           /* Constant on nabla^2 g */
  cross = f * c / loc4 * t4;          /* Constant on nabla g nabla f */
  f = f * (b - 1) / loc3 * 4.0;       /* Constant on nabla f nabla f */
  g = g *((c - 1) * t4 + 4.0*delta*delta / t2) / loc4;
  /* Constant on nabla g nabla g */

  /* First block of rows */
  loc3 = df[0]*f + dg[0]*cross;
  loc4 = dg[0]*g + df[0]*cross;

  J_A[0] = loc3*df[0] + loc4*dg[0];
  J_A[1] = loc3*df[1] + loc4*dg[1];
  J_A[2] = loc3*df[2] + loc4*dg[2];
  J_B[0] = loc3*df[3] + loc4*dg[3];
  J_B[1] = loc3*df[4] + loc4*dg[4];
  J_B[2] = loc3*df[5] + loc4*dg[5];
  J_C[0] = loc3*df[6] + loc4*dg[6];
  J_C[1] = loc3*df[7] + loc4*dg[7];
  J_C[2] = loc3*df[8] + loc4*dg[8];

  loc3 = df[1]*f + dg[1]*cross;
  loc4 = dg[1]*g + df[1]*cross;

  J_A[3] = loc3*df[1] + loc4*dg[1];
  J_A[4] = loc3*df[2] + loc4*dg[2];
  J_B[3] = loc3*df[3] + loc4*dg[3];
  J_B[4] = loc3*df[4] + loc4*dg[4];
  J_B[5] = loc3*df[5] + loc4*dg[5];
  J_C[3] = loc3*df[6] + loc4*dg[6];
  J_C[4] = loc3*df[7] + loc4*dg[7];
  J_C[5] = loc3*df[8] + loc4*dg[8];

  loc3 = df[2]*f + dg[2]*cross;
  loc4 = dg[2]*g + df[2]*cross;

  J_A[5] = loc3*df[2] + loc4*dg[2];
  J_B[6] = loc3*df[3] + loc4*dg[3];
  J_B[7] = loc3*df[4] + loc4*dg[4];
  J_B[8] = loc3*df[5] + loc4*dg[5];
  J_C[6] = loc3*df[6] + loc4*dg[6];
  J_C[7] = loc3*df[7] + loc4*dg[7];
  J_C[8] = loc3*df[8] + loc4*dg[8];

  /* First diagonal block */
  J_A[0] += loc0*(fmat[0] + ftmat[0] + aux[ 0]);
  J_A[1] += loc0*(          ftmat[1] + aux[ 1]);
  J_A[2] += loc0*(          ftmat[2] + aux[ 2]);

  J_A[3] += loc0*(fmat[0] + ftmat[3] + aux[ 9]);
  J_A[4] += loc0*(          ftmat[4] + aux[10]);

  J_A[5] += loc0*(fmat[0] + ftmat[5] + aux[17]);

  loc2 = invW[0][0]+invW[1][0];
  loc3 = invW[0][1]+invW[1][1];
  loc4 = invW[0][2]+invW[1][2];

  A[0]  = -J_A[0]*loc2 - J_A[1]*loc3 - J_A[2]*loc4;
  A[1]  =  J_A[0]*invW[0][0] + J_A[1]*invW[0][1] + J_A[2]*invW[0][2];
  A[2]  =  J_A[0]*invW[1][0] + J_A[1]*invW[1][1] + J_A[2]*invW[1][2];

  A[4]  = -J_A[1]*loc2 - J_A[3]*loc3 - J_A[4]*loc4;
  A[5]  =  J_A[1]*invW[0][0] + J_A[3]*invW[0][1] + J_A[4]*invW[0][2];
  A[6]  =  J_A[1]*invW[1][0] + J_A[3]*invW[1][1] + J_A[4]*invW[1][2];

  A[8]  = -J_A[2]*loc2 - J_A[4]*loc3 - J_A[5]*loc4;
  A[9]  =  J_A[2]*invW[0][0] + J_A[4]*invW[0][1] + J_A[5]*invW[0][2];
  A[10] =  J_A[2]*invW[1][0] + J_A[4]*invW[1][1] + J_A[5]*invW[1][2];

  h_obj[0][0][0] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][0][0] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][0][0] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];

  h_obj[3][0][0] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[4][0][0] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];

  h_obj[5][0][0] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];

  /* First off-diagonal block */
  J_B[0] += loc0*(fmat[1] + aux[3]);
  J_B[1] += loc0*aux[11];
  J_B[2] += loc0*aux[18];

  J_B[3] += loc0*aux[ 4];
  J_B[4] += loc0*(fmat[1] + aux[12]);
  J_B[5] += loc0*aux[19];

  J_B[6] += loc0*aux[ 5];
  J_B[7] += loc0*aux[13];
  J_B[8] += loc0*(fmat[1] + aux[20]);

  loc2 = matr[8]*loc1;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[7]*loc1;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[6]*loc1;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = invW[0][0]+invW[1][0];
  loc3 = invW[0][1]+invW[1][1];
  loc4 = invW[0][2]+invW[1][2];

  A[0]  = -J_B[0]*loc2 - J_B[1]*loc3 - J_B[2]*loc4;
  A[1]  =  J_B[0]*invW[0][0] + J_B[1]*invW[0][1] + J_B[2]*invW[0][2];
  A[2]  =  J_B[0]*invW[1][0] + J_B[1]*invW[1][1] + J_B[2]*invW[1][2];

  A[4]  = -J_B[3]*loc2 - J_B[4]*loc3 - J_B[5]*loc4;
  A[5]  =  J_B[3]*invW[0][0] + J_B[4]*invW[0][1] + J_B[5]*invW[0][2];
  A[6]  =  J_B[3]*invW[1][0] + J_B[4]*invW[1][1] + J_B[5]*invW[1][2];

  A[8]  = -J_B[6]*loc2 - J_B[7]*loc3 - J_B[8]*loc4;
  A[9]  =  J_B[6]*invW[0][0] + J_B[7]*invW[0][1] + J_B[8]*invW[0][2];
  A[10] =  J_B[6]*invW[1][0] + J_B[7]*invW[1][1] + J_B[8]*invW[1][2];

  h_obj[0][0][1] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][1][0] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][1][0] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];

  h_obj[1][0][1] = -A[1]*loc2 - A[5]*loc3 - A[9]*loc4;
  h_obj[3][0][1] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[4][1][0] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];

  h_obj[2][0][1] = -A[2]*loc2 - A[6]*loc3 - A[10]*loc4;
  h_obj[4][0][1] =  A[2]*invW[0][0] + A[6]*invW[0][1] + A[10]*invW[0][2];
  h_obj[5][0][1] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];

  /* Second off-diagonal block */
  J_C[0] += loc0*(fmat[2] + aux[6]);
  J_C[1] += loc0*aux[14];
  J_C[2] += loc0*aux[21];

  J_C[3] += loc0*aux[ 7];
  J_C[4] += loc0*(fmat[2] + aux[15]);
  J_C[5] += loc0*aux[22];

  J_C[6] += loc0*aux[ 8];
  J_C[7] += loc0*aux[16];
  J_C[8] += loc0*(fmat[2] + aux[23]);

  loc2 = matr[5]*loc1;
  J_C[1] -= loc2;
  J_C[3] += loc2;

  loc2 = matr[4]*loc1;
  J_C[2] += loc2;
  J_C[6] -= loc2;

  loc2 = matr[3]*loc1;
  J_C[5] -= loc2;
  J_C[7] += loc2;

  loc2 = invW[0][0]+invW[1][0];
  loc3 = invW[0][1]+invW[1][1];
  loc4 = invW[0][2]+invW[1][2];

  A[0]  = -J_C[0]*loc2 - J_C[1]*loc3 - J_C[2]*loc4;
  A[1]  =  J_C[0]*invW[0][0] + J_C[1]*invW[0][1] + J_C[2]*invW[0][2];
  A[2]  =  J_C[0]*invW[1][0] + J_C[1]*invW[1][1] + J_C[2]*invW[1][2];

  A[4]  = -J_C[3]*loc2 - J_C[4]*loc3 - J_C[5]*loc4;
  A[5]  =  J_C[3]*invW[0][0] + J_C[4]*invW[0][1] + J_C[5]*invW[0][2];
  A[6]  =  J_C[3]*invW[1][0] + J_C[4]*invW[1][1] + J_C[5]*invW[1][2];

  A[8]  = -J_C[6]*loc2 - J_C[7]*loc3 - J_C[8]*loc4;
  A[9]  =  J_C[6]*invW[0][0] + J_C[7]*invW[0][1] + J_C[8]*invW[0][2];
  A[10] =  J_C[6]*invW[1][0] + J_C[7]*invW[1][1] + J_C[8]*invW[1][2];

  h_obj[0][0][2] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][2][0] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][2][0] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];
    
  h_obj[1][0][2] = -A[1]*loc2 - A[5]*loc3 - A[9]*loc4;
  h_obj[3][0][2] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[4][2][0] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];
    
  h_obj[2][0][2] = -A[2]*loc2 - A[6]*loc3 - A[10]*loc4;
  h_obj[4][0][2] =  A[2]*invW[0][0] + A[6]*invW[0][1] + A[10]*invW[0][2];
  h_obj[5][0][2] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];

  /* Second block of rows */
  loc3 = df[3]*f + dg[3]*cross;
  loc4 = dg[3]*g + df[3]*cross;

  J_A[0] = loc3*df[3] + loc4*dg[3];
  J_A[1] = loc3*df[4] + loc4*dg[4];
  J_A[2] = loc3*df[5] + loc4*dg[5];
  J_B[0] = loc3*df[6] + loc4*dg[6];
  J_B[1] = loc3*df[7] + loc4*dg[7];
  J_B[2] = loc3*df[8] + loc4*dg[8];

  loc3 = df[4]*f + dg[4]*cross;
  loc4 = dg[4]*g + df[4]*cross;

  J_A[3] = loc3*df[4] + loc4*dg[4];
  J_A[4] = loc3*df[5] + loc4*dg[5];
  J_B[3] = loc3*df[6] + loc4*dg[6];
  J_B[4] = loc3*df[7] + loc4*dg[7];
  J_B[5] = loc3*df[8] + loc4*dg[8];

  loc3 = df[5]*f + dg[5]*cross;
  loc4 = dg[5]*g + df[5]*cross;

  J_A[5] = loc3*df[5] + loc4*dg[5];
  J_B[6] = loc3*df[6] + loc4*dg[6];
  J_B[7] = loc3*df[7] + loc4*dg[7];
  J_B[8] = loc3*df[8] + loc4*dg[8];

  /* Second diagonal block */
  J_A[0] += loc0*(fmat[3] + ftmat[0] + aux[24]);
  J_A[1] += loc0*(          ftmat[1] + aux[25]);
  J_A[2] += loc0*(          ftmat[2] + aux[26]);

  J_A[3] += loc0*(fmat[3] + ftmat[3] + aux[30]);
  J_A[4] += loc0*(          ftmat[4] + aux[31]);

  J_A[5] += loc0*(fmat[3] + ftmat[5] + aux[35]);

  loc2 = invW[0][0]+invW[1][0];
  loc3 = invW[0][1]+invW[1][1];
  loc4 = invW[0][2]+invW[1][2];

  A[0]  = -J_A[0]*loc2 - J_A[1]*loc3 - J_A[2]*loc4;
  A[1]  =  J_A[0]*invW[0][0] + J_A[1]*invW[0][1] + J_A[2]*invW[0][2];
  A[2]  =  J_A[0]*invW[1][0] + J_A[1]*invW[1][1] + J_A[2]*invW[1][2];

  A[4]  = -J_A[1]*loc2 - J_A[3]*loc3 - J_A[4]*loc4;
  A[5]  =  J_A[1]*invW[0][0] + J_A[3]*invW[0][1] + J_A[4]*invW[0][2];
  A[6]  =  J_A[1]*invW[1][0] + J_A[3]*invW[1][1] + J_A[4]*invW[1][2];

  A[8]  = -J_A[2]*loc2 - J_A[4]*loc3 - J_A[5]*loc4;
  A[9]  =  J_A[2]*invW[0][0] + J_A[4]*invW[0][1] + J_A[5]*invW[0][2];
  A[10] =  J_A[2]*invW[1][0] + J_A[4]*invW[1][1] + J_A[5]*invW[1][2];

  h_obj[0][1][1] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][1][1] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][1][1] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];

  h_obj[3][1][1] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[4][1][1] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];

  h_obj[5][1][1] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];

  /* Third off-diagonal block */
  J_B[0] += loc0*(fmat[4] + aux[27]);
  J_B[1] += loc0*aux[32];
  J_B[2] += loc0*aux[36];

  J_B[3] += loc0*aux[28];
  J_B[4] += loc0*(fmat[4] + aux[33]);
  J_B[5] += loc0*aux[37];

  J_B[6] += loc0*aux[29];
  J_B[7] += loc0*aux[34];
  J_B[8] += loc0*(fmat[4] + aux[38]);

  loc2 = matr[2]*loc1;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[1]*loc1;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[0]*loc1;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = invW[0][0]+invW[1][0];
  loc3 = invW[0][1]+invW[1][1];
  loc4 = invW[0][2]+invW[1][2];

  A[0]  = -J_B[0]*loc2 - J_B[1]*loc3 - J_B[2]*loc4;
  A[1]  =  J_B[0]*invW[0][0] + J_B[1]*invW[0][1] + J_B[2]*invW[0][2];
  A[2]  =  J_B[0]*invW[1][0] + J_B[1]*invW[1][1] + J_B[2]*invW[1][2];

  A[4]  = -J_B[3]*loc2 - J_B[4]*loc3 - J_B[5]*loc4;
  A[5]  =  J_B[3]*invW[0][0] + J_B[4]*invW[0][1] + J_B[5]*invW[0][2];
  A[6]  =  J_B[3]*invW[1][0] + J_B[4]*invW[1][1] + J_B[5]*invW[1][2];

  A[8]  = -J_B[6]*loc2 - J_B[7]*loc3 - J_B[8]*loc4;
  A[9]  =  J_B[6]*invW[0][0] + J_B[7]*invW[0][1] + J_B[8]*invW[0][2];
  A[10] =  J_B[6]*invW[1][0] + J_B[7]*invW[1][1] + J_B[8]*invW[1][2];

  h_obj[0][1][2] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][2][1] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][2][1] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];

  h_obj[1][1][2] = -A[1]*loc2 - A[5]*loc3 - A[9]*loc4;
  h_obj[3][1][2] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[4][2][1] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];

  h_obj[2][1][2] = -A[2]*loc2 - A[6]*loc3 - A[10]*loc4;
  h_obj[4][1][2] =  A[2]*invW[0][0] + A[6]*invW[0][1] + A[10]*invW[0][2];
  h_obj[5][1][2] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];

  /* Third block of rows */
  loc3 = df[6]*f + dg[6]*cross;
  loc4 = dg[6]*g + df[6]*cross;

  J_A[0] = loc3*df[6] + loc4*dg[6];
  J_A[1] = loc3*df[7] + loc4*dg[7];
  J_A[2] = loc3*df[8] + loc4*dg[8];

  loc3 = df[7]*f + dg[7]*cross;
  loc4 = dg[7]*g + df[7]*cross;

  J_A[3] = loc3*df[7] + loc4*dg[7];
  J_A[4] = loc3*df[8] + loc4*dg[8];

  loc3 = df[8]*f + dg[8]*cross;
  loc4 = dg[8]*g + df[8]*cross;

  J_A[5] = loc3*df[8] + loc4*dg[8];

  /* Third diagonal block */
  J_A[0] += loc0*(fmat[5] + ftmat[0] + aux[39]);
  J_A[1] += loc0*(          ftmat[1] + aux[40]);
  J_A[2] += loc0*(          ftmat[2] + aux[41]);

  J_A[3] += loc0*(fmat[5] + ftmat[3] + aux[42]);
  J_A[4] += loc0*(          ftmat[4] + aux[43]);

  J_A[5] += loc0*(fmat[5] + ftmat[5] + aux[44]);

  loc2 = invW[0][0]+invW[1][0];
  loc3 = invW[0][1]+invW[1][1];
  loc4 = invW[0][2]+invW[1][2];

  A[0]  = -J_A[0]*loc2 - J_A[1]*loc3 - J_A[2]*loc4;
  A[1]  =  J_A[0]*invW[0][0] + J_A[1]*invW[0][1] + J_A[2]*invW[0][2];
  A[2]  =  J_A[0]*invW[1][0] + J_A[1]*invW[1][1] + J_A[2]*invW[1][2];

  A[4]  = -J_A[1]*loc2 - J_A[3]*loc3 - J_A[4]*loc4;
  A[5]  =  J_A[1]*invW[0][0] + J_A[3]*invW[0][1] + J_A[4]*invW[0][2];
  A[6]  =  J_A[1]*invW[1][0] + J_A[3]*invW[1][1] + J_A[4]*invW[1][2];

  A[8]  = -J_A[2]*loc2 - J_A[4]*loc3 - J_A[5]*loc4;
  A[9]  =  J_A[2]*invW[0][0] + J_A[4]*invW[0][1] + J_A[5]*invW[0][2];
  A[10] =  J_A[2]*invW[1][0] + J_A[4]*invW[1][1] + J_A[5]*invW[1][2];

  h_obj[0][2][2] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][2][2] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][2][2] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];

  h_obj[3][2][2] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[4][2][2] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];

  h_obj[5][2][2] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];

  /* Complete diagonal blocks */
  h_obj[0].fill_lower_triangle();
  h_obj[3].fill_lower_triangle();
  h_obj[5].fill_lower_triangle();
  return true;
}

inline bool m_fcn_ridft3(double &obj, const Vector3D x[4], 
			 const Matrix3D &invW,
			 const double a, const Exponent b, const Exponent c,
			 const double delta)
{
  static double matr[9], f, t1, t2;
  static double fmat[6], g;

  /* Calculate M = A*inv(W). */
  f       = x[1][0] - x[0][0];
  g       = x[2][0] - x[0][0];
  t1      = x[3][0] - x[0][0];
  matr[0] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[1] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[2] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][1] - x[0][1];
  g       = x[2][1] - x[0][1];
  t1      = x[3][1] - x[0][1];
  matr[3] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[4] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[5] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][2] - x[0][2];
  g       = x[2][2] - x[0][2];
  t1      = x[3][2] - x[0][2];
  matr[6] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[7] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[8] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  /* Calculate det(M). */
  t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
       matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
       matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  t2 = sqrt(t1*t1 + 4.0*delta*delta);
  g = t1 + t2;
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] - 1.0;
  fmat[1] = matr[0]*matr[3] + matr[1]*matr[4] + matr[2]*matr[5];
  fmat[2] = matr[0]*matr[6] + matr[1]*matr[7] + matr[2]*matr[8];

  fmat[3] = matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] - 1.0;
  fmat[4] = matr[3]*matr[6] + matr[4]*matr[7] + matr[5]*matr[8];

  fmat[5] = matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  /* Calculate objective function. */
  obj = a * pow(f, b) * pow(g, c);
  return true;
}

inline bool g_fcn_ridft3(double &obj, Vector3D g_obj[4], const Vector3D x[4],
			 const Matrix3D &invW,
			 const double a, const Exponent b, const Exponent c,
			 const double delta)
{
  static double matr[9], f, t1, t2;
  static double fmat[6], g;
  static double adj_m[9], df[9], loc1, loc2, loc3;

  /* Calculate M = A*inv(W). */
  f       = x[1][0] - x[0][0];
  g       = x[2][0] - x[0][0];
  t1      = x[3][0] - x[0][0];
  matr[0] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[1] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[2] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][1] - x[0][1];
  g       = x[2][1] - x[0][1];
  t1      = x[3][1] - x[0][1];
  matr[3] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[4] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[5] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][2] - x[0][2];
  g       = x[2][2] - x[0][2];
  t1      = x[3][2] - x[0][2];
  matr[6] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[7] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[8] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  /* Calculate det(M). */
  loc1 = matr[4]*matr[8] - matr[5]*matr[7];
  loc2 = matr[5]*matr[6] - matr[3]*matr[8];
  loc3 = matr[3]*matr[7] - matr[4]*matr[6];
  t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;
  t2 = sqrt(t1*t1 + 4.0*delta*delta);
  g = t1 + t2;
  if (g < MSQ_MIN) { obj = g; return false; }

  /* Calculate norm(M). */
  fmat[0] = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] - 1.0;
  fmat[1] = matr[0]*matr[3] + matr[1]*matr[4] + matr[2]*matr[5];
  fmat[2] = matr[0]*matr[6] + matr[1]*matr[7] + matr[2]*matr[8];

  fmat[3] = matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] - 1.0;
  fmat[4] = matr[3]*matr[6] + matr[4]*matr[7] + matr[5]*matr[8];

  fmat[5] = matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  /* Compute objective function and constants.  Note that this */
  /* is now done by cases on $b$.  We only allow $b=1$, $b=2$, */
  /* or $b > 2$ in the computations becasue other values do    */
  /* not have defined Hessians when  $f = 0$.                  */
  
  if (1 == b) { 
    obj = a * f * pow(g, c);
    f = a * pow(g, c) * 2.0;            /* Constant on nabla f */
    g = c * obj / g * (1 + t1 / t2);    /* Constant on nabla g */
  }
  else if (2 == b) {
    obj = a * f * f * pow(g, c);
    f = a * f * pow(g, c) * 4.0;        /* Constant on nabla f */
    g = c * obj / g * (1 + t1 / t2);    /* Constant on nabla g */
  }
  else if (2 < b) {
    if (0 == f) {
      obj = 0;
      f = 0;
      g = 0;
    }
    else {
      obj = a * pow(f, b) * pow(g, c);
      f = b * obj / f * 2.0;            /* Constant on nabla f */
      g = c * obj / g * (1 + t1 / t2);  /* Constant on nabla g */
    }
  }
  else {
    printf("b = %5.4e not allowed in RI_DFT\n", (double)b);
    exit(-1);
  }

  df[0] = fmat[0]*matr[0] + fmat[1]*matr[3] + fmat[2]*matr[6];
  df[1] = fmat[0]*matr[1] + fmat[1]*matr[4] + fmat[2]*matr[7];
  df[2] = fmat[0]*matr[2] + fmat[1]*matr[5] + fmat[2]*matr[8];

  df[3] = fmat[1]*matr[0] + fmat[3]*matr[3] + fmat[4]*matr[6];
  df[4] = fmat[1]*matr[1] + fmat[3]*matr[4] + fmat[4]*matr[7];
  df[5] = fmat[1]*matr[2] + fmat[3]*matr[5] + fmat[4]*matr[8];

  df[6] = fmat[2]*matr[0] + fmat[4]*matr[3] + fmat[5]*matr[6];
  df[7] = fmat[2]*matr[1] + fmat[4]*matr[4] + fmat[5]*matr[7];
  df[8] = fmat[2]*matr[2] + fmat[4]*matr[5] + fmat[5]*matr[8];

  adj_m[0] = df[0]*f + loc1*g;
  adj_m[1] = df[1]*f + loc2*g;
  adj_m[2] = df[2]*f + loc3*g;

  loc1 = matr[0]*g;
  loc2 = matr[1]*g;
  loc3 = matr[2]*g;

  adj_m[3] = df[3]*f + loc3*matr[7] - loc2*matr[8];
  adj_m[4] = df[4]*f + loc1*matr[8] - loc3*matr[6];
  adj_m[5] = df[5]*f + loc2*matr[6] - loc1*matr[7];

  adj_m[6] = df[6]*f + loc2*matr[5] - loc3*matr[4];
  adj_m[7] = df[7]*f + loc3*matr[3] - loc1*matr[5];
  adj_m[8] = df[8]*f + loc1*matr[4] - loc2*matr[3];

  loc1 = invW[0][0]*adj_m[0];
  loc2 = invW[1][0]*adj_m[0];
  loc3 = invW[2][0]*adj_m[0];
  g_obj[0][0] = -loc1 - loc2 - loc3;
  g_obj[1][0] =  loc1;
  g_obj[2][0] =  loc2;
  g_obj[3][0] =  loc3;

  loc1 = invW[0][1]*adj_m[1];
  loc2 = invW[1][1]*adj_m[1];
  loc3 = invW[2][1]*adj_m[1];
  g_obj[0][0] -= loc1 + loc2 + loc3;
  g_obj[1][0] += loc1;
  g_obj[2][0] += loc2;
  g_obj[3][0] += loc3;

  loc1 = invW[0][2]*adj_m[2];
  loc2 = invW[1][2]*adj_m[2];
  loc3 = invW[2][2]*adj_m[2];
  g_obj[0][0] -= loc1 + loc2 + loc3;
  g_obj[1][0] += loc1;
  g_obj[2][0] += loc2;
  g_obj[3][0] += loc3;

  loc1 = invW[0][0]*adj_m[3];
  loc2 = invW[1][0]*adj_m[3];
  loc3 = invW[2][0]*adj_m[3];
  g_obj[0][1] = -loc1 - loc2 - loc3;
  g_obj[1][1] =  loc1;
  g_obj[2][1] =  loc2;
  g_obj[3][1] =  loc3;

  loc1 = invW[0][1]*adj_m[4];
  loc2 = invW[1][1]*adj_m[4];
  loc3 = invW[2][1]*adj_m[4];
  g_obj[0][1] -= loc1 + loc2 + loc3;
  g_obj[1][1] += loc1;
  g_obj[2][1] += loc2;
  g_obj[3][1] += loc3;

  loc1 = invW[0][2]*adj_m[5];
  loc2 = invW[1][2]*adj_m[5];
  loc3 = invW[2][2]*adj_m[5];
  g_obj[0][1] -= loc1 + loc2 + loc3;
  g_obj[1][1] += loc1;
  g_obj[2][1] += loc2;
  g_obj[3][1] += loc3;

  loc1 = invW[0][0]*adj_m[6];
  loc2 = invW[1][0]*adj_m[6];
  loc3 = invW[2][0]*adj_m[6];
  g_obj[0][2] = -loc1 - loc2 - loc3;
  g_obj[1][2] =  loc1;
  g_obj[2][2] =  loc2;
  g_obj[3][2] =  loc3;

  loc1 = invW[0][1]*adj_m[7];
  loc2 = invW[1][1]*adj_m[7];
  loc3 = invW[2][1]*adj_m[7];
  g_obj[0][2] -= loc1 + loc2 + loc3;
  g_obj[1][2] += loc1;
  g_obj[2][2] += loc2;
  g_obj[3][2] += loc3;

  loc1 = invW[0][2]*adj_m[8];
  loc2 = invW[1][2]*adj_m[8];
  loc3 = invW[2][2]*adj_m[8];
  g_obj[0][2] -= loc1 + loc2 + loc3;
  g_obj[1][2] += loc1;
  g_obj[2][2] += loc2;
  g_obj[3][2] += loc3;
  return true;
}

inline bool h_fcn_ridft3(double &obj, Vector3D g_obj[4], Matrix3D h_obj[10],
			 const Vector3D x[4], const Matrix3D &invW,
			 const double a, const Exponent b, const Exponent c,
			 const double delta)
{
  static double matr[9], f, t1, t2;
  static double fmat[6], ftmat[6], g, t3, t4;
  static double adj_m[9], df[9], dg[9], loc0, loc1, loc2, loc3, loc4;
  static double A[12], J_A[6], J_B[9], J_C[9], cross;
  static double aux[45];

  /* Calculate M = A*inv(W). */
  f       = x[1][0] - x[0][0];
  g       = x[2][0] - x[0][0];
  t1      = x[3][0] - x[0][0];
  matr[0] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[1] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[2] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][1] - x[0][1];
  g       = x[2][1] - x[0][1];
  t1      = x[3][1] - x[0][1];
  matr[3] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[4] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[5] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  f       = x[1][2] - x[0][2];
  g       = x[2][2] - x[0][2];
  t1      = x[3][2] - x[0][2];
  matr[6] = f*invW[0][0] + g*invW[1][0] + t1*invW[2][0];
  matr[7] = f*invW[0][1] + g*invW[1][1] + t1*invW[2][1];
  matr[8] = f*invW[0][2] + g*invW[1][2] + t1*invW[2][2];

  /* Calculate products for M*M' */
  aux[ 0] = matr[0]*matr[0];
  aux[ 1] = matr[0]*matr[1];
  aux[ 2] = matr[0]*matr[2];
  aux[ 3] = matr[0]*matr[3];
  aux[ 4] = matr[0]*matr[4];
  aux[ 5] = matr[0]*matr[5];
  aux[ 6] = matr[0]*matr[6];
  aux[ 7] = matr[0]*matr[7];
  aux[ 8] = matr[0]*matr[8];
  aux[ 9] = matr[1]*matr[1];
  aux[10] = matr[1]*matr[2];
  aux[11] = matr[1]*matr[3];
  aux[12] = matr[1]*matr[4];
  aux[13] = matr[1]*matr[5];
  aux[14] = matr[1]*matr[6];
  aux[15] = matr[1]*matr[7];
  aux[16] = matr[1]*matr[8];
  aux[17] = matr[2]*matr[2];
  aux[18] = matr[2]*matr[3];
  aux[19] = matr[2]*matr[4];
  aux[20] = matr[2]*matr[5];
  aux[21] = matr[2]*matr[6];
  aux[22] = matr[2]*matr[7];
  aux[23] = matr[2]*matr[8];
  aux[24] = matr[3]*matr[3];
  aux[25] = matr[3]*matr[4];
  aux[26] = matr[3]*matr[5];
  aux[27] = matr[3]*matr[6];
  aux[28] = matr[3]*matr[7];
  aux[29] = matr[3]*matr[8];
  aux[30] = matr[4]*matr[4];
  aux[31] = matr[4]*matr[5];
  aux[32] = matr[4]*matr[6];
  aux[33] = matr[4]*matr[7];
  aux[34] = matr[4]*matr[8];
  aux[35] = matr[5]*matr[5];
  aux[36] = matr[5]*matr[6];
  aux[37] = matr[5]*matr[7];
  aux[38] = matr[5]*matr[8];
  aux[39] = matr[6]*matr[6];
  aux[40] = matr[6]*matr[7];
  aux[41] = matr[6]*matr[8];
  aux[42] = matr[7]*matr[7];
  aux[43] = matr[7]*matr[8];
  aux[44] = matr[8]*matr[8];

  /* Calculate det(M). */
  dg[0] = aux[34] - aux[37];
  dg[1] = aux[36] - aux[29];
  dg[2] = aux[28] - aux[32];
  t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];
  t2 = t1*t1 + 4.0*delta*delta;
  t3 = sqrt(t2);
  g = t1 + t3;
  if (g < MSQ_MIN) { obj = g; return false; }

  fmat[0] = aux[ 0] + aux[ 9] + aux[17] - 1.0;
  fmat[1] = aux[ 3] + aux[12] + aux[20];
  fmat[2] = aux[ 6] + aux[15] + aux[23];

  fmat[3] = aux[24] + aux[30] + aux[35] - 1.0;
  fmat[4] = aux[27] + aux[33] + aux[38];

  fmat[5] = aux[39] + aux[42] + aux[44] - 1.0;

  f = fmat[0]*fmat[0] + 2.0*fmat[1]*fmat[1] + 2.0*fmat[2]*fmat[2] +
                            fmat[3]*fmat[3] + 2.0*fmat[4]*fmat[4] +
                                                  fmat[5]*fmat[5];

  loc3 = f;
  loc4 = g;

  /* Compute objective function and constants.  Note that this */
  /* is now done by cases on $b$.  We only allow $b=1$, $b=2$, */
  /* or $b > 2$ in the computations becasue other values do    */
  /* not have defined Hessians when  $f = 0$.                  */

  t4 = 1 + t1 / t3;
  if (1 == b) {
    obj = a * f * pow(g, c);
    f = a * pow(g, c) * 2.0;            /* Constant on nabla f */
    g = c * obj / g * t4;               /* Constant on nabla g */
  }
  else if (2 == b) {
    obj = a * f * f * pow(g, c);
    f = a * f * pow(g, c) * 4.0;        /* Constant on nabla f */
    g = c * obj / g * t4;               /* Constant on nabla g */
  }
  else if (2 < b) {
    if (0 == f) {
      obj = 0;
      f = 0;
      g = 0;
    }
    else {
      obj = a * pow(f, b) * pow(g, c);
      f = b * obj / f * 2.0;            /* Constant on nabla f */
      g = c * obj / g * t4;             /* Constant on nabla g */
    }
  }
  else {
    printf("b = %5.4e not allowed in RI_DFT\n", (double)b);
    exit(-1);
  }

  df[0] = fmat[0]*matr[0] + fmat[1]*matr[3] + fmat[2]*matr[6];
  df[1] = fmat[0]*matr[1] + fmat[1]*matr[4] + fmat[2]*matr[7];
  df[2] = fmat[0]*matr[2] + fmat[1]*matr[5] + fmat[2]*matr[8];

  df[3] = fmat[1]*matr[0] + fmat[3]*matr[3] + fmat[4]*matr[6];
  df[4] = fmat[1]*matr[1] + fmat[3]*matr[4] + fmat[4]*matr[7];
  df[5] = fmat[1]*matr[2] + fmat[3]*matr[5] + fmat[4]*matr[8];

  df[6] = fmat[2]*matr[0] + fmat[4]*matr[3] + fmat[5]*matr[6];
  df[7] = fmat[2]*matr[1] + fmat[4]*matr[4] + fmat[5]*matr[7];
  df[8] = fmat[2]*matr[2] + fmat[4]*matr[5] + fmat[5]*matr[8];

  dg[3] = aux[22] - aux[16];
  dg[4] = aux[ 8] - aux[21];
  dg[5] = aux[14] - aux[ 7];

  dg[6] = aux[13] - aux[19];
  dg[7] = aux[18] - aux[ 5];
  dg[8] = aux[ 4] - aux[11];

  adj_m[0] = df[0]*f + dg[0]*g;
  adj_m[1] = df[1]*f + dg[1]*g;
  adj_m[2] = df[2]*f + dg[2]*g;
  adj_m[3] = df[3]*f + dg[3]*g;
  adj_m[4] = df[4]*f + dg[4]*g;
  adj_m[5] = df[5]*f + dg[5]*g;
  adj_m[6] = df[6]*f + dg[6]*g;
  adj_m[7] = df[7]*f + dg[7]*g;
  adj_m[8] = df[8]*f + dg[8]*g;

  loc0 = invW[0][0]*adj_m[0];
  loc1 = invW[1][0]*adj_m[0];
  loc2 = invW[2][0]*adj_m[0];
  g_obj[0][0] = -loc0 - loc1 - loc2;
  g_obj[1][0] =  loc0;
  g_obj[2][0] =  loc1;
  g_obj[3][0] =  loc2;

  loc0 = invW[0][1]*adj_m[1];
  loc1 = invW[1][1]*adj_m[1];
  loc2 = invW[2][1]*adj_m[1];
  g_obj[0][0] -= loc0 + loc1 + loc2;
  g_obj[1][0] += loc0;
  g_obj[2][0] += loc1;
  g_obj[3][0] += loc2;

  loc0 = invW[0][2]*adj_m[2];
  loc1 = invW[1][2]*adj_m[2];
  loc2 = invW[2][2]*adj_m[2];
  g_obj[0][0] -= loc0 + loc1 + loc2;
  g_obj[1][0] += loc0;
  g_obj[2][0] += loc1;
  g_obj[3][0] += loc2;

  loc0 = invW[0][0]*adj_m[3];
  loc1 = invW[1][0]*adj_m[3];
  loc2 = invW[2][0]*adj_m[3];
  g_obj[0][1] = -loc0 - loc1 - loc2;
  g_obj[1][1] =  loc0;
  g_obj[2][1] =  loc1;
  g_obj[3][1] =  loc2;

  loc0 = invW[0][1]*adj_m[4];
  loc1 = invW[1][1]*adj_m[4];
  loc2 = invW[2][1]*adj_m[4];
  g_obj[0][1] -= loc0 + loc1 + loc2;
  g_obj[1][1] += loc0;
  g_obj[2][1] += loc1;
  g_obj[3][1] += loc2;

  loc0 = invW[0][2]*adj_m[5];
  loc1 = invW[1][2]*adj_m[5];
  loc2 = invW[2][2]*adj_m[5];
  g_obj[0][1] -= loc0 + loc1 + loc2;
  g_obj[1][1] += loc0;
  g_obj[2][1] += loc1;
  g_obj[3][1] += loc2;

  loc0 = invW[0][0]*adj_m[6];
  loc1 = invW[1][0]*adj_m[6];
  loc2 = invW[2][0]*adj_m[6];
  g_obj[0][2] = -loc0 - loc1 - loc2;
  g_obj[1][2] =  loc0;
  g_obj[2][2] =  loc1;
  g_obj[3][2] =  loc2;

  loc0 = invW[0][1]*adj_m[7];
  loc1 = invW[1][1]*adj_m[7];
  loc2 = invW[2][1]*adj_m[7];
  g_obj[0][2] -= loc0 + loc1 + loc2;
  g_obj[1][2] += loc0;
  g_obj[2][2] += loc1;
  g_obj[3][2] += loc2;

  loc0 = invW[0][2]*adj_m[8];
  loc1 = invW[1][2]*adj_m[8];
  loc2 = invW[2][2]*adj_m[8];
  g_obj[0][2] -= loc0 + loc1 + loc2;
  g_obj[1][2] += loc0;
  g_obj[2][2] += loc1;
  g_obj[3][2] += loc2;

  /* Start of the Hessian evaluation */

  ftmat[0] = aux[ 0] + aux[24] + aux[39];
  ftmat[1] = aux[ 1] + aux[25] + aux[40];
  ftmat[2] = aux[ 2] + aux[26] + aux[41];

  ftmat[3] = aux[ 9] + aux[30] + aux[42];
  ftmat[4] = aux[10] + aux[31] + aux[43];

  ftmat[5] = aux[17] + aux[35] + aux[44];

  loc0 = f;                           /* Constant on nabla^2 f */
  loc1 = g;                           /* Constant on nabla^2 g */
  if (1 == b) {
    cross = f * c / loc4 * t4;          /* Constant on nabla g nabla f */
    f = 0;                              /* Constant on nabla f nabla f */
    g = g *((c - 1) * t4 + 4.0*delta*delta / t2) / loc4;
                                        /* Constant on nabla g nabla g */
  }
  else if (2 == b) {
    cross = f * c / loc4 * t4;          /* Constant on nabla g nabla f */
    f = a * pow(loc4, c) * 8.0;         /* Constant on nabla f nabla f */
    g = g *((c - 1) * t4 + 4.0*delta*delta / t2) / loc4;
                                        /* Constant on nabla g nabla g */
  }
  else if (2 < b) {
    if (0 ==loc3) {
      cross = 0;
      f = 0;
      g = 0;
    }
    else {
      cross = f * c / loc4 * t4;        /* Constant on nabla g nabla f */
      f = f * (b - 1) / loc3 * 2.0;     /* Constant on nabla f nabla f */
      g = g *((c - 1) * t4 + 4.0*delta*delta / t2) / loc4;
                                        /* Constant on nabla g nabla g */
    }
  }
  else {
    cross = 0;
    f = 0;
    g = 0;
  }

  /* First block of rows */
  loc3 = df[0]*f + dg[0]*cross;
  loc4 = dg[0]*g + df[0]*cross;

  J_A[0] = loc3*df[0] + loc4*dg[0];
  J_A[1] = loc3*df[1] + loc4*dg[1];
  J_A[2] = loc3*df[2] + loc4*dg[2];
  J_B[0] = loc3*df[3] + loc4*dg[3];
  J_B[1] = loc3*df[4] + loc4*dg[4];
  J_B[2] = loc3*df[5] + loc4*dg[5];
  J_C[0] = loc3*df[6] + loc4*dg[6];
  J_C[1] = loc3*df[7] + loc4*dg[7];
  J_C[2] = loc3*df[8] + loc4*dg[8];

  loc3 = df[1]*f + dg[1]*cross;
  loc4 = dg[1]*g + df[1]*cross;

  J_A[3] = loc3*df[1] + loc4*dg[1];
  J_A[4] = loc3*df[2] + loc4*dg[2];
  J_B[3] = loc3*df[3] + loc4*dg[3];
  J_B[4] = loc3*df[4] + loc4*dg[4];
  J_B[5] = loc3*df[5] + loc4*dg[5];
  J_C[3] = loc3*df[6] + loc4*dg[6];
  J_C[4] = loc3*df[7] + loc4*dg[7];
  J_C[5] = loc3*df[8] + loc4*dg[8];

  loc3 = df[2]*f + dg[2]*cross;
  loc4 = dg[2]*g + df[2]*cross;

  J_A[5] = loc3*df[2] + loc4*dg[2];
  J_B[6] = loc3*df[3] + loc4*dg[3];
  J_B[7] = loc3*df[4] + loc4*dg[4];
  J_B[8] = loc3*df[5] + loc4*dg[5];
  J_C[6] = loc3*df[6] + loc4*dg[6];
  J_C[7] = loc3*df[7] + loc4*dg[7];
  J_C[8] = loc3*df[8] + loc4*dg[8];

  /* First diagonal block */
  J_A[0] += loc0*(fmat[0] + ftmat[0] + aux[ 0]);
  J_A[1] += loc0*(          ftmat[1] + aux[ 1]);
  J_A[2] += loc0*(          ftmat[2] + aux[ 2]);

  J_A[3] += loc0*(fmat[0] + ftmat[3] + aux[ 9]);
  J_A[4] += loc0*(          ftmat[4] + aux[10]);

  J_A[5] += loc0*(fmat[0] + ftmat[5] + aux[17]);

  loc2 = invW[0][0]+invW[1][0]+invW[2][0];
  loc3 = invW[0][1]+invW[1][1]+invW[2][1];
  loc4 = invW[0][2]+invW[1][2]+invW[2][2];

  A[0]  = -J_A[0]*loc2 - J_A[1]*loc3 - J_A[2]*loc4;
  A[1]  =  J_A[0]*invW[0][0] + J_A[1]*invW[0][1] + J_A[2]*invW[0][2];
  A[2]  =  J_A[0]*invW[1][0] + J_A[1]*invW[1][1] + J_A[2]*invW[1][2];
  A[3]  =  J_A[0]*invW[2][0] + J_A[1]*invW[2][1] + J_A[2]*invW[2][2];

  A[4]  = -J_A[1]*loc2 - J_A[3]*loc3 - J_A[4]*loc4;
  A[5]  =  J_A[1]*invW[0][0] + J_A[3]*invW[0][1] + J_A[4]*invW[0][2];
  A[6]  =  J_A[1]*invW[1][0] + J_A[3]*invW[1][1] + J_A[4]*invW[1][2];
  A[7]  =  J_A[1]*invW[2][0] + J_A[3]*invW[2][1] + J_A[4]*invW[2][2];

  A[8]  = -J_A[2]*loc2 - J_A[4]*loc3 - J_A[5]*loc4;
  A[9]  =  J_A[2]*invW[0][0] + J_A[4]*invW[0][1] + J_A[5]*invW[0][2];
  A[10] =  J_A[2]*invW[1][0] + J_A[4]*invW[1][1] + J_A[5]*invW[1][2];
  A[11] =  J_A[2]*invW[2][0] + J_A[4]*invW[2][1] + J_A[5]*invW[2][2];

  h_obj[0][0][0] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][0][0] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][0][0] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];
  h_obj[3][0][0] =  A[0]*invW[2][0] + A[4]*invW[2][1] + A[8]*invW[2][2];

  h_obj[4][0][0] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[5][0][0] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];
  h_obj[6][0][0] =  A[1]*invW[2][0] + A[5]*invW[2][1] + A[9]*invW[2][2];

  h_obj[7][0][0] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];
  h_obj[8][0][0] =  A[2]*invW[2][0] + A[6]*invW[2][1] + A[10]*invW[2][2];

  h_obj[9][0][0] =  A[3]*invW[2][0] + A[7]*invW[2][1] + A[11]*invW[2][2];

  /* First off-diagonal block */
  J_B[0] += loc0*(fmat[1] + aux[3]);
  J_B[1] += loc0*aux[11];
  J_B[2] += loc0*aux[18];

  J_B[3] += loc0*aux[ 4];
  J_B[4] += loc0*(fmat[1] + aux[12]);
  J_B[5] += loc0*aux[19];

  J_B[6] += loc0*aux[ 5];
  J_B[7] += loc0*aux[13];
  J_B[8] += loc0*(fmat[1] + aux[20]);

  loc2 = matr[8]*loc1;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[7]*loc1;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[6]*loc1;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = invW[0][0]+invW[1][0]+invW[2][0];
  loc3 = invW[0][1]+invW[1][1]+invW[2][1];
  loc4 = invW[0][2]+invW[1][2]+invW[2][2];

  A[0]  = -J_B[0]*loc2 - J_B[1]*loc3 - J_B[2]*loc4;
  A[1]  =  J_B[0]*invW[0][0] + J_B[1]*invW[0][1] + J_B[2]*invW[0][2];
  A[2]  =  J_B[0]*invW[1][0] + J_B[1]*invW[1][1] + J_B[2]*invW[1][2];
  A[3]  =  J_B[0]*invW[2][0] + J_B[1]*invW[2][1] + J_B[2]*invW[2][2];

  A[4]  = -J_B[3]*loc2 - J_B[4]*loc3 - J_B[5]*loc4;
  A[5]  =  J_B[3]*invW[0][0] + J_B[4]*invW[0][1] + J_B[5]*invW[0][2];
  A[6]  =  J_B[3]*invW[1][0] + J_B[4]*invW[1][1] + J_B[5]*invW[1][2];
  A[7]  =  J_B[3]*invW[2][0] + J_B[4]*invW[2][1] + J_B[5]*invW[2][2];

  A[8]  = -J_B[6]*loc2 - J_B[7]*loc3 - J_B[8]*loc4;
  A[9]  =  J_B[6]*invW[0][0] + J_B[7]*invW[0][1] + J_B[8]*invW[0][2];
  A[10] =  J_B[6]*invW[1][0] + J_B[7]*invW[1][1] + J_B[8]*invW[1][2];
  A[11] =  J_B[6]*invW[2][0] + J_B[7]*invW[2][1] + J_B[8]*invW[2][2];

  h_obj[0][0][1] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][1][0] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][1][0] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];
  h_obj[3][1][0] =  A[0]*invW[2][0] + A[4]*invW[2][1] + A[8]*invW[2][2];

  h_obj[1][0][1] = -A[1]*loc2 - A[5]*loc3 - A[9]*loc4;
  h_obj[4][0][1] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[5][1][0] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];
  h_obj[6][1][0] =  A[1]*invW[2][0] + A[5]*invW[2][1] + A[9]*invW[2][2];

  h_obj[2][0][1] = -A[2]*loc2 - A[6]*loc3 - A[10]*loc4;
  h_obj[5][0][1] =  A[2]*invW[0][0] + A[6]*invW[0][1] + A[10]*invW[0][2];
  h_obj[7][0][1] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];
  h_obj[8][1][0] =  A[2]*invW[2][0] + A[6]*invW[2][1] + A[10]*invW[2][2];

  h_obj[3][0][1] = -A[3]*loc2 - A[7]*loc3 - A[11]*loc4;
  h_obj[6][0][1] =  A[3]*invW[0][0] + A[7]*invW[0][1] + A[11]*invW[0][2];
  h_obj[8][0][1] =  A[3]*invW[1][0] + A[7]*invW[1][1] + A[11]*invW[1][2];
  h_obj[9][0][1] =  A[3]*invW[2][0] + A[7]*invW[2][1] + A[11]*invW[2][2];

  /* Second off-diagonal block */
  J_C[0] += loc0*(fmat[2] + aux[6]);
  J_C[1] += loc0*aux[14];
  J_C[2] += loc0*aux[21];

  J_C[3] += loc0*aux[ 7];
  J_C[4] += loc0*(fmat[2] + aux[15]);
  J_C[5] += loc0*aux[22];

  J_C[6] += loc0*aux[ 8];
  J_C[7] += loc0*aux[16];
  J_C[8] += loc0*(fmat[2] + aux[23]);

  loc2 = matr[5]*loc1;
  J_C[1] -= loc2;
  J_C[3] += loc2;

  loc2 = matr[4]*loc1;
  J_C[2] += loc2;
  J_C[6] -= loc2;

  loc2 = matr[3]*loc1;
  J_C[5] -= loc2;
  J_C[7] += loc2;

  loc2 = invW[0][0]+invW[1][0]+invW[2][0];
  loc3 = invW[0][1]+invW[1][1]+invW[2][1];
  loc4 = invW[0][2]+invW[1][2]+invW[2][2];

  A[0]  = -J_C[0]*loc2 - J_C[1]*loc3 - J_C[2]*loc4;
  A[1]  =  J_C[0]*invW[0][0] + J_C[1]*invW[0][1] + J_C[2]*invW[0][2];
  A[2]  =  J_C[0]*invW[1][0] + J_C[1]*invW[1][1] + J_C[2]*invW[1][2];
  A[3]  =  J_C[0]*invW[2][0] + J_C[1]*invW[2][1] + J_C[2]*invW[2][2];

  A[4]  = -J_C[3]*loc2 - J_C[4]*loc3 - J_C[5]*loc4;
  A[5]  =  J_C[3]*invW[0][0] + J_C[4]*invW[0][1] + J_C[5]*invW[0][2];
  A[6]  =  J_C[3]*invW[1][0] + J_C[4]*invW[1][1] + J_C[5]*invW[1][2];
  A[7]  =  J_C[3]*invW[2][0] + J_C[4]*invW[2][1] + J_C[5]*invW[2][2];

  A[8]  = -J_C[6]*loc2 - J_C[7]*loc3 - J_C[8]*loc4;
  A[9]  =  J_C[6]*invW[0][0] + J_C[7]*invW[0][1] + J_C[8]*invW[0][2];
  A[10] =  J_C[6]*invW[1][0] + J_C[7]*invW[1][1] + J_C[8]*invW[1][2];
  A[11] =  J_C[6]*invW[2][0] + J_C[7]*invW[2][1] + J_C[8]*invW[2][2];

  h_obj[0][0][2] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][2][0] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][2][0] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];
  h_obj[3][2][0] =  A[0]*invW[2][0] + A[4]*invW[2][1] + A[8]*invW[2][2];
    
  h_obj[1][0][2] = -A[1]*loc2 - A[5]*loc3 - A[9]*loc4;
  h_obj[4][0][2] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[5][2][0] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];
  h_obj[6][2][0] =  A[1]*invW[2][0] + A[5]*invW[2][1] + A[9]*invW[2][2];
    
  h_obj[2][0][2] = -A[2]*loc2 - A[6]*loc3 - A[10]*loc4;
  h_obj[5][0][2] =  A[2]*invW[0][0] + A[6]*invW[0][1] + A[10]*invW[0][2];
  h_obj[7][0][2] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];
  h_obj[8][2][0] =  A[2]*invW[2][0] + A[6]*invW[2][1] + A[10]*invW[2][2];

  h_obj[3][0][2] = -A[3]*loc2 - A[7]*loc3 - A[11]*loc4;
  h_obj[6][0][2] =  A[3]*invW[0][0] + A[7]*invW[0][1] + A[11]*invW[0][2];
  h_obj[8][0][2] =  A[3]*invW[1][0] + A[7]*invW[1][1] + A[11]*invW[1][2];
  h_obj[9][0][2] =  A[3]*invW[2][0] + A[7]*invW[2][1] + A[11]*invW[2][2];

  /* Second block of rows */
  loc3 = df[3]*f + dg[3]*cross;
  loc4 = dg[3]*g + df[3]*cross;

  J_A[0] = loc3*df[3] + loc4*dg[3];
  J_A[1] = loc3*df[4] + loc4*dg[4];
  J_A[2] = loc3*df[5] + loc4*dg[5];
  J_B[0] = loc3*df[6] + loc4*dg[6];
  J_B[1] = loc3*df[7] + loc4*dg[7];
  J_B[2] = loc3*df[8] + loc4*dg[8];

  loc3 = df[4]*f + dg[4]*cross;
  loc4 = dg[4]*g + df[4]*cross;

  J_A[3] = loc3*df[4] + loc4*dg[4];
  J_A[4] = loc3*df[5] + loc4*dg[5];
  J_B[3] = loc3*df[6] + loc4*dg[6];
  J_B[4] = loc3*df[7] + loc4*dg[7];
  J_B[5] = loc3*df[8] + loc4*dg[8];

  loc3 = df[5]*f + dg[5]*cross;
  loc4 = dg[5]*g + df[5]*cross;

  J_A[5] = loc3*df[5] + loc4*dg[5];
  J_B[6] = loc3*df[6] + loc4*dg[6];
  J_B[7] = loc3*df[7] + loc4*dg[7];
  J_B[8] = loc3*df[8] + loc4*dg[8];

  /* Second diagonal block */
  J_A[0] += loc0*(fmat[3] + ftmat[0] + aux[24]);
  J_A[1] += loc0*(          ftmat[1] + aux[25]);
  J_A[2] += loc0*(          ftmat[2] + aux[26]);

  J_A[3] += loc0*(fmat[3] + ftmat[3] + aux[30]);
  J_A[4] += loc0*(          ftmat[4] + aux[31]);

  J_A[5] += loc0*(fmat[3] + ftmat[5] + aux[35]);

  loc2 = invW[0][0]+invW[1][0]+invW[2][0];
  loc3 = invW[0][1]+invW[1][1]+invW[2][1];
  loc4 = invW[0][2]+invW[1][2]+invW[2][2];

  A[0]  = -J_A[0]*loc2 - J_A[1]*loc3 - J_A[2]*loc4;
  A[1]  =  J_A[0]*invW[0][0] + J_A[1]*invW[0][1] + J_A[2]*invW[0][2];
  A[2]  =  J_A[0]*invW[1][0] + J_A[1]*invW[1][1] + J_A[2]*invW[1][2];
  A[3]  =  J_A[0]*invW[2][0] + J_A[1]*invW[2][1] + J_A[2]*invW[2][2];

  A[4]  = -J_A[1]*loc2 - J_A[3]*loc3 - J_A[4]*loc4;
  A[5]  =  J_A[1]*invW[0][0] + J_A[3]*invW[0][1] + J_A[4]*invW[0][2];
  A[6]  =  J_A[1]*invW[1][0] + J_A[3]*invW[1][1] + J_A[4]*invW[1][2];
  A[7]  =  J_A[1]*invW[2][0] + J_A[3]*invW[2][1] + J_A[4]*invW[2][2];

  A[8]  = -J_A[2]*loc2 - J_A[4]*loc3 - J_A[5]*loc4;
  A[9]  =  J_A[2]*invW[0][0] + J_A[4]*invW[0][1] + J_A[5]*invW[0][2];
  A[10] =  J_A[2]*invW[1][0] + J_A[4]*invW[1][1] + J_A[5]*invW[1][2];
  A[11] =  J_A[2]*invW[2][0] + J_A[4]*invW[2][1] + J_A[5]*invW[2][2];

  h_obj[0][1][1] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][1][1] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][1][1] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];
  h_obj[3][1][1] =  A[0]*invW[2][0] + A[4]*invW[2][1] + A[8]*invW[2][2];

  h_obj[4][1][1] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[5][1][1] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];
  h_obj[6][1][1] =  A[1]*invW[2][0] + A[5]*invW[2][1] + A[9]*invW[2][2];

  h_obj[7][1][1] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];
  h_obj[8][1][1] =  A[2]*invW[2][0] + A[6]*invW[2][1] + A[10]*invW[2][2];

  h_obj[9][1][1] =  A[3]*invW[2][0] + A[7]*invW[2][1] + A[11]*invW[2][2];

  /* Third off-diagonal block */
  J_B[0] += loc0*(fmat[4] + aux[27]);
  J_B[1] += loc0*aux[32];
  J_B[2] += loc0*aux[36];

  J_B[3] += loc0*aux[28];
  J_B[4] += loc0*(fmat[4] + aux[33]);
  J_B[5] += loc0*aux[37];

  J_B[6] += loc0*aux[29];
  J_B[7] += loc0*aux[34];
  J_B[8] += loc0*(fmat[4] + aux[38]);

  loc2 = matr[2]*loc1;
  J_B[1] += loc2;
  J_B[3] -= loc2;

  loc2 = matr[1]*loc1;
  J_B[2] -= loc2;
  J_B[6] += loc2;

  loc2 = matr[0]*loc1;
  J_B[5] += loc2;
  J_B[7] -= loc2;

  loc2 = invW[0][0]+invW[1][0]+invW[2][0];
  loc3 = invW[0][1]+invW[1][1]+invW[2][1];
  loc4 = invW[0][2]+invW[1][2]+invW[2][2];

  A[0]  = -J_B[0]*loc2 - J_B[1]*loc3 - J_B[2]*loc4;
  A[1]  =  J_B[0]*invW[0][0] + J_B[1]*invW[0][1] + J_B[2]*invW[0][2];
  A[2]  =  J_B[0]*invW[1][0] + J_B[1]*invW[1][1] + J_B[2]*invW[1][2];
  A[3]  =  J_B[0]*invW[2][0] + J_B[1]*invW[2][1] + J_B[2]*invW[2][2];

  A[4]  = -J_B[3]*loc2 - J_B[4]*loc3 - J_B[5]*loc4;
  A[5]  =  J_B[3]*invW[0][0] + J_B[4]*invW[0][1] + J_B[5]*invW[0][2];
  A[6]  =  J_B[3]*invW[1][0] + J_B[4]*invW[1][1] + J_B[5]*invW[1][2];
  A[7]  =  J_B[3]*invW[2][0] + J_B[4]*invW[2][1] + J_B[5]*invW[2][2];

  A[8]  = -J_B[6]*loc2 - J_B[7]*loc3 - J_B[8]*loc4;
  A[9]  =  J_B[6]*invW[0][0] + J_B[7]*invW[0][1] + J_B[8]*invW[0][2];
  A[10] =  J_B[6]*invW[1][0] + J_B[7]*invW[1][1] + J_B[8]*invW[1][2];
  A[11] =  J_B[6]*invW[2][0] + J_B[7]*invW[2][1] + J_B[8]*invW[2][2];

  h_obj[0][1][2] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][2][1] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][2][1] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];
  h_obj[3][2][1] =  A[0]*invW[2][0] + A[4]*invW[2][1] + A[8]*invW[2][2];

  h_obj[1][1][2] = -A[1]*loc2 - A[5]*loc3 - A[9]*loc4;
  h_obj[4][1][2] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[5][2][1] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];
  h_obj[6][2][1] =  A[1]*invW[2][0] + A[5]*invW[2][1] + A[9]*invW[2][2];

  h_obj[2][1][2] = -A[2]*loc2 - A[6]*loc3 - A[10]*loc4;
  h_obj[5][1][2] =  A[2]*invW[0][0] + A[6]*invW[0][1] + A[10]*invW[0][2];
  h_obj[7][1][2] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];
  h_obj[8][2][1] =  A[2]*invW[2][0] + A[6]*invW[2][1] + A[10]*invW[2][2];

  h_obj[3][1][2] = -A[3]*loc2 - A[7]*loc3 - A[11]*loc4;
  h_obj[6][1][2] =  A[3]*invW[0][0] + A[7]*invW[0][1] + A[11]*invW[0][2];
  h_obj[8][1][2] =  A[3]*invW[1][0] + A[7]*invW[1][1] + A[11]*invW[1][2];
  h_obj[9][1][2] =  A[3]*invW[2][0] + A[7]*invW[2][1] + A[11]*invW[2][2];

  /* Third block of rows */
  loc3 = df[6]*f + dg[6]*cross;
  loc4 = dg[6]*g + df[6]*cross;

  J_A[0] = loc3*df[6] + loc4*dg[6];
  J_A[1] = loc3*df[7] + loc4*dg[7];
  J_A[2] = loc3*df[8] + loc4*dg[8];

  loc3 = df[7]*f + dg[7]*cross;
  loc4 = dg[7]*g + df[7]*cross;

  J_A[3] = loc3*df[7] + loc4*dg[7];
  J_A[4] = loc3*df[8] + loc4*dg[8];

  loc3 = df[8]*f + dg[8]*cross;
  loc4 = dg[8]*g + df[8]*cross;

  J_A[5] = loc3*df[8] + loc4*dg[8];

  /* Third diagonal block */
  J_A[0] += loc0*(fmat[5] + ftmat[0] + aux[39]);
  J_A[1] += loc0*(          ftmat[1] + aux[40]);
  J_A[2] += loc0*(          ftmat[2] + aux[41]);

  J_A[3] += loc0*(fmat[5] + ftmat[3] + aux[42]);
  J_A[4] += loc0*(          ftmat[4] + aux[43]);

  J_A[5] += loc0*(fmat[5] + ftmat[5] + aux[44]);

  loc2 = invW[0][0]+invW[1][0]+invW[2][0];
  loc3 = invW[0][1]+invW[1][1]+invW[2][1];
  loc4 = invW[0][2]+invW[1][2]+invW[2][2];

  A[0]  = -J_A[0]*loc2 - J_A[1]*loc3 - J_A[2]*loc4;
  A[1]  =  J_A[0]*invW[0][0] + J_A[1]*invW[0][1] + J_A[2]*invW[0][2];
  A[2]  =  J_A[0]*invW[1][0] + J_A[1]*invW[1][1] + J_A[2]*invW[1][2];
  A[3]  =  J_A[0]*invW[2][0] + J_A[1]*invW[2][1] + J_A[2]*invW[2][2];

  A[4]  = -J_A[1]*loc2 - J_A[3]*loc3 - J_A[4]*loc4;
  A[5]  =  J_A[1]*invW[0][0] + J_A[3]*invW[0][1] + J_A[4]*invW[0][2];
  A[6]  =  J_A[1]*invW[1][0] + J_A[3]*invW[1][1] + J_A[4]*invW[1][2];
  A[7]  =  J_A[1]*invW[2][0] + J_A[3]*invW[2][1] + J_A[4]*invW[2][2];

  A[8]  = -J_A[2]*loc2 - J_A[4]*loc3 - J_A[5]*loc4;
  A[9]  =  J_A[2]*invW[0][0] + J_A[4]*invW[0][1] + J_A[5]*invW[0][2];
  A[10] =  J_A[2]*invW[1][0] + J_A[4]*invW[1][1] + J_A[5]*invW[1][2];
  A[11] =  J_A[2]*invW[2][0] + J_A[4]*invW[2][1] + J_A[5]*invW[2][2];

  h_obj[0][2][2] = -A[0]*loc2 - A[4]*loc3 - A[8]*loc4;
  h_obj[1][2][2] =  A[0]*invW[0][0] + A[4]*invW[0][1] + A[8]*invW[0][2];
  h_obj[2][2][2] =  A[0]*invW[1][0] + A[4]*invW[1][1] + A[8]*invW[1][2];
  h_obj[3][2][2] =  A[0]*invW[2][0] + A[4]*invW[2][1] + A[8]*invW[2][2];

  h_obj[4][2][2] =  A[1]*invW[0][0] + A[5]*invW[0][1] + A[9]*invW[0][2];
  h_obj[5][2][2] =  A[1]*invW[1][0] + A[5]*invW[1][1] + A[9]*invW[1][2];
  h_obj[6][2][2] =  A[1]*invW[2][0] + A[5]*invW[2][1] + A[9]*invW[2][2];

  h_obj[7][2][2] =  A[2]*invW[1][0] + A[6]*invW[1][1] + A[10]*invW[1][2];
  h_obj[8][2][2] =  A[2]*invW[2][0] + A[6]*invW[2][1] + A[10]*invW[2][2];

  h_obj[9][2][2] =  A[3]*invW[2][0] + A[7]*invW[2][1] + A[11]*invW[2][2];

  /* Complete diagonal blocks */
  h_obj[0].fill_lower_triangle();
  h_obj[4].fill_lower_triangle();
  h_obj[7].fill_lower_triangle();
  h_obj[9].fill_lower_triangle();
  return true;
}

#if 0
bool RI_DFT::evaluate_element(PatchData& pd,
                              MsqMeshEntity* element,
                              double& value, MsqError &err)
{
  Matrix3D T[MSQ_MAX_NUM_VERT_PER_ENT];
  double c_k[MSQ_MAX_NUM_VERT_PER_ENT];
  double dft[MSQ_MAX_NUM_VERT_PER_ENT];
  bool return_flag;
  double h, tau;
    
  size_t num_T = element->vertex_count();
  compute_T_matrices(*element, pd, T, num_T, c_k, err); MSQ_ERRZERO(err);

  const double id[] = {1., 0., 0.,  0., 1., 0.,  0., 0., 1.};
  const Matrix3D I(id);
  Matrix3D TT;
  for (size_t i=0; i<num_T; ++i) {
    tau = det(T[i]);
    TT = transpose(T[i]);
    TT = TT * T[i];
    TT -= I; 
    dft[i] = .5 * Frobenius_2(TT);
    return_flag = get_barrier_function(pd, tau, h, err); MSQ_ERRZERO(err);
    dft[i] /= pow(h, 4.0/3.0);
  }
    
  value = weighted_average_metrics(c_k, dft, num_T, err); MSQ_ERRZERO(err);
    
  return return_flag;
}
#endif

bool RI_DFT::evaluate_element(PatchData& pd,
			      MsqMeshEntity* e,
			      double& m, 
			      MsqError &err)
{
  // Only works with the weighted average
  MsqError   mErr;
  MsqVertex *vertices = pd.get_vertex_array(err); MSQ_ERRZERO(err);

  EntityTopology topo = e->get_element_type();

  const size_t nv = e->vertex_count();
  const size_t *v_i = e->get_vertex_index_array();

  size_t idx = pd.get_element_index(e);
  const TargetMatrix *W = pd.targetMatrices.get_element_corner_tags(&pd, idx, err );
  MSQ_ERRZERO(err);

  // Initialize constants for the metric
  const double delta = pd.get_barrier_delta(err); MSQ_ERRZERO(err);

  const int tetInd[4][4] = {{0, 1, 2, 3}, {1, 2, 0, 3},
			    {2, 0, 1, 3}, {3, 2, 1, 0}};
  const int hexInd[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
			    {2, 3, 1, 6}, {3, 0, 2, 7},
			    {4, 7, 5, 0}, {5, 4, 6, 1},
			    {6, 5, 7, 2}, {7, 6, 4, 3}};

  // Variables used for computing the metric
  double   mMetric;		// Metric value
  bool     mValid;		// Validity of the metric
  int      i, j;

  m = 0.0;
  switch(topo) {
  case TRIANGLE:

    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }

      pd.get_domain_normal_at_vertex(v_i[tetInd[i][0]], true, mNormal, mErr); 

      if (mErr || mNormal.length() == 0) {
        mNormal = (vertices[v_i[tetInd[i][1]]] - vertices[v_i[tetInd[i][0]]]) *
                  (vertices[v_i[tetInd[i][2]]] - vertices[v_i[tetInd[i][0]]]);
        mNormal.normalize();
      }
      mNormal *= MSQ_3RT_2_OVER_6RT_3;

      inv(invW, W[i]);

      mValid = m_fcn_ridft2(mMetric, mCoords, mNormal, invW, a, b, c, delta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;
    }

    m *= MSQ_ONE_THIRD;
    break;

  case QUADRILATERAL:

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      pd.get_domain_normal_at_vertex(v_i[hexInd[i][0]], true, mNormal, mErr);
      if (mErr || mNormal.length() == 0) {
        mNormal = (vertices[v_i[hexInd[i][1]]] - vertices[v_i[hexInd[i][0]]]) *
                  (vertices[v_i[hexInd[i][2]]] - vertices[v_i[hexInd[i][0]]]);
        mNormal.normalize();
      }

      inv(invW, W[i]);

      mValid = m_fcn_ridft2(mMetric, mCoords, mNormal, invW, a, b, c, delta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;
    }

    m *= 0.25;
    break;

  case TETRAHEDRON:

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }

      inv(invW, W[i]);

      mValid = m_fcn_ridft3(mMetric, mCoords, invW, a, b, c, delta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;
    }

    m *= 0.25;
    break;

  case HEXAHEDRON:

    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      inv(invW, W[i]);

      mValid = m_fcn_ridft3(mMetric, mCoords, invW, a, b, c, delta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;
    }

    m *= 0.125;
    break;

  default:
    MSQ_SETERR(err)("element type not implemented.", MsqError::NOT_IMPLEMENTED);
    return false;
  }

  return true;
}

bool RI_DFT::compute_element_analytical_gradient(PatchData &pd,
						 MsqMeshEntity *e,
						 MsqVertex *fv[], 
						 Vector3D g[],
						 int nfv, 
						 double &m,
						 MsqError &err)
{
  // Only works with the weighted average

  MsqVertex *vertices = pd.get_vertex_array(err); MSQ_ERRZERO(err);

  EntityTopology topo = e->get_element_type();

  const size_t nv = e->vertex_count();
  const size_t *v_i = e->get_vertex_index_array();

  size_t idx = pd.get_element_index(e);
  const TargetMatrix *W = pd.targetMatrices.get_element_corner_tags(&pd, idx, err );
  MSQ_ERRZERO(err);

  // Initialize constants for the metric
  const double delta = pd.get_barrier_delta(err); MSQ_ERRZERO(err);

  const int tetInd[4][4] = {{0, 1, 2, 3}, {1, 2, 0, 3},
			    {2, 0, 1, 3}, {3, 2, 1, 0}};
  const int hexInd[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
			    {2, 3, 1, 6}, {3, 0, 2, 7},
			    {4, 7, 5, 0}, {5, 4, 6, 1},
			    {6, 5, 7, 2}, {7, 6, 4, 3}};

  // Variables used for computing the metric
  double   mMetric;		// Metric value
  bool     mValid;		// Validity of the metric
  int      i, j;

  m = 0.0;
  switch(topo) {
  case TRIANGLE:
    assert(3 == nv);

#ifndef ANALYTIC
    mValid = compute_element_numerical_gradient(pd, e, fv, g, nfv, m, err);
    return !MSQ_CHKERR(err) && mValid;
#else

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constaint, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    for (i = 0; i < 3; ++i) {
      mAccGrads[i] = 0.0;
    }

    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }
      
      pd.get_domain_normal_at_vertex(v_i[tetInd[i][0]], true, mNormal, err);
      MSQ_ERRZERO(err);
      mNormal *= MSQ_3RT_2_OVER_6RT_3);

      inv(invW, W[i]);

      mValid = g_fcn_ridft2(mMetric, mGrads, mCoords, mNormal,
			    invW, a, b, c, delta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;

      for (j = 0; j < 3; ++j) {
	mAccGrads[tetInd[i][j]] += W[i].get_cK() * mGrads[j];
      }
    }

    m *= MSQ_ONE_THIRD;
    for (i = 0; i < 3; ++i) {
      mAccGrads[i] *= MSQ_ONE_THIRD;
    }

    // This is not very efficient, but is one way to select correct gradients.
    // For gradients, info is returned only for free vertices, in the order
    // of fv[].

    for (i = 0; i < 3; ++i) {
      for (j = 0; j < nfv; ++j) {
        if (vertices + v_i[i] == fv[j]) {
          g[j] = mAccGrads[i];
        }
      }
    }
#endif

    break;

  case QUADRILATERAL:

#ifndef ANALYTIC
    mValid = compute_element_numerical_gradient(pd, e, fv, g, nfv, m, err);
    return !MSQ_CHKERR(err) && mValid;
#else

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constaint, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    for (i = 0; i < 4; ++i) {
      mAccGrads[i] = 0.0;
    }

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      pd.get_domain_normal_at_vertex(v_i[hexInd[i][0]], true, mNormal, err);
      MSQ_ERRZERO(err);

      inv(invW, W[i]);

      mValid = g_fcn_ridft2(mMetric, mGrads, mCoords, mNormal,
			    invW, a, b, c, delta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;

      for (j = 0; j < 3; ++j) {
	mAccGrads[hexInd[i][j]] += W[i].get_cK() * mGrads[j];
      }
    }

    m *= 0.25;
    for (i = 0; i < 4; ++i) {
      mAccGrads[i] *= 0.25;
    }

    // This is not very efficient, but is one way to select correct gradients
    // For gradients, info is returned only for free vertices, in the order 
    // of fv[].

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < nfv; ++j) {
        if (vertices + v_i[i] == fv[j]) {
          g[j] = mAccGrads[i];
        }
      }
    }
#endif

    break;

  case TETRAHEDRON:

    for (i = 0; i < 4; ++i) {
      mAccGrads[i] = 0.0;
    }

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }

      inv(invW, W[i]);

      mValid = g_fcn_ridft3(mMetric, mGrads, mCoords, invW, a, b, c, delta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;

      for (j = 0; j < 4; ++j) {
	mAccGrads[tetInd[i][j]] += W[i].get_cK() * mGrads[j];
      }
    }

    m *= 0.25;
    for (i = 0; i < 4; ++i) {
      mAccGrads[i] *= 0.25;
    }

    // This is not very efficient, but is one way to select correct gradients.
    // For gradients, info is returned only for free vertices, in the order
    // of fv[].

    for (i = 0; i < 4; ++i) {
      for (size_t k = 0; k < nv; ++k) {
        if (vertices + v_i[i] == fv[k]) {
          g[k] = mAccGrads[i];
        }
      }
    }
    break;

  case HEXAHEDRON:

    for (i = 0; i < 8; ++i) {
      mAccGrads[i] = 0.0;
    }

    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      inv(invW, W[i]);

      mValid = g_fcn_ridft3(mMetric, mGrads, mCoords, invW, a, b, c, delta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;

      for (j = 0; j < 4; ++j) {
	mAccGrads[hexInd[i][j]] += W[i].get_cK() * mGrads[j];
      }
    }

    m *= 0.125;
    for (i = 0; i < 8; ++i) {
      mAccGrads[i] *= 0.125;
    }

    // This is not very efficient, but is one way to select correct gradients
    // For gradients, info is returned only for free vertices, in the order 
    // of fv[].

    for (i = 0; i < 8; ++i) {
      for (j = 0; j < nfv; ++j) {
        if (vertices + v_i[i] == fv[j]) {
          g[j] = mAccGrads[i];
        }
      }
    }
    break;

  default:
    MSQ_SETERR(err)("element type not implemented.",MsqError::NOT_IMPLEMENTED);
    return false;
  }

  return true;
}

bool RI_DFT::compute_element_analytical_hessian(PatchData &pd,
						MsqMeshEntity *e,
						MsqVertex *fv[], 
						Vector3D g[],
						Matrix3D h[],
						int nfv, 
						double &m,
						MsqError &err)
{
  // Only works with the weighted average

  MsqVertex *vertices = pd.get_vertex_array(err);  MSQ_ERRZERO(err);

  EntityTopology topo = e->get_element_type();

  const size_t nv = e->vertex_count();
  const size_t *v_i = e->get_vertex_index_array();

  size_t idx = pd.get_element_index(e);
  const TargetMatrix *W = pd.targetMatrices.get_element_corner_tags(&pd, idx, err );
  MSQ_ERRZERO(err);

  // Initialize constants for the metric
  const double delta = pd.get_barrier_delta(err); MSQ_ERRZERO(err);

  const int tetInd[4][4] = {{0, 1, 2, 3}, {1, 2, 0, 3},
			    {2, 0, 1, 3}, {3, 2, 1, 0}};
  const int hexInd[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
			    {2, 3, 1, 6}, {3, 0, 2, 7},
			    {4, 7, 5, 0}, {5, 4, 6, 1},
			    {6, 5, 7, 2}, {7, 6, 4, 3}};

  // Variables used for computing the metric
  double   mMetric;		// Metric value
  bool     mValid;		// Validity of the metric
  int      i, j, k, l;
  int      row, col, loc;

  m = 0.0;
  switch(topo) {
  case TRIANGLE:
    assert(3 == nv);

#ifndef ANALYTIC
    mValid = compute_element_numerical_hessian(pd, e, fv, g, h, nfv, m, err);
    return !MSQ_CHKERR(err) && mValid;
#else

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constaint, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    // Zero out the hessian and gradient vector
    for (i = 0; i < 3; ++i) {
      g[i] = 0.0;
    }

    for (i = 0; i < 6; ++i) {
      h[i].zero();
    }

    // Compute the metric and sum them together
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }

      pd.get_domain_normal_at_vertex(v_i[tetInd[i][0]], true, mNormal, err);
      MSQ_ERRZERO(err);
      mNormal *= MSQ_3RT_2_OVER_6RT_3;

      inv(invW, W[i]);

      mValid = h_fcn_ridft2(mMetric, mGrads, mHessians, mCoords, mNormal,
			    invW, a, b, c, delta);
      if (!mValid) return false;

      m += W[i].get_cK() * mMetric;

      for (j = 0; j < 3; ++j) {
	g[tetInd[i][j]] += W[i].get_cK() * mGrads[j];
      }

      l = 0;
      for (j = 0; j < 3; ++j) {
	for (k = j; k < 3; ++k) {
	  row = tetInd[i][j];
	  col = tetInd[i][k];

	  if (row <= col) {
	    loc = 3*row - (row*(row+1)/2) + col;
	    h[loc] += W[i].get_cK() * mHessians[l];
	  }
	  else {
	    loc = 3*col - (col*(col+1)/2) + row;
	    h[loc] += W[i].get_cK() * transpose(mHessians[l]);
	  }
	  ++l;
	}
      }
    }

    m *= MSQ_ONE_THIRD;
    for (i = 0; i < 3; ++i) {
      g[i] *= MSQ_ONE_THIRD;
    }

    for (i = 0; i < 6; ++i) {
      h[i] *= MSQ_ONE_THIRD;
    }

    // zero out fixed elements of g
    j = 0;
    for (i = 0; i < 3; ++i) {
      if (vertices + v_i[i] == fv[j]) {
	// if free vertex, see next
        ++j;
      }
      else {
	// else zero gradient entry and hessian entries.
        g[i] = 0.;

        switch(i) {
        case 0:
          h[0].zero(); h[1].zero(); h[2].zero();
          break;
	  
        case 1:
          h[1].zero(); h[3].zero(); h[4].zero();
          break;
	  
        case 2:
          h[2].zero(); h[4].zero(); h[5].zero();
          break;
        }
      }
    }
#endif

    break;

  case QUADRILATERAL:

#ifndef ANALYTIC
    mValid = compute_element_numerical_hessian(pd, e, fv, g, h, nfv, m, err);
    return !MSQ_CHKERR(err) && mValid;
#else

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constaint, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    // Zero out the hessian and gradient vector
    for (i = 0; i < 4; ++i) {
      g[i] = 0.0;
    }

    for (i = 0; i < 10; ++i) {
      h[i].zero();
    }

    // Compute the metric and sum them together
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      pd.get_domain_normal_at_vertex(v_i[hexInd[i][0]], true, mNormal, err);
      MSQ_ERRZERO(err);

      inv(invW, W[i]);

      mValid = h_fcn_ridft2(mMetric, mGrads, mHessians, mCoords, mNormal,
			    invW, a, b, c, delta);
      if (!mValid) return false;

      m += W[i].get_cK() * mMetric;

      for (j = 0; j < 3; ++j) {
	g[hexInd[i][j]] += W[i].get_cK() * mGrads[j];
      }

      l = 0;
      for (j = 0; j < 3; ++j) {
	for (k = j; k < 3; ++k) {
	  row = hexInd[i][j];
	  col = hexInd[i][k];

	  if (row <= col) {
	    loc = 4*row - (row*(row+1)/2) + col;
	    h[loc] += W[i].get_cK() * mHessians[l];
	  }
	  else {
	    loc = 4*col - (col*(col+1)/2) + row;
	    h[loc] += W[i].get_cK() * transpose(mHessians[l]);
	  }
	  ++l;
	}
      }
    }

    m *= 0.25;
    for (i = 0; i < 4; ++i) {
      g[i] *= 0.25;
    }

    for (i = 0; i < 10; ++i) {
      h[i] *= 0.25;
    }

    // zero out fixed elements of gradient and Hessian
    j = 0;
    for (i = 0; i < 4; ++i) {
      if (vertices + v_i[i] == fv[j]) {
	// if free vertex, see next
        ++j;
      }
      else {
	// else zero gradient entry and hessian entries.
        g[i] = 0.;

        switch(i) {
        case 0:
          h[0].zero();   h[1].zero();   h[2].zero();   h[3].zero();
          break;
          
        case 1:
          h[1].zero();   h[4].zero();   h[5].zero();   h[6].zero();
          break;
          
        case 2:
          h[2].zero();   h[5].zero();   h[7].zero();  h[8].zero();
          break;
          
        case 3:
          h[3].zero();   h[6].zero();   h[8].zero();  h[9].zero();
          break;
        }
      }
    }
#endif

    break;

  case TETRAHEDRON:

    // Zero out the hessian and gradient vector
    for (i = 0; i < 4; ++i) {
      g[i] = 0.0;
    }

    for (i = 0; i < 10; ++i) {
      h[i].zero();
    }

    // Compute the metric and sum them together
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }

      inv(invW, W[i]);

      mValid = h_fcn_ridft3(mMetric, mGrads, mHessians, mCoords, 
			   invW, a, b, c, delta);
      if (!mValid) return false;

      m += W[i].get_cK() * mMetric;

      for (j = 0; j < 4; ++j) {
	g[tetInd[i][j]] += W[i].get_cK() * mGrads[j];
      }

      l = 0;
      for (j = 0; j < 4; ++j) {
	for (k = j; k < 4; ++k) {
	  row = tetInd[i][j];
	  col = tetInd[i][k];

	  if (row <= col) {
	    loc = 4*row - (row*(row+1)/2) + col;
	    h[loc] += W[i].get_cK() * mHessians[l];
	  }
	  else {
	    loc = 4*col - (col*(col+1)/2) + row;
	    h[loc] += W[i].get_cK() * transpose(mHessians[l]);
	  }
	  ++l;
	}
      }
    }

    m *= 0.25;
    for (i = 0; i < 4; ++i) {
      g[i] *= 0.25;
    }

    for (i = 0; i < 10; ++i) {
      h[i] *= 0.25;
    }

    // zero out fixed elements of g
    j = 0;
    for (i = 0; i < 4; ++i) {
      if (vertices + v_i[i] == fv[j]) {
	// if free vertex, see next
        ++j;
      }
      else {
	// else zero gradient entry and hessian entries.
        g[i] = 0.;

        switch(i) {
        case 0:
          h[0].zero(); h[1].zero(); h[2].zero(); h[3].zero();
          break;
	  
        case 1:
          h[1].zero(); h[4].zero(); h[5].zero(); h[6].zero();
          break;
	  
        case 2:
          h[2].zero(); h[5].zero(); h[7].zero(); h[8].zero();
          break;

        case 3:
          h[3].zero(); h[6].zero(); h[8].zero(); h[9].zero();
          break;
        }
      }
    }
    break;

  case HEXAHEDRON:

    // Zero out the hessian and gradient vector
    for (i = 0; i < 8; ++i) {
      g[i] = 0.0;
    }

    for (i = 0; i < 36; ++i) {
      h[i].zero();
    }

    // Compute the metric and sum them together
    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      inv(invW, W[i]);

      mValid = h_fcn_ridft3(mMetric, mGrads, mHessians, mCoords, 
			   invW, a, b, c, delta);
      if (!mValid) return false;

      m += W[i].get_cK() * mMetric;

      for (j = 0; j < 4; ++j) {
	g[hexInd[i][j]] += W[i].get_cK() * mGrads[j];
      }

      l = 0;
      for (j = 0; j < 4; ++j) {
	for (k = j; k < 4; ++k) {
	  row = hexInd[i][j];
	  col = hexInd[i][k];

	  if (row <= col) {
	    loc = 8*row - (row*(row+1)/2) + col;
	    h[loc] += W[i].get_cK() * mHessians[l];
	  }
	  else {
	    loc = 8*col - (col*(col+1)/2) + row;
	    h[loc] += W[i].get_cK() * transpose(mHessians[l]);
	  }
	  ++l;
	}
      }
    }

    m *= 0.125;
    for (i = 0; i < 8; ++i) {
      g[i] *= 0.125;
    }

    for (i = 0; i < 36; ++i) {
      h[i] *= 0.125;
    }

    // zero out fixed elements of gradient and Hessian
    j = 0;
    for (i = 0; i < 8; ++i) {
      if (vertices + v_i[i] == fv[j]) {
	// if free vertex, see next
        ++j;
      }
      else {
	// else zero gradient entry and hessian entries.
        g[i] = 0.;

        switch(i) {
        case 0:
          h[0].zero();   h[1].zero();   h[2].zero();   h[3].zero();
          h[4].zero();   h[5].zero();   h[6].zero();   h[7].zero();
          break;
          
        case 1:
          h[1].zero();   h[8].zero();   h[9].zero();   h[10].zero();
          h[11].zero();  h[12].zero();  h[13].zero();  h[14].zero();
          break;
          
        case 2:
          h[2].zero();   h[9].zero();   h[15].zero();  h[16].zero();
          h[17].zero();  h[18].zero();  h[19].zero();  h[20].zero();
          break;
          
        case 3:
          h[3].zero();   h[10].zero();  h[16].zero();  h[21].zero();
          h[22].zero();  h[23].zero();  h[24].zero();  h[25].zero();
          break;
          
        case 4:
          h[4].zero();   h[11].zero();  h[17].zero();  h[22].zero();
          h[26].zero();  h[27].zero();  h[28].zero();  h[29].zero();
          break;
          
        case 5:
          h[5].zero();   h[12].zero();  h[18].zero();  h[23].zero();
          h[27].zero();  h[30].zero();  h[31].zero();  h[32].zero();
          break;
          
        case 6:
          h[6].zero();   h[13].zero();  h[19].zero();  h[24].zero();
          h[28].zero();  h[31].zero();  h[33].zero();  h[34].zero();
          break;
          
        case 7:
          h[7].zero();   h[14].zero();  h[20].zero();  h[25].zero();
          h[29].zero();  h[32].zero();  h[34].zero();  h[35].zero();
          break;
        }
      }
    }
    break;

  default:
    MSQ_SETERR(err)("element type not implemented.",MsqError::NOT_IMPLEMENTED);
    return false;
  }

  return true;
}
