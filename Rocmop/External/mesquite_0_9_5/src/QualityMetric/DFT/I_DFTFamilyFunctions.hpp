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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file I_DFTFamilyFunctions.hpp

Header defining function, gradient, and Hessian for generalized
distance from target metrics.  This code replaces specialized versions
for the mu1 distance from target metric and the and inverse mean ratio
metric.  These metrics are recovered by setting the parameters to
specific values.  Note that the best efficiency will be obtained by
setting the parameter values as constants and eliminating dead code.
Further savings can be obtained by specializing the code for
particular weight matrices, such as the right tetrahedron and regular
tetrahedron.

The generalized distance from target metric is of the form:

                  alpha*|| A*inv(W) - beta*I ||_F^2
   ------------------------------------------------------------------
   0.5^gamma*(det(A*inv(W)) + sqrt(det(A*inv(W))^2 + 4*delta^2))^gamma

where gamma is a positive constant bounded above by one, and delta and
beta are nonnegative constants.  Delta is zero, the denominator becomes
max(det(A*inv(W)),0)^gamma, which is undefined for inverted elements.

Instead of passing inv(W) to the routine, Q and inv(R) from the QR 
factorization of W are passed, where Q and R have been formed so that 
transpose(Q)*Q = I, det(Q) = 1, and R and inv(R) are upper triangular 
matrices.  The metric is then equivalent to:

                  alpha*|| A*inv(R) - beta*Q ||_F^2
   ------------------------------------------------------------------
   0.5^gamma*(det(A*inv(R)) + sqrt(det(A*inv(R))^2 + 4*delta^2))^gamma

The function evaluation is cheaper to perform in this case because computing 
A*inv(R) requires 18 fewer operations than A*inv(W), while the modification 
for beta*Q adds 15 extra operations, resulting in a net reduction of 3 
operations.  Further savings are obtained in the analytic gradient and 
Hessian evaluations.  Moreover, the savings are magnified when beta has a 
value of zero or one.

Special cases of primary interest are:

   beta = 0, delta = 0:	inverse mean-ratio metric
   beta = 1, delta > 0:	distance from target metric (mu1)

The reason for the ordering of the constants is that we use default values 
in the function definitions (gamma= 1 or 2/3 depending on the element type, 
delta=0, beta=0) and rely upon the compiler to do specialization to gain 
better performance for these special cases.  Note that the defaults
recover the inverse mean-ratio metric.

\author Todd Munson
\date   2004-12-17
 */

#ifndef I_DFTFamilyFunctions_hpp
#define I_DFTFamilyFunctions_hpp

#include <math.h>
#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"

namespace Mesquite 
{
    //NOTE (mbrewer):  variables in the following discription may not be
    //consistent with the variables used above.
  /***************************************************************************/
  /* Gradient calculation courtesy of Paul Hovland.  The code was  modified  */
  /* to reduce the number of flops and intermediate variables, and improve   */
  /* the locality of reference.                                              */
  /***************************************************************************/
  /* The Hessian calculation is computes everyting in (x,y,z) blocks and     */
  /* stores only the upper triangular blocks.  The results are put into      */
  /* (c1,c2,c3,c4) order prior to returning the Hessian matrix.              */
  /***************************************************************************/
  /* The form of the function, gradient, and Hessian follows:                */
  /*   o(x) = a * pow(f(A(x)), b) * pow(g(A(x)), c)                          */
  /* where A(x) is the incidence matrix generated by:                        */
  /*           [x1-x0 x2-x0 x3-x0]                                           */
  /*    A(x) = [y1-y0 y2-y0 y3-y0] * inv(W)                                  */
  /*           [z1-z0 z2-z0 z3-z0]                                           */
  /* and f() is the squared Frobenius norm of A(x), and g() is the           */
  /* determinant of A(x).                                                    */
  /*                                                                         */
  /* The gradient is calculated as follows:                                  */
  /*   gamma := a*b*pow(f(A(x)),b-1)*pow(g(A(x)),c)                          */
  /*   delta  := a*c*pow(f(A(x)),b)*pow(g(A(x)),c-1)                         */
  /*                                                                         */
  /*   do/dx = (gamma * (df/dA) + delta * (dg/dA)) (dA/dx)                   */
  /*                                                                         */
  /*   (df/dA)_i = 2*A_i                                                     */
  /*   (dg/dA)_i = A_j*A_k - A_l*A_m for some {j,k,l,m}                      */
  /*                                                                         */
  /*   d^2o/dx^2 = (dA/dx)' * ((d gamma/dA) * (df/dA) +                      */
  /*                           (d  delta/dA) * (dg/dA)                       */
  /*                                  gamma * (d^2f/dA^2)                    */
  /*                                   delta * (d^2g/dA^2)) * (dA/dx)        */
  /*                                                                         */
  /*   Note: since A(x) is a linear function, there are no terms involving   */
  /*   d^2A/dx^2 since this matrix is zero.                                  */
  /*                                                                         */
  /*   gamma := a*b*c*pow(f(A(x)),b-1)*pow(g(A(x)),c-1)                      */
  /*   beta  := a*c*(c-1)*pow(f(A(x)),b)*pow(g(A(x)),c-2)                    */
  /*   psi   := a*b*(b-1)*pow(f(A(x)),b-2)*pow(g(A(x)),c)                    */
  /*                                                                         */
  /*   d^2o/dx^2 = (dA/dx)' * (gamma*((dg/dA)'*(df/dA) + (df/dA)'*(dg/dA)) + */
  /*                            beta* (dg/dA)'*(dg/dA) +                     */
  /*                             psi* (df/dA)'*(df/dA) +                     */
  /*                           gamma*(d^2f/dA^2) +                           */
  /*                           delta*(d^2g/dA^2)) * (dA/dx)                  */
  /*                                                                         */
  /*   Note: (df/dA) and (dg/dA) are row vectors and we only calculate the   */
  /*   upper triangular part of the inner matrices.                          */
  /***************************************************************************/

  /***************************************************************************/
  /* Triangular and quadrilateral versions of the metric.                    */
  /***************************************************************************/

  inline bool m_gdft_2(double &obj, 
                       const Vector3D x[3], const Vector3D &n,
		       const Matrix3D &invR,	/* upper triangular          */
		       const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
		       const double alpha = 0.5,/* constant                  */
		       const Exponent& gamma = 1.0,/* planar elements           */
		       const double delta = 0.0,/* max in denominator        */
		       const double beta  = 0.0)/* no modification           */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
         matr[3]*(matr[2]*matr[7] - matr[1]*matr[8]) +
         matr[6]*(matr[1]*matr[5] - matr[2]*matr[4]);

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    obj = alpha * pow(2.0, gamma) * f / pow(g, gamma);
    return true;
  }

  inline bool g_gdft_2(double &obj, Vector3D g_obj[3], 
                       const Vector3D x[3], const Vector3D &n,
		       const Matrix3D &invR,	/* upper triangular          */
		       const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
		       const double alpha = 0.5,/* constant                  */
		       const Exponent& gamma = 1.0,/* simplicial elements       */
		       const double delta = 0.0,/* max in denominator        */
		       const double beta  = 0.0)/* no modification           */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double adjm[9], loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Calculate adjoint matrix */
    adjm[0] = f*matd[0] + g*loc1;
    adjm[1] = f*matd[1] + g*loc2;
    adjm[2] = f*matd[2] + g*loc3;
    
    loc1 = g*matr[0];
    loc2 = g*matr[1];
    loc3 = g*matr[2];

    adjm[3] = f*matd[3] + loc3*matr[7] - loc2*matr[8];
    adjm[4] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
    adjm[5] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

    adjm[6] = f*matd[6] + loc2*matr[5] - loc3*matr[4];
    adjm[7] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
    adjm[8] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

    /* Construct gradients */
    g_obj[1][0] = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[2][0] =                    invR[1][1]*adjm[1]+invR[1][2]*adjm[2];
    g_obj[0][0] = -g_obj[1][0] - g_obj[2][0];

    g_obj[1][1] = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[2][1] =                    invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    g_obj[0][1] = -g_obj[1][1] - g_obj[2][1];

    g_obj[1][2] = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    g_obj[2][2] =                    invR[1][1]*adjm[7]+invR[1][2]*adjm[8];
    g_obj[0][2] = -g_obj[1][2] - g_obj[2][2];
    return true;
  }

  inline bool h_gdft_2(double &obj, Vector3D g_obj[3], Matrix3D h_obj[6],
		       const Vector3D x[3], const Vector3D &n,
		       const Matrix3D &invR,	/* upper triangular          */
		       const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
		       const double alpha = 0.5,/* constant                  */
		       const Exponent& gamma = 1.0,/* simplicial elements       */
		       const double delta = 0.0,/* max in denominator        */
		       const double beta  = 0.0)/* no modification           */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double adjm[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
    double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];
    double A[9];

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
    dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
    dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Calculate adjoint matrix */
    adjm[0] = dobj_df*matd[0] + dobj_dg*dg[0];
    adjm[1] = dobj_df*matd[1] + dobj_dg*dg[1];
    adjm[2] = dobj_df*matd[2] + dobj_dg*dg[2];
    adjm[3] = dobj_df*matd[3] + dobj_dg*dg[3];
    adjm[4] = dobj_df*matd[4] + dobj_dg*dg[4];
    adjm[5] = dobj_df*matd[5] + dobj_dg*dg[5];
    adjm[6] = dobj_df*matd[6] + dobj_dg*dg[6];
    adjm[7] = dobj_df*matd[7] + dobj_dg*dg[7];
    adjm[8] = dobj_df*matd[8] + dobj_dg*dg[8];

    /* Construct gradients */
    g_obj[1][0] = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[2][0] =                    invR[1][1]*adjm[1]+invR[1][2]*adjm[2];
    g_obj[0][0] = -g_obj[1][0] - g_obj[2][0];

    g_obj[1][1] = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[2][1] =                    invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    g_obj[0][1] = -g_obj[1][1] - g_obj[2][1];

    g_obj[1][2] = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    g_obj[2][2] =                    invR[1][1]*adjm[7]+invR[1][2]*adjm[8];
    g_obj[0][2] = -g_obj[1][2] - g_obj[2][2];

    /* Start of the Hessian evaluation */
    adjm[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
    adjm[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
    adjm[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
    adjm[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
    adjm[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
    adjm[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
    adjm[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
    adjm[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
    adjm[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[0] + matd[0];
    J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
    J_A[1] = dg[0]*matd[1] + loc1*dg[1];
    J_A[2] = dg[0]*matd[2] + loc1*dg[2];
    J_B[0] = dg[0]*matd[3] + loc1*dg[3];
    J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adjm[8];
    J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adjm[7];
    J_C[0] = dg[0]*matd[6] + loc1*dg[6];
    J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adjm[5];
    J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adjm[4];

    loc1 = dobj_dgdg*dg[1] + matd[1];
    J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
    J_A[4] = dg[1]*matd[2] + loc1*dg[2];
    J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adjm[8];
    J_B[4] = dg[1]*matd[4] + loc1*dg[4];
    J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adjm[6];
    J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adjm[5];
    J_C[4] = dg[1]*matd[7] + loc1*dg[7];
    J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adjm[3];

    loc1 = dobj_dgdg*dg[2] + matd[2];
    J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
    J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adjm[7];
    J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adjm[6];
    J_B[8] = dg[2]*matd[5] + loc1*dg[5];
    J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adjm[4];
    J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adjm[3];
    J_C[8] = dg[2]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[3] + matd[3];
    J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
    J_D[1] = dg[3]*matd[4] + loc1*dg[4];
    J_D[2] = dg[3]*matd[5] + loc1*dg[5];
    J_E[0] = dg[3]*matd[6] + loc1*dg[6];
    J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adjm[2];
    J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adjm[1];
  
    loc1 = dobj_dgdg*dg[4] + matd[4];
    J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
    J_D[4] = dg[4]*matd[5] + loc1*dg[5];
    J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adjm[2];
    J_E[4] = dg[4]*matd[7] + loc1*dg[7];
    J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adjm[0];
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
    J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adjm[1];
    J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adjm[0];
    J_E[8] = dg[5]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[6] + matd[6];
    J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
    J_F[1] = dg[6]*matd[7] + loc1*dg[7];
    J_F[2] = dg[6]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[7] + matd[7];
    J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
    J_F[4] = dg[7]*matd[8] + loc1*dg[8];
  
    J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

    /* Assemble matrix products */

    /* dx / dx */
    A[1]  =  J_A[0]*invR[0][0] + J_A[1]*invR[0][1] + J_A[2]*invR[0][2];
    A[2]  =                      J_A[1]*invR[1][1] + J_A[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_A[1]*invR[0][0] + J_A[3]*invR[0][1] + J_A[4]*invR[0][2];
    A[5]  =                      J_A[3]*invR[1][1] + J_A[4]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_A[2]*invR[0][0] + J_A[4]*invR[0][1] + J_A[5]*invR[0][2];
    A[8]  =                      J_A[4]*invR[1][1] + J_A[5]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[1][0][0] =  A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[2][0][0] =                    A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[0][0][0] = -h_obj[1][0][0] - h_obj[2][0][0];

    h_obj[3][0][0] =  A[1]*invR[0][0] + A[4]*invR[0][1] + A[7]*invR[0][2];
    h_obj[4][0][0] =                    A[4]*invR[1][1] + A[7]*invR[1][2];

    h_obj[5][0][0] =                    A[5]*invR[1][1] + A[8]*invR[1][2];

    /* dx / dy */
    A[1]  =  J_B[0]*invR[0][0] + J_B[1]*invR[0][1] + J_B[2]*invR[0][2];
    A[2]  =                      J_B[1]*invR[1][1] + J_B[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_B[3]*invR[0][0] + J_B[4]*invR[0][1] + J_B[5]*invR[0][2];
    A[5]  =                      J_B[4]*invR[1][1] + J_B[5]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_B[6]*invR[0][0] + J_B[7]*invR[0][1] + J_B[8]*invR[0][2];
    A[8]  =                      J_B[7]*invR[1][1] + J_B[8]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[1][1][0] = A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[3][0][1] = A[1]*invR[0][0] + A[4]*invR[0][1] + A[7]*invR[0][2];
    h_obj[4][0][1] = A[2]*invR[0][0] + A[5]*invR[0][1] + A[8]*invR[0][2];

    h_obj[2][1][0] = A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[4][1][0] = A[4]*invR[1][1] + A[7]*invR[1][2];
    h_obj[5][0][1] = A[5]*invR[1][1] + A[8]*invR[1][2];

    h_obj[0][0][1] = -h_obj[1][1][0] - h_obj[2][1][0];
    h_obj[1][0][1] = -h_obj[3][0][1] - h_obj[4][1][0];
    h_obj[2][0][1] = -h_obj[4][0][1] - h_obj[5][0][1];

    /* dx / dz */
    A[1]  =  J_C[0]*invR[0][0] + J_C[1]*invR[0][1] + J_C[2]*invR[0][2];
    A[2]  =                      J_C[1]*invR[1][1] + J_C[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_C[3]*invR[0][0] + J_C[4]*invR[0][1] + J_C[5]*invR[0][2];
    A[5]  =                      J_C[4]*invR[1][1] + J_C[5]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_C[6]*invR[0][0] + J_C[7]*invR[0][1] + J_C[8]*invR[0][2];
    A[8]  =                      J_C[7]*invR[1][1] + J_C[8]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[1][2][0] = A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[3][0][2] = A[1]*invR[0][0] + A[4]*invR[0][1] + A[7]*invR[0][2];
    h_obj[4][0][2] = A[2]*invR[0][0] + A[5]*invR[0][1] + A[8]*invR[0][2];

    h_obj[2][2][0] = A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[4][2][0] = A[4]*invR[1][1] + A[7]*invR[1][2];
    h_obj[5][0][2] = A[5]*invR[1][1] + A[8]*invR[1][2];

    h_obj[0][0][2] = -h_obj[1][2][0] - h_obj[2][2][0];
    h_obj[1][0][2] = -h_obj[3][0][2] - h_obj[4][2][0];
    h_obj[2][0][2] = -h_obj[4][0][2] - h_obj[5][0][2];

    /* dy / dy */
    A[1]  =  J_D[0]*invR[0][0] + J_D[1]*invR[0][1] + J_D[2]*invR[0][2];
    A[2]  =                      J_D[1]*invR[1][1] + J_D[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_D[1]*invR[0][0] + J_D[3]*invR[0][1] + J_D[4]*invR[0][2];
    A[5]  =                      J_D[3]*invR[1][1] + J_D[4]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_D[2]*invR[0][0] + J_D[4]*invR[0][1] + J_D[5]*invR[0][2];
    A[8]  =                      J_D[4]*invR[1][1] + J_D[5]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[1][1][1] =  A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[2][1][1] =                    A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[0][1][1] = -h_obj[1][1][1] - h_obj[2][1][1];

    h_obj[3][1][1] =  A[1]*invR[0][0] + A[4]*invR[0][1] + A[7]*invR[0][2];
    h_obj[4][1][1] =                    A[4]*invR[1][1] + A[7]*invR[1][2];

    h_obj[5][1][1] =                    A[5]*invR[1][1] + A[8]*invR[1][2];

    /* dy / dz */
    A[1]  =  J_E[0]*invR[0][0] + J_E[1]*invR[0][1] + J_E[2]*invR[0][2];
    A[2]  =                      J_E[1]*invR[1][1] + J_E[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_E[3]*invR[0][0] + J_E[4]*invR[0][1] + J_E[5]*invR[0][2];
    A[5]  =                      J_E[4]*invR[1][1] + J_E[5]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_E[6]*invR[0][0] + J_E[7]*invR[0][1] + J_E[8]*invR[0][2];
    A[8]  =                      J_E[7]*invR[1][1] + J_E[8]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[1][2][1] = A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[3][1][2] = A[1]*invR[0][0] + A[4]*invR[0][1] + A[7]*invR[0][2];
    h_obj[4][1][2] = A[2]*invR[0][0] + A[5]*invR[0][1] + A[8]*invR[0][2];

    h_obj[2][2][1] = A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[4][2][1] = A[4]*invR[1][1] + A[7]*invR[1][2];
    h_obj[5][1][2] = A[5]*invR[1][1] + A[8]*invR[1][2];

    h_obj[0][1][2] = -h_obj[1][2][1] - h_obj[2][2][1];
    h_obj[1][1][2] = -h_obj[3][1][2] - h_obj[4][2][1];
    h_obj[2][1][2] = -h_obj[4][1][2] - h_obj[5][1][2];

    /* dz / dz */
    A[1]  =  J_F[0]*invR[0][0] + J_F[1]*invR[0][1] + J_F[2]*invR[0][2];
    A[2]  =                      J_F[1]*invR[1][1] + J_F[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_F[1]*invR[0][0] + J_F[3]*invR[0][1] + J_F[4]*invR[0][2];
    A[5]  =                      J_F[3]*invR[1][1] + J_F[4]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_F[2]*invR[0][0] + J_F[4]*invR[0][1] + J_F[5]*invR[0][2];
    A[8]  =                      J_F[4]*invR[1][1] + J_F[5]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[1][2][2] =  A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[2][2][2] =                    A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[0][2][2] = -h_obj[1][2][2] - h_obj[2][2][2];

    h_obj[3][2][2] =  A[1]*invR[0][0] + A[4]*invR[0][1] + A[7]*invR[0][2];
    h_obj[4][2][2] =                    A[4]*invR[1][1] + A[7]*invR[1][2];

    h_obj[5][2][2] =                    A[5]*invR[1][1] + A[8]*invR[1][2];

    /* Complete diagonal blocks */
    h_obj[0].fill_lower_triangle();
    h_obj[3].fill_lower_triangle();
    h_obj[5].fill_lower_triangle();
    return true;
  }

  /***************************************************************************/
  /* Tetrahedral and hexahedral versions of the metric.                      */
  /***************************************************************************/

  inline bool m_gdft_3(double &obj, 
		       const Vector3D x[4],
		       const Matrix3D &invR,	/* upper triangular          */
		       const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
		       const double alpha = 1.0/3.0, /* constant             */
		       const Exponent& gamma = 2.0/3.0, /* simplicial elements  */
		       const double delta = 0.0,/* max in denominator        */
		       const double beta  = 0.0)/* no modification           */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    t1 = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
         matr[3]*(matr[2]*matr[7] - matr[1]*matr[8]) +
         matr[6]*(matr[1]*matr[5] - matr[2]*matr[4]);

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    obj = alpha * pow(2.0, gamma) * f / pow(g, gamma);
    return true;
  }

  inline bool g_gdft_3(double &obj, Vector3D g_obj[4], 
                       const Vector3D x[4],
		       const Matrix3D &invR,	/* upper triangular          */
		       const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
		       const double alpha = 1.0/3.0, /* constant             */
		       const Exponent& gamma = 2.0/3.0, /* simplicial elements  */
		       const double delta = 0.0,/* max in denominator        */
		       const double beta  = 0.0)/* no modification           */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double adjm[9], loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Calculate adjoint matrix */
    adjm[0] = f*matd[0] + g*loc1;
    adjm[1] = f*matd[1] + g*loc2;
    adjm[2] = f*matd[2] + g*loc3;

    loc1 = g*matr[0];
    loc2 = g*matr[1];
    loc3 = g*matr[2];

    adjm[3] = f*matd[3] + loc3*matr[7] - loc2*matr[8];
    adjm[4] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
    adjm[5] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

    adjm[6] = f*matd[6] + loc2*matr[5] - loc3*matr[4];
    adjm[7] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
    adjm[8] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

    /* Construct gradients */
    g_obj[1][0] = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[2][0] =                    invR[1][1]*adjm[1]+invR[1][2]*adjm[2];
    g_obj[3][0] =                                       invR[2][2]*adjm[2];
    g_obj[0][0] = -g_obj[1][0] - g_obj[2][0] - g_obj[3][0];

    g_obj[1][1] = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[2][1] =                    invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    g_obj[3][1] =                                       invR[2][2]*adjm[5];
    g_obj[0][1] = -g_obj[1][1] - g_obj[2][1] - g_obj[3][1];

    g_obj[1][2] = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    g_obj[2][2] =                    invR[1][1]*adjm[7]+invR[1][2]*adjm[8];
    g_obj[3][2] =                                       invR[2][2]*adjm[8];
    g_obj[0][2] = -g_obj[1][2] - g_obj[2][2] - g_obj[3][2];
    return true;
  }

  inline bool h_gdft_3(double &obj, Vector3D g_obj[4], Matrix3D h_obj[10],
		       const Vector3D x[4],
		       const Matrix3D &invR,	/* upper triangular          */
		       const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
		       const double alpha = 1.0/3.0, /* constant             */
		       const Exponent& gamma = 2.0/3.0, /* simplicial elements  */
		       const double delta = 0.0,/* max in denominator        */
		       const double beta  = 0.0)/* no modification           */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double adjm[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
    double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];
    double A[12];

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
    dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
    dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Calculate adjoint matrix */
    adjm[0] = dobj_df*matd[0] + dobj_dg*dg[0];
    adjm[1] = dobj_df*matd[1] + dobj_dg*dg[1];
    adjm[2] = dobj_df*matd[2] + dobj_dg*dg[2];
    adjm[3] = dobj_df*matd[3] + dobj_dg*dg[3];
    adjm[4] = dobj_df*matd[4] + dobj_dg*dg[4];
    adjm[5] = dobj_df*matd[5] + dobj_dg*dg[5];
    adjm[6] = dobj_df*matd[6] + dobj_dg*dg[6];
    adjm[7] = dobj_df*matd[7] + dobj_dg*dg[7];
    adjm[8] = dobj_df*matd[8] + dobj_dg*dg[8];

    /* Construct gradients */
    g_obj[1][0] = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[2][0] =                    invR[1][1]*adjm[1]+invR[1][2]*adjm[2];
    g_obj[3][0] =                                       invR[2][2]*adjm[2];
    g_obj[0][0] = -g_obj[1][0] - g_obj[2][0] - g_obj[3][0];

    g_obj[1][1] = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[2][1] =                    invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    g_obj[3][1] =                                       invR[2][2]*adjm[5];
    g_obj[0][1] = -g_obj[1][1] - g_obj[2][1] - g_obj[3][1];

    g_obj[1][2] = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    g_obj[2][2] =                    invR[1][1]*adjm[7]+invR[1][2]*adjm[8];
    g_obj[3][2] =                                       invR[2][2]*adjm[8];
    g_obj[0][2] = -g_obj[1][2] - g_obj[2][2] - g_obj[3][2];

    /* Start of the Hessian evaluation */
    adjm[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
    adjm[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
    adjm[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
    adjm[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
    adjm[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
    adjm[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
    adjm[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
    adjm[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
    adjm[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[0] + matd[0];
    J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
    J_A[1] = dg[0]*matd[1] + loc1*dg[1];
    J_A[2] = dg[0]*matd[2] + loc1*dg[2];
    J_B[0] = dg[0]*matd[3] + loc1*dg[3];
    J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adjm[8];
    J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adjm[7];
    J_C[0] = dg[0]*matd[6] + loc1*dg[6];
    J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adjm[5];
    J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adjm[4];

    loc1 = dobj_dgdg*dg[1] + matd[1];
    J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
    J_A[4] = dg[1]*matd[2] + loc1*dg[2];
    J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adjm[8];
    J_B[4] = dg[1]*matd[4] + loc1*dg[4];
    J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adjm[6];
    J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adjm[5];
    J_C[4] = dg[1]*matd[7] + loc1*dg[7];
    J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adjm[3];

    loc1 = dobj_dgdg*dg[2] + matd[2];
    J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
    J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adjm[7];
    J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adjm[6];
    J_B[8] = dg[2]*matd[5] + loc1*dg[5];
    J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adjm[4];
    J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adjm[3];
    J_C[8] = dg[2]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[3] + matd[3];
    J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
    J_D[1] = dg[3]*matd[4] + loc1*dg[4];
    J_D[2] = dg[3]*matd[5] + loc1*dg[5];
    J_E[0] = dg[3]*matd[6] + loc1*dg[6];
    J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adjm[2];
    J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adjm[1];
  
    loc1 = dobj_dgdg*dg[4] + matd[4];
    J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
    J_D[4] = dg[4]*matd[5] + loc1*dg[5];
    J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adjm[2];
    J_E[4] = dg[4]*matd[7] + loc1*dg[7];
    J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adjm[0];
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
    J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adjm[1];
    J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adjm[0];
    J_E[8] = dg[5]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[6] + matd[6];
    J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
    J_F[1] = dg[6]*matd[7] + loc1*dg[7];
    J_F[2] = dg[6]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[7] + matd[7];
    J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
    J_F[4] = dg[7]*matd[8] + loc1*dg[8];
  
    J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

    /* Assemble matrix products */

    /* dx / dx */
    A[1]  =  J_A[0]*invR[0][0] + J_A[1]*invR[0][1] + J_A[2]*invR[0][2];
    A[2]  =                      J_A[1]*invR[1][1] + J_A[2]*invR[1][2];
    A[3]  =                                          J_A[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_A[1]*invR[0][0] + J_A[3]*invR[0][1] + J_A[4]*invR[0][2];
    A[6]  =                      J_A[3]*invR[1][1] + J_A[4]*invR[1][2];
    A[7]  =                                          J_A[4]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_A[2]*invR[0][0] + J_A[4]*invR[0][1] + J_A[5]*invR[0][2];
    A[10] =                      J_A[4]*invR[1][1] + J_A[5]*invR[1][2];
    A[11] =                                          J_A[5]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[1][0][0] =  A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[2][0][0] =                    A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[3][0][0] =                                      A[8]*invR[2][2];
    h_obj[0][0][0] = -h_obj[1][0][0] - h_obj[2][0][0] - h_obj[3][0][0];

    h_obj[4][0][0] =  A[1]*invR[0][0] + A[5]*invR[0][1] + A[9]*invR[0][2];
    h_obj[5][0][0] =                    A[5]*invR[1][1] + A[9]*invR[1][2];
    h_obj[6][0][0] =                                      A[9]*invR[2][2];

    h_obj[7][0][0] =                    A[6]*invR[1][1] + A[10]*invR[1][2];
    h_obj[8][0][0] =                                      A[10]*invR[2][2];

    h_obj[9][0][0] =                                      A[11]*invR[2][2];

    /* dx / dy */
    A[1]  =  J_B[0]*invR[0][0] + J_B[1]*invR[0][1] + J_B[2]*invR[0][2];
    A[2]  =                      J_B[1]*invR[1][1] + J_B[2]*invR[1][2];
    A[3]  =                                          J_B[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_B[3]*invR[0][0] + J_B[4]*invR[0][1] + J_B[5]*invR[0][2];
    A[6]  =                      J_B[4]*invR[1][1] + J_B[5]*invR[1][2];
    A[7]  =                                          J_B[5]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_B[6]*invR[0][0] + J_B[7]*invR[0][1] + J_B[8]*invR[0][2];
    A[10] =                      J_B[7]*invR[1][1] + J_B[8]*invR[1][2];
    A[11] =                                          J_B[8]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[1][1][0] = A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[4][0][1] = A[1]*invR[0][0] + A[5]*invR[0][1] + A[9]*invR[0][2];
    h_obj[5][0][1] = A[2]*invR[0][0] + A[6]*invR[0][1] + A[10]*invR[0][2];
    h_obj[6][0][1] = A[3]*invR[0][0] + A[7]*invR[0][1] + A[11]*invR[0][2];

    h_obj[2][1][0] = A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[5][1][0] = A[5]*invR[1][1] + A[9]*invR[1][2];
    h_obj[7][0][1] = A[6]*invR[1][1] + A[10]*invR[1][2];
    h_obj[8][0][1] = A[7]*invR[1][1] + A[11]*invR[1][2];

    h_obj[3][1][0] = A[8]*invR[2][2];
    h_obj[6][1][0] = A[9]*invR[2][2];
    h_obj[8][1][0] = A[10]*invR[2][2];
    h_obj[9][0][1] = A[11]*invR[2][2];

    h_obj[0][0][1] = -h_obj[1][1][0] - h_obj[2][1][0] - h_obj[3][1][0];
    h_obj[1][0][1] = -h_obj[4][0][1] - h_obj[5][1][0] - h_obj[6][1][0];
    h_obj[2][0][1] = -h_obj[5][0][1] - h_obj[7][0][1] - h_obj[8][1][0];
    h_obj[3][0][1] = -h_obj[6][0][1] - h_obj[8][0][1] - h_obj[9][0][1];

    /* dx / dz */
    A[1]  =  J_C[0]*invR[0][0] + J_C[1]*invR[0][1] + J_C[2]*invR[0][2];
    A[2]  =                      J_C[1]*invR[1][1] + J_C[2]*invR[1][2];
    A[3]  =                                          J_C[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_C[3]*invR[0][0] + J_C[4]*invR[0][1] + J_C[5]*invR[0][2];
    A[6]  =                      J_C[4]*invR[1][1] + J_C[5]*invR[1][2];
    A[7]  =                                          J_C[5]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_C[6]*invR[0][0] + J_C[7]*invR[0][1] + J_C[8]*invR[0][2];
    A[10] =                      J_C[7]*invR[1][1] + J_C[8]*invR[1][2];
    A[11] =                                          J_C[8]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[1][2][0] = A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[4][0][2] = A[1]*invR[0][0] + A[5]*invR[0][1] + A[9]*invR[0][2];
    h_obj[5][0][2] = A[2]*invR[0][0] + A[6]*invR[0][1] + A[10]*invR[0][2];
    h_obj[6][0][2] = A[3]*invR[0][0] + A[7]*invR[0][1] + A[11]*invR[0][2];

    h_obj[2][2][0] = A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[5][2][0] = A[5]*invR[1][1] + A[9]*invR[1][2];
    h_obj[7][0][2] = A[6]*invR[1][1] + A[10]*invR[1][2];
    h_obj[8][0][2] = A[7]*invR[1][1] + A[11]*invR[1][2];

    h_obj[3][2][0] = A[8]*invR[2][2];
    h_obj[6][2][0] = A[9]*invR[2][2];
    h_obj[8][2][0] = A[10]*invR[2][2];
    h_obj[9][0][2] = A[11]*invR[2][2];

    h_obj[0][0][2] = -h_obj[1][2][0] - h_obj[2][2][0] - h_obj[3][2][0];
    h_obj[1][0][2] = -h_obj[4][0][2] - h_obj[5][2][0] - h_obj[6][2][0];
    h_obj[2][0][2] = -h_obj[5][0][2] - h_obj[7][0][2] - h_obj[8][2][0];
    h_obj[3][0][2] = -h_obj[6][0][2] - h_obj[8][0][2] - h_obj[9][0][2];

    /* dy / dy */
    A[1]  =  J_D[0]*invR[0][0] + J_D[1]*invR[0][1] + J_D[2]*invR[0][2];
    A[2]  =                      J_D[1]*invR[1][1] + J_D[2]*invR[1][2];
    A[3]  =                                          J_D[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_D[1]*invR[0][0] + J_D[3]*invR[0][1] + J_D[4]*invR[0][2];
    A[6]  =                      J_D[3]*invR[1][1] + J_D[4]*invR[1][2];
    A[7]  =                                          J_D[4]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_D[2]*invR[0][0] + J_D[4]*invR[0][1] + J_D[5]*invR[0][2];
    A[10] =                      J_D[4]*invR[1][1] + J_D[5]*invR[1][2];
    A[11] =                                          J_D[5]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[1][1][1] =  A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[2][1][1] =                    A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[3][1][1] =                                      A[8]*invR[2][2];
    h_obj[0][1][1] = -h_obj[1][1][1] - h_obj[2][1][1] - h_obj[3][1][1];

    h_obj[4][1][1] =  A[1]*invR[0][0] + A[5]*invR[0][1] + A[9]*invR[0][2];
    h_obj[5][1][1] =                    A[5]*invR[1][1] + A[9]*invR[1][2];
    h_obj[6][1][1] =                                      A[9]*invR[2][2];

    h_obj[7][1][1] =                    A[6]*invR[1][1] + A[10]*invR[1][2];
    h_obj[8][1][1] =                                      A[10]*invR[2][2];

    h_obj[9][1][1] =                                      A[11]*invR[2][2];

    /* dy / dz */
    A[1]  =  J_E[0]*invR[0][0] + J_E[1]*invR[0][1] + J_E[2]*invR[0][2];
    A[2]  =                      J_E[1]*invR[1][1] + J_E[2]*invR[1][2];
    A[3]  =                                          J_E[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_E[3]*invR[0][0] + J_E[4]*invR[0][1] + J_E[5]*invR[0][2];
    A[6]  =                      J_E[4]*invR[1][1] + J_E[5]*invR[1][2];
    A[7]  =                                          J_E[5]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_E[6]*invR[0][0] + J_E[7]*invR[0][1] + J_E[8]*invR[0][2];
    A[10] =                      J_E[7]*invR[1][1] + J_E[8]*invR[1][2];
    A[11] =                                          J_E[8]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[1][2][1] = A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[4][1][2] = A[1]*invR[0][0] + A[5]*invR[0][1] + A[9]*invR[0][2];
    h_obj[5][1][2] = A[2]*invR[0][0] + A[6]*invR[0][1] + A[10]*invR[0][2];
    h_obj[6][1][2] = A[3]*invR[0][0] + A[7]*invR[0][1] + A[11]*invR[0][2];

    h_obj[2][2][1] = A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[5][2][1] = A[5]*invR[1][1] + A[9]*invR[1][2];
    h_obj[7][1][2] = A[6]*invR[1][1] + A[10]*invR[1][2];
    h_obj[8][1][2] = A[7]*invR[1][1] + A[11]*invR[1][2];

    h_obj[3][2][1] = A[8]*invR[2][2];
    h_obj[6][2][1] = A[9]*invR[2][2];
    h_obj[8][2][1] = A[10]*invR[2][2];
    h_obj[9][1][2] = A[11]*invR[2][2];

    h_obj[0][1][2] = -h_obj[1][2][1] - h_obj[2][2][1] - h_obj[3][2][1];
    h_obj[1][1][2] = -h_obj[4][1][2] - h_obj[5][2][1] - h_obj[6][2][1];
    h_obj[2][1][2] = -h_obj[5][1][2] - h_obj[7][1][2] - h_obj[8][2][1];
    h_obj[3][1][2] = -h_obj[6][1][2] - h_obj[8][1][2] - h_obj[9][1][2];

    /* dz / dz */
    A[1]  =  J_F[0]*invR[0][0] + J_F[1]*invR[0][1] + J_F[2]*invR[0][2];
    A[2]  =                      J_F[1]*invR[1][1] + J_F[2]*invR[1][2];
    A[3]  =                                          J_F[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_F[1]*invR[0][0] + J_F[3]*invR[0][1] + J_F[4]*invR[0][2];
    A[6]  =                      J_F[3]*invR[1][1] + J_F[4]*invR[1][2];
    A[7]  =                                          J_F[4]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_F[2]*invR[0][0] + J_F[4]*invR[0][1] + J_F[5]*invR[0][2];
    A[10] =                      J_F[4]*invR[1][1] + J_F[5]*invR[1][2];
    A[11] =                                          J_F[5]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[1][2][2] =  A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[2][2][2] =                    A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[3][2][2] =                                      A[8]*invR[2][2];
    h_obj[0][2][2] = -h_obj[1][2][2] - h_obj[2][2][2] - h_obj[3][2][2];

    h_obj[4][2][2] =  A[1]*invR[0][0] + A[5]*invR[0][1] + A[9]*invR[0][2];
    h_obj[5][2][2] =                    A[5]*invR[1][1] + A[9]*invR[1][2];
    h_obj[6][2][2] =                                      A[9]*invR[2][2];

    h_obj[7][2][2] =                    A[6]*invR[1][1] + A[10]*invR[1][2];
    h_obj[8][2][2] =                                      A[10]*invR[2][2];

    h_obj[9][2][2] =                                      A[11]*invR[2][2];

    /* Complete diagonal blocks */
    h_obj[0].fill_lower_triangle();
    h_obj[4].fill_lower_triangle();
    h_obj[7].fill_lower_triangle();
    h_obj[9].fill_lower_triangle();
    return true;
  }

  /***************************************************************************/
  /* The following code computes the gradient and Hessian with respect to a  */
  /* particular vertex.  These specializations are required to get good      */
  /* performance out of a local code.  The fastest versions would need to    */
  /* do store a seperate version of invR and Q for each vertex, which is not */
  /* an efficient use of the storage.                                        */
  /***************************************************************************/

  inline bool g_gdft_2_v0(double &obj, Vector3D &g_obj, 
			  const Vector3D x[3], const Vector3D &n,
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 0.5,/* constant               */
			  const Exponent& gamma = 1.0,/* simplicial elements    */
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double adjm[9], loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Calculate adjoint matrix */
    adjm[0] = f*matd[0] + g*loc1;
    adjm[1] = f*matd[1] + g*loc2;
    adjm[2] = f*matd[2] + g*loc3;
    
    loc1 = g*matr[0];
    loc2 = g*matr[1];
    loc3 = g*matr[2];

    adjm[3] = f*matd[3] + loc3*matr[7] - loc2*matr[8];
    adjm[4] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
    adjm[5] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

    adjm[6] = f*matd[6] + loc2*matr[5] - loc3*matr[4];
    adjm[7] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
    adjm[8] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

    /* Construct gradients */
    g_obj[0]  = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[0] +=                    invR[1][1]*adjm[1]+invR[1][2]*adjm[2];
    g_obj[0]  = -g_obj[0];

    g_obj[1]  = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[1] +=                    invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    g_obj[1]  = -g_obj[1];

    g_obj[2]  = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    g_obj[2] +=                    invR[1][1]*adjm[7]+invR[1][2]*adjm[8];
    g_obj[2]  = -g_obj[2];
    return true;
  }

  inline bool g_gdft_2_v1(double &obj, Vector3D &g_obj, 
			  const Vector3D x[3], const Vector3D &n,
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 0.5,/* constant               */
			  const Exponent& gamma = 1.0,/* simplicial elements    */
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double adjm[9], loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Calculate adjoint matrix */
    adjm[0] = f*matd[0] + g*loc1;
    adjm[1] = f*matd[1] + g*loc2;
    adjm[2] = f*matd[2] + g*loc3;
    
    loc1 = g*matr[0];
    loc2 = g*matr[1];
    loc3 = g*matr[2];

    adjm[3] = f*matd[3] + loc3*matr[7] - loc2*matr[8];
    adjm[4] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
    adjm[5] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

    adjm[6] = f*matd[6] + loc2*matr[5] - loc3*matr[4];
    adjm[7] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
    adjm[8] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

    /* Construct gradients */
    g_obj[0] = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[1] = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[2] = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    return true;
  }

  inline bool g_gdft_2_v2(double &obj, Vector3D &g_obj, 
			  const Vector3D x[3], const Vector3D &n,
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 0.5,/* constant               */
			  const Exponent& gamma = 1.0,/* simplicial elements    */
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double adjm[6], loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Calculate adjoint matrix */
    adjm[0] = f*matd[1] + g*loc2;
    adjm[1] = f*matd[2] + g*loc3;
    
    loc1 = g*matr[0];
    loc2 = g*matr[1];
    loc3 = g*matr[2];

    adjm[2] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
    adjm[3] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

    adjm[4] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
    adjm[5] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

    /* Construct gradients */
    g_obj[0] = invR[1][1]*adjm[0]+invR[1][2]*adjm[1];
    g_obj[1] = invR[1][1]*adjm[2]+invR[1][2]*adjm[3];
    g_obj[2] = invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    return true;
  }

  inline bool h_gdft_2_v0(double &obj, Vector3D &g_obj, Matrix3D &h_obj,
			  const Vector3D x[3], const Vector3D &n,
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 0.5,/* constant               */
			  const Exponent& gamma = 1.0,/* simplicial elements    */
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double adjm[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
    double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];
    double A[9];

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
    dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
    dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Calculate adjoint matrix */
    adjm[0] = dobj_df*matd[0] + dobj_dg*dg[0];
    adjm[1] = dobj_df*matd[1] + dobj_dg*dg[1];
    adjm[2] = dobj_df*matd[2] + dobj_dg*dg[2];
    adjm[3] = dobj_df*matd[3] + dobj_dg*dg[3];
    adjm[4] = dobj_df*matd[4] + dobj_dg*dg[4];
    adjm[5] = dobj_df*matd[5] + dobj_dg*dg[5];
    adjm[6] = dobj_df*matd[6] + dobj_dg*dg[6];
    adjm[7] = dobj_df*matd[7] + dobj_dg*dg[7];
    adjm[8] = dobj_df*matd[8] + dobj_dg*dg[8];

    /* Construct gradients */
    g_obj[0]  = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[0] +=                    invR[1][1]*adjm[1]+invR[1][2]*adjm[2];
    g_obj[0]  = -g_obj[0];

    g_obj[1]  = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[1] +=                    invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    g_obj[1]  = -g_obj[1];

    g_obj[2]  = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    g_obj[2] +=                    invR[1][1]*adjm[7]+invR[1][2]*adjm[8];
    g_obj[2]  = -g_obj[2];

    /* Start of the Hessian evaluation */
    adjm[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
    adjm[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
    adjm[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
    adjm[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
    adjm[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
    adjm[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
    adjm[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
    adjm[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
    adjm[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[0] + matd[0];
    J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
    J_A[1] = dg[0]*matd[1] + loc1*dg[1];
    J_A[2] = dg[0]*matd[2] + loc1*dg[2];
    J_B[0] = dg[0]*matd[3] + loc1*dg[3];
    J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adjm[8];
    J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adjm[7];
    J_C[0] = dg[0]*matd[6] + loc1*dg[6];
    J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adjm[5];
    J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adjm[4];

    loc1 = dobj_dgdg*dg[1] + matd[1];
    J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
    J_A[4] = dg[1]*matd[2] + loc1*dg[2];
    J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adjm[8];
    J_B[4] = dg[1]*matd[4] + loc1*dg[4];
    J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adjm[6];
    J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adjm[5];
    J_C[4] = dg[1]*matd[7] + loc1*dg[7];
    J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adjm[3];

    loc1 = dobj_dgdg*dg[2] + matd[2];
    J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
    J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adjm[7];
    J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adjm[6];
    J_B[8] = dg[2]*matd[5] + loc1*dg[5];
    J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adjm[4];
    J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adjm[3];
    J_C[8] = dg[2]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[3] + matd[3];
    J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
    J_D[1] = dg[3]*matd[4] + loc1*dg[4];
    J_D[2] = dg[3]*matd[5] + loc1*dg[5];
    J_E[0] = dg[3]*matd[6] + loc1*dg[6];
    J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adjm[2];
    J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adjm[1];
  
    loc1 = dobj_dgdg*dg[4] + matd[4];
    J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
    J_D[4] = dg[4]*matd[5] + loc1*dg[5];
    J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adjm[2];
    J_E[4] = dg[4]*matd[7] + loc1*dg[7];
    J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adjm[0];
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
    J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adjm[1];
    J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adjm[0];
    J_E[8] = dg[5]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[6] + matd[6];
    J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
    J_F[1] = dg[6]*matd[7] + loc1*dg[7];
    J_F[2] = dg[6]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[7] + matd[7];
    J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
    J_F[4] = dg[7]*matd[8] + loc1*dg[8];
  
    J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

    /* Assemble matrix products */

    /* dx / dx */
    A[1]  =  J_A[0]*invR[0][0] + J_A[1]*invR[0][1] + J_A[2]*invR[0][2];
    A[2]  =                      J_A[1]*invR[1][1] + J_A[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_A[1]*invR[0][0] + J_A[3]*invR[0][1] + J_A[4]*invR[0][2];
    A[5]  =                      J_A[3]*invR[1][1] + J_A[4]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_A[2]*invR[0][0] + J_A[4]*invR[0][1] + J_A[5]*invR[0][2];
    A[8]  =                      J_A[4]*invR[1][1] + J_A[5]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[0][0]  =  A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[0][0] +=                    A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[0][0]  = -h_obj[0][0];

    /* dx / dy */
    A[1]  =  J_B[0]*invR[0][0] + J_B[1]*invR[0][1] + J_B[2]*invR[0][2];
    A[2]  =                      J_B[1]*invR[1][1] + J_B[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_B[3]*invR[0][0] + J_B[4]*invR[0][1] + J_B[5]*invR[0][2];
    A[5]  =                      J_B[4]*invR[1][1] + J_B[5]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_B[6]*invR[0][0] + J_B[7]*invR[0][1] + J_B[8]*invR[0][2];
    A[8]  =                      J_B[7]*invR[1][1] + J_B[8]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[0][1]  = A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[0][1] +=                   A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[0][1]  = -h_obj[0][1];

    /* dx / dz */
    A[1]  =  J_C[0]*invR[0][0] + J_C[1]*invR[0][1] + J_C[2]*invR[0][2];
    A[2]  =                      J_C[1]*invR[1][1] + J_C[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_C[3]*invR[0][0] + J_C[4]*invR[0][1] + J_C[5]*invR[0][2];
    A[5]  =                      J_C[4]*invR[1][1] + J_C[5]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_C[6]*invR[0][0] + J_C[7]*invR[0][1] + J_C[8]*invR[0][2];
    A[8]  =                      J_C[7]*invR[1][1] + J_C[8]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[0][2]  = A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[0][2] +=                   A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[0][2]  = -h_obj[0][2];

    /* dy / dy */
    A[1]  =  J_D[0]*invR[0][0] + J_D[1]*invR[0][1] + J_D[2]*invR[0][2];
    A[2]  =                      J_D[1]*invR[1][1] + J_D[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_D[1]*invR[0][0] + J_D[3]*invR[0][1] + J_D[4]*invR[0][2];
    A[5]  =                      J_D[3]*invR[1][1] + J_D[4]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_D[2]*invR[0][0] + J_D[4]*invR[0][1] + J_D[5]*invR[0][2];
    A[8]  =                      J_D[4]*invR[1][1] + J_D[5]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[1][1]  =  A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[1][1] +=                    A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[1][1]  = -h_obj[1][1];

    /* dy / dz */
    A[1]  =  J_E[0]*invR[0][0] + J_E[1]*invR[0][1] + J_E[2]*invR[0][2];
    A[2]  =                      J_E[1]*invR[1][1] + J_E[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_E[3]*invR[0][0] + J_E[4]*invR[0][1] + J_E[5]*invR[0][2];
    A[5]  =                      J_E[4]*invR[1][1] + J_E[5]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_E[6]*invR[0][0] + J_E[7]*invR[0][1] + J_E[8]*invR[0][2];
    A[8]  =                      J_E[7]*invR[1][1] + J_E[8]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[1][2]  = A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[1][2] +=                   A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[1][2] = -h_obj[1][2];

    /* dz / dz */
    A[1]  =  J_F[0]*invR[0][0] + J_F[1]*invR[0][1] + J_F[2]*invR[0][2];
    A[2]  =                      J_F[1]*invR[1][1] + J_F[2]*invR[1][2];
    A[0]  = -A[1] - A[2];

    A[4]  =  J_F[1]*invR[0][0] + J_F[3]*invR[0][1] + J_F[4]*invR[0][2];
    A[5]  =                      J_F[3]*invR[1][1] + J_F[4]*invR[1][2];
    A[3]  = -A[4] - A[5];

    A[7]  =  J_F[2]*invR[0][0] + J_F[4]*invR[0][1] + J_F[5]*invR[0][2];
    A[8]  =                      J_F[4]*invR[1][1] + J_F[5]*invR[1][2];
    A[6]  = -A[7] - A[8];

    h_obj[2][2]  =  A[0]*invR[0][0] + A[3]*invR[0][1] + A[6]*invR[0][2];
    h_obj[2][2] +=                    A[3]*invR[1][1] + A[6]*invR[1][2];
    h_obj[2][2]  = -h_obj[2][2];

    /* Complete diagonal blocks */
    h_obj.fill_lower_triangle();
    return true;
  }

  inline bool h_gdft_2_v1(double &obj, Vector3D &g_obj, Matrix3D &h_obj,
			  const Vector3D x[3], const Vector3D &n,
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 0.5,/* constant               */
			  const Exponent& gamma = 1.0,/* simplicial elements    */
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double adjm[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
    double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];
    double A[3];

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
    dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
    dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Calculate adjoint matrix */
    adjm[0] = dobj_df*matd[0] + dobj_dg*dg[0];
    adjm[1] = dobj_df*matd[1] + dobj_dg*dg[1];
    adjm[2] = dobj_df*matd[2] + dobj_dg*dg[2];
    adjm[3] = dobj_df*matd[3] + dobj_dg*dg[3];
    adjm[4] = dobj_df*matd[4] + dobj_dg*dg[4];
    adjm[5] = dobj_df*matd[5] + dobj_dg*dg[5];
    adjm[6] = dobj_df*matd[6] + dobj_dg*dg[6];
    adjm[7] = dobj_df*matd[7] + dobj_dg*dg[7];
    adjm[8] = dobj_df*matd[8] + dobj_dg*dg[8];

    /* Construct gradients */
    g_obj[0] = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[1] = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[2] = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];

    /* Start of the Hessian evaluation */
    adjm[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
    adjm[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
    adjm[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
    adjm[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
    adjm[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
    adjm[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
    adjm[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
    adjm[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
    adjm[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[0] + matd[0];
    J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
    J_A[1] = dg[0]*matd[1] + loc1*dg[1];
    J_A[2] = dg[0]*matd[2] + loc1*dg[2];
    J_B[0] = dg[0]*matd[3] + loc1*dg[3];
    J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adjm[8];
    J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adjm[7];
    J_C[0] = dg[0]*matd[6] + loc1*dg[6];
    J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adjm[5];
    J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adjm[4];

    loc1 = dobj_dgdg*dg[1] + matd[1];
    J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
    J_A[4] = dg[1]*matd[2] + loc1*dg[2];
    J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adjm[8];
    J_B[4] = dg[1]*matd[4] + loc1*dg[4];
    J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adjm[6];
    J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adjm[5];
    J_C[4] = dg[1]*matd[7] + loc1*dg[7];
    J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adjm[3];

    loc1 = dobj_dgdg*dg[2] + matd[2];
    J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
    J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adjm[7];
    J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adjm[6];
    J_B[8] = dg[2]*matd[5] + loc1*dg[5];
    J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adjm[4];
    J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adjm[3];
    J_C[8] = dg[2]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[3] + matd[3];
    J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
    J_D[1] = dg[3]*matd[4] + loc1*dg[4];
    J_D[2] = dg[3]*matd[5] + loc1*dg[5];
    J_E[0] = dg[3]*matd[6] + loc1*dg[6];
    J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adjm[2];
    J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adjm[1];
  
    loc1 = dobj_dgdg*dg[4] + matd[4];
    J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
    J_D[4] = dg[4]*matd[5] + loc1*dg[5];
    J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adjm[2];
    J_E[4] = dg[4]*matd[7] + loc1*dg[7];
    J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adjm[0];
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
    J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adjm[1];
    J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adjm[0];
    J_E[8] = dg[5]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[6] + matd[6];
    J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
    J_F[1] = dg[6]*matd[7] + loc1*dg[7];
    J_F[2] = dg[6]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[7] + matd[7];
    J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
    J_F[4] = dg[7]*matd[8] + loc1*dg[8];
  
    J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

    /* Assemble matrix products */

    /* dx / dx */
    A[0]  =  J_A[0]*invR[0][0] + J_A[1]*invR[0][1] + J_A[2]*invR[0][2];
    A[1]  =  J_A[1]*invR[0][0] + J_A[3]*invR[0][1] + J_A[4]*invR[0][2];
    A[2]  =  J_A[2]*invR[0][0] + J_A[4]*invR[0][1] + J_A[5]*invR[0][2];
    h_obj[0][0] =  A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dx / dy */
    A[0]  =  J_B[0]*invR[0][0] + J_B[1]*invR[0][1] + J_B[2]*invR[0][2];
    A[1]  =  J_B[3]*invR[0][0] + J_B[4]*invR[0][1] + J_B[5]*invR[0][2];
    A[2]  =  J_B[6]*invR[0][0] + J_B[7]*invR[0][1] + J_B[8]*invR[0][2];
    h_obj[0][1] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dx / dz */
    A[0]  =  J_C[0]*invR[0][0] + J_C[1]*invR[0][1] + J_C[2]*invR[0][2];
    A[1]  =  J_C[3]*invR[0][0] + J_C[4]*invR[0][1] + J_C[5]*invR[0][2];
    A[2]  =  J_C[6]*invR[0][0] + J_C[7]*invR[0][1] + J_C[8]*invR[0][2];
    h_obj[0][2] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dy / dy */
    A[0]  =  J_D[0]*invR[0][0] + J_D[1]*invR[0][1] + J_D[2]*invR[0][2];
    A[1]  =  J_D[1]*invR[0][0] + J_D[3]*invR[0][1] + J_D[4]*invR[0][2];
    A[2]  =  J_D[2]*invR[0][0] + J_D[4]*invR[0][1] + J_D[5]*invR[0][2];
    h_obj[1][1] =  A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dy / dz */
    A[0]  =  J_E[0]*invR[0][0] + J_E[1]*invR[0][1] + J_E[2]*invR[0][2];
    A[1]  =  J_E[3]*invR[0][0] + J_E[4]*invR[0][1] + J_E[5]*invR[0][2];
    A[2]  =  J_E[6]*invR[0][0] + J_E[7]*invR[0][1] + J_E[8]*invR[0][2];
    h_obj[1][2] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dz / dz */
    A[0]  =  J_F[0]*invR[0][0] + J_F[1]*invR[0][1] + J_F[2]*invR[0][2];
    A[1]  =  J_F[1]*invR[0][0] + J_F[3]*invR[0][1] + J_F[4]*invR[0][2];
    A[2]  =  J_F[2]*invR[0][0] + J_F[4]*invR[0][1] + J_F[5]*invR[0][2];
    h_obj[2][2] =  A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* Complete diagonal blocks */
    h_obj.fill_lower_triangle();
    return true;
  }

  inline bool h_gdft_2_v2(double &obj, Vector3D &g_obj, Matrix3D &h_obj,
			  const Vector3D x[3], const Vector3D &n,
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 0.5,/* constant               */
			  const Exponent& gamma = 1.0,/* simplicial elements    */
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double adjm[6], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
    double J_A[3], J_B[4], J_C[4], J_D[3], J_E[4], J_F[3];
    double A[2];

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + n[0]*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + n[1]*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + n[2]*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Calculate adjoint matrix */
    adjm[0] = dobj_df*matd[1] + dobj_dg*dg[1];
    adjm[1] = dobj_df*matd[2] + dobj_dg*dg[2];
    adjm[2] = dobj_df*matd[4] + dobj_dg*dg[4];
    adjm[3] = dobj_df*matd[5] + dobj_dg*dg[5];
    adjm[4] = dobj_df*matd[7] + dobj_dg*dg[7];
    adjm[5] = dobj_df*matd[8] + dobj_dg*dg[8];

    /* Construct gradients */
    g_obj[0] = invR[1][1]*adjm[0]+invR[1][2]*adjm[1];
    g_obj[1] = invR[1][1]*adjm[2]+invR[1][2]*adjm[3];
    g_obj[2] = invR[1][1]*adjm[4]+invR[1][2]*adjm[5];

    /* Start of the Hessian evaluation */
    adjm[0] = dobj_dg*matr[0];
    adjm[1] = dobj_dg*matr[3];
    adjm[2] = dobj_dg*matr[6];

    matd[1] *= dobj_dfdg;
    matd[2] *= dobj_dfdg;
    matd[4] *= dobj_dfdg;
    matd[5] *= dobj_dfdg;
    matd[7] *= dobj_dfdg;
    matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[1] + matd[1];
    J_A[0] = dobj_df + dg[1]*(matd[1] + loc1);
    J_A[1] = dg[1]*matd[2] + loc1*dg[2];
    J_B[0] = dg[1]*matd[4] + loc1*dg[4];
    J_B[1] = dg[1]*matd[5] + loc1*dg[5] + adjm[2];
    J_C[0] = dg[1]*matd[7] + loc1*dg[7];
    J_C[1] = dg[1]*matd[8] + loc1*dg[8] - adjm[1];

    loc1 = dobj_dgdg*dg[2] + matd[2];
    J_A[2] = dobj_df + dg[2]*(matd[2] + loc1);
    J_B[2] = dg[2]*matd[4] + loc1*dg[4] - adjm[2];
    J_B[3] = dg[2]*matd[5] + loc1*dg[5];
    J_C[2] = dg[2]*matd[7] + loc1*dg[7] + adjm[1];
    J_C[3] = dg[2]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[4] + matd[4];
    J_D[0] = dobj_df + dg[4]*(matd[4] + loc1);
    J_D[1] = dg[4]*matd[5] + loc1*dg[5];
    J_E[0] = dg[4]*matd[7] + loc1*dg[7];
    J_E[1] = dg[4]*matd[8] + loc1*dg[8] + adjm[0];
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    J_D[2] = dobj_df + dg[5]*(matd[5] + loc1);
    J_E[2] = dg[5]*matd[7] + loc1*dg[7] - adjm[0];
    J_E[3] = dg[5]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[7] + matd[7];
    J_F[0] = dobj_df + dg[7]*(matd[7] + loc1);
    J_F[1] = dg[7]*matd[8] + loc1*dg[8];
  
    J_F[2] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

    /* Assemble matrix products */

    /* dx / dx */
    A[0]  = J_A[0]*invR[1][1] + J_A[1]*invR[1][2];
    A[1]  = J_A[1]*invR[1][1] + J_A[2]*invR[1][2];
    h_obj[0][0] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dx / dy */
    A[0]  = J_B[0]*invR[1][1] + J_B[1]*invR[1][2];
    A[1]  = J_B[2]*invR[1][1] + J_B[3]*invR[1][2];
    h_obj[0][1] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dx / dz */
    A[0]  = J_C[0]*invR[1][1] + J_C[1]*invR[1][2];
    A[1]  = J_C[2]*invR[1][1] + J_C[3]*invR[1][2];
    h_obj[0][2] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dy / dy */
    A[0]  = J_D[0]*invR[1][1] + J_D[1]*invR[1][2];
    A[1]  = J_D[1]*invR[1][1] + J_D[2]*invR[1][2];
    h_obj[1][1] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dy / dz */
    A[0]  = J_E[0]*invR[1][1] + J_E[1]*invR[1][2];
    A[1]  = J_E[2]*invR[1][1] + J_E[3]*invR[1][2];
    h_obj[1][2] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dz / dz */
    A[0]  = J_F[0]*invR[1][1] + J_F[1]*invR[1][2];
    A[1]  = J_F[1]*invR[1][1] + J_F[2]*invR[1][2];
    h_obj[2][2] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* Complete diagonal blocks */
    h_obj.fill_lower_triangle();
    return true;
  }

  inline bool g_gdft_3_v0(double &obj, Vector3D &g_obj, 
			  const Vector3D x[4],
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 1.0/3.0,/* constant           */
			  const Exponent& gamma = 2.0/3.0,/* simplicial elements*/
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double adjm[9], loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Calculate adjoint matrix */
    adjm[0] = f*matd[0] + g*loc1;
    adjm[1] = f*matd[1] + g*loc2;
    adjm[2] = f*matd[2] + g*loc3;

    loc1 = g*matr[0];
    loc2 = g*matr[1];
    loc3 = g*matr[2];

    adjm[3] = f*matd[3] + loc3*matr[7] - loc2*matr[8];
    adjm[4] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
    adjm[5] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

    adjm[6] = f*matd[6] + loc2*matr[5] - loc3*matr[4];
    adjm[7] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
    adjm[8] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

    /* Construct gradients */
    g_obj[0]  = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[0] +=                    invR[1][1]*adjm[1]+invR[1][2]*adjm[2];
    g_obj[0] +=                                       invR[2][2]*adjm[2];
    g_obj[0]  = -g_obj[0];

    g_obj[1]  = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[1] +=                    invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    g_obj[1] +=                                       invR[2][2]*adjm[5];
    g_obj[1]  = -g_obj[1];

    g_obj[2]  = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    g_obj[2] +=                    invR[1][1]*adjm[7]+invR[1][2]*adjm[8];
    g_obj[2] +=                                       invR[2][2]*adjm[8];
    g_obj[2]  = -g_obj[2];
    return true;
  }

  inline bool g_gdft_3_v1(double &obj, Vector3D &g_obj,
			  const Vector3D x[4],
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 1.0/3.0,/* constant           */
			  const Exponent& gamma = 2.0/3.0,/* simplicial elements*/
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double adjm[9], loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Calculate adjoint matrix */
    adjm[0] = f*matd[0] + g*loc1;
    adjm[1] = f*matd[1] + g*loc2;
    adjm[2] = f*matd[2] + g*loc3;

    loc1 = g*matr[0];
    loc2 = g*matr[1];
    loc3 = g*matr[2];

    adjm[3] = f*matd[3] + loc3*matr[7] - loc2*matr[8];
    adjm[4] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
    adjm[5] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

    adjm[6] = f*matd[6] + loc2*matr[5] - loc3*matr[4];
    adjm[7] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
    adjm[8] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

    /* Construct gradients */
    g_obj[0] = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[1] = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[2] = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    return true;
  }

  inline bool g_gdft_3_v2(double &obj, Vector3D &g_obj, 
			  const Vector3D x[4],
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 1.0/3.0,/* constant           */
			  const Exponent& gamma = 2.0/3.0,/* simplicial elements*/
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double adjm[6], loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Calculate adjoint matrix */
    adjm[0] = f*matd[1] + g*loc2;
    adjm[1] = f*matd[2] + g*loc3;

    loc1 = g*matr[0];
    loc2 = g*matr[1];
    loc3 = g*matr[2];

    adjm[2] = f*matd[4] + loc1*matr[8] - loc3*matr[6];
    adjm[3] = f*matd[5] + loc2*matr[6] - loc1*matr[7];

    adjm[4] = f*matd[7] + loc3*matr[3] - loc1*matr[5];
    adjm[5] = f*matd[8] + loc1*matr[4] - loc2*matr[3];

    /* Construct gradients */
    g_obj[0] = invR[1][1]*adjm[0]+invR[1][2]*adjm[1];
    g_obj[1] = invR[1][1]*adjm[2]+invR[1][2]*adjm[3];
    g_obj[2] = invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    return true;
  }

  inline bool g_gdft_3_v3(double &obj, Vector3D &g_obj, 
			  const Vector3D x[4],
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 1.0/3.0,/* constant           */
			  const Exponent& gamma = 2.0/3.0,/* simplicial elements*/
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g;
    double loc1, loc2, loc3, loc4;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    loc1 = matr[4]*matr[8] - matr[5]*matr[7];
    loc2 = matr[5]*matr[6] - matr[3]*matr[8];
    loc3 = matr[3]*matr[7] - matr[4]*matr[6];
    t1 = matr[0]*loc1 + matr[1]*loc2 + matr[2]*loc3;

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = sqrt(t1*t1 + 4.0*delta*delta);
    g = t1 + t2;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc4 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc4;

    /* Calculate the derivative of the objective function. */
    f = 2.0 * loc4;
    g = -gamma * obj / t2;

    /* Construct gradients */
    loc1 = g*matr[0];
    loc2 = g*matr[1];

    g_obj[0] = invR[2][2]*(f*matd[2] + g*loc3);
    g_obj[1] = invR[2][2]*(f*matd[5] + loc2*matr[6] - loc1*matr[7]);
    g_obj[2] = invR[2][2]*(f*matd[8] + loc1*matr[4] - loc2*matr[3]);
    return true;
  }

  inline bool h_gdft_3_v0(double &obj, Vector3D &g_obj, Matrix3D &h_obj,
			  const Vector3D x[4],
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 1.0/3.0,/* constant           */
			  const Exponent& gamma = 2.0/3.0,/* simplicial elements*/
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double adjm[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
    double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];
    double A[12];

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
    dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
    dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Calculate adjoint matrix */
    adjm[0] = dobj_df*matd[0] + dobj_dg*dg[0];
    adjm[1] = dobj_df*matd[1] + dobj_dg*dg[1];
    adjm[2] = dobj_df*matd[2] + dobj_dg*dg[2];
    adjm[3] = dobj_df*matd[3] + dobj_dg*dg[3];
    adjm[4] = dobj_df*matd[4] + dobj_dg*dg[4];
    adjm[5] = dobj_df*matd[5] + dobj_dg*dg[5];
    adjm[6] = dobj_df*matd[6] + dobj_dg*dg[6];
    adjm[7] = dobj_df*matd[7] + dobj_dg*dg[7];
    adjm[8] = dobj_df*matd[8] + dobj_dg*dg[8];

    /* Construct gradients */
    g_obj[0]  = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[0] +=                    invR[1][1]*adjm[1]+invR[1][2]*adjm[2];
    g_obj[0] +=                                       invR[2][2]*adjm[2];
    g_obj[0]  = -g_obj[0];

    g_obj[1]  = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[1] +=                    invR[1][1]*adjm[4]+invR[1][2]*adjm[5];
    g_obj[1] +=                                       invR[2][2]*adjm[5];
    g_obj[1]  = -g_obj[1];

    g_obj[2]  = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];
    g_obj[2] +=                    invR[1][1]*adjm[7]+invR[1][2]*adjm[8];
    g_obj[2] +=                                       invR[2][2]*adjm[8];
    g_obj[2]  = -g_obj[2];

    /* Start of the Hessian evaluation */
    adjm[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
    adjm[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
    adjm[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
    adjm[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
    adjm[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
    adjm[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
    adjm[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
    adjm[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
    adjm[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[0] + matd[0];
    J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
    J_A[1] = dg[0]*matd[1] + loc1*dg[1];
    J_A[2] = dg[0]*matd[2] + loc1*dg[2];
    J_B[0] = dg[0]*matd[3] + loc1*dg[3];
    J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adjm[8];
    J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adjm[7];
    J_C[0] = dg[0]*matd[6] + loc1*dg[6];
    J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adjm[5];
    J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adjm[4];

    loc1 = dobj_dgdg*dg[1] + matd[1];
    J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
    J_A[4] = dg[1]*matd[2] + loc1*dg[2];
    J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adjm[8];
    J_B[4] = dg[1]*matd[4] + loc1*dg[4];
    J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adjm[6];
    J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adjm[5];
    J_C[4] = dg[1]*matd[7] + loc1*dg[7];
    J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adjm[3];

    loc1 = dobj_dgdg*dg[2] + matd[2];
    J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
    J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adjm[7];
    J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adjm[6];
    J_B[8] = dg[2]*matd[5] + loc1*dg[5];
    J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adjm[4];
    J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adjm[3];
    J_C[8] = dg[2]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[3] + matd[3];
    J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
    J_D[1] = dg[3]*matd[4] + loc1*dg[4];
    J_D[2] = dg[3]*matd[5] + loc1*dg[5];
    J_E[0] = dg[3]*matd[6] + loc1*dg[6];
    J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adjm[2];
    J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adjm[1];
  
    loc1 = dobj_dgdg*dg[4] + matd[4];
    J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
    J_D[4] = dg[4]*matd[5] + loc1*dg[5];
    J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adjm[2];
    J_E[4] = dg[4]*matd[7] + loc1*dg[7];
    J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adjm[0];
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
    J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adjm[1];
    J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adjm[0];
    J_E[8] = dg[5]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[6] + matd[6];
    J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
    J_F[1] = dg[6]*matd[7] + loc1*dg[7];
    J_F[2] = dg[6]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[7] + matd[7];
    J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
    J_F[4] = dg[7]*matd[8] + loc1*dg[8];
  
    J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

    /* Assemble matrix products */

    /* dx / dx */
    A[1]  =  J_A[0]*invR[0][0] + J_A[1]*invR[0][1] + J_A[2]*invR[0][2];
    A[2]  =                      J_A[1]*invR[1][1] + J_A[2]*invR[1][2];
    A[3]  =                                          J_A[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_A[1]*invR[0][0] + J_A[3]*invR[0][1] + J_A[4]*invR[0][2];
    A[6]  =                      J_A[3]*invR[1][1] + J_A[4]*invR[1][2];
    A[7]  =                                          J_A[4]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_A[2]*invR[0][0] + J_A[4]*invR[0][1] + J_A[5]*invR[0][2];
    A[10] =                      J_A[4]*invR[1][1] + J_A[5]*invR[1][2];
    A[11] =                                          J_A[5]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[0][0]  =  A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[0][0] +=                    A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[0][0] +=                                      A[8]*invR[2][2];
    h_obj[0][0]  = -h_obj[0][0];

    /* dx / dy */
    A[1]  =  J_B[0]*invR[0][0] + J_B[1]*invR[0][1] + J_B[2]*invR[0][2];
    A[2]  =                      J_B[1]*invR[1][1] + J_B[2]*invR[1][2];
    A[3]  =                                          J_B[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_B[3]*invR[0][0] + J_B[4]*invR[0][1] + J_B[5]*invR[0][2];
    A[6]  =                      J_B[4]*invR[1][1] + J_B[5]*invR[1][2];
    A[7]  =                                          J_B[5]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_B[6]*invR[0][0] + J_B[7]*invR[0][1] + J_B[8]*invR[0][2];
    A[10] =                      J_B[7]*invR[1][1] + J_B[8]*invR[1][2];
    A[11] =                                          J_B[8]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[0][1]  = A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[0][1] +=                   A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[0][1] +=                                     A[8]*invR[2][2];
    h_obj[0][1]  = -h_obj[0][1];

    /* dx / dz */
    A[1]  =  J_C[0]*invR[0][0] + J_C[1]*invR[0][1] + J_C[2]*invR[0][2];
    A[2]  =                      J_C[1]*invR[1][1] + J_C[2]*invR[1][2];
    A[3]  =                                          J_C[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_C[3]*invR[0][0] + J_C[4]*invR[0][1] + J_C[5]*invR[0][2];
    A[6]  =                      J_C[4]*invR[1][1] + J_C[5]*invR[1][2];
    A[7]  =                                          J_C[5]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_C[6]*invR[0][0] + J_C[7]*invR[0][1] + J_C[8]*invR[0][2];
    A[10] =                      J_C[7]*invR[1][1] + J_C[8]*invR[1][2];
    A[11] =                                          J_C[8]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[0][2]  = A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[0][2] +=                   A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[0][2] +=                                     A[8]*invR[2][2];
    h_obj[0][2]  = -h_obj[0][2];

    /* dy / dy */
    A[1]  =  J_D[0]*invR[0][0] + J_D[1]*invR[0][1] + J_D[2]*invR[0][2];
    A[2]  =                      J_D[1]*invR[1][1] + J_D[2]*invR[1][2];
    A[3]  =                                          J_D[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_D[1]*invR[0][0] + J_D[3]*invR[0][1] + J_D[4]*invR[0][2];
    A[6]  =                      J_D[3]*invR[1][1] + J_D[4]*invR[1][2];
    A[7]  =                                          J_D[4]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_D[2]*invR[0][0] + J_D[4]*invR[0][1] + J_D[5]*invR[0][2];
    A[10] =                      J_D[4]*invR[1][1] + J_D[5]*invR[1][2];
    A[11] =                                          J_D[5]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[1][1]  =  A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[1][1] +=                    A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[1][1] +=                                      A[8]*invR[2][2];
    h_obj[1][1]  = -h_obj[1][1];

    /* dy / dz */
    A[1]  =  J_E[0]*invR[0][0] + J_E[1]*invR[0][1] + J_E[2]*invR[0][2];
    A[2]  =                      J_E[1]*invR[1][1] + J_E[2]*invR[1][2];
    A[3]  =                                          J_E[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_E[3]*invR[0][0] + J_E[4]*invR[0][1] + J_E[5]*invR[0][2];
    A[6]  =                      J_E[4]*invR[1][1] + J_E[5]*invR[1][2];
    A[7]  =                                          J_E[5]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_E[6]*invR[0][0] + J_E[7]*invR[0][1] + J_E[8]*invR[0][2];
    A[10] =                      J_E[7]*invR[1][1] + J_E[8]*invR[1][2];
    A[11] =                                          J_E[8]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[1][2]  = A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[1][2] += A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[1][2] += A[8]*invR[2][2];
    h_obj[1][2]  = -h_obj[1][2];

    /* dz / dz */
    A[1]  =  J_F[0]*invR[0][0] + J_F[1]*invR[0][1] + J_F[2]*invR[0][2];
    A[2]  =                      J_F[1]*invR[1][1] + J_F[2]*invR[1][2];
    A[3]  =                                          J_F[2]*invR[2][2];
    A[0]  = -A[1] - A[2] - A[3];

    A[5]  =  J_F[1]*invR[0][0] + J_F[3]*invR[0][1] + J_F[4]*invR[0][2];
    A[6]  =                      J_F[3]*invR[1][1] + J_F[4]*invR[1][2];
    A[7]  =                                          J_F[4]*invR[2][2];
    A[4]  = -A[5] - A[6] - A[7];

    A[9]  =  J_F[2]*invR[0][0] + J_F[4]*invR[0][1] + J_F[5]*invR[0][2];
    A[10] =                      J_F[4]*invR[1][1] + J_F[5]*invR[1][2];
    A[11] =                                          J_F[5]*invR[2][2];
    A[8]  = -A[9] - A[10] - A[11];

    h_obj[2][2]  =  A[0]*invR[0][0] + A[4]*invR[0][1] + A[8]*invR[0][2];
    h_obj[2][2] +=                    A[4]*invR[1][1] + A[8]*invR[1][2];
    h_obj[2][2] +=                                      A[8]*invR[2][2];
    h_obj[2][2]  = -h_obj[2][2];

    /* Complete diagonal blocks */
    h_obj.fill_lower_triangle();
    return true;
  }

  inline bool h_gdft_3_v1(double &obj, Vector3D &g_obj, Matrix3D &h_obj,
			  const Vector3D x[4],
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 1.0/3.0,/* constant           */
			  const Exponent& gamma = 2.0/3.0,/* simplicial elements*/
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double adjm[9], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
    double J_A[6], J_B[10], J_C[10], J_D[6], J_E[10], J_F[6];
    double A[12];

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[3] = matr[2]*matr[7] - matr[1]*matr[8];
    dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[6] = matr[1]*matr[5] - matr[2]*matr[4];
    dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Calculate adjoint matrix */
    adjm[0] = dobj_df*matd[0] + dobj_dg*dg[0];
    adjm[1] = dobj_df*matd[1] + dobj_dg*dg[1];
    adjm[2] = dobj_df*matd[2] + dobj_dg*dg[2];
    adjm[3] = dobj_df*matd[3] + dobj_dg*dg[3];
    adjm[4] = dobj_df*matd[4] + dobj_dg*dg[4];
    adjm[5] = dobj_df*matd[5] + dobj_dg*dg[5];
    adjm[6] = dobj_df*matd[6] + dobj_dg*dg[6];
    adjm[7] = dobj_df*matd[7] + dobj_dg*dg[7];
    adjm[8] = dobj_df*matd[8] + dobj_dg*dg[8];

    /* Construct gradients */
    g_obj[0] = invR[0][0]*adjm[0]+invR[0][1]*adjm[1]+invR[0][2]*adjm[2];
    g_obj[1] = invR[0][0]*adjm[3]+invR[0][1]*adjm[4]+invR[0][2]*adjm[5];
    g_obj[2] = invR[0][0]*adjm[6]+invR[0][1]*adjm[7]+invR[0][2]*adjm[8];

    /* Start of the Hessian evaluation */
    adjm[0] = dobj_dg*matr[0]; matd[0] *= dobj_dfdg;
    adjm[1] = dobj_dg*matr[1]; matd[1] *= dobj_dfdg;
    adjm[2] = dobj_dg*matr[2]; matd[2] *= dobj_dfdg;
    adjm[3] = dobj_dg*matr[3]; matd[3] *= dobj_dfdg;
    adjm[4] = dobj_dg*matr[4]; matd[4] *= dobj_dfdg;
    adjm[5] = dobj_dg*matr[5]; matd[5] *= dobj_dfdg;
    adjm[6] = dobj_dg*matr[6]; matd[6] *= dobj_dfdg;
    adjm[7] = dobj_dg*matr[7]; matd[7] *= dobj_dfdg;
    adjm[8] = dobj_dg*matr[8]; matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[0] + matd[0];
    J_A[0] = dobj_df + dg[0]*(matd[0] + loc1);
    J_A[1] = dg[0]*matd[1] + loc1*dg[1];
    J_A[2] = dg[0]*matd[2] + loc1*dg[2];
    J_B[0] = dg[0]*matd[3] + loc1*dg[3];
    J_B[1] = dg[0]*matd[4] + loc1*dg[4] + adjm[8];
    J_B[2] = dg[0]*matd[5] + loc1*dg[5] - adjm[7];
    J_C[0] = dg[0]*matd[6] + loc1*dg[6];
    J_C[1] = dg[0]*matd[7] + loc1*dg[7] - adjm[5];
    J_C[2] = dg[0]*matd[8] + loc1*dg[8] + adjm[4];

    loc1 = dobj_dgdg*dg[1] + matd[1];
    J_A[3] = dobj_df + dg[1]*(matd[1] + loc1);
    J_A[4] = dg[1]*matd[2] + loc1*dg[2];
    J_B[3] = dg[1]*matd[3] + loc1*dg[3] - adjm[8];
    J_B[4] = dg[1]*matd[4] + loc1*dg[4];
    J_B[5] = dg[1]*matd[5] + loc1*dg[5] + adjm[6];
    J_C[3] = dg[1]*matd[6] + loc1*dg[6] + adjm[5];
    J_C[4] = dg[1]*matd[7] + loc1*dg[7];
    J_C[5] = dg[1]*matd[8] + loc1*dg[8] - adjm[3];

    loc1 = dobj_dgdg*dg[2] + matd[2];
    J_A[5] = dobj_df + dg[2]*(matd[2] + loc1);
    J_B[6] = dg[2]*matd[3] + loc1*dg[3] + adjm[7];
    J_B[7] = dg[2]*matd[4] + loc1*dg[4] - adjm[6];
    J_B[8] = dg[2]*matd[5] + loc1*dg[5];
    J_C[6] = dg[2]*matd[6] + loc1*dg[6] - adjm[4];
    J_C[7] = dg[2]*matd[7] + loc1*dg[7] + adjm[3];
    J_C[8] = dg[2]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[3] + matd[3];
    J_D[0] = dobj_df + dg[3]*(matd[3] + loc1);
    J_D[1] = dg[3]*matd[4] + loc1*dg[4];
    J_D[2] = dg[3]*matd[5] + loc1*dg[5];
    J_E[0] = dg[3]*matd[6] + loc1*dg[6];
    J_E[1] = dg[3]*matd[7] + loc1*dg[7] + adjm[2];
    J_E[2] = dg[3]*matd[8] + loc1*dg[8] - adjm[1];
  
    loc1 = dobj_dgdg*dg[4] + matd[4];
    J_D[3] = dobj_df + dg[4]*(matd[4] + loc1);
    J_D[4] = dg[4]*matd[5] + loc1*dg[5];
    J_E[3] = dg[4]*matd[6] + loc1*dg[6] - adjm[2];
    J_E[4] = dg[4]*matd[7] + loc1*dg[7];
    J_E[5] = dg[4]*matd[8] + loc1*dg[8] + adjm[0];
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    J_D[5] = dobj_df + dg[5]*(matd[5] + loc1);
    J_E[6] = dg[5]*matd[6] + loc1*dg[6] + adjm[1];
    J_E[7] = dg[5]*matd[7] + loc1*dg[7] - adjm[0];
    J_E[8] = dg[5]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[6] + matd[6];
    J_F[0] = dobj_df + dg[6]*(matd[6] + loc1);
    J_F[1] = dg[6]*matd[7] + loc1*dg[7];
    J_F[2] = dg[6]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[7] + matd[7];
    J_F[3] = dobj_df + dg[7]*(matd[7] + loc1);
    J_F[4] = dg[7]*matd[8] + loc1*dg[8];
  
    J_F[5] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

    /* Assemble matrix products */

    /* dx / dx */
    A[0] = J_A[0]*invR[0][0] + J_A[1]*invR[0][1] + J_A[2]*invR[0][2];
    A[1] = J_A[1]*invR[0][0] + J_A[3]*invR[0][1] + J_A[4]*invR[0][2];
    A[2] = J_A[2]*invR[0][0] + J_A[4]*invR[0][1] + J_A[5]*invR[0][2];
    h_obj[0][0] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dx / dy */
    A[0] = J_B[0]*invR[0][0] + J_B[1]*invR[0][1] + J_B[2]*invR[0][2];
    A[1] = J_B[3]*invR[0][0] + J_B[4]*invR[0][1] + J_B[5]*invR[0][2];
    A[2] = J_B[6]*invR[0][0] + J_B[7]*invR[0][1] + J_B[8]*invR[0][2];
    h_obj[0][1] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dx / dz */
    A[0] = J_C[0]*invR[0][0] + J_C[1]*invR[0][1] + J_C[2]*invR[0][2];
    A[1] = J_C[3]*invR[0][0] + J_C[4]*invR[0][1] + J_C[5]*invR[0][2];
    A[2] = J_C[6]*invR[0][0] + J_C[7]*invR[0][1] + J_C[8]*invR[0][2];
    h_obj[0][2] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dy / dy */
    A[0] = J_D[0]*invR[0][0] + J_D[1]*invR[0][1] + J_D[2]*invR[0][2];
    A[1] = J_D[1]*invR[0][0] + J_D[3]*invR[0][1] + J_D[4]*invR[0][2];
    A[2] = J_D[2]*invR[0][0] + J_D[4]*invR[0][1] + J_D[5]*invR[0][2];
    h_obj[1][1] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dy / dz */
    A[0] = J_E[0]*invR[0][0] + J_E[1]*invR[0][1] + J_E[2]*invR[0][2];
    A[1] = J_E[3]*invR[0][0] + J_E[4]*invR[0][1] + J_E[5]*invR[0][2];
    A[2] = J_E[6]*invR[0][0] + J_E[7]*invR[0][1] + J_E[8]*invR[0][2];
    h_obj[1][2] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* dz / dz */
    A[0] = J_F[0]*invR[0][0] + J_F[1]*invR[0][1] + J_F[2]*invR[0][2];
    A[1] = J_F[1]*invR[0][0] + J_F[3]*invR[0][1] + J_F[4]*invR[0][2];
    A[2] = J_F[2]*invR[0][0] + J_F[4]*invR[0][1] + J_F[5]*invR[0][2];
    h_obj[2][2] = A[0]*invR[0][0] + A[1]*invR[0][1] + A[2]*invR[0][2];

    /* Complete diagonal blocks */
    h_obj.fill_lower_triangle();
    return true;
  }

  inline bool h_gdft_3_v2(double &obj, Vector3D &g_obj, Matrix3D &h_obj,
			  const Vector3D x[4],
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 1.0/3.0,/* constant           */
			  const Exponent& gamma = 2.0/3.0,/* simplicial elements*/
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double adjm[6], dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;
    double J_A[3], J_B[4], J_C[4], J_D[3], J_E[4], J_F[3];
    double A[2];

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[4] = matr[0]*matr[8] - matr[2]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[7] = matr[2]*matr[3] - matr[0]*matr[5];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Calculate adjoint matrix */
    adjm[0] = dobj_df*matd[1] + dobj_dg*dg[1];
    adjm[1] = dobj_df*matd[2] + dobj_dg*dg[2];
    adjm[2] = dobj_df*matd[4] + dobj_dg*dg[4];
    adjm[3] = dobj_df*matd[5] + dobj_dg*dg[5];
    adjm[4] = dobj_df*matd[7] + dobj_dg*dg[7];
    adjm[5] = dobj_df*matd[8] + dobj_dg*dg[8];

    /* Construct gradients */
    g_obj[0] = invR[1][1]*adjm[0]+invR[1][2]*adjm[1];
    g_obj[1] = invR[1][1]*adjm[2]+invR[1][2]*adjm[3];
    g_obj[2] = invR[1][1]*adjm[4]+invR[1][2]*adjm[5];

    /* Start of the Hessian evaluation */
    adjm[0] = dobj_dg*matr[0];
    adjm[1] = dobj_dg*matr[3];
    adjm[2] = dobj_dg*matr[6];

    matd[1] *= dobj_dfdg;
    matd[2] *= dobj_dfdg;
    matd[4] *= dobj_dfdg;
    matd[5] *= dobj_dfdg;
    matd[7] *= dobj_dfdg;
    matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[1] + matd[1];
    J_A[0] = dobj_df + dg[1]*(matd[1] + loc1);
    J_A[1] = dg[1]*matd[2] + loc1*dg[2];
    J_B[0] = dg[1]*matd[4] + loc1*dg[4];
    J_B[1] = dg[1]*matd[5] + loc1*dg[5] + adjm[2];
    J_C[0] = dg[1]*matd[7] + loc1*dg[7];
    J_C[1] = dg[1]*matd[8] + loc1*dg[8] - adjm[1];

    loc1 = dobj_dgdg*dg[2] + matd[2];
    J_A[2] = dobj_df + dg[2]*(matd[2] + loc1);
    J_B[2] = dg[2]*matd[4] + loc1*dg[4] - adjm[2];
    J_B[3] = dg[2]*matd[5] + loc1*dg[5];
    J_C[2] = dg[2]*matd[7] + loc1*dg[7] + adjm[1];
    J_C[3] = dg[2]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[4] + matd[4];
    J_D[0] = dobj_df + dg[4]*(matd[4] + loc1);
    J_D[1] = dg[4]*matd[5] + loc1*dg[5];
    J_E[0] = dg[4]*matd[7] + loc1*dg[7];
    J_E[1] = dg[4]*matd[8] + loc1*dg[8] + adjm[0];
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    J_D[2] = dobj_df + dg[5]*(matd[5] + loc1);
    J_E[2] = dg[5]*matd[7] + loc1*dg[7] - adjm[0];
    J_E[3] = dg[5]*matd[8] + loc1*dg[8];
  
    loc1 = dobj_dgdg*dg[7] + matd[7];
    J_F[0] = dobj_df + dg[7]*(matd[7] + loc1);
    J_F[1] = dg[7]*matd[8] + loc1*dg[8];
  
    J_F[2] = dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]);

    /* Assemble matrix products */

    /* dx / dx */
    A[0] = J_A[0]*invR[1][1] + J_A[1]*invR[1][2];
    A[1] = J_A[1]*invR[1][1] + J_A[2]*invR[1][2];
    h_obj[0][0] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dx / dy */
    A[0] = J_B[0]*invR[1][1] + J_B[1]*invR[1][2];
    A[1] = J_B[2]*invR[1][1] + J_B[3]*invR[1][2];
    h_obj[0][1] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dx / dz */
    A[0] = J_C[0]*invR[1][1] + J_C[1]*invR[1][2];
    A[1] = J_C[2]*invR[1][1] + J_C[3]*invR[1][2];
    h_obj[0][2] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dy / dy */
    A[0] = J_D[0]*invR[1][1] + J_D[1]*invR[1][2];
    A[1] = J_D[1]*invR[1][1] + J_D[2]*invR[1][2];
    h_obj[1][1] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dy / dz */
    A[0] = J_E[0]*invR[1][1] + J_E[1]*invR[1][2];
    A[1] = J_E[2]*invR[1][1] + J_E[3]*invR[1][2];
    h_obj[1][2] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* dz / dz */
    A[0] = J_F[0]*invR[1][1] + J_F[1]*invR[1][2];
    A[1] = J_F[1]*invR[1][1] + J_F[2]*invR[1][2];
    h_obj[2][2] = A[0]*invR[1][1] + A[1]*invR[1][2];

    /* Complete diagonal blocks */
    h_obj.fill_lower_triangle();
    return true;
  }

  inline bool h_gdft_3_v3(double &obj, Vector3D &g_obj, Matrix3D &h_obj,
			  const Vector3D x[4],
			  const Matrix3D &invR,	/* upper triangular          */
			  const Matrix3D &Q, 	/* orthogonal, det(Q) = 1    */
			  const double alpha = 1.0/3.0,/* constant           */
			  const Exponent& gamma = 2.0/3.0,/* simplicial elements*/
			  const double delta = 0.0,/* max in denominator     */
			  const double beta  = 0.0)/* no modification        */
  {
    double matr[9], f, t1, t2;
    double matd[9], g, t3, loc1;
    double dg[9], dobj_df, dobj_dfdg, dobj_dg, dobj_dgdg;

    /* Calculate M = A*inv(R). */
    f       = x[1][0] - x[0][0];
    g       = x[2][0] - x[0][0];
    t1      = x[3][0] - x[0][0];
    matr[0] = f*invR[0][0];
    matr[1] = f*invR[0][1] + g*invR[1][1];
    matr[2] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][1] - x[0][1];
    g       = x[2][1] - x[0][1];
    t1      = x[3][1] - x[0][1];
    matr[3] = f*invR[0][0];
    matr[4] = f*invR[0][1] + g*invR[1][1];
    matr[5] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    f       = x[1][2] - x[0][2];
    g       = x[2][2] - x[0][2];
    t1      = x[3][2] - x[0][2];
    matr[6] = f*invR[0][0];
    matr[7] = f*invR[0][1] + g*invR[1][1];
    matr[8] = f*invR[0][2] + g*invR[1][2] + t1*invR[2][2];

    /* Calculate det(M). */
    dg[0] = matr[4]*matr[8] - matr[5]*matr[7];
    dg[1] = matr[5]*matr[6] - matr[3]*matr[8];
    dg[2] = matr[3]*matr[7] - matr[4]*matr[6];
    dg[5] = matr[1]*matr[6] - matr[0]*matr[7];
    dg[8] = matr[0]*matr[4] - matr[1]*matr[3];

    t1 = matr[0]*dg[0] + matr[1]*dg[1] + matr[2]*dg[2];

    if ((0.0 == delta) && (t1 < MSQ_MIN)) { obj = t1; return false; }

    /* Calculate sqrt(det(M)^2 + 4*delta^2) and denominator. */
    t2 = t1*t1 + 4.0*delta*delta;
    t3 = sqrt(t2);
    g = t1 + t3;
    
    /* Calculate N = M - beta*Q. */
    matd[0] = matr[0] - beta*Q[0][0];
    matd[1] = matr[1] - beta*Q[0][1];
    matd[2] = matr[2] - beta*Q[0][2];
    matd[3] = matr[3] - beta*Q[1][0];
    matd[4] = matr[4] - beta*Q[1][1];
    matd[5] = matr[5] - beta*Q[1][2];
    matd[6] = matr[6] - beta*Q[2][0];
    matd[7] = matr[7] - beta*Q[2][1];
    matd[8] = matr[8] - beta*Q[2][2];

    /* Calculate norm(N) */
    f = matd[0]*matd[0] + matd[1]*matd[1] + matd[2]*matd[2] +
        matd[3]*matd[3] + matd[4]*matd[4] + matd[5]*matd[5] +
        matd[6]*matd[6] + matd[7]*matd[7] + matd[8]*matd[8];

    /* Calculate objective function. */
    loc1 = alpha * pow(2.0, gamma) / pow(g, gamma);
    obj = f * loc1;

    /* Calculate the derivative of the objective function. */
    t3 = 1.0 / t3;
    dobj_df = 2.0 * loc1;
    dobj_dg = -gamma * obj * t3;
    dobj_dfdg = -gamma * dobj_df * t3;
    dobj_dgdg = dobj_dg * ((-gamma - 1.0)*t3 + 4.0*delta*delta/(t2*g));

    /* Construct gradients */
    g_obj[0] = invR[2][2]*(dobj_df*matd[2] + dobj_dg*dg[2]);
    g_obj[1] = invR[2][2]*(dobj_df*matd[5] + dobj_dg*dg[5]);
    g_obj[2] = invR[2][2]*(dobj_df*matd[8] + dobj_dg*dg[8]);

    /* Start of the Hessian evaluation */
    t1 = invR[2][2]*invR[2][2];
    matd[2] *= dobj_dfdg;
    matd[5] *= dobj_dfdg;
    matd[8] *= dobj_dfdg;

    /* Blocks for the Hessian construction */
    loc1 = dobj_dgdg*dg[2] + matd[2];
    h_obj[0][0] = t1*(dobj_df + dg[2]*(matd[2] + loc1));
    h_obj[0][1] = t1*(dg[2]*matd[5] + loc1*dg[5]);
    h_obj[0][2] = t1*(dg[2]*matd[8] + loc1*dg[8]);
  
    loc1 = dobj_dgdg*dg[5] + matd[5];
    h_obj[1][1] = t1*(dobj_df + dg[5]*(matd[5] + loc1));
    h_obj[1][2] = t1*(dg[5]*matd[8] + loc1*dg[8]);
  
    h_obj[2][2] = t1*(dobj_df + dg[8]*(2.0*matd[8] + dobj_dgdg*dg[8]));

    /* Complete diagonal blocks */
    h_obj.fill_lower_triangle();
    return true;
  }
}
#endif


