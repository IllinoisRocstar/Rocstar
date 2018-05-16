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
  \file   sRI_DFT.cpp
  \brief  

  \author Thomas Leurent

  \date   2004-04-12
*/

#include "sRI_DFT.hpp"

using namespace Mesquite;


bool sRI_DFT::evaluate_element(PatchData& pd,
                              MsqMeshEntity* element,
                              double& value, MsqError &err)
{
  Matrix3D T[MSQ_MAX_NUM_VERT_PER_ENT];
  double c_k[MSQ_MAX_NUM_VERT_PER_ENT];
  double dft[MSQ_MAX_NUM_VERT_PER_ENT];
  bool return_flag = false;
  double h, tau, s;
    
  size_t num_T = element->vertex_count();
  compute_T_matrices(*element, pd, T, num_T, c_k, err); MSQ_ERRZERO(err);

  const double id[] = {1., 0., 0.,  0., 1., 0.,  0., 0., 1.};
  const Matrix3D I(id);
  Matrix3D TT;
  for (size_t i=0; i<num_T; ++i) {
    tau = det(T[i]);
    s = pow(tau, MSQ_ONE_THIRD); // for 3D (tet & hex)
    TT = transpose(T[i]);
    TT = TT * T[i];
    TT -= (s*s)*I; 
    dft[i] = .5 * Frobenius_2(TT);
    return_flag = get_barrier_function(pd, tau, h, err); MSQ_ERRZERO(err);
    dft[i] /= pow(h, 4.0/3.0);
  }
    
  value = weighted_average_metrics(c_k, dft, num_T, err); MSQ_ERRZERO(err);
    
  return return_flag;
}



