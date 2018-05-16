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
  \file   IdealWeightInverseMeanRatio.cpp
  \brief  

  \author Michael Brewer
  \author Todd Munson
  \author Thomas Leurent

  \date   2002-06-9
*/
#include "IdealWeightInverseMeanRatio.hpp"
#include "MeanRatioFunctions.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
   using std::vector;
#endif

#include <math.h>

namespace Mesquite {

IdealWeightInverseMeanRatio::IdealWeightInverseMeanRatio( MsqError& err, double pow_dbl )
{
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_ERRRTN(err);

  set_gradient_type(ANALYTICAL_GRADIENT);
  set_hessian_type(ANALYTICAL_HESSIAN);
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  set_name("Inverse Mean Ratio");

  set_metric_power(pow_dbl, err);  MSQ_ERRRTN(err);
}

       //! Sets the power value in the metric computation.
void IdealWeightInverseMeanRatio::set_metric_power(double pow_dbl, MsqError& err)
{
  if(fabs(pow_dbl)<=MSQ_MIN){
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    return;
  }
  if(pow_dbl<0)
    set_negate_flag(-1);
  else
    set_negate_flag(1);
  a2Con=pow(.5,pow_dbl);
  b2Con=pow_dbl;
  c2Con=-pow_dbl;
  a3Con=pow(1.0/3.0,pow_dbl);
  b3Con=pow_dbl;
  c3Con=-2.0*pow_dbl/3.0;
}


bool IdealWeightInverseMeanRatio::evaluate_element(PatchData &pd,
                                              MsqMeshEntity *e,
                                              double &m,
                                              MsqError &err)
{
  EntityTopology topo = e->get_element_type();

  MsqVertex *vertices = pd.get_vertex_array(err);  MSQ_ERRZERO(err);
  const size_t *v_i = e->get_vertex_index_array();

  Vector3D n;			// Surface normal for 2D objects

  // Hex element descriptions
  static const int locs_hex[8][4] = {{0, 1, 3, 4},
				     {1, 2, 0, 5},
				     {2, 3, 1, 6},
				     {3, 0, 2, 7},
				     {4, 7, 5, 0},
				     {5, 4, 6, 1},
				     {6, 5, 7, 2},
				     {7, 6, 4, 3}};

  const Vector3D d_con(1.0, 1.0, 1.0);

  int i;

  m = 0.0;
  bool metric_valid = false;
  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n = n / n.length();		// Need unit normal
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    metric_valid = m_fcn_2e(m, mCoords, n, a2Con, b2Con, c2Con);
    if (!metric_valid) return false;
    break;
    
  case QUADRILATERAL:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n = n / n.length();	// Need unit normal
    for (i = 0; i < 4; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      metric_valid = m_fcn_2i(mMetrics[i], mCoords, n, 
			      a2Con, b2Con, c2Con, d_con);
      if (!metric_valid) return false;
    }
    m = average_metrics(mMetrics, 4, err);
    break;

  case TETRAHEDRON:
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    mCoords[3] = vertices[v_i[3]];
    metric_valid = m_fcn_3e(m, mCoords, a3Con, b3Con, c3Con);
    if (!metric_valid) return false;
    break;

  case HEXAHEDRON:
    for (i = 0; i < 8; ++i) {
      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      mCoords[3] = vertices[v_i[locs_hex[i][3]]];
      metric_valid = m_fcn_3i(mMetrics[i], mCoords, 
			      a3Con, b3Con, c3Con, d_con);
      if (!metric_valid) return false;
    }
    m = average_metrics(mMetrics, 8, err); MSQ_ERRZERO(err);
    break;

  default:
    break;
  } // end switch over element type
  return true;
}

bool IdealWeightInverseMeanRatio::compute_element_analytical_gradient(PatchData &pd,
								 MsqMeshEntity *e,
								 MsqVertex *v[], 
								 Vector3D g[],
								 int nv, 
								 double &m,
                         MsqError &err)
{
  EntityTopology topo = e->get_element_type();

  if (((topo == QUADRILATERAL) || (topo == HEXAHEDRON)) && 
      ((avgMethod == MINIMUM) || (avgMethod == MAXIMUM))) {
    MSQ_DBGOUT(1) <<
      "Minimum and maximum not continuously differentiable.\n"
      "Element of subdifferential will be returned.\n";
  }

  MsqVertex *vertices = pd.get_vertex_array(err);  MSQ_ERRZERO(err);
  const size_t *v_i = e->get_vertex_index_array();

  Vector3D n;			// Surface normal for 2D objects

  // Hex element descriptions
  static const int locs_hex[8][4] = {{0, 1, 3, 4},
				     {1, 2, 0, 5},
				     {2, 3, 1, 6},
				     {3, 0, 2, 7},
				     {4, 7, 5, 0},
				     {5, 4, 6, 1},
				     {6, 5, 7, 2},
				     {7, 6, 4, 3}};

  const Vector3D d_con(1.0, 1.0, 1.0);

  int i, j;

  bool metric_valid = false;
  int vert_per_elem = 0;

  m = 0.0;

  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n /= n.length();		// Need unit normal
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    if (!g_fcn_2e(m, mAccumGrad, mCoords, n, a2Con, b2Con, c2Con)) return false;

    vert_per_elem = 3;
    break;

  case QUADRILATERAL:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n /= n.length();	// Need unit normal
    for (i = 0; i < 4; ++i) {
      mAccumGrad[i] = 0.0;

      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      if (!g_fcn_2i(mMetrics[i], mGradients+3*i, mCoords, n,
		    a2Con, b2Con, c2Con, d_con)) return false;
    }

    m = average_metric_and_weights( mMetrics, 4, err ); MSQ_ERRZERO(err);
    for (i = 0; i < 4; ++i) {
      mAccumGrad[locs_hex[i][0]] += mMetrics[i]*mGradients[3*i+0];
      mAccumGrad[locs_hex[i][1]] += mMetrics[i]*mGradients[3*i+1];
      mAccumGrad[locs_hex[i][2]] += mMetrics[i]*mGradients[3*i+2];
    }

    vert_per_elem = 4;
    break;

  case TETRAHEDRON:
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    mCoords[3] = vertices[v_i[3]];
    metric_valid = g_fcn_3e(m, mAccumGrad, mCoords, a3Con, b3Con, c3Con);
    if (!metric_valid) return false;

    vert_per_elem = 4;
    break;

  case HEXAHEDRON:
    for (i = 0; i < 8; ++i) {
      mAccumGrad[i] = 0.0;

      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      mCoords[3] = vertices[v_i[locs_hex[i][3]]];
      if (!g_fcn_3i(mMetrics[i], mGradients+4*i, mCoords, 
		    a3Con, b3Con, c3Con, d_con)) return false;
    }

    m = average_metric_and_weights( mMetrics, 8, err ); MSQ_ERRZERO(err);
    for (i = 0; i < 8; ++i) {
      mAccumGrad[locs_hex[i][0]] += mMetrics[i]*mGradients[4*i+0];
      mAccumGrad[locs_hex[i][1]] += mMetrics[i]*mGradients[4*i+1];
      mAccumGrad[locs_hex[i][2]] += mMetrics[i]*mGradients[4*i+2];
      mAccumGrad[locs_hex[i][3]] += mMetrics[i]*mGradients[4*i+3];
    }
 
    vert_per_elem = 8;
    break;

  default:
    break;
  } // end switch over element type


  // This is not very efficient, but is one way to select correct gradients
  // For gradients, info is returned only for free vertices, in the order of v[].
  for (i = 0; i < vert_per_elem; ++i) {
    for (j = 0; j < nv; ++j) {
      if (vertices + v_i[i] == v[j]) {
        g[j] = mAccumGrad[i];
      }
    }
  }
  

  return true;
}


bool IdealWeightInverseMeanRatio::compute_element_analytical_hessian(PatchData &pd,
								MsqMeshEntity *e,
								MsqVertex *fv[], 
								Vector3D g[],
								Matrix3D h[],
                                                                int nfv, 
								double &m,
								MsqError &err)
{
  EntityTopology topo = e->get_element_type();

  if (((topo == QUADRILATERAL) || (topo == HEXAHEDRON)) && 
      ((avgMethod == MINIMUM) || (avgMethod == MAXIMUM))) {
    MSQ_DBGOUT(1) <<
      "Minimum and maximum not continuously differentiable.\n"
      "Element of subdifferential will be returned.\n"
      "Who knows what the Hessian is?\n" ;
  }

  MsqVertex *vertices = pd.get_vertex_array(err);  MSQ_ERRZERO(err);
  const size_t *v_i = e->get_vertex_index_array();


  Vector3D n;			// Surface normal for 2D objects
  Matrix3D outer;
  double   nm, t=0;

  // Hex element descriptions
  static const int locs_hex[8][4] = {{0, 1, 3, 4},  
				     {1, 2, 0, 5},
				     {2, 3, 1, 6},
				     {3, 0, 2, 7},
				     {4, 7, 5, 0},
				     {5, 4, 6, 1},
				     {6, 5, 7, 2},
				     {7, 6, 4, 3}};

  const Vector3D d_con(1.0, 1.0, 1.0);

  int i, j, k, l, ind;
  int r, c, loc;

  bool metric_valid = false;

  m = 0.0;

  switch(topo) {
  case TRIANGLE:
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n /= n.length();		// Need unit normal
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    if (!h_fcn_2e(m, g, h, mCoords, n, a2Con, b2Con, c2Con)) return false;

    // zero out fixed elements of g
    j = 0;
    for (i = 0; i < 3; ++i) {
      // if free vertex, see next
      if (j<nfv && vertices + v_i[i] == fv[j] )
        ++j;
      // else zero gradient and Hessian entries
      else {
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
        }
      }
    }
    break;

  case QUADRILATERAL:
    for (i=0; i < 10; ++i) {
      h[i].zero();
    }
    
    pd.get_domain_normal_at_element(e, n, err); MSQ_ERRZERO(err);
    n /= n.length();	// Need unit normal
    for (i = 0; i < 4; ++i) {
      g[i] = 0.0;

      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      if (!h_fcn_2i(mMetrics[i], mGradients+3*i, mHessians+6*i, mCoords, n,
		    a2Con, b2Con, c2Con, d_con)) return false;
    }

    switch(avgMethod) {
    case MINIMUM:
      MSQ_SETERR(err)("MINIMUM averaging method does not work.",MsqError::INVALID_STATE);
      return false;

    case MAXIMUM:
      MSQ_SETERR(err)("MAXIMUM averaging method does not work.",MsqError::INVALID_STATE);
      return false;

    case SUM:
      m = 0;
      for (i = 0; i < 4; ++i) {
        m += mMetrics[i];
      }

      l = 0;
      for (i = 0; i < 4; ++i) {
        g[locs_hex[i][0]] += mGradients[3*i+0];
        g[locs_hex[i][1]] += mGradients[3*i+1];
        g[locs_hex[i][2]] += mGradients[3*i+2];

        for (j = 0; j < 3; ++j) {
          for (k = j; k < 3; ++k) {
            r = locs_hex[i][j];
            c = locs_hex[i][k];

            if (r <= c) {
              loc = 4*r - (r*(r+1)/2) + c;
              h[loc] += mHessians[l];
            } 
            else {
              loc = 4*c - (c*(c+1)/2) + r;
              h[loc] += transpose(mHessians[l]);
            }
            ++l;
          }
        }
      }
      break;

    case SUM_SQUARED:
      m = 0;
      for (i = 0; i < 4; ++i) {
        m += (mMetrics[i]*mMetrics[i]);
	mMetrics[i] *= 2;
      }

      l = 0;
      for (i = 0; i < 4; ++i) {
        g[locs_hex[i][0]] += mMetrics[i]*mGradients[3*i+0];
        g[locs_hex[i][1]] += mMetrics[i]*mGradients[3*i+1];
        g[locs_hex[i][2]] += mMetrics[i]*mGradients[3*i+2];

        for (j = 0; j < 3; ++j) {
          for (k = j; k < 3; ++k) {
	    outer = 2.0*outer.outer_product(mGradients[3*i+j], 
					    mGradients[3*i+k]);
	    
            r = locs_hex[i][j];
            c = locs_hex[i][k];

            if (r <= c) {
              loc = 4*r - (r*(r+1)/2) + c;
              h[loc] += mMetrics[i]*mHessians[l] + outer;
            } 
            else {
              loc = 4*c - (c*(c+1)/2) + r;
              h[loc] += transpose(mMetrics[i]*mHessians[l] + outer);
            }
            ++l;
          }
        }
      }
      break;

    case LINEAR:
      m = 0;
      for (i = 0; i < 4; ++i) {
        m += mMetrics[i];
      }
      m *= 0.25;

      l = 0;
      for (i = 0; i < 4; ++i) {
        g[locs_hex[i][0]] += 0.25*mGradients[3*i+0];
        g[locs_hex[i][1]] += 0.25*mGradients[3*i+1];
        g[locs_hex[i][2]] += 0.25*mGradients[3*i+2];

        for (j = 0; j < 3; ++j) {
          for (k = j; k < 3; ++k) {
            r = locs_hex[i][j];
            c = locs_hex[i][k];

            if (r <= c) {
              loc = 4*r - (r*(r+1)/2) + c;
              h[loc] += mHessians[l];
            } 
            else {
              loc = 4*c - (c*(c+1)/2) + r;
              h[loc] += transpose(mHessians[l]);
            }
            ++l;
          }
        }
      }

      for (i=0; i<10; ++i) {
        h[i] *= 0.25;
      }
      break;

    case GEOMETRIC:
      MSQ_SETERR(err)("GEOMETRIC averaging method does not work.",MsqError::INVALID_STATE);
      return false;

    default:
      switch(avgMethod) {
      case RMS:
	t = 2.0;
	break;

      case HARMONIC:
	t = -1.0;
	break;

      case HMS:
	t = -2.0;
	break;

      default:
        MSQ_SETERR(err)("averaging method not available.",MsqError::NOT_IMPLEMENTED);
        break;
      }

      m = 0;
      for (i = 0; i < 4; ++i) {
	nm = pow(mMetrics[i], t);
	m += nm;

	g_factor[i] = 0.25*t*nm / mMetrics[i];
	h_factor[i] = (t-1)*g_factor[i] / mMetrics[i];
      }

      nm = 0.25 * m;

      l = 0;
      for (i = 0; i < 4; ++i) {
        g[locs_hex[i][0]] += g_factor[i]*mGradients[3*i+0];
	g[locs_hex[i][1]] += g_factor[i]*mGradients[3*i+1];
	g[locs_hex[i][2]] += g_factor[i]*mGradients[3*i+2];

        for (j = 0; j < 3; ++j) {
          for (k = j; k < 3; ++k) {
	    outer = h_factor[i] * outer.outer_product(mGradients[3*i+j], 
						      mGradients[3*i+k]);
	    
            r = locs_hex[i][j];
            c = locs_hex[i][k];

            if (r <= c) {
              loc = 4*r - (r*(r+1)/2) + c;
              h[loc] += g_factor[i]*mHessians[l] + outer;
            } 
            else {
              loc = 4*c - (c*(c+1)/2) + r;
              h[loc] += transpose(g_factor[i]*mHessians[l] + outer);
            }
            ++l;
          }
        }
      }

      m = pow(nm, 1.0 / t);
      g_factor[0] = m / (t*nm);
      h_factor[0] = (1.0 / t - 1)*g_factor[0] / nm;

      l = 0;
      for (i = 0; i < 4; ++i) {
	for (j = i; j < 4; ++j) {
	  outer = outer.outer_product(g[i], g[j]);
	  h[l] = g_factor[0]*h[l] + h_factor[0]*outer;
	  ++l;
	}
	g[i] *= g_factor[0];
      }
      break;
    }

    // zero out fixed elements of gradient and Hessian
    ind = 0;
    for (i=0; i<4; ++i) {
      // if free vertex, see next
      if (ind<nfv && vertices+v_i[i] == fv[ind] )
        ++ind;
      // else zero gradient entry and hessian entries.
      else {
        g[i] = 0.;
        switch(i) {
        case 0:
          h[0].zero();   h[1].zero();   h[2].zero();   h[3].zero();
          break;
          
        case 1:
          h[1].zero();   h[4].zero();   h[5].zero();   h[6].zero();
          break;
          
        case 2:
          h[2].zero();   h[5].zero();   h[7].zero();   h[8].zero();
          break;
          
        case 3:
          h[3].zero();   h[6].zero();   h[8].zero();   h[9].zero();
          break;
        }
      }
    }
    break;

  case TETRAHEDRON:
    mCoords[0] = vertices[v_i[0]];
    mCoords[1] = vertices[v_i[1]];
    mCoords[2] = vertices[v_i[2]];
    mCoords[3] = vertices[v_i[3]];
    metric_valid = h_fcn_3e(m, g, h, mCoords, a3Con, b3Con, c3Con);
    if (!metric_valid) return false;

    // zero out fixed elements of g
    j = 0;
    for (i = 0; i < 4; ++i) {
      // if free vertex, see next
      if (j<nfv && vertices + v_i[i] == fv[j] )
        ++j;
      // else zero gradient entry
      else {
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
    for (i=0; i<36; ++i)
      h[i].zero();

    for (i = 0; i < 8; ++i) {
      g[i] = 0.0;

      mCoords[0] = vertices[v_i[locs_hex[i][0]]];
      mCoords[1] = vertices[v_i[locs_hex[i][1]]];
      mCoords[2] = vertices[v_i[locs_hex[i][2]]];
      mCoords[3] = vertices[v_i[locs_hex[i][3]]];
      if (!h_fcn_3i(mMetrics[i], mGradients+4*i, mHessians+10*i, mCoords,
		    a3Con, b3Con, c3Con, d_con)) return false;
    }

    switch(avgMethod) {
    case MINIMUM:
      MSQ_SETERR(err)("MINIMUM averaging method does not work.",MsqError::NOT_IMPLEMENTED);
      return false;

    case MAXIMUM:
      MSQ_SETERR(err)("MAXIMUM averaging method does not work.",MsqError::NOT_IMPLEMENTED);
      return false;

    case SUM:
      m = 0;
      for (i = 0; i < 8; ++i) {
        m += mMetrics[i];
      }

      l = 0;
      for (i = 0; i < 8; ++i) {
        g[locs_hex[i][0]] += mGradients[4*i+0];
        g[locs_hex[i][1]] += mGradients[4*i+1];
        g[locs_hex[i][2]] += mGradients[4*i+2];
        g[locs_hex[i][3]] += mGradients[4*i+3];

        for (j = 0; j < 4; ++j) {
          for (k = j; k < 4; ++k) {
            r = locs_hex[i][j];
            c = locs_hex[i][k];

            if (r <= c) {
              loc = 8*r - (r*(r+1)/2) + c;
              h[loc] += mHessians[l];
            } 
            else {
              loc = 8*c - (c*(c+1)/2) + r;
              h[loc] += transpose(mHessians[l]);
            }
            ++l;
          }
        }
      }
      break;

    case SUM_SQUARED:
      m = 0;
      for (i = 0; i < 8; ++i) {
        m += (mMetrics[i]*mMetrics[i]);
	mMetrics[i] *= 2.0;
      }

      l = 0;
      for (i = 0; i < 8; ++i) {
        g[locs_hex[i][0]] += mMetrics[i]*mGradients[4*i+0];
        g[locs_hex[i][1]] += mMetrics[i]*mGradients[4*i+1];
        g[locs_hex[i][2]] += mMetrics[i]*mGradients[4*i+2];
        g[locs_hex[i][3]] += mMetrics[i]*mGradients[4*i+3];

        for (j = 0; j < 4; ++j) {
          for (k = j; k < 4; ++k) {
	    outer = 2.0*outer.outer_product(mGradients[4*i+j], 
					    mGradients[4*i+k]);

            r = locs_hex[i][j];
            c = locs_hex[i][k];

            if (r <= c) {
              loc = 8*r - (r*(r+1)/2) + c;
              h[loc] += mMetrics[i]*mHessians[l] + outer;
            } 
            else {
              loc = 8*c - (c*(c+1)/2) + r;
              h[loc] += transpose(mMetrics[i]*mHessians[l] + outer);
            }
            ++l;
          }
        }
      }
      break;

    case LINEAR:
      m = 0;
      for (i = 0; i < 8; ++i) {
        m += mMetrics[i];
      }
      m *= 0.125;

      l = 0;
      for (i = 0; i < 8; ++i) {
        g[locs_hex[i][0]] += 0.125*mGradients[4*i+0];
        g[locs_hex[i][1]] += 0.125*mGradients[4*i+1];
        g[locs_hex[i][2]] += 0.125*mGradients[4*i+2];
        g[locs_hex[i][3]] += 0.125*mGradients[4*i+3];

        for (j = 0; j < 4; ++j) {
          for (k = j; k < 4; ++k) {
            r = locs_hex[i][j];
            c = locs_hex[i][k];

            if (r <= c) {
              loc = 8*r - (r*(r+1)/2) + c;
              h[loc] += mHessians[l];
            } 
            else {
              loc = 8*c - (c*(c+1)/2) + r;
              h[loc] += transpose(mHessians[l]);
            }
            ++l;
          }
        }
      }

      for (i=0; i<36; ++i)
        h[i] *= 0.125;
      break;

    case GEOMETRIC:
      MSQ_SETERR(err)("GEOMETRIC averaging method does not work.",MsqError::NOT_IMPLEMENTED);
      return false;

    default:
      switch(avgMethod) {
      case RMS:
	t = 2.0;
	break;

      case HARMONIC:
	t = -1.0;
	break;

      case HMS:
	t = -2.0;
	break;

      default:
        MSQ_SETERR(err)("averaging method not available.",MsqError::NOT_IMPLEMENTED);
        break;
      }

      m = 0;
      for (i = 0; i < 8; ++i) {
	nm = pow(mMetrics[i], t);
	m += nm;

	g_factor[i] = 0.125*t*nm / mMetrics[i];
	h_factor[i] = (t-1)*g_factor[i] / mMetrics[i];
      }

      nm = 0.125 * m;

      l = 0;
      for (i = 0; i < 8; ++i) {
        g[locs_hex[i][0]] += g_factor[i]*mGradients[4*i+0];
        g[locs_hex[i][1]] += g_factor[i]*mGradients[4*i+1];
        g[locs_hex[i][2]] += g_factor[i]*mGradients[4*i+2];
        g[locs_hex[i][3]] += g_factor[i]*mGradients[4*i+3];

        for (j = 0; j < 4; ++j) {
          for (k = j; k < 4; ++k) {
	    outer = h_factor[i]*outer.outer_product(mGradients[4*i+j], 
						    mGradients[4*i+k]);

            r = locs_hex[i][j];
            c = locs_hex[i][k];

            if (r <= c) {
              loc = 8*r - (r*(r+1)/2) + c;
              h[loc] += g_factor[i]*mHessians[l] + outer;
            } 
            else {
              loc = 8*c - (c*(c+1)/2) + r;
              h[loc] += transpose(g_factor[i]*mHessians[l] + outer);
            }
            ++l;
          }
        }
      }

      m = pow(nm, 1.0 / t);
      g_factor[0] = m / (t*nm);
      h_factor[0] = (1.0 / t - 1)*g_factor[0] / nm;

      l = 0;
      for (i = 0; i < 8; ++i) {
	for (j = i; j < 8; ++j) {
	  outer = outer.outer_product(g[i], g[j]);
	  h[l] = g_factor[0]*h[l] + h_factor[0]*outer;
	  ++l;
	}
	g[i] *= g_factor[0];
      }
      break;
    }

    // zero out fixed elements of gradient and Hessian
    ind = 0;
    for (i=0; i<8; ++i) {
      // if free vertex, see next
      if (ind<nfv && vertices+v_i[i] == fv[ind] )
        ++ind;
      // else zero gradient entry and hessian entries.
      else {
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
    break;
  } // end switch over element type

  return true;
}

} // namespace Mesquite
