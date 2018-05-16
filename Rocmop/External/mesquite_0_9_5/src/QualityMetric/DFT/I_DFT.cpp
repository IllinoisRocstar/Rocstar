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
  \file   I_DFT.cpp
  \brief  

  \author Thomas Leurent

  \date   2004-04-12
*/

#include "I_DFT.hpp"
#include "I_DFTFamilyFunctions.hpp"
#include "TargetMatrix.hpp"

using namespace Mesquite;
   

bool I_DFT::evaluate_element(PatchData& pd,
			     MsqMeshEntity* e,
			     double& m, 
			     MsqError &err)
{
  // Only works with the weighted average

  MsqError  mErr;
  MsqVertex *vertices = pd.get_vertex_array(err); MSQ_ERRZERO(err);

  EntityTopology topo = e->get_element_type();

  const size_t nv = e->vertex_count();
  const size_t *v_i = e->get_vertex_index_array();

  size_t idx = pd.get_element_index(e);
  const TargetMatrix *W = pd.targetMatrices.get_element_corner_tags(&pd, idx, err );
  MSQ_ERRZERO(err);

  // Initialize constants for the metric
  const double delta = useBarrierDelta ? pd.get_barrier_delta(err) :
    (mGamma ? 0 : 1);
  MSQ_ERRZERO(err);

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

    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }
      
      mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_2(mMetric, mCoords, mNormals[i], invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;
      
    }

    m *= MSQ_ONE_THIRD;
    break;

  case QUADRILATERAL:
    assert(4 == nv);

    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 3; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_2(mMetric, mCoords, mNormals[i], invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;
    }

    m *= 0.25;
    break;

  case TETRAHEDRON:
    assert(4 == nv);

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[tetInd[i][j]]];
      }

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;
    }

    m *= 0.25;
    break;

  case HEXAHEDRON:
    assert(8 == nv);

    for (i = 0; i < 8; ++i) {
      for (j = 0; j < 4; ++j) {
	mCoords[j] = vertices[v_i[hexInd[i][j]]];
      }

      QR(mQ, mR, W[i]);
      inv(invR, mR);
      mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			mAlpha, mGamma, delta, mBeta);
      
      if (!mValid) return false;
      m += W[i].get_cK() * mMetric;
    }

    m *= 0.125;
    break;

  default:
    MSQ_SETERR(err)("element type not implemented.",MsqError::NOT_IMPLEMENTED);
    return false;
  }

  return true;
}

bool I_DFT::compute_element_analytical_gradient(PatchData &pd,
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
  const double delta = useBarrierDelta ? pd.get_barrier_delta(err) :
    (mGamma ? 0 : 1);
    //const double delta = useBarrierDelta ? pd.get_barrier_delta(err) : 0;
  MSQ_ERRZERO(err);

  const int tetInd[4][4] = {{0, 1, 2, 3}, {1, 2, 0, 3},
			    {2, 0, 1, 3}, {3, 2, 1, 0}};
  const int hexInd[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
			    {2, 3, 1, 6}, {3, 0, 2, 7},
			    {4, 7, 5, 0}, {5, 4, 6, 1},
			    {6, 5, 7, 2}, {7, 6, 4, 3}};

  // Variables used for computing the metric
  double   mMetric;		// Metric value
  bool     mValid;		// Validity of the metric
  int      i, j, mVert;

  m = 0.0;
  switch(topo) {
  case TRIANGLE:
    assert(3 == nv);
    
    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constant, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.
    
    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 3; ++i) {
	mVert = -1;
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	  if (vertices + v_i[tetInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}

	if (mVert >= 0) {
	  mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_2_v0(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_2_v1(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_2_v2(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	  g[0] += W[i].get_cK() * mGrads[0];
	}
	else {
	  // For triangles, the free vertex must appear in every element.
	  // Therefore, there these accumulations should not get used.

	  mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  mValid = m_gdft_2(mMetric, mCoords, mNormals[i],
			    invR, mQ, mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	}
      }

      m *= MSQ_ONE_THIRD;
      g[0] *= MSQ_ONE_THIRD;
    }
    else {
      for (i = 0; i < 3; ++i) {
	mAccGrads[i] = 0.0;
      }
    
      for (i = 0; i < 3; ++i) {
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	}
      
	mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
      
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_2(mMetric, mGrads, mCoords, mNormals[i], invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);

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

      // This is not very efficient, but is one way to select correct 
      // gradients.  For gradients, info is returned only for free 
      // vertices, in the order of fv[].

      for (i = 0; i < 3; ++i) {
	for (j = 0; j < nfv; ++j) {
	  if (vertices + v_i[i] == fv[j]) {
	    g[j] = mAccGrads[i];
	  }
	}
      }
    }
    break;

  case QUADRILATERAL:
    assert(4 == nv);
    
    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constant, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	  if (vertices + v_i[hexInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}

	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_2_v0(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_2_v1(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_2_v2(mMetric, mGrads[0], mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	  g[0] += W[i].get_cK() * mGrads[0];
	}
	else {
	  // For quadrilaterals, the free vertex only appears in three 
	  // elements.  Therefore, there these accumulations are needed 
	  // to get the true local objective function.  Note: this code 
          // can be commented out for local codes to improve performance 
          // because you are unable to change the contributions from the 
	  // elements where the free vertex does not appear.  (If the 
	  // weight matrices change, then the code needs to be modified.)
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  mValid = m_gdft_2(mMetric, mCoords, mNormals[i],
			    invR, mQ, mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	}
      }

      m *= 0.25;
      g[0] *= 0.25;
    }
    else {
      for (i = 0; i < 4; ++i) {
	mAccGrads[i] = 0.0;
      }

      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	}

	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_2(mMetric, mGrads, mCoords, mNormals[i], invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);

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
    }
    break;

  case TETRAHEDRON:
    assert(4 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	  if (vertices + v_i[tetInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_3_v0(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_3_v1(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = g_gdft_3_v2(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_3_v3(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	  g[0] += W[i].get_cK() * mGrads[0];
	}
	else {
	  // For tetrahedrons, the free vertex must appear in every element.
	  // Therefore, there these accumulations should not get used.

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	}
      }

      m *= 0.25;
      g[0] *= 0.25;
    }
    else {
      for (i = 0; i < 4; ++i) {
	mAccGrads[i] = 0.0;
      }
      
      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_3(mMetric, mGrads, mCoords, invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);
	
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
      
      // This is not very efficient, but is one way to select correct 
      // gradients.  For gradients, info is returned only for free vertices, 
      // in the order of fv[].

      for (i = 0; i < 4; ++i) {
	for (j = 0; j < nfv; ++j) {
	  if (vertices + v_i[i] == fv[j]) {
	    g[j] = mAccGrads[i];
	  }
	}
      }
    }
    break;

  case HEXAHEDRON:
    assert(8 == nv);

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the gradient.

      g[0] = 0.0;
      for (i = 0; i < 8; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	  if (vertices + v_i[hexInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}

	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = g_gdft_3_v0(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = g_gdft_3_v1(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = g_gdft_3_v2(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = g_gdft_3_v3(mMetric, mGrads[0], mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }
	  
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	  g[0] += W[i].get_cK() * mGrads[0];
	}
	else {
	  // For hexahedrons, the free vertex only appears in four elements.
	  // Therefore, there these accumulations are needed to get the
	  // true local objective function.  Note: this code can be commented 
	  // out for local codes to improve performance because you are 
	  // unable to change the contributions from the elements where the 
	  // free vertex does not appear.  (If the weight matrices change, 
	  // then the code needs to be modified.)

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	}
      }
      
      m *= 0.125;
      g[0] *= 0.125;
    }
    else {
      for (i = 0; i < 8; ++i) {
	mAccGrads[i] = 0.0;
      }
      
      for (i = 0; i < 8; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = g_gdft_3(mMetric, mGrads, mCoords, invR, mQ, 
			  mAlpha, mGamma, delta, mBeta);
	
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
      
      // This is not very efficient, but is one way to select correct 
      // gradients.  For gradients, info is returned only for free 
      // vertices, in the order of fv[].
      
      for (i = 0; i < 8; ++i) {
	for (j = 0; j < nfv; ++j) {
	  if (vertices + v_i[i] == fv[j]) {
	    g[j] = mAccGrads[i];
	  }
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

bool I_DFT::compute_element_analytical_hessian(PatchData &pd,
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
  const double delta = useBarrierDelta ? pd.get_barrier_delta(err) :
    (mGamma ? 0 : 1);  
    //const double delta = useBarrierDelta ? pd.get_barrier_delta(err) : 0;
  MSQ_ERRZERO(err);

  const int tetInd[4][4] = {{0, 1, 2, 3}, {1, 2, 0, 3},
			    {2, 0, 1, 3}, {3, 2, 1, 0}};
  const int hexInd[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
			    {2, 3, 1, 6}, {3, 0, 2, 7},
			    {4, 7, 5, 0}, {5, 4, 6, 1},
			    {6, 5, 7, 2}, {7, 6, 4, 3}};

  // Variables used for computing the metric
  double   mMetric;		// Metric value
  bool     mValid;		// Validity of the metric
  int      i, j, k, l, mVert;
  int      row, col, loc;

  m = 0.0;
  switch(topo) {
  case TRIANGLE:
    assert(3 == nv);
    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constant, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    // Zero out the hessian and gradient vector
    for (i = 0; i < 3; ++i) {
      g[i] = 0.0;
    }

    for (i = 0; i < 6; ++i) {
      h[i].zero();
    }

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 3; ++i) {
	mVert = -1;
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	  if (vertices + v_i[tetInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_2_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_2_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_2_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	  mG += W[i].get_cK() * mGrads[0];
	  mH += W[i].get_cK() * mHessians[0];
	}
	else {
	  // For triangles, the free vertex must appear in every element.
	  // Therefore, there these accumulations should not get used.

	  mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  mValid = m_gdft_2(mMetric, mCoords, mNormals[i],
			    invR, mQ, mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	}
      }

      m *= MSQ_ONE_THIRD;
      mG *= MSQ_ONE_THIRD;
      mH *= MSQ_ONE_THIRD;

      for (i = 0; i < 3; ++i) {
	if (vertices + v_i[i] == fv[0]) {
	  // free vertex, see next
	  g[i] = mG;
	  switch(i) {
	  case 0:
	    h[0] = mH;
	    break;

	  case 1:
	    h[3] = mH;
	    break;

	  default:
	    h[5] = mH;
	    break;
	  }
	  break;
	}
      }
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 3; ++i) {
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	}

	mNormals[i] *= MSQ_3RT_2_OVER_6RT_3;
      
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_2(mMetric, mGrads, mHessians, mCoords, mNormals[i], 
			  invR, mQ, mAlpha, mGamma, delta, mBeta);

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
    }
    break;

  case QUADRILATERAL:
    assert(4 == nv);
    e->compute_corner_normals( mNormals, pd, err ); MSQ_ERRZERO(err);

    // The following analytic calculation only works correctly if the
    // normal is constant.  If the normal is not constant, you need
    // to get the gradient of the normal with respect to the vertex
    // positions to obtain the correct values.

    // Zero out the hessian and gradient vector
    for (i = 0; i < 4; ++i) {
      g[i] = 0.0;
    }

    for (i = 0; i < 10; ++i) {
      h[i].zero();
    }

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	  if (vertices + v_i[hexInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_2_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_2_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_2_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, mNormals[i],
				 invR, mQ, mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	  mG += W[i].get_cK() * mGrads[0];
	  mH += W[i].get_cK() * mHessians[0];
	}
	else {
	  // For quadrilaterals, the free vertex only appears in three 
	  // elements.  Therefore, there these accumulations are needed 
	  // to get the true local objective function.  Note: this code 
          // can be commented out for local codes to improve performance 
          // because you are unable to change the contributions from the 
	  // elements where the free vertex does not appear.  (If the 
	  // weight matrices change, then the code needs to be modified.)
	  
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);

	  mValid = m_gdft_2(mMetric, mCoords, mNormals[i],
			    invR, mQ, mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	}
      }

      m *= 0.25;
      mG *= 0.25;
      mH *= 0.25;

      for (i = 0; i < 4; ++i) {
	if (vertices + v_i[i] == fv[0]) {
	  // free vertex, see next
	  g[i] = mG;
	  switch(i) {
	  case 0:
	    h[0] = mH;
	    break;

	  case 1:
	    h[4] = mH;
	    break;

	  case 2:
	    h[7] = mH;
	    break;

	  default:
	    h[9] = mH;
	    break;
	  }
	  break;
	}
      }
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 3; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	}

	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_2(mMetric, mGrads, mHessians, mCoords, mNormals[i], 
			  invR, mQ, mAlpha, mGamma, delta, mBeta);

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
    }
    break;

  case TETRAHEDRON:
    assert(4 == nv);

    // Zero out the hessian and gradient vector
    for (i = 0; i < 4; ++i) {
      g[i] = 0.0;
    }

    for (i = 0; i < 10; ++i) {
      h[i].zero();
    }

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 4; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	  if (vertices + v_i[tetInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_3_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_3_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = h_gdft_3_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_3_v3(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	  mG += W[i].get_cK() * mGrads[0];
	  mH += W[i].get_cK() * mHessians[0];
	}
	else {
	  // For tetrahedrons, the free vertex must appear in every element.
	  // Therefore, there these accumulations should not get used.

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	}
      }

      m *= 0.25;
      mG *= 0.25;
      mH *= 0.25;

      for (i = 0; i < 4; ++i) {
	if (vertices + v_i[i] == fv[0]) {
	  // free vertex, see next
	  g[i] = mG;
	  switch(i) {
	  case 0:
	    h[0] = mH;
	    break;

	  case 1:
	    h[4] = mH;
	    break;

	  case 2:
	    h[7] = mH;
	    break;

	  default:
	    h[9] = mH;
	    break;
	  }
	  break;
	}
      }
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 4; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[tetInd[i][j]]];
	}
	
	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_3(mMetric, mGrads, mHessians, mCoords, invR, mQ,
			  mAlpha, mGamma, delta, mBeta);
	
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
    }
    break;

  case HEXAHEDRON:
    assert(8 == nv);

    // Zero out the hessian and gradient vector
    for (i = 0; i < 8; ++i) {
      g[i] = 0.0;
    }

    for (i = 0; i < 36; ++i) {
      h[i].zero();
    }

    if (1 == nfv) {
      // One free vertex; use the specialized code for computing the 
      // gradient and Hessian.

      Vector3D mG;
      Matrix3D mH;

      mG = 0.0;
      mH.zero();

      for (i = 0; i < 8; ++i) {
	mVert = -1;
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	  if (vertices + v_i[hexInd[i][j]] == fv[0]) {
	    mVert = j;
	  }
	}
	
	if (mVert >= 0) {
	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  
	  switch(mVert) {
	  case 0:
	    mValid = h_gdft_3_v0(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 1:
	    mValid = h_gdft_3_v1(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  case 2:
	    mValid = h_gdft_3_v2(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	    
	  default:
	    mValid = h_gdft_3_v3(mMetric, mGrads[0], mHessians[0],
				 mCoords, invR, mQ, 
				 mAlpha, mGamma, delta, mBeta);
	    break;
	  }

	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	  mG += W[i].get_cK() * mGrads[0];
	  mH += W[i].get_cK() * mHessians[0];
	}
	else {
	  // For hexahedrons, the free vertex only appears in four elements.
	  // Therefore, there these accumulations are needed to get the
	  // true local objective function.  Note: this code can be commented 
	  // out for local codes to improve performance because you are 
	  // unable to change the contributions from the elements where the 
	  // free vertex does not appear.  (If the weight matrices change, 
	  // then the code needs to be modified.)

	  QR(mQ, mR, W[i]);
	  inv(invR, mR);
	  mValid = m_gdft_3(mMetric, mCoords, invR, mQ, 
			    mAlpha, mGamma, delta, mBeta);
	  if (!mValid) return false;
	  m += W[i].get_cK() * mMetric;
	}
      }

      m *= 0.125;
      mG *= 0.125;
      mH *= 0.125;

      for (i = 0; i < 8; ++i) {
	if (vertices + v_i[i] == fv[0]) {
	  // free vertex, see next
	  g[i] = mG;
	  switch(i) {
	  case 0:
	    h[0] = mH;
	    break;

	  case 1:
	    h[8] = mH;
	    break;

	  case 2:
	    h[15] = mH;
	    break;

	  case 3:
	    h[21] = mH;
	    break;

	  case 4:
	    h[26] = mH;
	    break;

	  case 5:
	    h[30] = mH;
	    break;

	  case 6:
	    h[33] = mH;
	    break;

	  default:
	    h[35] = mH;
	    break;
	  }
	  break;
	}
      }
    }
    else {
      // Compute the metric and sum them together
      for (i = 0; i < 8; ++i) {
	for (j = 0; j < 4; ++j) {
	  mCoords[j] = vertices[v_i[hexInd[i][j]]];
	}

	QR(mQ, mR, W[i]);
	inv(invR, mR);
	mValid = h_gdft_3(mMetric, mGrads, mHessians, mCoords, invR, mQ,
			  mAlpha, mGamma, delta, mBeta);
      
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
    }
    break;

  default:
    MSQ_SETERR(err)("element type not implemented.",MsqError::NOT_IMPLEMENTED);
    return false;
  }
  return true;
}
