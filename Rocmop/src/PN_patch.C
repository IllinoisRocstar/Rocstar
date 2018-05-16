/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/
#include "PN_patch.h"
#include <cmath>
using namespace std;

MOP_BEGIN_NAMESPACE


// Only need at most three coords, b_n = 1-(sum_{i=1}^{n-1}b_i)
Vector_3<double> PN_project(std::vector<int> v_ids,
			    Vector_3<double> bcoords,
			    const Vector_3<double> *pnts,
			    const Vector_3<double> *vnrms,
			    const Vector_3<double> *evects,
			    const Vector_3<double> *evals,
			    const Vector_3<double> *bs,
			    const int * tranks,
			    bool is_ridge[4]){
  Vector_3<double> ans(0.0,0.0,0.0);
  if(v_ids.size() == 3){
    
    Vector_3<double> cp[10]; // PN triangle control points
    cp[0] = pnts[v_ids[0]-1];
    cp[1] = pnts[v_ids[1]-1];
    cp[2] = pnts[v_ids[2]-1];
    cp[3] = (2.0*cp[0]+cp[2])/3.0;
    cp[4] = (2.0*cp[0]+cp[1])/3.0;
    cp[5] = (2.0*cp[1]+cp[0])/3.0;
    cp[6] = (2.0*cp[1]+cp[2])/3.0;
    cp[7] = (2.0*cp[2]+cp[1])/3.0;
    cp[8] = (2.0*cp[2]+cp[0])/3.0;
    
    int r_ind[3] = {4,6,8};
    int l_ind[3] = {3,5,7};
    
    for(int k =0; k < 3; ++k){
      int v_ind = k, t_ind = (4+k)%3;
      int v_id = v_ids[v_ind], t_id = v_ids[t_ind];
      project_edge(v_id,t_id, v_ids[(5+k)%3],// Nodal ids
		   cp[v_ind], // PN coordinates
		   cp[t_ind], 
		   cp[r_ind[v_id]],
		   cp[l_ind[t_id]],
		   pnts, // Nodal coords
		   vnrms,  // Normals
		   evects, // Eigenvectors
		   evals,
		   bs,
		   tranks, // Tangent ranks
		   is_ridge[k]); // Is the edge on a ridge?
      
      project_edge(t_id,v_id,v_ids[(5+k)%3], // Nodal ids
		   cp[t_ind], // PN coordinates
		   cp[v_ind], 
		   cp[l_ind[t_id]],
		   cp[r_ind[v_id]],
		   pnts, // Nodal coords
		   vnrms,  // Normals
		   evects, // Eigenvectors
		   evals,
		   bs,
		   tranks, // Tangent ranks
		   is_ridge[k]);  // Is the edge on a ridge?
    }
    Vector_3<double> E = (cp[3] + cp[4] + cp[5] +
			  cp[6]+cp[7]+cp[8])/6.0;
    Vector_3<double> V = (cp[0] + cp[1] + cp[2])/3.0;
    cp[9] = E+(E-V)/2.0;
    double w = bcoords[0], u = bcoords[1];
    double v = 1.0-(w+u);
    double w2 = w*w, u2 = u*u, v2=v*v;
    
    ans =  ( cp[0]*w2*w + cp[1]*u2*u + cp[2]*v2*v +
	     3.0*(cp[3]*w2*v + cp[4]*w2*u + cp[5]*u2*w +
		  cp[6]*u2*v + cp[7]*v2*u + cp[8]*v2*w) + 
	     6.0*cp[9]*w*u*v);    
  }  
  return ans;
}


void project_edge(int v_id, int t_id, int id3,// Nodal ids
		  Vector_3<double> v_crd, // PN coordinates
		  Vector_3<double> t_crd, 
		  Vector_3<double> & p_crd,
		  Vector_3<double> & p2_crd,
		  const Vector_3<double> * pnts, // Nodal coords
		  const Vector_3<double> * vnrms,  // Normals
		  const Vector_3<double> * evects, // Eigenvectors
		  const Vector_3<double> *evals,
		  const Vector_3<double> * bs,
		  const int* tranks, // Tangent ranks
		  bool is_ridge) // Is the edge on a ridge?
{
  // convert from ids to indices
  v_id -=1;
  t_id -=1;
  Vector_3<double> vp_vect = (v_crd-p_crd);
  Vector_3<double> f(0.0,0.0,0.0);
  switch(tranks[v_id]) {
  case 2:{ // V is smooth    
    f = (vp_vect)*evects[3*v_id]*evects[3*v_id];
    break;
  }
  case 1:{ // V is a ridge
    if(is_ridge) // edge VT is on the ridge
      f = vp_vect - vp_vect*evects[3*v_id+2]*evects[3*v_id+2];
    else{ // edge VT is not part of the ridge
      Vector_3<double> n_os = 
	one_sided_normal(v_id+1,t_id+1,id3,
			 evects,
			 pnts,
			 evals,
			 vnrms,
			 bs);
      // one_sided_normal(v_id,n_os);
      f = vp_vect*n_os*n_os;
    }
    break;
  }
  case 0:{ // V is a corner, do nothing.
    return;
  }
  default :
    COM_assertion_msg(0, "Invalid tangent space size");
  }
  if(tranks[t_id]==0){// Other point is a corner, handle it here.
    Vector_3<double> tangent = p_crd - v_crd;
    tangent.normalize();
    p2_crd += f - 2.0*f*tangent*tangent;
  }
  p_crd += f;
}

Vector_3<double> one_sided_normal(int id1, int id2, int id3,
				  const Vector_3<double> *evects,
				  const Vector_3<double> *pnts,
				  const Vector_3<double> *evals,
				  const Vector_3<double> *vnorms,
				  const Vector_3<double> *bs){
  int i = ((evects[3*id1-3]*bs[id1-1]/-evals[id1-1][0]) >
    (evects[3*id1-2]*bs[id1-1]/-evals[id1-1][1])) ? 0 : 1;
  int j = 1-i;
  Vector_3<double> os = std::sqrt(evals[id1-1][i]) * vnorms[id1-1];
  Vector_3<double> fnormal = Vector_3<double>::cross_product(pnts[id2-1]-pnts[id1-1],
							     pnts[id3-1]-pnts[id1-1]);
  Vector_3<double> y = Vector_3<double>::cross_product(vnorms[id1-1],evects[3*id1-1]);
  if(fnormal*y >= 0)
    os += std::sqrt(evals[id1-1][j])*y;
  else
    os -= std::sqrt(evals[id1-1][j])*y;
  os.normalize();
  return os;
}

MOP_END_NAMESPACE






