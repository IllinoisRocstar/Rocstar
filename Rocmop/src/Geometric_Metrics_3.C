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
// $Id: Geometric_Metrics_3.C,v 1.6 2008/12/06 08:45:24 mtcampbe Exp $

/*! \file Geometric_Metrics_3.C 
    \brief Implementation of 3D geometric metric quality measures.

    Implementation of geometric metrics:
    minimum and maximum dihedral angles,
    circumradius over inradius (R/r), and
    circumradius over shortest edge length (R/l)
    The latter two measures are defined only for simplicial elements.
*/

#include "geometry.h"
#include "Geometric_Metrics_3.h"
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

MOP_BEGIN_NAMESPACE

using namespace std;

void Geo_Metric_Base_3::initialize(Vector_3<double> n[], int type){
  _type = type;
  if (_type == COM::Connectivity::TET4){
    v.resize(7);
    v[0]= n[2]-n[1]; v[1]= n[1]-n[2]; v[2]= n[2]-n[0];
    v[3]= n[1]-n[0]; v[4]= n[3]-n[1]; v[5]= n[3]-n[0];
    v[6]= n[3]-n[2];
  }
  else if (_type == COM::Connectivity::HEX8){
    v.resize(8);
    v[0]= n[1]-n[0]; v[1]= n[3]-n[0]; v[2]= n[1]-n[2]; v[3]= n[3]-n[2]; 
    v[4]= n[5]-n[1]; v[5]= n[7]-n[3]; v[6]= n[5]-n[4]; v[7]= n[7]-n[4];
  }
}

void Geo_Metric_Base_3::initialize(Element_node_enumerator &ene){
  _type = ene.type();
  if(ene.type() == COM::Connectivity::TET4){
    v.resize(7);
    const std::string nc("nc");
    Element_node_vectors_k_const<double> n;
    n.set(ene.pane()->dataitem(nc),ene);
    v[0][0]=n(2,0)-n(1,0); v[0][1]=n(2,1)-n(1,1); v[0][2]=n(2,2)-n(1,2);
    v[1][0]=n(1,0)-n(2,0); v[1][1]=n(1,1)-n(2,1); v[1][2]=n(1,2)-n(2,2);    
    v[2][0]=n(2,0)-n(0,0); v[2][1]=n(2,1)-n(0,1); v[2][2]=n(2,2)-n(0,2);
    v[3][0]=n(1,0)-n(0,0); v[3][1]=n(1,1)-n(0,1); v[3][2]=n(1,2)-n(0,2);
    v[4][0]=n(3,0)-n(1,0); v[4][1]=n(3,1)-n(1,1); v[4][2]=n(3,2)-n(1,2);
    v[5][0]=n(3,0)-n(0,0); v[5][1]=n(3,1)-n(0,1); v[5][2]=n(3,2)-n(0,2);
    v[6][0]=n(3,0)-n(2,0); v[6][1]=n(3,1)-n(2,1); v[6][2]=n(3,2)-n(2,2);
  }
  else if (ene.type() == COM::Connectivity::HEX8){
    v.resize(8);
    const std::string nc("nc");
    Element_node_vectors_k_const<double> n;
    n.set(ene.pane()->dataitem(nc),ene);
    v[0][0]=n(1,0)-n(0,0); v[0][1]=n(1,1)-n(0,1); v[0][2]=n(1,2)-n(0,2);
    v[1][0]=n(3,0)-n(0,0); v[1][1]=n(3,1)-n(0,1); v[1][2]=n(3,2)-n(0,2);    
    v[2][0]=n(1,0)-n(2,0); v[2][1]=n(1,1)-n(2,1); v[2][2]=n(1,2)-n(2,2);
    v[3][0]=n(3,0)-n(2,0); v[3][1]=n(3,1)-n(2,1); v[3][2]=n(3,2)-n(2,2);
    v[4][0]=n(5,0)-n(1,0); v[4][1]=n(5,1)-n(1,1); v[4][2]=n(5,2)-n(1,2);
    v[5][0]=n(7,0)-n(3,0); v[5][1]=n(7,1)-n(3,1); v[5][2]=n(7,2)-n(3,2);
    v[6][0]=n(5,0)-n(4,0); v[6][1]=n(5,1)-n(4,1); v[6][2]=n(5,2)-n(4,2);
    v[7][0]=n(7,0)-n(4,0); v[7][1]=n(7,1)-n(4,1); v[7][2]=n(7,2)-n(4,2);
  }
  else if (ene.type() == COM::Connectivity::PYRIMID5 ||
	   ene.type() == COM::Connectivity::PRISM6){
    _type = ene.type();
    _ene_pane = ene.pane();
    _ene_i    = ene.id();
  }
  else COM_assertion_msg(0,"Element type not supported for 3D geometric metrics");
}

void Geo_Metric_Base_3::compute_angles(double& min, double& max) const{
  if (_type == COM::Connectivity::TET4){
    double a12,a13,a14,a23,a24,a34;
    Vector_3<double> n1,n2,n3,n4;
    unit_normal(v[2],v[3],n1); unit_normal(v[5],v[2],n2);
    unit_normal(v[3],v[5],n3); unit_normal(v[0],v[4],n4);
    a12=f_angle(n1,n2); a13=f_angle(n1,n3); a14=f_angle(n1,n4);
    a23=f_angle(n2,n3); a24=f_angle(n2,n4); a34=f_angle(n3,n4);
    min_max6(a12,a13,a14,a23,a24,a34,min,max);
  }
  else if (_type == COM::Connectivity::HEX8){
    double a12,a14,a25,a34,a26,a46,a13,a15,a23,a45,a36,a56;
    Vector_3<double> n1,n2,n3,n4,n5,n6;

    unit_normal(v[1],v[0],n1); unit_normal(v[0],v[4],n2); 
    unit_normal(v[4],v[2],n3); unit_normal(v[3],v[5],n4); 
    unit_normal(v[5],v[1],n5); unit_normal(v[6],v[7],n6);

    a12= f_angle(n1,n2); a13= f_angle(n1,n3); a14= f_angle(n1,n4);
    a15= f_angle(n1,n5); a25= f_angle(n2,n5); a23= f_angle(n2,n3);
    a34= f_angle(n3,n4); a45= f_angle(n4,n5); a26= f_angle(n2,n6); 
    a36= f_angle(n3,n6); a46= f_angle(n4,n6); a56= f_angle(n5,n6);

    min_max12(a12,a15,a34,a36,a13,a25,a45,a46,a14,a23,a26,a56,min,max);
  }
  else {

    const int face_node_lists_pyra[5][8] =
      { {0,3,2,1,8,7,6,5}, {0,1,4,5,10,9,-1,-1}, {1,2,4,6,11,10,-1,-1},
	{2,3,4,7,12,11,-1,-1}, {3,0,4,8,9,12,-1,-1}};
    
    const int face_node_lists_pris[5][9] =
      { {0,1,4,3,6,10,12,9,15}, {1,2,5,4,7,11,13,10,16}, {2,0,3,5,8,9,14,11,17},
	{0,2,1,8,7,6,-1,-1,-1}, {3,4,5,12,13,14,-1,-1,-1}};

    // Create an element node enumerator for handling this element
    Element_node_enumerator ene(_ene_pane, 
				_ene_i);
    // Loop over faces,
    //     a) determining which share an edge and
    //     b) face normals

    std::vector<Vector_3<double> > fnormals(ene.size_of_faces());
    Element_node_vectors_k_const<double> n;
    n.set(ene.pane()->dataitem("nc"),ene);
    // map from edges to faces
    std::map<std::pair<int,int>,int> e2f;
    std::map<std::pair<int,int>,int>::iterator e2f_it;
    // set of adjacent faces
    std::set<std::pair<int,int> > afs;
    std::set<std::pair<int,int> >::iterator afs_it;
    for(int i=0, ni = ene.size_of_faces(); i<ni; ++i){

          Facet_node_enumerator fne( &ene, i);

      const int *fn_list = (_type == COM::Connectivity::PYRIMID5) ? 
	face_node_lists_pyra[i] : face_node_lists_pris[i];

      int i0 = fn_list[0];
      int i1 = fn_list[1];
      int i2 = fn_list[2];

      // Find face normal
      Vector_3<double> e1,e2;
      e1[0] = n(i1,0)-n(i0,0);
      e1[1] = n(i1,1)-n(i0,1);
      e1[2] = n(i1,2)-n(i0,2);
      e2[0] = n(i2,0)-n(i1,0);
      e2[1] = n(i2,1)-n(i1,1);
      e2[2] = n(i2,2)-n(i1,2);
      unit_normal(e1,e2,fnormals[i]);

      for(int j=0, nj = fne.size_of_edges();j<nj;++j){

	// Loop over the edges of this face.
	// If we there is a matching edge, store the adjacent
	// faces.
	int node_id1 = fne[j],node_id2 = fne[(j+1)%nj];
	if(node_id1>node_id2)
	  std::swap(node_id1, node_id2);
	e2f_it = e2f.find(std::make_pair<int,int>(node_id1,node_id2));
	if(e2f_it != e2f.end()){
	  afs.insert(std::make_pair<int,int>(e2f_it->second,i));
	  e2f.erase(e2f_it);
	}
	else{
	  e2f.insert(std::make_pair<std::pair<int,int>,int>(
		     std::make_pair<int,int>(node_id1,node_id2)
		     ,i));
	}	
      }
      min = 180.0; max = 0.0;
      for(afs_it = afs.begin(); afs_it != afs.end(); ++afs_it){
	int ind1 = afs_it->first;
	int ind2 = afs_it->second;
	double angle = f_angle(fnormals[ind1],
			       fnormals[ind2]);
	min = (angle < min) ? angle : min;
	max = (angle > max) ? angle : max;
      }
    }
  }
}

void Geo_Metric_Base_3::compute_aspects(double& R, double& r, double& l) const {

  // This metric is only defined for tetrahedrons
  if (_type == COM::Connectivity::TET4){

    // Use formulas from mathworld to find circumradius (R),
    // inradius (r), and shortest edge length (l)

    double a1 = tri_area(v[0],v[2],v[3]);
    double a2 = tri_area(v[0],v[4],v[6]);
    double a3 = tri_area(v[2],v[5],v[6]);
    double a4 = tri_area(v[3],v[4],v[5]);
    double a_sum = a1+a2+a3+a4;


    double len_a = edge_length(v[0]);
    double len_b = edge_length(v[2]);
    double len_c = edge_length(v[3]);
    double len_ap = edge_length(v[5]);
    double len_bp = edge_length(v[4]);
    double len_cp = edge_length(v[6]);

    double p1 = len_a * len_ap;
    double p2 = len_b * len_bp;
    double p3 = len_c * len_cp;
    double p_prod = sqrt((p1+p2+p3)*(p1+p2-p3)*(p1-p2+p3)*(p2+p3-p1));

    double vol = tet_vol(v[5],v[4],v[6]);

    R= p_prod/(24.0*vol);
    r = (3.0*vol)/a_sum;
    l = len_a;
    if( len_b < l ){ l = len_b;}
    if( len_c < l ){ l = len_c;}
    if( len_ap < l ){ l = len_ap;}
    if( len_bp < l ){ l = len_bp;}
    if( len_cp < l ){ l = len_cp;}
  }
  else{
    R = 1.0;
    r = 1.0;
    l = 1.0;
  }
}

void Angle_Metric_3::compute(double atts[]) const {
    compute_angles(atts[0],atts[1]);
}

double Angle_Metric_3::maxValue() const { 
return 180.0; 
}

double Angle_Metric_3::minValue() const { 
return 0.0; 
}

void Aspect_Metric_3::compute(double atts[]) const {
  if (_type == COM::Connectivity::HEX8) {
    atts[0] = -1;
    atts[1] = -1;
  }
  else{
    double R, r, l;
    compute_aspects(R,r,l);
    atts[0] = (r/R)*3.0;
    atts[1] = (l/R)*.61237243569579447;
    if(atts[0] > 1.0)
      atts[0] = 1.0;
    if(atts[1] > 1.0)
      atts[0] = 1.0;    
  }
}

double Aspect_Metric_3::maxValue() const { 
  if (_type == COM::Connectivity::TET4) {
    return 1.0; 
  }
  else { return -1.0; }
}

double Aspect_Metric_3::minValue() const { 
  if (_type == COM::Connectivity::TET4) {
    return 0.0;
  }
  else { return -1.0; }
}

MOP_END_NAMESPACE






