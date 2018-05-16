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
#ifndef __PN_PATCH_H__
#define __PN_PATCH_H__

#include "mopbasic.h"

MOP_BEGIN_NAMESPACE

// Only need at most three coords, b_n = 1-(sum_{i=1}^{n-1}b_i)
Vector_3<double> PN_project(std::vector<int *> v_ids,
			    Vector_3<double> bcoords,
			    const Vector_3<double> *pnts,
			    const Vector_3<double> *vnrms,
			    const Vector_3<double> *evects,
			    const Vector_3<double> *evals,
			    const Vector_3<double> *bs,
			    const int * tranks,
			    std::vector<bool> is_ridge[4]);

void project_edge(int v_id, int t_id, int id3,// Nodal ids
		  Vector_3<double> v_crd, // PN coordinates
		  Vector_3<double> t_crd, 
		  Vector_3<double> & p_crd,
		  Vector_3<double> & p2_crd,
		  const Vector_3<double> * pnts, // Nodal coords
		  const Vector_3<double> * vnrms,  // Normals
		  const Vector_3<double> * evects, // Eigenvectors
		  const Vector_3<double> * evals,
		  const Vector_3<double> * bs,
		  const int* tranks, // Tangent ranks
		  bool is_ridge); // Is the edge on a ridge?

Vector_3<double> one_sided_normal(int id1, int id2, int id3,
				  const Vector_3<double> *evects,
				  const Vector_3<double> *pnts,
				  const Vector_3<double> *evals,
				  const Vector_3<double> *vnorms,
				  const Vector_3<double> *bs);

MOP_END_NAMESPACE

#endif






