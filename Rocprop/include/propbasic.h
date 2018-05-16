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
// $Id: propbasic.h,v 1.6 2008/12/06 08:45:27 mtcampbe Exp $

#ifndef __PROP_BASIC_H
#define __PROP_BASIC_H

#define PROP_BEGIN_NAMESPACE  namespace PROP {
#define PROP_END_NAMESPACE    }
#define USE_PROP_NAMESPACE    using namespace PROP;

#include "Manifold_2.h"

PROP_BEGIN_NAMESPACE

using SURF::Node;
using SURF::Halfedge;
using SURF::Edge_ID;

using SURF::get_normal;
using SURF::get_tangent;

enum { SMOOTHER_NONE, SMOOTHER_ANISOTROPIC, SMOOTHER_LAPLACIAN,
       SMOOTHER_ANGLE_WEIGHTED_CENTROID};

typedef SURF::Window_manifold_2 Manifold;
typedef SURF::Point_3<double>   Point_3;
typedef SURF::Point_2<double>   Point_2;
typedef SURF::Vector_3<double>  Vector_3;
typedef SURF::Vector_2<double>  Vector_2;

PROP_END_NAMESPACE

#endif






