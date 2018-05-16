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
// $Id: mopbasic.h,v 1.6 2008/12/06 08:45:24 mtcampbe Exp $

#ifndef __MOP_BASIC_H_
#define __MOP_BASIC_H_

#define MOP_BEGIN_NAMESPACE  namespace MOP {
#define MOP_END_NAMESPACE    }
#define USE_MOP_NAMESPACE    using namespace MOP;

#include "Manifold_2.h"
#include "Element_accessors.hpp"

#include <cstdlib>

MOP_BEGIN_NAMESPACE

using COM::Element_node_enumerator;
using COM::Element_node_enumerator_str_2;
using COM::Element_node_enumerator_uns;
using COM::Facet_node_enumerator;
using COM::Element_vectors_k_const;
using COM::Element_vectors_k;
using COM::Element_node_vectors_k_const;
using COM::Element_node_vectors_k;

using SURF::Origin;
using SURF::Null_vector;
using SURF::Vector_3;
using SURF::Point_3;
using SURF::Vector_2;
using SURF::Point_2;

using SURF::Window_manifold_2;
using SURF::Node;
using SURF::Halfedge;
using SURF::Edge_ID;

typedef double Real;

MOP_END_NAMESPACE

#endif






