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
  \file   LPtoPTemplate.cpp
  \brief  

  This Objective Function is evaluated using an L P norm to the pth power.
  total=(sum (x_i)^pVal)
  \author Michael Brewer
  \author Thomas Leurent
  \date   2002-01-23
*/
#include <math.h>
#include "LPtoPTemplate.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqTimer.hpp"
#include "MsqHessian.hpp"
#include "MsqDebug.hpp"

using  namespace Mesquite;  

LPtoPTemplate::LPtoPTemplate(QualityMetric *qualitymetric, short Pinput, MsqError &err){
  set_quality_metric(qualitymetric);
  pVal=Pinput;
  if(pVal<1){
    MSQ_SETERR(err)("P_VALUE must be greater than 0.", MsqError::INVALID_ARG);
    return;
  }
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
    //set_use_local_gradient(true);
  set_negate_flag(qualitymetric->get_negate_flag());
  dividingByN=false;
}

//Michael:  need to clean up here
LPtoPTemplate::~LPtoPTemplate(){

}

bool LPtoPTemplate::concrete_evaluate(PatchData &pd, double &fval,
                                      MsqError &err){
  size_t index=0;
  MsqMeshEntity* elems=pd.get_element_array(err);
  bool obj_bool=true;
    //double check for pVal=0;
  if(pVal==0){
    MSQ_SETERR(err)("pVal equal zero not allowed.  L_0 is not a valid norm.",
                    MsqError::INVALID_STATE);
    return false;
  }
  
    //Michael:  this may not do what we want
    //Set currentQM to be the first quality metric* in the list
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM==NULL)
    currentQM=get_quality_metric_list().front();
  if(currentQM==NULL) {
    MSQ_SETERR(err)("NULL QualityMetric pointer in LPtoPTemplate",
                    MsqError::INVALID_STATE);
    return false;
  }
  size_t num_elements=pd.num_elements();
  size_t num_vertices=pd.num_vertices();
  size_t total_num=0;
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED)   
    total_num=num_elements;
  else if (currentQM->get_metric_type()==QualityMetric::VERTEX_BASED)
    total_num=num_vertices;
  else {
    MSQ_SETERR(err)("Make sure MetricType is initialised in concrete "
                    "QualityMetric constructor.", MsqError::INVALID_STATE);
    return false;
  }
  
  msq_std::vector<double> metric_values(total_num);
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED)
  {
    for (index=0; index<num_elements;index++)
    {
        //if invalid return false after clean-up
      obj_bool=currentQM->evaluate_element(pd, (&elems[index]),
                                           metric_values[index], err);
      if(MSQ_CHKERR(err) || !obj_bool){
        fval=0.0;
        return false;
      }
      
      metric_values[index]=fabs(metric_values[index]);
      MSQ_DBGOUT(3) << "      o  Quality metric value for element "
                    << index << "\t: " << metric_values[index] << "\n";
    }
  }
  else if(currentQM->get_metric_type()==QualityMetric::VERTEX_BASED)
  {
    MsqVertex* vertices=pd.get_vertex_array(err); MSQ_ERRZERO(err);
    for (index=0; index<num_vertices;index++)
    {
        //evaluate metric for this vertex
      obj_bool=currentQM->evaluate_vertex(pd, (&vertices[index]),
                                          metric_values[index], err);
        //if invalid return false after clean-up
      if(MSQ_CHKERR(err) || !obj_bool){
        fval=0.0;
        return false;
      }
      
      metric_values[index]=fabs(metric_values[index]);
    }
  }
  fval=compute_function(&metric_values[0], total_num, err);
  return !MSQ_CHKERR(err);
}

/* virtual function reimplemented from QualityMetric. No doxygen doc needed. */
bool LPtoPTemplate::compute_analytical_gradient(PatchData &pd,
                                              Vector3D *const &grad,
                                              double &OF_val,
                                              MsqError &err, size_t array_size)
{
  MSQ_FUNCTION_TIMER( "LPtoPTemplate::compute_analytical_gradient" );
  
    //initialize the scaling value
  double scaling_value=1.0;
 
  size_t num_elements=pd.num_elements();
  size_t num_vertices=pd.num_vertices();
  if( num_vertices!=array_size && array_size>0)
  {
    MSQ_SETERR(err)("Incorrect array size.", MsqError::INVALID_ARG);
    return false;
  }
  
  MsqMeshEntity* elems=pd.get_element_array(err); MSQ_ERRZERO(err);
  MsqVertex* vertices=pd.get_vertex_array(err); MSQ_ERRZERO(err);
  bool qm_bool=true;
  double QM_val;
  OF_val = 0.;
  size_t i;
  int p1;
  
  //Set currentQM to be quality metric (possibly composite) associated with the objective function
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM==NULL) {
    MSQ_SETERR(err)("LPtoPTemplate has NULL QualityMetric pointer.",MsqError::INVALID_STATE);
    return false;
  }
  enum QualityMetric::MetricType qm_type=currentQM->get_metric_type();
  
  if (qm_type!=QualityMetric::ELEMENT_BASED &&
      qm_type!=QualityMetric::VERTEX_BASED) {
    MSQ_SETERR(err)("Make sure MetricType is initialised"
                    "in concrete QualityMetric constructor.",
                    MsqError::INVALID_STATE);
    return false;
  }


  // zeros out objective function gradient
  for (i=0; i<num_vertices; ++i)
    grad[i] =0;
  
  // Computes objective function gradient for an element based metric
  if(qm_type==QualityMetric::ELEMENT_BASED){
      //if scaling, divid by num_elements
    if(dividingByN){
      if(num_elements<=0) {
        MSQ_SETERR(err)("The number of elements should not be zero.",MsqError::INVALID_MESH);
        return false;
      }
      
      scaling_value/=num_elements;
    }
    size_t e, ve;
    size_t nfve; // num free vtx in element
    size_t nve; // num vtx in element
    MsqVertex* ele_free_vtces[MSQ_MAX_NUM_VERT_PER_ENT];
    const size_t *ele_vtces_ind;

    // loops over all elements.
    for (e=0; e<num_elements; ++e) {
      // stores the pointers to the free vertices within the element
      // (using pointer arithmetic).
      nfve = 0;
      nve = elems[e].vertex_count();
      ele_vtces_ind = elems[e].get_vertex_index_array();
      for (ve=0; ve<nve; ++ve) {
        if (vertices[ele_vtces_ind[ve]].is_free_vertex()) {
          ele_free_vtces[nfve] = vertices + ele_vtces_ind[ve];
          ++nfve;
        }
      }

      // Computes q and grad(q)
      Vector3D grad_vec[MSQ_MAX_NUM_VERT_PER_ENT];
      qm_bool = currentQM->compute_element_gradient(
                                     pd, &elems[e],
                                     ele_free_vtces,
                                     grad_vec, nfve, QM_val, err);
      if(MSQ_CHKERR(err) || !qm_bool) return false;

      // computes p*|Q(e)|^{p-1}
      QM_val = fabs(QM_val);
      double QM_pow=1.0;
      double factor;
      if (pVal==1) factor=1;
      else {
        QM_pow=QM_val;
        for (p1=1; p1<pVal-1; ++p1)
          QM_pow*=QM_val;
        factor = QM_pow * pVal;
      }
        //this scales the gradient
      factor *= (scaling_value * get_negate_flag());

        // computes Objective Function value \sum_{i=1}^{N_e} |q_i|^P
        // possibly scaled by 1/num.
      OF_val += (scaling_value * QM_pow * QM_val);
      
      // For each free vertex in the element ... 
      for (i=0; i<nfve; ++i) {
        // ... computes p*q^{p-1}*grad(q) ...
        grad_vec[i] *= factor;
        // ... and accumulates it in the objective function gradient.
        grad[pd.get_vertex_index(ele_free_vtces[i])] += grad_vec[i];
      }
    }
  }
  
  // Computes objective function gradient for a vertex based metric  
  else if (qm_type==QualityMetric::VERTEX_BASED){
      //if scaling, divide by the number of vertices
    if(dividingByN){
      if(num_elements<=0) {
        MSQ_SETERR(err)("The number of vertices should not be zero.",MsqError::INVALID_MESH);
        return false;
      }
      
      scaling_value/=num_vertices;
    }
    //vector for storing indices of vertex's connected elems
    msq_std::vector<size_t> vert_on_vert_ind;
    //position in pd's vertex array
    size_t vert_count=0;
    //position in vertex array
    size_t vert_pos=0;
    //loop over the free vertex indices to find the gradient...
    size_t vfv_array_length=10;//holds the current legth of vert_free_vtces
    msq_std::vector<MsqVertex*> vert_free_vtces(vfv_array_length);
    msq_std::vector<Vector3D> grad_vec(vfv_array_length); 
    for(vert_count=0; vert_count<num_vertices; ++vert_count){
      //For now we compute the metric for attached vertices and this
      //vertex, the above line gives us the attached vertices.  Now,
      //we must add this vertex.
      pd.get_adjacent_vertex_indices(vert_count,
                                     vert_on_vert_ind,err);
      vert_on_vert_ind.push_back(vert_count);
      size_t vert_num_vtces = vert_on_vert_ind.size();
      
      // dynamic memory management if arrays are too small.
      if(vert_num_vtces > vfv_array_length){
        vfv_array_length=vert_num_vtces+5;
        vert_free_vtces.resize(vfv_array_length);
        grad_vec.resize(vfv_array_length);
      }
      
      size_t vert_num_free_vtces=0;
      //loop over the vertices connected to this one (vert_count)
      //and put the free ones into vert_free_vtces
      while(!vert_on_vert_ind.empty()){
        vert_pos=(vert_on_vert_ind.back());
        //clear the vector as we go
        vert_on_vert_ind.pop_back();
        //if the vertex is free, add it to ver_free_vtces
        if(vertices[vert_pos].is_free_vertex()){
          vert_free_vtces[vert_num_free_vtces]=&vertices[vert_pos];
          ++vert_num_free_vtces ;
        }
      }
      
      qm_bool=currentQM->compute_vertex_gradient(pd,
                                                 vertices[vert_count],
                                                 &vert_free_vtces[0],
                                                 &grad_vec[0],
                                                 vert_num_free_vtces,
                                                 QM_val, err);
      if(MSQ_CHKERR(err) || !qm_bool){
        return false;
      }
       // computes p*|Q(e)|^{p-1}
      QM_val = fabs(QM_val);
      double QM_pow = 1.0;
      double factor;
      if (pVal==1) factor=1;
      else {
        QM_pow=QM_val;
        for (p1=1; p1<pVal-1; ++p1)
          QM_pow*=QM_val;
        factor = QM_pow * pVal;
      }
        //this scales the gradient
      factor *= (scaling_value * get_negate_flag());

      // computes Objective Function value \sum_{i=1}^{N_v} |q_i|^P
        // possibly scaled by 1/num
      OF_val += (scaling_value * QM_pow * QM_val);
      // For each free vertex around the vertex (and the vertex itself if free) ... 
      for (i=0; i < vert_num_free_vtces ; ++i) {
        // ... computes p*q^{p-1}*grad(q) ...
        grad_vec[i] *= factor;
        // ... and accumulates it in the objective function gradient.
        grad[pd.get_vertex_index(vert_free_vtces[i])] += grad_vec[i];
      }
    }
  }

  OF_val *= get_negate_flag();
  
  return true;  

}
  
	
/*! \fn LPtoPTemplate::compute_analytical_hessian(PatchData &pd, MsqHessian &hessian, MsqError &err)

    For each element, each entry to be accumulated in the Hessian for
    this objective function (\f$ \sum_{e \in E} Q(e)^p \f$ where \f$ E \f$
    is the set of all elements in the patch) has the form:
    \f$ pQ(e)^{p-1} \nabla^2 Q(e) + p(p-1)Q(e)^{p-2} \nabla Q(e) [\nabla Q(e)]^T \f$.

    For \f$ p=2 \f$, this simplifies to
    \f$ 2Q(e) \nabla^2 Q(e) + 2 \nabla Q(e) [\nabla Q(e)]^T \f$.

    For \f$ p=1 \f$, this simplifies to \f$ \nabla^2 Q(e) \f$.

    The \f$ p=1 \f$ simplified version is implemented directly
    to speed up computation. 

    This function does not support vertex-based metrics.
    
    \param pd The PatchData object for which the objective function
           hessian is computed.
    \param hessian: this object must have been previously initialized.
*/
bool LPtoPTemplate::compute_analytical_hessian(PatchData &pd,
                                               MsqHessian &hessian,
                                               Vector3D *const &grad,
                                               double &OF_val,
                                               MsqError &err)
{
  double scaling_value=1.0;
  
  MSQ_FUNCTION_TIMER( "LPtoPTemplate::compute_analytical_hessian" );

  MsqMeshEntity* elements = pd.get_element_array(err); MSQ_ERRZERO(err);
  MsqVertex* vertices = pd.get_vertex_array(err); MSQ_ERRZERO(err);
  size_t num_elems = pd.num_elements();
    //if scaling divide by the number of elements.
  if(dividingByN){
    if(num_elems<=0) {
      MSQ_SETERR(err)("LPtoP is attempting to divide by zero in analytical Hessian.",
                      MsqError::INVALID_MESH);
      return false;
    }
    scaling_value/=num_elems;
  }
  
  size_t num_vertices = pd.num_vertices();
  Matrix3D elem_hessian[MSQ_MAX_NUM_VERT_PER_ENT*(MSQ_MAX_NUM_VERT_PER_ENT+1)/2];
  Matrix3D elem_outer_product;
  Vector3D grad_vec[MSQ_MAX_NUM_VERT_PER_ENT];
  double QM_val;
  double fac1, fac2;
  Matrix3D grad_outprod;
  bool qm_bool;
  QualityMetric* currentQM = get_quality_metric();
  
  MsqVertex* ele_free_vtces[MSQ_MAX_NUM_VERT_PER_ENT];
  short i;
  for (i=0; i<MSQ_MAX_NUM_VERT_PER_ENT; ++i) ele_free_vtces[i]=NULL;
  
  const size_t* vtx_indices;
    
  size_t e, v;
  size_t nfve; // number of free vertices in element
  short j,n;

  hessian.zero_out();
  for (v=0; v<num_vertices; ++v) grad[v] = 0.;
  OF_val = 0.;
  
  // Loops over all elements in the patch.
  for (e=0; e<num_elems; ++e) {
    short nve = elements[e].vertex_count();
    
    // Gets a list of free vertices in the element.
    vtx_indices = elements[e].get_vertex_index_array();
    nfve=0;
    for (i=0; i<nve; ++i) {
      if ( vertices[vtx_indices[i]].is_free_vertex() ) {
        ele_free_vtces[nfve] = vertices + vtx_indices[i];
        ++nfve;
      }
    }
    
    // Computes \nabla^2 Q(e). Only the free vertices will have non-zero entries. 
    qm_bool = currentQM->compute_element_hessian(pd,
                                    elements+e, ele_free_vtces,
                                    grad_vec, elem_hessian,
                                    nfve, QM_val, err);
    if (MSQ_CHKERR(err) || !qm_bool) return false;


    // **** Computes Hessian ****
    double QM_pow=1.;
    if (pVal == 1) {
      n=0;
      for (i=0; i<nve; ++i) {
        for (j=i; j<nve; ++j) {
            //negate if necessary
          elem_hessian[n] *= (scaling_value * get_negate_flag());
          ++n;
        }
      }
      hessian.accumulate_entries(pd, e, elem_hessian, err);
      fac1 = 1;
    }
    else if (pVal >= 2) {
      // Computes the coefficients:
      QM_val = fabs(QM_val);
      QM_pow = 1;
      for (i=0; i<pVal-2; ++i)
        QM_pow *= QM_val;
      // 1 - computes p(p-1)Q(e)^{p-2}
      fac2 = pVal* (pVal-1) * QM_pow;
      // 2 - computes  pQ(e)^{p-1}
      QM_pow *= QM_val;
      fac1 = pVal * QM_pow;

        //fac1 *= get_negate_flag();
        //fac2 *= get_negate_flag();

      n=0;
      for (i=0; i<nve; ++i) {
        for (j=i; j<nve; ++j) {
          if ( vertices[vtx_indices[i]].is_free_vertex() &&
               vertices[vtx_indices[j]].is_free_vertex() ) {
            // Computes \nabla Q(e) [\nabla Q(e)]^T 
            elem_outer_product.outer_product(grad_vec[i], grad_vec[j]);

	    elem_outer_product *= fac2;
	    elem_hessian[n] *= fac1;
	    elem_hessian[n] += elem_outer_product;
          } else {
            // elem_outer_product is nul 
            elem_hessian[n] *= fac1;
          }
            //scale the hessian by the scaling factor
          elem_hessian[n] *= (scaling_value * get_negate_flag());
          ++n;
        }
      }

      hessian.accumulate_entries(pd, e, elem_hessian, err); MSQ_ERRZERO(err);

    } else {
      MSQ_SETERR(err)(" invalid P value.", MsqError::INVALID_STATE);
      return false;
    }


    // **** Computes Gradient ****

    // For each vertex in the element ... 
    for (i=0; i<nve; ++i) {
      if ( vertices[vtx_indices[i]].is_free_vertex() ) {
        // ... computes p*q^{p-1}*grad(q) ...
        grad_vec[i] *= fac1*get_negate_flag();
        // ... and accumulates it in the objective function gradient.
          //also scale the gradient by the scaling factor
        grad[vtx_indices[i]] += (scaling_value * grad_vec[i]);
      }
    }
    
    // **** computes Objective Function value \sum_{i=1}^{N_e} |q_i|^P ****
      //and scale by 1/num if necessary
    OF_val += (scaling_value * QM_pow * QM_val);
    
  }

  OF_val *= get_negate_flag();
  
  return true;
}
