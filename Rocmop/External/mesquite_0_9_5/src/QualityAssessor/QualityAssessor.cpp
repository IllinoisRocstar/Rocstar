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
  \file   QualityAssessor.cpp
  \brief  Member function of the Mesquite::QualityAssessor class

  \author Thomas Leurent
  \date   2002-05-23
*/

#include "QualityAssessor.hpp"
#include "PatchData.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqVertex.hpp"
#include "MsqDebug.hpp"
#include "MeshSet.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <vector.h>
#else
#  include <list>
#  include <vector>
   using namespace std;
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <iostream.h>
#  include <iomanip.h>
#else
#  include <iostream>
#  include <iomanip>
   using namespace std;
#endif

namespace Mesquite {

const int DEFAULT_HISTOGRAM_INTERVALS = 10;

QualityAssessor::QualityAssessor(string name) :
  qualityAssessorName(name),
  outputStream( cout ),
  printSummary( true ),
  stoppingMetric( assessList.end() ),
  stoppingFunction( NO_FUNCTION )
{ 
  MsqError err;
  set_patch_type( PatchData::GLOBAL_PATCH, err, 0 );
}

QualityAssessor::QualityAssessor(ostream& stream, string name) :
  qualityAssessorName(name),
  outputStream( stream ),
  printSummary( true ),
  stoppingMetric( assessList.end() ),
  stoppingFunction( NO_FUNCTION )
{ 
  MsqError err;
  set_patch_type( PatchData::GLOBAL_PATCH, err, 0 );
}

QualityAssessor::QualityAssessor( QualityMetric* metric,
                                  QAFunction function,
                                  MsqError& err,
                                  string name ) :
  qualityAssessorName(name),
  outputStream( cout ),
  printSummary( true ),
  stoppingMetric( assessList.end() ),
  stoppingFunction( (QAFunction)0 )
{ 
  set_patch_type( PatchData::GLOBAL_PATCH, err, 0 );
  add_quality_assessment( metric, function, err );
  set_stopping_assessment( metric, function, err );
}

QualityAssessor::QualityAssessor( QualityMetric* metric,
                                  QAFunction function,
                                  ostream& stream, 
                                  MsqError& err,
                                  string name ) :
  qualityAssessorName(name),
  outputStream( stream ),
  printSummary( true ),
  stoppingMetric( assessList.end() ),
  stoppingFunction( (QAFunction)0 )
{ 
  set_patch_type( PatchData::GLOBAL_PATCH, err, 0 );
  add_quality_assessment( metric, function, err );
  set_stopping_assessment( metric, function, err );
}

QualityAssessor::~QualityAssessor()
  { }

string QualityAssessor::get_QAFunction_name(
                              enum QualityAssessor::QAFunction fun)
{
  switch(fun){
    case(AVERAGE):
      return "Average   ";
    case(HISTOGRAM):
      return "Histogram of metric values: ";
    case(MAXIMUM):
      return "Maximum   ";
    case(MINIMUM):
      return "Minimum   ";
    case(RMS):
      return "RMS       ";
    case(STDDEV):
      return "Stan. Dev.";
    default:
      return "DEFAULT   ";
  };
}

double QualityAssessor::Assessor::get_average() const
{
  return count ? sum/count : 0;
}

double QualityAssessor::Assessor::get_rms() const 
{
  return count ? sqrt(sqrSum/count) : 0;
}

double QualityAssessor::Assessor::get_stddev() const
{
  double sqr = sqrSum/count - sum*sum/((double)count*count);
  return sqr < 0 ? 0 : sqrt(sqr);
}


/*!
    Several QualityMetric objects can be added to a single QualityAssessor
    object.  This allows to perform several quality assessments over a
    single mesh sweep.
    \param qm is the QualityMetric that will be used to evaluate the mesh
    quality
    \param func is the wrapper function used over the QualityMetric
    (min, max, etc..)
 */
void QualityAssessor::add_quality_assessment(QualityMetric* qm,
                                             int func,
                                             MsqError &/*err*/)
{ 
  msq_std::list<Assessor>::iterator iter;
  
  iter = find_or_add( qm );
  iter->funcFlags |= func;
  if (func&HISTOGRAM)
    iter->histogram.resize(DEFAULT_HISTOGRAM_INTERVALS+2);
}

list<QualityAssessor::Assessor>::iterator QualityAssessor::find_or_add( QualityMetric* qm )
{
  list<Assessor>::iterator iter;
  
    // If metric is already in list, find it
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->qualMetric == qm )
      break;
  
    // If metric not found in list, add it
  if (iter == assessList.end())
  {
    if (qm->get_metric_type() == QualityMetric::VERTEX_BASED)
    {
      assessList.push_back( Assessor(qm) );
      iter = --assessList.end();
    }
    else
    {
      assessList.push_front( Assessor(qm) );
      iter = assessList.begin();
    }
  }
  
  return iter;
}


/*!Sets which QualityMetric and QAFunction
combination is used to determine the value return from assess_mesh_quality().
It first ensures that the inputed QAFunction was not HISTOGRAM.  It then
calls add_quality_assessment with the given QualityMetric and QAFunction,
to ensure that this combination will be computed.  Finally, it sets
the stoppingMetric pointer and the stoppingFunction data members.
\param qm Pointer to QualityMetric.     
\param func (QAFUNCTION) Wrapper function for qm (e.g. MINIMUM, MAXIMUM,...).
    */
void QualityAssessor::set_stopping_assessment(QualityMetric* qm,
                                              QAFunction func,
                                              MsqError &err)
{
  if(func==HISTOGRAM){
    MSQ_SETERR(err)("HISTOGRAM DOES NOT GIVE A VALID RETURN VALUE", MsqError::INVALID_ARG);
    return;
  }
  else if (func == NO_FUNCTION) {
    MSQ_SETERR(err)("No function specified for stopping assessment", MsqError::INVALID_ARG);
    return;
  }
  
  stoppingMetric = find_or_add( qm );
  stoppingFunction = func;
}


/*! 
Checks first to see if the QualityMetric, qm, has been added to this
QualityAssessor, and if it has not, adds it.  It then adds HISTOGRAM as a
QAFunciton for that metric.  It then sets the minimum and maximum values
for the histogram.
\param qm Pointer to the QualityMetric to be used in histogram.
\param min_val (double) Minimum range of histogram.
\param max_val (double) Maximum range of histogram.
\param intervals Number of histogram intervals
    */
void QualityAssessor::add_histogram_assessment( QualityMetric* qm,
                                                double min_val, 
                                                double max_val,
                                                int intervals,
                                                MsqError &err )
{
  if (min_val >= max_val || intervals < 1) {
    MSQ_SETERR(err)("Invalid histogram range.", MsqError::INVALID_ARG );
    return;
  }
  
  list<Assessor>::iterator assessor = find_or_add( qm );
  assessor->funcFlags |= QualityAssessor::HISTOGRAM;
  assessor->histMin = min_val;
  assessor->histMax = max_val;
  assessor->histogram.resize( intervals + 2 );
} 



/*! 
  Computes the quality data for a given
  MeshSet, ms. What quality information is calculated, depends
  on what has been requested through the use of the QualityAssessor
  constructor, add_quality_assessment(), and set_stopping_assessment().
  The resulting data is printed in a table unless disable_printing_results()
  has been called.  The double returned depends on the QualityMetric
  and QAFunction "return" combination, which can be set using
  set_stopping_assessemnt().
  \param ms (const MeshSet &) MeshSet used for quality assessment.
 */
double QualityAssessor::loop_over_mesh(MeshSet &ms, MsqError& err)
{
    // Clear out any previous data
  reset_data();
  
    // Check for any metrics for which a histogram is to be 
    // calculated and for which the user has not specified 
    // minimum and maximum values.  
    // Element-based metrics are first in list, followed
    // by vertex-based metrics.  Find first vertex-based
    // metric also such that element metrics go from
    // assessList.begin() to elem_end and vertex metrics
    // go from elem_end to assessList.end()
  list<Assessor>::iterator elem_end = assessList.end();
  bool need_second_pass_for_elements = false;
  bool need_second_pass_for_vertices = false;
  list<Assessor>::iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
  {
    if (iter->get_metric()->get_metric_type() == QualityMetric::VERTEX_BASED)
      break;

    if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
      need_second_pass_for_elements = true;
  }
  elem_end = iter;
  for ( ; iter != assessList.end(); ++iter)
  {
    if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
      need_second_pass_for_vertices = true;
  }
  
  list<Assessor> histogramList;
  
    // Do element-based metrics
  if (assessList.begin() != elem_end)
  {
    bool first_pass = false;
    do { // might need to loop twice to calculate histograms
      first_pass = !first_pass;
     
      PatchData* pd;
      PatchData local_patch;
      no_culling_method();
      bool more_mesh = true;
      if (get_global_patch() == 0) {
        pd = &local_patch;
        more_mesh=ms.get_next_patch(*pd, this, err);  MSQ_ERRZERO(err);
      }
      else {
        pd = get_global_patch();
      }
      
        //until there are no more patches
        //there is another get_next_patch at
        //the end of this loop
      while (more_mesh)
      {
        for (unsigned i = 0; i < pd->num_elements(); ++i)
        {
          for (iter = assessList.begin(); iter != elem_end; ++iter)
          {
              // If first pass, get values for all metrics
            if (first_pass)
            {
              double value;
              bool valid = iter->get_metric()->evaluate_element( *pd, 
                                                           &pd->element_by_index(i),
                                                           value, err );
                                                           MSQ_ERRZERO(err);
              
              iter->add_value(value);
              if (!valid) 
                iter->add_invalid_value();
            }
              // If second pass, only do metrics for which the
              // histogram hasn't been calculated yet.
            else if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
            {
              double value;
              iter->get_metric()->evaluate_element( *pd, 
                                                    &pd->element_by_index(i),
                                                    value, err );
                                                    MSQ_ERRZERO(err);
              
              iter->add_hist_value(value);
            }
          }
        }
        
           // If dealing with local patches, get next element group (PatchData object)
        if (get_patch_type() != PatchData::GLOBAL_PATCH)
          more_mesh = ms.get_next_patch(*pd,this, err); MSQ_ERRZERO(err);
          //Michael:: Since we are doing global right now:
          //Remove this when no longer doing global
        more_mesh=false;
      }
  
        // Fix up any histogram ranges which were calculated
      for (iter = assessList.begin(); iter != elem_end; ++iter)
        if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
          if (first_pass)
            iter->calculate_histogram_range();
// Uncomment the following to have the QA keep the first
// calculated histogram range for all subsequent iterations.
//          else
//            iter->haveHistRange = true;
    
    } while (first_pass && need_second_pass_for_elements);
  }
      
    
      // Do vertex-based metrics
  if (assessList.end() != elem_end)
  {
    bool first_pass = false;
    do { // might need to loop twice to calculate histograms
      first_pass = !first_pass;
     
        //construct the patch we will send to get_next_patch
      PatchData* pd;
      PatchData local_patch;
      no_culling_method();
      bool more_mesh = true;
      if (get_global_patch() == 0) {
        pd = &local_patch; 
        more_mesh=ms.get_next_patch(*pd, this, err);  MSQ_ERRZERO(err);
      }
      else {
        pd = get_global_patch();
      }
      
        //until there are no more patches
        //there is another get_next_patch at
        //the end of this loop
      while (more_mesh)
      {
        for (unsigned i = 0; i < pd->num_vertices(); ++i)
        {
          for (iter = elem_end; iter != assessList.end(); ++iter)
          {
              // If first pass, get values for all metrics
            if (first_pass)
            {
              double value;
              bool valid = iter->get_metric()->evaluate_vertex( *pd, 
                                                           &pd->vertex_by_index(i),
                                                           value, err );
                                                           MSQ_ERRZERO(err);
              
              iter->add_value(value);
              if (!valid)
                iter->add_invalid_value();
            }
              // If second pass, only do metrics for which the
              // histogram hasn't been calculated yet.
            else if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
            {
              double value;
              iter->get_metric()->evaluate_vertex( *pd, 
                                                   &pd->vertex_by_index(i),
                                                   value, err );
                                                   MSQ_ERRZERO(err);
              
              iter->add_hist_value(value);
            }
          }
        }
        
        if (get_patch_type() != PatchData::GLOBAL_PATCH)
          more_mesh = ms.get_next_patch(*pd,this, err); MSQ_ERRZERO(err);
          //Michael:: Since we are doing global right now:
          //Remove this when no longer doing global
        more_mesh=false;
      }
  
        // Fix up any histogram ranges which were calculated
      for (iter = elem_end; iter != assessList.end(); ++iter)
        if (iter->funcFlags&HISTOGRAM && !iter->haveHistRange)
          if (first_pass)
            iter->calculate_histogram_range();
// Uncomment the following to have the QA keep the first
// calculated histogram range for all subsequent iterations.
//          else
//            iter->haveHistRange = true;
    
    } while (first_pass && need_second_pass_for_vertices);
  }  
  
  
    // Print results, if requested
  if (printSummary)
    print_summary( this->outputStream );
  
    // If no stopping function, just return zero
  if (!stoppingFunction)
    return 0.0;
  
    // Otherwise return requested value
  if      (stoppingFunction & STDDEV)
    return stoppingMetric->get_stddev();
  else if (stoppingFunction & AVERAGE)
    return stoppingMetric->get_average();
  else if (stoppingFunction & MAXIMUM)
    return stoppingMetric->get_maximum();
  else if (stoppingFunction & MINIMUM)
    return stoppingMetric->get_minimum();
  else if (stoppingFunction & RMS)
    return stoppingMetric->get_rms();
  else 
    MSQ_SETERR(err)("Invalid stopping function for QualityAssessor", 
                    MsqError::INVALID_STATE);
  
  return 0.0;
}

bool QualityAssessor::invalid_elements( ) const
{
  bool result = false;
  msq_std::list<Assessor>::const_iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->get_invalid_element_count())
      result = true;
  return result;
}

void QualityAssessor::reset_data() 
{
  msq_std::list<Assessor>::iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    iter->reset_data();
}

QualityAssessor::Assessor::Assessor( QualityMetric* metric )
  : qualMetric(metric),
    funcFlags(0),
    haveHistRange(false),
    histMin(1.0),
    histMax(0.0)
{
  reset_data();
}

const QualityAssessor::Assessor* QualityAssessor::get_results( QualityMetric* metric ) const
{
  msq_std::list<Assessor>::const_iterator iter;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->get_metric() == metric)
      return &*iter;
  return 0;
}


void QualityAssessor::Assessor:: get_histogram( double& lower_bound_out,
                                                double& upper_bound_out,
                                                msq_std::vector<int>& counts_out,
                                                MsqError& err ) const 
{
  if ( !(funcFlags & QualityAssessor::HISTOGRAM) )
  {
    MSQ_SETERR(err)("No histogram calculated.", MsqError::INVALID_STATE);
    return;
  }

  if (haveHistRange) {
    lower_bound_out = histMin;
    upper_bound_out = histMax;
  }
  else {
    lower_bound_out = minimum;
    upper_bound_out = maximum;
  }
  
  counts_out = histogram;
}

void QualityAssessor::Assessor::reset_data()
{
  count = 0;
  sum = 0;
  maximum = -HUGE_VAL;
  minimum = HUGE_VAL;
  sqrSum = 0;
  numInvalid = 0;
  memset( &histogram[0], 0, sizeof(int)*histogram.size() );
}

void QualityAssessor::Assessor::add_value( double metric_value )
{
  sum += metric_value;
  sqrSum += metric_value*metric_value;
  if (metric_value > maximum)
    maximum = metric_value;
  if (metric_value < minimum)
    minimum = metric_value;
    // Only add value to histogram data from this function if
    // the user has specified the range.  If user has not 
    // specified the range, QualityAssessor will call add_hist_value()
    // directly once the range has been calculated.
  if (funcFlags & QualityAssessor::HISTOGRAM && haveHistRange)
    add_hist_value( metric_value );
  
  ++count;
}

void QualityAssessor::Assessor::add_invalid_value()
{
  ++numInvalid;
}

void QualityAssessor::Assessor::add_hist_value( double metric_value )
{
    // Width of one interval in histogram
  double step = (histMax - histMin) / (histogram.size()-2);
    
    // First and last values in array are counts of values
    // outside the user-specified range of the histogram
    // (below and above, respectively.)
  if (metric_value < histMin)
    ++histogram[0];
  else if (metric_value > histMax)
    ++histogram[histogram.size()-1];
  else
  {
      // Calculate which interval the value is in.  Add one
      // because first entry is for values below user-specifed
      // minimum value for histogram.
    unsigned cell;
    if (step > DBL_EPSILON)
      cell = 1+(unsigned)((metric_value - histMin) / step);
    else
      cell = 1;
      
      // If value exactly equals maximum value, put in last
      // valid interval, not the count of values above the
      // maximum.
    if (cell + 1 == histogram.size())
      --cell;
      // Add value to interval.
    ++histogram[cell];
  }
}

void QualityAssessor::Assessor::calculate_histogram_range()
{
  double step = (maximum - minimum) / (histogram.size() - 2);
  double size = pow( 10, ceil(log10(step)) );
  histMin = size * floor( minimum / size );
  histMax = size *  ceil( maximum / size );
}  

void QualityAssessor::print_summary( msq_stdio::ostream& stream ) const
{
  const int NAMEW = 19;  // Width of name column in table output
  const int NUMW = 12;   // Width of value columns in table output
  
    // Print title
  stream << msq_stdio::endl 
         << "************** " 
         << qualityAssessorName
         << " Summary **************"
         << msq_stdio::endl
         << msq_stdio::endl;
         
    // Get union of function flags, and list any metrics with invalid values
  msq_std::list<Assessor>::const_iterator iter;
  unsigned flags = 0;
  int invalid_count = 0;
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
  {
    flags |= (iter->funcFlags & ~HISTOGRAM);
    
    if (iter->get_invalid_element_count())
    {
      ++invalid_count;
      
      stream << "  " << iter->get_invalid_element_count()
             << " OF " << iter->get_count()
             << " VALUES ARE INVALID FOR " 
             << iter->get_metric()->get_name() 
             << msq_stdio::endl << msq_stdio::endl;
    }
  }
  
  if (0 == invalid_count) {
    stream << "  No invalid values for any metric." 
           << msq_stdio::endl << msq_stdio::endl;
  }
  
    // If printing any values
  if (flags) 
  {
      // Print table header line
    stream << msq_stdio::setw(NAMEW) << "metric";
    if (flags & MINIMUM)
      stream << msq_stdio::setw(NUMW) << "minimum";
    if (flags & AVERAGE)
      stream << msq_stdio::setw(NUMW) << "average";
    if (flags & RMS)
      stream << msq_stdio::setw(NUMW) << "rms";
    if (flags & MAXIMUM)
      stream << msq_stdio::setw(NUMW) << "maximum";
    if (flags & STDDEV)
      stream << msq_stdio::setw(NUMW) << "std.dev.";
    stream << msq_stdio::endl;

      // Print out values for each assessor
    for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    {
        // If no output (other than histogram) for this metric, skip it
      if (!(iter->funcFlags & ~HISTOGRAM))
        continue;
      
        // Name column
      stream << msq_stdio::setw(NAMEW) << iter->get_metric()->get_name();

        // Value columns
      if (flags & MINIMUM)
      {
        if (iter->funcFlags & MINIMUM)
          stream << msq_stdio::setw(NUMW) << iter->get_minimum();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      if (flags & AVERAGE)
      {
        if (iter->funcFlags & AVERAGE)
          stream << msq_stdio::setw(NUMW) << iter->get_average();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      if (flags & RMS)
      {
        if (iter->funcFlags & RMS)
          stream << msq_stdio::setw(NUMW) << iter->get_rms();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      if (flags & MAXIMUM)
      {
        if (iter->funcFlags & MAXIMUM)
          stream << msq_stdio::setw(NUMW) << iter->get_maximum();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      if (flags & STDDEV)
      {
        if (iter->funcFlags & STDDEV)
          stream << msq_stdio::setw(NUMW) << iter->get_stddev();
        else
          stream << msq_stdio::setw(NUMW) << " ";
      }
      stream << msq_stdio::endl;
    } // for (assessList)
  } // if (flags)
  
  for (iter = assessList.begin(); iter != assessList.end(); ++iter)
    if (iter->funcFlags & HISTOGRAM)
      iter->print_histogram( stream );
}


void QualityAssessor::Assessor::print_histogram( msq_stdio::ostream& stream ) const
{
  // Portability notes:
  //  Use log10 rather than log10f because the float variations require
  //  including platform-dependent headers on some platforms.  
  //  Explicitly cast log10 argument to double because some platforms
  //  have overloaded float and double variations in C++ making an 
  //  implicit cast from an integer ambiguous.
  
  const char GRAPH_CHAR = '=';  // Character used to create bar graphs
  const int FLOATW = 12;        // Width of floating-point output
  const int GRAPHW = 50;        // Width of bar graph
  
    // range is either user-specified (histMin & histMax) or
    // calculated (minimum & maximum)
  double min, max;
  //if (haveHistRange) {
    min = histMin;
    max = histMax;
  //}
  //else {
  //  min = minimum;
  //  max = maximum;
  //}
    // Witdh of one interval of histogram
  double step = (max - min) / (histogram.size()-2);
  
    // Find maximum value for an interval of the histogram
  unsigned i;
  int max_interval = 1;
  for (i = 0; i < histogram.size(); ++i)
    if (histogram[i] > max_interval)
      max_interval = histogram[i];
  
  if (0 == max_interval)
    return; // no data 
  
    // Calculate width of field containing counts for 
    // histogram intervals (log10(max_interval)).
  int num_width = 1;
  for (int temp = max_interval; temp > 0; temp /= 10)
    ++num_width;

    // Create an array of bar graph characters for use in output
  char graph_chars[GRAPHW+1];
  memset( graph_chars, GRAPH_CHAR, sizeof(graph_chars) );
  
    // Check if bar-graph should be linear or log10 plot
    // Do log plot if standard deviation is less that 1.5
    // histogram intervals.
  bool log_plot = false;
  double stddev = get_stddev();
  if (stddev > 0 && stddev < 2.0*step)
  {
    int new_interval = (int)(log10((double)(1+max_interval)));
    if (new_interval > 0) {
      log_plot = true;
      max_interval = new_interval;
    }
  }

  
    // Write title
  stream << msq_stdio::endl << "   " << get_metric()->get_name() << " histogram:";
  if (log_plot)
    stream << " (log10 plot)";
  stream << msq_stdio::endl;

  
    // For each interval of histogram
  for (i = 0; i < histogram.size(); ++i)
  {
      // First value is the count of the number of values that
      // were below the minimum value of the histogram.
    if (0 == i)
    {
      if (0 == histogram[i])
        continue;
      stream << setw(FLOATW) << "under min";
    }
      // Last value is the count of the number of values that
      // were above the maximum value of the histogram.
    else if (i+1 == histogram.size())
    {
      if (0 == histogram[i])
        continue;
      stream << setw(FLOATW) << "over max";
    }
      // Anything else is a valid interval of the histogram.
      // Print the lower bound for each interval.
    else
    {
      stream << "   " << setw(FLOATW) << min + (i-1)*step;
    }
    
      // Print interval count.
    stream << ": " << setw(num_width) << histogram[i] << ": ";
    
      // Print bar graph
    
      // First calculate the number of characters to output
    int num_graph;
    if (log_plot)
      num_graph = GRAPHW * (int)log10((double)(1+histogram[i])) / max_interval;
    else
      num_graph = GRAPHW * histogram[i] / max_interval;
      
      // print num_graph characters using array of fill characters.
    graph_chars[num_graph] = '\0';
    stream << graph_chars << msq_stdio::endl;
    graph_chars[num_graph] = GRAPH_CHAR;
  }
  
  stream << msq_stdio::endl;
}
 


} //namespace Mesquite
