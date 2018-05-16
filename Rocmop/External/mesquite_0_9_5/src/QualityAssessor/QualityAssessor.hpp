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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file QualityAssessor.hpp

Header file for the Mesquite::QualityAssessor class

  \author Thomas Leurent
  \date   2002-05-01
  \author Jason Kraftcheck
  \date   2005-03-09
 */


#ifndef MSQ_QUALITYASSESSOR_HPP
#define MSQ_QUALITYASSESSOR_HPP

#include <math.h>

#include "Mesquite.hpp"
#include "PatchDataUser.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <string.h>
#else
#  include <list>
#  include <string>
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <ostream.h>
#else
#  include <iosfwd>
#endif


namespace Mesquite 
{

   class QualityMetric;
   class MsqError;
   class MeshSet;

  /*! \class QualityAssessor

      \brief A QualityAssessor instance can be inserted into an 
      \ref InstructionQueue to calculate and summarize registered
      \ref QualityMetrics for the mesh.

      The relevant quality assessments are set by the user or
      automatically (default) by Mesquite when an InstructionQueue
      object is used.  If the mesh has been changed (improved),
      it is often useful to reuse the same QualityAssessor object
      to reassess the mesh quality.
      
      The \ref QAFunction flags passed for each metric control the output
      for that metric.  If no \ref QAFuction flags are passed, no results
      are printed for the metric except the count of invalid 
      vertices/elements reported by that metric, if any.
      
      The "stopping assessor" and "stopping function", if set,
      determinte the value reported to Mesquite for the overall
      run of of the QualityAssessor.
      
      All summary data except the histogram is accumulated for all
      registered metrics, and can be accessed by the calling application.
      Histogram data is accumulated only if the \ref HISTOGRAM 
      \ref QUFunction output is requested for the metric (which is 
      impled by calling \ref add_histogram_assessment.
  */
  class QualityAssessor : public PatchDataUser
  {
  public:
    
    
    /*! \enum QAFunction
      type of function used in conjunction with QualityMetric to compute mesh quality */ 
    enum QAFunction {
       NO_FUNCTION = 0,
       AVERAGE=1,
       HISTOGRAM=2,
       MAXIMUM=4, 
       MINIMUM=8,
       RMS=16,
       STDDEV=32,
       ALL_MEASURES=255
    };
    
    static msq_std::string get_QAFunction_name(enum QualityAssessor::QAFunction);
    
    //! Constructor - output to std::cout
    QualityAssessor( msq_std::string name = "QualityAssessor" );
    
    //! Constructor - specified output stream 
    QualityAssessor( msq_stdio::ostream& output_stream,
                     msq_std::string name = "QualityAssessor" );
                     
    //! Constructor - initial stopping assessement and specified output stream
    QualityAssessor( QualityMetric* metric, QAFunction function,
                     msq_stdio::ostream& output_stream,
                     MsqError& err,
                     msq_std::string name = "QualityAssessor" );

                     
    //! Constructor - initial stopping assessement
    QualityAssessor( QualityMetric* metric, QAFunction function,
                     MsqError& err,
                     msq_std::string name = "QualityAssessor" );

      //!Destructor
    ~QualityAssessor();
    
      //! Provides a name to the QualityAssessor (use it for default name in constructor).
    void set_name(msq_std::string name) { qualityAssessorName = name; };
      //! Retrieves the QualityAssessor name. A default name should be set in the constructor.
    virtual msq_std::string get_name()  { return qualityAssessorName; }

    virtual AlgorithmType get_algorithm_type() { return QUALITY_ASSESSOR; }
    
      //! Adds a quality metric and a wrapper function (min, max, ...).
    void add_quality_assessment( QualityMetric* qm, 
                                 int function_flags, 
                                 MsqError &err);

      /*!Sets the QualityMetric and QAFunction combination that will
        be returned when loop_over_mesh is called.
      */
    void set_stopping_assessment( QualityMetric* qm, 
                                  QAFunction func,
                                  MsqError &err );

      //! Add a quality metric for which the histogram is to 
      //! be calculated, and set histogram parameters.
    void add_histogram_assessment( QualityMetric* qm, 
                                   double min, 
                                   double max,
                                   int intervals,
                                   MsqError& err );
    
      //! Does one sweep over the mesh and assess the quality with the metrics previously added.
    virtual double loop_over_mesh(MeshSet &ms, MsqError &err);

      //! Do not print results of assessment.
    void disable_printing_results()
       {
         printSummary = false;
       }
      
      //! Print accumulated summary data to specified stream. 
    void print_summary( msq_stdio::ostream& stream ) const;
    
      //! True if any metric evaluated to an invalid value
      //! for any element
    bool invalid_elements() const;
       
      //! Reset calculated data 
    void reset_data();
    
    /** \brief Per-metric QualityAssessor data
     *
     * The Assessor class holds QualityAssessor data for
     * each metric added by the calling application, including
     * a pointer to the metric instance, \ref QAFunction flags
     * dictating what is to be calculated and output, histogram
     * parameters, and the variables used to accumulate results as
     * the \ref QualityAssessor is running.  It also provides 
     * methods to access the calculated data once the QualityAssessor
     * pass is completed.
     */
    class Assessor
    {
      public:
      
        Assessor( QualityMetric* metric );
        
        double get_average() const ;
        double get_maximum() const { return maximum; }
        double get_minimum() const { return minimum; }
        double get_rms()     const ;
        double get_stddev()  const ;
        int get_count() const { return count; }
        
        int get_invalid_element_count() const { return numInvalid; }
        
        /** Get historgram of data, if calculated.
         *\param lower_bound_out  The lower bound of the histogram
         *\param upper_bound_out  The upper bound of the histogram
         *\param counts_out       An array of counts of elements where
         *              the first entry is the number of elements for
         *              which the metric is below the lower bound, the
         *              last entry is the number of elements above the
         *              upper bound, and all other values are the counts
         *              for histogram intervals between the lower and
         *              upper bounds.
         */
        void get_histogram( double& lower_bound_out,
                            double& upper_bound_out,
                            msq_std::vector<int>& counts_out,
                            MsqError& err ) const;
                            
        /** Reset all calculated data */
        void reset_data();
       
        /** Print the histogram */
        void print_histogram( msq_stdio::ostream& ) const;

        /** Get the QualityMetric */
        QualityMetric* get_metric() const { return qualMetric; }
        
        /** Add a value to the running counts */
        void add_value( double metric_value );
        
        /** Add a value to the hisogram data */
        void add_hist_value( double metric_value );
        
        /** Note invalid result */
        void add_invalid_value() ;
        
        /** If range of histogram has not yet been determined,
          * calculate it from the min/max values.
          */
        void calculate_histogram_range();
        
      private:
      
        friend class QualityAssessor;
        
        QualityMetric *const qualMetric; //< The quality metric
        unsigned funcFlags;             //< What to calculate
        
        unsigned long count;  //< The total number of times the metric was evaluated
        
        double sum;       //< The sum of the metric over all elements
        double maximum;   //< The maximum of the metric
        double minimum;   //< The minimum value of the metric
        double sqrSum;    //< The sum of the square of the metric values
        unsigned long numInvalid;  //< Count of invalid metric values
        
        /** The histogram counts, where the first and last values are
         * counts of values below the lower bound and above the upper
         * bound, respectively.  The remaining values are the histogram
         * counts.
         */
        bool haveHistRange;
        double histMin;   //< Lower bound of histogram
        double histMax;   //< Upper bound of histogram
        msq_std::vector<int> histogram;
    };    
        
    /** \brief Request summary data for a specific QualityMetric 
     * This method allows the application to request the summary
     * data for a metric it has registered with the QualityAssessor.
     * If the passed QualityMetric has not been registered with the
     * QualityAssessor instance, NULL is returned.
     */
    const Assessor* get_results( QualityMetric* metric ) const;
    
    /** \brief Get list of all summary data.
     *  Return a const reference to the internal list of 
     *  calculated data.
     */
   const msq_std::list<Assessor>& get_all_results() const
      { return assessList; }
      
  private:
  
    /** Find an Assessor corresponding to the passed
     *  QualityMetric, or create it if is not found in
     *  the list.
     */
    msq_std::list<Assessor>::iterator find_or_add( QualityMetric* qm );
   
    /** Name */
    msq_std::string qualityAssessorName;  
    
    /** List of quality metrics and corresponding data */
    msq_std::list<Assessor> assessList;
   
    /** Stream to which to write summary of metric data */
    msq_stdio::ostream& outputStream;
    /** Disable printing */
    bool printSummary;
    
      /** Metric in \ref assessList to use as return value for loop_over_mesh */
    msq_std::list<Assessor>::iterator stoppingMetric;
      /** Value to use as return value for loop_over_mesh */
    QAFunction stoppingFunction;
  };

  
} //namespace


#endif // QualityAssessor_hpp
