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

/*! \file InstructionQueue.hpp

Header file for the Mesquite::InstructionQueue class

  \author Thomas Leurent
  \date   2002-05-01
 */


#ifndef MSQ_INSTRUCTION_QUEUE_HPP
#define MSQ_INSTRUCTION_QUEUE_HPP

#include "Mesquite.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#else
#  include <list>
#endif

namespace Mesquite {

  class MsqError;
  class QualityImprover;
  class QualityAssessor;
  class MeshSet;
  class PatchDataUser;
  class PatchData;
  class TargetCalculator;

  /*! \class InstructionQueue
    \brief An InstructionQueue object gathers Mesquite Instructions and ensures
           that the instruction queue is coherent for mesh improvement and/or
           mesh quality assessment purposes.

           The user can instantiate several InstructionQueue objects to be used
           with various MeshSet objects.
           
           The most commonly used functions are:
           -# add_preconditioner(...)
           -# add_quality_assessor(...)
           -# set_master_quality_improver(...)
           -# run_instructions(...)
  */
  class InstructionQueue
  {

  public:
    InstructionQueue();

    virtual ~InstructionQueue() {};
    
    void add_target_calculator( TargetCalculator* tc, MsqError& err );
    
    void add_preconditioner(QualityImprover* instr, MsqError &err);
    void remove_preconditioner(size_t index, MsqError &err);
    void insert_preconditioner(QualityImprover* instr, size_t index, MsqError &err);
    
    void add_quality_assessor(QualityAssessor* instr, MsqError &err);
    void remove_quality_assessor(size_t index, MsqError &err);
    void insert_quality_assessor(QualityAssessor* instr, size_t index, MsqError &err);
    
    void set_master_quality_improver(QualityImprover* instr, MsqError &err);
    
    void disable_automatic_quality_assessment()
       { autoQualAssess = false; }
    void enable_automatic_quality_assessment()
       { autoQualAssess = true; }
    
    void disable_automatic_midnode_adjustment()
       { autoAdjMidNodes = false; }
    void enable_automatic_midnode_adjustment()
       { autoAdjMidNodes = true; }

      //!This function is virtual so that it may be redefined in the
      //! wraper classes.
    virtual void run_instructions(MeshSet &msc, MsqError &err);
    void clear();

      /**\brief Generate SIGFPE whenever a floating point exception occurs
       *
       * Generate a FPE signal when overflow, divbyzero, etc. occur
       * during floating-point arithmatic.  This is intended for debugging
       * purposes only, as enabling this will typically result in a 
       * crash when such arithmatic errors occur.
       *
       * If this option is enabled, Mesquite will attempt to set 
       * platform-specific flags such that a SIGFPE is generated for
       * floating point errors while the instruction queue is running.
       * If this option ins disabled, Mesquite will not change the
       * flags.  There is no option to explicitly disable such flags
       * because that is the default behavior on most platforms, and
       * presumably if the application has enabled such flags it makes
       * little sense to disable them while Mesquite is running.
       *
       * This functionality may not be supported on all platforms.  If
       * it is not supported, this option has no effect.
       */
    void trap_floating_point_exception( bool enable )
      { trapFPE = enable; }
    bool trap_floating_point_exception() const
      { return trapFPE; }
  
    
  protected:
    
  private:
    msq_std::list<PatchDataUser*>::iterator clear_master(MsqError &err);

    msq_std::list<PatchDataUser*> instructions;

    bool autoQualAssess;
    bool autoAdjMidNodes;
    
    size_t nbPreConditionners;
    bool isMasterSet;
    size_t masterInstrIndex; //!< 0-based. Keeping an index instead of an iterator
                             //!< in case list is reallocated

    PatchData* globalPatch; //!< Used to prevent reallocating a global patch
                            //!< for successive global algorithms.    
                            
    bool trapFPE;
  };


} //namespace


#endif // InstructionQueue_hpp
