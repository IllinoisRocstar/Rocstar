! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
! ******************************************************************************
!
! Purpose: Display minimum and maximum values of state vector for given local
!   domain.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: N/A.
!
! Notes: 
!   1. Depending on where this routine is called, one may not be able to 
!      switch to or from conservative variables because the arrays dv and gv 
!      are not allocated. To avoid this problem, this routine checks whether 
!      these arrays are allocated, and only carries out the conversion to 
!      primitive variables if that is the case. 
!   2. IMPORTANT: At present, the conversion from and back to conserved 
!      variables is commented out because it breaks perfect binary restart. 
!
! ******************************************************************************
!
! $Id: RFLU_PrintFlowInfo.F90,v 1.18 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PrintFlowInfo(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region 
  
  USE RFLU_ModConvertCv
  
  USE ModInterfaces, ONLY: RFLU_GetCvLoc, & 
                           RFLU_PrintLocInfo   

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Parameters
! ==============================================================================   

  TYPE(t_region), POINTER :: pRegion 

! ==============================================================================
! Local variables
! ==============================================================================
 
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: cvMixtPres,cvMixtXVel,cvMixtYVel,cvMixtZVel,errorFlag,uppLim
  INTEGER :: dummy(1) 
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: loc
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_global), POINTER :: global
  
  RCSIdentString = '$RCSfile: RFLU_PrintFlowInfo.F90,v $ $Revision: 1.18 $'
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PrintFlowInfo',&
  'RFLU_PrintFlowInfo.F90')

  IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing flow information...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal 
    IF ( global%flowType == FLOW_UNSTEADY ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                          global%currentTime
    END IF ! global%flowType
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pCv   => pRegion%mixt%cv
  pGrid => pRegion%grid

! ******************************************************************************
! Set upper limit: useful for checking of cell/dummy cell values
! ******************************************************************************

  uppLim = pGrid%nCells ! Only internal cells
!  uppLim = pGrid%nCellsTot ! Also dummy cells

! ******************************************************************************
! Convert to primitive state vector (only if dv and gv are allocated)
! ******************************************************************************

! TEMPORARY: Converting to another state vector format for printing purposes
! breaks perfect binary restart, because the cv values are overwritten, and on
! converting back to the original conserved variuables, it is not guaranteed 
! that the values are identical to the original ones. Statement further down
! to convert back is also commented out.
!  IF ( (ASSOCIATED(pRegion%mixt%dv) .EQV. .TRUE.) .AND. & 
!       (ASSOCIATED(pRegion%mixt%gv) .EQV. .TRUE.) ) THEN 
!    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)
!  END IF ! ASSOCIATED
! END TEMPORARY

! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************

  ALLOCATE(loc(pRegion%mixtInput%nCv,MIN_VAL:MAX_VAL),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'loc')
  END IF ! global%error

! ******************************************************************************
! Find locations of extrema: NOTE Asinine coding needed because of poor 
! FORTRAN interface for MINLOC and MAXLOC functions... 
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel ) 

! ==============================================================================
!   Incompressible fluid model
! ==============================================================================    

    CASE ( FLUID_MODEL_INCOMP ) 
      cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)
      cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)
      cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)
      cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)                  
        
      dummy = MINLOC(pCv(cvMixtXVel,1:uppLim))
      loc(cvMixtXVel,MIN_VAL) = dummy(1)

      dummy = MINLOC(pCv(cvMixtYVel,1:uppLim))
      loc(cvMixtYVel,MIN_VAL) = dummy(1)

      dummy = MINLOC(pCv(cvMixtZVel,1:uppLim))
      loc(cvMixtZVel,MIN_VAL) = dummy(1)

      dummy = MINLOC(pCv(cvMixtPres,1:uppLim))
      loc(cvMixtPres,MIN_VAL) = dummy(1)      
    
      dummy = MAXLOC(pCv(cvMixtXVel,1:uppLim))
      loc(cvMixtXVel,MAX_VAL) = dummy(1)

      dummy = MAXLOC(pCv(cvMixtYVel,1:uppLim))
      loc(cvMixtYVel,MAX_VAL) = dummy(1)

      dummy = MAXLOC(pCv(cvMixtZVel,1:uppLim))
      loc(cvMixtZVel,MAX_VAL) = dummy(1)

      dummy = MAXLOC(pCv(cvMixtPres,1:uppLim))
      loc(cvMixtPres,MAX_VAL) = dummy(1)    
        
! ==============================================================================
!   Compressible fluid model
! ==============================================================================    
    
    CASE ( FLUID_MODEL_COMP )
      dummy = MINLOC(pCv(CV_MIXT_DENS,1:uppLim))
      loc(CV_MIXT_DENS,MIN_VAL) = dummy(1)

      dummy = MINLOC(pCv(CV_MIXT_XMOM,1:uppLim))
      loc(CV_MIXT_XMOM,MIN_VAL) = dummy(1)

      dummy = MINLOC(pCv(CV_MIXT_YMOM,1:uppLim))
      loc(CV_MIXT_YMOM,MIN_VAL) = dummy(1)

      dummy = MINLOC(pCv(CV_MIXT_ZMOM,1:uppLim))
      loc(CV_MIXT_ZMOM,MIN_VAL) = dummy(1)

      dummy = MINLOC(pCv(CV_MIXT_ENER,1:uppLim))
      loc(CV_MIXT_ENER,MIN_VAL) = dummy(1)


      dummy = MAXLOC(pCv(CV_MIXT_DENS,1:uppLim))
      loc(CV_MIXT_DENS,MAX_VAL) = dummy(1)

      dummy = MAXLOC(pCv(CV_MIXT_XMOM,1:uppLim))
      loc(CV_MIXT_XMOM,MAX_VAL) = dummy(1)

      dummy = MAXLOC(pCv(CV_MIXT_YMOM,1:uppLim))
      loc(CV_MIXT_YMOM,MAX_VAL) = dummy(1)

      dummy = MAXLOC(pCv(CV_MIXT_ZMOM,1:uppLim))
      loc(CV_MIXT_ZMOM,MAX_VAL) = dummy(1)

      dummy = MAXLOC(pCv(CV_MIXT_ENER,1:uppLim))
      loc(CV_MIXT_ENER,MAX_VAL) = dummy(1)

! ==============================================================================
!   Default
! ==============================================================================    

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel
  
  
! ******************************************************************************
! Print locations of extrema
! ******************************************************************************
  
  SELECT CASE ( pRegion%mixtInput%fluidModel ) 

! ==============================================================================
!   Incompressible fluid model
! ==============================================================================    

    CASE ( FLUID_MODEL_INCOMP ) 
      cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)
      cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)
      cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)
      cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)   
 
      IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                      SOLVER_NAME,'X-velocity (m/s):', & 
                      MINVAL(pCv(cvMixtXVel,1:uppLim)), & 
                      MAXVAL(pCv(cvMixtXVel,1:uppLim)), & 
                      loc(cvMixtXVel,MIN_VAL),loc(cvMixtXVel,MAX_VAL)
        WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                      SOLVER_NAME,'Y-velocity (m/s):', & 
                      MINVAL(pCv(cvMixtYVel,1:uppLim)), & 
                      MAXVAL(pCv(cvMixtYVel,1:uppLim)), & 
                      loc(cvMixtYVel,MIN_VAL),loc(cvMixtYVel,MAX_VAL)
        WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                      SOLVER_NAME,'Z-velocity (m/s):', & 
                      MINVAL(pCv(cvMixtZVel,1:uppLim)), & 
                      MAXVAL(pCv(cvMixtZVel,1:uppLim)), & 
                      loc(cvMixtZVel,MIN_VAL),loc(cvMixtZVel,MAX_VAL)
        WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                      SOLVER_NAME,'Pressure (N/m^2):', & 
                      MINVAL(pCv(cvMixtPres,1:uppLim)), & 
                      MAXVAL(pCv(cvMixtPres,1:uppLim)), & 
                      loc(cvMixtPres,MIN_VAL),loc(cvMixtPres,MAX_VAL) 
      END IF !global%verbLevel 
  
! ==============================================================================
!   Compressible fluid model
! ==============================================================================    
    
    CASE ( FLUID_MODEL_COMP )
  
! ------------------------------------------------------------------------------
!     Conservative state vector
! ------------------------------------------------------------------------------
  
      SELECT CASE ( pRegion%mixt%cvState ) 
        CASE ( CV_MIXT_STATE_CONS )
          IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'Density (kg/m^3):    ', & 
                          MINVAL(pCv(CV_MIXT_DENS,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_DENS,1:uppLim)), & 
                          loc(CV_MIXT_DENS,MIN_VAL),loc(CV_MIXT_DENS,MAX_VAL)
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'X-momentum (kg/m^2s):', & 
                          MINVAL(pCv(CV_MIXT_XMOM,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_XMOM,1:uppLim)), & 
                          loc(CV_MIXT_XMOM,MIN_VAL),loc(CV_MIXT_XMOM,MAX_VAL)
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'Y-momentum (kg/m^2s):', & 
                          MINVAL(pCv(CV_MIXT_YMOM,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_YMOM,1:uppLim)), & 
                          loc(CV_MIXT_YMOM,MIN_VAL),loc(CV_MIXT_YMOM,MAX_VAL)
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'Z-momentum (kg/m^2s):', & 
                          MINVAL(pCv(CV_MIXT_ZMOM,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_ZMOM,1:uppLim)), & 
                          loc(CV_MIXT_ZMOM,MIN_VAL),loc(CV_MIXT_ZMOM,MAX_VAL)
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'Energy (N/m^2):      ', & 
                          MINVAL(pCv(CV_MIXT_ENER,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_ENER,1:uppLim)), & 
                          loc(CV_MIXT_ENER,MIN_VAL),loc(CV_MIXT_ENER,MAX_VAL)
          END IF !global%verbLevel

! ------------------------------------------------------------------------------
!       Primitive state vector
! ------------------------------------------------------------------------------

        CASE ( CV_MIXT_STATE_PRIM ) 
          IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'Density (kg/m^3):', & 
                          MINVAL(pCv(CV_MIXT_DENS,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_DENS,1:uppLim)), & 
                          loc(CV_MIXT_DENS,MIN_VAL),loc(CV_MIXT_DENS,MAX_VAL)
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'X-velocity (m/s):', & 
                          MINVAL(pCv(CV_MIXT_XVEL,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_XVEL,1:uppLim)), & 
                          loc(CV_MIXT_XVEL,MIN_VAL),loc(CV_MIXT_XVEL,MAX_VAL)
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'Y-velocity (m/s):', & 
                          MINVAL(pCv(CV_MIXT_YVEL,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_YVEL,1:uppLim)), & 
                          loc(CV_MIXT_YVEL,MIN_VAL),loc(CV_MIXT_YVEL,MAX_VAL)
            WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                          SOLVER_NAME,'Z-velocity (m/s):', & 
                          MINVAL(pCv(CV_MIXT_ZVEL,1:uppLim)), & 
                          MAXVAL(pCv(CV_MIXT_ZVEL,1:uppLim)), & 
                          loc(CV_MIXT_ZVEL,MIN_VAL),loc(CV_MIXT_ZVEL,MAX_VAL)
          END IF !global%verbLevel

          SELECT CASE ( pRegion%mixt%cvState )

! --------- Pressure as last entry ---------------------------------------------

            CASE ( CV_MIXT_STATE_DUVWP ) 
              IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
                WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                              SOLVER_NAME,'Pressure (N/m^2):', & 
                              MINVAL(pCv(CV_MIXT_PRES,1:uppLim)), & 
                              MAXVAL(pCv(CV_MIXT_PRES,1:uppLim)), & 
                              loc(CV_MIXT_PRES,MIN_VAL), &
                              loc(CV_MIXT_PRES,MAX_VAL)    
              END IF !global%verbLevel

! --------- Temperature as last entry ------------------------------------------

            CASE ( CV_MIXT_STATE_DUVWT ) 
              IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
                WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                              SOLVER_NAME,'Temperature (K): ', & 
                              MINVAL(pCv(CV_MIXT_TEMP,1:uppLim)), & 
                              MAXVAL(pCv(CV_MIXT_TEMP,1:uppLim)), & 
                              loc(CV_MIXT_TEMP,MIN_VAL), &
                              loc(CV_MIXT_TEMP,MAX_VAL)     
              END IF !global%verbLevel
            CASE DEFAULT  
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixt%cvState         
        CASE DEFAULT  
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)                    
      END SELECT ! pRegion%mixt%cvState  

! ==============================================================================
!   Default
! ==============================================================================    

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel
    
! ******************************************************************************
! Print out locations of cells at which extrema occur and deallocate temporary
! memory.
! ******************************************************************************
     
  IF ( global%verbLevel >= VERBOSE_MED ) THEN   
    CALL RFLU_PrintLocInfo(pRegion,loc,pRegion%mixtInput%nCv, & 
                           LOCINFO_MODE_SILENT,OUTPUT_MODE_MASTER_ONLY)
  END IF ! global%verbLevel
  
  DEALLOCATE(loc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'loc')
  END IF ! global%error
    
! ******************************************************************************
! Convert back to conservative state vector and  Check that cv has the correct 
! state: Defensive programming - DO NOT ADD ANY STATEMENTS AFTER THESE, THEY 
! MIGHT SCREW UP ROCFLU
! ******************************************************************************
 
! TEMPORARY: See comment above 
!  IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
!    CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)
!  END IF ! pRegion%mixt%cvState
! END TEMPORARY  
     
! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing flow information done.'
  END IF ! global%verbLevel    
   
  CALL DeregisterFunction(global)  
  
END SUBROUTINE RFLU_PrintFlowInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintFlowInfo.F90,v $
! Revision 1.18  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.17  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.16  2006/03/26 20:21:37  haselbac
! Removed error traps for GL model
!
! Revision 1.15  2004/11/14 19:40:50  haselbac
! Added printing of info for incompressible fluid model, cosmetics
!                                
! Revision 1.14  2004/01/22 16:04:48  haselbac                    
! Changed declaration to eliminate warning on ALC                 
!
! Revision 1.13  2003/12/04 03:23:49  haselbac                    
! Fixed bug in comment, cosmetic changes                          
!
! Revision 1.12  2003/06/04 22:03:12  haselbac                    
! Added argument to printing of locations                         
!
! Revision 1.11  2003/03/15 17:01:56  haselbac                    
! Increased precision, adapted subroutine call                    
!
! Revision 1.10  2003/01/28 15:37:52  haselbac                    
! Commented out calls to change state vector - see comments       
!
! Revision 1.9  2002/10/17 14:11:02  haselbac                     
! Added writing out of iRegionGlobal                              
!
! Revision 1.8  2002/10/16 21:12:06  haselbac                     
! Only print time for unsteady flows                              
!
! Revision 1.7  2002/10/12 14:50:05  haselbac                     
! Added WRITE statement for currentTime (for GENX use)            
!
! Revision 1.6  2002/09/09 14:15:01  haselbac                     
! global now under regions, distinguish between cons and prim cv  
!
! Revision 1.5  2002/07/25 14:43:12  haselbac                     
! Added upper limit, useful for checking of dummy cells           
!
! Revision 1.4  2002/06/17 13:31:22  haselbac                     
! Prefixed SOLVER_NAME to all screen output                       
!
! Revision 1.3  2002/05/04 19:10:49  haselbac                     
! Coded differently, did not compile on Turing                    
!
! Revision 1.2  2002/05/04 16:13:31  haselbac                     
! Added call to RFLU_PrintLocInfo                                 
!
! Revision 1.1  2002/04/11 18:42:25  haselbac                     
! Initial revision                                                
!
! ******************************************************************************







