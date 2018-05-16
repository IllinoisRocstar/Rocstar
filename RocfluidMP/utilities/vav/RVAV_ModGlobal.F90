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
!******************************************************************************
!
! Purpose: define globalRVAV variables and globalRVAV input data.
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: RVAV_ModGlobal.F90,v 1.4 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

MODULE RVAV_ModGlobal

  USE ModDataTypes
  IMPLICIT NONE

  TYPE t_compare

! RVAV_Comparisons ----------------------------------------------------------

    INTEGER :: operationS1
    INTEGER :: VariableIndexS1
    INTEGER :: blockS1
    INTEGER :: ibegS1, iendS1, ijumpS1, &
               jbegS1, jendS1, jjumpS1, &
               kbegS1, kendS1, kjumpS1

    INTEGER :: operationS2
    INTEGER :: VariableIndexS2
    INTEGER :: blockS2
    INTEGER :: ibegS2, iendS2, ijumpS2, &
               jbegS2, jendS2, jjumpS2, &
               kbegS2, kendS2, kjumpS2

  END TYPE t_compare

! global variable -------------------------------------------------------------

  TYPE t_RVAV_Global
  
    CHARACTER(CHRLEN) :: caseName
    INTEGER            :: nComparisons
    
! - the following fields correspond to stream1
    INTEGER                  :: FlowTypeS1
    INTEGER                  :: FileTypeS1
    INTEGER                  :: GridFormatS1
    INTEGER                  :: SolutFormatS1
    INTEGER                  :: iCOffS1,ijCOffS1
    
! - the following fields correspond to stream2
    INTEGER                  :: FlowTypeS2
    INTEGER                  :: FileTypeS2
    INTEGER                  :: GridFormatS2
    INTEGER                  :: SolutFormatS2
    INTEGER                  :: nRegionsS2
    INTEGER                  :: SimilarityTypeS2
    INTEGER                  :: iCOffS2,ijCOffS2
    
! - these variables hold the extracted values from streams 1 & 2
    REAL(RFREAL), POINTER    :: evS1(:,:,:), evS2(:,:,:) 
                                                        
! - the number of comparisons to be made in both streams and the fields therein

    TYPE(t_compare), POINTER :: RVAVcompare(:) !as many as the no of comparisons
    
  END TYPE t_RVAV_Global

  TYPE(t_RVAV_global) :: globalRVAV

END MODULE RVAV_ModGlobal

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ModGlobal.F90,v $
! Revision 1.4  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2002/06/15 17:52:28  f-najjar
! Bug fix for RVAV_extractVariables
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!
!******************************************************************************






