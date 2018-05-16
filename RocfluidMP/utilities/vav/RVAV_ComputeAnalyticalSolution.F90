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
! Purpose: Generate analytical data stream (S2) for various flow fields.
!
! Description: currently supported formats are:
!              - ASCII 
!              - Binary
!
! Input: case name from the list of arguments
!
! Output: File for Analytical solution
!
! Notes: 
!
!******************************************************************************
!
! $Id: RVAV_ComputeAnalyticalSolution.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ComputeAnalyticalSolution( similarityType,regionsS1 )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset, &
        RFLO_GetDimensPhys, RFLO_CopyGeometryDummy, RFLO_CalcCellCentroids
  USE ModMPI
  USE ModParameters
  USE RVAV_ModGlobal
  USE RVAV_ModParameters
  USE RVAV_ModInterfaces, ONLY : RVAV_BlasiusSolution, RVAV_GAMMBumpSolution, &
                                 RVAV_ProudmanCulickSolution
  IMPLICIT NONE

#include "Indexing.h"
  
! ... parameters
  INTEGER :: similarityType
  TYPE (t_region), POINTER :: regionsS1(:)
  
! ... loop variables
  INTEGER :: i, j, k 
  
! ... local variables
  CHARACTER(CHRLEN)    :: msg
  CHARACTER(CHRLEN+23) :: fname

  TYPE(t_grid)   , POINTER :: grid
  TYPE(t_mixt)   , POINTER :: mixt
  TYPE(t_global) , POINTER :: global

  INTEGER :: iReg, iLev  
  INTEGER :: ipc, jpc, kpc, ibc, iec
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iCOff, ijCOff
  INTEGER :: iCellOffset, ijCellOffset
  INTEGER :: ibpc,iepc
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iVars, iNodes, jNodes, kNodes, nVari, nVars, errorFlag

  REAL(RFREAL), POINTER :: cofg(:,:)
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: cvSA
  
!******************************************************************************

  global => regionsS1(1)%global
  
  CALL RegisterFunction( global, 'RVAV_ComputeAnalyticalSolution',&
  'RVAV_ComputeAnalyticalSolution.F90' )
   
! open file in which stream2 is stored
          
  WRITE(fname,'(A,A)') TRIM(globalRVAV%casename)//'_s2.anl'
  OPEN (IF_RVAV_FILE_S2,file=fname,form='formatted',status='unknown',&
        iostat=errorFlag)

  global%error = errorFlag
  IF (global%error /= 0) &
  CALL ErrorStop( global, ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )  

! write number of regions 
  
  WRITE(IF_RVAV_FILE_S2,*,err=10) global%nRegions 
  
! calculate cell centroids ----------------------------------------------------
      
  IF (global%verbLevel /= VERBOSE_NONE) &
    WRITE(STDOUT,'(/,A)') 'Compute Cell Centroids...'
        
  DO iReg=1, global%nRegions

    iLev =  regionsS1(iReg)%currLevel
    grid => regionsS1(iReg)%levels(iLev)%grid
    mixt => regionsS1(iReg)%levels(iLev)%mixt

! - get cell and node dimensions

    CALL RFLO_GetDimensDummy( regionsS1(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regionsS1(iReg),iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)
      
    CALL RFLO_GetDimensPhys( regionsS1(iReg),iLev,ipcbeg,ipcend, &
                             jpcbeg,jpcend,kpcbeg,kpcend )
    CALL RFLO_GetCellOffset( regionsS1(iReg),iLev,iCOff,ijCOff )
    
    ibpc = IndIJK(ipcbeg,jpcbeg,kpcbeg,iCOff,ijCOff)
    iepc = IndIJK(ipcend,jdcend,kpcend,iCOff,ijCOff)
      
    ALLOCATE( grid%cofg(3,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
        
    IF (global%verbLevel /= VERBOSE_NONE) &
       WRITE(STDOUT,'(/,A,2X,I5)') ' nGridLevels = ',regionsS1(iReg)%nGridLevels

    CALL RFLO_CopyGeometryDummy( regionsS1(iReg) )
    CALL RFLO_CalcCellCentroids( regionsS1(iReg) )
      
  ENDDO   ! iReg

! start case selection for analytical solutions -------------------------------

  SELECT CASE (similarityType)
                                
! generates Input DataStream for C0
      
    CASE (RVAV_CULICK)
                
      IF (global%verbLevel /= VERBOSE_NONE) &
        WRITE(STDOUT,'(/,A)') 'Generates Input Deck for C0 ...'
        
      CALL RVAV_ProudmanCulickSolution( fname,regionsS1 )
         
! generates Input DataStream for GAMM Bump
           
    CASE (RVAV_GAMMBUMP)
        
      IF (global%verbLevel /= VERBOSE_NONE) &
        WRITE(STDOUT,'(/,A)') 'Generates Input Deck for GAMM Bump ...'
  
      CALL RVAV_GAMMBumpSolution( fname,regionsS1 )
          
! generates Input DataStream for Laminar Blasius
           
    CASE (RVAV_BLASIUS) 
        
      IF (global%verbLevel /= VERBOSE_NONE) &
        WRITE(STDOUT,'(/,A)') 'Generates Input Deck for Laminar Blasius ...'

      CALL RVAV_BlasiusSolution( fname,regionsS1 )
          
  END SELECT ! similarityType  

! close file ------------------------------------------------------------------
      
  CLOSE(IF_RVAV_FILE_S2,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
  CALL ErrorStop( global, ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  
  GOTO 999

10 CONTINUE
   CALL ErrorStop( global, ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) ) 
      
! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )
  
1000 FORMAT('Region ',I5,', grid level= ',I2,'.')

END SUBROUTINE RVAV_ComputeAnalyticalSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ComputeAnalyticalSolution.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:46:35  fnajjar
! Initial revision after changing case
!
! Revision 1.15  2003/11/20 16:40:41  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.11  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.10  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.9  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.8  2002/08/15 19:48:06  jblazek
! Implemented grid deformation capability.
!
! Revision 1.7  2002/07/31 02:34:57  f-najjar
! Split Analytical Solutions into individual routines
!
! Revision 1.6  2002/06/24 16:05:05  f-najjar
! Cleanup for file writing
!
! Revision 1.5  2002/06/24 15:48:09  f-najjar
! Included Blasius Profile Computation
!
! Revision 1.4  2002/06/21 14:11:15  f-najjar
! Clean up for Consistency
!
! Revision 1.3  2002/06/19 20:29:02  f-najjar
! Included GAMM Bump Computation
!
! Revision 1.2  2002/06/19 14:40:14  f-najjar
! Included verbLevel calls for cleanup
!
! Revision 1.1  2002/06/18 03:19:15  f-najjar
! Initial Import
!
!******************************************************************************







