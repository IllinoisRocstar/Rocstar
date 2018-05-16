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
! Purpose: Generate Proudman-Culick Analytical Solution.
!
! Description: 
!
! Input: number of variables and array
!
! Output: array for Analytical solution
!
! Notes: 
!
!******************************************************************************
!
! $Id: RVAV_ProudmanCulickSolution.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ProudmanCulickSolution( fname,regionsS1 )
  
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE ModMPI
  USE ModParameters
  USE RVAV_ModGlobal
  USE RVAV_ModParameters
  IMPLICIT NONE

#include "Indexing.h"
  
! ... parameter
  CHARACTER(*) :: fname  
  TYPE (t_region), POINTER :: regionsS1(:)

! ... loop variables
  INTEGER :: i, j, k 
      
! ... local variables
  CHARACTER(CHRLEN) :: msg
  
  REAL(RFREAL), PARAMETER :: massFluxC0 = 2.42_RFREAL
  REAL(RFREAL), PARAMETER :: tempWallC0 = 303.0_RFREAL
  REAL(RFREAL), PARAMETER :: pressureHeadEnd = 1.53451E+05_RFREAL
  REAL(RFREAL), PARAMETER :: heightC0 = 0.02

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
  INTEGER :: iVars, iNodes, jNodes, kNodes, nVari, nVars
  INTEGER :: ijkCell, errorFlag
  
  REAL(RFREAL) :: pi
  REAL(RFREAL), POINTER :: cofg(:,:)
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: cvSA
  
!******************************************************************************

  global => regionsS1(1)%global
  
  CALL RegisterFunction( global, 'RVAV_ProudmanCulickSolution',&
  'RVAV_ProudmanCulickSolution.F90' )
  
  pi =  global%pi

! generate Input DataStream for Proudman Culick

  DO iReg=1, global%nRegions

    iLev =  regionsS1(iReg)%currLevel
    grid => regionsS1(iReg)%levels(iLev)%grid
    mixt => regionsS1(iReg)%levels(iLev)%mixt

    cofg => grid%cofg
    
! - get cell dimensions

    CALL RFLO_GetDimensDummy( regionsS1(iReg),iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( regionsS1(iReg),iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)   
     
! - get physical dimensions
      
    CALL RFLO_GetDimensPhys( regionsS1(iReg),iLev,ipcbeg,ipcend, &
                             jpcbeg,jpcend,kpcbeg,kpcend )
    CALL RFLO_GetCellOffset( regionsS1(iReg),iLev,iCOff,ijCOff )

! - allocate work array
              
    ALLOCATE( cvSA(3,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
    
! - initialize variables 
           
    cvSA = 0.0_RFREAL
    ipcbeg = 1
    ipcend = 1
    nVari  = 3 

    DO k = kpcbeg, kpcend
    DO j = jpcbeg, jpcend 
    DO i = ipcbeg, ipcend 
      ijkCell = IndIJK(i,j,k,iCOff,ijCOff)
	
      cvSA(1,ijkCell) = SIN(0.5_RFREAL*pi*cofg(YCOORD,ijkCell)/heightC0)
      cvSA(2,ijkCell) = COS(0.5_RFREAL*pi*cofg(YCOORD,ijkCell)/heightC0)
      cvSA(3,ijkCell) = pressureHeadEnd
    END DO ! i 
    END DO ! j 
    END DO ! k

! Write stream2 data

    iNodes = ipcend-ipcbeg+1
    jNodes = jpcend-jpcbeg+1
    kNodes = kpcend-kpcbeg+1
    nVars  = ZCOORD+nVari

    WRITE(IF_RVAV_FILE_S2,*,err=10) iNodes,jNodes,kNodes,nVars 

    WRITE(IF_RVAV_FILE_S2,*,err=10) &
    (((cofg(XCOORD,IndIJK(i,j,k,iCOff,ijCOff)), i=ipcbeg,ipcend), &
                                                j=jpcbeg,jpcend), &
                                                k=kpcbeg,kpcend), &
    (((cofg(YCOORD,IndIJK(i,j,k,iCOff,ijCOff)), i=ipcbeg,ipcend), &
                                                j=jpcbeg,jpcend), &
                                                k=kpcbeg,kpcend), &
    (((cofg(ZCOORD,IndIJK(i,j,k,iCOff,ijCOff)), i=ipcbeg,ipcend), &
                                                j=jpcbeg,jpcend), &
                                                k=kpcbeg,kpcend), &
    ((((cvSA(iVars,IndIJK(i,j,k,iCOff,ijCOff)), i=ipcbeg,ipcend), &
                                                j=jpcbeg,jpcend), &
                                                k=kpcbeg,kpcend), &
                                                iVars=1,nVari)
 
    DEALLOCATE( cvSA,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )   
      
  END DO ! iReg  

  GOTO 999
  
10  CONTINUE
    CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )
      
! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )
  
END SUBROUTINE RVAV_ProudmanCulickSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ProudmanCulickSolution.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:46:35  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.4  2002/10/21 18:56:05  f-najjar
! Bug fix for Character Length of fname to be consistent with RVAV_ModInterfaces
!
! Revision 1.3  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.2  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.1  2002/07/31 02:36:43  f-najjar
! Initial Import
!
!
!******************************************************************************







