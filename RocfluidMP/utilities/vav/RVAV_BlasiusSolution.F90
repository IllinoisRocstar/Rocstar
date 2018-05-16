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
! Purpose: Generate Blasius Analytical Solution.
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
! $Id: RVAV_BlasiusSolution.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_BlasiusSolution( fname, regionsS1)

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

  INTEGER, PARAMETER        :: nFine=131073 
  INTEGER                   :: ijkN0,ijkN1,ijkNR0, jTilde
  
  REAL(RFREAL), PARAMETER   :: simCoef=1.905_RFREAL
  REAL(RFREAL)              :: yFin(nFine),varFin(2,nFine),varVal(2),dumvar
  REAL(RFREAL)              :: refLen,refNu,refUvel,axialDist,deltaStar,simCoord

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

  REAL(RFREAL), POINTER :: cofg(:,:)
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: cvSA
  
!******************************************************************************

  global => regionsS1(1)%global
  
  CALL RegisterFunction( global, 'RVAV_BlasiusSolution',&
  'RVAV_BlasiusSolution.F90' )
    
! generates Input DataStream for Blasius

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
              
    ALLOCATE( cvSA(2,ibc:iec) ,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )

! - initialize variables 
           
    cvSA   = 0.0_RFREAL  
    yFin   = 0.0_RFREAL
    varFin = 0.0_RFREAL
    nVari  = 2 
                
! - open file containing analytical Blasius solution
          
    WRITE(fname,'(A,A)') TRIM(globalRVAV%casename)//'_finegrid.anl'
    OPEN(IF_RVAV_FILE_BL,file=fname,form='formatted',status='unknown') 

! - read in reference quantities and similarity data; rho, u stored in varFin(1:2,:)

    READ(IF_RVAV_FILE_BL,*,err=10,end=10) refLen
    READ(IF_RVAV_FILE_BL,*,err=10,end=10) refNu
    READ(IF_RVAV_FILE_BL,*,err=10,end=10) refUvel

    DO j = 1, nFine
      READ(IF_RVAV_FILE_BL,*,err=10,end=10) yFin(j),varFin(1,j),varFin(2,j),dumvar 
    ENDDO ! j

! - normalize solution component rho and u-velocity by their freestream values

    DO j = 1, nFine
      varFin(:,j) = varFin(:,j)/varFin(:,nFine) 
    ENDDO ! j

! - close file

    CLOSE(IF_RVAV_FILE_BL)

!-  interpolate analytical solution to cell centroid

    ijkNR0 = IndIJK(ipcbeg,jpcbeg,kpcbeg,iCOff,ijCOff)
    
    DO i = ipcbeg, ipcend
      ijkN0     = IndIJK(i     ,jpcbeg,kpcbeg,iCOff,ijCOff)
      axialDist = cofg(XCOORD,ijkN0)-cofg(XCOORD,ijkNR0)
      deltaStar = max(simCoef*sqrt(refNu*axialDist/refUvel),1.0E-08_RFREAL)

      jTilde=1 

      DO j = jpcbeg, jpcend 
        ijkN1    = IndIJK(i  ,j   ,kpcbeg,iCOff,ijCOff)
        simCoord = cofg(YCOORD,ijkN1)/deltaStar

! - find the smallest jtilde such that simCoord is strictly lower than yFin(jtilde)

        DO WHILE ((simCoord >= yFin(jTilde)) .AND. (jTilde < nFine))
          jTilde=jTilde+1
        END DO

! - interpolate between values stored in 'Fin' arrays varFin

        CALL LinIntpol(nFine, jTilde, nVari, simCoord, yFin, varFin, varVal)  

! - store the interpolated variables into cvSA
 
         cvSA(1:nVari,ijkN1)=varVal(1:nVari)
       
        END DO ! j
          
    END DO ! i 

! - copy the solution to other k-layers

    DO k = kpcbeg+1, kpcend
    DO j = jpcbeg, jpcend 
    DO i = ipcbeg, ipcend 
      ijkN0    = IndIJK(i  ,j  ,k     ,iCOff,ijCOff)
      ijkNR0   = IndIJK(i  ,j  ,kpcbeg,iCOff,ijCOff)
      cvSA(1:nVari,ijkN0) = cvSA(1:nVari,ijkNR0)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
    
! Write stream2 data

    iNodes = ipcend-ipcbeg+1
    jNodes = jpcend-jpcbeg+1
    kNodes = kpcend-kpcbeg+1
    nVars  = ZCOORD+nVari
          
    WRITE(IF_RVAV_FILE_S2,*,err=10) iNodes,jNodes,kNodes,nVars 


    WRITE(IF_RVAV_FILE_S2,*,err=10) &
    (((cofg(XCOORD,IndIJK(i,j,k,iCOff,ijCOff)),   i=ipcbeg,ipcend), &
                                                  j=jpcbeg,jpcend), &
                                                  k=kpcbeg,kpcend), &
    (((cofg(YCOORD,IndIJK(i,j,k,iCOff,ijCOff)),   i=ipcbeg,ipcend), &
                                                  j=jpcbeg,jpcend), &
                                                  k=kpcbeg,kpcend), &
    (((cofg(ZCOORD,IndIJK(i,j,k,iCOff,ijCOff)),   i=ipcbeg,ipcend), &
                                                  j=jpcbeg,jpcend), &
                                                  k=kpcbeg,kpcend), &
    ((((cvSA(iVars,IndIJK(i,j,k,iCOff,ijCOff)),   i=ipcbeg,ipcend), &
                                                  j=jpcbeg,jpcend), &
                                                  k=kpcbeg,kpcend), &
                                                  iVars=1,nVari)
 
    DEALLOCATE( cvSA,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ )   
      
  END DO ! iReg  
          
  GOTO 999
  
10  CONTINUE
    CALL ErrorStop( global, ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )  
  
! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )
  
!----------------------------------------------------------------

CONTAINS 

  SUBROUTINE LinIntpol( nFine, jTilde, nComp, crsCoord, &
                        finCoord, varFin, varVal)
  
! ... input variables
  INTEGER      :: nFine, jTilde, nComp  
  REAL(RFREAL) :: crsCoord
  REAL(RFREAL) :: finCoord(nFine), varFin(nComp,nFine)
  
! ... output variables
  REAL(RFREAL) :: varVal(nComp)
  
! ... local variables
  INTEGER      :: l
  REAL(RFREAL) :: ratio

! - perform simple linear interpolation between the fine mesh data

  ratio=(crsCoord-finCoord(jTilde-1))/(finCoord(jTilde)-finCoord(jTilde-1))
  DO l=1, nComp
    varVal(l) = varFin(l,jTilde-1) + ratio*(varFin(l,jTilde)-varFin(l,jTilde-1))
  ENDDO

  END SUBROUTINE LinIntpol  

END SUBROUTINE RVAV_BlasiusSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_BlasiusSolution.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:18  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/05/15 02:57:07  jblazek
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
! Revision 1.1  2002/07/31 02:36:00  f-najjar
! Initial Import
!
!
!******************************************************************************







