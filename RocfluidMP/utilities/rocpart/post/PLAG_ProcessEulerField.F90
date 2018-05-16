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
! Purpose: process Lagragian particle field onto Eulerian-based field data.
!
! Description: none.
!
! Input: regions  = pointer to all regions
!        iReg     = region number
!        nPclsSum = sum of number of particles in all regions
!
! Output: to file.
!
! Notes: 
!   1.   The output is collected in one file, but the regions are processed
!        separately to save memory.
!   2.   The routine assumes that the superparticle loading is fixed for
!        all the particles.
!   3.   Need to fix when particles breakup in the nozzle region.
!   4.   No data is computed in dummy cells, copying is performed for 
!        Tecplot plotting purpose.
!   5.   Raw cell-based data is saved independently for further 
!        post-processing.
!
! ******************************************************************************
!
! $Id: PLAG_ProcessEulerField.F90,v 1.3 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_ProcessEulerField( regions, iReg, nPclsSum )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModGrid,       ONLY : t_grid
  USE ModPartLag,    ONLY : t_plag
  USE PLAG_ModInterfacesPost, ONLY : RFLO_CopyGeometryDummy,   &
                                     RFLO_GenerateCoarseGrids, &
				     RFLO_GetCellOffset,       &
				     RFLO_GetDimensPhysNodes,  &
                                     RFLO_GetDimensDummy,      &
                                     RFLO_GetDimensDummyNodes, &
				     RFLO_GetNodeOffset,       &
				     RFLO_GetCellOffset,       &
				     RFLO_GetDimensPhys,       &
                                     RFLO_ReadGridRegion, Aver

  USE ModParameters
  USE PLAG_ModParameters
  
  IMPLICIT NONE

#include "Indexing.h"

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: iReg, nPclsSum

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================  
! Locals
! ==============================================================================   

  CHARACTER(CHRLEN)   :: RCSIdentString
  CHARACTER(CHRLEN+4) :: fname,fnameTec
  CHARACTER(256)      :: varStr,varStrTec
  CHARACTER(16)       :: compStr
  
  INTEGER, PARAMETER :: IF_PLOT_PLAGEUL = IF_PLOT +200
  INTEGER, PARAMETER :: IF_PLOT_PLAGEUL_TEC = IF_PLOT_PLAGEUL +1

  INTEGER :: iLev, ipc, jpc, kpc, ibc, iec, ibn, ien
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: errorFlag
  INTEGER :: i,j,k, ijkN
  INTEGER :: cell(8),node(8)
  INTEGER :: ic,iCont,iPcl
  INTEGER :: ijkCell,ijkDum
  INTEGER :: nCont,nCvEulerTot,nDensTot,nDiamTot,nMassTot,nPlagCell
  INTEGER, POINTER,      DIMENSION(:)   :: cvMass
  INTEGER, POINTER,      DIMENSION(:,:) :: aiv

  REAL(RFREAL) :: currentTime,diamMicron,massL,massSqrL,nPlagCellInv
  REAL(RFREAL) :: densNum,densNumSqr,diamL,diamL3,diamL4
  REAL(RFREAL) :: xMom,yMom,zMom,temp
  REAL(RFREAL) :: xMomSqr,yMomSqr,zMomSqr,tempSqr
  REAL(RFREAL) :: xcell,ycell,zcell
  REAL(RFREAL) :: u,v,w,uSqr,vSqr,wSqr
  REAL(RFREAL), ALLOCATABLE, DIMENSION(:) :: comp,compSqr
  REAL(RFREAL), POINTER, DIMENSION(:,:)   :: arv,cv,dv
  REAL(RFREAL), POINTER, DIMENSION(:,:)   :: cvEuler,cvSqrEuler,diam,mass,&
                                             massSqr,numDens
  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid)  , POINTER :: grid
  TYPE(t_plag)  , POINTER :: pPlag
  
!*******************************************************************************

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ProcessEulerField.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_ProcessEulerField',&
  'PLAG_ProcessEulerField.F90' )

  IF (.NOT. global%plagUsed) GOTO 999

! ******************************************************************************
! Obtain region grid size
! ******************************************************************************

  DO iLev=2,regions(iReg)%nGridLevels
    ipc = regions(iReg)%levels(iLev-1)%grid%ipc
    jpc = regions(iReg)%levels(iLev-1)%grid%jpc
    kpc = regions(iReg)%levels(iLev-1)%grid%kpc
    regions(iReg)%levels(iLev)%grid%ipc = ipc/2
    regions(iReg)%levels(iLev)%grid%jpc = jpc/2
    regions(iReg)%levels(iLev)%grid%kpc = kpc/2
  ENDDO ! iLev

! ******************************************************************************
! Allocate memory for the grid (all grid levels)
! ******************************************************************************

  DO iLev=1,regions(iReg)%nGridLevels  
    CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
    ibn  =  IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien  =  IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

    grid => regions(iReg)%levels(iLev)%grid
    ALLOCATE( grid%xyz(3,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ENDDO ! iLev

! ******************************************************************************
! Set current grid level and get dummy cell dimensions
! ******************************************************************************

  iLev =  regions(iReg)%currLevel
  grid => regions(iReg)%levels(iLev)%grid

  CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
  ibn = IndIJK(ipnbeg,jpnbeg,kpnbeg,iNOff,ijNOff)
  ien = IndIJK(ipnend,jpnend,kpnend,iNOff,ijNOff)

  CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

! ******************************************************************************
! Save currentTime as global value is clobbered by RFLO_ReadGridRegion
! ******************************************************************************

  currentTime = global%currentTime 

! ******************************************************************************
! Read grid
! ******************************************************************************

  CALL RFLO_ReadGridRegion( iReg,regions )
  CALL RFLO_GenerateCoarseGrids( regions(iReg) )
  CALL RFLO_CopyGeometryDummy( regions(iReg) )

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************
 
  pPlag => regions(iReg)%levels(iLev)%plag

  cv  => pPlag%cv
  dv  => pPlag%dv
  aiv => pPlag%aiv
  arv => pPlag%arv

  cvMass => pPlag%cvPlagMass
  
  nCont = regions(iReg)%plagInput%nCont

  nCvEulerTot = 4
  nDensTot    = 2          ! Includes square of number density 
  nDiamTot    = 2          ! Includes 3rd and 4th power
  nMassTot    = nCont+1       

! ******************************************************************************
! Allocate local arrays
! ******************************************************************************

  ALLOCATE(comp(nCont),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE(compSqr(nCont),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )  

  ALLOCATE(cvEuler(nCvEulerTot,ibc:iec),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE(cvSqrEuler(nCvEulerTot,ibc:iec),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  
  ALLOCATE(diam(nDiamTot,ibc:iec),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE(mass(nMassTot,ibc:iec),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE(massSqr(nMassTot,ibc:iec),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE(numDens(nDensTot,ibc:iec),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ******************************************************************************
! Initialize local arrays
! ******************************************************************************

  DO ic = ibc, iec
    diam(:,ic)       = 0.0_RFREAL
    mass(:,ic)       = 0.0_RFREAL
    massSqr(:,ic)    = 0.0_RFREAL
    cvEuler(:,ic)    = 0.0_RFREAL
    cvSqrEuler(:,ic) = 0.0_RFREAL 
    numDens(:,ic)    = 0.0_RFREAL
  END DO

  comp(:)    = 0.0_RFREAL
  compSqr(:) = 0.0_RFREAL

! ******************************************************************************
! Compute Eulerian-based field data
! ******************************************************************************

  DO iPcl = 1, pPlag%nPcls
    ic = aiv(AIV_PLAG_ICELLS,iPcl)
    diamMicron = dv(DV_PLAG_DIAM,iPcl)*1.0E+06_RFREAL 
    massL = SUM(cv(cvMass(:),iPcl))  

    numDens(1,ic) = numDens(1,ic) +1.0_RFREAL
    numDens(2,ic) = numDens(2,ic) +numDens(1,ic)**2.0_RFREAL
    
    DO iCont = 1, nCont
      mass(iCont,ic)    = mass(iCont,ic)    +cv(cvMass(iCont),iPcl)
      massSqr(iCont,ic) = massSqr(iCont,ic) +cv(cvMass(iCont),iPcl)**2
    END DO ! iCont
    
    mass(nCont+1,ic)    = mass(nCont+1,ic) +massL
    massSqr(nCont+1,ic) = mass(nCont+1,ic) +massL**2.0_RFREAL  

    diam(1,ic) = diam(1,ic) +diamMicron**3.0_RFREAL
    diam(2,ic) = diam(2,ic) +diamMicron**4.0_RFREAL

    cvEuler(1,ic) = cvEuler(1,ic) +cv(CV_PLAG_XMOM,iPcl)
    cvEuler(2,ic) = cvEuler(2,ic) +cv(CV_PLAG_YMOM,iPcl)
    cvEuler(3,ic) = cvEuler(3,ic) +cv(CV_PLAG_ZMOM,iPcl)   
    cvEuler(4,ic) = cvEuler(4,ic) +dv(DV_PLAG_TEMP,iPcl)

    cvSqrEuler(1,ic) = cvSqrEuler(1,ic) +cv(CV_PLAG_XMOM,iPcl)**2.0_RFREAL
    cvSqrEuler(2,ic) = cvSqrEuler(2,ic) +cv(CV_PLAG_YMOM,iPcl)**2.0_RFREAL
    cvSqrEuler(3,ic) = cvSqrEuler(3,ic) +cv(CV_PLAG_ZMOM,iPcl)**2.0_RFREAL   
    cvSqrEuler(4,ic) = cvSqrEuler(4,ic) +dv(DV_PLAG_TEMP,iPcl)**2.0_RFREAL 
  END DO ! iPcl

! ******************************************************************************
! Compute Eulerian-based derived field data
!   Scale field by cell-based number of particles
! ******************************************************************************

  DO ic = ibc,iec
    nPlagCell = numDens(1,ic)
   
    IF ( nPlagCell > 0 ) THEN
      nPlagCellInv = 1.0/REAL(nPlagCell,KIND=RFREAL)
      
      diam(:,ic)       = diam(:,ic)       *nPlagCellInv
      numDens(2,ic)    = numDens(2,ic)    *nPlagCellInv
      mass(:,ic)       = mass(:,ic)       *nPlagCellInv
      massSqr(:,ic)    = massSqr(:,ic)    *nPlagCellInv
      diam(:,ic)       = diam(:,ic)       *nPlagCellInv
      cvEuler(:,ic)    = cvEuler(:,ic)    *nPlagCellInv
      cvSqrEuler(:,ic) = cvSqrEuler(:,ic) *nPlagCellInv
     ENDIF ! numPlagCell 
   END DO ! ic  

   PRINT*,' PLAG_ProcessEulerField: iReg numDenSum,nDimPlag = ', &
    iReg,SUM(numDens(1,:)),pPlag%nPcls

! ******************************************************************************
! Copy field to dummy cells
! ******************************************************************************

! ==============================================================================
! Left boundary
! ==============================================================================

  i = ipcbeg
  DO k=kpcbeg,kpcend
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j
  ENDDO ! k

! ------------------------------------------------------------------------------
! Edges
! ------------------------------------------------------------------------------

  i = ipcbeg
  
  j = jpcbeg
  DO k=kpcbeg,kpcend
    ijkDum  = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k
  
  j = jpcend
  DO k=kpcbeg,kpcend
    ijkDum  = IndIJK(i-1,j+1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

  k = kpcbeg
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j

  k = kpcend
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i-1,j  ,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j 

! ------------------------------------------------------------------------------
! Corners
! ------------------------------------------------------------------------------

  i = ipcbeg
  j = jpcbeg
  k = kpcbeg

    ijkDum  = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  
  i = ipcbeg
  j = jpcend
  k = kpcbeg

    ijkDum  = IndIJK(i-1,j+1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  
  i = ipcbeg
  j = jpcbeg
  k = kpcend

    ijkDum  = IndIJK(i-1,j-1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  i = ipcbeg
  j = jpcend
  k = kpcend

    ijkDum  = IndIJK(i-1,j+1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

! ==============================================================================
! Right boundary
! ==============================================================================

  i = ipcend
  DO k=kpcbeg,kpcend
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i+1,j  ,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j
  ENDDO ! k

! ------------------------------------------------------------------------------
! Edges
! ------------------------------------------------------------------------------

  i = ipcend
  j = jpcbeg
  DO k=kpcbeg,kpcend
    ijkDum  = IndIJK(i+1,j-1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k
  
  i = ipcend
  j = jpcend
  DO k=kpcbeg,kpcend
    ijkDum  = IndIJK(i+1,j+1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

  i = ipcend
  k = kpcbeg
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i+1,j  ,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j

  i = ipcend
  k = kpcend
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i+1,j  ,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j 

! ------------------------------------------------------------------------------
! Corners
! ------------------------------------------------------------------------------

  i = ipcend
  j = jpcbeg
  k = kpcbeg

    ijkDum  = IndIJK(i+1,j-1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  
  i = ipcend
  j = jpcend
  k = kpcbeg

    ijkDum  = IndIJK(i+1,j+1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  
  i = ipcend
  j = jpcbeg
  k = kpcend

    ijkDum  = IndIJK(i+1,j-1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  i = ipcend
  j = jpcend
  k = kpcend

    ijkDum  = IndIJK(i+1,j-1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

! ==============================================================================
! Bottom Boundary
! ==============================================================================

  j = jpcbeg
  DO k=kpcbeg,kpcend
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i ,j-1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! i
  ENDDO ! k

! ------------------------------------------------------------------------------
! Edges
! ------------------------------------------------------------------------------

  j = jpcbeg
  
  i = ipcbeg
  DO k=kpcbeg,kpcend
    ijkDum  = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

  i = ipcend
  DO k=kpcbeg,kpcend
    ijkDum  = IndIJK(i+1,j-1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

  k = kpcbeg
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

  k = kpcend
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i  ,j-1,k+1 ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

! ------------------------------------------------------------------------------
! Corners
! ------------------------------------------------------------------------------

  j = jpcbeg
  
  i = ipcbeg
  k = kpcbeg
    ijkDum  = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  i = ipcend
  k = kpcend
    ijkDum  = IndIJK(i+1,j-1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  k = kpcbeg
  i = ipcbeg
    ijkDum  = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  k = kpcend
  i = ipcend
    ijkDum  = IndIJK(i+1,j-1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

! ==============================================================================
! Top Boundary
! ==============================================================================

  j = jpcend
  DO k=kpcbeg,kpcend
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i ,j+1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! i
  ENDDO ! k

! ------------------------------------------------------------------------------
! Edges
! ------------------------------------------------------------------------------

  j = jpcend
  
  i = ipcbeg
  DO k=kpcbeg,kpcend
    ijkDum  = IndIJK(i-1,j+1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

  i = ipcend
  DO k=kpcbeg,kpcend
    ijkDum  = IndIJK(i+1,j+1,k  ,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

  k = kpcbeg
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i  ,j+1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

  k = kpcend
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i  ,j+1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! k

! ------------------------------------------------------------------------------
! Corners
! ------------------------------------------------------------------------------

  j = jpcend
  
  i = ipcbeg
  k = kpcbeg
    ijkDum  = IndIJK(i-1,j+1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  i = ipcend
  k = kpcend
    ijkDum  = IndIJK(i+1,j+1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  k = kpcbeg
  i = ipcbeg
    ijkDum  = IndIJK(i-1,j+1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  k = kpcend
  i = ipcend
    ijkDum  = IndIJK(i+1,j+1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

! ==============================================================================
! Front Boundary
! ==============================================================================

  k = kpcbeg
  DO j=jpcbeg,jpcend
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i ,j  ,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! i
  ENDDO ! j

! ------------------------------------------------------------------------------
! Edges
! ------------------------------------------------------------------------------
  
  k = kpcbeg
  j = jpcbeg
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! i

  j = jpcend
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i ,j+1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! i

  i = ipcbeg
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j

  i = ipcend
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i+1,j  ,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j

! ------------------------------------------------------------------------------
! Corners
! ------------------------------------------------------------------------------
  
  k = kpcbeg
  
  j = jpcbeg
  i = ipcbeg
    ijkDum  = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  
  j = jpcbeg
  i = ipcend
    ijkDum  = IndIJK(i+1,j-1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  j = jpcend
  i = ipcbeg
    ijkDum  = IndIJK(i-1,j+1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  j = jpcend
  i = ipcend
    ijkDum  = IndIJK(i+1,j+1,k-1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

! ==============================================================================
! Back Boundary
! ==============================================================================

  k = kpcend
  DO j=jpcbeg,jpcend
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i ,j  ,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! i
  ENDDO ! j

! ------------------------------------------------------------------------------
! Edges
! ------------------------------------------------------------------------------

  k = kpcend
  j = jpcbeg
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i ,j-1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! i
  
  j = jpcend
  DO i=ipcbeg,ipcend
    ijkDum  = IndIJK(i ,j+1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! i
  
  i = ipcbeg
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i-1,j  ,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j
  
  i = ipcend
  DO j=jpcbeg,jpcend
    ijkDum  = IndIJK(i+1,j  ,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  ENDDO ! j     

! ------------------------------------------------------------------------------
! Corners
! ------------------------------------------------------------------------------

  k = kpcend
  
  j = jpcbeg
  i = ipcbeg
    ijkDum  = IndIJK(i-1,j-1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  j = jpcbeg
  i = ipcend
    ijkDum  = IndIJK(i+1,j-1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 
  
  j = jpcend
  i = ipcbeg
    ijkDum  = IndIJK(i-1,j+1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell) 

  j = jpcend
  i = ipcend
    ijkDum  = IndIJK(i+1,j+1,k+1,iCOff,ijCOff)
    ijkCell = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)

    diam(:,ijkDum)       = diam(:,ijkCell)     
    numDens(:,ijkDum)    = numDens(:,ijkCell)  
    mass(:,ijkDum)       = mass(:,ijkCell)     
    massSqr(:,ijkDum)    = massSqr(:,ijkCell)  
    diam(:,ijkDum)       = diam(:,ijkCell)     
    cvEuler(:,ijkDum)    = cvEuler(:,ijkCell)   
    cvSqrEuler(:,ijkDum) = cvSqrEuler(:,ijkCell)  
 
! ******************************************************************************
! Write data to Tecplot file
! ******************************************************************************

! ==============================================================================
!  Open file and write the header
! ==============================================================================

   IF (iReg == 1) THEN
     WRITE(fnameTec,'(A,ES11.5,A)') &
     TRIM(global%casename)//'.plag_eulerVarsTec_',currentTime,'.dat'
     OPEN(IF_PLOT_PLAGEUL_TEC,FILE=fnameTec,STATUS='unknown',FORM='formatted',&
	  iostat=errorFlag)
     global%error = errorFlag
     IF (global%error /= 0) &
       CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fnameTec )     
     
     WRITE(IF_PLOT_PLAGEUL_TEC,1005,err=10) &
       TRIM(global%casename),global%currentTime,nPclsSum

     varStrTec = 'x y z nDens nDensSqr diam diam3 diam4 mass u v w temp massSqr uSqr vSqr wSqr tempSqr'

     DO iCont=1,nCont
       SELECT CASE (iCont)
         CASE ( 0: 9)
          WRITE(compStr,'(A,I1)') 'comp_', iCont
         CASE (10:99)
          WRITE(compStr,'(A,I2)') 'comp_', iCont
         CASE DEFAULT
          WRITE(compStr,'(A)') 'comp_?'
       END SELECT ! iCont
       varStrTec = TRIM(varStrTec)//' '//TRIM(compStr)
     ENDDO ! iCont

     DO iCont=1,nCont
       SELECT CASE (iCont)
         CASE ( 0: 9)
          WRITE(compStr,'(A,I1)') 'compSqr_', iCont
         CASE (10:99)
          WRITE(compStr,'(A,I2)') 'compSqr_', iCont
         CASE DEFAULT
          WRITE(compStr,'(A)') 'compSqr_?'
       END SELECT ! iCont
       varStrTec = TRIM(varStrTec)//' '//TRIM(compStr)
     ENDDO ! iCont     

     WRITE(IF_PLOT_PLAGEUL_TEC,1010,err=10) TRIM(varStrTec)
   ENDIF   ! iReg=1

! ==============================================================================
!  Write zone header for plotting data based on nodes
! ==============================================================================

   WRITE(IF_PLOT_PLAGEUL_TEC,1015) iReg,ipnend-ipnbeg+1,jpnend-jpnbeg+1,&
                                   kpnend-kpnbeg+1

! ==============================================================================
! Write plotting data
! ==============================================================================

  DO k=kpnbeg,kpnend
  DO j=jpnbeg,jpnend
  DO i=ipnbeg,ipnend
    ijkN = IndIJK(i,j,k,iNOff,ijNOff)

    cell(1) = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
    cell(2) = IndIJK(i+1,j  ,k  ,iCOff,ijCOff)
    cell(3) = IndIJK(i  ,j+1,k  ,iCOff,ijCOff)
    cell(4) = IndIJK(i+1,j+1,k  ,iCOff,ijCOff)
    cell(5) = IndIJK(i  ,j  ,k+1,iCOff,ijCOff)
    cell(6) = IndIJK(i+1,j  ,k+1,iCOff,ijCOff)
    cell(7) = IndIJK(i  ,j+1,k+1,iCOff,ijCOff)
    cell(8) = IndIJK(i+1,j+1,k+1,iCOff,ijCOff)

    diamL3 = Aver(cell,1,diam)
    diamL4 = Aver(cell,2,diam)

    IF( diamL3 > 0.0_RFREAL) THEN
      diamL = diamL4/diamL3
    ELSE
      diamL = 0.0_RFREAL
    ENDIF ! diamL3

    densNum    = Aver(cell,1,numDens)
    densNumSqr = Aver(cell,2,numDens)

    massL      = Aver(cell,nCont+1,mass)
    massSqrL   = Aver(cell,nCont+1,massSqr)

    xMom    = Aver(cell,1,cvEuler)
    yMom    = Aver(cell,2,cvEuler)
    zMom    = Aver(cell,3,cvEuler)
    temp    = Aver(cell,4,cvEuler)

    xMomSqr    = Aver(cell,1,cvSqrEuler)
    yMomSqr    = Aver(cell,2,cvSqrEuler)
    zMomSqr    = Aver(cell,3,cvSqrEuler)
    tempSqr    = Aver(cell,4,cvSqrEuler)

    IF ( massL > 0.0_RFREAL ) THEN
      u = xMom/massL
      v = yMom/massL
      w = zMom/massL
    ELSE
      u = 0.0_RFREAL; v = 0.0_RFREAL; w = 0.0_RFREAL;
    END IF ! massL
    
    IF ( massSqrL > 0.0_RFREAL ) THEN
      uSqr = xMomSqr/massSqrL
      vSqr = yMomSqr/massSqrL
      wSqr = zMomSqr/massSqrL
    ELSE
      uSqr = 0.0_RFREAL; vSqr = 0.0_RFREAL; wSqr = 0.0_RFREAL;
    END IF ! massSqrL
    
    DO iCont = 1, nCont
      comp(iCont)    = Aver(cell,iCont,mass)
      compSqr(iCont) = Aver(cell,iCont,mass)
      
      IF ( massL > 0.0_RFREAL ) THEN
       comp(iCont) = comp(iCont)/massL
      ELSE
       comp(iCont) = 0.0_RFREAL
      END IF ! massL
      
      IF ( massSqrL > 0.0_RFREAL) THEN
        compSqr(iCont) = compSqr(iCont)/massSqrL
      ELSE
        compSqr(iCont) = 0.0_RFREAL
      END IF ! massSqrL
    END DO ! iCont

!!    ic = cell(1)
!    ic = ijkN
!        
!
!    diamL3 = diam(1,ic)
!    diamL4 = diam(2,ic)
!    IF( diamL3 > 0.0_RFREAL) THEN
!      diamL = diamL4/diamL3
!    ELSE
!      diamL = 0.0_RFREAL
!    ENDIF ! diamL3
!
!    densNum    = numDens(1,ic)
!    densNumSqr = numDens(2,ic)
!    massL      = mass(nCont+1,ic)
!    massSqrL   = massSqr(nCont+1,ic)
!
!    xMom  = cvEuler(1,ic)
!    yMom  = cvEuler(2,ic)
!    zmom  = cvEuler(3,ic)
!    temp  = cvEuler(4,ic)
!
!    xMomSqr = cvSqrEuler(1,ic)
!    yMomSqr = cvSqrEuler(2,ic)
!    zMomSqr = cvSqrEuler(3,ic)
!    tempSqr = cvSqrEuler(4,ic)    
!
   
    WRITE(IF_PLOT_PLAGEUL_TEC,1020,err=10) grid%xyz(XCOORD,ijkN),            &
                                           grid%xyz(YCOORD,ijkN),            &
                                           grid%xyz(ZCOORD,ijkN),            &
				           densNum,densNumSqr,               &
                                           diamL,diamL3,diamL4,              &
				           massL,u,v,w,temp,                 &
				           massSqrL,uSqr,vSqr,wSqr,tempSqr,  &
                                          (comp(iCont),iCont=1,nCont),       &
                                          (compSqr(iCont),iCont=1,nCont)

  ENDDO ! i
  ENDDO ! j
  ENDDO ! k

! ******************************************************************************
! Write raw data to formatted file
! ******************************************************************************

! ==============================================================================
!  Open file and write the header
! ==============================================================================

   IF (iReg == 1) THEN
     WRITE(fname,'(A,ES11.5,A)') &
     TRIM(global%casename)//'.plag_eulerVars_',currentTime,'.dat'
     OPEN(IF_PLOT_PLAGEUL,FILE=fname,STATUS='unknown',FORM='formatted',&
	  iostat=errorFlag)
     global%error = errorFlag
     IF (global%error /= 0) &
       CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

     WRITE(IF_PLOT_PLAGEUL,1005,err=10) &
       TRIM(global%casename),global%currentTime,nPclsSum

     varStr = 'x y z nDens nDensSqr diam diam3 diam4 mass xmom ymom zmom temp massSqr xmomSqr ymomSqr zmomSqr tempSqr'

     DO iCont=1,nCont
       SELECT CASE (iCont)
         CASE ( 0: 9)
          WRITE(compStr,'(A,I1)') 'mass_', iCont
         CASE (10:99)
          WRITE(compStr,'(A,I2)') 'mass_', iCont
         CASE DEFAULT
          WRITE(compStr,'(A)') 'mass_?'
       END SELECT ! iCont
       varStr = TRIM(varStr)//' '//TRIM(compStr)
     ENDDO ! iCont

     DO iCont=1,nCont
       SELECT CASE (iCont)
         CASE ( 0: 9)
          WRITE(compStr,'(A,I1)') 'massSqr_', iCont
         CASE (10:99)
          WRITE(compStr,'(A,I2)') 'massSqr_', iCont
         CASE DEFAULT
          WRITE(compStr,'(A)') 'massSqr_?'
       END SELECT ! iCont
       varStr = TRIM(varStr)//' '//TRIM(compStr)
     ENDDO ! iCont     

     WRITE(IF_PLOT_PLAGEUL,1010,err=10) TRIM(varStr)
   ENDIF   ! iReg=1

! ==============================================================================
!  Write zone header for raw data based on cells
! ==============================================================================

   WRITE(IF_PLOT_PLAGEUL,1015) iReg,ipcend-ipcbeg+1,jpcend-jpcbeg+1,&
                               kpcend-kpcbeg+1

! ==============================================================================
! Write cell-based data
! ==============================================================================

  DO k=kpcbeg,kpcend
  DO j=jpcbeg,jpcend
  DO i=ipcbeg,ipcend
    ijkN = IndIJK(i,j,k,iNOff,ijNOff)

    node(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
    node(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
    node(3) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
    node(4) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
    node(5) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
    node(6) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
    node(7) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)
    node(8) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
    
    xcell = Aver(node,1,grid%xyz)
    ycell = Aver(node,2,grid%xyz)
    zcell = Aver(node,3,grid%xyz)
    
    cell(1) = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
    ic      = cell(1)
        
    diamL  = 0.0_RFREAL 
    diamL3 = diam(1,ic)
    diamL4 = diam(2,ic)
    IF( diamL3 > 0.0_RFREAL) THEN
      diamL = diamL4/diamL3
    ENDIF ! diamL3

    densNum    = numDens(1,ic)
    densNumSqr = numDens(2,ic)
    massL      = mass(nCont+1,ic)
    massSqrL   = massSqr(nCont+1,ic)

    xMom  = cvEuler(1,ic)
    yMom  = cvEuler(2,ic)
    zmom  = cvEuler(3,ic)
    temp  = cvEuler(4,ic)

    xMomSqr = cvSqrEuler(1,ic)
    yMomSqr = cvSqrEuler(2,ic)
    zMomSqr = cvSqrEuler(3,ic)
    tempSqr = cvSqrEuler(4,ic)

    WRITE(IF_PLOT_PLAGEUL,1020,err=10) xcell,ycell,zcell,               &
				       densNum,densNumSqr,              &
                                       diamL,diamL3,diamL4,             &
				       massL,xMom,yMom,zMom,temp,       &
				       massSqrL,xmomSqr,yMomSqr,        &
				       zMomSqr,tempSqr,                 &
                                      (mass(iCont,ic),iCont=1,nCont),   &
                                      (massSqr(iCont,ic),iCont=1,nCont)

  ENDDO ! i
  ENDDO ! j
  ENDDO ! k

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DO iLev=1,regions(iReg)%nGridLevels
    grid => regions(iReg)%levels(iLev)%grid
    DEALLOCATE( grid%xyz,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
  ENDDO ! iLev

  iLev = regions(iReg)%currLevel

! ==============================================================================
!   Dealocate local arrays
! ==============================================================================

  DEALLOCATE(comp,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE(compSqr,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE(cvEuler,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE(cvSqrEuler,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE(diam,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE(mass,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE(massSqr,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE(numDens,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! ******************************************************************************
! Finalize
! ******************************************************************************
 
  IF (iReg == global%nRegions) THEN
    CLOSE(IF_PLOT_PLAGEUL,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )

    CLOSE(IF_PLOT_PLAGEUL_TEC,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fnameTec )
  ENDIF

  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,fname )

999  CONTINUE
  CALL DeregisterFunction( global )

! ******************************************************************************
! Formats
! ******************************************************************************

1005 FORMAT('TITLE="',A,'. Time: ',ES11.5,'. nPclsSum: ',I8,'."')
1010 FORMAT('VARIABLES= ',A)
1015 FORMAT('ZONE T="',I5.5,'", I=',I6,', J=',I6,', K=',I6,', F=POINT')
1020 FORMAT(999(1X,ES13.6))

END SUBROUTINE PLAG_ProcessEulerField

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ProcessEulerField.F90,v $
! Revision 1.3  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/17 22:15:23  fnajjar
! Initial import
!
!******************************************************************************







