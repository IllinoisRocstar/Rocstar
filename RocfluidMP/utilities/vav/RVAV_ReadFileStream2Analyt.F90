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
! Purpose: read File for Stream 2 Data based on Analytical Solutions.
!
! Description: none.
!
! Input: 
!
! Output: Memory location with data for Stream 2
!
! ISSUE: Where do you read the Analytical and Experimental Data on Nodes or Cells
! 
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_ReadFileStream2Analyt.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ReadFileStream2Analyt ( global, regionsS2 )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE ModMPI
  USE ModParameters
  USE RVAV_ModParameters
  USE RVAV_ModGlobal
  USE RVAV_ModDataStruct
  IMPLICIT NONE

#include "Indexing.h"
  
! ... parameter variables
  TYPE (t_region), POINTER :: regionsS2(:)
  TYPE(t_global)  , POINTER  :: global
  
! ... loop variables
  INTEGER :: iReg, iLev

! ... local variables
  CHARACTER(CHRLEN) :: msg, fname
  
  INTEGER :: i, j, k
  INTEGER :: iNodes, jNodes, kNodes, nVars, iVars
  INTEGER :: ipc, jpc, kpc, ibc, iec, ibn, ien
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: nRegionsS2
  INTEGER :: ijkN0, ijkC0
  INTEGER :: nGridLevels, ivar, errorFlag

  REAL(RFREAL), DIMENSION(3)               :: xyzMin, xyzMax
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE  ::  cvMin,  cvMax
  REAL(RFREAL), POINTER                    ::  cv(:,:),xyz(:,:)  
  
  TYPE(t_grid)    , POINTER  :: grid
  TYPE(t_mixt)    , POINTER  :: mixt  
    
!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_ReadFileStream2Analyt',&
  'RVAV_ReadFileStream2Analyt.F90' )
	   
  IF ( global%verbLevel/=VERBOSE_NONE ) &	   
    WRITE(STDOUT,'(/,A)') 'Reading grid and solution from Stream 2 - ANALYTICAL...' 

  global%casename = TRIM(globalRVAV%casename)//'_s2'

! open file

  IF ( globalRVAV%fileTypeS2 == FILE_ANALYTICAL ) THEN
    WRITE(fname,'(A,I5.5,A,1PE11.5)') TRIM(global%casename)//'.anl'
  ENDIF ! fileTypeS2

  OPEN(IF_RVAV_FILE_S2,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global, ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! read number of Regions 

  READ(IF_RVAV_FILE_S2,*,err=10,end=10) nRegionsS2
    
  globalRVAV%nRegionsS2 = nRegionsS2 
     
  ALLOCATE( regionsS2(globalRVAV%nRegionsS2),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
    
  DO iReg=1,globalRVAV%nRegionsS2

    WRITE(STDOUT,'(A,I5.5)') '  - region ',iReg

    regionsS2(iReg)%startLevel = global%startLevel
    regionsS2(iReg)%currLevel  = global%startLevel

    iLev =  regionsS2(iReg)%currLevel
    regionsS2(iReg)%nGridLevels = regionsS2(iReg)%startLevel
    nGridLevels = regionsS2(iReg)%nGridLevels
      
    IF ( global%verbLevel/=VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(/,A,I5)') 'regionsS2(iReg)%startLevel=',regionsS2(iReg)%startLevel
      WRITE(STDOUT,'(A,I5)')   'regionsS2(iReg)%currLevel =',regionsS2(iReg)%currLevel
      WRITE(STDOUT,'(A,I5)')   'iLev = ',iLev
    END IF ! verbLevel
      
    ALLOCATE( regionsS2(iReg)%levels(nGridLevels),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )

  ENDDO ! iReg

  DO iReg=1,globalRVAV%nRegionsS2
          
    regionsS2(iReg)%startLevel = global%startLevel
    regionsS2(iReg)%currLevel  = global%startLevel 

    iLev =  regionsS2(iReg)%currLevel
      
    IF ( global%verbLevel/=VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,I5.5)') '  - region ',iReg
      WRITE(STDOUT,'(/,A,I5)') 'regionsS2(iReg)%startLevel=',regionsS2(iReg)%startLevel
      WRITE(STDOUT,'(A,I5)')   'regionsS2(iReg)%currLevel =',regionsS2(iReg)%currLevel
      WRITE(STDOUT,'(A,I5)')   'iLev = ',iLev
    END IF ! verbLevel

    grid => regionsS2(iReg)%levels(iLev)%grid
    mixt => regionsS2(iReg)%levels(iLev)%mixt
             
! - read dimensions

    READ(IF_RVAV_FILE_S2,*,err=10,end=10) iNodes,jNodes,kNodes,nVars 
            
    IF ( global%verbLevel/=VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3(I5,2X))')   'nRegionsS2 = ',nRegionsS2
      WRITE(STDOUT,'(A,3(I5,2X))')   'iNodes,jNodes,kNodes = ',iNodes,jNodes,kNodes
      WRITE(STDOUT,'(A,1(I5,2X))')   'nVars = ',nVars
    END IF ! verbLevel

! - find array extent
      
    ipnbeg = 1
    jpnbeg = 1
    kpnbeg = 1
      
    ipnend = iNodes
    jpnend = jNodes
    kpnend = kNodes
 
    iNOff  = ipnend - ipnbeg + 1
    ijNOff = iNOff*(jpnend-jpnbeg+1)
      
    ibn = IndIJK(ipnbeg,jpnbeg,kpnbeg,iNOff,ijNOff)
    ien = IndIJK(ipnend,jpnend,kpnend,iNOff,ijNOff) 
                              
    globalRVAV%iCOffS2  = iNOff
    globalRVAV%ijCOffS2 = ijNOff
      
! - allocate memory

    ALLOCATE( grid%xyz(3,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )

    ALLOCATE( mixt%cv(nVars-3,ibn:ien) ,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
    
! - initialize variables 

    grid%xyz = 0.0_RFREAL
    mixt%cv  = 0.0_RFREAL

! - create new pointers for reading data 

    xyz => grid%xyz
     cv => mixt%cv
        
! - read stream2 data 

    READ(IF_RVAV_FILE_S2,*,err=10,end=10) &
    (((xyz(XCOORD,IndIJK(i,j,k,iNOff,ijNOff)),    i=ipnbeg,ipnend), &
                                                  j=jpnbeg,jpnend), &
                                                  k=kpnbeg,kpnend), &
    (((xyz(YCOORD,IndIJK(i,j,k,iNOff,ijNOff)),    i=ipnbeg,ipnend), &
                                                  j=jpnbeg,jpnend), &
                                                  k=kpnbeg,kpnend), &
    (((xyz(ZCOORD,IndIJK(i,j,k,iNOff,ijNOff)),    i=ipnbeg,ipnend), &
                                                  j=jpnbeg,jpnend), &
                                                  k=kpnbeg,kpnend), &
    ((((cv(iVars, IndIJK(i,j,k,iNOff,ijNOff)),    i=ipnbeg,ipnend), &
                                                  j=jpnbeg,jpnend), &
                                                  k=kpnbeg,kpnend), &
                                              iVars=1,nVars-ZCOORD)

! - Write Min-Max Values

    IF ( global%verbLevel/=VERBOSE_NONE ) THEN
      
      WRITE(STDOUT,'(/,A,2(I5,2X))') 'ibn-ien ',ibn,ien
      WRITE(STDOUT,'(A,3(I5,2X))')   'ipnbeg,jpnbeg,kpnbeg',ipnbeg,jpnbeg,kpnbeg
      WRITE(STDOUT,'(A,3(I5,2X))')   'ipnend,jpnend,kpnend',ipnend,jpnend,kpnend
      WRITE(STDOUT,'(A,2(I5,2X))')   'iNOff,ijNOff',iNOff,ijNOff
      WRITE(STDOUT,'(A,2(I5,2X))')   'iCOffS2,ijCOffS2',globalRVAV%iCOffS2,globalRVAV%ijCOffS2
        
      xyzMin = +1.0E+30_RFREAL
      xyzMax = -1.0E+30_RFREAL

      DO i = ibn, ien
        xyzMin(1) = MIN(grid%xyz(1,i), xyzMin(1))
        xyzMin(2) = MIN(grid%xyz(2,i), xyzMin(2))
        xyzMin(3) = MIN(grid%xyz(3,i), xyzMin(3))

        xyzMax(1) = MAX(grid%xyz(1,i), xyzMax(1))
        xyzMax(2) = MAX(grid%xyz(2,i), xyzMax(2))
        xyzMax(3) = MAX(grid%xyz(3,i), xyzMax(3))
      END DO ! i

      DO k = kpnbeg, kpnbeg
      DO j = jpnbeg, jpnbeg
      DO i = ipnbeg, ipnend
        ijkN0  = IndIJK(i  ,j,k,iNOff,ijNOff)
        WRITE(STDOUT,'(I5,3(3X,E12.5))') i,grid%xyz(1,ijkN0),grid%xyz(2,ijkN0),grid%xyz(3,ijkN0)
      END DO ! i 
      END DO ! j
      END DO ! k

      WRITE(STDOUT,'(/,A,2(E12.5,2X))') 'Min-Max of X     ', xyzMin(1),xyzMax(1)
      WRITE(STDOUT,'(A,2(E12.5,2X))')   'Min-Max of Y     ', xyzMin(2),xyzMax(2)
      WRITE(STDOUT,'(A,2(E12.5,2X))')   'Min-Max of Z     ', xyzMin(3),xyzMax(3)

      ALLOCATE(cvMin(nVars-ZCOORD),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )
      ALLOCATE(cvMax(nVars-ZCOORD),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ )

      cvMin = +1.0E+30_RFREAL
      cvMax = -1.0E+30_RFREAL

      DO i = ibn, ien
        DO ivar = 1, nVars-ZCOORD
          cvMin(ivar) = MIN(mixt%cv(ivar,i), cvMin(ivar))
          cvMax(ivar) = MAX(mixt%cv(ivar,i), cvMax(ivar))
        END DO ! ivar 
      END DO ! i

        WRITE(STDOUT,'(/,A)') 'Distribution of Stream2 in I-Direction '
        DO k = kpnbeg, kpnbeg
        DO j = jpnbeg, jpnbeg
        DO i = ipnbeg, ipnend
          ijkN0  = IndIJK(i  ,j,k,iNOff,ijNOff)
          WRITE(STDOUT,'(I5,4(3X,1PE15.7))') i,(mixt%cv(iVar,ijkN0),iVar=1,nVars-ZCOORD)
        END DO ! i
        END DO ! j 
        END DO ! k

        WRITE(STDOUT,'(/,A)') 'Distribution of Stream2 in J-Direction '
        DO k = kpnbeg, kpnbeg
        DO j = jpnbeg, jpnend
        DO i = ipnbeg, ipnbeg
          ijkN0  = IndIJK(i  ,j,k,iNOff,ijNOff)
          WRITE(STDOUT,'(I5,4(3X,1PE15.7))') j,(mixt%cv(iVar,ijkN0),iVar=1,nVars-ZCOORD)
        END DO ! i
        END DO ! j
        END DO ! k

        DO ivar = 1, nVars-ZCOORD
          WRITE(STDOUT,'(/,A,I5)') 'ivar = ',ivar
          WRITE(STDOUT,'(/,A,2(E12.5,2X))')'Min-Max of CV ', cvMin(ivar),cvMax(ivar)
        END DO ! ivar

        DEALLOCATE(cvMin,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
        DEALLOCATE(cvMax,stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
	
      END IF ! verbLevel
      
    END DO ! iReg 
  
  GOTO 999
  
10  CONTINUE
  CALL ErrorStop( global, ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) ) 
   
! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )
  
1000 FORMAT('Region ',I5,', grid level= ',I2,'.')

END SUBROUTINE RVAV_readFileStream2Analyt

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ReadFileStream2Analyt.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:46:35  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.3  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.2  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.1  2002/07/16 22:34:16  f-najjar
! Initial Import
!
!
!******************************************************************************







