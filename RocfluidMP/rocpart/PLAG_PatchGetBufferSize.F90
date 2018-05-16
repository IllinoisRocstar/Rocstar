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
! Purpose: obtains exact buffer size for patches. 
!
! Description: none.
!
! Input: region = current region.
!
! Output: regions(iReg)%levels%patch%buffPlag%nBuffSize = buffer size.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_PatchGetBufferSize.F90,v 1.4 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_PatchGetBufferSize( region )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_buffer_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset,   RFLO_GetNodeOffset

  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, idum, iPatch, iPcls, j,k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, ibeg, iC, iCOff, idir, iend, ijCOff, ijNOff,      &
             ijkCPatch, ijkCPlag, ijkDPatch,iLev, inode, iNOff,        &
             jbeg, jC, jdir, jend, jnode, kbeg, kC, kdir, kend, knode, & 
             lbound, n1, n2, nDumCells, nPatches, nPcls, nPclsBuffTot, nOff 
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv

  
  TYPE(t_patch),       POINTER :: pPatch
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_buffer_plag), POINTER :: pBuffPlag
  TYPE(t_global),      POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_PatchGetBufferSize.F90,v $ $Revision: 1.4 $'

  global => region%global
    
  CALL RegisterFunction( global, 'PLAG_PatchGetBufferSize',&
  'PLAG_PatchGetBufferSize.F90' )

! Get dimensions --------------------------------------------------------------

  iLev         = region%currLevel
  nPatches     = region%nPatches
  nDumCells    = region%nDumCells
  
  nPcls        = region%levels(iLev)%plag%nPcls 
  nPclsBuffTot = region%plagInput%nPclsBuffTot
  
! Set pointers ----------------------------------------------------------------
   
  pPlag   => region%levels(iLev)%plag   
  pAiv    => pPlag%aiv

! Loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    pPatch  => region%levels(iLev)%patches(iPatch)

    bcType = pPatch%bcType
    
! - Get patch dimensions ------------------------------------------------------

    lbound = pPatch%lbound

    CALL RFLO_GetPatchIndices( region,pPatch,iLev, &
                               ibeg,iend,jbeg,jend,kbeg,kend )
    CALL RFLO_GetPatchDirection( pPatch,idir,jdir,kdir )
    CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    nOff        = ABS(pPatch%l1end-pPatch%l1beg) + 1

    inode = 0
    jnode = 0
    knode = 0
    IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
      inode = -idir
      jnode = -jdir
      knode = -kdir
    ENDIF ! lbound
  
! - Set pointers and initialize buffer sizes ----------------------------------

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
        (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN 
    
      pBuffPlag => pPatch%bufferPlag     
    
      pBuffPlag%nBuffSize  = 0
    
      pBuffPlag%nSendBuffI = 0
      pBuffPlag%nSendBuffR = 0
          
      pBuffPlag%nRecvBuffI = 0
      pBuffPlag%nRecvBuffR = 0

!- Bypass loop for null number of particles -----------------------------------

      IF ( nPcls == 0 ) GOTO 1999
    
! - Loop over particles -------------------------------------------------------

      DO iPcls = 1, nPcls
        ijkCPlag = pAiv(AIV_PLAG_ICELLS,iPcls)
        iC       = pAiv(AIV_PLAG_INDEXI,iPcls)
        jC       = pAiv(AIV_PLAG_INDEXJ,iPcls)
        kC       = pAiv(AIV_PLAG_INDEXK,iPcls)

! -- Loop over patch and dummy cells ------------------------------------------

        DO idum=1,nDumCells
          DO k=kbeg,kend
            DO j=jbeg,jend
              DO i=ibeg,iend
                ijkCPatch  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir, iCOff,ijCOff)
                ijkDPatch  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
                IF ( ijkCPlag == ijkDPatch ) THEN
                  pBuffPlag%nBuffSize = pBuffPlag%nBuffSize + 1
                END IF ! ijkCPlag
      
              END DO ! i
            END DO ! j
          END DO ! k
        END DO  ! idum 

      END DO ! iPcls

1999  CONTINUE 
    ENDIF  ! bcType         
  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_PatchGetBufferSize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PatchGetBufferSize.F90,v $
! Revision 1.4  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:56  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/04/09 23:14:00  fnajjar
! Moved bypass of null particle size after setting buffer sizes to zero
!
! Revision 1.6  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.5  2003/01/16 20:58:00  f-najjar
! Reincluded missing CVS comments
!
! Revision 1.4  2003/01/16 20:24:32  f-najjar
! Removed iRegionGlobal
!
! Revision 1.3  2003/01/13 19:52:10  f-najjar
! Remove WRITE statement
!
! Revision 1.2  2003/01/13 19:31:44  f-najjar
! Included I8 in FORMAT call
!
! Revision 1.1  2002/10/25 14:18:35  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







