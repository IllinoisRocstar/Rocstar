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
! Purpose: initialize flow field in a region.
!
! Description: none.
!
! Input: iReg = region number
!        iLev = grid level
!
! Output: region = region data (mixt%cv)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_InitializeFlowField.F90,v 1.8 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE InitializeFlowField( iLev,region )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE PREP_ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModMixture, ONLY : t_mixt, t_mixt_input
  USE ModParameters
  USE PREP_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iLev

  TYPE(t_region), TARGET :: region

! ... loop variables
  INTEGER :: i,j,k, ijk, ijkN(8)

! ... local variables
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ibc, iec, ndum, ijkNBeg, ijkNEnd
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend, errorFlag

  REAL(RFREAL) :: gam1, qvel, xc, xb, xe, velX

  TYPE(t_mixt), POINTER :: mixt
  TYPE(t_mixt_input)    :: input

!******************************************************************************

  CALL RegisterFunction( region%global,'InitializeFlowField',&
  'PREP_InitializeFlowField.F90' )

! start

  mixt  => region%levels(iLev)%mixt
  input =  region%mixtInput
  ndum  =  region%nDumCells

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  ALLOCATE( mixt%cv(5,ibc:iec),stat=errorFlag )
  region%global%error = errorFlag
  IF (region%global%error /= 0) &
    CALL ErrorStop( region%global,ERR_ALLOCATE,__LINE__ )

  IF (input%prepIniCase==INITFLO_UNIFORM) THEN
    IF (input%iniXsplit > 0.5_RFREAL*HUGE(1._RFREAL)) THEN
      DO ijk = ibc,iec
        gam1 = region%global%refGamma - 1._RFREAL
        qvel = 0.5_RFREAL*(input%iniVelX*input%iniVelX+ &
                           input%iniVelY*input%iniVelY+ &
                           input%iniVelZ*input%iniVelZ)

        mixt%cv(CV_MIXT_DENS,ijk) = input%iniDens
        mixt%cv(CV_MIXT_XMOM,ijk) = input%iniVelX*input%iniDens
        mixt%cv(CV_MIXT_YMOM,ijk) = input%iniVelY*input%iniDens
        mixt%cv(CV_MIXT_ZMOM,ijk) = input%iniVelZ*input%iniDens
        mixt%cv(CV_MIXT_ENER,ijk) = input%iniPress/gam1 + input%iniDens*qvel
      ENDDO
    ELSE
      DO k=kdcbeg,kdcend
        DO j=jdcbeg,jdcend
          DO i=idcbeg,idcend

            ijk     = IndIJK(i,j,k,iCOff,ijCOff)
            ijkN(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
            ijkN(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
            ijkN(3) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
            ijkN(4) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
            ijkN(5) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
            ijkN(6) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
            ijkN(7) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
            ijkN(8) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)

            gam1 = region%global%refGamma - 1._RFREAL
            xc   = 0.125_RFREAL* &
                  (region%levels(iLev)%grid%xyz(XCOORD,ijkN(1)) + &
                   region%levels(iLev)%grid%xyz(XCOORD,ijkN(2)) + &
                   region%levels(iLev)%grid%xyz(XCOORD,ijkN(3)) + &
                   region%levels(iLev)%grid%xyz(XCOORD,ijkN(4)) + &
                   region%levels(iLev)%grid%xyz(XCOORD,ijkN(5)) + &
                   region%levels(iLev)%grid%xyz(XCOORD,ijkN(6)) + &
                   region%levels(iLev)%grid%xyz(XCOORD,ijkN(7)) + &
                   region%levels(iLev)%grid%xyz(XCOORD,ijkN(8)))
            IF (xc <= input%iniXsplit) THEN 
              qvel = 0.5_RFREAL*(input%iniVelX*input%iniVelX+ &
                                 input%iniVelY*input%iniVelY+ &
                                 input%iniVelZ*input%iniVelZ)

              mixt%cv(CV_MIXT_DENS,ijk) = input%iniDens
              mixt%cv(CV_MIXT_XMOM,ijk) = input%iniVelX*input%iniDens
              mixt%cv(CV_MIXT_YMOM,ijk) = input%iniVelY*input%iniDens
              mixt%cv(CV_MIXT_ZMOM,ijk) = input%iniVelZ*input%iniDens
              mixt%cv(CV_MIXT_ENER,ijk) = input%iniPress/gam1+input%iniDens*qvel
            ELSE
              qvel = 0.5_RFREAL*(input%iniVelX2*input%iniVelX2+ &
                                 input%iniVelY2*input%iniVelY2+ &
                                 input%iniVelZ2*input%iniVelZ2)

              mixt%cv(CV_MIXT_DENS,ijk) = input%iniDens2
              mixt%cv(CV_MIXT_XMOM,ijk) = input%iniVelX2*input%iniDens2
              mixt%cv(CV_MIXT_YMOM,ijk) = input%iniVelY2*input%iniDens2
              mixt%cv(CV_MIXT_ZMOM,ijk) = input%iniVelZ2*input%iniDens2
              mixt%cv(CV_MIXT_ENER,ijk) = input%iniPress2/gam1+input%iniDens2*qvel
            ENDIF  ! xc
          ENDDO    ! i
        ENDDO      ! j
      ENDDO        ! k
    ENDIF          ! iniXsplit
  ENDIF            ! prepIniCase

  IF (input%prepIniCase==INITFLO_PISTON_EXPAN) THEN
    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend

          ijk     = IndIJK(i,j,k,iCOff,ijCOff)
          ijkN(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          ijkN(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          ijkN(3) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
          ijkN(4) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
          ijkN(5) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
          ijkN(6) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
          ijkN(7) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
          ijkN(8) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)
          ijkNBeg = IndIJK(1            ,j  ,k  ,iNOff,ijNOff)
          ijkNEnd = IndIJK(idcend-ndum+1,j  ,k  ,iNOff,ijNOff)

          gam1 = region%global%refGamma - 1._RFREAL
          xc   = 0.125_RFREAL* &
                (region%levels(iLev)%grid%xyz(XCOORD,ijkN(1)) + &
                 region%levels(iLev)%grid%xyz(XCOORD,ijkN(2)) + &
                 region%levels(iLev)%grid%xyz(XCOORD,ijkN(3)) + &
                 region%levels(iLev)%grid%xyz(XCOORD,ijkN(4)) + &
                 region%levels(iLev)%grid%xyz(XCOORD,ijkN(5)) + &
                 region%levels(iLev)%grid%xyz(XCOORD,ijkN(6)) + &
                 region%levels(iLev)%grid%xyz(XCOORD,ijkN(7)) + &
                 region%levels(iLev)%grid%xyz(XCOORD,ijkN(8)))
          xb   = region%levels(iLev)%grid%xyz(XCOORD,ijkNBeg)
          xe   = region%levels(iLev)%grid%xyz(XCOORD,ijkNEnd)

!          velX = (xc-xb)/(xe-xb)*input%iniVelX           ! 1-block
          velX = (xc-input%iniXsplit)/ &                  ! multi-block
                 region%global%refLength*input%iniVelX 

          qvel = 0.5_RFREAL*(velX*velX+ &
                             input%iniVelY*input%iniVelY+ &
                             input%iniVelZ*input%iniVelZ)

          mixt%cv(CV_MIXT_DENS,ijk) = input%iniDens
          mixt%cv(CV_MIXT_XMOM,ijk) =          velX*input%iniDens
          mixt%cv(CV_MIXT_YMOM,ijk) = input%iniVelY*input%iniDens
          mixt%cv(CV_MIXT_ZMOM,ijk) = input%iniVelZ*input%iniDens
          mixt%cv(CV_MIXT_ENER,ijk) = input%iniPress/gam1+input%iniDens*qvel
        ENDDO  ! i
      ENDDO  ! j
    ENDDO  ! k
  ENDIF  ! iniCase

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE InitializeFlowField

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_InitializeFlowField.F90,v $
! Revision 1.8  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2005/09/21 20:36:52  wasistho
! added xsplit parameter in pistonExpand case
!
! Revision 1.5  2005/09/21 19:04:15  wasistho
! modified initial condition for expanding piston case
!
! Revision 1.4  2005/09/09 03:30:45  wasistho
! added Expanding-Piston case kernel
!
! Revision 1.3  2005/08/03 17:51:31  wasistho
! added domain-splitting feature to initial condition
!
! Revision 1.2  2004/12/03 03:29:08  wasistho
! rflo_modinterfacesprep to prep_modinterfaces
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.9  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.8  2003/04/23 22:51:43  olawlor
! Added TARGET dataitem to region variable, since
! pointers are cached into it.
!
! Revision 1.7  2003/03/20 22:27:56  haselbac
! Renamed ModInterfaces
!
! Revision 1.6  2003/03/20 19:44:22  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.5  2003/03/20 19:35:43  haselbac
! Modified RegFun call to avoid probs with long 'PREP_InitializeFlowField.F90' names
!
! Revision 1.4  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.3  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/02/21 23:25:07  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/02 15:57:08  jblazek
! Added flow initialization.
!
!******************************************************************************








