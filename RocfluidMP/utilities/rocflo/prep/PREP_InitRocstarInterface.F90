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
! Purpose: register grid and flow variables with GenX for pre-processing
!
! Description: none.
!
! Input: iReg       = current region number
!        region     = data of current region
!        wins, winv = name of fluid surface and volume windows
!
! Output: to Roccom.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_InitGenxInterface.F90,v 1.15 2009/08/27 14:04:53 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_InitGenxInterfacePrep( iReg,region,wins,winv )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE PREP_ModInterfaces, ONLY : RFLO_GetDimensDummyNodes
#ifdef PEUL
!  USE PREP_ModInterfaces, ONLY : PREP_PEUL_InitGenxInterface
#endif
#ifdef PLAG
!  USE PREP_ModInterfaces, ONLY : PREP_PLAG_InitGenxInterface
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  INTEGER :: iReg
  TYPE(t_region) :: region
  CHARACTER(CHRLEN) :: wins, winv

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  INTEGER :: iLev, bcType, pid, icount, errorFlag, ilb
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER, POINTER :: dims(:), iConstrType
  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_InitGenxInterfacePrep', &
                         'PREP_InitRocstarInterface.F90' )

! open data windows and register variables ------------------------------------

! input surface data

  CALL COM_new_dataitem( TRIM(wins)//'.bcflag'  ,'p',COM_INTEGER,1,'' )
  CALL COM_new_dataitem( TRIM(wins)//'.cnstr_type','p',COM_INTEGER,1,'' )

! input volume data

  CALL COM_new_dataitem( TRIM(winv)//'.rhof' ,'e',COM_DOUBLE,1,&
                          'kg/(m^3)')
  CALL COM_new_dataitem( TRIM(winv)//'.rhovf','e',COM_DOUBLE,3,&
                          'kg/(m^2 s)')
  CALL COM_new_dataitem( TRIM(winv)//'.rhoEf','e',COM_DOUBLE,1,&
                          '(J/kg)')

! store pointers to variables -------------------------------------------------

  ALLOCATE( dims(3),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  iLev   = region%currLevel
  icount = 0

! surface data

  DO iPatch=1,region%nPatches
    patch  => region%levels(iLev)%patches(iPatch)
    bcType =  patch%bcType
!    IF(bcType .NE. BC_SYMMETRY) THEN
    
    IF (patch%bcCoupled == BC_EXTERNAL) THEN        ! data from outside
      icount  = icount + 1
      pid     = iReg*REGOFF + icount

! --- burning pane?

      CALL COM_set_size(TRIM(wins)//'.bcflag',pid,1)
      CALL COM_set_array( TRIM(wins)//'.bcflag',pid,patch%bcFlag )

      IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
        patch%bcFlag(1) = 1   ! ignitable patch
      ELSE
        patch%bcFlag(1) = 0   ! non-ignitable patch
      ENDIF

! --- surface grid

      dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
      dims(2) = ABS(patch%l2end-patch%l2beg) + 2

      CALL COM_set_array_const( TRIM(wins)//'.:st2:',pid,dims )
      CALL COM_set_array( TRIM(wins)//'.nc',pid,patch%surfCoord )

! --- constraint type

      CALL COM_set_size(TRIM(wins)//'.cnstr_type',pid,1)
      CALL COM_allocate_array( TRIM(wins)//'.cnstr_type',pid,iConstrType )

      iConstrType = 0
      IF (bcType==BC_SYMMETRY_FREE) THEN
         iConstrType = 0
      ELSEIF (bcType==BC_SYMMETRY_FIXED) THEN
         iConstrType = 2
      ELSEIF (bcType==BC_SYMMETRY_XSLIDE) THEN
         iConstrType = 120
      ELSEIF (bcType==BC_SYMMETRY_YSLIDE) THEN
         iConstrType = 121
      ELSEIF (bcType==BC_SYMMETRY_ZSLIDE) THEN
         iConstrType = 122
      ELSEIF (bcType==BC_SYMMETRY_XYSLIDE) THEN
         iConstrType = -122
      ELSEIF (bcType==BC_SYMMETRY_XZSLIDE) THEN
         iConstrType = -121
      ELSEIF (bcType==BC_SYMMETRY_YZSLIDE) THEN
         iConstrType = -120
      ENDIF
      IF (bcType==BC_SLIPWALL_FREE) THEN
         iConstrType = 0
      ELSEIF (bcType==BC_SLIPWALL_FIXED) THEN
         iConstrType = 2
      ELSEIF (bcType==BC_SLIPWALL_XSLIDE) THEN
         iConstrType = 120
      ELSEIF (bcType==BC_SLIPWALL_YSLIDE) THEN
         iConstrType = 121
      ELSEIF (bcType==BC_SLIPWALL_ZSLIDE) THEN
         iConstrType = 122
      ELSEIF (bcType==BC_SLIPWALL_XYSLIDE) THEN
         iConstrType = -122
      ELSEIF (bcType==BC_SLIPWALL_XZSLIDE) THEN
         iConstrType = -121
      ELSEIF (bcType==BC_SLIPWALL_YZSLIDE) THEN
         iConstrType = -120
      ENDIF
      IF (bcType==BC_NOSLIPWALL_FREE) THEN
         iConstrType = 0
      ELSEIF (bcType==BC_NOSLIPWALL_FIXED) THEN
         iConstrType = 2
      ELSEIF (bcType==BC_NOSLIPWALL_XSLIDE) THEN
         iConstrType = 120
      ELSEIF (bcType==BC_NOSLIPWALL_YSLIDE) THEN
         iConstrType = 121
      ELSEIF (bcType==BC_NOSLIPWALL_ZSLIDE) THEN
         iConstrType = 122
      ELSEIF (bcType==BC_NOSLIPWALL_XYSLIDE) THEN
         iConstrType = -122
      ELSEIF (bcType==BC_NOSLIPWALL_XZSLIDE) THEN
         iConstrType = -121
      ELSEIF (bcType==BC_NOSLIPWALL_YZSLIDE) THEN
         iConstrType = -120
      ENDIF
      IF (bcType==BC_OUTFLOW_FREE) THEN
         iConstrType = 0
      ELSEIF (bcType==BC_OUTFLOW_FIXED) THEN
         iConstrType = 2
      ELSEIF (bcType==BC_OUTFLOW_XSLIDE) THEN
         iConstrType = 120
      ELSEIF (bcType==BC_OUTFLOW_YSLIDE) THEN
         iConstrType = 121
      ELSEIF (bcType==BC_OUTFLOW_ZSLIDE) THEN
         iConstrType = 122
      ELSEIF (bcType==BC_OUTFLOW_XYSLIDE) THEN
         iConstrType = -122
      ELSEIF (bcType==BC_OUTFLOW_XZSLIDE) THEN
         iConstrType = -121
      ELSEIF (bcType==BC_OUTFLOW_YZSLIDE) THEN
         iConstrType = -120
      ENDIF
      
   ELSE     ! internal BC
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
           (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) .OR. &
           (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
           (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) .OR. &
           (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE)) THEN
         
      ELSE
         icount  = icount + 1
         pid     = iReg*REGOFF + icount
         
         CALL COM_set_size(TRIM(wins)//'.bcflag',pid,1)
         CALL COM_set_array( TRIM(wins)//'.bcflag',pid,patch%bcFlag )
         
         patch%bcFlag(1) = 2   ! non-interacting patch
         
! ----- surface grid
         
         dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
         dims(2) = ABS(patch%l2end-patch%l2beg) + 2
         
         CALL COM_set_array_const( TRIM(wins)//'.:st2:',pid,dims )
         CALL COM_set_array( TRIM(wins)//'.nc',pid,patch%surfCoord )
         
! ----- constrain type
         
         CALL COM_set_size(TRIM(wins)//'.cnstr_type',pid,1)
         CALL COM_allocate_array( TRIM(wins)//'.cnstr_type',pid,iConstrType )
         
         iConstrType = 2
         IF (bcType==BC_SYMMETRY_FREE) THEN
            iConstrType = 0
         ELSEIF (bcType==BC_SYMMETRY_FIXED) THEN
            iConstrType = 2
         ELSEIF (bcType==BC_SYMMETRY_XSLIDE) THEN
            iConstrType = 120
         ELSEIF (bcType==BC_SYMMETRY_YSLIDE) THEN
            iConstrType = 121
         ELSEIF (bcType==BC_SYMMETRY_ZSLIDE) THEN
            iConstrType = 122
         ELSEIF (bcType==BC_SYMMETRY_XYSLIDE) THEN
            iConstrType = -122
         ELSEIF (bcType==BC_SYMMETRY_XZSLIDE) THEN
            iConstrType = -121
         ELSEIF (bcType==BC_SYMMETRY_YZSLIDE) THEN
            iConstrType = -120
         ENDIF
         IF (bcType==BC_SLIPWALL_FREE) THEN
            iConstrType = 0
         ELSEIF (bcType==BC_SLIPWALL_FIXED) THEN
            iConstrType = 2
         ELSEIF (bcType==BC_SLIPWALL_XSLIDE) THEN
            iConstrType = 120
         ELSEIF (bcType==BC_SLIPWALL_YSLIDE) THEN
            iConstrType = 121
         ELSEIF (bcType==BC_SLIPWALL_ZSLIDE) THEN
            iConstrType = 122
         ELSEIF (bcType==BC_SLIPWALL_XYSLIDE) THEN
            iConstrType = -122
         ELSEIF (bcType==BC_SLIPWALL_XZSLIDE) THEN
            iConstrType = -121
         ELSEIF (bcType==BC_SLIPWALL_YZSLIDE) THEN
            iConstrType = -120
         ENDIF
         IF (bcType==BC_NOSLIPWALL_FREE) THEN
            iConstrType = 0
         ELSEIF (bcType==BC_NOSLIPWALL_FIXED) THEN
            iConstrType = 2
         ELSEIF (bcType==BC_NOSLIPWALL_XSLIDE) THEN
            iConstrType = 120
         ELSEIF (bcType==BC_NOSLIPWALL_YSLIDE) THEN
            iConstrType = 121
         ELSEIF (bcType==BC_NOSLIPWALL_ZSLIDE) THEN
            iConstrType = 122
         ELSEIF (bcType==BC_NOSLIPWALL_XYSLIDE) THEN
            iConstrType = -122
         ELSEIF (bcType==BC_NOSLIPWALL_XZSLIDE) THEN
            iConstrType = -121
         ELSEIF (bcType==BC_NOSLIPWALL_YZSLIDE) THEN
            iConstrType = -120
         ENDIF
         IF (bcType==BC_OUTFLOW_FREE) THEN
            iConstrType = 0
         ELSEIF (bcType==BC_OUTFLOW_FIXED) THEN
            iConstrType = 2
         ELSEIF (bcType==BC_OUTFLOW_XSLIDE) THEN
            iConstrType = 120
         ELSEIF (bcType==BC_OUTFLOW_YSLIDE) THEN
            iConstrType = 121
         ELSEIF (bcType==BC_OUTFLOW_ZSLIDE) THEN
            iConstrType = 122
         ELSEIF (bcType==BC_OUTFLOW_XYSLIDE) THEN
            iConstrType = -122
         ELSEIF (bcType==BC_OUTFLOW_XZSLIDE) THEN
            iConstrType = -121
         ELSEIF (bcType==BC_OUTFLOW_YZSLIDE) THEN
            iConstrType = -120
         ENDIF
         
      ENDIF  ! block-interfaces
      
   ENDIF    ! external/internal BC
! ENDIF
  ENDDO      ! iPatch

! volume data

  pid = iReg*REGOFF

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  dims(1) = idnend - idnbeg + 1
  dims(2) = jdnend - jdnbeg + 1
  dims(3) = kdnend - kdnbeg + 1

  CALL COM_set_size( TRIM(winv)//".:st3:",pid,3, region%nDumCells)
  CALL COM_set_array_const( TRIM(winv)//".:st3:",pid,dims )
  CALL COM_set_array( TRIM(winv)//'.nc',pid, &
                      region%levels(iLev)%grid%xyz )

! set cv

  ilb = LBOUND(region%levels(iLev)%mixt%cv,2)

  CALL COM_set_array( TRIM(winv)//'.rhof',pid, &
       region%levels(iLev)%mixt%cv(1,ilb),5)
  CALL COM_set_array( TRIM(winv)//'.1-rhovf',pid, &
       region%levels(iLev)%mixt%cv(2,ilb),5)
  CALL COM_set_array( TRIM(winv)//'.2-rhovf',pid, &
       region%levels(iLev)%mixt%cv(3,ilb),5)
  CALL COM_set_array( TRIM(winv)//'.3-rhovf',pid, &
       region%levels(iLev)%mixt%cv(4,ilb),5)
  CALL COM_set_array( TRIM(winv)//'.rhoEf',pid, &
       region%levels(iLev)%mixt%cv(5,ilb),5)

! deallocate local arrays

  DEALLOCATE( dims )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! Genx initialization of physical modules -------------------------------------

#ifdef PEUL
!  IF (global%peulUsed) CALL PEUL_initGenxInterfacePrep( iReg,region,wins,winv )
#endif
#ifdef PLAG
!  IF (global%plagUsed) CALL PLAG_initGenxInterfacePrep( iReg,region,wins,winv )
#endif

! finalize --------------------------------------------------------------------

  CALL COM_window_init_done( TRIM(wins))
  CALL COM_window_init_done( TRIM(winv))

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_InitGenxInterfacePrep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_InitGenxInterface.F90,v $
! Revision 1.15  2009/08/27 14:04:53  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.14  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.13  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.12  2006/08/28 11:42:13  rfiedler
! Add grid motion constraint types for outflow BC.
!
! Revision 1.11  2006/08/24 15:05:29  rfiedler
! Support XYSLIDE, XZSLIDE, and YZSLIDE instead of TANGEN.  Use numbers not ICHAR.
!
! Revision 1.10  2005/06/19 05:33:07  wasistho
! shift index rocprop slipwalls and change default iConstrType
!
! Revision 1.9  2005/06/17 03:09:45  wasistho
! relocated cnstr_type kernel inside both external internal loops
!
! Revision 1.8  2005/06/17 02:23:47  jiao
! Fixed bug in registering cnstr_type.
!
! Revision 1.7  2005/06/16 22:34:05  wasistho
! constr_type to cnstr_type
!
! Revision 1.6  2005/06/13 04:37:25  wasistho
! registered constr_type for both internal external
!
! Revision 1.5  2005/06/13 02:17:15  wasistho
! registered constr_type
!
! Revision 1.4  2005/05/10 15:04:28  wasistho
! exclude block interfaces in registration of internal surfaces
!
! Revision 1.3  2005/04/18 18:12:08  wasistho
! registered non-interacting patches
!
! Revision 1.2  2004/12/03 03:29:02  wasistho
! rflo_modinterfacesprep to prep_modinterfaces
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.2  2004/06/30 04:07:06  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.1  2004/06/30 00:06:05  wasistho
! initial import for GEN3
!
!
!******************************************************************************







