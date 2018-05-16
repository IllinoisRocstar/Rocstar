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
! Purpose: convert vapor energy of particles to gas energy
!
! Description: none.
!
! Input: region  = current region.
!
! Output: region%levels(iLev)%plag%cv
!         region%levels(iLev)%mixt%cv
!
! Notes:
!
!   Need to investigate release of energy into dummy cells
!
!******************************************************************************
!
! $Id: INRT_VaporEnergyConversion.F90,v 1.4 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_VaporEnergyConversion( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModMixture,    ONLY : t_mixt
  USE ModPartLag,    ONLY : t_plag
  USE ModInteract,   ONLY : t_inrt_interact
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

#ifdef PLAG
  USE PLAG_ModParameters
#endif

#ifdef RFLO
  USE ModInterfaces, ONLY: RFLO_GetDimensDummy,RFLO_GetCellOffset

#include "Indexing.h"
#endif

  USE ModInterfaces, ONLY: MixtPerf_G_CpR, MixtPerf_R_M, MixtureProperties
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls, ic

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif
  INTEGER :: indPlagVapor, errorFlag, ibc, iec, indCp, indMol, nPcls, iCell

  INTEGER, POINTER, DIMENSION(:,:) :: pPlagAiv

  REAL(RFREAL) :: releaseTemp, gasTemp, keepFrac, rgas, gamma, hcap, deltaEner

  REAL(RFREAL), POINTER, DIMENSION(:)   :: vol
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pMixtCv, pMixtDv, pMixtGv
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pPlagCv, pPlagArv

  REAL(RFREAL), ALLOCATABLE :: vaporTot(:)

  TYPE(t_inrt_interact), POINTER :: pInrtBurn
  TYPE(t_mixt),          POINTER :: pMixt
  TYPE(t_plag),          POINTER :: pPlag
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_VaporEnergyConversion.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_VaporEnergyConversion',&
  'INRT_VaporEnergyConversion.F90' )

#ifdef PLAG
! begin -----------------------------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)
  pMixt => region%levels(iLev)%mixt
  pPlag => region%levels(iLev)%plag
  vol   => region%levels(iLev)%grid%vol
  indCp  = pMixt%indCp
  indMol = pMixt%indMol
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot
  pMixt => region%mixt
  pPlag => region%plag
  vol   => region%grid%vol
  indCp  = region%mixtInput%indCp
  indMol = region%mixtInput%indMol
#endif

! Check if there are any particles --------------------------------------------

  nPcls = 0
  IF (global%plagUsed) nPcls = pPlag%nPcls

  IF (nPcls < 1) GO TO 9

! Check if particle vapor energy is active ------------------------------------

  indPlagVapor = region%inrtInput%indPlagVapor
  IF (region%inrtInput%globActiveness(indPlagVapor) /= INRT_ACT_ACTIVE) GO TO 9

! Set pointers and values -----------------------------------------------------

  pMixtCv  => pMixt%cv
  pMixtDv  => pMixt%dv
  pMixtGv  => pMixt%gv

  pPlagCv  => pPlag%cv
  pPlagAiv => pPlag%aiv
  pPlagArv => pPlag%arv

  pInrtBurn => region%inrtInput%inrts(INRT_TYPE_BURNING)

! temperature below which to release vapor energy, which is set to be the same
! as the temperature above which burning shunts energy to vapor

  releaseTemp = pInrtBurn%data(INRT_DAT_BURNING_VAPOR_TEMP)

  IF (releaseTemp < 500._RFREAL .OR. releaseTemp > 10000._RFREAL) THEN
    CALL ErrorStop( global,ERR_INRT_BADVAL,__LINE__ )
  END IF ! releaseTemp

! need mixture temperature and spht, so compute all dv and gv for mixture

  IF ( region%mixtInput%gasModel == GAS_MODEL_TCPERF ) THEN ! cp, Mol=const.
    CALL MixtureProperties(region,ibc,iec,.FALSE.)
  ELSE
    CALL MixtureProperties(region,ibc,iec,.TRUE.)
  END IF ! region%mixtInput%gasModel

! allocate and zero temporary array -------------------------------------------

  ALLOCATE( vaporTot(ibc:iec),stat=errorFlag )
  errorFlag = global%error
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  vaporTot(:) = 0._RFREAL

! sum up vapor energy in cells ------------------------------------------------

  DO iPcls = 1,nPcls

    iCell = pPlagAiv(AIV_PLAG_ICELLS,iPcls)

    vaporTot(iCell) = vaporTot(iCell) + pPlagCv( CV_PLAG_ENERVAPOR,iPcls) * &
                                        pPlagArv(ARV_PLAG_SPLOAD,  iPcls)
  END DO ! iPcls

  DO ic = ibc,iec

    keepFrac = -1._RFREAL
    gasTemp  = pMixtDv(DV_MIXT_TEMP,ic)

    IF (vaporTot(ic) > 0._RFREAL .AND. gasTemp < releaseTemp) THEN

! --- normalize vaporTot by cell volume

! --- Note: normalization is necessary because cv values are not per volume
! --- for particles, but are per volume for continuua, in contrast to rhs
! --- values, which are not per unit volume either for particles or continuua.

      vaporTot(ic) = vaporTot(ic) / vol(ic)

! --- Compute heat capacity of gas (this is based on how revision 1.12 of
! --- perfgasDependentVars.F90 computes these quantities)

      rgas  = MixtPerf_R_M(  pMixtGv(GV_MIXT_MOL,ic*indMol))
      gamma = MixtPerf_G_CpR(pMixtGv(GV_MIXT_CP, ic*indCp ),rgas)

      hcap  = (rgas/(gamma - 1._RFREAL)) * pMixtCv(CV_MIXT_DENS,ic)

! --- Compute difference between current energy and that which would
! --- result from raising the gas temperature to releaseTemp

      deltaEner = hcap*(releaseTemp - gasTemp)

! --- Augment gas energy

      IF (vaporTot(ic) > deltaEner) THEN

        pMixtCv(CV_MIXT_ENER,ic) = pMixtCv(CV_MIXT_ENER,ic) + deltaEner

        keepFrac = 1._RFREAL - deltaEner / vaporTot(ic)

      ELSE

        pMixtCv(CV_MIXT_ENER,ic) = pMixtCv(CV_MIXT_ENER,ic) + vaporTot(ic)

        keepFrac = 0._RFREAL

      END IF ! vaporTot(ic) > deltaEner

    ELSE

      keepFrac = 1._RFREAL

    END IF ! vaporTot(ic) > 0._RFREAL .AND. gasTemp < releaseTemp

    IF (keepFrac < 0._RFREAL .OR. keepFrac > 1._RFREAL) THEN
      CALL ErrorStop( global,ERR_INVALID_VALUE,__LINE__ )
    END IF ! keepFrac

! - save value of keepFrac in the vaporTot array

    vaporTot(ic) = keepFrac

  END DO ! ic

! Loop over all the particles -------------------------------------------------

  DO iPcls = 1,nPcls

    iCell = pPlagAiv(AIV_PLAG_ICELLS,iPcls)

    pPlagCv(CV_PLAG_ENERVAPOR,iPcls) = pPlagCv(CV_PLAG_ENERVAPOR,iPcls) * &
                                       vaporTot(iCell)
  END DO ! iPcls

! Deallocate temporary array --------------------------------------------------

  DEALLOCATE( vaporTot,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! finalize --------------------------------------------------------------------

9 CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE INRT_VaporEnergyConversion

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_VaporEnergyConversion.F90,v $
! Revision 1.4  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/10/31 21:09:37  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/12/01 21:56:49  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2004/03/25 21:15:33  jferry
! added MixtureProperties to ModInterfaces list
!
! Revision 1.2  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
!******************************************************************************







