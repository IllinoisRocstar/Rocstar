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
! Purpose: compute gas dynamic properties of the mixture.
!
! Description: none.
!
! Input: region    = mixture variables and index pointers
!        inBeg     = first value to update
!        inEnd     = last value to update
!        gasUpdate = should the gas variables be also updated?
!
! Output: mixt%dv = dependent variables (p, T, c)
!         mixt%tv = transport variables (viscosity, heat conductivity)
!         mixt%gv = gas variables (cp, molecular mass)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: MixtureProperties.F90,v 1.4 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE MixtureProperties( region,inBeg,inEnd,gasUpdate )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : PerfgasGasVars, PerfgasDependentVars, &
                            PerfgasTransportVars
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: inBeg, inEnd

  LOGICAL :: gasUpdate

  TYPE(t_region) :: region

! ... local variables
  INTEGER :: iLev, indCp, indMol, gasModel, viscModel

  LOGICAL :: computeTv

  REAL(RFREAL)          :: prLam, refTemp, refVisc, suthCoef
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), tv(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'MixtureProperties',&
  'MixtureProperties.F90' )

! set local quantities to avoid duplicate calls -------------------------------

#ifdef RFLO
  iLev      =  region%currLevel
  indCp     =  region%levels(iLev)%mixt%indCp
  indMol    =  region%levels(iLev)%mixt%indMol
  prLam     =  region%levels(iLev)%mixt%prLam
  cv        => region%levels(iLev)%mixt%cv
  dv        => region%levels(iLev)%mixt%dv
  gv        => region%levels(iLev)%mixt%gv
  tv        => region%levels(iLev)%mixt%tv
#endif
#ifdef RFLU
  indCp     =  region%mixtInput%indCp
  indMol    =  region%mixtInput%indMol
  prLam     =  region%mixtInput%prLam
  cv        => region%mixt%cv
  dv        => region%mixt%dv
  gv        => region%mixt%gv
  tv        => region%mixt%tv
#endif

  gasModel = region%mixtInput%gasModel
  computeTv = region%mixtInput%computeTv

  viscModel = region%mixtInput%viscModel
  refVisc   = region%mixtInput%refVisc
  refTemp   = region%mixtInput%refTemp
  suthCoef  = region%mixtInput%suthCoef

! perfect gas -----------------------------------------------------------------

  IF (gasModel == GAS_MODEL_TCPERF) THEN

! - update gas variables

    IF (gasUpdate) CALL PerfgasGasVars( inBeg,inEnd,indCp,indMol, &
                                        region%global%refCp, &
                                        region%global%refGamma,cv,gv )

! - update dependent variables

    CALL PerfgasDependentVars( inBeg,inEnd,indCp,indMol,cv,gv,dv )

! - update transport variables

    IF (computeTv) THEN
      CALL PerfgasTransportVars( inBeg,inEnd,indCp,indMol,viscModel,prLam, &
                                 refVisc,refTemp,suthCoef,cv,dv,gv,tv )
    ENDIF

! -----------------------------------------------------------------------------

  ELSE
    CALL ErrorStop( region%global,ERR_OPTION_TYPE,__LINE__,'Species model?' )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE MixtureProperties

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtureProperties.F90,v $
! Revision 1.4  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/10/31 21:09:34  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/12/01 16:49:41  haselbac
! Initial revision after changing case
!
! Revision 1.13  2004/03/03 23:55:39  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.12  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.8  2003/04/10 23:26:48  fnajjar
! Modified calling sequence for perfgasTransportVars with viscosity model
!
! Revision 1.7  2002/09/09 13:59:13  haselbac
! mixtInput now under regions
!
! Revision 1.6  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/08/16 21:33:47  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.4  2002/07/05 23:20:46  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.3  2002/05/04 16:33:34  haselbac
! Added RFLU statements
!
! Revision 1.2  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/10 00:02:06  jblazek
! Added calculation of mixture properties.
!
!******************************************************************************







