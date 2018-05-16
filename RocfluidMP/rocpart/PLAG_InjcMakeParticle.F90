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
! Purpose: set the next particle diameter and superparticle loading
!          for the multiphase injection algorithm.
!
! Description: none.
!
! Input: region     = current region
!        injcDiamDist  = injection model type
!        diam       = particle diameter
!        spLoad     = superparticle loading
!
! Output: diam and spLoad 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InjcMakeParticle.F90,v 1.5 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InjcMakeParticle( region, injcDiamDist, diam, spLoad )

  USE ModDataTypes  
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag_input
  USE ModRandom, ONLY     : t_rand_data, Rand1Uniform, Rand1LogNormal,Rand1ImposedPDF
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER        :: injcDiamDist
  REAL(RFREAL)   :: diam, spLoad

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  
  REAL(RFREAL) :: diamMeanLog, injcDiamMean, injcStdDev
  REAL(RFREAL) :: alpha, diamPeakLog, injcDiamMax, injcDiamMin, &
                  injcDiamPeak, power,valMax
  INTEGER      :: locMaxPdf
  
  TYPE(t_plag_input), POINTER :: plagInput
  TYPE(t_global),     POINTER :: global
  REAL(RFREAL),       POINTER :: pdfvalues(:,:) 
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcMakeParticle.F90,v $ $Revision: 1.5 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InjcMakeParticle',&
  'PLAG_InjcMakeParticle.F90' )

! Select injection model ------------------------------------------------------

  SELECT CASE(injcDiamDist)

!- Compute particle diameter based on Log normal distribution -----------------
  
    CASE (PLAG_INJC_LOGNORM)
      injcDiamMean = region%plagInput%injcDiamMean
      injcStdDev   = region%plagInput%injcStdDev 
         
      diamMeanLog  = LOG(injcDiamMean)
      diam = Rand1LogNormal(diamMeanLog,injcStdDev,region%randData)
      spLoad = region%plagInput%spLoad

!- Compute particle diameter based on skewed Log distribution -----------------
  
    CASE (PLAG_INJC_LOGSKWD)
      injcDiamPeak = region%plagInput%injcDiamMean
      injcDiamMin  = region%plagInput%injcDiamMin
      injcDiamMax  = region%plagInput%injcDiamMax
      injcStdDev   = region%plagInput%injcStdDev 
      
      power        = 2.0_RFREAL
      diamPeakLog  = LOG(injcDiamPeak)
           
      alpha = diamPeakLog + injcStdDev*injcStdDev*power/ &
          ( (injcDiamMax/injcDiamPeak)**power - 1.0_RFREAL)

      diam   = PLAG_rand1LogSkewed(injcDiamPeak, injcStdDev, injcDiamMin, &
                               injcDiamMax, alpha, power, region%randData ) 
      spLoad = region%plagInput%spLoad

   CASE (PLAG_INJC_PDF)
      pdfvalues => region%plagInput%PDF%pdfvalues
      locMaxPdf =  region%plagInput%PDF%locMax
      valMax    =  region%plagInput%PDF%valMax

      diam   = rand1ImposedPDF(region%randData,pdfvalues,locMaxPdf,valMax)
      spLoad = region%plagInput%spLoad
               
    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )
       
  END SELECT
     
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )
    
CONTAINS

!******************************************************************************
  REAL(RFREAL) FUNCTION PLAG_rand1LogSkewed(dPeak,sDev,dMin,dMax,alpha, &
       power,rdata)
!******************************************************************************
! Skewed logarithmic distribution
! dPeak = peak
! sdev  = standard deviation
! dMin  = minimum
! dMax  = maximum
! alpha = median
! power = shape parameter

  REAL(RFREAL), INTENT(IN)  :: dMax,dPeak,dMin,sdev,alpha,power
  TYPE(t_rand_data), INTENT(INOUT) :: rdata

  INTEGER, PARAMETER :: iTerMax = 1000
  INTEGER   :: iTer

  REAL(RFREAL) :: dLogNorm, fSkewedRatio, xDistRand

  DO iTer = 1, iTerMax 
    dLogNorm = Rand1LogNormal(alpha,sdev,rdata)    

! decide to either accept or reject dLogNorm ----------------------------------

    IF ( dLogNorm >= dMin .AND. dLogNorm <= dMax ) THEN
      xDistRand = Rand1Uniform(rdata)

      fSkewedRatio = 1.0_RFREAL - ( (1.0_RFREAL-(dLogNorm/dMax)**power) &
                                  * (1.0_RFREAL-(dMin/dLogNorm)**power) )
      IF ( xDistRand > fSkewedRatio ) THEN
        PLAG_rand1LogSkewed = dLogNorm
        GOTO 8
      END IF ! xDistRand
      
    END IF ! dLogNorm

  END DO  ! iTer

  WRITE(STDOUT,'(A)') SOLVER_NAME//'### WARNING: PLAG_rand1LogSkewed failed!'
  WRITE(STDOUT,'(A)') SOLVER_NAME//'### Setting random value to peak'
  
  PLAG_rand1LogSkewed = dPeak
   
8 CONTINUE

  END FUNCTION PLAG_rand1LogSkewed
!******************************************************************************

END SUBROUTINE PLAG_injcMakeParticle

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcMakeParticle.F90,v $
! Revision 1.5  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/05/31 21:36:38  fnajjar
! Bug Fix to add spload for PLAG_INJC_PDF model
!
! Revision 1.2  2005/04/25 18:39:10  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.1  2004/12/01 20:57:41  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/06/16 23:06:33  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.3  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.2  2003/09/17 21:05:13  fnajjar
! Added infrastructure for skewed Log distribution in injection model
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************







