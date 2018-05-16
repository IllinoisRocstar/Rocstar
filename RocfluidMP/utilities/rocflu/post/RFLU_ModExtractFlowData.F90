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
! Purpose: Collect routines to extract data from flow solution.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes:
!   1. These routines are hardcoded to extract data for particular flows on 
!      particular grids. This means that one CANNOT use these routines for any 
!      grid. 
!
! ******************************************************************************
!
! $Id: RFLU_ModExtractFlowData.F90,v 1.24 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2004-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModExtractFlowData

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region

#ifdef PLAG
  USE PLAG_ModParameters
#endif

  USE RFLU_ModExtractFlowDataUtils
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: &
    RCSIdentString = '$RCSfile: RFLU_ModExtractFlowData.F90,v $ $Revision: 1.24 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_ExtractFlowData

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_ExtractFlowDataBlasius, & 
             RFLU_ExtractFlowDataNSCBC, &
             RFLU_ExtractFlowDataLineFarf, &
             RFLU_ExtractFlowDataProudman, &
             RFLU_ExtractFlowDataSod, &
             RFLU_ExtractFlowDataSTG2D, &   
             RFLU_ExtractFlowDataSkews, &
             RFLU_WriteMeshBump

! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS








! ******************************************************************************
!
! Purpose: Extract data from flow solution.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowData(pRegion)
                 
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowData', &
                        'RFLU_ModExtractFlowData.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extracting data from flow solution...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel
  
! ******************************************************************************
! Initialize flow field based on user input
! ******************************************************************************

  SELECT CASE ( global%casename )
  
! ==============================================================================
!   Incompressible laminar flat plate
! ==============================================================================  
  
    CASE ( "lfpbli-16.48x8.8x3"   ,"lfpbli-32.96x16.16x3", & 
           "lfpbli-64.192x32.32x3","lfpbli-128.384x32.32x3"  )
      CALL RFLU_ExtractFlowDataBlasius(pRegion)
    CASE( "lfpblim-64x16x1" , "lfpblim-128x32x1", & 
          "lfpblim-256x64x1", "lfpblim-512x128x1" ) 
      CALL RFLU_ExtractFlowDataBlasius(pRegion)

! ==============================================================================
!   ONERA C0
! ==============================================================================

    CASE ( "onera_c0_2d_100x50" )
      CALL RFLU_ExtractFlowDataProudman(pRegion)

! ==============================================================================
!   Skews diffracting shock
! ==============================================================================

    CASE ( "skews_ms2p0","skews_ms3p0","skews_ms4p0" )
      CALL RFLU_ExtractFlowDataSkews(pRegion)
 
! ==============================================================================
!   Shock tubes
! ==============================================================================

    CASE ( "st_sod1","st_sod1_mp2","st_sod2","st_sod2_mp2" )
      CALL RFLU_ExtractFlowDataSod(pRegion)
    CASE ( "stg1d" )
      CALL RFLU_ExtractFlowDataSTG1D(pRegion)      
    CASE ( "stg2d" )
      CALL RFLU_ExtractFlowDataSTG2D(pRegion)      

! ==============================================================================
!   Sommerfeld shock-particle interaction
! ==============================================================================

    CASE ( "somm_spi" )
      CALL RFLU_ExtractFlowDataSommSPI(pRegion)
 
! ==============================================================================
!   NSCBC test cases 'nscbcX'
! ==============================================================================

    CASE ( "nscbc1","nscbc2","nscbc3","nscbc4","nscbc5","nscbc6","nscbc7", & 
           "nscbc8" )
      CALL RFLU_ExtractFlowDataNSCBC(pRegion)

! ==============================================================================
!   NSCBC farfield 
! ==============================================================================

    CASE ( "farf" )
      CALL RFLU_ExtractFlowDataLineFarf(pRegion)
 
! ==============================================================================
!   Bump test case 'bumpq10'
! ==============================================================================

    CASE ( "bumpq10" )
      CALL RFLU_WriteMeshBump(pRegion)
 
! ==============================================================================
!   Default - due to input error or missing CALL in this routine
! ==============================================================================  
            
    CASE DEFAULT 
      global%warnCounter = global%warnCounter + 1
    
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,2(1X,A))') SOLVER_NAME,'*** WARNING ***', & 
              'Extraction of data not available.', & 
              'Returning to calling procedure.'                                   
      END IF ! global%verbLevel
  END SELECT ! global%casename

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extracting data from flow solution done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowData









! ******************************************************************************
!
! Purpose: Extract data for Blasius flow
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Many assumptions are made, most of which are documented directly in the 
!      code below. Others are listed in the following: 
!   2. Plate is assumed to lie on x-z plane, with the plate leading edge at x=0,
!      and plate being defined by plane y=0.
!   3. Plate length is assumed to be unity, so that x-stations at which 
!      velocity profiles are extracted are equally spaced in [0,1].
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataBlasius(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1,iFileName2
  INTEGER, PARAMETER :: NPOINTS = 18
  INTEGER :: errorFlag,i,icg,ifl,ile,ite,ivp,j,jbl,jfs,k,nvp
  INTEGER, DIMENSION(:), ALLOCATABLE :: iloc
  REAL(RFREAL) :: cf,cfTheory,delta0,delta0Theory,delta1,delta1Theory,delta2, &
                  delta2Theory,delta3,delta3Theory,dist1,dist2,dy,eta,inter, &
                  ir,lRef,r,ReRef,Rex,ru,rv,slope,term1,term2,u,uNorm,uRef, &
                  v,vNorm,x,y
  REAL(RFREAL), DIMENSION(NPOINTS), PARAMETER :: etaBlas = & 
    (/0.0000_RFREAL,0.1000_RFREAL,0.2000_RFREAL,0.3000_RFREAL,0.4000_RFREAL, &
      0.5000_RFREAL,0.6000_RFREAL,0.8000_RFREAL,1.0000_RFREAL,1.2000_RFREAL, &
      1.4000_RFREAL,1.6000_RFREAL,1.8000_RFREAL,2.0000_RFREAL,2.5000_RFREAl, &
      3.0000_RFREAL,3.5000_RFREAL,4.0000_RFREAL/)                  
  REAL(RFREAL), DIMENSION(NPOINTS), PARAMETER :: uNormBlas = & 
    (/0.0000_RFREAL,0.0664_RFREAL,0.1328_RFREAL,0.1989_RFREAL,0.2647_RFREAL, &
      0.3298_RFREAL,0.3938_RFREAL,0.5168_RFREAL,0.6298_RFREAL,0.7290_RFREAL, &
      0.8115_RFREAL,0.8761_RFREAL,0.9233_RFREAL,0.9555_RFREAL,0.9916_RFREAl, &
      0.9990_RFREAL,0.9999_RFREAL,1.0000_RFREAL/)
  REAL(RFREAL), DIMENSION(NPOINTS), PARAMETER :: vNormBlas = & 
    (/0.0000_RFREAL,0.0033_RFREAL,0.0133_RFREAL,0.0298_RFREAL,0.0528_RFREAL, &
      0.0821_RFREAL,0.1173_RFREAL,0.2033_RFREAL,0.3048_RFREAL,0.4136_RFREAL, &
      0.5206_RFREAL,0.6172_RFREAL,0.6972_RFREAL,0.7581_RFREAL,0.8373_RFREAl, &
      0.8482_RFREAL,0.8600_RFREAL,0.8604_RFREAL/)      
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: yu
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataBlasius', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for Blasius flow...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

  nvp = 9 ! Number of equally-spaced stations at which profiles are extracted

! ******************************************************************************
! Get data about freestream conditions
! ******************************************************************************

  uRef  = global%refVelocity
  ReRef = global%refREnum
  lRef  = global%refLength 

! ******************************************************************************
! Get data about dimensions
!   ile         i-index of first cell on plate (leading edge)
!   ite         i-index of last cell on plate (trailing edge)
!   jbl         j-index of (nominally) last cell in boundary layer
!   jfs         j-index of cell abutting freestream boundary
!   k           k-index of middle layer of cells
! ******************************************************************************

  SELECT CASE ( global%casename ) 
    CASE ( "lfpbli-16.48x8.8x3" ) 
      ile =  17 
      ite =  64
      jbl =   8
      jfs =  16
      k   =   2
    CASE ( "lfpbli-32.96x16.16x3" ) 
      ile =  33 
      ite = 128
      jbl =  16
      jfs =  32
      k   =   2  
    CASE ( "lfpbli-64.192x32.32x3" ) 
      ile =  65 
      ite = 256
      jbl =  32
      jfs =  64
      k   =   2 
    CASE ( "lfpbli-128.384x64.64x3" ) 
      ile = 129 
      ite = 512
      jbl =  64
      jfs = 128
      k   =   2   
    CASE ( "lfpblim-64x16x1" ) 
      ile =  17
      ite =  64
      jbl =   8
      jfs =  16
      k   =   1 
    CASE ( "lfpblim-128x32x1" ) 
      ile =  33
      ite = 128
      jbl =  16
      jfs =  32
      k   =   1      
    CASE ( "lfpblim-256x64x1" ) 
      ile =  65
      ite = 256
      jbl =  32
      jfs =  64
      k   =   1
    CASE ( "lfpblim-512x128x1" ) 
      ile = 129
      ite = 512
      jbl =  64
      jfs = 128
      k   =   1                               
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%casename

! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************

  ALLOCATE(yu(2,0:jfs),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'yu')
  END IF ! global%error

  ALLOCATE(iloc(nvp),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'iloc')
  END IF ! global%error  

! ******************************************************************************
! Find i-indices of x-stations at which data about velocity profiles is to be
! extracted. NOTE assume that the cell preceding a given cell in the stream-
! wise direction has i-index smaller by one. 
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
      'Determining stations at which velocity data is extracted...'
  END IF ! global%verbLevel 
   
  i = ile 

  DO ivp = 1,nvp
    x = ivp*0.1_RFREAL
  
    emptyLoop: DO 
      icg = i + (k-1)*ite*jfs ! NOTE no j-term because always zero here 

      dist1 = pGrid%cofg(XCOORD,icg  ) - x
      dist2 = pGrid%cofg(XCOORD,icg-1) - x 
            
      IF ( (dist1 > 0.0_RFREAL) .AND. (dist2 < 0.0_RFREAL) ) THEN 
        IF ( ABS(dist1) < ABS(dist2) ) THEN 
          iloc(ivp) = i
          icg = i + (k-1)*ite*jfs
        ELSE 
          iloc(ivp) = i-1  
          icg = i - 1 + (k-1)*ite*jfs            
        END IF ! ABS

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,5X,A,1X,I2,1X,A,1X,I3,A,1X,E13.6)') SOLVER_NAME, &
            'Station',ivp,'located at i=',iloc(ivp),', x=',pGrid%cofg(XCOORD,icg)
        END IF ! global%verbLevel
        
        EXIT emptyLoop  
      ELSE 
        i = i + 1
      END IF ! dist1      
    END DO emptyLoop
  END DO ! ivp

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
      'Determining stations at which velocity data is extracted done.'
  END IF ! global%verbLevel

! ******************************************************************************
! Extract velocity profiles
! ******************************************************************************

  DO ivp = 1,nvp
    i = iloc(ivp)

    yu(1,0) = 0.0_RFREAL
    yu(2,0) = 0.0_RFREAL

! ==============================================================================
!   Open file
! ==============================================================================
      
    WRITE(iFileName2,'(A,I2.2,A)') 'blasius-vel',ivp,'.dat'

    OPEN(IF_EXTR_DATA2,FILE=iFileName2,FORM='FORMATTED',STATUS='UNKNOWN', &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName2))
    END IF ! global%error 

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Writing velocity-profile data to file: '// & 
                               TRIM(iFileName2)
    END IF ! global%verbLevel              

! ==============================================================================
!   Normalize velocity and write to file at selected stations
! ==============================================================================
            
    DO j = 1,jfs
      icg = i + (j-1)*ite + (k-1)*ite*jfs
      
      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)
      
      Rex = ReRef*x/lRef
      eta = y/(2.0_RFREAL*x)*SQRT(Rex)
      
      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      rv = pRegion%mixt%cv(CV_MIXT_YMOM,icg)      
      
      ir = 1.0_RFREAL/r
      u  = ir*ru
      v  = ir*rv
      
      uNorm = u/uRef
      vNorm = v/uRef
      
      yu(1,j) = y
      yu(2,j) = uNorm
      
      WRITE(IF_EXTR_DATA2,'(3(1X,E13.6))') eta,uNorm,vNorm*SQRT(Rex)
    END DO ! j

! ==============================================================================
!   Close file
! ==============================================================================
      
    CLOSE(IF_EXTR_DATA2,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '// & 
                     TRIM(iFileName2))
    END IF ! global%error       
  END DO ! ivp

! ******************************************************************************
! Extract thicknesses
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  iFileName1 = 'blasius-thick.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                             'Writing thickness data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

! ==============================================================================
! Extract data
! ==============================================================================

  DO i = ile,ite
    icg = i + (k-1)*ite*jfs ! NOTE no j-term because always zero here 
      
    x = pGrid%cofg(XCOORD,icg)
      
    Rex = ReRef*x/lRef
 
! ------------------------------------------------------------------------------
!   Compute normalized velocity profile
! ------------------------------------------------------------------------------  
  
    yu(1,0) = 0.0_RFREAL
    yu(2,0) = 0.0_RFREAL  
  
    DO j = 1,jfs
      icg = i + (j-1)*ite + (k-1)*ite*jfs
      
      y = pGrid%cofg(YCOORD,icg)
      
      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      
      ir = 1.0_RFREAL/r
      u  = ir*ru
      
      uNorm = u/uRef
      
      yu(1,j) = y
      yu(2,j) = uNorm
    END DO ! j  
  
! ------------------------------------------------------------------------------  
!   Determine boundary layer thickness (indicated by normalized streamwise 
!   velocity of 0.99).
! ------------------------------------------------------------------------------  

    delta0       = CRAZY_VALUE_INT
    delta0Theory = 5.0_RFREAL*x/SQRT(Rex)
        
    DO j = 1,jfs
      IF ( yu(2,j-1) < 0.99_RFREAL .AND. yu(2,j) >= 0.99_RFREAL ) THEN 
        slope = (yu(2,j  )         - yu(2,j-1)          )/(yu(1,j) - yu(1,j-1)) 
        inter = (yu(2,j-1)*yu(1,j) - yu(2,j  )*yu(1,j-1))/(yu(1,j) - yu(1,j-1)) 
        delta0 = (0.99_RFREAL - inter)/slope
      END IF ! yu            
    END DO ! j    
        
! ------------------------------------------------------------------------------  
!   Integrate normalized velocity to get displacement, momentum, and energy 
!   thicknesses and write thicknesses to file. For the moment, use simple 
!   trapezoidal rule. 
! ------------------------------------------------------------------------------  

    delta1 = 0.0_RFREAL
    delta2 = 0.0_RFREAL
    delta3 = 0.0_RFREAL
    
    delta1Theory = 1.720_RFREAL*x/SQRT(Rex)
    delta2Theory = 0.664_RFREAL*x/SQRT(Rex)
    delta3Theory = 1.044_RFREAL*x/SQRT(Rex)            
    
    DO j = 1,jbl
      dy = yu(1,j) - yu(1,j-1)
                       
      term1 = yu(2,j-1)                 
      term2 = yu(2,j  )                       
                                    
      delta1 = delta1 & 
             + 0.5_RFREAL*dy*(        (1.0_RFREAL-      term1) & 
                              +       (1.0_RFREAL-      term2))
      delta2 = delta2 & 
             + 0.5_RFREAL*dy*(  term1*(1.0_RFREAL-      term1) & 
                              + term2*(1.0_RFREAL-      term2))
      delta3 = delta3 & 
             + 0.5_RFREAL*dy*(  term1*(1.0_RFREAL-term1*term1) & 
                              + term2*(1.0_RFREAL-term2*term2))                                             
    END DO ! j

    WRITE(IF_EXTR_DATA1,'(9(1X,E13.6))') Rex,delta0/x, &
                                         delta1/delta0, &
                                         delta2/delta0, &
                                         delta3/delta0, & 
                                         delta0Theory/x, &
                                         delta1Theory/delta0Theory, & 
                                         delta2Theory/delta0Theory, &
                                         delta3Theory/delta0Theory                                         
  END DO ! i

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(yu,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'yu')
  END IF ! global%error

  DEALLOCATE(iloc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'iloc')
  END IF ! global%error 

! ******************************************************************************
! Write file with skin-friction data. NOTE assume that flat plate is always on
! patch 2. 
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  iFileName1 = 'blasius-cf.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                             'Writing skin-friction data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

! ==============================================================================
! Write data
! ==============================================================================

  pPatch => pRegion%patches(2) ! NOTE assumption

  DO ifl = 1+(k-1)*(ite-ile+1),k*(ite-ile+1)
    x = pPatch%fc(XCOORD,ifl)
        
    Rex = ReRef*x/lRef
  
    cfTheory = 0.664_RFREAL/SQRT(Rex)  
    cf       = pPatch%cf(XCOORD,ifl)
  
    WRITE(IF_EXTR_DATA1,'(3(1X,E13.6))') Rex,cfTheory,cf
  END DO ! ifl

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error    

! ******************************************************************************
! Write file with exact velocity profile data for comparison purposes
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  iFileName1 = 'blasius-vel-exact.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ==============================================================================
! Write data
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                             'Writing exact velocity-profile data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

  DO j = 1,NPOINTS
    WRITE(IF_EXTR_DATA1,'(3(1X,E11.4))') etaBlas(j),uNormBlas(j),vNormBlas(j)
  END DO ! j

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract data for Blasius flow done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataBlasius









! ******************************************************************************
!
! Purpose: Extract Line data for NSCBC case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataNSCBC(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: nExtract
  INTEGER :: errorFlag,icg,ix
  REAL(RFREAL) :: a,p,r,u,v,w,M,xx,yy,angle
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataNSCBC', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Line data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pPatch => pRegion%patches(1)

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line_',global%currentTime,'.plt'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to  '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Compute number of cells
! ******************************************************************************

  nExtract = pPatch%nBFaces

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************
  DO ix = 1,nExtract
    icg  = pPatch%bf2c(ix)

    xx = pGrid%cofg(XCOORD,icg)
    yy = pGrid%cofg(YCOORD,icg)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg)

    WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') xx,yy,r,u,v,p
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract Line data for NSCBC case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataNSCBC








! ******************************************************************************
!
! Purpose: Extract Surf data for Farfield 
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataLineFarf(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1,iFileName2,iFileName3,iFileName4
  INTEGER :: nExtract
  INTEGER :: errorFlag,icg,ix,icg1,icg2,icg3,icg4
  REAL(RFREAL) :: a,p,r,u,v,w,M,xx,yy,angle,radius,distance
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataLineFarf', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Surf data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pPatch => pRegion%patches(1)

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line1_',global%currentTime,'.plt'

  WRITE(iFileName2,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line2_',global%currentTime,'.plt'

  WRITE(iFileName3,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line3_',global%currentTime,'.plt'

  WRITE(iFileName4,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line4_',global%currentTime,'.plt'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

  OPEN(IF_EXTR_DATA2,FILE=iFileName2,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName2))
  END IF ! global%error

  OPEN(50,FILE=iFileName3,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName3))
  END IF ! global%error

  OPEN(51,FILE=iFileName4,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName4))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Compute number of cells
! ******************************************************************************

  nExtract = 95 
  radius   = 0.006125_RFREAL

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************

  DO ix = 1,nExtract
    icg   = pPatch%bf2c(ix)
    icg1  =   1 + (ix-1)*258 
    icg2  =  65 + (ix-1)*258 
    icg3  = 130 + (ix-1)*258 
    icg4  = 195 + (ix-1)*258 

    xx = pGrid%cofg(XCOORD,icg1)
    yy = pGrid%cofg(YCOORD,icg1)
    distance = SQRT(xx*xx + yy*yy)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg1)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg1)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg1)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg1)
    WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') distance-radius,r,u,v,p 

    xx = pGrid%cofg(XCOORD,icg2)
    yy = pGrid%cofg(YCOORD,icg2)
    distance = SQRT(xx*xx + yy*yy)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg2)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg2)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg2)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg2)
    WRITE(IF_EXTR_DATA2,'(6(1X,E23.16))') distance-radius,r,u,v,p 

    xx = pGrid%cofg(XCOORD,icg3)
    yy = pGrid%cofg(YCOORD,icg3)
    distance = SQRT(xx*xx + yy*yy)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg3)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg3)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg3)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg3)
    WRITE(50,'(6(1X,E23.16))') distance-radius,r,u,v,p 

    xx = pGrid%cofg(XCOORD,icg4)
    yy = pGrid%cofg(YCOORD,icg4)
    distance = SQRT(xx*xx + yy*yy)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg4)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg4)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg4)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg4)
    WRITE(51,'(6(1X,E23.16))') distance-radius,r,u,v,p 
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

  CLOSE(IF_EXTR_DATA2,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName2))
  END IF ! global%error

  CLOSE(50,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName3))
  END IF ! global%error

  CLOSE(51,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName4))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract Surface data for Farf case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataLineFarf










! ******************************************************************************
!
! Purpose: Extract data for Proudman flow
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Many assumptions are made, most of which are documented directly in the 
!      code below, otherwise also assume that case satisfies same restrictions
!      as those for defining exact solution.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataProudman(pRegion)

  USE RFLU_ModExactFlow, ONLY: RFLU_ComputeExactFlowProudman
  USE RFLU_ModFlowHardCode, ONLY: RFLU_GetParamsHardCodeProudman

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,i,icg,ivp,j,nvp,nx,ny
  INTEGER, DIMENSION(:), ALLOCATABLE :: iloc
  REAL(RFREAL) :: dInc,height,ir,mInj,vInj,p,pTot,r,ru,rv,u,uNorm,v,vNorm,w,x,y
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataProudman', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for Blasius flow...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

  nvp = 9

! ******************************************************************************
! Get data about flow and geometry
! ******************************************************************************

  CALL RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)

  height = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))

! ******************************************************************************
! Get data about dimensions
! ******************************************************************************

  SELECT CASE ( global%casename ) 
    CASE ( "onera_c0_2d_100x50" ) 
      nx = 100 
      ny =  50            
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%casename

! ******************************************************************************
! Extract data at selected stations
! ******************************************************************************

  DO ivp = 1,nvp
    i = nx/(nvp+1)*ivp

! ==============================================================================
!   Open file
! ==============================================================================
      
    WRITE(iFileName1,'(A,I2.2,A)') 'onera_c0-vel',ivp,'.dat'

    OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
    END IF ! global%error 

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Writing velocity-profile data to file: '// & 
                               TRIM(iFileName1)
    END IF ! global%verbLevel              

! ==============================================================================
!   Normalize velocity and write to file at selected stations. NOTE assume grid 
!   is uniform in x-direction so that stations are equally spaced.
! ==============================================================================
            
    DO j = 1,ny
      icg = i + (j-1)*nx ! NOTE assumption about cell numbering

      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)

      IF ( global%verbLevel > VERBOSE_NONE .AND. j == 1 ) THEN
        WRITE(STDOUT,'(A,5X,A,1X,I2,1X,A,1X,E13.6)') SOLVER_NAME,'Station', &
                                                     ivp,'located at x=',x
      END IF ! global%verbLevel      
      
      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      rv = pRegion%mixt%cv(CV_MIXT_YMOM,icg)      
      
      ir = 1.0_RFREAL/r
      u  = ir*ru
      v  = ir*rv
      
      uNorm = -u/(0.5_RFREAL*global%pi*x/height*vInj)
      vNorm =  v/vInj
      
      WRITE(IF_EXTR_DATA1,'(3(1X,E13.6))') (height-y)/height,uNorm,vNorm
    END DO ! j

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '// & 
                     TRIM(iFileName1))
    END IF ! global%error
  END DO ! ivp       

! ******************************************************************************
! Write file with exact velocity profile data for comparison purposes
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  iFileName1 = 'onera_c0-vel-exact.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ==============================================================================
! Write data
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                             'Writing exact velocity-profile data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

  DO j = 1,ny
    icg = 1 + (j-1)*nx ! NOTE set i to 1 as only need y-coordinate variation

    x = pGrid%cofg(XCOORD,icg)
    y = pGrid%cofg(YCOORD,icg)

    CALL RFLU_ComputeExactFlowProudman(global,x,y,height,dInc,vInj,pTot, &
                                       r,u,v,w,p)

    uNorm = -u/(0.5_RFREAL*global%pi*x/height*vInj)
    vNorm =  v/vInj

    WRITE(IF_EXTR_DATA1,'(3(1X,E13.6))') (height-y)/height,uNorm,vNorm
  END DO ! j

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract data for Proudman flow done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataProudman






! ******************************************************************************
!
! Purpose: Extract data for Sod shocktube case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie along x-axis.
!   2. Assume cross-section to have 3x3 cells, so can find number of cells
!      along x-axis by dividing total number of cells by 9.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSod(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgBeg,icgEnd,nCellsX
  REAL(RFREAL) :: a,cp,mw,p,r,T,u,v
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSod', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for Sod shocktube case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file for data 
! ******************************************************************************

  iFileName1 = 'sodst.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

! ******************************************************************************
! Find cell indices for extraction
! ******************************************************************************

  nCellsX = pGrid%nCellsTot/9 ! NOTE integer division

  icgBeg = 4*nCellsX + 1
  icgEnd = 5*nCellsX

! ******************************************************************************
! Extract along line of cells 
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%gasModel )
    CASE ( GAS_MODEL_TCPERF ) 
      DO icg = icgBeg,icgEnd
        r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        u = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
        p = pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T = pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a = pRegion%mixt%dv(DV_MIXT_SOUN,icg)

        WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a
      END DO ! icg    
    CASE ( GAS_MODEL_MIXT_TCPERF,GAS_MODEL_MIXT_PSEUDO ) 
      DO icg = icgBeg,icgEnd
        r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        u  = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
        p  = pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T  = pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a  = pRegion%mixt%dv(DV_MIXT_SOUN,icg)
        mw = pRegion%mixt%gv(GV_MIXT_MOL ,icg)
        cp = pRegion%mixt%gv(GV_MIXT_CP  ,icg)

        WRITE(IF_EXTR_DATA1,'(8(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a,mw,cp
      END DO ! icg        
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%gasModel

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extracting data for shocktube done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSod





! ******************************************************************************
!
! Purpose: Extract data for Sommerfeld shock-particle-interaction case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie along x-axis.
!   2. Assume cross-section to have 1 cell layer in z.
!   3. Assume patches 1 and 2 to be for x=xlow and x=xhigh boundaries.
!   4. Assumptions 2 and 3 mean that number of y-layers can be determined by 
!      number of faces on patch 1 or 2.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSommSPI(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgBeg,icgEnd,nCellLayersY,nCellsX
  REAL(RFREAL) :: a,cp,idx,ir,lx,mw,p,r,T,u,v,xs
  TYPE(t_grid), POINTER :: pGrid
#ifdef PLAG
  LOGICAL :: plagFlag
  INTEGER :: iLocTp,iLocUp
  REAL(RFREAL) :: Tp,up
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSommSPI', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extracting data for Sommerfeld SPI case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Find number of cell layers in y, cell indices, and grid spacing for extraction
! ******************************************************************************

  nCellLayersY = pRegion%patches(1)%nBFaces
  nCellsX      = pGrid%nCellsTot/nCellLayersY ! NOTE integer division

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,I5)') SOLVER_NAME,'Number of cell layers in y:', &
                                   nCellLayersY
    WRITE(STDOUT,'(A,3X,A,1X,I5)') SOLVER_NAME,'Number of cells in x:      ', &
                                   nCellsX
  END IF ! global%verbLevel
   
  icgBeg = (nCellLayersY-1)*nCellsX + 1
  icgEnd = nCellLayersY*nCellsX

  lx  = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nCellsTot)) &
      - MINVAL(pGrid%xyz(XCOORD,1:pGrid%nCellsTot))
  idx = REAL(nCellsX,KIND=RFREAL)/lx

#ifdef PLAG
  plagFlag = .FALSE.

  IF ( (global%plagUsed .EQV. .TRUE.) .AND. & 
       (global%postLag2EulFlag .EQV. .TRUE.) ) THEN 
    iLocUp = pRegion%plot%pv2pvi(PV_PLAG_XVEL)
    iLocTp = pRegion%plot%pv2pvi(PV_PLAG_TEMP) 

    IF ( (iLocUp /= CRAZY_VALUE_INT) .AND. &
         (iLocTp /= CRAZY_VALUE_INT) ) THEN 
      plagFlag = .TRUE.
    END IF ! iLocUp
  END IF ! global%plagUsed
#endif

! ******************************************************************************
! Write usual data
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  WRITE(iFileName1,'(A,1PE11.5)') 'somm_spi.dat_',global%currentTime  

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  
      
! ==============================================================================
! Extract usual data along line of cells 
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing solution data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel

  SELECT CASE ( pRegion%mixtInput%gasModel )
    CASE ( GAS_MODEL_TCPERF )             
      DO icg = icgBeg,icgEnd
        r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r

        u = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
        p =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)
#ifndef PLAG
        WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a                                                                                                    
#else
        IF ( plagFlag .EQV. .TRUE. ) THEN
          up = pRegion%plot%pv(iLocUp,icg)
          Tp = pRegion%plot%pv(iLocTp,icg)

          WRITE(IF_EXTR_DATA1,'(8(1X,E23.16))') pGrid%cofg(XCOORD,icg), &
                                                r,u,p,T,a,up,Tp
        ELSE 
          WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), &
                                                r,u,p,T,a
        END IF ! plagFlag
#endif
      END DO ! icg    
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%gasModel
  
! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error 
    
! ******************************************************************************
! Find shock position and write to file
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  WRITE(iFileName1,'(A,1PE11.5)') 'somm_spi_xs.dat_',global%currentTime

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error   

! ==============================================================================
! Extract data and write to file
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing shock position to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel
  
  CALL RFLU_ExtractShockLocation1D(pRegion,icgBeg,icgEnd,nCellsX,xs)

  WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') global%currentTime,xs

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error                                                                         

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Extracting data for Sommerfeld SPI case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSommSPI






! ******************************************************************************
!
! Purpose: Extract data for generic 1d shocktube case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie along x-axis.
!   2. Assume have 2 species if running with mixture gas models.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSTG1D(pRegion)

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgShock
  REAL(RFREAL) :: a,cp,ir,mw,p,r,T,u,v,xs
  TYPE(t_grid), POINTER :: pGrid
#ifdef SPEC
  INTEGER :: iSpec,iSpecEEv
  REAL(RFREAL) :: Y1,Y2  
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv   
  TYPE(t_spec_type), POINTER :: pSpecType  
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSTG1D', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for generic 1d shocktube case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Write usual data
! ******************************************************************************

! ==============================================================================
! Open file 
! ==============================================================================

  WRITE(iFileName1,'(A,1PE11.5)') 'stg1d.dat_',global%currentTime

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  
  
! ==============================================================================
! Extract usual data along line of cells 
! ==============================================================================

  SELECT CASE ( pRegion%mixtInput%gasModel )
    CASE ( GAS_MODEL_TCPERF ) 
      DO icg = 1,pGrid%nCellsTot
        r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r

        u = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
        p =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)

        WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a
      END DO ! icg    
    CASE ( GAS_MODEL_MIXT_TCPERF,GAS_MODEL_MIXT_PSEUDO ) 
      DO icg = 1,pGrid%nCellsTot
        r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r
 
        u  = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
        p  =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T  =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a  =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)
        mw =    pRegion%mixt%gv(GV_MIXT_MOL ,icg)
        cp =    pRegion%mixt%gv(GV_MIXT_CP  ,icg)
#ifdef SPEC
        Y1 = ir*pRegion%spec%cv(1,icg)
        Y2 = ir*pRegion%spec%cv(2,icg)

        WRITE(IF_EXTR_DATA1,'(10(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                               r,u,p,T,a,Y1,Y2,mw,cp
#endif
      END DO ! icg        
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%gasModel

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! Find shock position and write to file
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  WRITE(iFileName1,'(A,1PE11.5)') 'stg1d_xs.dat_',global%currentTime

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error   

! ==============================================================================
! Extract data and write to file
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing shock position to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel
  
  CALL RFLU_ExtractShockLocation1D(pRegion,1,pGrid%nCells,pGrid%nCells,xs)

  WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') global%currentTime,xs 

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error         

#ifdef SPEC  
! ******************************************************************************
! Write EE data
! ******************************************************************************

  IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 

! ==============================================================================
!   Loop over species
! ==============================================================================

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)
     
! ------------------------------------------------------------------------------
!     Write EE data if species evolved with EE approach
! ------------------------------------------------------------------------------

      IF ( pSpecType%velocityMethod == SPEC_METHV_EQEUL ) THEN
        iSpecEEv = pSpecType%iSpec2iSpecEEv
  
        pEEv => pRegion%spec%eev

! ----- Open file --------------------------------------------------------------

        WRITE(iFileName1,'(A,I2.2,A)') 'stg1d',iSpecEEv,'.dat'

        OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
             IOSTAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '// &
                         TRIM(iFileName1))
        END IF ! global%error  

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I2,1X,A)') SOLVER_NAME, & 
                'Writing eev data for species',iSpec,'to file: '// &
                TRIM(iFileName1)
        END IF ! global%verbLevel          
        
! ----- Write data -------------------------------------------------------------        
        
        DO icg = 1,pGrid%nCellsTot
          WRITE(IF_EXTR_DATA1,'(5(1X,E23.16))') & 
                pGrid%cofg(XCOORD,icg), & 
!                pEEv(EEV_SPEC_XVEL,iSpecEEv,icg), & 
!                pEEv(EEV_SPEC_YVEL,iSpecEEv,icg), &
!                pEEv(EEV_SPEC_ZVEL,iSpecEEv,icg), &  
!                pEEv(EEV_SPEC_TEMP,iSpecEEv,icg)
                pEEv(1,iSpecEEv,icg), & 
                pEEv(2,iSpecEEv,icg), &
                pEEv(3,iSpecEEv,icg), &  
                pEEv(4,iSpecEEv,icg)                
        END DO ! icg 

! ----- Close file -------------------------------------------------------------
        
        CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '// & 
                         TRIM(iFileName1))
        END IF ! global%error                    
      END IF ! pSpecType%velocityMethod
    END DO ! iSpecEE
  END IF ! pRegion%specInput%nSpeciesEE
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Extracting data for generic 1d shocktube case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSTG1D




! ******************************************************************************
!
! Purpose: Extract data for generic 2d shocktube case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie along x-axis.
!   2. Assume cross-section to have 3 cells, so can find number of cells
!      along x-axis by dividing total number of cells by 3.
!   3. Assume have 2 species if running with mixture gas models.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSTG2D(pRegion)

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgBeg,icgEnd,nCellsX
  REAL(RFREAL) :: a,cp,ir,mw,p,r,T,u,v
  TYPE(t_grid), POINTER :: pGrid
#ifdef SPEC
  INTEGER :: iSpec,iSpecEEv
  REAL(RFREAL) :: Y1,Y2  
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv   
  TYPE(t_spec_type), POINTER :: pSpecType  
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSTG2D', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for generic 2d shocktube case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Find cell indices for extraction
! ******************************************************************************

  nCellsX = pGrid%nCellsTot/3 ! NOTE integer division

  icgBeg =   nCellsX + 1
  icgEnd = 2*nCellsX

! ******************************************************************************
! Write usual data
! ******************************************************************************

! ==============================================================================
! Open file 
! ==============================================================================

  iFileName1 = 'stg2d.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  
  
! ==============================================================================
! Extract usual data along line of cells 
! ==============================================================================

  SELECT CASE ( pRegion%mixtInput%gasModel )
    CASE ( GAS_MODEL_TCPERF ) 
      DO icg = icgBeg,icgEnd
        r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r

        u = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
        p =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)

        WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a
      END DO ! icg    
    CASE ( GAS_MODEL_MIXT_TCPERF,GAS_MODEL_MIXT_PSEUDO ) 
      DO icg = icgBeg,icgEnd
        r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r
 
        u  = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
        p  =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T  =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a  =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)
        mw =    pRegion%mixt%gv(GV_MIXT_MOL ,icg)
        cp =    pRegion%mixt%gv(GV_MIXT_CP  ,icg)
#ifdef SPEC
        Y1 = ir*pRegion%spec%cv(1,icg)
        Y2 = ir*pRegion%spec%cv(2,icg)

        WRITE(IF_EXTR_DATA1,'(10(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                               r,u,p,T,a,Y1,Y2,mw,cp
#endif
      END DO ! icg        
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%gasModel

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

#ifdef SPEC  
! ******************************************************************************
! Write EE data
! ******************************************************************************

  IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 

! ==============================================================================
!   Loop over species
! ==============================================================================

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)
     
! ------------------------------------------------------------------------------
!     Write EE data if species evolved with EE approach
! ------------------------------------------------------------------------------

      IF ( pSpecType%velocityMethod == SPEC_METHV_EQEUL ) THEN
        iSpecEEv = pSpecType%iSpec2iSpecEEv
  
        pEEv => pRegion%spec%eev

! ----- Open file --------------------------------------------------------------

        WRITE(iFileName1,'(A,I2.2,A)') 'stg2d',iSpecEEv,'.dat'

        OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
             IOSTAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '// &
                         TRIM(iFileName1))
        END IF ! global%error  

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I2,1X,A)') SOLVER_NAME, & 
                'Writing eev data for species',iSpec,'to file: '// &
                TRIM(iFileName1)
        END IF ! global%verbLevel          
        
! ----- Write data -------------------------------------------------------------        
        
        DO icg = icgBeg,icgEnd
          WRITE(IF_EXTR_DATA1,'(5(1X,E23.16))') & 
                pGrid%cofg(XCOORD,icg), & 
!                pEEv(EEV_SPEC_XVEL,iSpecEEv,icg), & 
!                pEEv(EEV_SPEC_YVEL,iSpecEEv,icg), &
!                pEEv(EEV_SPEC_ZVEL,iSpecEEv,icg), &  
!                pEEv(EEV_SPEC_TEMP,iSpecEEv,icg)
                pEEv(1,iSpecEEv,icg), & 
                pEEv(2,iSpecEEv,icg), &
                pEEv(3,iSpecEEv,icg), &  
                pEEv(4,iSpecEEv,icg)                
        END DO ! icg 

! ----- Close file -------------------------------------------------------------
        
        CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '// & 
                         TRIM(iFileName1))
        END IF ! global%error                    
      END IF ! pSpecType%velocityMethod
    END DO ! iSpecEE
  END IF ! pRegion%specInput%nSpeciesEE
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Extracting data for generic 2d shocktube case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSTG2D







! ******************************************************************************
!
! Purpose: Extract data for Skews shock diffraction case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie in z=constant plane.
!   2. Assume domain to be divided into two pieces, a chamber (into which the
!      shock diffracts) and a tube (from which the shock emanates).
!   3. Assume that the chamber is bounded by (xMin,xMax) x (yMin,yMax) and tube
!      by (xMin,0) x (0,yMax).
!   4. Assume that grid is uniform and that cells in tube are numbered in
!      y-direction fastest. 
!   5. All these assumptions are satisfied by grids generated by gg_skews.f90
!      grid generator...
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSkews(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgOffs,ix,iy,nx,ny
  REAL(RFREAL) :: a,h,p,r,T,u,v,xMax,xMin,yMax,yMin
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSkews', &
                        'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extracting data for Skews case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file for data 
! ******************************************************************************

  iFileName1 = 'skews.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

! ******************************************************************************
! Find cell from which to start extracting
! ******************************************************************************

! ==============================================================================
! Compute spacing from coordinate extrema
! ==============================================================================

  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))

  h = SQRT((xMax*(yMax-yMin)-xMin*yMax)/DBLE(pGrid%nCells))

! ==============================================================================
! Compute number of cells in duct portion of domain 
! ==============================================================================

  nx = INT(-xMin/h+0.5_RFREAL)
  ny = INT( yMax/h+0.5_RFREAL)

! ==============================================================================
! Compute number of cells in chamber portion of domain and use as offset 
! ==============================================================================

  icgOffs = INT(xMax/h+0.5_RFREAL)*INT((yMax-yMin)/h+0.5_RFREAL)  

! ******************************************************************************
! Extract along line of cells halfway up duct portion of domain
! ******************************************************************************

  iy = nx/2 + 1

  DO ix = 1,nx
    icg = icgOffs + iy + (ix - 1)*ny

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg)
    T = pRegion%mixt%dv(DV_MIXT_TEMP,icg)
    a = pRegion%mixt%dv(DV_MIXT_SOUN,icg)

    WRITE(IF_EXTR_DATA1,'(7(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                          pGrid%cofg(YCOORD,icg), &
                                          r,u,p,T,a
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extracting data for Skews case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSkews








! ******************************************************************************
!
! Purpose: Extract Mesh boundaries for Bump case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_WriteMeshBump(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: iPatch,iv1,iv2
  INTEGER :: errorFlag,icg,ix
  REAL(RFREAL) :: a,p,r,u,v,w,M,xx,yy,angle
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_WriteMeshBump',&
  'RFLU_ModExtractFlowData.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract boundary grid points data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '_bdry','.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to  '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************

!  DO iPatch = 1,pGrid%nPatches
  DO iPatch = 1,4
    pPatch => pRegion%patches(iPatch)
     
    WRITE(IF_EXTR_DATA1,*) pPatch%nBFaces+1 
     
    IF ( iPatch < 3 ) THEN
      ix = 1
      iv1 = pPatch%bv(pPatch%bf2v(1,ix))
      xx = pGrid%xyz(XCOORD,iv1) 
      yy = pGrid%xyz(YCOORD,iv1) 
      WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') xx,yy

      DO ix = 1,pPatch%nBFaces
        iv2 = pPatch%bv(pPatch%bf2v(2,ix))
        xx = pGrid%xyz(XCOORD,iv2) 
        yy = pGrid%xyz(YCOORD,iv2) 
        WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') xx,yy
      END DO ! ix
    ELSE
      ix = pPatch%nBFaces
      iv1 = pPatch%bv(pPatch%bf2v(1,ix))
      xx = pGrid%xyz(XCOORD,iv1) 
      yy = pGrid%xyz(YCOORD,iv1) 
      WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') xx,yy

      DO ix = pPatch%nBFaces,1,-1
        iv2 = pPatch%bv(pPatch%bf2v(2,ix))
        xx = pGrid%xyz(XCOORD,iv2) 
        yy = pGrid%xyz(YCOORD,iv2) 
        WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') xx,yy
      END DO ! ix
    END IF
  END DO ! iPatch

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract boundary grid points for Bump case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WriteMeshBump







END MODULE RFLU_ModExtractFlowData

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModExtractFlowData.F90,v $
! Revision 1.24  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.23  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.22  2007/04/12 00:29:56  haselbac
! Updates to take into account changes to shock location routine
!
! Revision 1.21  2007/04/05 01:14:28  haselbac
! Added USE of new extract module, added xs extraction for stg1d
!
! Revision 1.20  2007/03/27 00:47:27  haselbac
! Added extraction of particle data to Sommerfeld case
!
! Revision 1.19  2007/03/19 21:43:28  haselbac
! Fixed typo in somm_spi write statement
!
! Revision 1.18  2007/03/02 17:56:33  haselbac
! Updated SommSPI data extraction to allow for 1d/2d cases
!
! Revision 1.17  2007/02/27 13:23:40  haselbac
! Added stg1d case
!
! Revision 1.16  2007/02/17 20:57:20  haselbac
! Added shock-position extraction to Sommerfeld case
!
! Revision 1.15  2007/02/16 20:01:29  haselbac
! Added code for somm_spi case
!
! Revision 1.14  2006/08/19 15:41:17  mparmar
! Added data extractions for nscbc[1-8], farf, bumpq10
!
! Revision 1.13  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.12  2006/01/04 15:48:54  haselbac
! Changed computation of thicknesses for lfp
!
! Revision 1.11  2005/11/27 02:00:22  haselbac
! Added extraction for EEv, bug fix
!
! Revision 1.10  2005/11/17 14:50:48  haselbac
! Added extraction of mass fractions for STG
!
! Revision 1.9  2005/11/14 17:05:24  haselbac
! Added support for pseudo-gas model
!
! Revision 1.8  2005/11/11 17:21:04  haselbac
! Added extraction for generic 2d shocktube
!
! Revision 1.7  2005/11/10 02:51:30  haselbac
! Added support for MP sod cases, clean-up
!
! Revision 1.6  2005/10/09 15:11:40  haselbac
! Added ONERA C0 data extraction
!
! Revision 1.5  2005/07/19 19:19:25  haselbac
! Added new lfp cases, generalized and fixed extraction for lfp
!
! Revision 1.4  2005/06/14 01:09:50  haselbac
! Adapted to new case names for Skews case
!
! Revision 1.3  2005/03/18 23:09:34  haselbac
! Added routine to extract data for Skews case
!
! Revision 1.2  2004/10/27 18:08:51  haselbac
! Cosmetics only
!
! Revision 1.1  2004/10/26 15:19:11  haselbac
! Initial revision
!
! ******************************************************************************

















