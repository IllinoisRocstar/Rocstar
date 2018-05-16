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
! Purpose: Suite of routines to compute exact flow solutions.
!
! Description: None.
!
! Notes: 
!   1. Collected routines in single module because computation of exact 
!      solutions is needed in at least three places: initialization of 
!      solution, setting of boundary profiles, and computation of errors. 
!
! ******************************************************************************
!
! $Id: RFLU_ModExactFlow.F90,v 1.10 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModExactFlow

  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
      
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_ComputeExactFlowCulick, & 
            RFLU_ComputeExactFlowPAcoust, &
            RFLU_ComputeExactFlowProudman, & 
            RFLU_ComputeExactFlowRingleb, & 
            RFLU_ComputeExactFlowSsVortex, & 
            RFLU_SetExactFlowLinear, & 
            RFLU_SetExactFlowTrig

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModExactFlow.F90,v $ $Revision: 1.10 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Compute exact solution for Culick flow.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   z		z-coordinate
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. Assume cylindrical domain with axis along x-coordinate direction.
!   2. Assume radius equal to 0.01.
!   3. Assume injection velocity equal to unity.
!   4. Total pressure is not uniform but assumed so for now.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowCulick(global,x,y,z,d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: x,y,z
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: r,ro,theta,vr

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    r     = SQRT(y*y + z*z)
    ro    = 0.01_RFREAL
    theta = ATAN2(z,y)

    d  = 1.0_RFREAL
    
    vr = -SIN(0.5_RFREAL*global%pi*(r/ro)**2)/(r/ro)
    u  = global%pi*x/ro*COS(0.5_RFREAL*global%pi*(r/ro)**2)    
    v  = vr*COS(theta)
    w  = vr*SIN(theta)                               

    p = 1.0E+5 ! TEMPORARY

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowCulick






! ******************************************************************************
!
! Purpose: Compute exact solution for pipe acoustics.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   z           z-coordinate
!   t           time
!   L           Length of pipe
!   ro          Outer radius of pipe
!   iBc         Type of boundary condition in axial direction (0 = zero
!               pressure disturbance, 1 = zero pressure disturbance gradient)
!   im          Index of mode in circumferential direction
!   in          Index of mode in axial direction
!   iq          Index of mode in radial direction
!   etaqm       q-th root of m-th order Bessel function of first kind
!   omega       Frequency
!   dTot        Total density
!   pTot        Total pressure
!   aTot        Total speed of sound
!   const       Constant (appears in theoretical solution, but is not equal 
!               to amplitude)
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. Exact solution assumes that z-axis corresponds to axis of symmetry of 
!      pipe.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowPAcoust(global,x,y,z,t,L,ro,iBc,im,in,iq, &
                                          etaqm,omega,dTot,pTot,aTot,const, & 
                                          d,u,v,w,p)

    USE RFLU_ModBessel

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iBc,im,in,iq
    REAL(RFREAL), INTENT(IN) :: aTot,const,dTot,etaqm,L,omega,pTot,ro,t,x,y,z
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: cimt,cinz,cot,dummyReal,Jm,Jmd,r,simt,sinz,sot,term, & 
                    theta,ur,ut,uz

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    r = SQRT(x*x + y*y)
    theta = ATAN2(y,x)

    CALL RFLU_JYNDD(im,etaqm*r/ro,Jm,Jmd,dummyReal,dummyReal,dummyReal, & 
                    dummyReal)
    
    cimt = COS(im*theta)
    simt = SIN(im*theta)
    
    cinz = COS(in*global%pi*z/L)
    sinz = SIN(in*global%pi*z/L)    
    
    cot  = COS(omega*t)
    sot  = SIN(omega*t)    
    
    term = const/(dTot*omega)
        
    IF ( iBc == 0 ) THEN     
      p = const*Jm*cimt*sinz*cot
    ELSE IF ( iBc == 1 ) THEN 
      p = const*Jm*cimt*cinz*cot
    ELSE 
      p = 0.0_RFREAL
    END IF ! iBc
    
    d = dTot + p/aTot*aTot
    p = pTot + p

    IF ( iBc == 0 ) THEN   
      ur = -   term  *Jmd*(etaqm/ro)*cimt*sinz                 *sot
      ut =  im*term/r*Jm            *simt*sinz                 *sot
      uz = -   term  *Jm            *cimt*cinz*(in*global%pi/L)*sot    
    ELSE IF ( iBc == 1 ) THEN
      ur = -   term  *Jmd*(etaqm/ro)*cimt*cinz                 *sot
      ut =  im*term/r*Jm            *simt*cinz                 *sot
      uz =     term  *Jm            *cimt*sinz*(in*global%pi/L)*sot
    ELSE 
      ur = 0.0_RFREAL
      ut = 0.0_RFREAL
      uz = 0.0_RFREAL      
    END IF ! iBc

    u = ur*COS(theta) - ut*SIN(theta)
    v = ur*SIN(theta) + ut*COS(theta)
    w = uz

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowPAcoust






! ******************************************************************************
!
! Purpose: Compute exact solution for Proudman-Culick flow.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   height      Height of domain
!   dInc        Density
!   vInj        Injection velocity
!   pTot        Total pressure at (x,y) = (0,0)
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. Total pressure is not uniform!
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowProudman(global,x,y,height,dInc,vInj,pTot, &
                                           d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: dInc,height,pTot,vInj,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    u = -0.5_RFREAL*global%pi*x/height*vInj*COS(0.5_RFREAL*global%pi*y/height)
    v =                                vInj*SIN(0.5_RFREAL*global%pi*y/height)

    w = 0.0_RFREAL

    d = dInc
    p = pTot - 0.25_RFREAL*vInj**2*(1.0_RFREAL &
                                  + 0.5_RFREAL*(global%pi*x/height)**2 & 
                                  -         COS(global%pi*y/height))

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowProudman






! ******************************************************************************
!
! Purpose: Compute exact solution for Ringleb flow.
!
! Description: The solution to the Ringleb flow is given in terms of qBar 
!   (normalized velocity magnitude) and the k (streamline constant). It is
!   not possible to determine these as an explicit function of x and y. 
!   Instead, given x and y, qBar and k are determined iteratively. Once 
!   qBar and k are known, the solution can be computed in terms of the 
!   primitive variables. 
!
! Input: 
!   x           x-coordinate
!   y           y-coordinate
!   rGas        Gas constant
!   pTot        Total pressure
!   tTot        Total temperature
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. This routine assumes a perfect gas and that the ratio of specific heats
!      is equal to 1.4.
!   2. Depending on how the geometry is defined (which values of qBar and k 
!      were chosen), it may become necessary to alter the initial guess for 
!      qBar.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowRingleb(x,y,rGas,pTot,tTot,d,u,v,w,p)

    USE ModInterfaces, ONLY: MixtPerf_C_GRT,MixtPerf_D_PRT

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: pTot,rGas,tTot,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w

! ==============================================================================
! Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: RCSIdentString
    INTEGER :: cntr,cntrMax
    REAL(RFREAL) :: a,aBar,aTot,alpha,daBardqBar,dFdqBar,dJdaBar,dJdqBar, &
                    dqBar,drhoBardqBar,F,J,k,qBar,qBarMax,rhoBar,t,tol

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

! ==============================================================================
!   Determine solution in terms of qBar and k with Newton-Raphson method
! ==============================================================================

    cntr    = 0
    cntrMax = 100 ! Maximum number of Newton-Raphson steps

    qBar    = 0.3_RFREAL       ! Initial guess for qBar
    qBarMax = SQRT(5.0_RFREAL) ! Physically possible maximum value for qBar 

    tol = 1.0E-16_RFREAL ! Convergence tolerance

! -------------------------------------------------------------------------------
!   Loop over iterations
! -------------------------------------------------------------------------------

    DO 
      cntr = cntr + 1

      aBar   = SQRT(1.0_RFREAL - 0.2_RFREAL*qBar**2)                        
      rhoBar = aBar**5.0_RFREAL

      J      = 1.0_RFREAL/aBar                 & 
             + 1.0_RFREAL/(3.0_RFREAL*aBar**3) & 
             + 1.0_RFREAL/(5.0_RFREAL*aBar**5) &
             - 0.5_RFREAL*LOG((1.0_RFREAL+aBar)/(1.0_RFREAL-aBar))              

! --- Compute function value ---------------------------------------------------                

      F = (x - 0.5_RFREAL*J)**2 + y**2 - 1.0_RFREAL/(4.0_RFREAL*rhoBar**2*qBar**4)                

! --- Compute derivative -------------------------------------------------------        

      drhoBardqBar = -qBar*aBar**3      
      dJdaBar = -1.0_RFREAL/aBar**2 & 
                -1.0_RFREAL/aBar**4 & 
                -1.0_RFREAL/aBar**6 &
                -1.0_RFREAL/((1.0_RFREAL+aBar)*(1.0_RFREAL-aBar))
      daBardqBar = -0.2_RFREAL*qBar/aBar
      dJdqBar = dJdaBar*daBardqBar        
      dFdqBar = -(x - 0.5_RFREAL*J)*dJdqBar & 
              + 0.5_RFREAL/(rhoBar**3*qBar**4)*drhoBardqBar & 
              + 1.0_RFREAL/(rhoBar**2*qBar**5)

! --- Update -------------------------------------------------------------------        

      dqBar = -F/dFdqBar         

      qBar = qBar + dqBar                                                                                 

! --- Limit to physically possible maximum -------------------------------------

      qBar = MIN(qBar,qBarMax)                   

! --- Convergence check and update k -------------------------------------------        

      IF ( (ABS(dqBar) < tol ) .OR. (cntr == cntrMax) ) THEN
        k = SQRT(2.0_RFREAL/(1.0_RFREAL/qBar**2 & 
                           - 2.0_RFREAL*rhoBar*(x - 0.5_RFREAL*J)))

        EXIT
      END IF ! ABS(dqBar)                   
    END DO ! <empty>

! ==============================================================================
!   Determine solution in terms of primitive variables
! ==============================================================================

    p    = pTot*rhoBar**1.4_RFREAL
    t    = tTot*rhoBar**0.4_RFREAL          
    d    = MixtPerf_D_PRT(p,rGas,t)      
    aTot = MixtPerf_C_GRT(1.4_RFREAL,rGas,tTot)
    a    = aTot*aBar

    alpha = ASIN(MIN(1.0_RFREAL,qBar/k))  
    u     = aTot*qBar*COS(alpha)*SIGN(1.0_RFREAL,y)
    v     = aTot*qBar*SIN(alpha) 
    w     = 0.0_RFREAL

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowRingleb








! ******************************************************************************
!
! Purpose: Compute exact solution for supersonic vortex.
!
! Description: None.
!
! Input: 
!   x           x-coordinate
!   y           y-coordinate
!   gGas        ratio of specific heats
!   rGas        Gas constant
!   ri          Inner radius
!   Mi          Mach number at inner radius
!   pTot        Total pressure
!   tTot        Total temperature
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. This routine assumes a perfect gas.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowSsVortex(x,y,gGas,rGas,ri,Mi,pTot,tTot, &
                                           d,u,v,w,p)

    USE ModInterfaces, ONLY: MixtPerf_C_GRT, & 
                             MixtPerf_D_DoGMa, &
                             MixtPerf_D_PRT, & 
                             MixtPerf_P_GMaPo, &
                             MixtPerf_P_DDoGPo, &
                             MixtPerf_T_DPR 

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: gGas,Mi,pTot,rGas,ri,tTot,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: ai,alpha,di,dTot,pi,r,term,ti,Vi

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    dTot = MixtPerf_D_PRT(pTot,rGas,tTot)

    di = MixtPerf_D_DoGMa(dTot,gGas,Mi)      
    pi = MixtPerf_P_GMaPo(gGas,Mi,pTot)
    ti = MixtPerf_T_DPR(di,pi,rGas) 
    ai = MixtPerf_C_GRT(gGas,rGas,ti)
    Vi = Mi*ai

    r     = SQRT(x*x + y*y)
    alpha = ATAN(y/x) 

    term = 1.0_RFREAL - (ri/r)**2
    term = 1.0_RFREAL + 0.5_RFREAL*(gGas - 1.0_RFREAL)*Mi*Mi*term
    d    = di*term**(1.0_RFREAL/(gGas - 1.0_RFREAL))
    p    = MixtPerf_P_DDoGPo(d,dTot,gGas,pTot)

    u =  Vi*ri/r*SIN(alpha)
    v = -Vi*ri/r*COS(alpha)
    w =  0.0_RFREAL

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowSsVortex







! ******************************************************************************
!
! Purpose: Set linear variable behavior.
!
! Description: None.
!
! Input: 
!   x           x-coordinate
!   y           y-coordinate
!   z           z-coordinate
!   iVar        Variable index
!
! Output: 
!   var         Variable
!   gx          x-component of gradient
!   gy          y-component of gradient
!   gz          z-component of gradient
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetExactFlowLinear(x,y,z,iVar,var,gx,gy,gz)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar
    REAL(RFREAL), INTENT(IN) :: x,y,z
    REAL(RFREAL), INTENT(OUT) :: gx,gy,gz,var

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: a,b,c

! ******************************************************************************
!   Start, set exact solution
! ******************************************************************************

    a = 1 + (iVar-1)*(ZCOORD-XCOORD)
    b = a + 1
    c = a + 2

    var = a*x + b*y + c*z
    
    gx = a
    gy = b
    gz = c 

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_SetExactFlowLinear








! ******************************************************************************
!
! Purpose: Set trigonometric variable behavior.
!
! Description: None.
!
! Input: 
!   global      Pointer to globa data
!   nx          Wave number for x-direction
!   ny          Wave number for y-direction
!   nz          Wave number for z-direction
!   x           x-coordinate
!   y           y-coordinate
!   z           z-coordinate
!   iVar        Variable index
!
! Output: 
!   var         Variable
!   gx          x-component of gradient
!   gy          y-component of gradient
!   gz          z-component of gradient
!
! Notes: 
!   1. At present, do not make use of iVar.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,iVar,var,gx,gy,gz)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar
    REAL(RFREAL), INTENT(IN) :: nx,ny,nz,x,y,z
    REAL(RFREAL), INTENT(OUT) :: gx,gy,gz,var
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: a,b,c

! ******************************************************************************
!   Start, set exact solution
! ******************************************************************************

    a = nx*global%pi
    b = ny*global%pi
    c = nz*global%pi

    var = COS(a*x)*SIN(b*y)*COS(c*z)
    
    gx = -a*SIN(a*x)*SIN(b*y)*COS(c*z)
    gy =  b*COS(a*x)*COS(b*y)*COS(c*z)
    gz = -c*COS(a*x)*SIN(b*y)*SIN(c*z) 

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_SetExactFlowTrig






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModExactFlow


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModExactFlow.F90,v $
! Revision 1.10  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/04/13 18:07:34  haselbac
! Added routine for Culick flow
!
! Revision 1.7  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.6  2006/01/06 22:10:26  haselbac
! Added routines for linear and trig solution for grad testing
!
! Revision 1.5  2005/04/29 12:54:36  haselbac
! Adapted routine for exact pipe acoust solution to accept time argument
!
! Revision 1.4  2005/04/20 14:41:53  haselbac
! Bug fix and extensions for pipe acoustics case
!
! Revision 1.3  2005/03/31 17:02:44  haselbac
! Fixed bug in initialization of pressure for ONERA C0 case
!
! Revision 1.2  2005/03/15 20:44:44  haselbac
! Added routine to compute exact solution for pipe acoustics
!
! Revision 1.1  2004/07/06 15:14:27  haselbac
! Initial revision
!
! ******************************************************************************






