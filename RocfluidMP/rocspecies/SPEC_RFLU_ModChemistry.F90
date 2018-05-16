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
! Purpose: Simple chemistry model for burning crack application.
!
! Description: None.
!
! Notes: 
!   1. Based on collection of routines developed by Luca Massa.
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_ModChemistry.F90,v 1.6 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE SPEC_RFLU_ModChemistry

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
    
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: SPEC_RFLU_IntegrateChemSrcTerm

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  INTEGER, PARAMETER :: MAXEQN = 2,MAXRCT = 1
  INTEGER :: nreaction,ntemp_spec,maxnwtn
  REAL(RFREAL) :: eps,epy,beta,cpGas
  REAL(RFREAL) :: damkholer(MAXRCT),mx(MAXEQN),mn(MAXEQN),nn(MAXRCT), &
                  pre(MAXRCT),qg(MAXRCT),thetay(MAXRCT),VEC(MAXEQN)
  REAL(RFREAL) :: MAT(MAXEQN,MAXEQN),yf(MAXEQN,MAXRCT)
  REAL(RFREAL) :: ratey(MAXRCT)    

! ******************************************************************************
! Routines
! ******************************************************************************
  
  CONTAINS
  



! ******************************************************************************
!
! Purpose: Integrate source term for mixture and species using ODE solver.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Numerical method is due to Luca Massa and Mark Short.
!
! ******************************************************************************
  
  SUBROUTINE SPEC_RFLU_IntegrateChemSrcTerm(pRegion,CALLFLAG)
  
    USE ModInterfaces, ONLY: MixtPerf_Eo_GRTUVW, &
                             MixtPerf_G_CpR, &
                             MixtPerf_R_M, &
                             MixtPerf_D_PRT

  
! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: errorString,RCSIdentString
    INTEGER :: icg,indCp,indMol,iSpec,nout,ncumul,ibc,iec,CALLFLAG
    REAL(RFREAL) :: fact,gGas,ir,molMass,p,r,rGas,t,u,v,w,rvcumul,x,y,aaa,bbb
    REAL(RFREAL), DIMENSION(MAXEQN) :: f,f_old
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec,pDvMixt,pGvMixt,Mrhs,Srhs
    REAL(RFREAL), DIMENSION(:), POINTER :: vol
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    RCSIdentString = '$RCSfile: SPEC_RFLU_ModChemistry.F90,v $ $Revision: 1.6 $'

    global => pRegion%global

    CALL RegisterFunction(global,'SPEC_RFLU_IntegrateChemSrcTerm',&
  'SPEC_RFLU_ModChemistry.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
    
    pCvMixt => pRegion%mixt%cv
    pDvMixt => pRegion%mixt%dv
    pGvMixt => pRegion%mixt%gv    
    pCvSpec => pRegion%spec%cv    
    Mrhs    => pRegion%mixt%rhs   
    Srhs    => pRegion%spec%rhs
    vol     => pRegion%grid%vol

    indCp  = pRegion%mixtInput%indCp
    indMol = pRegion%mixtInput%indMol

    fact = 1.0_RFREAL
    nout = 100

    CALL set_rate(pRegion) ! Set rate constants

! ******************************************************************************
!   Checks 
! ******************************************************************************

    IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
      CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
    END IF ! pRegion%mixt%cvState

    IF ( pRegion%specInput%nSpecies /= MAXEQN-1 ) THEN
      WRITE(errorString,'(I3)') MAXEQN 
      CALL ErrorStop(global,ERR_SPEC_MAXEQN,__LINE__,errorString)
    END IF ! pRegion%specInput%nSpecies
    
! ******************************************************************************
!   Loop limits and sum initialization
! ******************************************************************************

    ibc = 1
    iec = merge(pGrid%nCellsTot,pGrid%nCells,CALLFLAG == 0)
    ncumul = 0
    rvcumul = 0.0_RFREAL

! ******************************************************************************
!   Integrate source term
! ******************************************************************************
  
    DO icg =  ibc,iec
      r =  pCvMixt(CV_MIXT_DENS,icg)  
      ir = 1.0_RFREAL/r   
          
      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)        
      t = pDvMixt(DV_MIXT_TEMP,icg)         
      p = pDvMixt(DV_MIXT_PRES,icg)
       
      molMass = pGvMixt(GV_MIXT_MOL,indMol*icg)
      cpGas   = pGvMixt(GV_MIXT_CP ,indCp *icg) 
      
      rGas = MixtPerf_R_M(molMass)
      gGas = MixtPerf_G_CpR(cpGas,rGas) 
      
! ==============================================================================
!     Call integration routine
! ==============================================================================            
      
      !CALL rate_split(f,f_old,r,p,fact,global%dtMin,nout)
      pre(1) = damkholer(1)*p**nn(1)

! ==============================================================================
!     Update state vectors. NOTE have conserved state vector for the mixture, 
!     so need to update total internal energy because it is affected by the 
!     temperature. The density is assumed to stay constant during this update. 
!     The dependent variables will be updated outside this routine.
! ==============================================================================  

      if(CALLFLAG == 1) then
         ratey(1) = pre(1) * pCvSpec(1,icg) * EXP(-thetay(1)/t) / r * global%dtMin
         pCvMixt(CV_MIXT_ENER,icg) = pCvMixt(CV_MIXT_ENER,icg) + ratey(1)*yf(1,1)
         pCvSpec(           1,icg) = pCvSpec(           1,icg) + ratey(1)*yf(2,1)
      elseif(CALLFLAG == 0) then
         ratey(1) = pre(1) * pCvSpec(1,icg) * EXP(-thetay(1)/t) / r * vol(icg)

         Mrhs(CV_MIXT_ENER,icg) = Mrhs(CV_MIXT_ENER,icg) - ratey(1)*yf(1,1)
         Srhs(1,icg) = Srhs(1,icg) - ratey(1)*yf(2,1) 
!         if (abs(pRegion%grid%cofg(YCOORD,icg)) > 4.5200d-3) Mrhs(CV_MIXT_ENER,icg) = 0_RFREAL
!         if (abs(pRegion%grid%cofg(YCOORD,icg)) > 4.5200d-3) Srhs(:,icg) = 0_RFREAL
      endif
              
    END DO ! icg
      
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  
  
  END SUBROUTINE SPEC_RFLU_IntegrateChemSrcTerm
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Set rate constants. 
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! *******************************************************************************
  
  SUBROUTINE set_rate(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: n

! ******************************************************************************
!   Start
! ******************************************************************************

    maxnwtn = 1
    eps = 1.0E-8_RFREAL
    epy = 2.0E-4_RFREAL

    nreaction = 1
    ntemp_spec = 1+pRegion%specInput%nSpecies

    beta = 1.0_RFREAL

    thetay(1) =  15000_RFREAL 

    qg(1) = 2511600.0_RFREAL

    nn(1) = 0.0_RFREAL

    damkholer(1) = 1174775.7803_RFREAL


!!!one reaction, one species
    yf(1,1) = qg(1)

    yf(2,1) = -1.0_RFREAL

!TEMPERATURE AND SPECIES BOUNDS !IN REACTION STEP INTEGRATION

    mx(1) =   6000.0_RFREAL
    mn(1) = 0.0_RFREAL
    
    DO n = 2,ntemp_spec
       mn(n) = -10.0_RFREAL
       mx(n) = 1.0_RFREAL - mn(n)
    END DO ! n

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE set_rate
  
  
  
  
! ******************************************************************************
!
! Purpose: Integrate source terms using ODE solver.
!
! Description: None.
!
! Input: 
!   old_f
!   rho
!   time_step
!   nout
!
! Output: 
!   f 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE rate_split(f,old_f,rho,press,fact,time_step,nout)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: fact,press,rho,time_step
    REAL(RFREAL), DIMENSION(MAXEQN) :: old_f,f

! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: outcome
    INTEGER :: eqn,i,icnvrgd,nout
    REAL(RFREAL) :: toll,dt_over_rho

! ******************************************************************************
!   Start
! ******************************************************************************

    dt_over_rho = time_step/rho*fact
 
    DO eqn = 1, ntemp_spec
      old_f(eqn) = f(eqn) 
    END DO ! eqn
    
! ******************************************************************************
!   Integrate
! ******************************************************************************
       
    CALL ODESOLVE(old_f(1:ntemp_spec),f(1:ntemp_spec),press, &
                  ntemp_spec,dt_over_rho,nout,outcome)

    IF ( outcome ) THEN
      WRITE(*,*)'DIVERGENCE IN RATE',old_f(1:ntemp_spec),f(1:ntemp_spec), &
                                     dt_over_rho
      f(1:ntemp_spec)=old_f(1:ntemp_spec)
      icnvrgd = 0
    ELSE
      icnvrgd = 1
    END IF ! outcome

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE rate_split




! ******************************************************************************
!
! Purpose: None.
!
! Description: None.
!
! Input: 
!   yin         Old solution
!   press       Pressure
!   neqin       Number of equations
!   dtout       Time step
!   nout        Number of subiterations
!
! Output: 
!   y           New solution
!   outcome     Success flag (TRUE if diverged!!!)
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE ODESOLVE(yin,y,press,neqin,dtout,nout,outcome)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    LOGICAL :: outcome
    INTEGER, INTENT(IN) :: neqin,nout
    REAL(RFREAL) :: dtout,press
    REAL(RFREAL) :: y(neqin),yin(neqin)

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icol,inwtn,iout,irow
    REAL(RFREAL) :: t,tout
    REAL(RFREAL) :: pp(MAXEQN),ytmp(MAXEQN)

! ******************************************************************************
!   Start
! ******************************************************************************

    dtout = dtout/REAL(nout,KIND=RFREAL)

    t    = 0.0_RFREAL
    tout = t
    y    = yin
    ytmp = yin
    pp   = yin

    DO iout = 1,nout
      DO inwtn = 1, maxnwtn
        call rate_stiff(neqin, tout, y, press, VEC)
        call drate_stiff(neqin, tout, y, press, 1,1,MAT,neqin)

        DO irow = 1,neqin
          VEC(irow) = ytmp(irow) - y(irow) + VEC(irow)*dtout
          
          DO icol = 1,neqin
            IF ( icol == irow ) THEN
              MAT(irow,irow) =  MAT(irow,irow)*dtout - 1.0_RFREAL
            ELSE
              MAT(irow,icol) = MAT(irow,icol)*dtout
            END IF ! icol
          END DO ! icol
        END DO ! irow 
        
        CALL LUSOLVE(MAT(1:neqin,1:neqin),VEC(1:neqin),neqin)
        y(1:neqin) = y(1:neqin) - VEC(1:neqin)
      END DO ! inwtn

      tout = tout + dtout
      ytmp = y

!      WRITE(23,*)iout,tout,y(1)
!      WRITE(6,*)iout,tout,y(1)
    END DO ! iout
    
    outcome = (.NOT. y(1) < mx(1)) .OR. (.NOT. y(1) > mn(1))

    IF ( outcome ) THEN
       DO irow = 1,neqin
          WRITE (STDOUT,*)irow,'TXY',yin(irow),y(irow)
       END DO ! irow
    END IF ! outcome

300 format(////' STATE FLAG =',i3,/,'FAILED INEGRATION')

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE ODESOLVE  
  
  
  
  

! ******************************************************************************
!
! Purpose: Set up rate for two step kinetics and for fine AP/binder mixture
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
! ***************************************************
!     ff(1)  = T    (gas phase temperature)
!     ff(2)  = F    (gas phase fuel Y)
!     ff(3)  = OX   (gas phase oxidizer X)
!     ff(4)  = Z    (gas phase intermediate Z)
! ********************************************************************

  SUBROUTINE rate_stiff(neqin,time_stiff,ff,press,rrate)
   
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: neqin
    REAL(RFREAL) :: press,time_stiff
    REAL(RFREAL) :: ff(neqin),rrate(neqin)

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: m,n
    REAL(RFREAL) :: ratey(MAXRCT)     
    REAL(RFREAL) :: dratey(MAXRCT,MAXEQN)  

! ******************************************************************************
!   Start
! ******************************************************************************

    pre(1) = damkholer(1)*press**nn(1)

    ratey(1) = ff(2)*EXP(-thetay(1)/ff(1))        

    ratey(1) = pre(1)*ratey(1)

!   yf(species,reaction), m species_eq, n=species_var

    DO m = 1,neqin     
      rrate(m) = SUM(yf(m,1:nreaction)*ratey(1:nreaction))
    END DO ! m

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE rate_stiff






! ******************************************************************************
!
! Purpose: None.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE drate_stiff(neqin,time_stiff,ff,press,ml,imu,drate,nrowpd)
   
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: neqin,ml,imu
    REAL(RFREAL) :: press,time_stiff
    REAL(RFREAL) :: ff(neqin),drate(nrowpd,neqin)
    REAL(RFREAL) :: dratey(MAXRCT,MAXEQN)

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: m,n,nrowpd
    REAL(RFREAL) :: tsqinv
    REAL(RFREAL) :: ratey(MAXRCT)

! ******************************************************************************
!   Start
! ******************************************************************************

    pre(1) = damkholer(1)*press**nn(1)

    tsqinv = 1.0_RFREAL/(ff(1)*ff(1))

    ratey(1) = ff(2)*EXP(-thetay(1)/ff(1))        

    ratey(1) = pre(1)*ratey(1)

    dratey(1,1) = thetay(1)*ratey(1)*tsqinv
    dratey(1,2) = ratey(1)/MERGE(eps,ff(2),ff(2)==0.0_RFREAL)

!   yf(species,reaction), m species_eq, n=species_var

    DO m = 1,neqin     
       DO n = 1,ntemp_spec
          drate(m,n) = SUM(yf(m,1:nreaction)*dratey(1:nreaction,n))
       END DO ! n
    END DO ! m

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE drate_stiff





! ******************************************************************************
!
! Purpose: Solve system of linear equations Ax = f by LU decomposition. 
!
! Description: None.
!
! Input: 
!   A           Matrix of size n by n
!   f           Vector of size n
!   n           Size of matrix and vector
!
! Output: 
!   f           Solution to Ax=f
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE LUSOLVE(A,f,n)

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: n
    REAL(RFREAL) :: f(n)
    REAL(RFREAL) :: A(n,n)

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: i,j
    REAL(RFREAL) :: d
    REAL(RFREAL), DIMENSION(MAXEQN) :: b,y
    REAL(RFREAL), DIMENSION(MAXEQN,MAXEQN) :: u,l

! ******************************************************************************
!   Start
! ******************************************************************************

    l = 1.0_RFREAL
    u = 0.0_RFREAL
    
    u(1,1)   = a(1,1)
    u(1,2:n) = a(1,2:n)
    l(2:n,1) = a(2:n,1)/a(1,1)

    DO i = 2,n-1
       u(i,i) = a(i,i) - SUM(l(i,1:i-1)*u(1:i-1,i))
       d = 1.0_RFREAL/u(i,i)

       DO j = i+1,n
          u(i,j) =  a(i,j) - SUM(l(i,1:i-1)*u(1:i-1,j))
          l(j,i) = (a(j,i) - SUM(l(j,1:i-1)*u(1:i-1,i)))*d
       END DO ! j
    END DO ! i

    i=n
    u(i,i) = a(i,i) - SUM(l(i,1:i-1)*u(1:i-1,i))

    y(1) = f(1)

    DO i = 2,n
       y(i) = f(i) - SUM(l(i,1:i-1)*y(1:i-1))
    END DO ! i

    f(n) = y(n)/u(n,n)
    DO i = n-1,1,-1
       f(i) = (y(i) - SUM(u(i,i+1:n)*f(i+1:n)))/u(i,i)
    END DO ! i

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE LUSOLVE




 
END MODULE SPEC_RFLU_ModChemistry

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_ModChemistry.F90,v $
! Revision 1.6  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.3  2005/06/08 16:06:23  haselbac
! Bug fix: Correcting incomplete check-in
!
! Revision 1.2  2005/06/06 14:24:11  haselbac
! Adapted to Lucas changes
!
! Revision 1.1  2004/04/01 21:22:41  haselbac
! Initial revision
!
! ******************************************************************************








