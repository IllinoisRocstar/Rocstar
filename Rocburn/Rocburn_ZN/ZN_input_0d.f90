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
!                                                                             
! ---------------------------------------------------------------------------  
!                                                                              
!   SUBROUTINE  : ZN_input_0d                                                  
!                                                                              
!   This subroutine read inputs for burning rate model Rocburn_1D_ZN
!
!   Authors          :  K. Tang
!
!   Creation Date    :  Sep. 3, 2002
!   
!   Modifications    :
!
!    No.     Date         Authors       Description
!
!                                                                              
! ---------------------------------------------------------------------------  
!                                                                              
!                                                                              
!   arguments   :                                                              
!                                                                              
!      G_ZN     : Global variables for Rocburn_1D_ZN
!      Indir     : directory for input data file
!                                                                              
! ---------------------------------------------------------------------------  
!
  SUBROUTINE ZN_input_0d(G_ZN, Indir)

    USE M_Rocburn_ZN_Global_Data

    IMPLICIT NONE
    INCLUDE 'mpif.h'

!
!   Global data for Rocburn_1D_ZN passed as a pointer
!

    TYPE (G_BURN_1D), POINTER :: G_ZN
!
!   arguments
!
    CHARACTER(*), INTENT(IN)  :: Indir

!
!
! ----------------------------------------------------------------
!   local variables
    CHARACTER(LEN=80) :: Infile
    CHARACTER(*), PARAMETER :: ControlFile = "RocburnZNControl.txt"

    INTEGER   :: ir, ioerr

!
!       read propellant thermophysical properties
!

    ir = 10


    if (Indir(LEN_TRIM(Indir):LEN_TRIM(Indir)) == '/') then
       Infile= TRIM(Indir) // ControlFile
    else
       Infile= TRIM(Indir) // ControlFile
    endif

    OPEN (unit=ir,file=Infile,status='old')

    READ(ir,*)  G_ZN%Model_combustion
    IF (G_ZN%rank==0) PRINT *,' ROCBURN_ZN: rank=',G_ZN%rank,          &
               ' ; input propellant thermophysical properties', &
                   G_ZN%Model_combustion
!    
!      Model_combustion = 1  :: WSB homogeneouse propellant 
!                               combustion model
!      Model_combustion = 2  :: ZN phenomenological combustion model 
!                               for composite and homogeneous 
!                               propellant
!      Model_combustion = 3  :: rb=a*P**n
!    

    READ(ir,*)   G_ZN%a_p
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  a_p =', G_ZN%a_p

    READ(ir,*)   G_ZN%n_p
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  n_p =', G_ZN%n_p

    READ(ir,*)   G_ZN%Ac
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  Ac= ',G_ZN%Ac

    READ(ir,*)   G_ZN%Bg
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  Bg= ',G_ZN%Bg

    READ(ir,*)   G_ZN%Ec
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  Ec= ',G_ZN%Ec

    READ(ir,*)   G_ZN%Qc
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  Qc= ',G_ZN%Qc

    READ(ir,*)   G_ZN%Qg
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  Qg= ',G_ZN%Qg

    READ(ir,*)   G_ZN%alfac
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  alfac= ',G_ZN%alfac

    READ(ir,*)   G_ZN%C
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  C= ',G_ZN%C

    READ(ir,*)   G_ZN%rhoc
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  rhoc= ',G_ZN%rhoc

    READ(ir,*)   G_ZN%lamg
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  lamg= ',G_ZN%lamg

    READ(ir,*)   G_ZN%MW
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  MW= ',G_ZN%MW

    READ(ir,*)   G_ZN%Ka
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  Ka= ',G_ZN%Ka

    READ(ir,*)   G_ZN%nxmax
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  nxmax= ',G_ZN%nxmax

    READ(ir,*)   G_ZN%nx
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  nx= ',G_ZN%nx

    READ(ir,*)   G_ZN%delt_max
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  delt= ',G_ZN%delt_max

    READ(ir,*)   G_ZN%igrid
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  igrid= ',G_ZN%igrid

    READ(ir,*)   G_ZN%xmax
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  xmax= ',G_ZN%xmax

    READ(ir,*)   G_ZN%beta
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  beta= ',G_ZN%beta

    READ(ir,*)   G_ZN%tol_Ts
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  tol_Ts= ',G_ZN%tol_Ts

    READ(ir,*)   G_ZN%itermax
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  itermax= ',G_ZN%itermax

    READ(ir,*)   G_ZN%Tf_adiabatic
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  Tf_adiabatic= ',G_ZN%Tf_adiabatic

    READ(ir,*)   G_ZN%To
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  To= ',G_ZN%To

    READ(ir,*)   G_ZN%ign_flag
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  ign_flag= ',G_ZN%ign_flag

    READ(ir,*)   G_ZN%To_cold
    IF (G_ZN%rank==0) WRITE(*,*) 'ROCBURN_ZN:  To_cold= ',G_ZN%To_cold

    CLOSE(ir)

    G_ZN%lamc=G_ZN%alfac*G_ZN%rhoc*G_ZN%C

!
!   initialize variables not allowed to be changed
!
   
    G_ZN%R=1.9872
    G_ZN%a_T=1.0
    G_ZN%n_T=1.0

    IF (G_ZN%rank==0) PRINT *,'ROCBURN_ZN: rank=',G_ZN%rank,       &
            ' done input', ControlFile

    RETURN 

  END SUBROUTINE ZN_input_0d






