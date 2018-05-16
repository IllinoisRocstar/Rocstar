! * *******************************************************************
! * Rocstar Simulation Suite                                          *
! * Version: Sandia Evaluation Distribution                           *
! * Licensed To: Sandia National Laboratories                         *
! * License Type: Evaluation                                          *
! * License Expiration: March 13, 2013                                *
!*********************************************************************
! * *******************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2012, IllinoisRocstar LLC. All rights reserved.         *
! *                                                                   *
! * The Rocstar Simulation Suite is the property of IllinoisRocstar   *
! * LLC. No use or distribution of this version of the Rocstar        *
! * Simulation Suite beyond the license provided through separate     *
! * contract is permitted.                                            *
! *                                                                   *
! * IllinoisRocstar LLC                                               *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *********************************************************************
! * *******************************************************************
! *  Initial open source Rocstar software developed by                *
!*     Center for Simulation of Advanced Rockets                     *
! *     University of Illinois at Urbana-Champaign                    *
! *     Urbana, IL                                                    *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
! * Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
! *********************************************************************
! * *******************************************************************
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
MODULE RocFracSubInterface
  
  USE ROCSTAR_RocFracInterp
  USE ROCSTAR_RocFrac 
  
  PUBLIC :: RocFracInterfaceInitial, RocFracInterfaceUpdate,readsdv
  
CONTAINS

  SUBROUTINE RocFracInterfaceInitial( glb,obtain_attr,surfIn)

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INCLUDE 'comf90.h'
    TYPE(ROCFRAC_GLOBAL) :: glb
    INTEGER :: MyId, NumProcs, Ierror
    CHARACTER(*), INTENT(IN)      :: surfIn
    INTEGER, INTENT(IN)           :: obtain_attr
    INTEGER :: npanes
    INTEGER, POINTER, DIMENSION(:) :: paneIDs
    
    CHARACTER*9 :: timeLevel
    CHARACTER*120 :: fracHDFFname, meshFile

    INTEGER, POINTER :: zero, one, two
    integer :: i,j, ierr
    integer :: NumElTypes2D 
    character(LEN=1), POINTER, DIMENSION(:) :: names
    character(LEN=4) :: ChrElType
    integer :: endPt, startPt, chrlngth
    integer :: pane, bcflag

    CALL MPI_COMM_RANK( glb%MPI_COMM_ROCFRAC, MyId, Ierror)
    CALL MPI_COMM_SIZE( glb%MPI_COMM_ROCFRAC, NumProcs, Ierror)

! Subroutine to register interface data.
    CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,i)
    IF(myid.eq.0 .AND. glb%debug_state) THEN
       WRITE(6,'(A)') 'Rocfrac: Calling RocFrac Register...'
    ENDIF
! Get propellant density from Rocburn ! fix if more then one material
    IF(myid.eq.0 .AND. glb%Verb.gt.1) THEN
       WRITE(6,'(A,e11.4)') 'Rocfrac: Propellant density rho(1) is',&
                       glb%rho(1)
    ENDIF
!!$ CALL COM_create_window(surWin)
    CALL COM_new_window(surWin)

! fix, should not this be allocated NumMaterials instead of 1 for more then one material
    CALL COM_new_dataitem(surWin//'.rhos', 'w', COM_DOUBLE, 1, 'kg/m^3')
    CALL COM_set_array( surWin//'.rhos', 0, glb%rho,1 )

    CALL COM_new_dataitem(surWin//'.u', 'n', COM_DOUBLE, 3, 'm')
    CALL COM_new_dataitem(surWin//'.vs', 'n', COM_DOUBLE, 3, 'm/s')
    CALL COM_new_dataitem(surWin//'.uhat', 'n', COM_DOUBLE, 3, 'm')


    CALL COM_new_dataitem(surWin//'.ts_alp', 'e', COM_DOUBLE, 1, 'Pa')
    CALL COM_new_dataitem(surWin//'.bv','n',COM_INTEGER, 1, '')
    CALL COM_new_dataitem( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
    CALL COM_new_dataitem( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')

    IF (glb%HeatTransSoln.eqv..true.) THEN
       CALL COM_new_dataitem(surWin//'.qs', 'e', COM_DOUBLE, 1, 'W/m^2')
       CALL COM_new_dataitem(surWin//'.Ts', 'n', COM_DOUBLE, 1, 'K')
    ENDIF
       
!!$ #OLD   CALL COM_new_dataitem(surWin//'.ts_alp', 'e', & 
!!$ #OLD        COM_DOUBLE_PRECISION, 3, 'Pa')

    IF ( glb%ALEenabled .eqv. .true.) CALL COM_new_dataitem( surWin//'.vbar_alp', 'n', COM_DOUBLE, 3, 'm/s')

    CALL COM_get_panes( surfIn, npanes, paneIDs)



    DO j = 1, npanes

       pane = paneIDs(j)
       
       CALL COM_copy_array( surfIn//'.bcflag', pane, bcflag)

!!!!!!!!!!!!!!!!!!!!!! Burning Surface !!!!!!!!!!!!!!!!!!
       IF(bcflag.EQ.1)THEN

          CALL COM_get_size( surfIn//".nc", pane, glb%InterfaceSFNumNodes)
          CALL COM_set_size( surWin//".nc", pane, glb%InterfaceSFNumNodes)
          CALL COM_set_size( surWin//'.bcflag', pane, 1)
          CALL COM_resize_array(surWin//'.bcflag', pane)
          
          CALL COM_get_connectivities(surfIn,pane,NumElTypes2D,names)
          
          startPt = 1
          DO i = 1, NumElTypes2D
             ! Search for the next dataitem name
             endPt = startPt
             chrlngth = 0
             DO WHILE (endPt .LE. UBOUND(names,1))
                IF (names(endPt) .NE. ' ') THEN
                   chrlngth = chrlngth + 1
                   ChrElType(chrlngth:chrlngth) = names(endPt)
                   endPt = endPt + 1
                ELSE
                   EXIT
                END IF
             END DO
             
             startPt = endPt + 1
             
             IF(ChrElType(1:chrlngth).EQ.':t6')THEN
                CALL COM_get_size( surfIn//".:t6", pane, glb%InterfaceSFNumElems)
                CALL COM_set_size( surWin//".:t6", pane, glb%InterfaceSFNumElems)
                glb%iElType2D = 6
             ELSE IF(ChrElType(1:chrlngth).EQ.':t3')THEN
                CALL COM_get_size( surfIn//".:t3", pane, glb%InterfaceSFNumElems)
                CALL COM_set_size( surWin//".:t3", pane, glb%InterfaceSFNumElems)
                glb%iElType2D = 3
             ELSE IF(ChrElType(1:chrlngth).EQ.':q4')THEN
                CALL COM_get_size( surfIn//".:q4", pane, glb%InterfaceSFNumElems)
                CALL COM_set_size( surWin//".:q4", pane, glb%InterfaceSFNumElems)
                glb%iElType2D = 4
             ELSE
                WRITE(0,'(A,A)') 'Rocfrac: Error: Surface mesh type',&
                               ' element not supported'
                WRITE(0,'(A,A)') 'Read in Element Type :: ',&
                                 ChrElType(1:chrlngth)
                CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
             ENDIF
             
          END DO

          CALL COM_free_buffer(names)  
! Fluids/Solids Interface Mesh 

          ALLOCATE(glb%InterfaceSFNodalCoors(1:3,1:glb%InterfaceSFNumNodes))
          ALLOCATE(glb%MapNodeSF(1:glb%InterfaceSFNumNodes))

          ALLOCATE(glb%InterfaceSFElemConn(1:glb%iElType2D,1:glb%InterfaceSFNumElems))
          ALLOCATE(glb%MapSFElVolEl(1:glb%InterfaceSFNumElems)) 

          
! -- The array's containing fluid-solid values

          ALLOCATE(glb%InterfaceSFNodalDisps(1:3,1:glb%InterfaceSFNumNodes))
          ALLOCATE(glb%InterfaceSFTotalNodalDisps(1:3,1:glb%InterfaceSFNumNodes))         
          ALLOCATE(glb%InterfaceSFNodalVels(1:3,1:glb%InterfaceSFNumNodes))
!#OLD    ALLOCATE(glb%InterfaceSFElemTract(1:3,1:glb%InterfaceSFNumElems))
          ALLOCATE(glb%InterfaceSFElemTract(1:glb%InterfaceSFNumElems))
!    IF(glb%ALEenabled)THEN
          ALLOCATE(glb%InterfaceSFVbar(1:3,1:glb%InterfaceSFNumNodes))

          IF(glb%HeatTransSoln)THEN
             ALLOCATE(glb%InterfaceSFHeatFlux(1:glb%InterfaceSFNumElems))
             glb%InterfaceSFHeatFlux(:) = glb%DummyFlux
             
             ALLOCATE(glb%InterfaceSFNodalTemp(1:glb%InterfaceSFNumNodes))
             glb%InterfaceSFNodalTemp(:) = glb%Temperature0
          ENDIF
!    ENDIF

          IF(glb%ipstatic.eqv..true.)THEN
             ALLOCATE(glb%pstatic(1:3,1:glb%InterfaceSFNumElems))
             glb%pstatic(1:3,1:glb%InterfaceSFNumElems) = 0.d0
          ENDIF
          glb%InterfaceSFNodalDisps(:,:) = 0.d0
          glb%InterfaceSFTotalNodalDisps(:,:) = 0.d0
          glb%InterfaceSFNodalVels(:,:) = 0.d0
          glb%InterfaceSFElemTract(:) = 0.d0
!#OLD    glb%InterfaceSFElemTract(:,:) = 0.d0
          glb%InterfaceSFVbar = 0.d0

          IF ( glb%InterfaceSFNumNodes > 0) THEN ! Fluid-solid interface

! Register Coordinates of 2D Fluid-solid interface

!!$    CALL COM_init_mesh( surWin//'.nc', MyId+1, glb%InterfaceSFNodalCoors, glb%InterfaceSFNumNodes)


!       print*,'glb%InterfaceSFNumNodes', glb%InterfaceSFNumNodes

             CALL COM_set_array(     surWin//'.nc', pane, glb%InterfaceSFNodalCoors,3 )
!
! Registering 2D Element Connectivity of Fluid-solid interface
!
             IF(glb%iElType2D.EQ.3)THEN
!!$      CALL COM_init_mesh( surWin//'.t3', MyId+1, glb%InterfaceSFElemConn, glb%InterfaceSFNumElems)

          
                CALL COM_set_array( surWin//'.:t3', pane, glb%InterfaceSFElemConn,3)

             ELSE IF(glb%iElType2D.EQ.6)THEN
!!$       CALL COM_init_mesh( surWin//'.t6', MyId+1, glb%InterfaceSFElemConn, glb%InterfaceSFNumElems)
                CALL COM_set_array( surWin//'.:t6', pane, glb%InterfaceSFElemConn,6)

             ELSE IF(glb%iElType2D.EQ.4)THEN
!!$       CALL COM_init_mesh( surWin//'.q4', MyId+1, glb%InterfaceSFElemConn, glb%InterfaceSFNumElems)

                CALL COM_set_array( surWin//'.:q4', pane, glb%InterfaceSFElemConn,4)
                
             ENDIF

             CALL COM_set_array( surWin//'.bv', pane, glb%MapNodeSF,1)
             CALL COM_set_array(surWin//'.bf2c', pane, glb%MapSFElVolEl, 1)

             CALL COM_set_array( surWin//'.u',pane, glb%InterfaceSFNodalDisps,3)
             CALL COM_set_array( surWin//'.vs', pane, glb%InterfaceSFNodalVels,3)
             CALL COM_set_array( surWin//'.uhat', pane, glb%InterfaceSFTotalNodalDisps,3)
             
             CALL COM_set_array( surWin//'.ts_alp', pane, glb%InterfaceSFElemTract, 1)

             IF (glb%ALEenabled) CALL COM_set_array( surWin//'.vbar_alp', pane, glb%InterfaceSFVbar,3)

             IF (glb%HeatTransSoln)THEN
                CALL COM_set_array( surWin//'.qs', pane, glb%InterfaceSFHeatFlux, 1)
                CALL COM_set_array( surWin//'.Ts', pane, glb%InterfaceSFNodalTemp, 1)
             ENDIF

          ENDIF

!!!!!!!!!!!!!!!!!!!!!! Non-Interacting Surface !!!!!!!!!!!!!!!!!!
       ELSE IF(bcflag.EQ.2)THEN

          CALL COM_get_size( surfIn//".nc", pane, glb%InterfaceSNumNodes)
!          WRITE(*,*) 'Number of Nodes on noninteracting structures surface: ',glb%InterfaceSNumNodes,myid
          CALL COM_set_size( surWin//".nc", pane, glb%InterfaceSNumNodes)

          CALL COM_set_size( surWin//'.bcflag', pane, 1)
          CALL COM_resize_array(surWin//'.bcflag', pane)

          CALL COM_get_connectivities(surfIn,pane,NumElTypes2D,names)
          startPt = 1
          DO i = 1, NumElTypes2D
      ! Search for the next dataitem name
             endPt = startPt
             chrlngth = 0
             DO WHILE (endPt .LE. UBOUND(names,1))
                IF (names(endPt) .NE. ' ') THEN
                   chrlngth = chrlngth + 1
                   ChrElType(chrlngth:chrlngth) = names(endPt)
                   endPt = endPt + 1
                ELSE
                   EXIT
                END IF
             END DO
             
             startPt = endPt + 1
      
             IF(ChrElType(1:chrlngth).EQ.':t6')THEN
                CALL COM_get_size( surfIn//".:t6", pane, glb%InterfaceSNumElems)
                CALL COM_set_size( surWin//".:t6", pane, glb%InterfaceSNumElems)
                glb%iElType2D = 6
             ELSE IF(ChrElType(1:chrlngth).EQ.':t3')THEN
                CALL COM_get_size( surfIn//".:t3", pane, glb%InterfaceSNumElems)
                CALL COM_set_size( surWin//".:t3", pane, glb%InterfaceSNumElems)
                glb%iElType2D = 3
             ELSE IF(ChrElType(1:chrlngth).EQ.':q4')THEN
                CALL COM_get_size( surfIn//".:q4", pane, glb%InterfaceSNumElems)
                CALL COM_set_size( surWin//".:q4", pane, glb%InterfaceSNumElems)
                glb%iElType2D = 4
             ELSE
                WRITE(0,'(A,A)') 'Rocfrac: Error: Surface mesh type',&
                               ' element not supported'
                WRITE(0,'(A,A)') 'Read in Element Type :: ',&
                                 ChrElType(1:chrlngth)
                CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
             ENDIF
             
          END DO
!
! Read No Solid/Fluid mesh
!
          CALL COM_free_buffer(names)

          ALLOCATE(glb%InterfaceSNodalCoors(1:3,1:glb%InterfaceSNumNodes))
          ALLOCATE(glb%MapNodeS(1:glb%InterfaceSNumNodes))

          ALLOCATE(glb%InterfaceSVbar(1:3,1:glb%InterfaceSNumNodes))
          glb%InterfaceSVbar = 0.d0
          
          ALLOCATE(glb%InterfaceSElemConn(1:glb%iElType2D,1:glb%InterfaceSNumElems))
          IF(glb%EnforceTractionS.eqv..true.)THEN
             ALLOCATE(glb%MapSElVolEl(1:glb%InterfaceSNumElems))
          ENDIF

          IF ( glb%InterfaceSNumNodes > 0) THEN ! Non-fluid-solid interface

! Register Coordinates of 2D Non-Fluid-solid interface

             CALL COM_set_array( surWin//'.nc', pane, glb%InterfaceSNodalCoors,3)
             IF(glb%iElType2D.EQ.3)THEN
                
                CALL COM_set_array( surWin//'.:t3', pane, glb%InterfaceSElemConn,3)
                
             ELSE IF(glb%iElType2D.EQ.6)THEN
                
                CALL COM_set_array( surWin//'.:t6', pane, glb%InterfaceSElemConn,6)
                
             ELSE IF(glb%iElType2D.EQ.4)THEN
                
                CALL COM_set_array( surWin//'.:q4', pane, glb%InterfaceSElemConn,4)
                
             ENDIF
       !fix, why was this numprocs+myid+1
!       IF(glb%ALEenabled) CALL COM_set_array( surWin//'.vbar_alp', NumProcs+pane, glb%InterfaceSVbar,3)


             IF(glb%EnforceTractionS)THEN
                CALL COM_set_array(surWin//'.bf2c', pane, glb%MapSElVolEl, 1)
             ENDIF

!!$       CALL COM_set_array( surWin//'.u',pane, glb%InterfaceSNodalDisps,3)
!!$       CALL COM_set_array( surWin//'.vs', pane, glb%InterfaceSNodalVels,3)
!!$       CALL COM_set_array( surWin//'.uhat', pane, glb%InterfaceSTotalNodalDisps,3)

             CALL COM_set_array( surWin//'.bv', pane, glb%MapNodeS,1)

             IF( glb%ALEenabled) CALL COM_set_array( surWin//'.vbar_alp', pane, glb%InterfaceSVbar,3)
       
          ENDIF


!!!!!!!!!!!!!!!!!!!!!! Non-Burning Interacting Surface !!!!!!!!!!!!!!!!!!
       else if ( bcflag.eq.0 )THEN

          CALL COM_get_size( surfIn//".nc", pane, glb%InterfaceSFnbNumNodes)
          CALL COM_set_size( surWin//".nc", pane, glb%InterfaceSFnbNumNodes)

          CALL COM_set_size( surWin//'.bcflag', pane, 1)
          CALL COM_resize_array(surWin//'.bcflag', pane)

          CALL COM_get_connectivities(surfIn,pane,NumElTypes2D,names)
          
          startPt = 1
          DO i = 1, NumElTypes2D
             ! Search for the next dataitem name
             endPt = startPt
             chrlngth = 0
             DO WHILE (endPt .LE. UBOUND(names,1))
                IF (names(endPt) .NE. ' ') THEN
                   chrlngth = chrlngth + 1
                   ChrElType(chrlngth:chrlngth) = names(endPt)
                   endPt = endPt + 1
                ELSE
                   EXIT
                END IF
             END DO
             
             startPt = endPt + 1
             
             IF(ChrElType(1:chrlngth).EQ.':t6')THEN
                CALL COM_get_size( surfIn//".:t6", pane, glb%InterfaceSFnbNumElems)
                CALL COM_set_size( surWin//".:t6", pane, glb%InterfaceSFnbNumElems)
                glb%iElType2D = 6
             ELSE IF(ChrElType(1:chrlngth).EQ.':t3')THEN
                CALL COM_get_size( surfIn//".:t3", pane, glb%InterfaceSFnbNumElems)
                CALL COM_set_size( surWin//".:t3", pane, glb%InterfaceSFnbNumElems)
                glb%iElType2D = 3
             ELSE IF(ChrElType(1:chrlngth).EQ.':q4')THEN
                CALL COM_get_size( surfIn//".:q4", pane, glb%InterfaceSFnbNumElems)
                CALL COM_set_size( surWin//".:q4", pane, glb%InterfaceSFnbNumElems)
                glb%iElType2D = 4
             ELSE
                WRITE(0,'(A,A)') 'Rocfrac: Error: Surface mesh type',&
                               ' element not supported'
                WRITE(0,'(A,A)') 'Read in Element Type :: ',&
                                 ChrElType(1:chrlngth)
                CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
             ENDIF
             
          END DO

          CALL COM_free_buffer(names)  
! Fluids/Solids Interface Mesh 

          ALLOCATE(glb%InterfaceSFnbNodalCoors(1:3,1:glb%InterfaceSFnbNumNodes))
          ALLOCATE(glb%MapNodeSFnb(1:glb%InterfaceSFnbNumNodes))

          ALLOCATE(glb%InterfaceSFnbElemConn(1:glb%iElType2D,1:glb%InterfaceSFnbNumElems))
          ALLOCATE(glb%MapSFnbElVolEl(1:glb%InterfaceSFnbNumElems)) 

          
! -- The array's containing fluid-solid values

          ALLOCATE(glb%InterfaceSFnbNodalDisps(1:3,1:glb%InterfaceSFnbNumNodes))
          ALLOCATE(glb%InterfaceSFnbTotalNodalDisps(1:3,1:glb%InterfaceSFnbNumNodes))         
          ALLOCATE(glb%InterfaceSFnbNodalVels(1:3,1:glb%InterfaceSFnbNumNodes))
!#OLD    ALLOCATE(glb%InterfaceSFElemTract(1:3,1:glb%InterfaceSFNumElems))
          ALLOCATE(glb%InterfaceSFnbElemTract(1:glb%InterfaceSFnbNumElems))
!    IF(glb%ALEenabled)THEN
          ALLOCATE(glb%InterfaceSFnbVbar(1:3,1:glb%InterfaceSFnbNumNodes))
!   ALLOCATE(glb%InterfaceSVbar(1:3,1:glb%InterfaceSNumNodes)) ! needed?
!    ENDIF

          IF(glb%ipstatic)THEN
             ALLOCATE(glb%pstatic(1:3,1:glb%InterfaceSFNumElems))
             glb%pstatic(1:3,1:glb%InterfaceSFNumElems) = 0.d0
             ALLOCATE(glb%pstaticnb(1:3,1:glb%InterfaceSFnbNumElems))
             glb%pstaticnb(1:3,1:glb%InterfaceSFnbNumElems) = 0.d0
          ENDIF
          glb%InterfaceSFnbNodalDisps(:,:) = 0.d0
          glb%InterfaceSFnbTotalNodalDisps(:,:) = 0.d0
          glb%InterfaceSFnbNodalVels(:,:) = 0.d0
          glb%InterfaceSFnbElemTract(:) = 0.d0
!#OLD    glb%InterfaceSFElemTract(:,:) = 0.d0
          glb%InterfaceSFnbVbar = 0.d0

          IF ( glb%InterfaceSFnbNumNodes > 0) THEN ! Fluid-solid interface

! Register Coordinates of 2D Fluid-solid interface

!!$    CALL COM_init_mesh( surWin//'.nc', MyId+1, glb%InterfaceSFNodalCoors, glb%InterfaceSFNumNodes)


!       print*,'glb%InterfaceSFNumNodes', glb%InterfaceSFNumNodes

             CALL COM_set_array(surWin//'.nc', pane, glb%InterfaceSFnbNodalCoors,3 )
!
! Registering 2D Element Connectivity of Fluid-solid interface
!
             IF(glb%iElType2D.EQ.3)THEN
!!$      CALL COM_init_mesh( surWin//'.t3', MyId+1, glb%InterfaceSFElemConn, glb%InterfaceSFNumElems)

          
                CALL COM_set_array( surWin//'.:t3', pane, glb%InterfaceSFnbElemConn,3)

             ELSE IF(glb%iElType2D.EQ.6)THEN
!!$       CALL COM_init_mesh( surWin//'.t6', MyId+1, glb%InterfaceSFElemConn, glb%InterfaceSFNumElems)
                CALL COM_set_array( surWin//'.:t6', pane, glb%InterfaceSFnbElemConn,6)

             ELSE IF(glb%iElType2D.EQ.4)THEN
!!$       CALL COM_init_mesh( surWin//'.q4', MyId+1, glb%InterfaceSFElemConn, glb%InterfaceSFNumElems)

                CALL COM_set_array( surWin//'.:q4', pane, glb%InterfaceSFnbElemConn,4)
                
             ENDIF

             CALL COM_set_array( surWin//'.bv', pane, glb%MapNodeSFnb,1)
             CALL COM_set_array(surWin//'.bf2c', pane, glb%MapSFnbElVolEl, 1)           
             
             CALL COM_set_array( surWin//'.u',pane, glb%InterfaceSFnbNodalDisps,3)
             CALL COM_set_array( surWin//'.vs', pane, glb%InterfaceSFnbNodalVels,3)
             CALL COM_set_array( surWin//'.uhat', pane, glb%InterfaceSFnbTotalNodalDisps,3)
             
             CALL COM_set_array( surWin//'.ts_alp', pane, glb%InterfaceSFnbElemTract, 1)

             IF ( glb%ALEenabled) CALL COM_set_array( surWin//'.vbar_alp', pane, glb%InterfaceSFnbVbar,3)

           endif
                
      ELSE
          WRITE(0,'(A,i4,A,i4)') 'Rocfrac: Error: Invalid bcflag', &
                                  bcflag,' on surface pane',pane
          STOP
      ENDIF


   enddo
!!$    IF ( glb%NumNdsBC > 0) THEN ! Non-fluid-solid interface
!!$       
!!$   !    CALL COM_init_mesh( surWin//'.nc', 2*NumProcs+MyId+1, glb%NumNdsBC)
!!$
!!$       CALL COM_set_size( surWin//'.nc',, MyId+1, glb%NumNdsBC)
!!$
!!$       CALL COM_set_array( surWin//'.vs',     2*NumProcs+MyId+1, glb%VeloBndry)
!!$       CALL COM_set_array( surWin//'.u',     2*NumProcs+MyId+1, glb%AccelBndry)
!!$    ENDIF

! No longer needed
!!$    CALL COM_new_dataitem( surWin//'.bcflag', 'p', COM_INTEGER, 1, '')
!!$    CALL COM_allocate_array( surWin//'.bcflag', MyId+1)
!!$
!!$    stop
!!$
!!$    IF ( glb%InterfaceSFNumNodes > 0) THEN ! Fluid-solid interface
!!$       CALL COM_get_array( surWin//'.bcflag', MyId+1, one)
!!$       one = 1
!!$    END IF
!!$    
!!$    IF ( glb%InterfaceSNumNodes > 0) THEN ! Non-fluid-solid interface
!!$       CALL COM_get_array( surWin//'.bcflag', NumProcs+MyId+1, two)
!!$       two = 2
!!$    ENDIF
! put into volume mesh

!!$    IF ( glb%NumNdsBC > 0) THEN ! Surface Boundary Condition flags
!!$       CALL COM_get_array( surWin//'.bcflag', 2*NumProcs+MyId+1, two)
!!$       two = 2
!!$    ENDIF


    CALL COM_window_init_done( surWin)


         
    CALL COM_free_buffer(paneIDs)  

    CALL COM_call_function( obtain_attr, 2, &
         COM_get_dataitem_handle_const( surfIn//".all"), &
         COM_get_dataitem_handle( surWin//".all"))

    CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,i)
    IF(myid.eq.0 .AND. glb%debug_state) THEN
       WRITE(6,'(A)') 'Rocfrac: Calling RocFrac Register... Done'
    ENDIF
11  FORMAT(A,'_',A,A1)
!---------------------------------------------------------------------------------------------------

  END SUBROUTINE RocFracInterfaceInitial


  SUBROUTINE RocFracInterfaceUpdate( glb, CurrentTime, &
       CurrentTimeStep, MAN_update_inbuff)
    IMPLICIT NONE

    TYPE(ROCFRAC_GLOBAL), POINTER :: glb
    REAL*8, INTENT(IN)            :: CurrentTime, CurrentTimeStep
    INTEGER, INTENT(IN)           :: MAN_update_inbuff
    
! Subroutine to gather current interface position computed by the solids code into
! the database used by the interpolation procedures.

    ! Local variables:

    REAL*8  :: SolidsDisp, SolidsCoor
    
    INTEGER :: iInterfaceNode, SolidsNodeNum
    
! Update ALL the surface mesh coordinates 

! Burning/Interacting surfaces

    DO iInterfaceNode = 1, glb%InterfaceSFNumNodes
       
       SolidsNodeNum = ABS(glb%MapNodeSF(iInterfaceNode))
       
       glb%InterfaceSFNodalCoors(1,iInterfaceNode) = glb%MeshCoor(1,SolidsNodeNum )
       glb%InterfaceSFNodalCoors(2,iInterfaceNode) = glb%MeshCoor(2,SolidsNodeNum )
       glb%InterfaceSFNodalCoors(3,iInterfaceNode) = glb%MeshCoor(3,SolidsNodeNum )
       
    END DO

! Non-Burning/Interacting surfaces

    DO iInterfaceNode = 1, glb%InterfaceSFnbNumNodes
       
       SolidsNodeNum = ABS(glb%MapNodeSFnb(iInterfaceNode))
       
       glb%InterfaceSFnbNodalCoors(1,iInterfaceNode) = glb%MeshCoor(1,SolidsNodeNum )
       glb%InterfaceSFnbNodalCoors(2,iInterfaceNode) = glb%MeshCoor(2,SolidsNodeNum )
       glb%InterfaceSFnbNodalCoors(3,iInterfaceNode) = glb%MeshCoor(3,SolidsNodeNum )
       
    END DO

! Non-Interacting surfaces

    DO iInterfaceNode = 1, glb%InterfaceSNumNodes
       
       SolidsNodeNum = ABS(glb%MapNodeS(iInterfaceNode))
       IF(SolidsNodeNum .ne. 0) THEN
       glb%InterfaceSNodalCoors(1,iInterfaceNode) = glb%MeshCoor(1,SolidsNodeNum )
       glb%InterfaceSNodalCoors(2,iInterfaceNode) = glb%MeshCoor(2,SolidsNodeNum )
       glb%InterfaceSNodalCoors(3,iInterfaceNode) = glb%MeshCoor(3,SolidsNodeNum )
       ENDIF
       
    END DO
!---------------------------------------------------------------------------------------------------

  END SUBROUTINE RocFracInterfaceUpdate


  SUBROUTINE readsdv(glb,myid)

    IMPLICIT NONE

    INCLUDE 'comf90.h'
    INCLUDE 'mpif.h'

    INTEGER :: myid

    TYPE(ROCFRAC_GLOBAL),POINTER  :: glb

    INTEGER :: pid


    INTEGER :: hdl_read, hdl_obtain, hdl_all

    CHARACTER(*), PARAMETER :: OverlayWin = "Overlay"
    CHARACTER(*), PARAMETER :: OverlayWin2 = "Overlay2"


    INTEGER :: comm_self
    CHARACTER(*), PARAMETER :: prefix1 = "A"
    CHARACTER(*), PARAMETER :: prefix2 = "B"

    CHARACTER(LEN=5) :: sdv_material
    CHARACTER(LEN=12) :: sdv_wname
    CHARACTER(LEN=27) :: fname1, fname2 
    CHARACTER(LEN=3) ::  ichr3
    INTEGER :: i,j

    
    INTEGER, POINTER, DIMENSION(:) :: MapFaceEl2Vol1a, FaceOfVolEL1a

    INTEGER :: iprocs,ios, iaux
! obtain function handle ------------------------------------------------------


!    IF(myid.eq.1) THEN

    pid = (myid+1)*100 + 3


    WRITE(ichr3,'(I3.3)') pid ! problem if over 999 processors, fix

    fname1 = 'Rocfrac/Rocin/A_'//ichr3//'_sdv.hdf'
    fname2 = 'Rocfrac/Rocin/B_'//ichr3//'_sdv.hdf'


    CALL SimIN_load_module( "SDV_IN")

    hdl_read   = COM_get_function_handle( 'SDV_IN.read_windows')
    hdl_obtain = COM_get_function_handle( 'SDV_IN.obtain_dataitem')




! Define the base-window and sdv-window names

    sdv_material = prefix1//'_sdv'
    sdv_wname = OverlayWin//sdv_material
    CALL COM_new_window( OverlayWin )

!  // Read the pane from the given file. Note that the file contains both
!  // the parent and the subdivided windows. Read only the subdivided one.
    comm_self = MPI_COMM_SELF
    CALL COM_call_function( hdl_read, 4, fname1, OverlayWin, &
         sdv_material, comm_self)
    hdl_all = COM_get_dataitem_handle( sdv_wname//".all")
    CALL COM_call_function( hdl_obtain, 3, hdl_all, hdl_all, pid)
!  // Obtain number of sub-nodes, sub-faces, nodes, and faces
  
    CALL COM_get_size( sdv_wname//'.sn_parent_fcID', pid, glb%nsubn1)
    
!    PRINT*,'Number of sub-nodes =',glb%nsubn1

    CALL COM_get_size( sdv_wname//'.:t3', pid, glb%nsubf1) 

!    PRINT*,'Number of sub-faces =', glb%nsubf1

    ALLOCATE(glb%Sthresh1(1:3, 1:glb%nsubf1)) 
    glb%Sthresh1(1:3, 1:glb%nsubf1) = glb%Sinit(1) ! fix, this should be read-in


!    // Obtain the connectivity

    ALLOCATE( glb%sd_subfaces1(1:3,1:glb%nsubf1) ) 

    CALL COM_copy_array( sdv_wname//".:t3", pid, glb%sd_subfaces1, 3)

!!$    DO i = 1, glb%nsubf1
!!$       PRINT*,glb%sd_subfaces1(1:3,i)
!!$    ENDDO

    ALLOCATE(glb%sd_coor1(1:3,1:glb%nsubn1))

    CALL COM_copy_array( sdv_wname//".nc", pid, glb%sd_coor1, 3)

!!$    DO i = 1, glb%nsubn1
!!$       PRINT*,glb%sd_coor1(1:3,i)
!!$  ENDDO
  
    ALLOCATE(glb%sd_subface_parents1(1:glb%nsubf1))

    CALL COM_copy_array( sdv_wname//".sf_parent", pid, glb%sd_subface_parents1)
    
!!$  DO i = 1, glb%nsubf1
!!$     PRINT*,'sd_subface_parents',glb%sd_subface_parents1(i)
!!$  ENDDO


    ALLOCATE(glb%sd_subface_nat_coors1(1:6,1:glb%nsubf1))



    CALL COM_copy_array( sdv_wname//".sf_ncs", pid, &
         glb%sd_subface_nat_coors1, 6)

!
!  // NOTE: The last argument (6) indicates that the local coordinates are 
!  // stored consecutively (xi1, eta1, xi2, eta2, xi3, eta3). Use the number 
!  // one (1) if the xi1 for all nodes are stored together and then xi2, etc..
 
!!$  DO i = 1, glb%nsubf1
!!$     PRINT*, i,glb%nsubf1,glb%sd_subface_nat_coors1(1:6,i)
!!$  END DO
    ALLOCATE(glb%sd_subface_counterparts1(1:glb%nsubf1) )


    CALL COM_copy_array( sdv_wname//".sf_cntpt_fcID", pid, glb%sd_subface_counterparts1)



!!$  DO i = 1, glb%nsubf1
!!$     PRINT*,glb%sd_subface_counterparts1(i)
!!$  ENDDO



!!$  PRINT*,'finished readsdv',myid
!!$  CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,i)
!!$  stop

!  // Delete the window created by Rocin.

!!$    CALL COM_delete_window( sdv_wname)
!!$
!!$    
!!$    CALL Rocin_unload_module( "SDV_IN")
!!$
!!$
!!$    CALL Rocin_load_module( "SDV_IN")

! Define the base-window and sdv-window names

    sdv_material = prefix2//'_sdv'
    sdv_wname = OverlayWin//sdv_material

!    CALL COM_new_window( OverlayWin )





!  // Read the pane from the given file. Note that the file contains both
!  // the parent and the subdivided windows. Read only the subdivided one.
    comm_self = MPI_COMM_SELF
    CALL COM_call_function( hdl_read, 4, fname2, OverlayWin, &
         sdv_material, comm_self)
    hdl_all = COM_get_dataitem_handle( sdv_wname//".all")
    CALL COM_call_function( hdl_obtain, 3, hdl_all, hdl_all, pid)
    CALL COM_get_size( sdv_wname//'.sn_parent_fcID', pid, glb%nsubn2)
    



!  PRINT*,'Number of sub-nodes =',glb%nsubn2
  
    CALL COM_get_size( sdv_wname//'.:t3', pid, glb%nsubf2)


  
!  PRINT*,'Number of sub-faces =', glb%nsubf2
  
    ALLOCATE(glb%sd_subface_counterparts2(1:glb%nsubf2) )
     
    CALL COM_copy_array( sdv_wname//".sf_cntpt_fcID", pid, glb%sd_subface_counterparts2)


!    // Obtain the connectivity

    ALLOCATE( glb%sd_subfaces2(1:3,1:glb%nsubf1) ) ! GETS STUCK HERE FOR more then ONE processor

!!$    PRINT*,'finished readsdv!!!!!!!!!',myid,glb%nsubn2,glb%nsubf2
!!$    CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,i)

    CALL COM_copy_array( sdv_wname//".:t3", pid, glb%sd_subfaces2, 3)

!!$  DO i = 1, glb%nsubf2
!!$     PRINT*,glb%sd_subfaces2(1:3,i)
!!$  ENDDO

    ALLOCATE(glb%sd_coor2(1:3,1:glb%nsubn2))
    CALL COM_copy_array( sdv_wname//".nc", pid, glb%sd_coor2, 3) 


!!$  DO i = 1, glb%nsubn2
!!$     PRINT*,glb%sd_coor2(1:3,i)
!!$  ENDDO



!  PRINT*,'ldfdfkjsdflkj'
    ALLOCATE(glb%sd_subface_nat_coors2(1:6,1:glb%nsubf2))
!  PRINT*,'dsflkjsdflkj'
    CALL COM_copy_array( sdv_wname//".sf_ncs", pid, &
         glb%sd_subface_nat_coors2, 6)
!
!  PRINT*,'lkjsdflkj'
 
!
!  // NOTE: The last argument (6) indicates that the local coordinates are 
!  // stored consecutively (xi1, eta1, xi2, eta2, xi3, eta3). Use the number 
!  // one (1) if the xi1 for all nodes are stored together and then xi2, etc..
 
!!$  DO i = 1, glb%nsubf2
!!$     PRINT*, i,glb%nsubf2,glb%sd_subface_nat_coors2(1:6,i)
!!$  END DO  



!!$  DO i = 1, glb%nsubf2
!!$     PRINT*,i,glb%nsubf2,glb%sd_subface_counterparts2(i)
!!$  ENDDO

!  PRINT*,'lkjsdflkj',glb%nsubf2,myid
  ALLOCATE(glb%sd_subface_parents2(1:glb%nsubf2))
!  PRINT*,'copy',myid
  CALL COM_copy_array( sdv_wname//".sf_parent", pid, glb%sd_subface_parents2)
!  PRINT*,'sdfsdflkjsdflkj',myid
!!$  DO i = 1, glb%nsubf2
!!$     PRINT*,'sd_subface_parents',glb%sd_subface_parents2(i)
!!$  ENDDO


!  // Delete the window created by Rocin.
  CALL COM_delete_window( sdv_wname)
!  PRINT*,'sdfsdflkjsdflkj',myid

! // Unload Rocin from Roccom.

  CALL SimIN_unload_module( "SDV_IN")



!  CALL COM_new_dataitem( surWin//'.bf2c', 'e', COM_INTEGER, 1, '')
!  CALL COM_allocate_array(surWin//'.bf2c', iNI, ElFlag_List, 1)
!  CALL COM_new_dataitem( surWin//'.faceOnCell', 'e', COM_INTEGER, 1, '')
!  CALL COM_allocate_array(surWin//'.faceOnCell', iNI, FaceOnCell, 1)
  

!  ENDIF

!  stop



END SUBROUTINE  readsdv
  
END MODULE RocFracSubInterface

