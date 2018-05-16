!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
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
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
     MODULE data_declarations

       IMPLICIT NONE

       INTEGER :: i,j,k,l,m,n,       &  ! Loop Counters
                  etype,             &  ! Element type (4 or 10 node)
                  ftype,             &  ! Face type (3 or 6 node)
                  numnp,             &  ! Number of nodes
                  numbndnp,          &  ! Number of boundary marker nodes
                  numel,             &  ! Number of elements
                  npfixed,           &  ! Number of original surface nodes (in .points file)
                  numsub,            &  ! Number of subregions assigned to elements as flags
                  numfaces,          &  ! Number of faces describing the surface mesh
                  numbndfaces,       &  ! Number of faces with markers .NE. '0'
                  numscale_bk,       &  ! Number of sister nodes on back face
                  numscale_ft,       &  ! Number of sister nodes on front face
                  numscale_np,       &  ! Total Number of sister node for both faces
                  mark,              &  ! Temp marker holder
                  scale,             &  ! The number of times for duplication
                  count,             &  ! Temporary counter
                  frontfile,         &  ! previous file ID number
                  backfile,          &  ! next file ID number
                  priv1, priv2          ! Temp holders for passing private variables between files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       INTEGER :: numprop,           &  ! ?
                  numbound,          &  ! Number of Boundary flags
                  numboundmesh,      &  ! Number of boundary mesh nodes
                  NdBCflag,          &  ! Nodal Boundary condition flag
                  numco,             &
                  iaux,              &
                  numclst,           &
                  numcstet,          &
                  numcohshared,      &
                  nproc_neigh,       &
                  num_border_coh,    &
                  num_border_cst,    &
                  numboundmid,       &
                  numboundend,       &
                  numboundfirst,     &
                  numbmeshmid,       &
                  numbmeshend,       &
                  numbmeshfirst  

      REAL*8  :: delta             ! ?
      REAL*8  :: maxZ, minZ        ! Variables used to hold max and min values for node coords

      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: &
                 coor,          &  ! Node coordinates array
                 meshcoor,      &  ! Wish I knew
                 xm                ! This one, too

      INTEGER, ALLOCATABLE, DIMENSION(:) :: &
                 nboundtype,    &  ! Nodal boundary marker array
                 matcstet,      &
                 ts_proportion     ! ?

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: &
                 lmcstet          ! Element connectivity


      REAL*8, ALLOCATABLE, DIMENSION(:) :: &
                 proportion

      INTEGER :: NBoundryEl3D

      TYPE id_struct
         INTEGER :: ID
         INTEGER :: NdBCflag
         INTEGER :: sister
      END TYPE id_struct

      TYPE(id_struct), ALLOCATABLE, DIMENSION(:) ::  &
                 id,            &
                 idmesh

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
       REAL*8  :: XYZ(3)                ! Temp coordinates holder
       REAL*8  :: Ztrans                ! The displacement by which the mesh nodes are moved

       CHARACTER*40 :: PREFIX           ! The generic filename prefix
       CHARACTER*4  :: scaleID          ! The ID number of the filename extension

       

       TYPE node_struct
          INTEGER :: marker             ! Marker (flag) assigned to this node
          INTEGER :: meshmark           ! Mesh motion boundary flag assigned to this node
          REAL*8  :: coord(3)           ! x,y,z coordinates for this node
       END TYPE node_struct

       TYPE element_struct
          INTEGER :: marker
          INTEGER, POINTER, DIMENSION(:) :: &
                     conn               ! Element connectivity array
       END TYPE element_struct

       TYPE face_struct
          INTEGER :: ID                 ! The face ID number
          INTEGER :: marker             ! Marker (flag) assigned to this face
          INTEGER, POINTER, DIMENSION(:) :: &
                     conn               ! Connectivity array for this face
          TYPE(face_struct), POINTER :: next
       END TYPE face_struct

       TYPE File_struct
          INTEGER :: comm(2)            ! The current files being communicated with 1=back 2=front
          INTEGER :: num_sister_np      ! The number of shared sister nodes between files
          INTEGER :: numboundnp         ! The number of boundary nodes accounted for
                                        !    This will change because of internal nodes on original
                                        !    exterior surfaces that are now interior surfaces
       END TYPE File_struct

       TYPE scale_struct
          INTEGER :: ID
          INTEGER :: marker
          INTEGER :: mesh
          REAL*8  :: coord(3)
          TYPE(scale_struct), POINTER :: next
          TYPE(scale_struct), POINTER :: sister
       END TYPE scale_struct

       TYPE elem_list
          INTEGER :: ID
          INTEGER :: conn(4)
          INTEGER :: flag
          TYPE(elem_list), POINTER :: next
          TYPE(elem_list), POINTER :: previous
       END TYPE elem_list

       TYPE surface_struct
          INTEGER :: marker
          INTEGER :: conn(3)
       END TYPE surface_struct

       TYPE surfscale_struct
          INTEGER :: ID
          INTEGER :: marker
          REAL*8  :: coord(3)
          TYPE(surfscale_struct), POINTER :: next
          TYPE(surfscale_struct), POINTER :: sister
       END TYPE surfscale_struct

       TYPE(surface_struct), ALLOCATABLE, DIMENSION(:) ::  &
                  surface
       TYPE(surfscale_struct), POINTER :: frontsnode, first_frontsnode, backsnode, first_backsnode


     ! User defined type declarations
     TYPE(element_struct),   ALLOCATABLE, DIMENSION(:) :: element
     TYPE(node_struct),      ALLOCATABLE, DIMENSION(:) :: node
     TYPE(node_struct),      ALLOCATABLE, DIMENSION(:) :: surfnode
     TYPE(face_struct), POINTER :: facelist, firstface
     TYPE(File_struct), ALLOCATABLE, DIMENSION(:) :: FileID
     TYPE(scale_struct), POINTER ::        &
                      frontnode,           &  ! Front face node list
                      first_front_node,    &
                      backnode,            &  ! Back face node list
                      first_back_node,     &  
                      BNPlist,             &  ! All other boundary nodes list
                      firstnode
     
     ! These linked lists will be used to keep track of what elements have points that 
     ! lie on the boundary
     TYPE(elem_list), POINTER :: elemlist, firstelem, &
                                 eboundlist, firstebound

     CONTAINS

!-----------------------------------------------------------------------------------------!
!            FUNCTION LENGTH: obtains the length of a passed in character string          !
!-----------------------------------------------------------------------------------------!
       INTEGER FUNCTION LENGTH(STRING)

       ! Determines the length of string
       ! Preconditions: STRING is defined
       ! Postconditions: The function result is the actual string length
       !                 and includes only the characters preceding blanks

       ! Declarations
         CHARACTER*(*) STRING    ! Input Dummy character type

       ! Local
         CHARACTER *1 BLANK      ! Used for conditional statement
         PARAMETER (BLANK = ' ')

         ! Length holder
         INTEGER NEXT

       ! Starts with the first character and finds the first nonblank
         DO 10 NEXT = LEN(STRING), 1, -1
            IF (STRING(NEXT:NEXT) .NE. BLANK) THEN
               LENGTH = NEXT
               RETURN
             END IF
     10  CONTINUE

       ! If all characters are blanks
         LENGTH = 0

       ! Exit function
         RETURN
       END FUNCTION LENGTH

     END MODULE data_declarations


