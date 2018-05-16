
!!****
!!
!!  NAME
!!     HeatLoad_Hex
!!
!!  FUNCTION
!!     This subroutine takes a heat flux (InterfaceHeatFlux)(W/m2) and divides it amoung the 
!!     nodes associated with the interface patch.
!!
!!  INPUTS
!!     coor -- Coordinates of the nodes
!!     ElConnVol -- The connectivity table
!!     dmat -- Material stiffness matrix
!!     numnp -- Number of nodes
!!     NumEl -- Number of elements
!!     MatType -- Mapping of materials to elements
!!     NumMatType -- Total number of materials
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE HeatLoad_Hex(Q_ex, numnp, InterfaceHeatFlux, InterfaceNumElems, InterfaceNumNodes, &
     InterfaceElemConn, MapNode, LwrBnd, UppBnd, coor)

  IMPLICIT NONE

  INTEGER :: LwrBnd,UppBnd

  INTEGER :: InterfaceNumElems,InterfaceNumNodes
  INTEGER, DIMENSION(1:UppBnd,1:InterfaceNumElems) :: InterfaceElemConn
  INTEGER, DIMENSION(1:InterfaceNumNodes) :: MapNode
  REAL*8,  DIMENSION(1:InterfaceNumElems) :: InterfaceHeatFlux

  INTEGER :: numnp          ! number of nodes

  REAL*8, DIMENSION(1:numnp) :: Q_ex ! external force

  INTEGER :: nx, ny, nz
  INTEGER :: i,j

  INTEGER,DIMENSION(1:4) :: Conn2D
  real*8 :: qmag
  REAL*8 :: ElemArea

  ! ... mesh coordinates
  REAL*8, DIMENSION(1:3,1:numnp) :: coor

  REAL*8, DIMENSION(1:3) :: veca, vecb, outwd_surface_normal, &
       half_normal_veca, half_normal_vecb

  DO i = 1, InterfaceNumElems

     ! ... Nodes of the Quad ( Vertex Nodes ).  Gives surface node numbers,
     ! ... not global node numbers
     Conn2D(1:4) = InterfaceElemConn(1:4,i)

     ! ... Uses the fact that the norm of the cross product vector
     ! ... is the area of the parallelogram they form.  The triangle they
     ! ... form has half that area.

     ! ... Vector between element interface surface nodes 2 and 1
     veca(1:3) =  coor(1:3,MapNode(Conn2D(2)))- coor(1:3,MapNode(Conn2D(1)))
     ! ... Vector between element interface surface nodes 2 and 1
     vecb(1:3) =  coor(1:3,MapNode(Conn2D(4)))- coor(1:3,MapNode(Conn2D(1)))

     ! ... The norm of the vector 'half_normal_veca' is the area of the parallelogram
     ! ... formed by element interface surface nodes 2-1-4.
     half_normal_veca(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
     half_normal_veca(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
     half_normal_veca(3) = veca(1)*vecb(2) - vecb(1)*veca(2)

     ! ... Vector between element interface surface nodes 4 and 3        
     veca(1:3) = coor(1:3,MapNode(Conn2D(4))) - coor(1:3,MapNode(Conn2D(3)))
     ! ... Vector between element interface surface nodes 2 and 3
     vecb(1:3) = coor(1:3,MapNode(Conn2D(2))) - coor(1:3,MapNode(Conn2D(3)))

     ! ... The norm of the vector 'half_normal_vecb' is the area of the parallelogram
     ! ... formed by element interface surface nodes 2-3-4.        
     half_normal_vecb(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
     half_normal_vecb(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
     half_normal_vecb(3) = veca(1)*vecb(2) - vecb(1)*veca(2)

     ! ... Outward surface normal is a vector normal to the element
     ! ... interface surface with magnitude of
     ! ... 2*(element interface surface area)
     outwd_surface_normal = half_normal_veca + half_normal_vecb 

     ! ... Element interface area
     ElemArea = 0.0d0
     do j = 1,3
        ElemArea = ElemArea + outwd_surface_normal(j) * outwd_surface_normal(j)
     end do
     ElemArea = SQRT(ElemArea) * 0.50d0

     ! ... qmag is the heat flux/(2*(element interface surface area)) magnitude into
     ! ... each of the four nodes on the element interface surface.
     ! ... Note: multiplied by factor of two contained in outward_surface_normal below.
     qmag = InterfaceHeatFlux(i)*0.250d0
     qmag = qmag * ElemArea
     
     ! ... Assembly into global external load vector Q_ex
     ! ... InterfaceElemConn contains interface surface
     ! ... node numbers corresponding to element i.
     ! ... (Not global element i, surface element i)
     Conn2D(1:4) = InterfaceElemConn(LwrBnd:UppBnd,i)

     ! ... MapNode gives the global node number of a given surface
     ! ... node number j.
     DO j = 1, 4
        Q_ex(MapNode(Conn2D(j))) = Q_ex(MapNode(Conn2D(j))) + qmag
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE HeatLoad_Hex
