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
      Program select_bbx
      
        implicit none

        INTEGER, PARAMETER ::  double = SELECTED_REAL_KIND(P=14,R=30)
        real(kind=double) :: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,eps=1.0e-06
        real(kind=double) :: G_xmin,G_xmax,G_ymin,G_ymax,G_zmin,G_zmax
        real(kind=double) :: xmn,xmx,ymn,ymx,zmn,zmx
        integer           :: idom,ios
        character(len=6)  :: string

!
! Select blocks that lie in a specified bounding box
!
! Reads output from hdf2bbx_d
!
      open(7,file='bbox',form='formatted',status='old',iostat=ios)
      if (ios /= 0) then
        write(*,*) 'I need file bbox from program hdf2bbx_d or reader'
        stop
      else
        open(10,file='BBOX',form='formatted',status='old',iostat=ios)
        read(10,*) string, G_xmin,G_xmax,G_ymin,G_ymax,G_zmin,G_zmax
        write(*,*) 'Global BBX is ',  &
                          G_xmin,G_xmax,G_ymin,G_ymax,G_zmin,G_zmax
      endif
!
! Read image bounding box
!
      write(*,*) 'Enter Xmin, Xmax, Ymin, Ymax, Zmin, Zmax'
      read (*,*) Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
!
! Read bounding box for each block, and write list of visible blocks 
!
      open(8,file='Bbox',form='formatted',status='unknown')
      open(9,file='Ibox',form='formatted',status='unknown')
!
      do
        read(7,*,iostat=ios) idom,xmn,xmx,ymn,ymx,zmn,zmx
        if (ios /= 0) exit
!
! Eliminate blocks outside the image bounding box
!
        if (xmn.gt.Xmax-eps .or. ymn.gt.Ymax-eps .or. zmn.gt.Zmax-eps) cycle
        if (xmx.lt.Xmin+eps .or. ymx.lt.Ymin+eps .or. zmx.lt.Zmin+eps) cycle
!
! This block is visible
!
        write(8,10) idom
        write(9,*) idom
 10     format('block_',i4.4)
!
      enddo
!
      end program select_bbx






