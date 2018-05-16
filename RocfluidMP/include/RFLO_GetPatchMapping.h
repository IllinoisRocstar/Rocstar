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
INTERFACE
  SUBROUTINE RFLO_GetPatchMapping( lb,lbs,l1SrcDir,l2SrcDir,align, &
                                   idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                                   ibeg,iend,jbeg,jend,kbeg,kend, &
                                   ibegSrc,iendSrc,jbegSrc,jendSrc, &
                                   kbegSrc,kendSrc,mapMat )
    INTEGER :: lb, lbs, l1SrcDir, l2SrcDir
    INTEGER :: idir, jdir, kdir, idirSrc, jdirSrc, kdirSrc
    INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
    INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc
    INTEGER :: mapMat(3,4)
    LOGICAL :: align
  END SUBROUTINE RFLO_GetPatchMapping
END INTERFACE






