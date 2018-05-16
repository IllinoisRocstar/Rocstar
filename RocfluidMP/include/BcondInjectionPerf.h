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
  SUBROUTINE BcondInjectionPerf( distrib,minj,tinj,rhoVrel,sxn,syn,szn, &
                                 cpgas,mm,p,rhob,rhoub,rhovb,rhowb,rhoeb,pb, &
                                 uinj,vinj,winj )
    USE ModDataTypes
    INTEGER,      INTENT(IN)  :: distrib
    REAL(RFREAL), INTENT(IN)  :: cpgas, minj, mm, sxn, syn, szn, tinj, &
                                 rhoVrel(3), p
    REAL(RFREAL), INTENT(OUT) :: rhob, rhoub, rhovb, rhowb, rhoeb, pb, &
                                 uinj, vinj, winj
  END SUBROUTINE BcondInjectionPerf
END INTERFACE






