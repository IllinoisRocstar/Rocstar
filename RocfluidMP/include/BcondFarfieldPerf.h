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
  SUBROUTINE BcondFarfieldPerf( machInf,alphaInf,betaInf,pInf,tInf, & 
                                sxn,syn,szn,cpgas,mol, & 
                                rho,rhou,rhov,rhow,rhoe,press, &
                                rhob,rhoub,rhovb,rhowb,rhoeb,pb )
    USE ModDataTypes
    REAL(RFREAL) :: machInf, alphaInf, betaInf, pInf, tInf
    REAL(RFREAL) :: rho, rhou, rhov, rhow, rhoe, press
    REAL(RFREAL) :: sxn, syn, szn, cpgas, mol
    REAL(RFREAL) :: rhob, rhoub, rhovb, rhowb, rhoeb, pb
  END SUBROUTINE BcondFarfieldPerf
END INTERFACE






