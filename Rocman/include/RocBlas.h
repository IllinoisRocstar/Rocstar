/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/

#ifndef _ROCBLAS_H_
#define _ROCBLAS_H_

class RocBlas
{
public:
  RocBlas() {}
  static void init();
  static void initHandles();
public:
  static int copy;
  static int add;
  static int sub;
  static int limit1;
  static int mul;
  static int div;
  static int neg;
  static int axpy;
  static int nrm2;
  static int copy_scalar;
  static int sub_scalar;
  static int axpy_scalar;
  static int div_scalar;
  static int mul_scalar;
  static int max_scalar_MPI;
  static int min_scalar_MPI;
  static int sum_scalar_MPI;
  static int nrm2_scalar_MPI;
  static int maxof_scalar;
};

#endif






