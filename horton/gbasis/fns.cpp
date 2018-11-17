// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2017 The HORTON Development Team
//
// This file is part of HORTON.
//
// HORTON is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// HORTON is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

// #define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include "horton/moments.h"
#include "horton/gbasis/boys.h"
#include "horton/gbasis/cartpure.h"
#include "horton/gbasis/fns.h"
#include "horton/gbasis/dir_gaussderive.h"


/*
    GB1GridFn
*/

GB1GridFn::GB1GridFn(long max_shell_type, long dim_work, long dim_output)
    : GBCalculator(max_shell_type), dim_work(dim_work), dim_output(dim_output),
      shell_type0(0), r0(NULL), point(NULL), i1p() {
  nwork = max_nbasis*dim_work;
  work_cart = new double[nwork];
  work_pure = new double[nwork];
}

void GB1GridFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type)) {
    throw std::domain_error("shell_type0 out of range.");
  }
  shell_type0 = _shell_type0;
  r0 = _r0;
  point = _point;
  // We make use of the fact that a floating point zero consists of
  // consecutive zero bytes.
  memset(work_cart, 0, nwork*sizeof(double));
  memset(work_pure, 0, nwork*sizeof(double));
}

void GB1GridFn::cart_to_pure() {
  /*
     The initial results are always stored in work_cart. The projection
     routine always outputs its result in work_pure. Once that is done,
     the pointers to both blocks are swapped such that the final result is
     always back in work_cart.
  */

  if (shell_type0 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type0,
      1,          // anterior
      dim_work);  // posterior
    swap_work();
  }
}


/*
    GB1ExpGridOrbitalFn
*/

void GB1ExpGridOrbitalFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  if (shell_type0 != 0) {
    poly_work[1] = point[0] - r0[0];
    poly_work[2] = point[1] - r0[1];
    poly_work[3] = point[2] - r0[2];
    offset = fill_cartesian_polynomials(poly_work+1, abs(shell_type0))+1;
  } else {
    offset = 0;
  }
}

void GB1ExpGridOrbitalFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The (Cartesian) basis function evaluated in `point` (see reset method) is added to
  // work_cart, where ibasis0 is an index for the primitive in the current Cartesian
  // shell.
  for (long ibasis0=get_shell_nbasis(abs(shell_type0))-1; ibasis0 >= 0; ibasis0--) {
    work_cart[ibasis0] += pre*scales0[ibasis0]*poly_work[ibasis0+offset];
  }
}

void GB1ExpGridOrbitalFn::compute_point_from_exp(double* work_basis, double* coeffs,
                                                 long nbasis, double* output) {
  for (long i=0; i < norb; i++) {
    long iorb = iorbs[i];
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
      // Just evaluate the contribution of each basis function to an orbital.
      // The values of the basis functions at `point` (see reset method) are available in
      // work_basis.
      output[i] += coeffs[ibasis*nfn + iorb]*work_basis[ibasis];
    }
  }
}


/*
 GB1ExpGridOrbGradientFn
 */

void GB1ExpGridOrbGradientFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher and lower polynomials are required because of first derivative.
  offset_h1 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+1)+1;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
}

void GB1ExpGridOrbGradientFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // For every primitive, work_cart contains four values computed at `point` (see reset
  // method): the basis function and its derivatives toward x, y and z.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function derived toward x
    work_cart[3*i1p.ibasis0+0] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[3*i1p.ibasis0+0] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
      work_cart[3*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[3*i1p.ibasis0+1] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[3*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[3*i1p.ibasis0+2] += i1p.n0[2]*pre0*poly_work[i1p.ibasis0-nnotx-1+offset_l1];
    } while (i1p.inc());
}

void GB1ExpGridOrbGradientFn::compute_point_from_exp(double* work_basis, double* coeffs,
                                                     long nbasis, double* output) {
  for (long i=0; i < norb; i++) {
    double g_x = 0, g_y = 0, g_z = 0;
    long iorb = iorbs[i];
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
      g_x += coeffs[ibasis*nfn + iorb]*work_basis[ibasis*3+0];
      g_y += coeffs[ibasis*nfn + iorb]*work_basis[ibasis*3+1];
      g_z += coeffs[ibasis*nfn + iorb]*work_basis[ibasis*3+2];
    }
    output[i*3+0] += g_x;
    output[i*3+1] += g_y;
    output[i*3+2] += g_z;
  }
}


/*
    GB1DMGridDensityFn
*/

void GB1DMGridDensityFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  /*XHC -_shell_type0 (angular momentum of the shell)
        -_r0 (coordinates of the center)
        -_point (coordinates of the grid point)
  */
  GB1GridFn::reset(_shell_type0, _r0, _point);
  //XHC
  poly_work[0] = 1.0;
  //XHC Only for non-s orbitals because the factor becomes 1 (i-A_i)
  if (shell_type0 != 0) {
    poly_work[1] = point[0] - r0[0];
    poly_work[2] = point[1] - r0[1];
    poly_work[3] = point[2] - r0[2];
    offset = fill_cartesian_polynomials(poly_work+1, abs(shell_type0))+1;
  } else {
    offset = 0;
  }
}

void GB1DMGridDensityFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The (Cartesian) basis function evaluated in `point` (see reset method) is added to
  // work_cart, where ibasis0 is an index for the primitive in the current Cartesian
  // shell.
  //printf("add\n");
  for (long ibasis0=get_shell_nbasis(abs(shell_type0))-1; ibasis0 >= 0; ibasis0--) {
    work_cart[ibasis0] += pre*scales0[ibasis0]*poly_work[ibasis0+offset];
    //printf("poly_work[%ld+%ld] = %f\n", ibasis0, offset, poly_work[ibasis0+offset]);
  }
}

void GB1DMGridDensityFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The values of the basis functions at `point` (see reset method) are available in
  // work_basis.

  // The epsilon argument allows one to skip part of the calculation where the density is
  // low.
  if (epsilon > 0) {
    double absmax_basis = 0.0;
    // compute the maximum basis function
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
      double tmp = fabs(work_basis[ibasis]);
      if (tmp > absmax_basis) absmax_basis = tmp;
    }
    // upper estimate of the density
    double rho_upper = 0.0;
    for (long ibasis=0; ibasis < nbasis; ibasis++) {
      rho_upper += fabs(work_basis[ibasis])*dmmaxrow[ibasis];
    }
    rho_upper *= nbasis*absmax_basis;

    // if the upper bound is too low, do not compute density.
    if (rho_upper < epsilon) return;

    // modify epsilon to avoid recomputation
    epsilon /= absmax_basis*nbasis*nbasis;
  }

  // Loop over all basis functions and add significant contributions
  double rho = 0.0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    // if the contribution of this loop is smaller than epsilon/nbasis, skip it.
    if (epsilon > 0) {
      if (fabs(work_basis[ibasis0])*dmmaxrow[ibasis0] < epsilon)
          continue;
    }
    double tmp = 0;
    // Loop for off-diagonal contributions of density matrix.
    for (long ibasis1=ibasis0-1; ibasis1 >= 0; ibasis1--) {
      tmp += work_basis[ibasis1]*dm[ibasis0*nbasis+ibasis1];
    }
    // Finally, also include diagonal contribution
    rho += (2*tmp +  // off-diagonal
            dm[ibasis0*(nbasis+1)]*work_basis[ibasis0])  // diagonal
           *work_basis[ibasis0];
  }
  *output += rho;
}

void GB1DMGridDensityFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The potential in `point` (see reset method) is given in `*pot`. It is the functional
  // derivative of the energy w.r.t. to the density in `point`. This gets transformed to
  // a contribution to the Fock matrix, i.e. the derivative of the energy w.r.t. the 1RDM
  // elements, using the chain rule.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp1 = (*pot)*work_basis[ibasis0];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double tmp2 = tmp1*work_basis[ibasis1];
      fock[ibasis0*nbasis+ibasis1] += tmp2;
      if (ibasis0 != ibasis1) fock[ibasis1*nbasis+ibasis0] += tmp2;
    }
  }
}



/*
    GB1DMGridGradientFn
*/

void GB1DMGridGradientFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  /*printf("reset\n");
  printf("poly_work[0] = %f\n", poly_work[0]);
  printf("poly_work[1] = %f\n", poly_work[1]);
  printf("poly_work[2] = %f\n", poly_work[2]);
  printf("poly_work[3] = %f\n", poly_work[3]);*/
  // One order higher and lower polynomials are required because of first derivative.
  offset_h1 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+1)+1;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
}

void GB1DMGridGradientFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // For every primitive, work_cart contains four values computed at `point` (see reset
  // method): the basis function and its derivatives toward x, y and z.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    //printf("exp = %f \n", exp(-alpha0*dist_sq(r0, point)));
    //printf("scale[%d] = %f \n", i1p.ibasis0, scales0[i1p.ibasis0]);
    //printf("pre0 = %f \n", pre0);
    double pre0_h = -pre0*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[4*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    //printf("before work_cart[%d] = %f \n", 4*i1p.ibasis0+1, work_cart[4*i1p.ibasis0+1]);
    //XHCwork_cart[4*i1p.ibasis0+1] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 0);
    //printf("poly_work = %f \n", poly_work[i1p.ibasis0+offset_h1]);
    work_cart[4*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[4*i1p.ibasis0+1] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];
    //printf("work_cart[%d] = %f \n", 4*i1p.ibasis0+1, work_cart[4*i1p.ibasis0+1]);


    // Basis function derived toward y
    //printf("before work_cart[%d] = %f \n", 4*i1p.ibasis0+2, work_cart[4*i1p.ibasis0+2]);
    //XHCwork_cart[4*i1p.ibasis0+2] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 0);
    work_cart[4*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[4*i1p.ibasis0+2] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];
    //printf("work_cart[%d] = %f \n", 4*i1p.ibasis0+2, work_cart[4*i1p.ibasis0+2]);

    // Basis function derived toward z
    //printf("before work_cart[%d] = %f \n", 4*i1p.ibasis0+3, work_cart[4*i1p.ibasis0+3]);
    //XHCwork_cart[4*i1p.ibasis0+3] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 1);
    work_cart[4*i1p.ibasis0+3] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[4*i1p.ibasis0+3] += i1p.n0[2]*pre0*poly_work[i1p.ibasis0-nnotx-1+offset_l1];
    //printf("work_cart[%d] = %f \n", 4*i1p.ibasis0+3, work_cart[4*i1p.ibasis0+3]);
  } while (i1p.inc());
}

void GB1DMGridGradientFn::compute_point_from_dm(double* work_basis, double* dm,
                                                long nbasis, double* output,
                                                double epsilon, double* dmmaxrow) {
  // The value of the basis function and its derivatives toward x, y and z (at `point`)
  // are available in work_basis. These are used, together with the 1RDM coefficients,
  // to compute the density gradient at `point`.
  double rho_x = 0, rho_y = 0, rho_z = 0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      row += work_basis[ibasis1*4]*dm[ibasis0*nbasis+ibasis1];
    }
    rho_x += row*work_basis[ibasis0*4+1];
    rho_y += row*work_basis[ibasis0*4+2];
    rho_z += row*work_basis[ibasis0*4+3];
  }
  output[0] += 2*rho_x;
  output[1] += 2*rho_y;
  output[2] += 2*rho_z;
}

void GB1DMGridGradientFn::compute_fock_from_pot(double* pot, double* work_basis,
                                                long nbasis, double* fock) {
  // The functional derivative of the energy w.r.t. to the density gradient is given in
  // the argument *pot (three components: x, y, and z). The chain rule is used to
  // transform these into a contribution to Fock matrix.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp0 = work_basis[ibasis0*4];
    double tmp1 = pot[0]*work_basis[ibasis0*4+1] +
                  pot[1]*work_basis[ibasis0*4+2] +
                  pot[2]*work_basis[ibasis0*4+3];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double result = tmp0*(pot[0]*work_basis[ibasis1*4+1] +
                            pot[1]*work_basis[ibasis1*4+2] +
                            pot[2]*work_basis[ibasis1*4+3]) +
                      tmp1*work_basis[ibasis1*4];
      fock[ibasis1*nbasis+ibasis0] += result;
      if (ibasis1 != ibasis0) {
        // Enforce symmetry
        fock[ibasis0*nbasis+ibasis1] += result;
      }
    }
  }
}


/*
    GB1DMGridGGAFn
*/

void GB1DMGridGGAFn::compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                           double* output, double epsilon,
                                           double* dmmaxrow) {
  // The value of the basis function and its derivatives toward x, y and z (at `point`)
  // are available in work_basis. These are used, together with the 1RDM coefficients,
  // to compute the density and its gradient at `point`.
  double rho = 0, rho_x = 0, rho_y = 0, rho_z = 0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      row += work_basis[ibasis1*4]*dm[ibasis0*nbasis+ibasis1];
    }
    rho += row*work_basis[ibasis0*4];
    rho_x += row*work_basis[ibasis0*4+1];
    rho_y += row*work_basis[ibasis0*4+2];
    rho_z += row*work_basis[ibasis0*4+3];
  }
  output[0] += rho;
  output[1] += 2*rho_x;
  output[2] += 2*rho_y;
  output[3] += 2*rho_z;
}

void GB1DMGridGGAFn::compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                           double* fock) {
  // The functional derivative of the energy w.r.t. to the density and its gradient are
  // given inthe argument *pot (four components). The chain rule is used to transform
  // these into a contribution to Fock matrix.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp0 = pot[0]*work_basis[ibasis0*4] +
                  pot[1]*work_basis[ibasis0*4+1] +
                  pot[2]*work_basis[ibasis0*4+2] +
                  pot[3]*work_basis[ibasis0*4+3];
    double tmp1 = pot[1]*work_basis[ibasis0*4];
    double tmp2 = pot[2]*work_basis[ibasis0*4];
    double tmp3 = pot[3]*work_basis[ibasis0*4];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double result = tmp0*work_basis[ibasis1*4] +
                      tmp1*work_basis[ibasis1*4+1] +
                      tmp2*work_basis[ibasis1*4+2] +
                      tmp3*work_basis[ibasis1*4+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      if (ibasis1 != ibasis0) {
        // Enforce symmetry
        fock[ibasis0*nbasis+ibasis1] += result;
      }
    }
  }
}


/*
    GB1DMGridKineticFn
*/

void GB1DMGridKineticFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h1 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+1)+1;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
}

void GB1DMGridKineticFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // For every primitive, work_cart contains three values computed at `point` (see reset
  // method): the derivatives of the basis function toward x, y and z.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function derived toward x
    work_cart[3*i1p.ibasis0] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[3*i1p.ibasis0] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
    work_cart[3*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[3*i1p.ibasis0+1] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[3*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[3*i1p.ibasis0+2] += i1p.n0[2]*pre0*poly_work[i1p.ibasis0-nnotx-1+offset_l1];
  } while (i1p.inc());
}

void GB1DMGridKineticFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The derivatives of the basis function w.r.t. x, y and z, at `point` (see reset
  // method) are given in work_basis. Together with the 1RDM coefficients, dm, these are
  // used to compute the kinetic energy density in `point`.
  double tau = 0.0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;
    for (long ibasis1=0; ibasis1 < ibasis0; ibasis1++) {
      tmp_x += work_basis[ibasis1*3  ]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*3+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*3+2]*dm[ibasis0*nbasis+ibasis1];
    }
    tmp_x += 0.5*work_basis[ibasis0*3  ]*dm[ibasis0*nbasis+ibasis0];
    tmp_y += 0.5*work_basis[ibasis0*3+1]*dm[ibasis0*nbasis+ibasis0];
    tmp_z += 0.5*work_basis[ibasis0*3+2]*dm[ibasis0*nbasis+ibasis0];
    tau += tmp_x*work_basis[ibasis0*3  ] +
           tmp_y*work_basis[ibasis0*3+1] +
           tmp_z*work_basis[ibasis0*3+2];
  }
  *output += tau;
}

void GB1DMGridKineticFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the kinetic energy density in `point`, is given
  // in *pot. This is used to compute a contribution to the Fock matrix from `point`.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp_x = 0.5*(*pot)*work_basis[ibasis0*3  ];
    double tmp_y = 0.5*(*pot)*work_basis[ibasis0*3+1];
    double tmp_z = 0.5*(*pot)*work_basis[ibasis0*3+2];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double result = tmp_x*work_basis[ibasis1*3  ] +
                      tmp_y*work_basis[ibasis1*3+1] +
                      tmp_z*work_basis[ibasis1*3+2];
      fock[ibasis1*nbasis+ibasis0] += result;
      if (ibasis0 != ibasis1)
        fock[ibasis0*nbasis+ibasis1] += result;
    }
  }
}


/*
    GB1DMGridHessianFn
*/

void GB1DMGridHessianFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h2 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+2)+1;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
//printf("1XHC l: %d, offset_h1: %d, offset: %d, offset_l1: %d\n",abs(shell_type0),offset_h1,offset,offset_l1);
//printf("2XHC offset_h2: %d, offset_l2: %d\n",offset_h2,offset_l2);
}

void GB1DMGridHessianFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The value of the basis functions in the Cartesian primitive shell are added to
  // work_basis, together with its first and second derivatives, all evaluated at `point`.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    double pre0_hh = -pre0_h*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[10*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[10*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[10*i1p.ibasis0+1] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
    work_cart[10*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[10*i1p.ibasis0+2] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[10*i1p.ibasis0+3] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[10*i1p.ibasis0+3] += i1p.n0[2]*pre0*poly_work[i1p.ibasis0-nnotx-1+offset_l1];

    // Basis function derived toward xx
    work_cart[10*i1p.ibasis0+4] += pre0_hh*poly_work[i1p.ibasis0+offset_h2];
    work_cart[10*i1p.ibasis0+4] += (2*i1p.n0[0]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[0] > 1)
      work_cart[10*i1p.ibasis0+4] += i1p.n0[0]*(i1p.n0[0]-1)*pre0*
                                      poly_work[i1p.ibasis0+offset_l2];

    // Basis function derived toward xy
    work_cart[10*i1p.ibasis0+5] += pre0_hh*poly_work[i1p.ibasis0+1+nnotx+offset_h2];
    if (i1p.n0[0] > 0)
      work_cart[10*i1p.ibasis0+5] += i1p.n0[0]*pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset];
    if (i1p.n0[1] > 0)
      work_cart[10*i1p.ibasis0+5] += i1p.n0[1]*pre0_h*poly_work[i1p.ibasis0-nnotx+offset];
    if ((i1p.n0[0] > 0) && (i1p.n0[1] > 0))
      work_cart[10*i1p.ibasis0+5] += i1p.n0[0]*i1p.n0[1]*pre0*
                                     poly_work[i1p.ibasis0-nnotx+offset_l2];

    // Basis function derived toward xz
    work_cart[10*i1p.ibasis0+6] += pre0_hh*poly_work[i1p.ibasis0+2+nnotx+offset_h2];
    if (i1p.n0[0] > 0)
      work_cart[10*i1p.ibasis0+6] += i1p.n0[0]*pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset];
    if (i1p.n0[2] > 0)
      work_cart[10*i1p.ibasis0+6] += i1p.n0[2]*pre0_h*poly_work[i1p.ibasis0-1-nnotx+offset];
    if ((i1p.n0[0] > 0) && (i1p.n0[2] > 0))
      work_cart[10*i1p.ibasis0+6] += i1p.n0[0]*i1p.n0[2]*pre0*
                                     poly_work[i1p.ibasis0-1-nnotx+offset_l2];

    // Basis function derived toward yy
    work_cart[10*i1p.ibasis0+7] += pre0_hh*poly_work[i1p.ibasis0+3+2*nnotx+offset_h2];
    work_cart[10*i1p.ibasis0+7] += (2*i1p.n0[1]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[1] > 1)
      work_cart[10*i1p.ibasis0+7] += i1p.n0[1]*(i1p.n0[1]-1)*pre0*
                                     poly_work[i1p.ibasis0+1-2*nnotx+offset_l2];

    // Basis function derived toward yz
    work_cart[10*i1p.ibasis0+8] += pre0_hh*poly_work[i1p.ibasis0+4+2*nnotx+offset_h2];
    if (i1p.n0[1] > 0)
      work_cart[10*i1p.ibasis0+8] += i1p.n0[1]*pre0_h*poly_work[i1p.ibasis0+1+offset];
    if (i1p.n0[2] > 0)
      work_cart[10*i1p.ibasis0+8] += i1p.n0[2]*pre0_h*poly_work[i1p.ibasis0-1+offset];
    if ((i1p.n0[1] > 0) && (i1p.n0[2] > 0))
      work_cart[10*i1p.ibasis0+8] += i1p.n0[1]*i1p.n0[2]*pre0*
                                     poly_work[i1p.ibasis0-2*nnotx+offset_l2];

    // Basis function derived toward zz
    work_cart[10*i1p.ibasis0+9] += pre0_hh*poly_work[i1p.ibasis0+5+2*nnotx+offset_h2];
    work_cart[10*i1p.ibasis0+9] += (2*i1p.n0[2]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[2] > 1)
      work_cart[10*i1p.ibasis0+9] += i1p.n0[2]*(i1p.n0[2]-1)*pre0*
                                     poly_work[i1p.ibasis0-1-2*nnotx+offset_l2];
  } while (i1p.inc());
}

void GB1DMGridHessianFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The density Hessian is computed in `point` using the results in work_basis (basis
  // function values and their first and second derivatives).
  double rho_xx = 0, rho_xy = 0, rho_xz = 0;
  double rho_yy = 0, rho_yz = 0, rho_zz = 0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      row += work_basis[ibasis1*10]*dm[ibasis0*nbasis+ibasis1];
      tmp_x += work_basis[ibasis1*10+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*10+2]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*10+3]*dm[ibasis0*nbasis+ibasis1];
    }
    rho_xx += row*work_basis[ibasis0*10+4] + tmp_x*work_basis[ibasis0*10+1];
    rho_xy += row*work_basis[ibasis0*10+5] + tmp_x*work_basis[ibasis0*10+2];
    rho_xz += row*work_basis[ibasis0*10+6] + tmp_x*work_basis[ibasis0*10+3];
    rho_yy += row*work_basis[ibasis0*10+7] + tmp_y*work_basis[ibasis0*10+2];
    rho_yz += row*work_basis[ibasis0*10+8] + tmp_y*work_basis[ibasis0*10+3];
    rho_zz += row*work_basis[ibasis0*10+9] + tmp_z*work_basis[ibasis0*10+3];
  }
  output[0] += 2*rho_xx;
  output[1] += 2*rho_xy;
  output[2] += 2*rho_xz;
  output[3] += 2*rho_yy;
  output[4] += 2*rho_yz;
  output[5] += 2*rho_zz;
}


void GB1DMGridHessianFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the density Hessian matrix elements, evaluated in
  // `point` is given in `*pot` (six elements). These are transformed with the chain rule
  // to a contribution to the Fock matrix from `point`.
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp = pot[0]*work_basis[ibasis0*10+4]
                +pot[1]*work_basis[ibasis0*10+5]
                +pot[2]*work_basis[ibasis0*10+6]
                +pot[3]*work_basis[ibasis0*10+7]
                +pot[4]*work_basis[ibasis0*10+8]
                +pot[5]*work_basis[ibasis0*10+9];
    double tmp_x = pot[0]*work_basis[ibasis0*10+1]
                  +pot[1]*work_basis[ibasis0*10+2]
                  +pot[2]*work_basis[ibasis0*10+3];
    double tmp_y = pot[1]*work_basis[ibasis0*10+1]
                  +pot[3]*work_basis[ibasis0*10+2]
                  +pot[4]*work_basis[ibasis0*10+3];
    double tmp_z = pot[2]*work_basis[ibasis0*10+1]
                  +pot[4]*work_basis[ibasis0*10+2]
                  +pot[5]*work_basis[ibasis0*10+3];
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      double result = tmp*work_basis[ibasis1*10]
                     +tmp_x*work_basis[ibasis1*10+1]
                     +tmp_y*work_basis[ibasis1*10+2]
                     +tmp_z*work_basis[ibasis1*10+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      fock[ibasis0*nbasis+ibasis1] += result;
    }
  }
}


/*
    GB1DMGridLaplacianFn
*/

void GB1DMGridLaplacianFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h2 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+2)+1;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
}

void GB1DMGridLaplacianFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The value of the basis functions in the Cartesian primitive shell are added to
  // work_basis, together with its first and second derivatives, all evaluated at `point`.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    //double pre0_hh = -pre0_h*2.0*alpha0;
    //long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[10*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[10*i1p.ibasis0+1] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 0);

    // Basis function derived toward y
    work_cart[10*i1p.ibasis0+2] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 0);

    // Basis function derived toward z
    work_cart[10*i1p.ibasis0+3] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 1);

    // Basis function derived toward xx
    work_cart[10*i1p.ibasis0+4] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 0);

    // Basis function derived toward xy
    work_cart[10*i1p.ibasis0+5] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 0);

    // Basis function derived toward xz
    work_cart[10*i1p.ibasis0+6] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 1);

    // Basis function derived toward yy
    work_cart[10*i1p.ibasis0+7] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 0);

    // Basis function derived toward yz
    work_cart[10*i1p.ibasis0+8] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 1);

    // Basis function derived toward zz
    work_cart[10*i1p.ibasis0+9] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 2);
  } while (i1p.inc());
}

void GB1DMGridLaplacianFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The density Laplacian is computed in `point` using the results in work_basis (basis
  // function values and their first and second derivatives).
  double rho_xx = 0;//rho_xy = 0, rho_xz = 0;
  double rho_yy = 0;//, rho_yz = 0, rho_zz = 0;
  double rho_zz = 0;
  //XHC Atomic Orbital b
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;
    //XHC Atomic Orbital a
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      //XHC row only has one value (scalar) and tmp has 3 (gradient is vector)
      row += work_basis[ibasis1*10]*dm[ibasis0*nbasis+ibasis1];
      tmp_x += work_basis[ibasis1*10+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*10+2]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*10+3]*dm[ibasis0*nbasis+ibasis1];
    }
    rho_xx += row*work_basis[ibasis0*10+4] + tmp_x*work_basis[ibasis0*10+1];
    //rho_xy += row*work_basis[ibasis0*10+5] + tmp_x*work_basis[ibasis0*10+2];
    //rho_xz += row*work_basis[ibasis0*10+6] + tmp_x*work_basis[ibasis0*10+3];
    rho_yy += row*work_basis[ibasis0*10+7] + tmp_y*work_basis[ibasis0*10+2];
    //rho_yz += row*work_basis[ibasis0*10+8] + tmp_y*work_basis[ibasis0*10+3];
    rho_zz += row*work_basis[ibasis0*10+9] + tmp_z*work_basis[ibasis0*10+3];
  }
  *output += 2*rho_xx;
  *output += 2*rho_yy;
  *output += 2*rho_zz;
  //output[1] += 2*rho_xy;
  //output[2] += 2*rho_xz;
  //output[3] += 2*rho_yy;
  //output[4] += 2*rho_yz;
  //output[5] += 2*rho_zz;
  //printf("Hola\n");
}


void GB1DMGridLaplacianFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the density Hessian matrix elements, evaluated in
  // `point` is given in `*pot` (six elements). These are transformed with the chain rule
  // to a contribution to the Fock matrix from `point`.
  /*for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp = pot[0]*work_basis[ibasis0*10+4]
                +pot[1]*work_basis[ibasis0*10+5]
                +pot[2]*work_basis[ibasis0*10+6]
                +pot[3]*work_basis[ibasis0*10+7]
                +pot[4]*work_basis[ibasis0*10+8]
                +pot[5]*work_basis[ibasis0*10+9];
    double tmp_x = pot[0]*work_basis[ibasis0*10+1]
                  +pot[1]*work_basis[ibasis0*10+2]
                  +pot[2]*work_basis[ibasis0*10+3];
    double tmp_y = pot[1]*work_basis[ibasis0*10+1]
                  +pot[3]*work_basis[ibasis0*10+2]
                  +pot[4]*work_basis[ibasis0*10+3];
    double tmp_z = pot[2]*work_basis[ibasis0*10+1]
                  +pot[4]*work_basis[ibasis0*10+2]
                  +pot[5]*work_basis[ibasis0*10+3];
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      double result = tmp*work_basis[ibasis1*10]
                     +tmp_x*work_basis[ibasis1*10+1]
                     +tmp_y*work_basis[ibasis1*10+2]
                     +tmp_z*work_basis[ibasis1*10+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      fock[ibasis0*nbasis+ibasis1] += result;
    }
  }*/
}


/*
    GB1DMGridGradofSqGradFn
*/

void GB1DMGridGradofSqGradFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h2 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+2)+1;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
}

void GB1DMGridGradofSqGradFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The value of the basis functions in the Cartesian primitive shell are added to
  // work_basis, together with its first and second derivatives, all evaluated at `point`.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    //double pre0_hh = -pre0_h*2.0*alpha0;
    //long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[10*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[10*i1p.ibasis0+1] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 0);

    // Basis function derived toward y
    work_cart[10*i1p.ibasis0+2] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 0);

    // Basis function derived toward z
    work_cart[10*i1p.ibasis0+3] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 1);

    // Basis function derived toward xx
    work_cart[10*i1p.ibasis0+4] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 0);

    // Basis function derived toward xy
    work_cart[10*i1p.ibasis0+5] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 0);

    // Basis function derived toward xz
    work_cart[10*i1p.ibasis0+6] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 1);

    // Basis function derived toward yy
    work_cart[10*i1p.ibasis0+7] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 0);

    // Basis function derived toward yz
    work_cart[10*i1p.ibasis0+8] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 1);

    // Basis function derived toward zz
    work_cart[10*i1p.ibasis0+9] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 2);
  } while (i1p.inc());
}

void GB1DMGridGradofSqGradFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The density Laplacian is computed in `point` using the results in work_basis (basis
  // function values and their first and second derivatives).
  double factor_x = 0, x_term = 0;
  double factor_y = 0, y_term = 0;
  double factor_z = 0, z_term = 0;
  double tmp_ab = 0, tmp_cd = 0;
  //XHC Atomic Orbital d
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    //double tmp = 0;
    //XHC Atomic Orbital c
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      tmp_cd += work_basis[ibasis1*10]*dm[ibasis0*nbasis+ibasis1];
    }
    factor_x += tmp_cd*work_basis[ibasis0*10+1];
    factor_y += tmp_cd*work_basis[ibasis0*10+2];
    factor_z += tmp_cd*work_basis[ibasis0*10+3];
  }
  //XHC Atomic Orbital b
  for (long ibasis2=0; ibasis2 < nbasis; ibasis2++) {
    //double tmp = 0, tmp_x = 0, tmp_y = 0, tmp_z = 0;
    double tmp_x = 0, tmp_y = 0, tmp_z = 0;
    //XHC Atomic Orbital a
    for (long ibasis3=0; ibasis3 < nbasis; ibasis3++) {
      tmp_ab += work_basis[ibasis3*10]*dm[ibasis2*nbasis+ibasis3];
      tmp_x += work_basis[ibasis3*10+1]*dm[ibasis2*nbasis+ibasis3];
      tmp_y += work_basis[ibasis3*10+2]*dm[ibasis2*nbasis+ibasis3];
      tmp_z += work_basis[ibasis3*10+3]*dm[ibasis2*nbasis+ibasis3];
    }
    x_term += factor_x*(tmp_ab*work_basis[ibasis2*10+4] +tmp_x*work_basis[ibasis2*10+1]);
    x_term += factor_y*(tmp_ab*work_basis[ibasis2*10+5] +tmp_x*work_basis[ibasis2*10+2]);
    x_term += factor_z*(tmp_ab*work_basis[ibasis2*10+6] +tmp_x*work_basis[ibasis2*10+3]);
    y_term += factor_x*(tmp_ab*work_basis[ibasis2*10+5] +tmp_y*work_basis[ibasis2*10+1]);
    y_term += factor_y*(tmp_ab*work_basis[ibasis2*10+7] +tmp_y*work_basis[ibasis2*10+2]);
    y_term += factor_z*(tmp_ab*work_basis[ibasis2*10+8] +tmp_y*work_basis[ibasis2*10+3]);
    z_term += factor_x*(tmp_ab*work_basis[ibasis2*10+6] +tmp_z*work_basis[ibasis2*10+1]);
    z_term += factor_y*(tmp_ab*work_basis[ibasis2*10+8] +tmp_z*work_basis[ibasis2*10+2]);
    z_term += factor_z*(tmp_ab*work_basis[ibasis2*10+9] +tmp_z*work_basis[ibasis2*10+3]);
  }
  output[0] += 8*x_term;
  output[1] += 8*y_term;
  output[2] += 8*z_term;
}


void GB1DMGridGradofSqGradFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the density Hessian matrix elements, evaluated in
  // `point` is given in `*pot` (six elements). These are transformed with the chain rule
  // to a contribution to the Fock matrix from `point`.
  /*for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp = pot[0]*work_basis[ibasis0*10+4]
                +pot[1]*work_basis[ibasis0*10+5]
                +pot[2]*work_basis[ibasis0*10+6]
                +pot[3]*work_basis[ibasis0*10+7]
                +pot[4]*work_basis[ibasis0*10+8]
                +pot[5]*work_basis[ibasis0*10+9];
    double tmp_x = pot[0]*work_basis[ibasis0*10+1]
                  +pot[1]*work_basis[ibasis0*10+2]
                  +pot[2]*work_basis[ibasis0*10+3];
    double tmp_y = pot[1]*work_basis[ibasis0*10+1]
                  +pot[3]*work_basis[ibasis0*10+2]
                  +pot[4]*work_basis[ibasis0*10+3];
    double tmp_z = pot[2]*work_basis[ibasis0*10+1]
                  +pot[4]*work_basis[ibasis0*10+2]
                  +pot[5]*work_basis[ibasis0*10+3];
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      double result = tmp*work_basis[ibasis1*10]
                     +tmp_x*work_basis[ibasis1*10+1]
                     +tmp_y*work_basis[ibasis1*10+2]
                     +tmp_z*work_basis[ibasis1*10+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      fock[ibasis0*nbasis+ibasis1] += result;
    }
  }*/
}


/*
    GB1DMGridLapofSqGradFn
*/

void GB1DMGridLapofSqGradFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h3 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+3)+1;
  offset_h2 = offset_h3 - ((abs(shell_type0)+3)*(abs(shell_type0)+4))/2;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
  offset_l3 = offset_l2 - ((abs(shell_type0)-2)*(abs(shell_type0-1)))/2;
}

void GB1DMGridLapofSqGradFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The value of the basis functions in the Cartesian primitive shell are added to
  // work_basis, together with its first and second derivatives, all evaluated at `point`.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    //double pre0_hh = -pre0_h*2.0*alpha0;
    //long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[20*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[20*i1p.ibasis0+1] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 0);

    // Basis function derived toward y
    work_cart[20*i1p.ibasis0+2] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 0);

    // Basis function derived toward z
    work_cart[20*i1p.ibasis0+3] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 1);

    // Basis function derived toward xx
    work_cart[20*i1p.ibasis0+4] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 0);

    // Basis function derived toward xy
    work_cart[20*i1p.ibasis0+5] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 0);

    // Basis function derived toward xz
    work_cart[20*i1p.ibasis0+6] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 1);

    // Basis function derived toward yy
    work_cart[20*i1p.ibasis0+7] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 0);

    // Basis function derived toward yz
    work_cart[20*i1p.ibasis0+8] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 1);

    // Basis function derived toward zz
    work_cart[20*i1p.ibasis0+9] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 2);

    // Basis function derived toward xxx
    work_cart[20*i1p.ibasis0+10] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 3, 0, 0);

    // Basis function derived toward xxy
    work_cart[20*i1p.ibasis0+11] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 1, 0);

    // Basis function derived toward xxz
    work_cart[20*i1p.ibasis0+12] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 1);

    // Basis function derived toward xyy
    work_cart[20*i1p.ibasis0+13] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 2, 0);

    // Basis function derived toward xyz
    work_cart[20*i1p.ibasis0+14] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 1);

    // Basis function derived toward xzz
    work_cart[20*i1p.ibasis0+15] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 2);

    // Basis function derived toward yyy
    work_cart[20*i1p.ibasis0+16] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 3, 0);

    // Basis function derived toward yyz
    work_cart[20*i1p.ibasis0+17] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 1);

    // Basis function derived toward yzz
    work_cart[20*i1p.ibasis0+18] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 2);

    // Basis function derived toward zzz
    work_cart[20*i1p.ibasis0+19] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 3);
  } while (i1p.inc());
}

void GB1DMGridLapofSqGradFn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The density Nabla3 is computed in `point` using the results in work_basis (basis
  // function values and their first and second derivatives).
  double rho_x = 0, rho_y = 0, rho_z = 0;
  //double rho_xx, rho_yy, rho_zz;
  double tmp_x = 0, tmp_y = 0, tmp_z = 0;
  //double nabla3_x_ao, nabla3_y_ao, nabla3_z_ao;
  double tmp = 0, tmp_0 = 0, tmp_x2 = 0;
  double tmp_y2 = 0, tmp_z2 = 0;
  double factor_0x, factor_xx, factor_xy, factor_xz;
  double factor_0y, factor_yx, factor_yy, factor_yz;
  double factor_0z, factor_zx, factor_zy, factor_zz;
  double factor_0_x2, factor_0_x3, factor_x_x, factor_x_x2;
  double factor_x2_x, factor_0_xy, factor_0_x2y, factor_x_y;
  double factor_x_xy, factor_x2_y, factor_0_xz, factor_0_x2z;
  double factor_x_z, factor_x_xz, factor_x2_z, factor_0_xy2;
  double factor_y_x, factor_y_xy, factor_y2_x, factor_z0;
  double lap_ao;
  //XHC Atomic Orbital d
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp = 0;
    //XHC Atomic Orbital c
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      //XHC row only has one value (scalar) and tmp has 3 (gradient is vector)
      tmp += work_basis[ibasis1*20]*dm[ibasis0*nbasis+ibasis1];
      tmp_x += work_basis[ibasis1*20+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*20+2]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*20+3]*dm[ibasis0*nbasis+ibasis1];
    }
    factor_0x += tmp*work_basis[ibasis0*20+1];
    factor_xx += tmp_x*work_basis[ibasis0*20+1] + tmp*work_basis[ibasis0*20+4];
    factor_xy += tmp_x*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+5];
    factor_xz += tmp_x*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+6];
    factor_0y += tmp*work_basis[ibasis0*20+2];
    factor_yx += tmp_y*work_basis[ibasis0*20+1] + tmp*work_basis[ibasis0*20+5];
    factor_yy += tmp_y*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+7];
    factor_yz += tmp_y*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+8];
    factor_0z += tmp*work_basis[ibasis0*20+3];
    factor_zx += tmp_z*work_basis[ibasis0*20+1] + tmp*work_basis[ibasis0*20+6];
    factor_zy += tmp_z*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+8];
    factor_zz += tmp_z*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+9];
  }
  //XHC Atomic Orbital b
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    //double row = 0;
    //XHC Atomic Orbital a
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      //XHC row only has one value (scalar) and tmp has 3 (gradient is vector)
      //XHC P_ab*AO_a
      tmp_0 += work_basis[ibasis1*20]*dm[ibasis0*nbasis+ibasis1];
      tmp_x += work_basis[ibasis1*20+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_x2 += work_basis[ibasis1*20+4]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*20+2]*dm[ibasis0*nbasis+ibasis1];
      tmp_y2 += work_basis[ibasis1*20+7]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*20+3]*dm[ibasis0*nbasis+ibasis1];
      tmp_z2 += work_basis[ibasis1*20+9]*dm[ibasis0*nbasis+ibasis1];
    } 
    //XHC First term P_ab*AO_a*nabla3_AO_b
    factor_0_x2 += tmp*work_basis[ibasis0*20+1];
    factor_0_x3 += tmp*work_basis[ibasis0*20+1];
    factor_x_x += tmp_x*work_basis[ibasis0*20+1] + tmp*work_basis[ibasis0*20+4];
    factor_x_x2 += 2.0*tmp_x*work_basis[ibasis0*20+1] + tmp*work_basis[ibasis0*20+4];
    factor_x2_x += tmp_x*work_basis[ibasis0*20+1] + tmp*work_basis[ibasis0*20+4];

    factor_0_xy += tmp_x*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+5];
    factor_0_x2y += tmp_x*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+5];
    factor_x_y += tmp_x*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+5];
    factor_x_xy += tmp_x*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+5];
    factor_x2_y += tmp_x*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+5];

    factor_0_xz += tmp_x*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+6];
    factor_0_x2z += tmp_x*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+6];
    factor_x_z += tmp_x*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+6];
    factor_x_xz += tmp_x*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+6];
    factor_x2_z += tmp_x*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+6];

    factor_0_xy += tmp*work_basis[ibasis0*20+2];
    factor_0_xy2 += tmp*work_basis[ibasis0*20+2];
    factor_y_x += tmp*work_basis[ibasis0*20+2];
    factor_y_xy += tmp*work_basis[ibasis0*20+2];
    factor_y2_x += tmp*work_basis[ibasis0*20+2];

    factor_yx += tmp_y*work_basis[ibasis0*20+1] + tmp*work_basis[ibasis0*20+5];
    factor_yy += tmp_y*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+7];
    factor_yz += tmp_y*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+8];
    factor_z0 += tmp*work_basis[ibasis0*20+3];
    factor_zx += tmp_z*work_basis[ibasis0*20+1] + tmp*work_basis[ibasis0*20+6];
    factor_zy += tmp_z*work_basis[ibasis0*20+2] + tmp*work_basis[ibasis0*20+8];
    factor_zz += tmp_z*work_basis[ibasis0*20+3] + tmp*work_basis[ibasis0*20+9];
    rho_x += tmp_x*lap_ao;
    rho_y += tmp_y*lap_ao;
    rho_z += tmp_z*lap_ao;
 
  }
  output[0] += 2*rho_x;
  output[1] += 2*rho_y;
  output[2] += 2*rho_z;
  //output[3] += 2*rho_yy;
  //output[4] += 2*rho_yz;
  //output[5] += 2*rho_zz;
  //printf("Hola\n");
}


void GB1DMGridLapofSqGradFn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the density Hessian matrix elements, evaluated in
  // `point` is given in `*pot` (six elements). These are transformed with the chain rule
  // to a contribution to the Fock matrix from `point`.
  /*for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp = pot[0]*work_basis[ibasis0*10+4]
                +pot[1]*work_basis[ibasis0*10+5]
                +pot[2]*work_basis[ibasis0*10+6]
                +pot[3]*work_basis[ibasis0*10+7]
                +pot[4]*work_basis[ibasis0*10+8]
                +pot[5]*work_basis[ibasis0*10+9];
    double tmp_x = pot[0]*work_basis[ibasis0*10+1]
                  +pot[1]*work_basis[ibasis0*10+2]
                  +pot[2]*work_basis[ibasis0*10+3];
    double tmp_y = pot[1]*work_basis[ibasis0*10+1]
                  +pot[3]*work_basis[ibasis0*10+2]
                  +pot[4]*work_basis[ibasis0*10+3];
    double tmp_z = pot[2]*work_basis[ibasis0*10+1]
                  +pot[4]*work_basis[ibasis0*10+2]
                  +pot[5]*work_basis[ibasis0*10+3];
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      double result = tmp*work_basis[ibasis1*10]
                     +tmp_x*work_basis[ibasis1*10+1]
                     +tmp_y*work_basis[ibasis1*10+2]
                     +tmp_z*work_basis[ibasis1*10+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      fock[ibasis0*nbasis+ibasis1] += result;
    }
  }*/
}


/*
    GB1DMGridNabla3Fn
*/

void GB1DMGridNabla3Fn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h3 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+3)+1;
  offset_h2 = offset_h3 - ((abs(shell_type0)+3)*(abs(shell_type0)+4))/2;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
  offset_l3 = offset_l2 - ((abs(shell_type0)-2)*(abs(shell_type0-1)))/2;
}

void GB1DMGridNabla3Fn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The value of the basis functions in the Cartesian primitive shell are added to
  // work_basis, together with its first and second derivatives, all evaluated at `point`.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    //double pre0_h = -pre0*2.0*alpha0;
    //double pre0_hh = -pre0_h*2.0*alpha0;
    //long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[20*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[20*i1p.ibasis0+1] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 0);

    // Basis function derived toward y
    work_cart[20*i1p.ibasis0+2] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 0);

    // Basis function derived toward z
    work_cart[20*i1p.ibasis0+3] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 1);

    // Basis function derived toward xx
    work_cart[20*i1p.ibasis0+4] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 0);

    // Basis function derived toward xy
    work_cart[20*i1p.ibasis0+5] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 0);

    // Basis function derived toward xz
    work_cart[20*i1p.ibasis0+6] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 1);

    // Basis function derived toward yy
    work_cart[20*i1p.ibasis0+7] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 0);

    // Basis function derived toward yz
    work_cart[20*i1p.ibasis0+8] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 1);

    // Basis function derived toward zz
    work_cart[20*i1p.ibasis0+9] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 2);

    // Basis function derived toward xxx
    work_cart[20*i1p.ibasis0+10] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 3, 0, 0);

    // Basis function derived toward xxy
    work_cart[20*i1p.ibasis0+11] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 1, 0);

    // Basis function derived toward xxz
    work_cart[20*i1p.ibasis0+12] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 1);

    // Basis function derived toward xyy
    work_cart[20*i1p.ibasis0+13] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 2, 0);

    // Basis function derived toward xyz
    work_cart[20*i1p.ibasis0+14] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 1);

    // Basis function derived toward xzz
    work_cart[20*i1p.ibasis0+15] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 2);

    // Basis function derived toward yyy
    work_cart[20*i1p.ibasis0+16] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 3, 0);

    // Basis function derived toward yyz
    work_cart[20*i1p.ibasis0+17] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 1);

    // Basis function derived toward yzz
    work_cart[20*i1p.ibasis0+18] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 2);

    // Basis function derived toward zzz
    work_cart[20*i1p.ibasis0+19] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 3);
  } while (i1p.inc());
}

void GB1DMGridNabla3Fn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The density Nabla3 is computed in `point` using the results in work_basis (basis
  // function values and their first and second derivatives).
  double rho_x = 0, rho_y = 0, rho_z = 0;
  double rho_xx, rho_yy, rho_zz;
  double nabla3_x_ao, nabla3_y_ao, nabla3_z_ao;
  double lap_ao;
  //XHC Atomic Orbital b
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;
    //XHC Atomic Orbital a
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      //XHC row only has one value (scalar) and tmp has 3 (gradient is vector)
      //XHC P_ab*AO_a
      row += work_basis[ibasis1*20]*dm[ibasis0*nbasis+ibasis1];
      //XHC P_ab*grad_AO_a
      tmp_x += work_basis[ibasis1*20+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*20+2]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*20+3]*dm[ibasis0*nbasis+ibasis1];
      //XHC Third term P_ab*grad(grad_AO_a*grad_AO_b)
      rho_xx = work_basis[ibasis0*20+4]*work_basis[ibasis1*20+1];
      rho_xx += work_basis[ibasis1*20+4]*work_basis[ibasis0*20+1];
      rho_xx += work_basis[ibasis0*20+5]*work_basis[ibasis1*20+2];
      rho_xx += work_basis[ibasis1*20+5]*work_basis[ibasis0*20+2];
      rho_xx += work_basis[ibasis0*20+6]*work_basis[ibasis1*20+3];
      rho_xx += work_basis[ibasis1*20+6]*work_basis[ibasis0*20+3];
      rho_x += rho_xx*dm[ibasis0*nbasis+ibasis1];
      rho_yy = work_basis[ibasis0*20+5]*work_basis[ibasis1*20+1];
      rho_yy += work_basis[ibasis1*20+5]*work_basis[ibasis0*20+1];
      rho_yy += work_basis[ibasis0*20+7]*work_basis[ibasis1*20+2];
      rho_yy += work_basis[ibasis1*20+7]*work_basis[ibasis0*20+2];
      rho_yy += work_basis[ibasis0*20+8]*work_basis[ibasis1*20+3];
      rho_yy += work_basis[ibasis1*20+8]*work_basis[ibasis0*20+3];
      rho_y += rho_yy*dm[ibasis0*nbasis+ibasis1];
      rho_zz = work_basis[ibasis0*20+6]*work_basis[ibasis1*20+1];
      rho_zz += work_basis[ibasis1*20+6]*work_basis[ibasis0*20+1];
      rho_zz += work_basis[ibasis0*20+8]*work_basis[ibasis1*20+2];
      rho_zz += work_basis[ibasis1*20+8]*work_basis[ibasis0*20+2];
      rho_zz += work_basis[ibasis0*20+9]*work_basis[ibasis1*20+3];
      rho_zz += work_basis[ibasis1*20+9]*work_basis[ibasis0*20+3];
      rho_z += rho_zz*dm[ibasis0*nbasis+ibasis1];
    } 
    //XHC First term P_ab*AO_a*nabla3_AO_b
    nabla3_x_ao = work_basis[ibasis0*20+10];
    nabla3_x_ao += work_basis[ibasis0*20+13];
    nabla3_x_ao += work_basis[ibasis0*20+15];
    rho_x += row*nabla3_x_ao;
    nabla3_y_ao = work_basis[ibasis0*20+11];
    nabla3_y_ao += work_basis[ibasis0*20+16];
    nabla3_y_ao += work_basis[ibasis0*20+18];
    rho_y += row*nabla3_y_ao;
    nabla3_z_ao = work_basis[ibasis0*20+12];
    nabla3_z_ao += work_basis[ibasis0*20+17];
    nabla3_z_ao += work_basis[ibasis0*20+19];
    rho_z += row*nabla3_z_ao;
    
    //XHC Second term P_ab*grad_AO_a*lap_AO_b
    lap_ao = work_basis[ibasis0*20+4] + work_basis[ibasis0*20+7] + work_basis[ibasis0*20+9];
    rho_x += tmp_x*lap_ao;
    rho_y += tmp_y*lap_ao;
    rho_z += tmp_z*lap_ao;
 
  }
  output[0] += 2*rho_x;
  output[1] += 2*rho_y;
  output[2] += 2*rho_z;
  //output[3] += 2*rho_yy;
  //output[4] += 2*rho_yz;
  //output[5] += 2*rho_zz;
  //printf("Hola\n");
}


void GB1DMGridNabla3Fn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the density Hessian matrix elements, evaluated in
  // `point` is given in `*pot` (six elements). These are transformed with the chain rule
  // to a contribution to the Fock matrix from `point`.
  /*for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp = pot[0]*work_basis[ibasis0*10+4]
                +pot[1]*work_basis[ibasis0*10+5]
                +pot[2]*work_basis[ibasis0*10+6]
                +pot[3]*work_basis[ibasis0*10+7]
                +pot[4]*work_basis[ibasis0*10+8]
                +pot[5]*work_basis[ibasis0*10+9];
    double tmp_x = pot[0]*work_basis[ibasis0*10+1]
                  +pot[1]*work_basis[ibasis0*10+2]
                  +pot[2]*work_basis[ibasis0*10+3];
    double tmp_y = pot[1]*work_basis[ibasis0*10+1]
                  +pot[3]*work_basis[ibasis0*10+2]
                  +pot[4]*work_basis[ibasis0*10+3];
    double tmp_z = pot[2]*work_basis[ibasis0*10+1]
                  +pot[4]*work_basis[ibasis0*10+2]
                  +pot[5]*work_basis[ibasis0*10+3];
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      double result = tmp*work_basis[ibasis1*10]
                     +tmp_x*work_basis[ibasis1*10+1]
                     +tmp_y*work_basis[ibasis1*10+2]
                     +tmp_z*work_basis[ibasis1*10+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      fock[ibasis0*nbasis+ibasis1] += result;
    }
  }*/
}


/*
    GB1DMGridNabla4Fn
*/

void GB1DMGridNabla4Fn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h4 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+4)+1;
  offset_h3 = offset_h4 - ((abs(shell_type0)+4)*(abs(shell_type0)+5))/2;
  offset_h2 = offset_h3 - ((abs(shell_type0)+3)*(abs(shell_type0)+4))/2;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
  offset_l3 = offset_l2 - ((abs(shell_type0)-2)*(abs(shell_type0-1)))/2;
  offset_l4 = offset_l3 - ((abs(shell_type0)-3)*(abs(shell_type0-2)))/2;
}

void GB1DMGridNabla4Fn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The value of the basis functions in the Cartesian primitive shell are added to
  // work_basis, together with its first and second derivatives, all evaluated at `point`.
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    //double pre0_h = -pre0*2.0*alpha0;
    //double pre0_hh = -pre0_h*2.0*alpha0;
    //long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[35*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[35*i1p.ibasis0+1] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 0);

    // Basis function derived toward y
    work_cart[35*i1p.ibasis0+2] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 0);

    // Basis function derived toward z
    work_cart[35*i1p.ibasis0+3] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 1);

    // Basis function derived toward xx
    work_cart[35*i1p.ibasis0+4] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 0);

    // Basis function derived toward xy
    work_cart[35*i1p.ibasis0+5] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 0);

    // Basis function derived toward xz
    work_cart[35*i1p.ibasis0+6] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 1);

    // Basis function derived toward yy
    work_cart[35*i1p.ibasis0+7] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 0);

    // Basis function derived toward yz
    work_cart[35*i1p.ibasis0+8] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 1);

    // Basis function derived toward zz
    work_cart[35*i1p.ibasis0+9] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 2);

    // Basis function derived toward xxx
    work_cart[35*i1p.ibasis0+10] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 3, 0, 0);

    // Basis function derived toward xxy
    work_cart[35*i1p.ibasis0+11] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 1, 0);

    // Basis function derived toward xxz
    work_cart[35*i1p.ibasis0+12] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 1);

    // Basis function derived toward xyy
    work_cart[35*i1p.ibasis0+13] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 2, 0);

    // Basis function derived toward xyz
    work_cart[35*i1p.ibasis0+14] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 1);

    // Basis function derived toward xzz
    work_cart[35*i1p.ibasis0+15] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 2);

    // Basis function derived toward yyy
    work_cart[35*i1p.ibasis0+16] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 3, 0);

    // Basis function derived toward yyz
    work_cart[35*i1p.ibasis0+17] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 1);

    // Basis function derived toward yzz
    work_cart[35*i1p.ibasis0+18] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 2);

    // Basis function derived toward zzz
    work_cart[35*i1p.ibasis0+19] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 3);

    // Basis function derived toward xxxx
    work_cart[35*i1p.ibasis0+20] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 4, 0, 0);

    // Basis function derived toward xxxy
    work_cart[35*i1p.ibasis0+21] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 3, 1, 0);

    // Basis function derived toward xxxz
    work_cart[35*i1p.ibasis0+22] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 3, 0, 1);

    // Basis function derived toward xxyy
    work_cart[35*i1p.ibasis0+23] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 2, 0);

    // Basis function derived toward xxyz
    work_cart[35*i1p.ibasis0+24] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 1, 1);

    // Basis function derived toward xxzz
    work_cart[35*i1p.ibasis0+25] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 2, 0, 2);

    // Basis function derived toward xyyy
    work_cart[35*i1p.ibasis0+26] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 3, 0);

    // Basis function derived toward xyyz
    work_cart[35*i1p.ibasis0+27] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 2, 1);

    // Basis function derived toward xyzz
    work_cart[35*i1p.ibasis0+28] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 1, 2);

    // Basis function derived toward xzzz
    work_cart[35*i1p.ibasis0+29] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 1, 0, 3);

    // Basis function derived toward yyyy
    work_cart[35*i1p.ibasis0+30] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 4, 0);

    // Basis function derived toward yyyz
    work_cart[35*i1p.ibasis0+31] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 3, 1);

    // Basis function derived toward yyzz
    work_cart[35*i1p.ibasis0+32] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 2, 2);

    // Basis function derived toward yzzz
    work_cart[35*i1p.ibasis0+33] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 1, 3);

    // Basis function derived toward zzzz
    work_cart[35*i1p.ibasis0+34] += dir_gaussian_derivative(r0, point, i1p.n0, alpha0, coeff*scales0[i1p.ibasis0], 0, 0, 4);
  } while (i1p.inc());
}

void GB1DMGridNabla4Fn::compute_point_from_dm(double* work_basis, double* dm,
                                               long nbasis, double* output,
                                               double epsilon, double* dmmaxrow) {
  // The density Nabla4 is computed in `point` using the results in work_basis (basis
  // function values and their first and second derivatives).
  double nabla4_ao, nabla4_rho = 0, lap_ao;
  double nabla3_x_ao, nabla3_y_ao, nabla3_z_ao;
  //XHC Atomic Orbital b
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp1 = 0;
    double tmp2_x = 0, tmp2_y = 0, tmp2_z = 0;
    double tmp3 = 0;
    double tmp4_x = 0, tmp4_y = 0, tmp4_z = 0;
    double tmp4_xx = 0, tmp4_xy= 0, tmp4_xz = 0;
    double tmp4_yy = 0, tmp4_yz = 0, tmp4_zz = 0;
    //XHC Atomic Orbital a
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      //XHC row only has one value (scalar) and tmp has 3 (gradient is vector)
      //XHC First term P_ab*AO_a
      tmp1 += work_basis[ibasis1*35]*dm[ibasis0*nbasis+ibasis1];
      //XHC Second term P_ab*grad_AO_a
      tmp2_x += work_basis[ibasis1*35+1]*dm[ibasis0*nbasis+ibasis1];
      tmp2_y += work_basis[ibasis1*35+2]*dm[ibasis0*nbasis+ibasis1];
      tmp2_z += work_basis[ibasis1*35+3]*dm[ibasis0*nbasis+ibasis1];
      //XHC Third term P_ab*lap_AO_a
      lap_ao = work_basis[ibasis1*35+4] + work_basis[ibasis1*35+7] + work_basis[ibasis1*35+9];
      tmp3 += lap_ao;
      //XHC Fourth term P_ab*lap(grad_AO_a*grad_AO_b)
      tmp4_x += work_basis[ibasis1*35+10];
      tmp4_x += work_basis[ibasis1*35+13];
      tmp4_x += work_basis[ibasis1*35+15];
      tmp4_y += work_basis[ibasis1*35+11];
      tmp4_y += work_basis[ibasis1*35+16];
      tmp4_y += work_basis[ibasis1*35+18];
      tmp4_z += work_basis[ibasis1*35+12];
      tmp4_z += work_basis[ibasis1*35+17];
      tmp4_z += work_basis[ibasis1*35+19];
      tmp4_xx += work_basis[ibasis1*35+4];
      tmp4_xy += 2*work_basis[ibasis1*35+5];
      tmp4_xz += 2*work_basis[ibasis1*35+6];
      tmp4_yy += 2*work_basis[ibasis1*35+7];
      tmp4_yz += 2*work_basis[ibasis1*35+8];
      tmp4_zz += work_basis[ibasis1*35+9];
    } 
 
    //XHC First term P_ab*AO_a*nabla4_AO_b
    nabla4_ao = work_basis[ibasis0*35+20];
    nabla4_ao += 2*work_basis[ibasis0*35+23];
    nabla4_ao += 2*work_basis[ibasis0*35+25];
    nabla4_ao += work_basis[ibasis0*35+30];
    nabla4_ao += 2*work_basis[ibasis0*35+32];
    nabla4_ao += work_basis[ibasis0*35+34];
    nabla4_rho += tmp1*nabla4_ao;

    //XHC Second term P_ab*grad_AO_a*nabla3_AO_b
    nabla3_x_ao = work_basis[ibasis0*35+10];
    nabla3_x_ao += work_basis[ibasis0*35+13];
    nabla3_x_ao += work_basis[ibasis0*35+15];
    nabla4_rho += 2*tmp2_x*nabla3_x_ao;
    nabla3_y_ao = work_basis[ibasis0*35+11];
    nabla3_y_ao += work_basis[ibasis0*35+16];
    nabla3_y_ao += work_basis[ibasis0*35+18];
    nabla4_rho += 2*tmp2_y*nabla3_y_ao;
    nabla3_z_ao = work_basis[ibasis0*35+12];
    nabla3_z_ao += work_basis[ibasis0*35+17];
    nabla3_z_ao += work_basis[ibasis0*35+19];
    nabla4_rho += 2*tmp2_z*nabla3_z_ao;
    
    //XHC Third term P_ab*lap_AO_a*lap_AO_b
    lap_ao = work_basis[ibasis0*20+4] + work_basis[ibasis0*20+7] + work_basis[ibasis0*20+9];
    nabla4_rho += tmp3*lap_ao;
 
    //XHC Fourth term P_ab*lap(grad_AO_a*grad_AO_b)
    nabla4_rho += 2*tmp4_x*work_basis[ibasis0*35+1];
    nabla4_rho += 2*tmp4_y*work_basis[ibasis0*35+2];
    nabla4_rho += 2*tmp4_z*work_basis[ibasis0*35+3];
    nabla4_rho += 2*tmp4_xx*work_basis[ibasis0*35+4];
    nabla4_rho += 2*tmp4_xy*work_basis[ibasis0*35+5];
    nabla4_rho += 2*tmp4_xz*work_basis[ibasis0*35+6];
    nabla4_rho += 2*tmp4_yy*work_basis[ibasis0*35+7];
    nabla4_rho += 2*tmp4_yz*work_basis[ibasis0*35+8];
    nabla4_rho += 2*tmp4_zz*work_basis[ibasis0*35+9];
  }
  *output += 2*nabla4_rho;
  //output[1] += 2*rho_xy;
  //output[2] += 2*rho_xz;
  //output[3] += 2*rho_yy;
  //output[4] += 2*rho_yz;
  //output[5] += 2*rho_zz;
  //printf("Hola\n");
}


void GB1DMGridNabla4Fn::compute_fock_from_pot(double* pot, double* work_basis,
                                               long nbasis, double* fock) {
  // The derivative of the energy w.r.t. the density Hessian matrix elements, evaluated in
  // `point` is given in `*pot` (six elements). These are transformed with the chain rule
  // to a contribution to the Fock matrix from `point`.
  /*for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp = pot[0]*work_basis[ibasis0*10+4]
                +pot[1]*work_basis[ibasis0*10+5]
                +pot[2]*work_basis[ibasis0*10+6]
                +pot[3]*work_basis[ibasis0*10+7]
                +pot[4]*work_basis[ibasis0*10+8]
                +pot[5]*work_basis[ibasis0*10+9];
    double tmp_x = pot[0]*work_basis[ibasis0*10+1]
                  +pot[1]*work_basis[ibasis0*10+2]
                  +pot[2]*work_basis[ibasis0*10+3];
    double tmp_y = pot[1]*work_basis[ibasis0*10+1]
                  +pot[3]*work_basis[ibasis0*10+2]
                  +pot[4]*work_basis[ibasis0*10+3];
    double tmp_z = pot[2]*work_basis[ibasis0*10+1]
                  +pot[4]*work_basis[ibasis0*10+2]
                  +pot[5]*work_basis[ibasis0*10+3];
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      double result = tmp*work_basis[ibasis1*10]
                     +tmp_x*work_basis[ibasis1*10+1]
                     +tmp_y*work_basis[ibasis1*10+2]
                     +tmp_z*work_basis[ibasis1*10+3];
      fock[ibasis1*nbasis+ibasis0] += result;
      fock[ibasis0*nbasis+ibasis1] += result;
    }
  }*/
}


/*
    GB1DMGridMGGAFn
*/

void GB1DMGridMGGAFn::reset(long _shell_type0, const double* _r0, const double* _point) {
  GB1GridFn::reset(_shell_type0, _r0, _point);
  poly_work[0] = 1.0;
  poly_work[1] = point[0] - r0[0];
  poly_work[2] = point[1] - r0[1];
  poly_work[3] = point[2] - r0[2];
  // One order higher polynomials are required because of first derivative.
  offset_h2 = fill_cartesian_polynomials(poly_work+1, abs(shell_type0)+2)+1;
  offset_h1 = offset_h2 - ((abs(shell_type0)+2)*(abs(shell_type0)+3))/2;
  offset = offset_h1 - ((abs(shell_type0)+1)*(abs(shell_type0)+2))/2;
  offset_l1 = offset - ((abs(shell_type0))*(abs(shell_type0)+1))/2;
  offset_l2 = offset_l1 - ((abs(shell_type0)-1)*(abs(shell_type0)))/2;
}

void GB1DMGridMGGAFn::add(double coeff, double alpha0, const double* scales0) {
  double pre = coeff*exp(-alpha0*dist_sq(r0, point));
  // The primitive basis functions in one Cartesian shell are added to work_cart,
  // including the first derivatives and the Laplacian. (Five elements in total per basis
  // function.)
  i1p.reset(abs(shell_type0));
  do {
    double pre0 = pre*scales0[i1p.ibasis0];
    double pre0_h = -pre0*2.0*alpha0;
    double pre0_hh = -pre0_h*2.0*alpha0;
    long nnotx = i1p.n0[1] + i1p.n0[2];

    // Basis function
    work_cart[5*i1p.ibasis0] += pre0*poly_work[i1p.ibasis0+offset];

    // Basis function derived toward x
    work_cart[5*i1p.ibasis0+1] += pre0_h*poly_work[i1p.ibasis0+offset_h1];
    if (i1p.n0[0] > 0)
      work_cart[5*i1p.ibasis0+1] += i1p.n0[0]*pre0*poly_work[i1p.ibasis0+offset_l1];

    // Basis function derived toward y
    work_cart[5*i1p.ibasis0+2] += pre0_h*poly_work[i1p.ibasis0+1+nnotx+offset_h1];
    if (i1p.n0[1] > 0)
      work_cart[5*i1p.ibasis0+2] += i1p.n0[1]*pre0*poly_work[i1p.ibasis0-nnotx+offset_l1];

    // Basis function derived toward z
    work_cart[5*i1p.ibasis0+3] += pre0_h*poly_work[i1p.ibasis0+2+nnotx+offset_h1];
    if (i1p.n0[2] > 0)
      work_cart[5*i1p.ibasis0+3] += i1p.n0[2]*pre0*
                                    poly_work[i1p.ibasis0-nnotx-1+offset_l1];

    // Laplacian of the Basis function
    work_cart[5*i1p.ibasis0+4] += pre0_hh*poly_work[i1p.ibasis0+offset_h2];
    work_cart[5*i1p.ibasis0+4] += (2*i1p.n0[0]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[0] > 1)
      work_cart[5*i1p.ibasis0+4] += i1p.n0[0]*(i1p.n0[0]-1)*
                                    pre0*poly_work[i1p.ibasis0+offset_l2];
    work_cart[5*i1p.ibasis0+4] += pre0_hh*poly_work[i1p.ibasis0+3+2*nnotx+offset_h2];
    work_cart[5*i1p.ibasis0+4] += (2*i1p.n0[1]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[1] > 1)
      work_cart[5*i1p.ibasis0+4] += i1p.n0[1]*(i1p.n0[1]-1)*pre0*
                                    poly_work[i1p.ibasis0+1-2*nnotx+offset_l2];
    work_cart[5*i1p.ibasis0+4] += pre0_hh*poly_work[i1p.ibasis0+5+2*nnotx+offset_h2];
    work_cart[5*i1p.ibasis0+4] += (2*i1p.n0[2]+1)*pre0_h*poly_work[i1p.ibasis0+offset];
    if (i1p.n0[2] > 1)
      work_cart[5*i1p.ibasis0+4] += i1p.n0[2]*(i1p.n0[2]-1)*pre0*
                                    poly_work[i1p.ibasis0-1-2*nnotx+offset_l2];
  } while (i1p.inc());
}

void GB1DMGridMGGAFn::compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                            double* output, double epsilon,
                                            double* dmmaxrow) {
  // The density, its gradient, the kinetic energy density and the Laplacian, are computed
  // in `point` (see reset method), using the results in work_basis and the 1RDM
  // coefficients in dm.
  double rho = 0, rho_x = 0, rho_y = 0, rho_z = 0;
  double lapl_part = 0, tau2 = 0.0;
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double row = 0;
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;
    for (long ibasis1=0; ibasis1 < nbasis; ibasis1++) {
      row += work_basis[ibasis1*5]*dm[ibasis0*nbasis+ibasis1];
      tmp_x += work_basis[ibasis1*5+1]*dm[ibasis0*nbasis+ibasis1];
      tmp_y += work_basis[ibasis1*5+2]*dm[ibasis0*nbasis+ibasis1];
      tmp_z += work_basis[ibasis1*5+3]*dm[ibasis0*nbasis+ibasis1];
    }
    rho += row*work_basis[ibasis0*5];
    rho_x += row*work_basis[ibasis0*5+1];
    rho_y += row*work_basis[ibasis0*5+2];
    rho_z += row*work_basis[ibasis0*5+3];
    lapl_part += row*work_basis[ibasis0*5+4];
    tau2 += tmp_x*work_basis[ibasis0*5+1] +
            tmp_y*work_basis[ibasis0*5+2] +
            tmp_z*work_basis[ibasis0*5+3];
  }
  output[0] += rho;
  output[1] += 2*rho_x;
  output[2] += 2*rho_y;
  output[3] += 2*rho_z;
  output[4] += 2*lapl_part + 2*tau2;
  output[5] += 0.5*tau2;
}

void GB1DMGridMGGAFn::compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                            double* fock) {
  // The functional derivative of the energy w.r.t. density, its gradient, the kinetic
  // energy density and the Laplacian, are given in *pot. These are transformed to a
  // contribution to the Fock matrix using the chain rule.
  double auxpot = 0.5*pot[5] + 2.0*pot[4];
  for (long ibasis0=0; ibasis0 < nbasis; ibasis0++) {
    double tmp0 = pot[0]*work_basis[ibasis0*5] +
                  pot[1]*work_basis[ibasis0*5+1] +
                  pot[2]*work_basis[ibasis0*5+2] +
                  pot[3]*work_basis[ibasis0*5+3] +
                  pot[4]*work_basis[ibasis0*5+4];
    double tmp1 = pot[1]*work_basis[ibasis0*5] + auxpot*work_basis[ibasis0*5+1];
    double tmp2 = pot[2]*work_basis[ibasis0*5] + auxpot*work_basis[ibasis0*5+2];
    double tmp3 = pot[3]*work_basis[ibasis0*5] + auxpot*work_basis[ibasis0*5+3];
    double tmp4 = pot[4]*work_basis[ibasis0*5];
    for (long ibasis1=0; ibasis1 <= ibasis0; ibasis1++) {
      double result = tmp0*work_basis[ibasis1*5] +
                      tmp1*work_basis[ibasis1*5+1] +
                      tmp2*work_basis[ibasis1*5+2] +
                      tmp3*work_basis[ibasis1*5+3] +
                      tmp4*work_basis[ibasis1*5+4];
      fock[ibasis1*nbasis+ibasis0] += result;
      if (ibasis1 != ibasis0) {
        // Enforce symmetry
        fock[ibasis0*nbasis+ibasis1] += result;
      }
    }
  }
}


/*
    GB2DMGridFn
*/

GB2DMGridFn::GB2DMGridFn(long max_shell_type)
    : GBCalculator(max_shell_type), shell_type0(0), shell_type1(0), r0(NULL), r1(NULL),
      point(NULL), i2p() {
  nwork = max_nbasis*max_nbasis;
  work_cart = new double[nwork];
  work_pure = new double[nwork];
}

void GB2DMGridFn::reset(long _shell_type0, long _shell_type1, const double* _r0,
                        const double* _r1, const double* _point) {
  if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type)) {
    throw std::domain_error("shell_type0 out of range.");
  }
  if ((_shell_type1 < -max_shell_type) || (_shell_type1 > max_shell_type)) {
    throw std::domain_error("shell_type0 out of range.");
  }
  shell_type0 = _shell_type0;
  shell_type1 = _shell_type1;
  r0 = _r0;
  r1 = _r1;
  point = _point;
  // We make use of the fact that a floating point zero consists of
  // consecutive zero bytes.
  memset(work_cart, 0, nwork*sizeof(double));
  memset(work_pure, 0, nwork*sizeof(double));
}

void GB2DMGridFn::cart_to_pure() {
  /*
     The initial results are always stored in work_cart. The projection
     routine always outputs its result in work_pure. Once that is done,
     the pointers to both blocks are swapped such that the final result is
     always back in work_cart.
  */

  // For now, this is just a copy from GB2Integral.
  // (This must be changed when electrical fields are implemented.)

  // Project along index 0 (rows)
  if (shell_type0 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type0,
      1,                                     // anterior
      get_shell_nbasis(abs(shell_type1)) );  // posterior
    swap_work();
  }

  // Project along index 1 (cols)
  if (shell_type1 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type1,
      get_shell_nbasis(shell_type0),  // anterior
      1);                             // posterior
    swap_work();
  }
}


/*
    GB2DMGridHartreeFn
*/

GB2DMGridHartreeFn::GB2DMGridHartreeFn(long max_shell_type) : GB2DMGridFn(max_shell_type) {
  work_g0 = new double[2*max_shell_type+1];
  work_g1 = new double[2*max_shell_type+1];
  work_g2 = new double[2*max_shell_type+1];
  work_boys = new double[2*max_shell_type+1];
}


GB2DMGridHartreeFn::~GB2DMGridHartreeFn() {
  delete[] work_g0;
  delete[] work_g1;
  delete[] work_g2;
  delete[] work_boys;
}


void GB2DMGridHartreeFn::add(double coeff, double alpha0, double alpha1,
                            const double* scales0, const double* scales1) {
  double pre, gamma, gamma_inv, arg;
  double gpt_center[3], pa[3], pb[3], pc[3];

  gamma = alpha0 + alpha1;
  gamma_inv = 1.0/gamma;
  pre = 2*M_PI*gamma_inv*coeff*exp(-alpha0*alpha1*gamma_inv*dist_sq(r0, r1));
  compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
  pa[0] = gpt_center[0] - r0[0];
  pa[1] = gpt_center[1] - r0[1];
  pa[2] = gpt_center[2] - r0[2];
  pb[0] = gpt_center[0] - r1[0];
  pb[1] = gpt_center[1] - r1[1];
  pb[2] = gpt_center[2] - r1[2];

  // thrid center for the current charge
  pc[0] = gpt_center[0] - point[0];
  pc[1] = gpt_center[1] - point[1];
  pc[2] = gpt_center[2] - point[2];

  // Fill the work array with the Boys function values
  arg = gamma*(pc[0]*pc[0] + pc[1]*pc[1] + pc[2]*pc[2]);
  for (long nu=abs(shell_type0)+abs(shell_type1); nu >= 0; nu--) {
    work_boys[nu] = boys_function(nu, arg);
  }

  // Iterate over all combinations of Cartesian exponents
  i2p.reset(abs(shell_type0), abs(shell_type1));
  do {
    // Fill the work arrays with the polynomials
    nuclear_attraction_helper(work_g0, i2p.n0[0], i2p.n1[0], pa[0], pb[0], pc[0], gamma_inv);
    nuclear_attraction_helper(work_g1, i2p.n0[1], i2p.n1[1], pa[1], pb[1], pc[1], gamma_inv);
    nuclear_attraction_helper(work_g2, i2p.n0[2], i2p.n1[2], pa[2], pb[2], pc[2], gamma_inv);

    // Take the product
    arg = 0;
    for (long i0 = i2p.n0[0]+i2p.n1[0]; i0 >= 0; i0--)
      for (long i1 = i2p.n0[1]+i2p.n1[1]; i1 >= 0; i1--)
        for (long i2 = i2p.n0[2]+i2p.n1[2]; i2 >= 0; i2--)
          arg += work_g0[i0]*work_g1[i1]*work_g2[i2]*work_boys[i0+i1+i2];

    // Finally add to the work array
    work_cart[i2p.offset] += pre*scales0[i2p.ibasis0]*scales1[i2p.ibasis1]*arg;
  } while (i2p.inc());
}
