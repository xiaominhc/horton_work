#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include "horton/gbasis/dir_gaussderive.h"

double dir_gaussian_derivative(const double* r, const double* atom, long int* n, double alpha0, const double scales0, int derx, int dery, int derz){
//double dir_gaussian_derivative(double* r, double* atom, long int* n, double alpha0, const double scales0, int derx, int dery, int derz){
    /*r"""
    This function takes the specified derivative of a gaussian function and returns the value of the derivative at a point of the grid.

    Parameters
    ----------
    r : np.ndarray
        Coordinates of the grid point: ``x``, ``y``, and ``z``.
    atom : np.ndarray
        Coordinates of the nucleus: ``A_x``, ``A_y``, and ``A_z``.
    n : np.ndarray
        Angular momentum of orbital: ``n_x``, ``n_y``, and ``n_z``.
    alpha : real
        Exponent of the gaussian function: ``\alpha_i``.
    scales0 : real
        Contraction coefficient: ``C_i``.
    order : integer
        Total order of the derivative.
    derx : integer
        Order of the derivative with respect to x.
    dery : integer
        Order of the derivative with respect to y.
    derz : integer
        Order of the derivative with respect to z.

    Returns
    -------
    result: real
        Value of the specified derivative of the gaussian function.
    """*/
    int order = derx + dery + derz;
 
    // Evaluate distance between point and nucleus
    double poly[4];
    /*poly[0] = poly_work[0];
    poly[1] = poly_work[1];
    poly[2] = poly_work[2];
    poly[3] = poly_work[3];*/
    poly[0] = 1.0;
    poly[1] = atom[0] - r[0];
    poly[2] = atom[1] - r[1];
    poly[3] = atom[2] - r[2];
    /*printf("reset\n");
    printf("poly_work[0] = %f\n", poly[0]);
    printf("poly_work[1] = %f\n", poly[1]);
    printf("poly_work[2] = %f\n", poly[2]);
    printf("poly_work[3] = %f\n", poly[3]);*/



    // Evaluate exponential part multiplied by the contraction coefficient
    double pre0 = exp(-alpha0*dist_sq(r, atom))*scales0;
    //printf("exp = %f\n", exp(-alpha0*dist_sq(r, atom)));
    //printf("scales = %f\n", scales0);
    //printf("pre0 = %f\n", pre0);
 
    // Evaluate exponential part for +1 order in derivative
    double pre0_h = -2.0*alpha0;
    
    // Number of elements needed for work_cart array = 2^0 + 2^1 +...+ 2^order
    // For each derivative 2^order gaussians are generated
    // Lower order dimensions are calculated for offset calculations
    int workdim[order+1];
    workdim[0] = 1; 
    for (int iorder=1; iorder<order+1; iorder++){
        workdim[iorder] = workdim[iorder-1] + pow(2,iorder);
        //print "workdim[", iorder, "]=", workdim[iorder]
    }
    //printf("wasabi");
   
    // Evaluate polynomial part of gaussian function
    double tmp;
    tmp = pow(poly[1], n[0]);
    tmp *= pow(poly[2], n[1]);
    tmp *= pow(poly[3], n[2]);

    /*"""
    Offsets are needed to point to the function to derive and where to store the
    new 2 gaussian functions in work_cart. The derivatives branch out like a tree. 
    Offset[0] points to the element in the middle of the tree (the original function). 
    
    Case A (order = 1):
 
              d/di(f(x))_high
            /     2
        f(x)
          1 \
              d/di(f(x))_low
                  0
        Total number of elements in work_cart = 3
        Offset[0] = 1
        Offset[1] points to the middle + 1 (higher angular momentum term of the first derivative).
        Offset[1] = 2
        Offset[2] points to the middle - 1 (lower angular momentum term of the first derivative).
        Offset[1] = 0

    Case A (order = 2):

                              d/di(d/dj(f(x)))_high2
                             /       6
              d/di(f(x))_high1
             /     5         \
            /                 d/di(d/dj(f(x)))_zero
        f(x)                         4
         3  \                 d/di(d/dj(f(x)))_zero
             \               /       2
              d/di(f(x))_low1
                   1         \
                              d/di(d/dj(f(x)))_low2
                                     0
        Offset[0] = 3
        Offset[1] = 5
        Offset[1] points to the middle + 2 (higher angular momentum term of the first derivative).
        Offset[2] = 1
        Offset[2] points to the middle - 2 (lower angular momentum term of the first derivative).
        Offset[3] = 6
        Offset[3] points to the middle + 2 + 1.
        Offset[4] = 4
        Offset[4] points to the middle + 2 - 1.
        Offset[5] = 2
        Offset[5] points to the middle - 2 + 1.
        Offset[6] = 0
        Offset[5] points to the middle - 2 - 1.
    """*/
    int offset[workdim[order]];
    offset[0] = (workdim[order] - 1)/2; // Pointer to the middle of the tree
    
    // Evaluate original gaussian function
    double work_cart[workdim[order]];
    //print "offset[0]=", offset[0]
    //print "workdim[", order, "]=", workdim[order]
    work_cart[offset[0]] = pre0*tmp;

    // Auxiliary angular momentum arrays that keep track of the modified n
    int tmp_nx[workdim[order]];
    int tmp_ny[workdim[order]];
    int tmp_nz[workdim[order]];
    tmp_nx[0] = n[0];
    tmp_ny[0] = n[1];
    tmp_nz[0] = n[2];
 
    // Iterate over derivatives of gaussian function
    int idx = 0;
    for (int i = 0; i < order; i++) {
      //print "wasabi i", i+1, " of order", order
      int nelem = workdim[order - i - 1];
      for (int j = 0; j < pow(2,i); j++) {
        //print "   "
        //print "wasabi j", j+1, " of workdim[i-1]+1", workdim[i]+1
        //print "wasabi idx", idx, "offset", offset[idx]
        //print "wasabi i=", i,"wasabi j=", j, "wasabi idx=", idx, "offsetidx=", offset[idx]
        // Derive toward x
        if (derx > i) {
          //print "entra a x"
          offset[2*idx+1] = offset[idx] + (nelem + 1)/2;
          tmp_nx[2*idx+1] = tmp_nx[idx] + 1;
          //print "offset[", 2*idx+1, "]=", offset[2*idx+1]
          offset[2*idx+2] = offset[idx] - (nelem + 1)/2;
          tmp_nx[2*idx+2] = tmp_nx[idx] - 1;
          //print "offset[", 2*idx+2, "]=", offset[2*idx+2]
          work_cart[offset[2*idx+1]] = pre0_h*poly[1]*work_cart[offset[idx]];
          if (poly[1] == 0.0) {
          work_cart[offset[2*idx+2]] = 0.0;
          }
          else {
          work_cart[offset[2*idx+2]] = tmp_nx[idx]/poly[1]*work_cart[offset[idx]];
          }
          tmp_ny[2*idx+1] = tmp_ny[idx];
          tmp_ny[2*idx+2] = tmp_ny[idx];
          tmp_nz[2*idx+1] = tmp_nz[idx];
          tmp_nz[2*idx+2] = tmp_nz[idx];
        }
        // Derive toward y
        else if (dery>(i-derx)) {
          //print "entra a y"
          offset[2*idx+1] =offset[idx] + (nelem + 1)/2;
          tmp_ny[2*idx+1] = tmp_ny[idx] + 1;
          //print "offset[", 2*idx+1, "]=", offset[2*idx+1]
          offset[2*idx+2] =offset[idx] - (nelem + 1)/2;
          tmp_ny[2*idx+2] = tmp_ny[idx] - 1;
          //print "offset[", 2*idx+2, "]=", offset[2*idx+2]
          work_cart[offset[2*idx+1]] = pre0_h*poly[2]*work_cart[offset[idx]];
          if (poly[2] == 0.0) {
          work_cart[offset[2*idx+2]] = 0.0;
          }
          else {
          work_cart[offset[2*idx+2]] = tmp_ny[idx]/poly[2]*work_cart[offset[idx]];  
          }
          tmp_nx[2*idx+1] = tmp_nx[idx];
          tmp_nx[2*idx+2] = tmp_nx[idx];
          tmp_nz[2*idx+1] = tmp_nz[idx];
          tmp_nz[2*idx+2] = tmp_nz[idx];
        }
        // Derive toward z
        else {
          //print "entra a z"
          offset[2*idx+1] =offset[idx] + (nelem + 1)/2;
          tmp_nz[2*idx+1] = tmp_nz[idx] + 1;
          //print "offset[", 2*idx+1, "]=", offset[2*idx+1]
          offset[2*idx+2] =offset[idx] - (nelem + 1)/2;
          tmp_nz[2*idx+2] = tmp_nz[idx] - 1;
          //print "offset[", 2*idx+2, "]=", offset[2*idx+2]
          work_cart[offset[2*idx+1]] = pre0_h*poly[3]*work_cart[offset[idx]];
          if (poly[3] == 0.0) {
          work_cart[offset[2*idx+2]] = 0.0;
          }
          else {
          work_cart[offset[2*idx+2]] = tmp_nz[idx]/poly[3]*work_cart[offset[idx]];
          }
          tmp_nx[2*idx+1] = tmp_nx[idx];
          tmp_nx[2*idx+2] = tmp_nx[idx];
          tmp_ny[2*idx+1] = tmp_ny[idx];
          tmp_ny[2*idx+2] = tmp_ny[idx];
        }
        //printf("work_cart_in[%d]=%f\n",offset[2*idx+1],work_cart[offset[2*idx+1]]);
        //printf("work_cart_in[%d]=%f\n",offset[2*idx+2],work_cart[offset[2*idx+2]]);
        idx += 1;
      }
    }

 
    // Add terms to obtain actual derivatives (Only the terms in the last branch contribute).
    double result = 0.0;
    for (int i = 0; i < pow(2,order); i++) {
      result +=work_cart[offset[workdim[order]-1-i]];
      //print "result", result, ",i", i, ",offset[", workdim[order]-1-i,"]=", offset[workdim[order]-1-i]

    }

    // If order is zero, then return the value of the function
    if (order==0) {
      result = work_cart[offset[0]];
    }
        
    return result;
}
