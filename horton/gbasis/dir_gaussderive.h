//#include <cstdio>
//#include <cmath>
#include "horton/gbasis/common.h"


double dir_gaussian_derivative(const double* r, const double* atom, long int* n, double alpha0, const double scales0, int derx, int dery, int derz);
//double dir_gaussian_derivative(double* poly_work, long int* n, double alpha, const double scales0, int derx, int dery, int derz);
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
