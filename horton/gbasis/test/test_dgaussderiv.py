# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


import numpy as np
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_gauss_derivative_s():
    # reference data computed with sympy.mpmath arbitrary precision libraray
   
    r = np.array([0.0, -0.1, 1.2]) 
    atom = np.array([0.0, 1.1, 0.2]) 
    n = np.array([0, 0, 0]) 
    alpha0 = 2.5 
    scales0 = 1.0
    result = np.array([
        # Function
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 0, 0),
        # First derivative of the function
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 1, 0, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 1, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 0, 1),
        # Second derivative of the function
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 2, 0, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 1, 1, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 1, 0, 1),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 2, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 1, 1),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 0, 2),
        # Third derivative of the function
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 3, 0, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 2, 1, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 2, 0, 1),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 1, 2, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 1, 1, 1),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 1, 0, 2),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 3, 0),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 2, 1),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 1, 2),
        dir_gaussian_derivative(r, atom, n, alpha0, scales0, 0, 0, 3)])
    check = np.array([
        # Function
        0.00224287,
        # First derivative of the function
        0.00000000, -0.0134572, 0.0112143,
        # Second derivative of the function
        -0.0112143, 0.00000000, 0.00000000, 0.0695289, -0.067286, 0.0448574,
        # Third derivative of the function
        0.00000000, 0.067286, -0.0560717, 0.0000000, 0.00000000,
        0.00000000, -0.282601, 0.347644, -0.269144, 0.112143])

    # relative errors for the bigger ones
    #print result
    #print check 
    big_mask = abs(check) > 1e-10
    errors = (result[big_mask] - check[big_mask])/check[big_mask]
    max_error = abs(errors).max()
    rms_error = np.sqrt((errors**2).mean())
    assert max_error < 1e-5
    assert rms_error < 1e-5

#    # absolute errors for the smaller ones
#    small_mask = abs(check) > 1e-10
#    errors = result[small_mask] - check[small_mask]
#    max_error = abs(errors).max()
#    rms_error = np.sqrt((errors**2).mean())
#    assert max_error < 1e-16
#    assert rms_error < 1e-16


#def test_boys_domain_error():
#    for m, t in (-1, 0.0), (get_max_shell_type()*4+1, 0.0), (5, -1):
#        with assert_raises(ValueError):
#            boys_function(m, t)


#def test_boys_array():
#    for mmax in xrange(get_max_shell_type()*4+1):
#        for t in np.random.uniform(0, 200, 500):
#            output = boys_function_array(mmax, t)
#            for m in xrange(mmax+1):
#                assert output[m] == boys_function(m, t)
