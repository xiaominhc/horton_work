# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


from scipy.special import erf
import numpy as np

from horton import *


def test_solve_poisson_becke_n2():
    mol = Molecule.from_file(context.get_fn('test/n2_hfs_sto3g.fchk'))
    lmaxmax = 4

    # compute hartree potential on a molecular grid
    molgrid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False, mode='keep')
    dm_full = mol.get_dm_full()
    reference = mol.obasis.compute_grid_hartree_dm(dm_full, molgrid.points)

    # construct the same potential numerically with Becke's method
    rho = mol.obasis.compute_grid_density_dm(dm_full, molgrid.points)
    begin = 0
    hds = []
    for i in xrange(mol.natom):
        atgrid = molgrid.subgrids[i]
        end = begin + atgrid.size
        becke_weights = molgrid.becke_weights[begin:end]
        density_decomposition = atgrid.get_spherical_decomposition(rho[begin:end], becke_weights, lmax=lmaxmax)
        hartree_decomposition = solve_poisson_becke(density_decomposition)
        hds.append(hartree_decomposition)
        begin = end

    # Evaluate the splines obtained with Becke's method on the molecular grid
    # Increasing angular momenta are used to check the convergence.
    last_error = None
    for lmax in xrange(0, lmaxmax+1):
        result = molgrid.zeros()
        for i in xrange(mol.natom):
            molgrid.eval_decomposition(hds[i][:(lmax+1)**2], mol.coordinates[i], result)
        potential_error = result - reference
        error = molgrid.integrate(potential_error, potential_error)**0.5
        if last_error is not None:
            assert error < last_error
        last_error = error
        if False:
            worst = molgrid.integrate(reference, reference)**0.5
            print 'lmax=%i  %12.4e  %12.4e' % (lmax, error, worst)
            for rho_low, rho_high in (0, 1e-8), (1e-8, 1e-4), (1e-4, 1e0), (1e0, 1e4), (1e4, 1e100):
                mask = ((rho >= rho_low) & (rho < rho_high)).astype(float)
                error = molgrid.integrate(potential_error, potential_error, mask)**0.5
                worst = molgrid.integrate(reference, reference, mask)**0.5
                print '%10.2e : %10.2e   |   %12.4e  %12.4e' % (rho_low, rho_high, error, worst)
            print
    assert error < 6e-2

    if False:
        # Plot stuff
        import matplotlib.pyplot as pt
        linegrid = LineGrid(mol.coordinates[0], mol.coordinates[1], 500, 1)
        rho = mol.obasis.compute_grid_density_dm(dm_full, linegrid.points)
        reference = mol.obasis.compute_grid_hartree_dm(dm_full, linegrid.points)
        for lmax in xrange(0, lmaxmax+1):
            result = linegrid.zeros()
            for i in xrange(mol.natom):
                linegrid.eval_decomposition(hds[i][:(lmax+1)**2], mol.coordinates[i], result)
            pt.clf()
            #pt.plot(linegrid.x, reference)
            #pt.plot(linegrid.x, result)
            pt.plot(linegrid.x, (result - reference))
            pt.ylim(-0.3, 0.3)
            pt.savefig('test_poisson_%i.png' % lmax)


def test_solve_poisson_becke_sa():
    sigma = 8.0
    rtf = ExpRTransform(4e-4, 4e1, 100)
    r = rtf.get_radii()
    rhoy = np.exp(-0.5*(r/sigma)**2)/sigma**3/(2*np.pi)**1.5
    rhod = np.exp(-0.5*(r/sigma)**2)/sigma**3/(2*np.pi)**1.5*(-r/sigma)/sigma
    rho = CubicSpline(rhoy, rhod, rtf)
    v = solve_poisson_becke([rho])[0]

    s2s = np.sqrt(2)*sigma
    soly = erf(r/s2s)/r
    sold = np.exp(-(r/s2s)**2)*2/np.sqrt(np.pi)/s2s/r - erf(r/s2s)/r**2

    if False:
        import matplotlib.pyplot as pt
        n = 10
        pt.clf()
        pt.plot(r[:n], soly[:n], label='exact')
        pt.plot(r[:n], v.y[:n], label='spline')
        pt.legend(loc=0)
        pt.savefig('denu.png')

    assert abs(v.y - soly).max()/abs(soly).max() < 1e-6
    assert abs(v.dx - sold).max()/abs(sold).max() < 1e-4