# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import os
import sys
from tools import vrrotvec2mat

class Lattice(object):
    """
    This class contains all the crystallographic information required for each
    atom type. Currently there are only two pre-configured atoms available in
    class i.e. 'Al' and 'Mg'. Up on need user can create a new instance of this
    class with the same attributes. The attributes of this class are:
    ...

    Attributes
    ----------
    elem_type: string
        Element of Interest

    pearson: string
        Pearson symbol for the lattice

    lat_params: dictionary
        Lattice parameters ('a', 'b', 'c', 'alpha', 'beta', 'gamma')

    l_p_po: numpy array
        Primitve basis of the lattice

    basis_atoms:
        Location of the basis atoms in the primitive lattice

    cryst_ptgrp: string
        Crystallographic point group of the lattice

    burgers_mag: float
        The smallest burgers vector in the lattice

    eam_file:
        eam_file name for atomistic simulations

    Methods
    -------
    str
    Method for printing the lattice class

    Notes
    --------
    Examples of elem_type
    elem_type = 'Mg'; \v
    elem_type = 'Al'; \v
    elem_type = 'Cu'; \v
    elem_type = 'Ni'; \v
    elem_type = 'cF_Id'; \v
    elem_type = 'cI_Id'; \v
    elem_type = 'cP_Id'; \v
    elem_type = 'hP_Id'; \v
    """

    def __init__(self, *args):
        nargs = len(args)
        if nargs == 0:
            elem_type = 'cF_Id'
        else:
            elem_type = args[0]

        # Simple cubic lattice
        # Ideal
        if elem_type == 'cP_Id':
            self.elem_type = 'cP_Id'
            self.pearson = 'cP'
            a = 1.0
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': np.pi/2,
                               'beta': np.pi/2, 'gamma': np.pi/2}
            b1x = a*np.array([1.0, 0.0, 0.0])
            b1y = a*np.array([0.0, 1.0, 0.0])
            b1z = a*np.array([0.0, 0.0, 1.0])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
            # Point group symmetry of the lattice
            self.cryst_ptgrp = 'Oh'
        # ------------------------------------------------------------------------------------------------------
        # Polonium
        if elem_type.lower() == 'po':
            self.elem_type = 'Po'
            self.pearson = 'cP'
            a = 3.35
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': np.pi/2,
                               'beta': np.pi/2, 'gamma': np.pi/2}
            b1x = a*np.array([1.0, 0.0, 0.0])
            b1y = a*np.array([0.0, 1.0, 0.0])
            b1z = a*np.array([0.0, 0.0, 1.0])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
            # Point group symmetry of the lattice
            self.cryst_ptgrp = 'Oh'
        # ------------------------------------------------------------------------------------------------------

        # BCC Lattices
        # Ideal
        if elem_type == 'cI_Id':
            self.elem_type = 'cI_Id'
            self.pearson = 'cI'
            a = 1.0
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': np.pi/2,
                               'beta': np.pi/2, 'gamma': np.pi/2}
            b1x = a*np.array([-0.5,  0.5,  0.5])
            b1y = a*np.array([ 0.5, -0.5,  0.5])
            b1z = a*np.array([ 0.5,  0.5, -0.5])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
            # Point group symmetry of the lattice
            self.cryst_ptgrp = 'Oh'
        # ------------------------------------------------------------------------------------------------------
        # $\alpha$-Fe
        if elem_type.lower() == 'fe_alpha':
            self.elem_type = 'Fe_alpha'
            self.pearson = 'cI'
            a = 2.855324
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': np.pi/2,
                               'beta': np.pi/2, 'gamma': np.pi/2}
            b1x = a*np.array([-0.5,  0.5,  0.5])
            b1y = a*np.array([ 0.5, -0.5,  0.5])
            b1z = a*np.array([ 0.5,  0.5, -0.5])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
            # Point group symmetry of the lattice
            self.cryst_ptgrp = 'Oh'
            self.burgers_mag = a * np.sqrt(3)/ 2
            self.eam_file = np.array(['alloy', 'Fe_2.eam.fs', -4.12243510073005, 'Fe'])
            self.basis_atoms = np.array([0, 0, 0])
        # ------------------------------------------------------------------------------------------------------
        # FCC Lattices
        if elem_type == 'cF_Id':
            self.elem_type = 'cF_Id'
            self.pearson = 'cF'
            a = 1.0
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': np.pi/2,
                               'beta': np.pi/2, 'gamma': np.pi/2}
            b1x = a*np.array([0.0, 0.5, 0.5])
            b1y = a*np.array([0.5, 0.0, 0.5])
            b1z = a*np.array([0.5, 0.5, 0.0])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
            # Point group symmetry of the lattice
            self.cryst_ptgrp = 'Oh'
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'al':
            self.elem_type = 'Al'
            self.pearson = 'cF'

            a = 4.05
            a_ang = np.pi/2
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': a_ang, 'beta': a_ang, 'gamma': a_ang}

            self.cryst_ptgrp = 'Oh'
            self.burgers_mag = a / np.sqrt(2)
            self.eam_file = np.array(['alloy', 'Al99.eam.alloy', -3.36000010077, 'Al'])

            b1x = a*np.array([0.0, 0.5, 0.5])
            b1y = a*np.array([0.5, 0.0, 0.5])
            b1z = a*np.array([0.5, 0.5, 0.0])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))

            self.basis_atoms = np.array([0, 0, 0])
        # ------------------------------------------------------------------------------------------------------
        ## potential used in  doi:10.1016/j.actamat.2009.04.015
        if elem_type.lower() == 'ni':
            self.elem_type = 'Ni'
            self.pearson = 'cF'

            a = 3.52
            a_ang = np.pi/2
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': a_ang, 'beta': a_ang, 'gamma': a_ang}

            self.cryst_ptgrp = 'Oh'
            self.burgers_mag = a / np.sqrt(2)
            self.eam_file = np.array(['alloy', 'ni1.set', -4.45000000526637, 'Ni'])

            b1x = a*np.array([0.0, 0.5, 0.5])
            b1y = a*np.array([0.5, 0.0, 0.5])
            b1z = a*np.array([0.5, 0.5, 0.0])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))

            self.basis_atoms = np.array([0, 0, 0])
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'cu':
            self.elem_type = 'Cu'
            self.pearson = 'cF'

            a = 3.615
            a_ang = np.pi/2
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': a_ang, 'beta': a_ang, 'gamma': a_ang}

            self.cryst_ptgrp = 'Oh'
            self.burgers_mag = a / np.sqrt(2)
            self.eam_file = np.array(['alloy', 'Cu01.eam.alloy', -3.53999996838, 'Cu'])
            # self.eam_file = np.array(['alloy', 'Cu01.eam.alloy', -3.54021833020492]) # 2007 version of Mishin's potential (2001)

            b1x = a*np.array([0.0, 0.5, 0.5])
            b1y = a*np.array([0.5, 0.0, 0.5])
            b1z = a*np.array([0.5, 0.5, 0.0])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))

            self.basis_atoms = np.array([0, 0, 0])

        # Hexagonal Lattices
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'hp_id':
            self.elem_type = 'hP_Id'
            self.pearson = 'hP'
            a = 1.0
            b = a
            CAratio = 1.633
            c = a*CAratio
            a_ang = np.pi/2
            b_ang = np.pi/2
            g_ang = 2*np.pi/3
            self.lat_params = {'a': a, 'b': b, 'c': c, 'alpha': a_ang, 'beta': b_ang, 'gamma': g_ang}
            self.cryst_ptgrp = 'D6h'
            self.burgers_mag = a
            b1x = a*np.array([1., 0., 0.])
            b1y = np.dot(vrrotvec2mat(np.array([0., 0., 1., g_ang])), b1x)
            b1z = c*np.array([0., 0., 1.])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
            self.basis_atoms = np.array([[0., 0., 0.], [1./3, 2./3, 1./2]])
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'hp_ca':
            if nargs < 2:
                raise Exception('Provide the c/a ratio')
            elif nargs == 2:
                    self.elem_type = 'hP_ca'
                    ca_ratio = args[1]
            self.pearson = 'hP'
            a = 1.0
            b = a
            c = a*ca_ratio
            a_ang = np.pi/2
            b_ang = np.pi/2
            g_ang = 2*np.pi/3
            self.lat_params = {'a': a, 'b': b, 'c': c, 'alpha': a_ang, 'beta': b_ang, 'gamma': g_ang}
            self.cryst_ptgrp = 'D6h'
            b1x = a*np.array([1., 0., 0.])
            b1y = np.dot(vrrotvec2mat(np.array([0., 0., 1., g_ang])), b1x)
            b1z = c*np.array([0., 0., 1.])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))

        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'mg':
            self.elem_type = 'Mg'
            self.pearson = 'hP'

            a = 3.181269601; b = a; CAratio = 1.632993162; c = a*CAratio
            a_ang = np.pi/2; b_ang = np.pi/2; g_ang = 2*np.pi/3;
            self.lat_params = {'a': a, 'b': b, 'c': c, 'alpha': a_ang, 'beta': b_ang, 'gamma': g_ang}

            self.cryst_ptgrp = 'D6h'
            self.burgers_mag = a
            self.eam_file = np.array(['fs', 'Mg.eam.fs'])

            b1x = a*np.array([1., 0., 0.])
            b1y = np.dot(vrrotvec2mat(np.array([0., 0., 1., g_ang])), b1x)
            b1z = c*np.array([0., 0., 1.])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))

            self.basis_atoms = np.array([[0., 0., 0.], [1./3, 2./3, 1./2]])
        # ------------------------------------------------------------------------------------------------------

        #### Tetragonal Lattices
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'tp_id':
            self.elem_type = 'tP_Id'
            self.pearson = 'tP'
            a = 1.0; b = a; CAratio = 1.2; c = a*CAratio
            a_ang = np.pi/2; b_ang = np.pi/2; g_ang = np.pi/2;
            self.lat_params = {'a': a, 'b': b, 'c': c, 'alpha': a_ang, 'beta': b_ang, 'gamma': g_ang}

            self.cryst_ptgrp = 'D4h'

            b1x = a*np.array([1., 0., 0.])
            b1y = a*np.array([0., 1., 0.])
            b1z = c*np.array([0., 0., 1.])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'tp_ca':
            if nargs < 2:
                raise Exception('Provide the c/a ratio')
            elif nargs == 2:
                    self.elem_type = 'tP_ca'
                    ca_ratio = args[1]
            self.pearson = 'tP'

            a = 1.0; b = a; c = a*ca_ratio
            a_ang = np.pi/2; b_ang = np.pi/2; g_ang = np.pi/2;
            self.lat_params = {'a': a, 'b': b, 'c': c, 'alpha': a_ang, 'beta': b_ang, 'gamma': g_ang}

            self.cryst_ptgrp = 'D4h'

            b1x = a*np.array([1., 0., 0.])
            b1y = a*np.array([0., 1., 0.])
            b1z = c*np.array([0., 0., 1.])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
        # ------------------------------------------------------------------------------------------------------

        # Rhombohedral lattices
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'hr_id':
            self.elem_type = 'hR_Id'
            self.pearson = 'hR'
            ## Lattice Parameters for Corundum
            a = 4.75
            c = 12.982
            ## Using the relation between hexagonal and rhombohedral lattices
            tau = 1/3 - (a**2/c**2)/2;
            cos_ang = tau/(1 - 2*tau);
            a_ang = np.arccos(cos_ang)
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': a_ang, 'beta': a_ang, 'gamma': a_ang}

            k1 = a*np.sqrt((2 - 2*cos_ang)/3.0)
            k2 = a*np.sqrt((1 + 2*cos_ang)/3.0)
            b1x = np.array([ np.sqrt(3)*k1/2, k1/2, k2])
            b1y = np.array([-np.sqrt(3)*k1/2, k1/2, k2])
            b1z = np.array([0, -k1, k2])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))

            self.cryst_ptgrp = 'D3d'
        # ------------------------------------------------------------------------------------------------------
        if elem_type.lower() == 'hr_ca':
            if nargs < 2:
                raise Exception('Provide the c/a ratio')
            elif nargs == 2:
                    self.elem_type = 'hR_ca'
                    ca_ratio = args[1]
            self.pearson = 'hR'

            a = 1
            c = a*ca_ratio
            ## Using the relation between hexagonal and rhombohedral lattices
            tau = 1./3. - (a**2/c**2)/2;
            cos_ang = tau/(1 - 2*tau);
            a_ang = np.arccos(cos_ang)
            self.lat_params = {'a': a, 'b': a, 'c': a, 'alpha': a_ang, 'beta': a_ang, 'gamma': a_ang}

            k1 = a*np.sqrt((2 - 2*cos_ang)/3.0)
            k2 = a*np.sqrt((1 + 2*cos_ang)/3.0)
            b1x = np.array([ np.sqrt(3)*k1/2, k1/2, k2])
            b1y = np.array([-np.sqrt(3)*k1/2, k1/2, k2])
            b1z = np.array([0, -k1, k2])
            self.l_p_po = np.column_stack((b1x, b1y, b1z))
            self.cryst_ptgrp = 'D3d'
        # ------------------------------------------------------------------------------------------------------

    def __str__(self):
        l1 = self
        str1 = 'Lattice:'
        str1 += 'Pearson Symbol: %s \n' %(l1.pearson)
        str1 += 'Lattice Parameters: \n a = %f \t b = %f \t c = %f \t alpha = %f \t beta = %f \t gamma = %f \n' \
            %(l1.lat_params['a'], l1.lat_params['b'], l1.lat_params['c'],
              l1.lat_params['alpha'], l1.lat_params['beta'], l1.lat_params['gamma'])
        str1 += 'Point Group: %s \n' %(l1.cryst_ptgrp)
        str1 += 'Primitive lattice (l_p_po): \n'
        str1 += str(l1.l_p_po)
        str1 += '\n'
        str1 += 'Crystal Point group: '
        str1 += self.cryst_ptgrp
        str1 += '\n'
        return str1
# ------------------------------------------------------------------------------------------------------
