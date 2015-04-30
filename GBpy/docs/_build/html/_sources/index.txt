
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GBpy's documentation!
================================

.. toctree::
   :maxdepth: 2
   
   docs/find_csl_dsc
   Geometry Tools
   Integer Manipulations
   CSL Utility Function
   CSL/DSC Computation
   Boundary Plane Basis
   Misorientation Fundamental Zones
   CSL/DSC Computation
   Boundary Plane Basis
   Misorientation Fundamental Zones
   Generate Symmetry Operators
   Tools
   ...
   
Tutorials
=================  
For a complete list of tutorials please click `here. <../../tutorials/index.html>`__

Lattice Class
=================
.. automodule:: lattice

.. autoclass:: Lattice
    :members:


Geometry Tools
=================
.. automodule:: quaternion

.. autoclass:: quaternion
    :members:
.. autofunction:: sph2vec


Integer Manipulations
=====================
.. automodule:: integer_manipulations
.. autofunction:: gcd_array
.. autofunction:: lcm_vec
.. autofunction:: lcm_array
.. autofunction:: int_check
.. autofunction:: rat
.. autofunction:: int_finder
.. autofunction:: int_mult


CSL Utility Function
=====================
.. automodule:: csl_utility_functions

.. autofunction:: csl_rotations
.. autofunction:: proper_ptgrp
.. autofunction:: largest_odd_factor
.. autofunction:: compute_inp_params
.. autofunction:: mesh_muvw
.. autofunction:: mesh_muvw_fz
.. autofunction:: check_fsig_int
.. autofunction:: eliminate_idrots
.. autofunction:: sigtype_muvw
.. autofunction:: eliminate_mults
.. autofunction:: check_sigma
.. autofunction:: gcd1d_arr
.. autofunction:: compute_tmat
.. autofunction:: disorient_sigmarots
.. autofunction:: check_sigma_rots

CSL/DSC Computation
===================
.. automodule:: find_csl_dsc
.. autofunction:: find_csl_dsc
.. autofunction:: sigma_calc
.. autofunction:: reciprocal_mat
.. autofunction:: csl_elem_div_thm_l1
.. autofunction:: csl_elem_div_thm_l2
.. autofunction:: csl_finder_smith
.. autofunction:: check_csl_finder_smith
.. autofunction:: dsc_finder
.. autofunction:: check_dsc_finder

Boundary Plane Basis
=======================
.. automodule:: bp_basis
.. autofunction:: check_int_mats
.. autofunction:: check_2d_csl
.. autofunction:: lbi_dioph_soln
.. autofunction:: compute_basis_vec
.. autofunction:: bp_basis
.. autofunction:: pl_density
.. autofunction:: csl_finder_2d
.. autofunction:: gb_2d_csl
.. autofunction:: bicryst_planar_den


Misorientation Fundamental Zones
=======================================
.. automodule:: misorient_fz

.. autofunction:: misorient_fz
.. autofunction:: check_cond


Generate Symmetry Operators
======================================
.. automodule:: generate_symm_ops
.. autofunction:: generate_symm_mats
.. autofunction:: generate_symm_quats
.. autofunction:: save_symm_pkl

Tools
======================================
.. automodule:: tools
.. autofunction:: vrrotmat2vec
.. autofunction:: vrrotvec2mat
.. autofunction:: unique_rows_tol
.. autofunction:: axang2quat
.. autofunction:: mat2quat
.. autofunction:: quat2mat
.. autofunction:: lll_reduction
.. autofunction:: eq
.. autofunction:: message_display
.. autofunction:: ehermite
.. autofunction:: smith_nf
.. autofunction:: extgcd
.. autoclass:: Col
    :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

