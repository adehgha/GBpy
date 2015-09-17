**GBpy** is a python package for finding the geometrical properties of a Bicrystal. It includes all the necessary tools for constructing a simulation box for grain boundary simulation.
                
.. image:: https://raw.githubusercontent.com/adehgha/GBpy/master/GBpy/docs/images/pic.png
        
* GBpy on `GitHub <https://github.com/adehgha/GBpy>`__
* GBpy `Documentation <https://cdn.rawgit.com/adehgha/GBpy/master/GBpy/docs/_build/html/index.html>`__
* GBpy `Tutorials <https://cdn.rawgit.com/adehgha/GBpy/master/GBpy/docs/tutorials/index.html>`__

.. image:: https://img.shields.io/pypi/v/GBpy.svg
    :target: https://pypi.python.org/pypi/GBpy
    :alt: Latest PyPI version

.. image:: https://img.shields.io/pypi/dm/GBpy.svg
    :target: https://pypi.python.org/pypi/GBpy
    :alt: Number of PyPI downloads

     
Functions:
==========
        
* ``GBpy.find_csl_dsc``, collection of functions for computing the CSL and the DSC lattices of a general bicrystal (general lattices), if the transformation **T** is given.
* ``GBpy.generate_symm_ops``, a function for generating various point group symmetry operations.
* ``GBpy.bp_basis``, collection of functions for calculating the basis vectors of a two-dimensional lattice of an interface.
* ``GBpy.quaternion``, collection of functions for quaternion operations.
* ``GBpy.misorient_fz``, function for finding the unique disorientations in fundamental zone of various crystal point groups.
* ``GBpy.integer_manipulations``, collection of many useful ineteger manipulation functions.
                
and many other useful tools. Please refer to the `documentation <https://cdn.rawgit.com/adehgha/GBpy/master/GBpy/docs/_build/html/index.html>`__ and `tutorials <https://cdn.rawgit.com/adehgha/GBpy/master/GBpy/docs/tutorials/index.html>`__ for detailed description and utility of functions.
                
Classes:
========
                
- ``lattice``: Includes all the crystallographic data required for an element used by the code.
- ``quaternion``: Quaternion construction and operations.
        
        
How to Use This Package:
========================
1.  **To install the stable version of GBpy:**      
    
    .. code-block:: console
                
        $ pip install GBpy
                                       
                
    *To install the development version of GBpy* Clone the repository:   
        
    .. code-block:: console
                
        $ git clone https://github.com/adehgha/GBpy.git   
             
    and run the setup script.                	

    .. code-block:: console     
           
        $ python setup.py install
                   
2.  **Import the package:** 
                
    .. code-block:: pycon
                
        >>> import GBpy
                          
3.  **Call the function by using:**
                
    .. code-block:: pycon
                
        >>> GBpy.<name_of_the_function>
                	
    * for example to find the 2D basis vectors of a plane with Miller indices of (h,k,l):
                
    .. code-block:: pycon
                
        >>> GBpy.bp_basis.bp_basis([h,k,l])
                
4.  **You can also use the tools provided in this package individually by importing the functions separately.** For example use :``from GBpy import <name_of_the_function> as <a_name>``.


                
Consult the `documentation <https://cdn.rawgit.com/adehgha/GBpy/master/GBpy/docs/_build/html/index.html>`__ for further details.
        
        
Prerequisites:
==============
                
1. install ``numpy`` from `here. <http://www.numpy.org/>`__
                
2. install ``scipy`` from `here. <http://www.scipy.org/>`__
                
3. install ``setuptools`` from `here. <https://pypi.python.org/pypi/setuptools>`__
                
Cite GBpy:
========================

"An Efficient Algorithm for Computing the Primitive Bases of a General Lattice Plane", A. Banadaki, S. Patala, *Journal of Applied Crystallography*, v. 48, 2015, `doi:10.1107/S1600576715004446. <http://scripts.iucr.org/cgi-bin/paper?S1600576715004446>`__"

                
Credits:
========
GBpy is written by:
                
* `Srikanth Patala <spatala@ncsu.edu>`__
* `Arash Dehghan Banadaki <adehgha@ncsu.edu>`__
* `Patala Research Group <http://research.mse.ncsu.edu/patala/>`__.
        
Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
