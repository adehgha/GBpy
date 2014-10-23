<<<<<<< HEAD
GBpy
======================
is a complete python package for constructing the geometrical properties of
a Bicrystal. It includes all the necessary tools for constructing a simulation box
for grain boundary simulation calculations.

- `csl_fidner_smith`, a state machine
- `dsc_finder`, a state superclass
- `symmetry_operations`, a whitespace-sensitive version of `StateMachine`
etc.......

Classes:

- `LatType()`: Includes all the crystallographic data required for an element used by the code.


How To Use This Module
======================
(See the individual classes, methods, and attributes for details.)

1. Import it: ``import GBpy`` or ``from GBpy import ...``.
   You will also need to ``import numpy``.

2. Find the appropriate misorientation `R_G1toG2_G1` for the particular 
   Sigma value of your interest::

       class Sigma_Rots(Sigma):

   Within the class definition:

   a) Include find the numerator and denominator of the element of the rotation matrix 
      form the second and third element of the Sigma_Rots instance.

          R_N = Sigma_Rots['Sigma_Rots'][0, cnt][0][1]
          R_D = Sigma_Rots['Sigma_Rots'][0, cnt][0][2]

3. Find the CSL and DSC lattices from in G1 frame::

          [L_CSL_G1, L_DSC_G1] = find_csl_dsc(L_G1_GO1, R_G1toG2_G1)
       
4. Convert the obtained basis vectors to orthogonal basis frame.

.
.
.

- $Id: GBpy ???? 08-14-2014 ??:??:??Z milde $
- Authors: Srikanth Patala <spatala@ncsu.edu>
-          Arash D Banadaki <adehgha@ncsu.edu>
- Copyright: This module has been placed in the public domain.\n
=======
GBpy
======================
is a complete python package for constructing the geometrical properties of
a Bicrystal. It includes all the necessary tools for constructing a simulation box
for grain boundary simulation calculations.

- `csl_fidner_smith`, a state machine
- `dsc_finder`, a state superclass
- `symmetry_operations`, a whitespace-sensitive version of `StateMachine`
etc.......

Classes:

- `LatType()`: Includes all the crystallographic data required for an element used by the code.


How To Use This Module
======================
(See the individual classes, methods, and attributes for details.)

1. Import it: ``import GBpy`` or ``from GBpy import ...``.
   You will also need to ``import numpy``.

2. Find the appropriate misorientation `R_G1toG2_G1` for the particular 
   Sigma value of your interest::

       class Sigma_Rots(Sigma):

   Within the class definition:

   a) Include find the numerator and denominator of the element of the rotation matrix 
      form the second and third element of the Sigma_Rots instance.

          R_N = Sigma_Rots['Sigma_Rots'][0, cnt][0][1]
          R_D = Sigma_Rots['Sigma_Rots'][0, cnt][0][2]

3. Find the CSL and DSC lattices from in G1 frame::

          [L_CSL_G1, L_DSC_G1] = find_csl_dsc(L_G1_GO1, R_G1toG2_G1)
       
4. Convert the obtained basis vectors to orthogonal basis frame.

.
.
.

- $Id: GBpy ???? 08-14-2014 ??:??:??Z milde $
- Authors: Srikanth Patala <spatala@ncsu.edu>
-          Arash D Banadaki <adehgha@ncsu.edu>
- Copyright: This module has been placed in the public domain.\n

>>>>>>> e99ecbcc3f522d28a256d3103c68641b11244dc3
