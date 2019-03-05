.. Linearized Free 3D Surface Wave Resistance Code documentation master file, created by
   sphinx-quickstart on Sat Jan 26 16:27:14 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to documentation for the 3D Linearized Free Surface Wave Resistance Code
===========================================================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Features
--------

- Linearized free surface boundary condition
- Fortran 90 implementation


Build
------------

Build the executable:

- make file set up to use gfortran
- assumes you have make tools.
- might work with g95 but not tested
- cd into the main directory for this project and run:

$    make


Run
------------

Run the included wigley hull example with

$  ./p3 fifi.dat test1.out .2 1

- ./p3 runs the executable
- fifi.dat selects the included panel input file of the wigely hull
- test1.out sets an output file
- .2 sets the Froude number to 0.2
-  1 selects the SimQuit matrix inversion routine

Results
------------
- Hydrodynamic results are viewable in .vtp files
- VTK Paraview is the recommended viewing engine.

Notes
------------

This version tested on:

- OSX
- Linux
