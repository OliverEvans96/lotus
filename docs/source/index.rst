.. _`index.rst`:

===================
Lotus Documentation
===================

Analyze molecular dynamics simulations of nanoscale wetting.

This is a C++ package which parses output from LAMMPS_,
and uses ROOT_ to determine geometrical quantities over time
such as droplet height, radius, and contact angle.

The source code is available on GitHub_.

.. _ROOT: https://root.cern.ch/

.. _LAMMPS: http://lammps.sandia.gov/

.. _GitHub: https://github.com/OliverEvans96/lotus

.. toctree::
   :caption: INSTALLATION
   :maxdepth: 2

   install

.. toctree::
   :caption: USAGE
   :maxdepth: 2

   overview
   output

.. toctree::
   :caption: API REFERENCE
   :maxdepth: 2
   :glob:

   api/*

.. toctree::
   :caption: MISCELLANEOUS
   :maxdepth: 2

   todo
   genindex
