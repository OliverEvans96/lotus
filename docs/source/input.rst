.. _`input.rst`:

Lotus Input File
================

The ``lotus`` binary takes a single argument:
the path to a lotus input file.

The input file format is a basic subset of YAML_
which has the following restrictions:

- Root node must be a mapping
- Keys must scalars or lists (no nested mappings).

The reason for these restrictions is that lotus
is required to be C++98 compliant to run on our
cluster which hasn't been updated since 2005.
Therefore, I had to manually write a parser
using libyaml_ rather than use the much friendlier
`yaml-cpp`_ which requires C++11.

A more complex parser could of course be written with
further effort, but so far this has been sufficient.

The following is an example YAML script from the test suite
(which won't run without a subset of our simulation data
which we can't share until we publish our results).

.. literalinclude:: ../../test/data/test_config.yaml
   :language: yaml
   :caption: ``test/data/test_config.yaml``

All public member variables of :cpp:class:`Options`
are valid YAML keys for this input file,
some of which are required while others are optional.

.. seealso::
   All options and their acceptable and default values
   are described in depth in the :cpp:class:`Options` documentation.

.. todo:: mention LAMMPS dumpfile/datafile requirements somewhere.

.. _YAML: http://yaml.org/
.. _libyaml: https://pyyaml.org/wiki/LibYAML
.. _`yaml-cpp`: https://github.com/jbeder/yaml-cpp/

