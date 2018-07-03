.. _`output.rst`:

Output Files
============

.. _`output-dir`:

Directory Structure
-------------------
The following is the output directory structure,
where ``outputDirectory`` is defined by :cpp:var:`Options::outLoc`.
See below for descriptions of the files
in the `Data Directory`_, `Image Directory`_, and `ROOT Directory`_.

.. code-block:: text

  --\ `outputDirectory`
    |
    |--\ data
    |  |
    |  |-- contactAngle.txt
    |  |-- monoEdge.txt
    |  |-- bulkEdge.txt
    |  |-- dropletHeight.txt
    |  |
    |  |--\ circlePoints
    |     |-- 001.txt
    |     |-- 002.txt
    |     |-- ...
    |
    |--\ img
        |--\ dens
        |  |-- 001.png
        |  |-- 002.png
        |  |-- ...
        |
        |--\ droplet
        |  |-- 001.png
        |  |-- 002.png
        |  |-- ...
        |
        |--\ tanh
           |-- 001.png
           |-- 002.png
           |-- ...


.. _`output-data`:

Data Directory
--------------
The ``data`` directory is produced in :doc:`api/Writers`,
and contains all numerical output.

* The ``.txt`` files in the ``data`` directory
  are produced by :cpp:class:`ScalarWriter`.
  Each of these files are single-column data,
  one row per :cpp:class:`Frame`.

* The subdirectories of ``data``
  are produced by :cpp:class:`ArrayWriter`.
  Each subdirectory (there is currently only ``circlePoints``)
  contains one ``.txt`` file per :cpp:class:`Frame`,
  each of which contains multiple columns,
  with headers in the first line of each column.
  For the ``circlePoints`` files, these are the
  :math:`r` and :math:`z` coordinates of the boundary points.


.. _`output-image`:

Image Directory
---------------
The ``img`` directory is produced in :doc:`api/Visualization`,
and contains all figures that are generated each frame.
These figures are produced by default,
but can be disabled via :cpp:var:`Options::saveImages`.

* ``dens``: See :ref:`dens-figure`, produced by :cpp:class:`DensFigure`.
* ``tanh``: See :ref:`tanh-figure`, produced by :cpp:class:`TanhFigure`.
* ``droplet``: See :ref:`droplet-figure`, produced by :cpp:class:`DropletFigure`.


.. _`output-root`:

ROOT Directory
--------------

The ``root`` directory contains the same figures as the ``img``
directory, but saved as `ROOT Macros`_ rather than raster images,
so that they can be manipulated and inspected.
This directory is not produced by default,
but may be enabled via :cpp:var:`Options::saveROOT`.

.. _`ROOT Macros`: https://root.cern.ch/root/htmldoc/guides/primer/ROOTPrimer.html#root-macros
