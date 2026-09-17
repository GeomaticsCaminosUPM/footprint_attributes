footprint_attributes
=====================

Automated computation of seismic behaviour modifiers from 2-D building
footprint polygons. Companion package to:

    Ureña-Pliego, M., Rodríguez-Saiz, J., Núñez-Álvarez, G.,
    Marchamalo-Sacristán, M., González-Rodrigo, B. (2026). *A methodology
    for the automated estimation of footprint-derived seismic behaviour
    modifiers in building exposure assessment.* Advanced Modeling and
    Simulation in Engineering Sciences, 13:3.
    `doi.org/10.1186/s40323-026-00323-y
    <https://link.springer.com/content/pdf/10.1186/s40323-026-00323-y.pdf>`_

Developed by the `Advanced Geomatics group (AGA)
<https://blogs.upm.es/aga/en/>`_ at the Universidad Politécnica de Madrid.
See :doc:`citation` for the full author list, ORCIDs, and funding.

Source code: `github.com/GeomaticsCaminosUPM/footprint_attributes
<https://github.com/GeomaticsCaminosUPM/footprint_attributes>`_

.. image:: ../figures/graphical_abstract.jpg
   :width: 80%
   :align: center

Interactive map
---------------

Every value below is computed by this package, from the raw footprint
geometry, for three real pilot regions (Guatemala Zona 10, San José Mata
Redonda, Santo Domingo Ensanche Quisquella). It starts in an auto-playing
tour -- relative position, then the shape index, then each of the 15 raw
shape metrics in turn -- and orbits the camera in 3D; interact with it
(drag/zoom/click) to take over, or see :doc:`examples` for the full
walkthrough, checkboxes, and screenshots.

.. raw:: html

   <iframe src="_static/maps/index.html" width="100%" height="640" style="border:1px solid #444;" loading="lazy"></iframe>

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   api/index
   formulas
   examples
   citation
   license

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
