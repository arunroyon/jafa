Usage
=====

Command line
------------

Single feature:

.. code-block:: bash

   jafa make-map \
     --cubes cube1.fits cube2.fits \
     --feature C60_18p9 \
     --outdir results/c60_18p9

Ratio map:

.. code-block:: bash

   jafa make-ratio \
     --cubes cube1.fits cube2.fits cube3.fits \
     --feature1 C60_7p0 \
     --feature2 C60_18p9 \
     --outdir results/c60_7_over_18p9

Custom feature:

.. code-block:: bash

   jafa make-map \
     --cubes cube1.fits \
     --wavelength 18.9 \
     --feature-window 18.75 19.15 \
     --anchors 18.0 18.2 19.2 19.5 \
     --outdir results/custom_18p9

Python API
----------

.. code-block:: python

   from jafa import FeatureMapper

   mapper = FeatureMapper(["cube_ch1.fits", "cube_ch4.fits"])
   mapper.make_feature_map("C60_18p9", outdir="results/c60")
   mapper.make_ratio_map("C60_7p0", "C60_18p9", outdir="results/ratio")

Lower-level functions are available for workflows that manage cubes in memory:

.. code-block:: python

   from jafa import load_cube, load_feature_database, make_feature_map

   cube = load_cube("cube_ch4.fits")
   features = load_feature_database()
   result = make_feature_map([cube], features["C60_18p9"], write_outputs=False)

