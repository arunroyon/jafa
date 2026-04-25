Installation
============

Core install
------------

The core package keeps import-time dependencies light:

.. code-block:: bash

   python -m pip install jafa

JAFA — JWST Aromatic Feature Analyzer uses ``jafa`` as the Python package
name and primary command-line executable.

Core runtime dependencies are ``numpy``, ``scipy``, ``astropy``, and
``PyYAML``.

JWST science workflow
---------------------

For JWST cube loading, morphology continua, reprojection, and plotting:

.. code-block:: bash

   python -m pip install "jafa[jwst,plot]"

For local development from this repository:

.. code-block:: bash

   python -m pip install -e ".[jwst,plot,test,docs,dev]"

Conda environment
-----------------

Compiled astronomy dependencies are often easiest through conda-forge:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate jafa
   python -m pip install -e .

The former ``jwst_feature_mapper`` import path and ``jwst-feature-map``
command are available as legacy aliases during the rename transition.
