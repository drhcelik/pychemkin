.. PyChemkin documentation master file, created by
   sphinx-quickstart on Thu Aug 22 17:14:45 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyChemkin documentation
=======================

Introduction
------------
**PyChemkin** provides Pythonic access to `Ansys Chemkin`_. It facilitates programmatic customization of Chemkin simulation workflow within the Python ecosystem and permits access to Chemkin property and rate utilities as well as selected reactor models:

- Process Chemkin-compatible gas-phase mechanisms
- Evaluate species and mixture thermodynamic and transport properties
- Compute reaction rate of progress and species rate of production (ROP)
- Combine gas mixtures isothermally or adiabatically
- Find the equilibrium state of a gas mixture
- Run gas-phase batch reactor models

Prerequisites
-------------
- Ansys Chemkin release 2025R1 or newer and a valid license
- Python 3.9 or newer
- Numpy 1.14.0 or newer
- PyYAML 6.0 or newer
- Jupyter Notebook for running tutorials

.. note::
   - This project is under active development.
   - Ansys `Chemkin documents`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. _Ansys Chemkin: https://www.ansys.com/products/fluids/ansys-chemkin-pro/

.. _Chemkin Documents: https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/prod_page.html?pn=Chemkin&pid=ChemkinPro&lang=en/