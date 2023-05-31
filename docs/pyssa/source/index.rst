.. PyMOL-Plugin documentation master file, created by
   sphinx-quickstart on Sat Jun 18 20:20:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PySSA!
=========================================================
This plugin combines prediction engine from DeepMind's AlphaFold through the ColabFold with the molecular visualization software PyMOL.

.. figure:: /assets/pyssa_logo.png
   :height: 350px
   :width: 350px

   Logo of PySSA.

The prediction is done locally with the ColabFold which is based on the AlphaFold of DeepMind.
The structure alignment is done with the align command from PyMOL and
the images are also taken with PyMOL. All tasks were done with this plugin.

*Warren's Philosophy*

   We have chosen a free and open-source approach for PyMOL because we believe this strategy will have the greatest positive impact on humanity.
   Visualization is key part of understanding the nature of life at the molecular level, and powerful visualization tools 
   need to be universally available to all students and scientists if we are to make rapid progress in biomedical research.


.. _official PyMOL: https://pymolwiki.org/index.php/Main_Page
.. _modules: https://github.com/schrodinger/pymol-open-source/tree/master/modules

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   getting_started_project
   getting_started_prediction
   getting_started_analysis
   getting_started_image
   settings
   file_organization
