# PySSA
<!-- [![DOI](https://zenodo.org/badge/220207097.svg)](https://zenodo.org/badge/latestdoi/220207097) -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/zielesny/PySSA/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/zielesny/PySSA)](https://GitHub.com/zielesny/PySSA/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/zielesny/PySSA.svg)](https://GitHub.com/zielesny/PySSA/graphs/contributors/)
[![GitHub release](https://img.shields.io/github/release/zielesny/PySSA.svg)](https://github.com/urban233/PySSA/releases/)

[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/pylint-dev/pylint)


![PySSA_logo](assets/images/title_logo.png)
# Python rich client for visual protein Sequence to Structure Analysis (PySSA)

## Contents of this document
* [Description](#Description)
* [Contents of this repository](#Contents-of-this-repository)
  * [Sources](#Sources)
  * [Tutorial](#Tutorial)
  * [Images](#Images)
* [Installation](#Installation)
    * [Windows](#Windows)
    * [Source code](#Source-code)
* [Dependencies](#Dependencies)
* [Citation](#Citation)
* [References and useful links](#References-and-useful-links)
* [Acknowledgements](#Acknowledgements)

## Description
PySSA ('Python rich client for visual protein Sequence to Structure Analysis') is an open software project that 
aims to combine PyMOL and ColabFold to enable the prediction and analysis of 
3D protein structures for the scientific end-user.
PySSA allows the creation of managed and shareable projects with defined workflows for 
the prediction and analysis of protein structures, 
which can be conveniently carried out by scientists without 
any special computer skills or programming knowledge on their local computers. 
PySSA can help make protein structure prediction accessible for research and 
development in protein chemistry and molecular biology, and for teaching. <br>
In addition to the prediction and analysis capabilities, PySSA has a more user-friendly interface 
for interacting with PyMOL, such as creating high-quality ray-tracing images in a few clicks.
<!-- The scientific article describing PySSA can be found here: <a href="doi"> Title </a> -->

## Contents of this repository
### Sources
There are five different Python packages that contain functionality for specific aspects of the architecture. 

- _pyssa_ 
  - The package contains Python modules, Qt .ui files and cascading stylesheets.
- _application_process_
  - The package contains a Python module for the _ApplicationProcessManager_ class.
- _auxiliary_pymol_ 
  - The package contains Python modules for the communication and integration of _Auxiliary PyMOL_.
- _pyssa_colabfold_
  - The package contains a modified version of the [_batch.py_](https://github.com/sokrypton/ColabFold/blob/main/colabfold/batch.py) module of the [ColabFold library](https://github.com/sokrypton/ColabFold). Moreover, the package contains the colabfold_run.py  and the service.py module for managing the ColabFold microservice.
- _pyssa_pymol_
  - The package contains Python modules for the communication and integration of _User PyMOL_.

### Documentation
The <a href="https://github.com/urban233/PySSA/tree/main/docs">"docs"</a> folder 
contains the end-user documentation in the form of markdown and HTML files. 
The subfolder <a href="https://github.com/urban233/PySSA/tree/main/docs/dev-notes">"dev-notes"</a>, contains development notes.

Use cases describing specific workflows can be found here (Available soon). 

### Assets
The <a href="https://github.com/urban233/PySSA/tree/main/assets">"assets"</a> folder consists of
the subfolder <a href="https://github.com/urban233/PySSA/tree/main/assets">"images"</a> which contains the PySSA logo.
If you are using PySSA for your own projects, you are welcome to give credit to PySSA by using the logo in your presentations, etc.

## Installation
PySSA is tested and available for Windows 10 and 11.
### Windows
For a convenient and user-friendly installation, 
the <a href="https://github.com/urban233/ComponentInstaller">"PySSA-Installer"</a> is available
(click <a href="https://github.com/urban233/ComponentInstaller">here</a> to
automatically download the setup.exe of the latest version).
Download the installer executable, start, and follow the instructions to install the PySSA-Installer.
After that, open the installer and install each component (WSL2, LocalColabfold and PySSA) one after the other.
To start PySSA, double-click the created shortcut on the desktop.
To uninstall PySSA, open the PySSA-Installer and uninstall the components which should get removed.

To install PySSA in an offline environment the offline setup needs to be downloaded (available [here]()).
This setup contains the PySSA-Installer together with all necessary files that are needed for installing ColabFold and PySSA without an internet connection. 
The WSL2 needs to be setup beforehand. 
To update ColabFold or PySSA it is necessary to download a new offline setup.

### Source code
PySSA should be build using the provided setup.py script with different 
build tasks.

To build the package that is used by the PySSA-Installer to install PySSA, run
`python setup.py create_win_package`.

To generate only the user docs, use `python setup.py make_docs`.

## Dependencies
**Managed by PySSA-Installer:**
* Windows Subsystem for Linux 2
  * WSL2
  * License: Microsoft Software License Terms
* Colabfold
  * [LocalColabfold](https://github.com/YoshitakaMo/localcolabfold)
  * License: MIT License
* [Mamba](https://github.com/mamba-org/mamba)
  * License: BSD 3-Clause "New" or "Revised" License
* [PyQt5](https://riverbankcomputing.com/software/pyqt/intro)
  * License: GNU General Public License (GPL)
* [NumPy](https://numpy.org/)
  * License: BSD 3-Clause "New" or "Revised" License
* [Pandas](https://github.com/pandas-dev/pandas)
  * License: BSD 3-Clause "New" or "Revised" License
* [Matplotlib](https://matplotlib.org/)
  * License: Python Software Foundation License (PSF)
* [Biopython](https://biopython.org/)
  * License: BSD 3-Clause License
* [PyMOL Open-Source](https://github.com/schrodinger/pymol-open-source)
  * License: [BSD-like license](https://github.com/schrodinger/pymol-open-source/blob/master/LICENSE)
* SQLite3
  * License: [Public Domain](https://www.sqlite.org/copyright.html)

## Architecture
![PySSA_logo](assets/images/architecture.png)

## Citation
TODO: Add citation here.

## References and useful links
TODO: Add any references or useful links here.

## Acknowledgements
**Developers:**
* Martin Urban
* Hannah Kullik

**End-user testers**
* Jonas Schaub
* Achim Zielesny

**Logo:**
* Martin Urban
* Hannah Kullik

**Initialization, conceptualization, and supervision:**
* Achim Zielesny and Angelika Loidl-Stahlhofen

**The PySSA project team would like to thank
the communities behind the open software libraries and especially Warren L. DeLano
for their amazing work.**

<!--
**PySSA was developed at:**
<br>
<br>Zielesny Research Group
<br>Westphalian University of Applied Sciences
<br>August-Schmidt-Ring 10
<br>D-45665 Recklinghausen Germany
--!>

