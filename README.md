# PySSA
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12667158.svg)](https://zenodo.org/doi/10.5281/zenodo.12667157)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/urban233/PySSA/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/urban233/PySSA)](https://GitHub.com/urban233/PySSA/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/urban233/PySSA.svg)](https://GitHub.com/urban233/PySSA/graphs/contributors/)
[![GitHub release](https://img.shields.io/github/release/urban233/PySSA.svg)](https://github.com/urban233/PySSA/releases/)
[![Software Article - JChemInf](https://img.shields.io/badge/Software_Article_Preprint-ChemRxiv-blue)](https://doi.org/10.26434/chemrxiv-2024-srx5d)
[![Application Example - Use Cases](https://img.shields.io/badge/Application%20Example-Use%20Cases-blue)](https://github.com/urban233/PySSA/wiki/Use-Cases)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)


<div id='header' align='center'>
  <img src='assets/images/title_logo.png' width='450'/>
</div>

# Python rich client for visual protein Sequence to Structure Analysis (PySSA)

## <img src='https://github.com/primer/octicons/blob/main/icons/download-24.svg' width='32'/> [Quick Installation](https://github.com/urban233/PySSA/wiki/Installation-for-Windows-Operating-System)
## <img src='https://github.com/primer/octicons/blob/main/icons/video-24.svg' width='32'/> [Use Cases](https://github.com/urban233/PySSA/wiki/Use-Cases)
A "Use Case" refers to a specific scenario that demonstrates how the software can be used to accomplish a particular task or objective.
There are three interconnected use cases that progressively demonstrate the capabilities of PySSA: the first use case illustrates how to use PySSA to validate the ColabFold prediction method, the second involves comparing the wild-type protein to a mutated protein to determine if the amino acids crucial for receptor interaction remain in the same position, and the third showcases how to generate a high-quality, ray-traced image of a specific amino acid.

## Contents of this document
* [Description](#Description)
* [Contents of this repository](#Contents-of-this-repository)
  * [Sources](#Sources)
  * [Documentation](#Documentation)
  * [Assets](#Assets)
* [Installation](#Installation)
    * [Windows](#Windows)
      * [Installation for Windows OS](#installation-for-windows-os)
      * [Offline environment](#offline-environment)
      * [Online installation (For experts)](#online-installation-for-experts)
    * [Source code](#Source-code)
* [Dependencies](#Dependencies)
* [Architecture](#Architecture)
* [Use cases](#Use-cases)
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

Use cases describing specific workflows can be found [here](https://github.com/urban233/PySSA/wiki/Use-Cases). 

### Assets
The <a href="https://github.com/urban233/PySSA/tree/main/assets">"assets"</a> folder consists of
the subfolder <a href="https://github.com/urban233/PySSA/tree/main/assets/images">"images"</a> which contains the PySSA logo.
If you are using PySSA for your own projects, you are welcome to give credit to PySSA by using the logo in your presentations, etc.

## Installation
PySSA is tested and available for Windows 10 and 11.
### Windows
For a convenient and user-friendly installation, the <a href="https://github.com/urban233/ComponentInstaller">"PySSA Component Installer"</a> is available.

**Important:**
* WSL2 **cannot** be uninstalled, once it is installed! Windows will integrate the WSL2 as a system component.
* Be aware that the computer needs to be **restarted** after installing WSL2.

#### Installation for Windows OS
A step-by-step guide can be found in the Wiki (click [here](https://github.com/urban233/PySSA/wiki/Installation-for-Windows-Operating-System) to go to the guide).

#### Offline Environment
A short guide on how to install PySSA in an offline environment can be found [here](https://github.com/urban233/PySSA/wiki/Advanced-Installation#offline-environment).

#### Online Installation (For experts)
A short guide on how to install PySSA by downloading ColabFold and PySSA during installation (i.e. online) can be found [here](https://github.com/urban233/PySSA/wiki/Advanced-Installation#online-installation-for-experts).

The PySSA Component Installer [user guide](https://github.com/urban233/ComponentInstaller/blob/v1.0.1/deployment/inno_setup/PySSA-Component-Installer-User-Guide.pdf) is available in the installer's menu under _Help_.
If WSL2 or ColabFold installation fails, consult the user guide and read the troubleshooting section.

### Source code
This is a Python project based on a virtual environment. 
To modify the source code, download or clone the repository 
and open it in an IDE that supports virtual environments (e.g. PyCharm).
Finally, run `pip install -r requirements_dev.txt` to set up the virtual environment used for development.

The project supports using a setup.py file to create a package that works with the PySSA Component Installer.
To build the package run `python setup.py create_win_package`.

The setup.py also supports building the documentation using the command:
`python setup.py make_docs`.

## Dependencies
**Managed by PySSA-Installer:**
* Windows Subsystem for Linux 2
  * WSL2
  * License: Microsoft Software License Terms
* Colabfold
  * [LocalColabfold](https://github.com/YoshitakaMo/localcolabfold)
  * License: MIT License
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
<div id='header' align='center'>
  <img src='assets/images/architecture.png' width='650'/>
</div>

## Use cases
Three distinct use cases are presented to illustrate the functionality of PySSA. 
Each use case builds upon the preceding one, 
and it is therefore recommended to view all three in order to gain 
a comprehensive understanding of the software's capabilities.

The use cases can be found under this wiki page: [Use Cases](https://github.com/urban233/PySSA/wiki/Use-Cases)

## Citation
You can cite this software or this repository as it is defined in the CITATION.cff file.

## References and useful links
**ColabFold**
* [ColabFold GitHub repository](https://github.com/sokrypton/ColabFold)
* [Localcolabfold GitHub repository](https://github.com/YoshitakaMo/localcolabfold)
* [Mirdita, M., Schütze, K., Moriwaki, Y. et al. ColabFold: making protein folding accessible to all. Nat Methods 19, 679–682 (2022). https://doi.org/10.1038/s41592-022-01488-1](https://doi.org/10.1038/s41592-022-01488-1)

**PyMOL**
* [Open-source GitHub repository](https://github.com/schrodinger/pymol-open-source)
* [Open-source Windows Python wheelfiles](https://github.com/cgohlke/pymol-open-source-wheels)

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

