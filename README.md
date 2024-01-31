# PySSA
<!-- [![DOI](https://zenodo.org/badge/220207097.svg)](https://zenodo.org/badge/latestdoi/220207097) -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/urban233/PySSA/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/urban233/PySSA)](https://GitHub.com/urban233/PySSA/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/urban233/PySSA.svg)](https://GitHub.com/urban233/PySSA/graphs/contributors/)
[![GitHub release](https://img.shields.io/github/release/urban233/PySSA.svg)](https://github.com/urban233/PySSA/releases/)

[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/pylint-dev/pylint)


![PySSA_logo](https://github.com/urban233/PySSA/blob/main/assets/images/graphical_abstract_pyssa.png)
# Python plugin for protein Sequence to Structure Analysis (PySSA)

## Contents of this document
* [Description](#Description)
* [Contents of this repository](#Contents-of-this-repository)
  * [Sources](#Sources)
  * [Tutorial](#Tutorial)
  * [Images](#Images)
* [Installation](#Installation)
    * [Application](#Application)
    * [Source code](#Source-code)
* [Dependencies](#Dependencies)
* [Citation](#Citation)
* [References and useful links](#References-and-useful-links)
* [Acknowledgements](#Acknowledgements)

## Description
PySSA ('Python plugin for protein Sequence to Structure Analysis') is an open software project that
enables protein structure predictions and analysis on a local computer. The GUI-Plugin offers
extensive graphical functions for visualising and evaluating the quality of structure prediction.
In addition, it is possible to create high quality images of protein structures in a few clicks.
<!-- The scientific article describing PySSA can be found here: <a href="doi"> Title </a> -->

## Contents of this repository
### Sources
The <a href="https://github.com/urban233/PySSA/tree/main/pyssa">"pyssa"</a> subfolder contains all source code packages.

### Documentation
The <a href="https://github.com/urban233/PySSA/tree/main/docs">"docs"</a> folder contains the documentation as PDF document,
as well as all html pages of the internal help.
Use cases can be found here. <!-- TODO: Add reference to use cases here -->

### Assets
The <a href="https://github.com/urban233/PySSA/tree/main/assets">"assets"</a> folder consists of
the subfolder <a href="https://github.com/urban233/PySSA/tree/main/assets">"images"</a> which contains the PySSA logo.
<!-- If you are using PySSA for your own projects, feel free to acknowledge it by using the logo in your presentations etc. -->

## Installation
Until now, PySSA is tested and available for Windows 10 and 11.
### Windows
For a convenient and user-friendly installation, the <a href="https://github.com/urban233/PySSAInstaller">"PySSA-Installer"</a> is available
(click <a href="https://github.com/urban233/PySSAInstaller">here</a> to
automatically download the installer .exe of the latest version).
Download the installer
executable, start, and follow the instructions to install the PySSA-Installer.
After that, open the installer and install each component (WSL2, LocalColabfold and PySSA)
one after the other.
To start PySSA, double-click the created shortcut on the desktop.
To uninstall PySSA, open the PySSA-Installer and uninstall the components which
should get removed.

### Source code
It is also possible to install PySSA without the help of the PySSA-Installer.
In order to install PySSA a WSL2 environment needs to be setup and a WSL2
distro with a LocalColabfold installation is necessary to be able to run
structure predictions.
Mamba can be used to create a conda environment with the dependencies PyQt5, NumPy,
Pandas, Matplotlib, Biopython and the open-source version of PyMOL.

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

## Citation

## References and useful links

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

**PySSA was developed at:**
<br>
<br>Zielesny Research Group
<br>Westphalian University of Applied Sciences
<br>August-Schmidt-Ring 10
<br>D-45665 Recklinghausen Germany

**The PySSA project team would like to thank
the communities behind the open software libraries and especially Warren L. DeLano
for their amazing work.**
