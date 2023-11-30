# PySSA

<img align="center" src="https://github.com/urban233/tmpPySSA/blob/main/pyssa_github.gif" alt="pyssa_intro" />

This is the temporary repo for the PySSA-PyMOL-Plugin.

### Code checked with: 

[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/pylint-dev/pylint)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

<!-- [![DOI](https://zenodo.org/badge/220207097.svg)](https://zenodo.org/badge/latestdoi/220207097) -->
DOI: 
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/urban233/PySSA/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/FelixBaensch/MORTAR.svg)](https://GitHub.com/urban233/PySSA/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/FelixBaensch/MORTAR.svg)](https://GitHub.com/urban233/PySSA/graphs/contributors/)
[![GitHub release](https://img.shields.io/github/release/FelixBaensch/MORTAR.svg)](https://github.com/urban233/PySSA/releases/)

![PySSA_logo](https://github.com/urban233/PySSA/blob/main/graphical_abstract_pyssa.tiff)
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
PySSA ('Python plugin for protein Sequence to Structure Analysis') is an open software project that supports workflows of
molecular fragmentation and substructure analysis. The Java/JavaFX rich-client application offers extensive graphical 
functions for visualising the fragmentation results of individual compounds or entire compound sets. With several 
views and analysis functions, MORTAR supports the interpretation of fragmentation results. In addition to three 
currently integrated methods for fragmentation and substructure analysis - 
<a href="https://github.com/zielesny/ErtlFunctionalGroupsFinder">ErtlFunctionalGroupsFinder</a>, 
<a href="https://github.com/JonasSchaub/SugarRemoval">Sugar Removal Utility</a>, 
and <a href="https://github.com/Julian-Z98/ScaffoldGenerator">Scaffold Generator</a> - MORTAR allows straightforward integration of 
additional fragmentation algorithms with automatic generation of settings 
menus. All cheminformatics functionalities are implemented based on the <a href="https://github.com/cdk/cdk">Chemistry Development Kit (CDK)</a>.<br>
The scientific article describing MORTAR can be found here: 
<a href="https://doi.org/10.1186/s13321-022-00674-9"> MORTAR: a rich client application for in silico molecule fragmentation (Baensch et al. 2023) </a>

## Contents of this repository
### Sources
The <a href="https://github.com/urban233/PySSA/tree/main/pyssa">"pyssa"</a> subfolder contains all source code packages.

### Tutorial
The <a href="https://github.com/FelixBaensch/MORTAR/tree/master/Tutorial">"Tutorial" folder</a> contains a PDF document 
with a detailed tutorial on how to install and use MORTAR, together with
a test data set. 

### Images
The <a href="https://github.com/FelixBaensch/MORTAR/tree/master/Images">"Images" folder</a> contains the MORTAR logo and 
icon as image files that were created by <a href="https://github.com/Kohulan">Kohulan Rajan</a>.
If you are using MORTAR for your own projects, feel free to acknowledge it by using the logo in your presentations etc.

## Installation
### Application
Pre-compiled and executable MORTAR distributions can be found attached to the 
<a href="https://github.com/FelixBaensch/MORTAR/releases">marked releases</a>.

<p>
<b>Windows:</b> A convenient Windows OS installer executable for MORTAR is available 
(click <a href="https://github.com/FelixBaensch/MORTAR/releases/download/v1.1.1.0/MORTAR_v1.1.1.0_WINx64_installer.exe">here</a> to 
automatically download the installer .exe of the latest version). Download the installer 
executable, start, and follow the instructions to install MORTAR. Note that the installation includes a full 
Java Runtime Environment (JRE). After installation, create a shortcut to an appropriate MORTAR start batch file on your 
Windows desktop. E.g. for MORTAR to use up to 4 gigabyte of RAM, copy a shortcut to batch file "MORTAR.bat" which is 
located in the MORTAR program folder (default "C:\Program Files\MORTAR\MORTARv1.1.1.0\bin" or the path specified at 
installation). To start MORTAR, double click the created shortcut. MORTAR can be uninstalled by the provided 
Uninstall.exe executable in the MORTAR program folder or standard Windows functions.
<br>
As an alternative to "MORTAR.bat", there is also the "MORTAR_20GB.bat" batch file available that allocates up to 
20 GB of RAM for MORTAR. If you want to configure your own heap space settings, open one of the provided batch files 
and adjust the line
</p>

<p><code>set DEFAULT_JVM_OPTS="-Xms4g" "-Xmx4g"</code></p>

with your chosen initially allocated memory (-Xms) and maximum value (-Xmx) accordingly.<br>

<p>Should this installation or the execution of the batch files not work for you, try the guidelines for Linux and MaxOS 
described below. As an alternative way, they should also work on Windows.
</p>

<p><b>Linux and macOS:</b> Every release has the executable Java ARchive (JAR) "MORTAR-fat-1.1.1.0.jar"
attached, which contains the packaged MORTAR code together with all dependencies 
(click <a href="https://github.com/FelixBaensch/MORTAR/releases/download/v1.1.1.0/MORTAR-fat-1.1.1.0.jar">here</a> to 
automatically download the JAR of the latest version). 
To run MORTAR (with up to 4 GB of RAM available, e.g.), 
execute the JAR from the command-line using</p>

<p><code>java -jar -Xms512m -Xmx4g [path to]MORTAR-fat-1.1.1.0.jar</code></p>

A JDK or JRE of version 17.0.4 or higher needs to be installed on your system and linked to the "java" command. 
Otherwise, replace "java" with the path to the java command of your JDK or JRE.<br>

<p>Please note that MORTAR only supports x64 (on all three platforms) and AArch64 (on macOS and Linux) architectures in general. 
For the latter, a special "fat JAR" named "MORTAR-fat-aarch64-1.1.1.0.jar" is available from the distributions attached to the releases and must be used 
(click <a href="https://github.com/FelixBaensch/MORTAR/releases/download/v1.1.1.0/MORTAR-fat-aarch64-1.1.1.0.jar">here</a> to 
automatically download the AArch64 JAR of the latest version).</p>
Also note that using the Windows Subsystem for Linux (WSL) is not recommended, since a lot of additional configurations 
have to be made there to run Java GUI applications.

### Source code
This is a Gradle project. In order to use the source code for your own software or do your own MORTAR build, download or 
clone the repository and open it in a Gradle-supporting IDE (e.g. IntelliJ) as a Gradle project and execute the 
build.gradle file. Gradle will then take care of installing all dependencies. A Java Development Kit (JDK) of version 17.0.4 
or higher must also be pre-installed and set as project JDK / project compiler.
The Gradle build process is configured to include a specific Java Runtime Environment (JRE) in the "install" folder.
For this to work, you need to create an "AdoptOpenJDK\jdk-17.0.4_8_jre\" folder and put the JRE with the specified version 
into it (i.e. sub-folders of "AdoptOpenJDK\jdk-17.0.4_8_jre\" need to be "bin", "conf", "legal", "lib", etc.).

## Dependencies
**Needs to be pre-installed:**
* Java Development Kit (JDK) version 17.0.4 or higher
    * [Adoptium Open JDK](https://adoptium.net) (as one possible source of the JDK)
* Gradle version 7.3
    * [Gradle Build Tool](https://gradle.org)

**Managed by Gradle:**
* JavaFX version 17.0.2
  * [Open JavaFX](https://openjfx.io)
  * GNU General Public License (GPL) Version 2
* Chemistry Development Kit (CDK) version 2.8
    * [Chemistry Development Kit on GitHub](https://cdk.github.io/)
    * License: GNU Lesser General Public License 2.1
* JUnit version 5.9.3
    * [JUnit 5](https://junit.org/junit5/)
    * License: Eclipse Public License 2.0
* LibrePDF OpenPDF version 1.3.26
  * [OpenPDF GitHub repository](https://github.com/LibrePDF/OpenPDF)
  * License: GNU Lesser General Public License 2.1
* Spotless version 6.19
  * [Spotless GitHub repository](https://github.com/diffplug/spotless)
  * License: Apache-2.0 license

## Citation
You can cite this software or this repository as it is defined in the CITATION.cff file. Also, please cite our scientific 
article (<a href= "https://doi.org/10.1186/s13321-022-00674-9"> MORTAR: a rich client application for in silico molecule 
fragmentation (Baensch et al. 2023)).

## References and useful links
* None

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

**The PySSA project team would like to thank the communities behind the open software libraries the application employs 
for their amazing work.**
