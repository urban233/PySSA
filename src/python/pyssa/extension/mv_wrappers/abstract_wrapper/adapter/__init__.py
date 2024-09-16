#
# PyDD - Python rich client for Drug Discovery
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
This module contains a description of the package "adapter".

The package contains the logic to connect PyDD with
a molecular visualization tool like ChimeraX or PyMOL.

In the case of ChimeraX, a bundle is needed to connect PyDD with ChimeraX.
The bundle should be structured as a PyQt gui tool.
The contents of this package will be used as bundle and
therefore be copied into the correct ChimeraX directory, if installed
with Toolshed.

In the case of PyMOL, a plugin is needed to connect PyDD with PyMOL.
The plugin should be structured as a PyQt gui plugin.
The contents of this package will be used as plugin and
therefore need to be copied into the plugin directory of PyMOL.
"""
