#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
import subprocess
from pyssa.util import constants

if __name__ == "__main__":
    print(constants.SETTINGS_DIR)
    tmp = constants.SETTINGS_DIR.replace("\\", "/")
    tmp2 = tmp.replace(":", "")
    tmp3 = tmp2.replace("C", "c")
    SETTINGS_DIR_UNIX_NOTATION = f"/mnt/{tmp3}"
    print(SETTINGS_DIR_UNIX_NOTATION)
    #subprocess.run(["podman", "machine", "start"])
    #subprocess.run(["podman", "run", "-itd", "--name", "localcolabfold-container",
    #                "localhost/localcolabfold-ubuntu2204:1.5.1.2", "exit"])

    #subprocess.run(["podman", "stop", "localcolabfold-container"])
    #subprocess.run(["podman", "container", "start", "localcolabfold-container"])
    #subprocess.run(["podman", "container", "list"])
    #subprocess.run(["podman", "container", "exec", "localcolabfold-container",
    #                "/home/ubuntu_colabfold/localcolabfold/colabfold-conda/bin/colabfold_batch", "--help"])
    # subprocess.run(["podman", "container", "exec", "busy_hugle",
    #                 "/home/ubuntu_colabfold/localcolabfold/colabfold-conda/bin/colabfold_batch",
    #                 "/home/ubuntu_colabfold/test/fasta", "/home/ubuntu_colabfold/test/pdb",
    #                 "--amber", "--templates"])
    #subprocess.run(["podman", "machine", "stop"])
