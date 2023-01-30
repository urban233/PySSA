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

import os
import json
import pathlib


class Protein:
    """This class stores one protein in a PyMOL compatible form

    Var:
        selection:
            a pymol selection string which needs to conform with the selection algebra
            from pymol
    """
    selection: str = None

    def __init__(self, molecule_object: str, import_data_dir: str = None,
                 export_data_dir: str = None):
        """Constructor.

        Args:
            molecule_object (str):
                name of the reference Protein in the pymol session
            import_data_dir (str):
                directory where the pdb files of both model and
                reference are stored
            export_data_dir (str, optional):
                directory where all results related to the Protein
                will be stored.
                All subdirectories like ``images``, ``alignment_files``
                and ``distances`` will be created automatically.
                The export_data_dir will then function as parent directory.

        Raises:
            NotADirectoryError: If directory not found.
            FileNotFoundError: If file not found.
        """
        self.molecule_object = molecule_object
        self.import_data_dir = import_data_dir
        self.export_data_dir = export_data_dir

        # argument test
        if import_data_dir is not None:
            if not os.path.exists(f"{import_data_dir}"):
                raise NotADirectoryError(f"The path {import_data_dir} was not "
                                         f"found.")
        if export_data_dir is not None:
            if not os.path.exists(f"{export_data_dir}"):
                raise NotADirectoryError(f"The path {export_data_dir} was not "
                                         f"found.")
        if import_data_dir is not None:
            if not os.path.exists(f"{import_data_dir}/{molecule_object}.pdb"):
                raise FileNotFoundError(f"The pdb file {molecule_object} was "
                                        f"not found under {import_data_dir}.")

    def set_selection(self, selection: str) -> None:
        """This function sets a selection for the Protein object.

        Args:
            selection (str):
                A pymol conform selection as a string.

        Example:
            This is a pymol conform selection::

                "pymol_6omn////CA"
        """
        self.selection = selection

def create_generic_dictionary_from_directory(path):
    i = 0
    tmp_generic_dict = {}
    for item in os.listdir(path):
        tmp_generic_dict.update({f"token_{i}": item})
        i += 1
    return tmp_generic_dict


if __name__ == '__main__':
    print(create_generic_dictionary_from_directory(pathlib.Path("C:\\Users\\martin\\github_repos\\tmpPySSA\\pyssa\\util")))
