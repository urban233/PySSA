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
#
import os
import pathlib
from PyQt5.QtCore import QObject, pyqtSignal, QThread

from pyssa.internal.data_structures import protein, protein_pair
from pyssa.io_pyssa import filesystem_io
from pyssa.io_pyssa.xml_pyssa import element_names, attribute_names
from pyssa.util import tools, constants
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project


def setup_worker_for_work(tmp_thread, tmp_worker, return_value_func):
    tmp_worker.moveToThread(tmp_thread)
    tmp_thread.started.connect(tmp_worker.run)
    tmp_worker.finished.connect(tmp_thread.quit)
    tmp_worker.finished.connect(tmp_worker.deleteLater)
    tmp_thread.finished.connect(tmp_thread.deleteLater)
    tmp_worker.return_value.connect(return_value_func)
    return tmp_thread


class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(tuple)
    workspace_path: str

    def __init__(self, workspace):
        super().__init__()
        self.workspace_path = workspace

    def run(self):
        protein_dict, protein_names = tools.scan_workspace_for_non_duplicate_proteins(pathlib.Path(self.workspace_path))
        self.return_value.emit((protein_dict, protein_names))
        self.finished.emit()


class CreateUseProjectWorker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(list)
    workspace_path: str
    proteins_to_copy: list

    def __init__(self, workspace, proteins_to_copy):
        super().__init__()
        self.workspace_path = workspace
        self.proteins_to_copy = proteins_to_copy

    def run(self):
        proteins_for_new_project = []
        protein_infos = (tools.scan_workspace_for_non_duplicate_proteins(pathlib.Path(self.workspace_path)))[0]
        for tmp_protein in self.proteins_to_copy:
            for tmp_protein_info in protein_infos:
                if tmp_protein_info.name == tmp_protein:
                    """Var: project_proteins is a list which contains all proteins from a single project"""
                    xml_deserializer = filesystem_io.XmlDeserializer(
                        pathlib.Path(f"{self.workspace_path}/{tmp_protein_info.project_name}.xml"))
                    for xml_protein in xml_deserializer.xml_root.iter(element_names.PROTEIN):
                        if xml_protein.attrib[attribute_names.ID] == tmp_protein_info.id:
                            basic_information = xml_protein.attrib
                            pdb_lines = []
                            session_data_base64 = ""
                            for tmp_data in xml_protein:
                                if tmp_data.tag == "pdb_data":
                                    for tmp_atom in tmp_data.findall("atom"):
                                        pdb_lines.append(tmp_atom.text)
                                elif tmp_data.tag == "session_data":
                                    session_data_base64 = tmp_data.attrib[attribute_names.PROTEIN_SESSION]
                                else:
                                    raise ValueError
                            tmp_protein_obj = protein.Protein(
                                molecule_object=basic_information[attribute_names.PROTEIN_MOLECULE_OBJECT],
                                pdb_xml_string=xml_protein)
                            tmp_protein_obj.set_all_attributes(basic_information, pdb_lines, session_data_base64)
            proteins_for_new_project.append(tmp_protein_obj)

        self.return_value.emit(proteins_for_new_project)
        self.finished.emit()


class LoadResultsWorker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(str)
    protein_pair_of_results: protein_pair.ProteinPair
    app_project_xml_filepath: str
    image_type: str

    def __init__(self, protein_pair_of_results, app_project_xml_filepath):
        super().__init__()
        self.protein_pair_of_results = protein_pair_of_results
        self.image_type = constants.IMAGES_NONE
        self.app_project_xml_filepath = app_project_xml_filepath

    def run(self):
        filesystem_io.XmlDeserializer(self.app_project_xml_filepath).deserialize_analysis_images(
            self.protein_pair_of_results.name, self.protein_pair_of_results.distance_analysis.analysis_results)
        if len(self.protein_pair_of_results.distance_analysis.analysis_results.structure_aln_image) != 0 and len(
                self.protein_pair_of_results.distance_analysis.analysis_results.interesting_regions_images) != 0:
            # if both image types were made during analysis
            self.protein_pair_of_results.distance_analysis.analysis_results.create_image_png_files_from_base64()
            self.image_type = constants.IMAGES_ALL
        elif len(self.protein_pair_of_results.distance_analysis.analysis_results.structure_aln_image) != 0 and len(
                self.protein_pair_of_results.distance_analysis.analysis_results.interesting_regions_images) == 0:
            # only struct align image were made
            self.protein_pair_of_results.distance_analysis.analysis_results.create_image_png_files_from_base64()
            self.image_type = constants.IMAGES_STRUCT_ALN_ONLY
        else:
            # no images were made
            self.image_type = constants.IMAGES_NONE

        self.return_value.emit(self.image_type)
        self.finished.emit()
