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
import fnmatch
import os
import json
import pathlib
import shutil
import logging
import ast
import numpy as np
from xml.etree import ElementTree
from pyssa.io_pyssa.xml_pyssa import element_names, attribute_names
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures import settings
from pyssa.internal.data_structures import sequence
from pyssa.internal.analysis_types import distance_analysis
from pyssa.internal.data_structures import results
from pyssa.io_pyssa import safeguard
from pyssa.util import constants
from pyssa.util import protein_util
from pyssa.util import project_util
from pyssa.io_pyssa import path_util

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class XmlDeserializer:
    """"""
    """
    
    """
    xml_root: ElementTree.Element

    def __init__(self, filepath):
        # Read the XML file
        xml_file = open(filepath, "r")
        xml_contents = xml_file.read()
        self.xml_root = ElementTree.fromstring(xml_contents)

    def create_all_proteins_from_xml(self):
        proteins = []
        for tmp_protein in self.xml_root.iter(element_names.PROTEIN):
            basic_information = tmp_protein.attrib
            pdb_lines = []
            session_data_base64 = ""
            for tmp_data in tmp_protein:
                if tmp_data.tag == "pdb_data":
                    for tmp_atom in tmp_data.findall("atom"):
                        pdb_lines.append(tmp_atom.text)
                elif tmp_data.tag == "session_data":
                    session_data_base64 = tmp_data.attrib[attribute_names.PROTEIN_SESSION]
                else:
                    raise ValueError
            tmp_protein_obj = protein.Protein(molecule_object=basic_information[attribute_names.PROTEIN_MOLECULE_OBJECT],
                                              pdb_xml_string=tmp_protein)
            tmp_protein_obj.set_all_attributes(basic_information, pdb_lines, session_data_base64)

            proteins.append(tmp_protein_obj)
        return proteins

    def create_all_protein_pairs(self, tmp_project: 'project.Project', app_settings: 'settings.Settings'):
        protein_pairs = []
        if len(self.xml_root.findall(f".//{element_names.PROTEIN_PAIR}")) == 0:
            return
        for tmp_protein_pair in self.xml_root.findall(f".//{element_names.PROTEIN_PAIR}"):
            basic_information = tmp_protein_pair.attrib
            prot_1_molecule_object = basic_information[attribute_names.PROTEIN_PAIR_PROT_1_MOLECULE_OBJECT]
            if prot_1_molecule_object.find("_1") != -1:
                org_prot_1_molecule_object = prot_1_molecule_object[:prot_1_molecule_object.find("_1")]
                protein_2: protein.Protein = tmp_project.search_protein(org_prot_1_molecule_object).duplicate_protein()
                protein_1: protein.Protein = protein_2.duplicate_protein()
                protein_1.set_molecule_object(f"{protein_1.get_molecule_object()}_1")
                protein_2.set_molecule_object(f"{protein_2.get_molecule_object()}_2")
            else:
                protein_1 = tmp_project.search_protein(basic_information[attribute_names.PROTEIN_PAIR_PROT_1_MOLECULE_OBJECT])
                protein_2 = tmp_project.search_protein(basic_information[attribute_names.PROTEIN_PAIR_PROT_2_MOLECULE_OBJECT])

            tmp_protein_pair_obj = protein_pair.ProteinPair(protein_1=protein_1, protein_2=protein_2)

            pymol_session = tmp_protein_pair.find(element_names.PROTEIN_PAIR_SESSION).attrib
            tag_distance_analysis = tmp_protein_pair.find(element_names.DISTANCE_ANALYSIS)
            distance_analysis_settings = tag_distance_analysis.attrib
            tag_results = tag_distance_analysis.find(element_names.DISTANCE_ANALYSIS_RESULTS)
            rmsd_aligned_aa = tag_results.attrib
            for tmp_tag_results_distance in tag_results.findall(f".//{element_names.DISTANCE_ANALYSIS_DISTANCE_RESULTS}"):
                indexes = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_INDEX_LIST).text
                prot_1_chains = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_PROT_1_CHAIN_LIST).text
                prot_1_positions = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_PROT_1_POSITION_LIST).text
                prot_1_residues = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_PROT_1_RESIDUE_LIST).text

                prot_2_chains = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_PROT_2_CHAIN_LIST).text
                prot_2_positions = tmp_tag_results_distance.find(
                    element_names.DISTANCE_ANALYSIS_PROT_2_POSITION_LIST).text
                prot_2_residues = tmp_tag_results_distance.find(
                    element_names.DISTANCE_ANALYSIS_PROT_2_RESIDUE_LIST).text

                distances = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_DISTANCES_LIST).text
                index_array: np.ndarray = np.array(ast.literal_eval(indexes))
                ref_chain_array: np.ndarray = np.array(ast.literal_eval(prot_1_chains))
                ref_pos_array: np.ndarray = np.array(ast.literal_eval(prot_1_positions))
                ref_resi_array: np.ndarray = np.array(ast.literal_eval(prot_1_residues))
                model_chain_array: np.ndarray = np.array(ast.literal_eval(prot_2_chains))
                model_pos_array: np.ndarray = np.array(ast.literal_eval(prot_2_positions))
                model_resi_array: np.ndarray = np.array(ast.literal_eval(prot_2_residues))
                distance_array: np.ndarray = np.array(ast.literal_eval(distances))

                result_hashtable: dict[str, np.ndarray] = {'index': index_array,
                                                           'ref_chain': ref_chain_array,
                                                           'ref_pos': ref_pos_array,
                                                           'ref_resi': ref_resi_array,
                                                           'model_chain': model_chain_array,
                                                           'model_pos': model_pos_array,
                                                           'model_resi': model_resi_array,
                                                           'distance': distance_array
                                                           }

            tmp_protein_pair_obj.distance_analysis = distance_analysis.DistanceAnalysis(tmp_protein_pair_obj, app_settings)
            tmp_protein_pair_obj.distance_analysis.analysis_results = results.DistanceAnalysisResults(distance_data=result_hashtable,
                                                                                                      pymol_session=pymol_session[attribute_names.DISTANCE_ANALYSIS_SESSION],
                                                                                                      rmsd=float(rmsd_aligned_aa[attribute_names.DISTANCE_ANALYSIS_RMSD]),
                                                                                                      aligned_aa=int(rmsd_aligned_aa[attribute_names.DISTANCE_ANALYSIS_ALIGNED_AA]))
            tmp_protein_pair_obj.distance_analysis.cutoff = float(distance_analysis_settings[attribute_names.DISTANCE_ANALYSIS_CUTOFF])
            tmp_protein_pair_obj.distance_analysis.cycles = distance_analysis_settings[attribute_names.DISTANCE_ANALYSIS_CYCLES]
            tmp_protein_pair_obj.distance_analysis.name = distance_analysis_settings[attribute_names.DISTANCE_ANALYSIS_NAME]
            tmp_protein_pair_obj.distance_analysis.rmsd_dict['rmsd'] = float(rmsd_aligned_aa[attribute_names.DISTANCE_ANALYSIS_RMSD])
            tmp_protein_pair_obj.distance_analysis.rmsd_dict['aligned_residues'] = int(rmsd_aligned_aa[attribute_names.DISTANCE_ANALYSIS_ALIGNED_AA])
            protein_pairs.append(tmp_protein_pair_obj)

        return protein_pairs

    def deserialize_project(self, app_settings):
        project_dict = {}
        for info in self.xml_root.iter(element_names.PROJECT_INFO):
            project_dict = info.attrib
        tmp_project = project.Project(project_dict[attribute_names.PROJECT_NAME], pathlib.Path(project_dict[attribute_names.PROJECT_WORKSPACE_PATH]))
        protein_objs = self.create_all_proteins_from_xml()
        for tmp_protein_obj in protein_objs:
            tmp_project.add_existing_protein(tmp_protein_obj)
        protein_pair_objs = self.create_all_protein_pairs(tmp_project, app_settings)
        if protein_pair_objs is None:
            return tmp_project
        else:
            for tmp_protein_pair_obj in protein_pair_objs:
                tmp_project.add_protein_pair(tmp_protein_pair_obj)
            return tmp_project

    def deserialize_analysis_images(self, protein_pair_name: str, analysis_results: 'results.DistanceAnalysisResults'):
        for tmp_protein_pair in self.xml_root.findall(f".//{element_names.PROTEIN_PAIR}"):
            if tmp_protein_pair.attrib["name"] == protein_pair_name:
                structure_aln_images = tmp_protein_pair.findall(f".//{element_names.DISTANCE_ANALYSIS_IMAGES}/{element_names.DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE}")
                analysis_results.structure_aln_image = (structure_aln_images[0].attrib[attribute_names.DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE_BASENAME], structure_aln_images[0].text)
                interesting_reg_images = tmp_protein_pair.findall(f".//{element_names.DISTANCE_ANALYSIS_IMAGES}/{element_names.DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS}")
                for tmp_image in interesting_reg_images:
                    analysis_results.interesting_regions_images.append((tmp_image.attrib[attribute_names.DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS_BASENAME], tmp_image.text))


class ObjectSerializer:
    object_dict: dict

    def __init__(self, object_to_serialize, filepath, filename):
        """Constructor

        Args:
            object_to_serialize:
                objects which gets serialized
            filepath:
                filepath of the json file, without the file name
            filename:
                name of the json file WITHOUT .json extension!

        """
        self.object = object_to_serialize
        self.filepath = filepath
        self.filename = filename

    def create_standard_object_dict(self):
        self.object_dict = self.object.__dict__
    
    def set_custom_object_dict(self, custom_dict: dict):
        self.object_dict = custom_dict
    
    def update_object_dict(self, dict_with_updates: dict):
        self.object_dict.update(dict_with_updates)
    
    def serialize_object(self):
        tmp_object_file = open(pathlib.Path(f"{self.filepath}/{self.filename}.json"), "w", encoding="utf-8")
        json.dump(self.object_dict, tmp_object_file, indent=4)
        tmp_object_file.close()


class ObjectDeserializer:
    
    def __init__(self, filepath, filename):
        tmp_object_file = open(pathlib.Path(f"{filepath}/{filename}.json"), "r", encoding="utf-8")
        self.object_dict = json.load(tmp_object_file)
    
    def deserialize_protein(self) -> 'protein.Protein':
        if self.object_dict.get("export_data_dir") == "None":
            update = {"export_data_dir": None}
            self.object_dict.update(update)
        tmp_protein = protein.Protein(molecule_object=self.object_dict.get("filename"), # important for duplicated proteins
                                      proteins_dirname=self.object_dict.get("proteins_dirname"),
                                      pdb_filepath=path_util.FilePath(self.object_dict.get("import_data_dir")))
        tmp_protein.molecule_object = self.object_dict.get("molecule_object")
        tmp_protein.pymol_selection.selection_string = self.object_dict.get("selection")
        tmp_protein.chains = protein_util.create_chains_from_list_of_tuples(self.object_dict.get("chains"))
        tmp_protein.pymol_session_file = self.object_dict.get("pymol_session_file")
        return tmp_protein
    
    def deserialize_protein_pair(self) -> 'protein_pair.ProteinPair':
        tmp_protein_pair = protein_pair.ProteinPair(self.object_dict.get("prot_1_molecule_object"),
                                                    self.object_dict.get("prot_2_molecule_object"),
                                                    pathlib.Path(self.object_dict.get("export_dirname")))
        tmp_protein_pair.pymol_session_file = self.object_dict.get("pymol_session_file")
        tmp_protein_pair.name = self.object_dict.get("name")

        if self.object_dict.get("prot_1_export_data_dir") == "None":
            export_data_dir = None
        else:
            export_data_dir = self.object_dict.get("prot_1_export_data_dir")

        tmp_protein_1 = protein.Protein(molecule_object=self.object_dict.get("prot_1_filename"), # important for duplicated proteins
                                        proteins_dirname=self.object_dict.get("proteins_dirname"),
                                        pdb_filepath=self.object_dict.get("prot_1_import_data_dir"))
        tmp_protein_1.molecule_object = self.object_dict.get("prot_1_molecule_object")
        tmp_protein_1.pymol_selection.selection_string = self.object_dict.get("prot_1_selection")
        tmp_protein_1.chains = protein_util.create_chains_from_list_of_tuples(self.object_dict.get("prot_1_chains"))
        tmp_protein_pair.ref_obj = tmp_protein_1

        if self.object_dict.get("prot_2_export_data_dir") == "None":
            export_data_dir = None
        else:
            export_data_dir = self.object_dict.get("prot_2_export_data_dir")
        tmp_protein_2 = protein.Protein(molecule_object=self.object_dict.get("prot_2_filename"),
                                        # important for duplicated proteins
                                        proteins_dirname=self.object_dict.get("proteins_dirname"),
                                        pdb_filepath=self.object_dict.get("prot_2_import_data_dir"))
        tmp_protein_2.molecule_object = self.object_dict.get("prot_2_molecule_object")
        tmp_protein_2.pymol_selection.selection_string = self.object_dict.get("prot_2_selection")
        tmp_protein_2.chains = protein_util.create_chains_from_list_of_tuples(self.object_dict.get("prot_2_chains"))
        tmp_protein_pair.model_obj = tmp_protein_2

        return tmp_protein_pair
    
    def deserialize_project(self, app_settings: 'settings.Settings') -> 'project.Project':
        tmp_project: project.Project = project.Project(self.object_dict.get("_project_name"), self.object_dict.get("_workspace"))
        tmp_project.set_folder_paths(self.object_dict.get("_folder_paths"))
        tmp_project.project_path = self.object_dict.get("project_path")
        for tmp_protein_json_filepath in project_util.get_all_protein_json_filepaths_from_project(tmp_project):
            tmp_project.add_existing_protein(
                ObjectDeserializer(
                    tmp_protein_json_filepath.get_dirname(),
                    tmp_protein_json_filepath.get_filename()
                ).deserialize_protein()
            )
        for tmp_protein_pair_json_filepath in project_util.get_all_protein_pair_json_filepaths_from_project(tmp_project):
            tmp_project.add_protein_pair(
                ObjectDeserializer(
                    tmp_protein_pair_json_filepath.get_dirname(),
                    tmp_protein_pair_json_filepath.get_filename()
                ).deserialize_protein_pair()
            )
        for tmp_distance_analysis_json_filepath in project_util.get_all_distance_analysis_json_filepaths_from_project(tmp_project):
            tmp_project.add_distance_analysis(
                ObjectDeserializer(
                    tmp_distance_analysis_json_filepath.get_dirname(),
                    tmp_distance_analysis_json_filepath.get_filename()
                ).deserialize_distance_analysis(app_settings=app_settings, app_project=tmp_project)
            )
        return tmp_project

    def deserialize_distance_analysis(self, app_settings: 'settings.Settings', app_project: 'project.Project') -> 'distance_analysis.DistanceAnalysis':
        tmp_protein_pair = app_project.get_specific_protein_pair(self.object_dict.get("protein_pair_for_analysis_name"))
        tmp_distance_analysis = distance_analysis.DistanceAnalysis(tmp_protein_pair, app_settings, self.object_dict.get("distance_analysis_dirname"))
        tmp_distance_analysis.cutoff = self.object_dict.get("cutoff")
        tmp_distance_analysis.cycles = self.object_dict.get("cycles")
        tmp_distance_analysis.export_dirname = self.object_dict.get("export_dirname")
        tmp_distance_analysis.pymol_session_filepath = self.object_dict.get("pymol_session_filepath")
        return tmp_distance_analysis

    def deserialize_settings(self) -> 'settings.Settings':
        tmp_settings: settings.Settings = settings.Settings(self.object_dict.get("dir_settings"), self.object_dict.get("filename"))
        if safeguard.Safeguard.check_filepath(self.object_dict.get("workspace_path")):
            tmp_settings.workspace_path = self.object_dict.get("workspace_path")
        else:
            raise ValueError
        if safeguard.Safeguard.check_filepath(self.object_dict.get("prediction_path")):
            tmp_settings.prediction_path = self.object_dict.get("prediction_path")
        else:
            raise ValueError
        if safeguard.Safeguard.check_if_number_is_positive(int(self.object_dict.get("cycles"))):
            tmp_settings.cycles = self.object_dict.get("cycles")
        else:
            raise ValueError
        if safeguard.Safeguard.check_if_number_is_positive(float(self.object_dict.get("cutoff"))):
            tmp_settings.cutoff = self.object_dict.get("cutoff")
        else:
            raise ValueError
        if int(self.object_dict.get("app_launch")) == 0 or int(self.object_dict.get("app_launch")) == 1:
            tmp_settings.app_launch = self.object_dict.get("app_launch")
        else:
            raise ValueError
        if not safeguard.Safeguard.check_filepath(self.object_dict.get("dir_settings")):
            raise ValueError
        if not self.object_dict.get("filename") == "settings.json":
            raise ValueError
        if int(self.object_dict.get("wsl_install")) == 0 or int(self.object_dict.get("wsl_install")) == 1:
            tmp_settings.wsl_install = int(self.object_dict.get("wsl_install"))
        else:
            raise ValueError
        if int(self.object_dict.get("local_colabfold")) == 0 or int(self.object_dict.get("local_colabfold")) == 1:
            tmp_settings.local_colabfold = int(self.object_dict.get("local_colabfold"))
        else:
            raise ValueError
        tmp_settings.wsl_username = self.object_dict.get("wsl_username")
        return tmp_settings

    def deserialize_sequence(self) -> 'sequence.Sequence':
        """This function deserializes the sequence object from a .json file.

        Returns:
            a ProteinSequence object
        """
        return sequence.Sequence(self.object_dict.get("name"), self.object_dict.get("sequence"))


class ProjectScanner:

    def __init__(self, project_obj: 'project.Project'):
        self.project = project_obj

    def scan_project_for_valid_proteins(self, list_view_project_proteins=None) -> list[str]:
        """This function scans the project pdb path and optionally fills a list view

        Args:
            list_view_project_proteins:

        Returns:

        """
        pattern = "*.pdb"
        protein_basenames = []
        # iterates over possible project directories
        if list_view_project_proteins is not None:
            for tmp_protein_dir in os.listdir(self.project.folder_paths["project"]):
                path = pathlib.Path(f"{self.project.folder_paths['project']}/{tmp_protein_dir}")
                if len(os.listdir(pathlib.Path(f"{path}/pdb"))) > 1:
                    logger.error("Too many pdb files in one directory.")
                    raise ValueError("Too many pdb files in one directory.")
                if fnmatch.fnmatch(str(os.listdir(pathlib.Path(f"{path}/pdb"))[0]), pattern) and list_view_project_proteins is not None:
                    list_view_project_proteins.addItem(str(os.listdir(pathlib.Path(f"{path}/pdb"))[0]))
                elif fnmatch.fnmatch(str(os.listdir(pathlib.Path(f"{path}/pdb"))[0]), pattern) and list_view_project_proteins is None:
                    protein_basenames.append(str(os.listdir(pathlib.Path(f"{path}/pdb"))[0]))
        return protein_basenames


class WorkspaceScanner:

    def __init__(self, workspace_path):
        self.workspace_path = workspace_path

    def scan_workspace_for_valid_projects(self, list_new_projects) -> list:
        workspace_projects: list[str] = os.listdir(self.workspace_path)
        valid_directories = []
        # iterates over possible project directories
        for directory in workspace_projects:
            try:
                directory_content = os.listdir(f"{self.workspace_path}/{directory}")
                # iterates over the content in a single project directory
                for content in directory_content:
                    if content == "project.json":
                        valid_directories.append(directory)
            except NotADirectoryError:
                print(f"This: {directory} is not a directory.")

        valid_directories.sort()
        for project in valid_directories:
            list_new_projects.addItem(project)
        return valid_directories

    def scan_workspace_for_non_duplicate_proteins(self) -> dict:
        """This function scans the workspace directory for protein structures and eliminates all duplicates

        Args:
            workspace_path:
                the path of the workspace
        """
        pdb_basenames_filepaths = []
        for tmp_project_name in os.listdir(self.workspace_path):
            for tmp_protein_name in os.listdir(pathlib.Path(f"{self.workspace_path}/{tmp_project_name}/proteins")):
                pdb_basename = os.listdir(pathlib.Path(f"{self.workspace_path}/{tmp_project_name}/proteins/{tmp_protein_name}/pdb"))
                if len(pdb_basename) > 1:
                    logger.error("Too many pdb files in one directory.")
                    raise ValueError("Too many pdb files in one directory.")
                pdb_filepath = path_util.FilePath(pathlib.Path(f"{self.workspace_path}/{tmp_project_name}/proteins/{tmp_protein_name}/pdb/{pdb_basename[0]}"))
                pdb_basenames_filepaths.append((pdb_basename[0], pdb_filepath))
        return dict(list(set(pdb_basenames_filepaths)))
        # """Var: workspace_proteins is a list which contains all proteins from all projects in the workspace"""
        # workspace_proteins = []
        # protein_names = []
        # protein_tuples_notation = []
        # for valid_project in valid_projects:
        #     # if valid_project != current_project_name: # I don't know why this if-statement should be important
        #     """Var: project_proteins is a list which contains all proteins from a single project"""
        #     project_proteins = tools.scan_project_for_valid_proteins(pathlib.Path(f"{self.workspace_path}/{valid_project}"),
        #                                                              list_widget)
        #     list_widget.clear()
        #     for protein in project_proteins:
        #         tmp_protein = protein_info.ProteinInfo(protein,
        #                                                pathlib.Path(f"{self.workspace_path}/{valid_project}/pdb/{protein}"))
        #         workspace_proteins.append(tmp_protein)
        #         if tmp_protein.name not in protein_names:
        #             protein_names.append(tmp_protein.name)
        # # this for-loop is necessary for the creation of the protein dictionary
        # for protein in workspace_proteins:
        #     protein_tuples_notation.append(protein.get_tuple_notation())
        # protein_dict = dict(protein_tuples_notation)
        # return protein_dict, protein_names


class FilesystemCleaner:

    def __int__(self):
        pass

    @staticmethod
    def clean_prediction_scratch_folder():
        shutil.rmtree(constants.PREDICTION_FASTA_DIR)
        os.mkdir(constants.PREDICTION_FASTA_DIR)
        shutil.rmtree(constants.PREDICTION_PDB_DIR)
        os.mkdir(constants.PREDICTION_PDB_DIR)
