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
"""Module for handling filesystem operations."""
import concurrent
import concurrent.futures
import collections
import copy
import itertools
from queue import Queue
from threading import Lock
import os
import json
import pathlib
import shutil
import logging
import ast
import collections
import numpy as np
from xml import sax
from xml.etree import ElementTree
from PyQt5 import QtCore
from Bio import SeqRecord
from pyssa.internal.data_structures.data_classes import basic_protein_info
from pyssa.io_pyssa.xml_pyssa import element_names, attribute_names
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures import protein, structure_analysis, chain, selection
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures import settings
from pyssa.internal.data_structures import sequence
from pyssa.internal.analysis_types import distance_analysis
from pyssa.internal.data_structures import results
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa import path_util
from pyssa.util import constants
from pyssa.util import protein_util
from pyssa.util import project_util

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class ProjectParserHandler(sax.ContentHandler):
    def __init__(self, a_project, the_settings):
        self.one_element_history: collections.deque = collections.deque()
        self.one_element_history.append("root")
        self.current_element = ""
        self.protein_data = {}
        self._project = a_project
        self._app_settings = the_settings
        self.pdb_data = []
        self.chains = []
        self._protein_in_pair: bool = False
        self.protein_pair_data = {}
        self._proteins_for_pair = []
        self.result = []
        self.result_obj = None
        self.distance_analysis = {}
        self.distance_analysis_obj = None
        self.current_content = ""

        index_list = np.array([])
        ref_chain_list = np.array([])
        ref_pos_list = np.array([])
        ref_resi_list = np.array([])
        model_chain_list = np.array([])
        model_pos_list = np.array([])
        model_resi_list = np.array([])
        distance_list = np.array([])
        self.distances_results = {
            "index": index_list,
            "ref_chain": ref_chain_list,
            "ref_pos": ref_pos_list,
            "ref_resi": ref_resi_list,
            "model_chain": model_chain_list,
            "model_pos": model_pos_list,
            "model_resi": model_resi_list,
            "distance": distance_list,
        }

    def startElement(self, name, attrs):
        self.current_element = name
        tmp_last_element = self.one_element_history[-1]
        if name == "project_info":
            self._project = project.Project(attrs.getValue("name"), pathlib.Path(attrs.getValue("workspace_path")))
        elif name == "sequence":
            if tmp_last_element == "sequences":
                self._project.sequences.append(
                    SeqRecord.SeqRecord(name=attrs.getValue("name"), seq=attrs.getValue("seq"))
                )
        elif name == "protein":
            self.protein_data = {
                "id": attrs.getValue("id"),
                "pymol_molecule_object": attrs.getValue("pymol_molecule_object"),
                "pymol_selection": attrs.getValue("pymol_selection")
            }
            if tmp_last_element == "proteins":
                self._protein_in_pair = False
            elif tmp_last_element == "proteins_of_pair":
                self._protein_in_pair = True
        elif name == "chain":
            tmp_sequence = sequence.Sequence(self.protein_data["pymol_molecule_object"], attrs.getValue("chain_sequence"))
            tmp_chain = chain.Chain(attrs.getValue("chain_letter"),
                                    tmp_sequence,
                                    attrs.getValue("chain_type"))
            self.chains.append(tmp_chain)
        elif name == "pymol_parameters":
            tmp_last_chain: "chain.Chain" = self.chains[-1]
            tmp_last_chain.pymol_parameters = {
                "chain_color": attrs.getValue("chain_color"),
                "chain_representation": attrs.getValue("chain_representation"),
            }
        elif name == "session_data":
            if tmp_last_element == "atom":
                self.protein_data["session"] = attrs.getValue("session")
            elif tmp_last_element == "distances":
                self.result.append(attrs.getValue("session"))
        elif name == "protein_pair":
            self.protein_pair_data = {
                "name": attrs.getValue("name")
            }
        elif name == "distance_analysis":
            self.distance_analysis = {
                "name": attrs.getValue("name"),
                "cutoff": attrs.getValue("cutoff"),
                "cycles": attrs.getValue("cycles")
            }
        elif name == "result":
            self.result.append(attrs.getValue("rmsd"))
            self.result.append(attrs.getValue("aligned_amino_acids"))

        self.one_element_history.append(name)
        self.one_element_history.popleft()

    def characters(self, content):
        if content == "":
            return
        if self.current_element == "atom":
            self.pdb_data.append(content)
        # Distance results
        elif self.current_element == "index":
            self.current_content += content
            print(f"index content: {content}")
        elif self.current_element == "ref_chain":
            self.current_content += content
        elif self.current_element == "ref_pos":
            self.current_content += content
        elif self.current_element == "ref_resi":
            self.current_content += content
        elif self.current_element == "model_chain":
            self.current_content += content
        elif self.current_element == "model_pos":
            self.current_content += content
        elif self.current_element == "model_resi":
            self.current_content += content
        elif self.current_element == "distance":
            self.current_content += content

    def endElement(self, name):
        tmp_last_element = self.one_element_history[-1]
        if name == "protein":
            # Creates new protein object
            tmp_protein = protein.Protein(self.protein_data["pymol_molecule_object"], pdb_data=self.pdb_data)
            tmp_protein.set_id(self.protein_data["id"])
            tmp_protein.chains = copy.deepcopy(self.chains)
            tmp_protein.pymol_session = self.protein_data["session"]
            tmp_protein.pymol_selection = selection.Selection(self.protein_data["pymol_molecule_object"])
            tmp_protein.pymol_selection.selection_string = self.protein_data["pymol_selection"]
            # Clears vars for next protein
            self.protein_data = {}
            self.chains = []
            if not self._protein_in_pair:
                # Appends protein to project for the protein model
                self._project.proteins.append(copy.deepcopy(tmp_protein))
            elif self._protein_in_pair:
                # Appends protein to project for the protein pair model
                self._proteins_for_pair.append(copy.deepcopy(tmp_protein))
        elif self.current_element == "index":
            if self.current_content != "":
                self.distances_results["index"] = np.array(self.current_content)
                print(f"After np.array conversion: {self.distances_results['index']}")
            self.current_content = ""
        elif self.current_element == "ref_chain":
            self.distances_results["ref_chain"] = np.array(ast.literal_eval(self.current_content))
            self.current_content = ""
        elif self.current_element == "ref_pos":
            self.distances_results["ref_pos"] = np.array(ast.literal_eval(self.current_content))
            self.current_content = ""
        elif self.current_element == "ref_resi":
            self.distances_results["ref_resi"] = np.array(ast.literal_eval(self.current_content))
            self.current_content = ""
        elif self.current_element == "model_chain":
            self.distances_results["model_chain"] = np.array(ast.literal_eval(self.current_content))
            self.current_content = ""
        elif self.current_element == "model_pos":
            self.distances_results["model_pos"] = np.array(ast.literal_eval(self.current_content))
            self.current_content = ""
        elif self.current_element == "model_resi":
            self.distances_results["model_resi"] = np.array(ast.literal_eval(self.current_content))
            self.current_content = ""
        elif self.current_element == "distance":
            self.distances_results["distance"] = np.array(ast.literal_eval(self.current_content))
            self.current_content = ""
        elif name == "result":
            print(self.distances_results)
            index_array: np.ndarray = np.array(self.distances_results["index"])
            ref_chain_array: np.ndarray = np.array(self.distances_results["ref_chain"])
            ref_pos_array: np.ndarray = np.array(self.distances_results["ref_pos"])
            ref_resi_array: np.ndarray = np.array(self.distances_results["ref_resi"])
            model_chain_array: np.ndarray = np.array(self.distances_results["model_chain"])
            model_pos_array: np.ndarray = np.array(self.distances_results["model_pos"])
            model_resi_array: np.ndarray = np.array(self.distances_results["model_resi"])
            distance_array: np.ndarray = np.array(self.distances_results["distance"])
            result_hashtable: dict[str, np.ndarray] = {
                "index": index_array,
                "ref_chain": ref_chain_array,
                "ref_pos": ref_pos_array,
                "ref_resi": ref_resi_array,
                "model_chain": model_chain_array,
                "model_pos": model_pos_array,
                "model_resi": model_resi_array,
                "distance": distance_array,
            }
            self.result_obj = results.DistanceAnalysisResults(result_hashtable,
                                                              self.result[2],
                                                              self.result[0],
                                                              self.result[1])
        elif name == "distance_analysis":
            self.distance_analysis_obj = structure_analysis.DistanceAnalysis(
                self._app_settings, distance_analysis_name=self.distance_analysis["name"]
            )
            self.distance_analysis_obj.analysis_results = self.result_obj
            self.distance_analysis_obj.cutoff = self.distance_analysis["cutoff"]
            self.distance_analysis_obj.cycles = self.distance_analysis["cycles"]
        elif name == "protein_pair":
            tmp_protein_pair = protein_pair.ProteinPair(self._proteins_for_pair[0], self._proteins_for_pair[1])
            tmp_protein_pair.name = self.protein_pair_data["name"]
            tmp_protein_pair.distance_analysis = self.distance_analysis_obj
            self._project.protein_pairs.append(tmp_protein_pair)
            self.protein_pair_data = {}
            self._proteins_for_pair = []
            self.result = []
            self.result_obj = None
            self.distance_analysis = {}
            self.distance_analysis_obj = None

    def get_project(self):
        return self._project

    # def __init__(self, a_project):
    #     self.current_element = ""
    #     self.protein_data = {}
    #     self.pdb_data = []
    #     self.proteins = []
    #     self._project = a_project
    #
    # def startElement(self, name, attrs):
    #     self.current_element = name
    #     if name == "protein":
    #         self.protein_data = {"id": attrs.getValue("id"), "pymol_molecule_object": attrs.getValue("pymol_molecule_object")}
    #     elif name == "protein_pair":
    #         pass
    #
    # def characters(self, content):
    #     if self.current_element == "atom":
    #         self.pdb_data.append(content)
    #         #self.protein_data["pdb_data"][-1] += content
    #
    # def endElement(self, name):
    #     if name == "protein":
    #         tmp_protein = protein.Protein(self.protein_data["id"], pdb_data=self.pdb_data)
    #         self.proteins.append(tmp_protein)
    #         self.protein_data = {}
    #         self._project.proteins.append(tmp_protein)


class ProjectParser:
    _project: "project.Project"

    def __init__(self, xml_reader: QtCore.QXmlStreamReader):
        self.xml_reader = xml_reader
        self._project = project.Project()

    def parse_complete_xml(self):
        while not self.xml_reader.atEnd() and not self.xml_reader.hasError():
            token = self.xml_reader.readNext()
            element_name = self.xml_reader.name()
            if element_name == "project_info":
                # Process project_info element
                self._project = project.Project(self.xml_reader.attributes().value("name"),
                                                pathlib.Path(self.xml_reader.attributes().value("workspace_path")))
            elif element_name == "protein":
                # Process protein element
                id = self.xml_reader.attributes().value("id")
                molecule_object = self.xml_reader.attributes().value("pymol_molecule_object")
                selection = self.xml_reader.attributes().value("pymol_selection")
                # Continue reading the elements within the protein element
                while self.xml_reader.readNextStartElement():
                    sub_element = self.xml_reader.name()
                    if sub_element == "pdb_data":
                        # Process pdb_data element
                        # Continue reading the elements within the pdb_data element
                        tmp_pdb_data: collections.deque = collections.deque()
                        while self.xml_reader.readNextStartElement():
                            atom_element = self.xml_reader.name()
                            if atom_element == "atom":
                                # Process atom element
                                tmp_pdb_data.append(self.xml_reader.readElementText())
                                # Process atom data as needed
                            else:
                                self.xml_reader.skipCurrentElement()  # Skip unknown elements
                    elif sub_element == "session_data":
                        # Process session_data element
                        session_data = self.xml_reader.readElementText()
                        # Process session_data as needed
                    else:
                        self.xml_reader.skipCurrentElement()  # Skip unknown elements
                tmp_protein = protein.Protein(molecule_object, pdb_data=list(tmp_pdb_data))
                print(tmp_protein.get_pdb_data())
            elif element_name == "protein_pair":
                # Process protein_pair element
                pair_name = self.xml_reader.attributes().value("name").toString()
                # Process other attributes or data as needed

            else:
                self.xml_reader.skipCurrentElement()  # Skip unknown elements

        if self.xml_reader.hasError():
            print("Error parsing XML:", self.xml_reader.errorString())

    def parse_project(self):
        while not self.xml_reader.atEnd():
            self.xml_reader.readNext()
            if self.xml_reader.isStartElement() and self.xml_reader.name() == "project":
                self._project = self.parse_project_info()
                self.parse_sequences()
                self.parse_proteins()
                #self.parse_protein_pairs()
        return self._project

    def parse_project_info(self):
        while not self.xml_reader.atEnd():
            self.xml_reader.readNext()
            if self.xml_reader.isStartElement() and self.xml_reader.name() == "project_info":
                # Process attributes
                self.xml_reader.attributes()
                project_infos = []
                for attribute in self.xml_reader.attributes():
                    attribute_name = attribute.name()
                    attribute_value = attribute.value()
                    project_infos.append((attribute_name, attribute_value))
                return project.Project(project_infos[0][1], pathlib.Path(project_infos[1][1]))

    def parse_sequences(self):
        while not self.xml_reader.atEnd():
            self.xml_reader.readNext()
            if self.xml_reader.isStartElement() and self.xml_reader.name() == "sequences":
                self.parse_sequence()

    def parse_sequence(self):
        sequence_data = {}
        while not self.xml_reader.atEnd():
            self.xml_reader.readNext()
            if self.xml_reader.isStartElement() and self.xml_reader.name() == "sequence":
                # Process attributes
                for attribute in self.xml_reader.attributes():
                    attribute_name = attribute.name()
                    attribute_value = attribute.value()
                    sequence_data[attribute_name] = attribute_value
                self._project.sequences.append(SeqRecord.SeqRecord(name=sequence_data["name"], seq=sequence_data["seq"]))

    def parse_proteins(self):
        while not self.xml_reader.atEnd():
            self.xml_reader.readNext()
            if self.xml_reader.isStartElement() and self.xml_reader.name() == "proteins":
                self.parse_protein()

    def parse_protein(self):
        protein_data = {}
        while not self.xml_reader.atEnd():
            self.xml_reader.readNext()
            if self.xml_reader.isStartElement() and self.xml_reader.name() == "protein":
                # Process attributes
                for attribute in self.xml_reader.attributes():
                    attribute_name = attribute.name()
                    attribute_value = attribute.value()
                    protein_data[attribute_name] = attribute_value
                    print(protein_data)

    def parse_protein_pairs(self):
        while not self.xml_reader.atEnd():
            self.xml_reader.readNext()
            if self.xml_reader.isStartElement() and self.xml_reader.name() == "protein_pairs":
                self.xml_reader.attributes()
            # Process individual protein pairs


class XmlDeserializer:
    """Class used to deserialize XML documents into the respective objects."""

    """
    the root element of the xml file
    """
    xml_root: ElementTree.Element

    def __init__(self, filepath: pathlib.Path) -> None:
        """Constructor.

        Args:
            filepath: the filepath to the xml file.
        """
        # Read the XML file
        xml_file = open(str(filepath), "r")
        xml_contents = xml_file.read()
        self.xml_root = ElementTree.fromstring(xml_contents)

    def create_all_proteins_from_xml(self) -> list["protein.Protein"]:
        """Creates a list of proteins from the XML."""
        proteins: list["protein.Protein"] = []
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
            tmp_protein_obj = protein.Protein(
                molecule_object=basic_information[attribute_names.PROTEIN_MOLECULE_OBJECT],
                pdb_xml_string=tmp_protein,
            )
            tmp_protein_obj.set_all_attributes(basic_information, pdb_lines, session_data_base64)

            proteins.append(tmp_protein_obj)
        return proteins

    def create_all_protein_pairs(
        self,
        a_project: "project.Project",
        app_settings: "settings.Settings",
    ) -> list["protein_pair.ProteinPair"]:
        """Creates a list of protein pairs from the XML.

        Args:
            a_project: a project of which the protein pairs should be used.
            app_settings: the settings of the app.
        """
        protein_pairs = []
        if len(self.xml_root.findall(f".//{element_names.PROTEIN_PAIR}")) == 0:
            return  # noqa: RET502 #TODO: fix this
        for tmp_protein_pair in self.xml_root.findall(f".//{element_names.PROTEIN_PAIR}"):
            basic_information = tmp_protein_pair.attrib
            prot_1_molecule_object = basic_information[attribute_names.PROTEIN_PAIR_PROT_1_MOLECULE_OBJECT]
            if prot_1_molecule_object.find("_1") != -1:
                org_prot_1_molecule_object = prot_1_molecule_object[: prot_1_molecule_object.find("_1")]
                protein_2: protein.Protein = a_project.search_protein(org_prot_1_molecule_object).duplicate_protein()
                protein_1: protein.Protein = protein_2.duplicate_protein()
                protein_1.set_molecule_object(f"{protein_1.get_molecule_object()}_1")
                protein_2.set_molecule_object(f"{protein_2.get_molecule_object()}_2")
            else:
                protein_1 = a_project.search_protein(
                    basic_information[attribute_names.PROTEIN_PAIR_PROT_1_MOLECULE_OBJECT],
                )
                protein_2 = a_project.search_protein(
                    basic_information[attribute_names.PROTEIN_PAIR_PROT_2_MOLECULE_OBJECT],
                )

            tmp_protein_pair_obj = protein_pair.ProteinPair(protein_1=protein_1, protein_2=protein_2)
            tmp_protein_pair_obj.name = tmp_protein_pair.attrib["name"]
            pymol_session = tmp_protein_pair.find(element_names.PROTEIN_PAIR_SESSION).attrib
            tag_distance_analysis = tmp_protein_pair.find(element_names.DISTANCE_ANALYSIS)
            tmp_protein_pair_obj.pymol_session = pymol_session[attribute_names.DISTANCE_ANALYSIS_SESSION]
            distance_analysis_settings = tag_distance_analysis.attrib
            tag_results = tag_distance_analysis.find(element_names.DISTANCE_ANALYSIS_RESULTS)
            rmsd_aligned_aa = tag_results.attrib
            for tmp_tag_results_distance in tag_results.findall(
                f".//{element_names.DISTANCE_ANALYSIS_DISTANCE_RESULTS}",
            ):
                indexes = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_INDEX_LIST).text
                prot_1_chains = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_PROT_1_CHAIN_LIST).text
                prot_1_positions = tmp_tag_results_distance.find(
                    element_names.DISTANCE_ANALYSIS_PROT_1_POSITION_LIST,
                ).text
                prot_1_residues = tmp_tag_results_distance.find(
                    element_names.DISTANCE_ANALYSIS_PROT_1_RESIDUE_LIST,
                ).text

                prot_2_chains = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_PROT_2_CHAIN_LIST).text
                prot_2_positions = tmp_tag_results_distance.find(
                    element_names.DISTANCE_ANALYSIS_PROT_2_POSITION_LIST,
                ).text
                prot_2_residues = tmp_tag_results_distance.find(
                    element_names.DISTANCE_ANALYSIS_PROT_2_RESIDUE_LIST,
                ).text

                distances = tmp_tag_results_distance.find(element_names.DISTANCE_ANALYSIS_DISTANCES_LIST).text
                index_array: np.ndarray = np.array(ast.literal_eval(indexes))
                ref_chain_array: np.ndarray = np.array(ast.literal_eval(prot_1_chains))
                ref_pos_array: np.ndarray = np.array(ast.literal_eval(prot_1_positions))
                ref_resi_array: np.ndarray = np.array(ast.literal_eval(prot_1_residues))
                model_chain_array: np.ndarray = np.array(ast.literal_eval(prot_2_chains))
                model_pos_array: np.ndarray = np.array(ast.literal_eval(prot_2_positions))
                model_resi_array: np.ndarray = np.array(ast.literal_eval(prot_2_residues))
                distance_array: np.ndarray = np.array(ast.literal_eval(distances))

                result_hashtable: dict[str, np.ndarray] = {
                    "index": index_array,
                    "ref_chain": ref_chain_array,
                    "ref_pos": ref_pos_array,
                    "ref_resi": ref_resi_array,
                    "model_chain": model_chain_array,
                    "model_pos": model_pos_array,
                    "model_resi": model_resi_array,
                    "distance": distance_array,
                }

            tmp_protein_pair_obj.distance_analysis = structure_analysis.DistanceAnalysis(
                protein_pair_name=tmp_protein_pair_obj.name,
                the_app_settings=app_settings,
            )
            tmp_protein_pair_obj.distance_analysis.analysis_results = results.DistanceAnalysisResults(
                distance_data=result_hashtable,
                pymol_session=pymol_session[attribute_names.DISTANCE_ANALYSIS_SESSION],
                rmsd=float(rmsd_aligned_aa[attribute_names.DISTANCE_ANALYSIS_RMSD]),
                aligned_aa=str(rmsd_aligned_aa[attribute_names.DISTANCE_ANALYSIS_ALIGNED_AA]),
            )
            tmp_protein_pair_obj.distance_analysis.cutoff = float(
                distance_analysis_settings[attribute_names.DISTANCE_ANALYSIS_CUTOFF],
            )
            tmp_protein_pair_obj.distance_analysis.cycles = distance_analysis_settings[
                attribute_names.DISTANCE_ANALYSIS_CYCLES
            ]
            tmp_protein_pair_obj.distance_analysis.name = distance_analysis_settings[
                attribute_names.DISTANCE_ANALYSIS_NAME
            ]
            tmp_protein_pair_obj.distance_analysis.analysis_results.rmsd = float(
                rmsd_aligned_aa[attribute_names.DISTANCE_ANALYSIS_RMSD],
            )
            tmp_protein_pair_obj.distance_analysis.analysis_results.aligned_aa = str(
                rmsd_aligned_aa[attribute_names.DISTANCE_ANALYSIS_ALIGNED_AA],
            )
            self.deserialize_analysis_images(
                tmp_protein_pair_obj.name,
                tmp_protein_pair_obj.distance_analysis.analysis_results,
            )
            protein_pairs.append(tmp_protein_pair_obj)

        return protein_pairs

    def create_all_sequences_from_xml(self):
        seq_records: list = []
        for tmp_seq in self.xml_root.iter(element_names.SEQUENCE):
            basic_information = tmp_seq.attrib
            tmp_seq_record_obj = SeqRecord.SeqRecord(basic_information[attribute_names.SEQUENCE_SEQ], name=basic_information[attribute_names.SEQUENCE_NAME])
            seq_records.append(tmp_seq_record_obj)
        return seq_records

    @staticmethod
    def deserialize_project(file_path, app_settings: "settings.Settings") -> "project.Project":
        """Deserialize the project from the XML."""
        file = QtCore.QFile(str(file_path))
        try:
            if not file.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text):
                print("Error: Cannot open file for reading")
        except Exception as e:
            print(f"Exception occurred: {e}")
        else:
            tmp_project = project.Project()
            handler = ProjectParserHandler(tmp_project, app_settings)
            parser = sax.make_parser()
            parser.setContentHandler(handler)
            parser.parse(file_path)
            file.close()
            return handler.get_project()

        # project_dict = {}
        # for info in self.xml_root.iter(element_names.PROJECT_INFO):
        #     project_dict = info.attrib
        # tmp_project = project.Project(
        #     project_dict[attribute_names.PROJECT_NAME],
        #     pathlib.Path(project_dict[attribute_names.PROJECT_WORKSPACE_PATH]),
        # )
        # tmp_seq_record_objs = self.create_all_sequences_from_xml()
        # if tmp_seq_record_objs is not None:
        #     for tmp_seq_record_obj in tmp_seq_record_objs:
        #         tmp_project.sequences.append(tmp_seq_record_obj)
        # protein_objs = self.create_all_proteins_from_xml()
        # for tmp_protein_obj in protein_objs:
        #     tmp_project.add_existing_protein(tmp_protein_obj)
        # protein_pair_objs = self.create_all_protein_pairs(tmp_project, app_settings)
        # if protein_pair_objs is not None:
        #     for tmp_protein_pair_obj in protein_pair_objs:
        #         tmp_project.add_protein_pair(tmp_protein_pair_obj)
        #
        # return tmp_project

    def deserialize_analysis_images(
        self,
        protein_pair_name: str,
        analysis_results: "results.DistanceAnalysisResults",
    ) -> None:
        """Deserialize the analysis images.

        Args:
            protein_pair_name: the name of the protein pair to use.
            analysis_results: the results of the distance analysis.
        """
        for tmp_protein_pair in self.xml_root.findall(f".//{element_names.PROTEIN_PAIR}"):
            if tmp_protein_pair.attrib["name"] == protein_pair_name:
                structure_aln_images = tmp_protein_pair.findall(
                    f".//{element_names.DISTANCE_ANALYSIS_IMAGES}/{element_names.DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE}",
                )
                try:
                    analysis_results.structure_aln_image = (
                        structure_aln_images[0].attrib[attribute_names.DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE_BASENAME],
                        structure_aln_images[0].text,
                    )
                except IndexError:
                    analysis_results.structure_aln_image = ()
                interesting_reg_images = tmp_protein_pair.findall(
                    f".//{element_names.DISTANCE_ANALYSIS_IMAGES}/{element_names.DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS}",
                )
                for tmp_image in interesting_reg_images:
                    analysis_results.interesting_regions_images.append(
                        (
                            tmp_image.attrib[attribute_names.DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS_BASENAME],
                            tmp_image.text,
                        ),
                    )


class ObjectSerializer:
    """Class which handles serialization of objects."""

    object_dict: dict

    def __init__(self, object_to_serialize, filepath: pathlib.Path, filename: str) -> None:  # noqa: ANN001
        """Constructor.

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

    def create_standard_object_dict(self) -> None:
        """Creates standard object dictionary for serialization."""
        self.object_dict = self.object.__dict__

    def set_custom_object_dict(self, custom_dict: dict) -> None:
        """Sets a custom dictionary for serialization."""
        self.object_dict = custom_dict

    def update_object_dict(self, dict_with_updates: dict) -> None:
        """Updates a dictionary for serialization."""
        self.object_dict.update(dict_with_updates)

    def serialize_object(self) -> None:
        """Serializes the object to a JSON file."""
        tmp_object_file = open(pathlib.Path(f"{self.filepath}/{self.filename}.json"), "w", encoding="utf-8")
        json.dump(self.object_dict, tmp_object_file, indent=4)
        tmp_object_file.close()


class ObjectDeserializer:
    """Class which deserializes objects."""

    def __init__(self, filepath: pathlib.Path, filename: str) -> None:
        """Constructor.

        Args:
            filepath: filepath of the json file, without the file name
            filename: name of the json file WITHOUT .json extension!
        """
        tmp_object_file = open(pathlib.Path(f"{filepath}/{filename}.json"), "r", encoding="utf-8")
        self.object_dict = json.load(tmp_object_file)

    def deserialize_protein(self) -> "protein.Protein":
        """Deserializes the protein object from a JSON file."""
        if self.object_dict.get("export_data_dir") == "None":
            update = {"export_data_dir": None}
            self.object_dict.update(update)
        tmp_protein = protein.Protein(
            molecule_object=self.object_dict.get("filename"),  # important for duplicated proteins
            proteins_dirname=self.object_dict.get("proteins_dirname"),
            pdb_filepath=path_util.FilePath(self.object_dict.get("import_data_dir")),
        )
        tmp_protein.molecule_object = self.object_dict.get("molecule_object")
        tmp_protein.pymol_selection.selection_string = self.object_dict.get("selection")
        tmp_protein.chains = protein_util.create_chains_from_list_of_tuples(self.object_dict.get("chains"))
        tmp_protein.pymol_session_file = self.object_dict.get("pymol_session_file")
        return tmp_protein

    def deserialize_protein_pair(self) -> "protein_pair.ProteinPair":
        """Deserialize a protein pair object from a JSON file."""
        tmp_protein_pair = protein_pair.ProteinPair(
            self.object_dict.get("prot_1_molecule_object"),
            self.object_dict.get("prot_2_molecule_object"),
            pathlib.Path(self.object_dict.get("export_dirname")),
        )
        tmp_protein_pair.pymol_session_file = self.object_dict.get("pymol_session_file")
        tmp_protein_pair.name = self.object_dict.get("name")

        if self.object_dict.get("prot_1_export_data_dir") == "None":
            export_data_dir = None
        else:
            export_data_dir = self.object_dict.get("prot_1_export_data_dir")

        tmp_protein_1 = protein.Protein(
            molecule_object=self.object_dict.get("prot_1_filename"),  # important for duplicated proteins
            proteins_dirname=self.object_dict.get("proteins_dirname"),
            pdb_filepath=self.object_dict.get("prot_1_import_data_dir"),
        )
        tmp_protein_1.molecule_object = self.object_dict.get("prot_1_molecule_object")
        tmp_protein_1.pymol_selection.selection_string = self.object_dict.get("prot_1_selection")
        tmp_protein_1.chains = protein_util.create_chains_from_list_of_tuples(self.object_dict.get("prot_1_chains"))
        tmp_protein_pair.ref_obj = tmp_protein_1

        if self.object_dict.get("prot_2_export_data_dir") == "None":
            export_data_dir = None
        else:
            export_data_dir = self.object_dict.get("prot_2_export_data_dir")  # noqa: F841
        tmp_protein_2 = protein.Protein(
            molecule_object=self.object_dict.get("prot_2_filename"),
            # important for duplicated proteins
            proteins_dirname=self.object_dict.get("proteins_dirname"),
            pdb_filepath=self.object_dict.get("prot_2_import_data_dir"),
        )
        tmp_protein_2.molecule_object = self.object_dict.get("prot_2_molecule_object")
        tmp_protein_2.pymol_selection.selection_string = self.object_dict.get("prot_2_selection")
        tmp_protein_2.chains = protein_util.create_chains_from_list_of_tuples(self.object_dict.get("prot_2_chains"))
        tmp_protein_pair.model_obj = tmp_protein_2

        return tmp_protein_pair

    def deserialize_project(self, app_settings: "settings.Settings") -> "project.Project":
        """Deserializes the project object from the JSON file."""
        tmp_project: "project.Project" = project.Project(
            self.object_dict.get("_project_name"),
            self.object_dict.get("_workspace"),
        )
        tmp_project.set_folder_paths(self.object_dict.get("_folder_paths"))
        tmp_project.project_path = self.object_dict.get("project_path")
        for tmp_protein_json_filepath in project_util.get_all_protein_json_filepaths_from_project(tmp_project):
            tmp_project.add_existing_protein(
                ObjectDeserializer(
                    tmp_protein_json_filepath.get_dirname(),
                    tmp_protein_json_filepath.get_filename(),
                ).deserialize_protein(),
            )
        for tmp_protein_pair_json_filepath in project_util.get_all_protein_pair_json_filepaths_from_project(
            tmp_project,
        ):
            tmp_project.add_protein_pair(
                ObjectDeserializer(
                    tmp_protein_pair_json_filepath.get_dirname(),
                    tmp_protein_pair_json_filepath.get_filename(),
                ).deserialize_protein_pair(),
            )
        for tmp_distance_analysis_json_filepath in project_util.get_all_distance_analysis_json_filepaths_from_project(
            tmp_project,
        ):
            tmp_project.add_distance_analysis(
                ObjectDeserializer(
                    tmp_distance_analysis_json_filepath.get_dirname(),
                    tmp_distance_analysis_json_filepath.get_filename(),
                ).deserialize_distance_analysis(app_settings=app_settings, app_project=tmp_project),
            )
        return tmp_project

    def deserialize_distance_analysis(
        self,
        app_settings: "settings.Settings",
        app_project: "project.Project",
    ) -> "distance_analysis.DistanceAnalysis":
        """Deserializes distance analysis object from a JSON file."""
        tmp_protein_pair = app_project.search_protein_pair(self.object_dict.get("protein_pair_for_analysis_name"))
        tmp_distance_analysis = distance_analysis.DistanceAnalysis(tmp_protein_pair, app_settings)
        tmp_distance_analysis.cutoff = self.object_dict.get("cutoff")
        tmp_distance_analysis.cycles = self.object_dict.get("cycles")
        tmp_distance_analysis.export_dirname = self.object_dict.get("export_dirname")
        tmp_distance_analysis.pymol_session_filepath = self.object_dict.get("pymol_session_filepath")
        return tmp_distance_analysis

    def deserialize_settings(self) -> "settings.Settings":
        """Deserializes settings object from a JSON file."""
        tmp_settings: settings.Settings = settings.Settings(
            self.object_dict.get("dir_settings"),
            self.object_dict.get("filename"),
        )
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
        if self.object_dict.get("filename") != "settings.json":
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

    def deserialize_sequence(self) -> "sequence.Sequence":
        """Deserializes sequence object from a JSON file."""
        return sequence.Sequence(self.object_dict.get("name"), self.object_dict.get("sequence"))


class WorkspaceScanner:
    """Class to scan the workspace."""

    def __init__(self, workspace_path: pathlib.Path) -> None:
        """Constructor.

        Args:
            workspace_path: the path to the workspace.
        """
        self.workspace_path = workspace_path

    def scan_workspace_for_valid_projects(self) -> list[str]:
        """Scans the workspace for valid projects."""
        directory_content = os.listdir(str(self.workspace_path))
        return [file for file in directory_content if file.endswith(".xml")]

    def scan_workspace_for_non_duplicate_proteins(self) -> np.ndarray:
        """Scans the workspace for non-duplicate proteins."""
        tmp_workspace_path: str = str(self.workspace_path)
        # List of XML file paths
        xml_files = [
            f"{tmp_workspace_path}/{tmp_project_file}"
            for tmp_project_file in os.listdir(tmp_workspace_path)
            if not os.path.isdir(pathlib.Path(f"{tmp_workspace_path}/{tmp_project_file}"))
        ]
        # Using ThreadPoolExecutor to read XML files in parallel
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Create a queue and a lock for thread safety
            result_queue = Queue()
            lock = Lock()

            # Map each XML file to the read_xml_file function in parallel
            futures = [
                executor.submit(self._get_proteins_of_project_xml, pathlib.Path(file_path), result_queue, lock)
                for file_path in xml_files
            ]

            # Wait for all tasks to complete
            concurrent.futures.wait(futures)

            # Retrieve results from the queue in the main thread
            all_proteins_of_workspace: collections.deque = collections.deque()
            while not result_queue.empty():
                result = result_queue.get()
                protein_name, protein_info = result
                all_proteins_of_workspace.append(result[1])
                # print(f"Name: {protein_name}, Info: {protein_info}")
            flatten_proteins = list(itertools.chain(*all_proteins_of_workspace))
            print(flatten_proteins)
            # Use a set to store unique elements
            unique_proteins_set = set()
            # Iterate through the list and add unique elements to the set
            for tmp_protein in flatten_proteins:
                unique_proteins_set.add(tmp_protein)

            return np.array(list(unique_proteins_set))

    def _get_proteins_of_project_xml(self, xml_filepath: pathlib.Path, result_queue, lock) -> None:  # noqa: ANN001
        """Gets the protein information from the project xml file.

        Args:
            xml_filepath: a path to the project xml file.
            result_queue: a queue to put the protein information into.
            lock: a lock to prevent the unauthorized access to the memory.
        """
        tmp_protein_names: collections.deque = collections.deque()
        tmp_protein_infos: collections.deque = collections.deque()
        xml_deserializer = XmlDeserializer(xml_filepath)
        project_name = str(xml_filepath.name).replace(".xml", "")
        for tmp_protein in xml_deserializer.xml_root.iter(element_names.PROTEIN):
            molecule_object = tmp_protein.attrib[attribute_names.PROTEIN_MOLECULE_OBJECT]
            if molecule_object not in tmp_protein_names:
                tmp_protein_names.append(molecule_object)
                tmp_protein_infos.append(
                    basic_protein_info.BasicProteinInfo(
                        molecule_object,
                        tmp_protein.attrib[attribute_names.ID],
                        project_name,
                    ),
                )
        with lock:
            result_queue.put((list(tmp_protein_names), list(tmp_protein_infos)))


class FilesystemCleaner:
    """Class for cleaning up the filesystem."""

    def __int__(self) -> None:
        """Empty constructor."""
        pass

    @staticmethod
    def clean_prediction_scratch_folder() -> None:
        """Deletes the scratch folder for fasta and pdb files and creates new ones."""
        shutil.rmtree(constants.PREDICTION_FASTA_DIR)
        os.mkdir(constants.PREDICTION_FASTA_DIR)
        shutil.rmtree(constants.PREDICTION_PDB_DIR)
        os.mkdir(constants.PREDICTION_PDB_DIR)
