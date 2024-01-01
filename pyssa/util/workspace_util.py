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
import concurrent.futures
import collections
import itertools
from queue import Queue
from threading import Lock

import numpy as np

from pyssa.internal.data_structures.data_classes import basic_protein_info
from pyssa.io_pyssa import filesystem_io
from pyssa.io_pyssa.xml_pyssa import element_names, attribute_names


def get_proteins_of_project_xml(xml_filepath: pathlib.Path, result_queue, lock):
    tmp_protein_names: collections.deque = collections.deque()
    tmp_protein_infos: collections.deque = collections.deque()
    xml_deserializer = filesystem_io.XmlDeserializer(xml_filepath)
    project_name = str(xml_filepath.name).replace(".xml", "")
    for tmp_protein in xml_deserializer.xml_root.iter(element_names.PROTEIN):
        molecule_object = tmp_protein.attrib[attribute_names.PROTEIN_MOLECULE_OBJECT]
        if molecule_object not in tmp_protein_names:
            tmp_protein_names.append(molecule_object)
            tmp_protein_infos.append(
                basic_protein_info.BasicProteinInfo(
                    molecule_object, tmp_protein.attrib[attribute_names.ID], project_name
                )
            )
    with lock:
        result_queue.put((list(tmp_protein_names), list(tmp_protein_infos)))


def scan_workspace_for_non_duplicate_proteins(the_workspace_path) -> np.ndarray:
    # List of XML file paths
    xml_files = [
        f"{the_workspace_path}/{tmp_project_file}"
        for tmp_project_file in os.listdir(the_workspace_path)
        if not os.path.isdir(pathlib.Path(f"{the_workspace_path}/{tmp_project_file}"))
    ]
    # Using ThreadPoolExecutor to read XML files in parallel
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Create a queue and a lock for thread safety
        result_queue = Queue()
        lock = Lock()

        # Map each XML file to the read_xml_file function in parallel
        futures = [
            executor.submit(get_proteins_of_project_xml, pathlib.Path(file_path), result_queue, lock)
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
        for protein in flatten_proteins:
            unique_proteins_set.add(protein)

        return np.array(list(unique_proteins_set))


def scan_workspace_for_valid_projects(the_workspace_path):
    directory_content = os.listdir(the_workspace_path)
    # directory_content.sort()
    xml_files = [file for file in directory_content if file.endswith(".xml")]
    return xml_files
