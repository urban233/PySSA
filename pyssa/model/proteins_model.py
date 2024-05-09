#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module contains the proteins model."""
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from auxiliary_pymol import auxiliary_pymol_client
from pyssa.internal.data_structures import protein, job
from pyssa.util import enums

__docformat__ = "google"


class ProteinsModel(QtGui.QStandardItemModel):

    def __init__(self):
        super().__init__()

    def build_model_from_scratch(self, the_protein_objects: list["protein.Protein"], the_main_socket, a_socket):
        """
        Builds a model from scratch using the given protein objects.

        Args:
            the_protein_objects (list["protein.Protein"]): A list of protein objects.
            the_main_socket: The main socket used for communication.
            a_socket: A socket used for communication.

        Example:
            protein_objects = [protein1, protein2, protein3]
            main_socket = "main_socket"
            a_socket = "a_socket"
            build_model_from_scratch(protein_objects, main_socket, a_socket)

        Note:
            This method builds a model by creating items for each protein, scenes, and chains. It retrieves the scenes of
            each protein using the given sockets and adds them as child items to the corresponding protein item. It also adds
            chains as child items to the corresponding protein item.

            The model is built using QStandardItem objects.
        """
        # TODO: multiprocessing code does not work (gets blocked if a prediction is running) and it is not really faster
        # pool_information = []
        # for tmp_protein in the_protein_objects:
        #     pool_information.append((tmp_protein.pymol_session, tmp_protein.get_molecule_object()))
        # # Use a context manager to create a multiprocessing pool
        # with multiprocessing.Pool(multiprocessing.cpu_count() - 2) as pool:
        #     # Map the function to the pool to run in parallel
        #     tmp_all_scenes_from_all_sessions = pool.starmap(auxiliary_pymol.AuxiliaryPyMOL.get_all_scenes_of_session_multi, pool_information)

        tmp_root_item = self.invisibleRootItem()
        i = 0
        for tmp_protein in the_protein_objects:
            tmp_job_description = job.GeneralPurposeJobDescription(enums.JobShortDescription.GET_ALL_SCENES_OF_SESSION)
            tmp_job_description.setup_dict(
                {enums.JobDescriptionKeys.PYMOL_SESSION.value: str(tmp_protein.pymol_session)}
            )
            tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
                the_main_socket, a_socket, tmp_job_description
            )
            tmp_all_scenes = tmp_reply["data"]
            # protein node (type = protein)
            tmp_protein_item = QtGui.QStandardItem(tmp_protein.get_molecule_object())
            tmp_protein_item.setData(tmp_protein, enums.ModelEnum.OBJECT_ROLE)
            tmp_protein_item.setData("protein", enums.ModelEnum.TYPE_ROLE)
            tmp_root_item.appendRow(tmp_protein_item)
            # Scenes node (type = header)
            tmp_scenes_item = QtGui.QStandardItem("Scenes")
            tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_item.appendRow(tmp_scenes_item)
            # scene nodes (type = scene)
            for tmp_scene in tmp_all_scenes:
                tmp_scene_item = QtGui.QStandardItem(tmp_scene)
                tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
                tmp_scenes_item.appendRow(tmp_scene_item)
            # Chains node (type = header)
            tmp_chains_item = QtGui.QStandardItem("Chains")
            tmp_chains_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_item.appendRow(tmp_chains_item)
            # chain nodes (type = chain)
            for tmp_chain in tmp_protein.chains:
                tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
                tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
                tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
                tmp_chain_item.setData(
                    tmp_chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value], enums.ModelEnum.CHAIN_COLOR_ROLE
                )
                tmp_chains_item.appendRow(tmp_chain_item)
            i += 1

    def check_if_scratch_scene_exists(self, a_model_index: QtCore.QModelIndex) -> bool:
        """
        Check if a scratch scene exists for the given model index.

        Args:
            a_model_index (QtCore.QModelIndex): The model index.

        Returns:
            bool: True if a scratch scene exists, False otherwise.

        Raises:
            ValueError: If the model index has an invalid type.
        """
        tmp_type = a_model_index.data(enums.ModelEnum.TYPE_ROLE)
        # on protein node
        if tmp_type == "protein":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).child(0, 0)
        # on header node (Scenes)
        elif tmp_type == "header" and a_model_index.data(Qt.DisplayRole) == "Scenes":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index)
        # on scene node
        elif tmp_type == "scene":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent()
        # on header node (Chains)
        elif tmp_type == "header" and a_model_index.data(Qt.DisplayRole) == "Chains":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().child(0, 0)
        # on chain node
        elif tmp_type == "chain":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().parent().child(0, 0)
        else:
            raise ValueError("Wrong type!")
        for i in range(tmp_scenes_header_item.rowCount()):
            if tmp_scenes_header_item.child(i, 0).data(Qt.DisplayRole) == "_scratch_":
                return True
            i += 1
        return False

    def add_scene(self, a_model_index: QtCore.QModelIndex, the_scene_item_to_add: QtGui.QStandardItem) -> None:
        """
        Add a scene item to the tree model.

        Args:
            a_model_index: The index of the tree model where the scene item should be added.
            the_scene_item_to_add: The scene item to be added.

        Returns:
            None

        Raises:
            ValueError: If the type of the model index is invalid.
        """
        tmp_type = a_model_index.data(enums.ModelEnum.TYPE_ROLE)
        # on protein node
        if tmp_type == "protein":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).child(0, 0)
        # on header node (Scenes)
        elif tmp_type == "header" and a_model_index.data(Qt.DisplayRole) == "Scenes":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index)
        # on scene node
        elif tmp_type == "scene":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent()
        # on header node (Chains)
        elif tmp_type == "header" and a_model_index.data(Qt.DisplayRole) == "Chains":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().child(0, 0)
        # on chain node
        elif tmp_type == "chain":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().parent().child(0, 0)
        else:
            raise ValueError("Wrong type!")
        tmp_scenes_header_item.appendRow(the_scene_item_to_add)

    def remove_scene(self, the_model_index_of_the_scene: QtCore.QModelIndex) -> None:
        """
        Removes a scene from the model.

        Args:
            the_model_index_of_the_scene: The index of the scene to be removed.

        Returns:
            None

        Raises:
            ValueError: If the type of the item at the given index is not "scene".

        """
        # <editor-fold desc="Checks">
        if self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE) != "scene":
            raise ValueError("Wrong type!")

        # </editor-fold>

        tmp_scene_item = self.itemFromIndex(the_model_index_of_the_scene)
        tmp_scenes_item = tmp_scene_item.parent()
        tmp_scenes_item.removeRow(tmp_scene_item.row())

    def add_protein(self, a_protein: "protein.Protein") -> None:
        """
        Add a protein to the data model.

        Args:
            a_protein (protein.Protein): The protein to add.

        Returns:
            None
        """
        tmp_protein_item = QtGui.QStandardItem(a_protein.get_molecule_object())
        tmp_protein_item.setData("protein", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_item.setData(a_protein, enums.ModelEnum.OBJECT_ROLE)
        self.appendRow(tmp_protein_item)
        # Scenes node (type = header)
        tmp_scenes_item = QtGui.QStandardItem("Scenes")
        tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_item.appendRow(tmp_scenes_item)
        # scene nodes (type = scene)
        # fixme: I don't know if only the base scene can be used
        # a_protein.load_protein_pymol_session()
        # for tmp_scene in pymol_io.get_all_scenes_from_pymol_session():
        #     tmp_scene_item = QtGui.QStandardItem(tmp_scene)
        #     tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
        #     tmp_scenes_item.appendRow(tmp_scene_item)
        tmp_scene_item = QtGui.QStandardItem("base")
        tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
        tmp_scenes_item.appendRow(tmp_scene_item)
        # Chains node (type = header)
        tmp_chains_item = QtGui.QStandardItem("Chains")
        tmp_chains_item.setData("header", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_item.appendRow(tmp_chains_item)
        # chain nodes (type = chain)
        for tmp_chain in a_protein.chains:
            tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
            tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
            tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
            tmp_chains_item.appendRow(tmp_chain_item)

    def remove_protein(self, the_model_index_of_the_scene: QtCore.QModelIndex):
        """
        Removes a protein from the scene.

        Args:
            the_model_index_of_the_scene (QtCore.QModelIndex): The index of the protein in the scene.

        Raises:
            ValueError: If the type of the item at the given index is not "protein".
        """
        # <editor-fold desc="Checks">
        if self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE) != "protein":
            raise ValueError("Wrong type!")

        # </editor-fold>

        tmp_protein_item = self.itemFromIndex(the_model_index_of_the_scene)
        self.removeRow(tmp_protein_item.row())


class TemporaryProteinsModel(ProteinsModel):

    def __init__(self):
        super().__init__()

    def build_model_from_scratch(self, the_protein_objects: list["protein.Protein"]):
        """
        Builds a model from scratch using the given protein objects.

        Args:
            the_protein_objects (list[protein.Protein]): A list of protein objects.

        """
        tmp_root_item = self.invisibleRootItem()
        for tmp_protein in the_protein_objects:
            # protein node (type = protein)
            tmp_protein_item = QtGui.QStandardItem(tmp_protein.get_molecule_object())
            tmp_protein_item.setData("protein", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_item.setData(tmp_protein, enums.ModelEnum.OBJECT_ROLE)
            self.appendRow(tmp_protein_item)
            # Scenes node (type = header)
            tmp_scenes_item = QtGui.QStandardItem("Scenes")
            tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_item.appendRow(tmp_scenes_item)
            # scene nodes (type = scene)
            # a_protein.load_protein_pymol_session()
            # for tmp_scene in pymol_io.get_all_scenes_from_pymol_session():
            tmp_scene_item = QtGui.QStandardItem("generic")
            tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
            tmp_scenes_item.appendRow(tmp_scene_item)
            # Chains node (type = header)
            tmp_chains_item = QtGui.QStandardItem("Chains")
            tmp_chains_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_item.appendRow(tmp_chains_item)
            # chain nodes (type = chain)
            for tmp_chain in tmp_protein.chains:
                tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
                tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
                tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
                tmp_chains_item.appendRow(tmp_chain_item)

    def add_temporary_protein(self, a_protein: "protein.Protein"):
        """
        Adds a temporary protein to the model.

        Args:
            a_protein: An object of type Protein.

        Returns:
            None
        """
        tmp_protein_item = QtGui.QStandardItem(a_protein.get_molecule_object())
        tmp_protein_item.setData("protein", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_item.setData(a_protein, enums.ModelEnum.OBJECT_ROLE)
        self.appendRow(tmp_protein_item)
        # Scenes node (type = header)
        tmp_scenes_item = QtGui.QStandardItem("Scenes")
        tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_item.appendRow(tmp_scenes_item)
        # scene nodes (type = scene)
        # a_protein.load_protein_pymol_session()
        # for tmp_scene in pymol_io.get_all_scenes_from_pymol_session():
        tmp_scene_item = QtGui.QStandardItem("generic")
        tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
        tmp_scenes_item.appendRow(tmp_scene_item)
        # Chains node (type = header)
        tmp_chains_item = QtGui.QStandardItem("Chains")
        tmp_chains_item.setData("header", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_item.appendRow(tmp_chains_item)
        # chain nodes (type = chain)
        for tmp_chain in a_protein.chains:
            tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
            tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
            tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
            tmp_chains_item.appendRow(tmp_chain_item)
