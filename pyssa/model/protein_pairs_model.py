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
"""Module contains the protein pairs model."""
import logging

import zmq
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from auxiliary_pymol import auxiliary_pymol_client
from pyssa.internal.data_structures import protein_pair, job
from pyssa.logging_pyssa import log_handlers
from pyssa.util import enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class ProteinPairsModel(QtGui.QStandardItemModel):
    """Contains the protein pairs of the project in form of a QStandardItemModel."""
    
    def __init__(self) -> None:
        """Constructor."""
        super().__init__()

    def build_model_from_scratch(self,
                                 the_protein_pair_objects: list["protein_pair.ProteinPair"],
                                 the_main_socket: zmq.Socket,
                                 a_socket: zmq.Socket) -> None:
        """Builds a model from scratch using the given protein objects.

        Args:
            the_protein_pair_objects (list[protein_pair.ProteinPair]): A list of protein objects.
            the_main_socket (zmq.Socket): The main socket used for communication.
            a_socket (zmq.Socket): A socket used for communication.
        
        Raises:
            exception.IllegalArgumentError: If any of the arguments are None.
        """
        # <editor-fold desc="Checks">
        if the_protein_pair_objects is None:
            logger.error("the_protein_pair_objects is None.")
            raise exception.IllegalArgumentError("the_protein_pair_objects is None.")
        if the_main_socket is None:
            logger.error("the_main_socket is None.")
            raise exception.IllegalArgumentError("the_main_socket is None.")
        if a_socket is None:
            logger.error("a_socket is None.")
            raise exception.IllegalArgumentError("a_socket is None.")

        # </editor-fold>
        
        # TODO: multiprocessing code does not work (gets blocked if a prediction is running) and it is not really faster
        # pool_information = []
        # for tmp_protein_pair in the_protein_pair_objects:
        #     pool_information.append((tmp_protein_pair.pymol_session, tmp_protein_pair.name))
        # # Use a context manager to create a multiprocessing pool
        # with multiprocessing.Pool(processes=multiprocessing.cpu_count() - 2) as pool:
        #     # Map the function to the pool to run in parallel
        #     tmp_all_scenes_from_all_sessions = pool.starmap(
        #         auxiliary_pymol.AuxiliaryPyMOL.get_all_scenes_of_session_multi, pool_information)

        tmp_root_item = self.invisibleRootItem()
        for tmp_protein_pair in the_protein_pair_objects:
            tmp_job_description = job.GeneralPurposeJobDescription(
                enums.JobShortDescription.GET_ALL_SCENES_OF_SESSION,
            )
            tmp_job_description.setup_dict(
                {enums.JobDescriptionKeys.PYMOL_SESSION.value: str(tmp_protein_pair.pymol_session)})
            tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
                the_main_socket, a_socket, tmp_job_description,
            )
            tmp_all_scenes = tmp_reply["data"]

            tmp_protein_pair_item = QtGui.QStandardItem(tmp_protein_pair.name)
            tmp_protein_pair_item.setData(tmp_protein_pair, enums.ModelEnum.OBJECT_ROLE)
            tmp_protein_pair_item.setData("protein_pair", enums.ModelEnum.TYPE_ROLE)
            tmp_scenes_item = QtGui.QStandardItem("Scenes")
            tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_pair_item.appendRow(tmp_scenes_item)
            for tmp_scene in tmp_all_scenes:
                tmp_scene_item = QtGui.QStandardItem(tmp_scene)
                tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
                tmp_scenes_item.appendRow(tmp_scene_item)
            # Create protein 1 item
            tmp_protein_item_1 = QtGui.QStandardItem(tmp_protein_pair.protein_1.get_molecule_object())
            tmp_protein_item_1.setData(tmp_protein_pair.protein_1, enums.ModelEnum.OBJECT_ROLE)
            tmp_protein_item_1.setData("protein", enums.ModelEnum.TYPE_ROLE)
            tmp_chains_item_1 = QtGui.QStandardItem("Chains")
            tmp_chains_item_1.setData("header", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_item_1.appendRow(tmp_chains_item_1)
            for tmp_chain in tmp_protein_pair.protein_1.chains:
                tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
                tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
                tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
                tmp_chains_item_1.appendRow(tmp_chain_item)
            # Create protein 2 item
            tmp_protein_item_2 = QtGui.QStandardItem(tmp_protein_pair.protein_2.get_molecule_object())
            tmp_protein_item_2.setData(tmp_protein_pair.protein_2, enums.ModelEnum.OBJECT_ROLE)
            tmp_protein_item_2.setData("protein", enums.ModelEnum.TYPE_ROLE)
            tmp_chains_item_2 = QtGui.QStandardItem("Chains")
            tmp_chains_item_2.setData("header", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_item_2.appendRow(tmp_chains_item_2)
            for tmp_chain in tmp_protein_pair.protein_2.chains:
                tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
                tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
                tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
                tmp_chains_item_2.appendRow(tmp_chain_item)
            tmp_protein_pair_item.appendRow(tmp_protein_item_1)
            tmp_protein_pair_item.appendRow(tmp_protein_item_2)
            tmp_root_item.appendRow(tmp_protein_pair_item)

    def check_if_scratch_scene_exists(self,
                                      a_model_index: QtCore.QModelIndex) -> bool:
        """Check if a scratch scene exists for the given model index.

        Args:
            a_model_index (QtCore.QModelIndex): The model index.

        Returns:
            bool: True if a scratch scene exists, False otherwise.

        Raises:
            exception.IllegalArgumentError: If a_model_index is None.
            ValueError: If the model index has an invalid type.
        """
        # <editor-fold desc="Checks">
        if a_model_index is None:
            logger.error("a_model_index is None.")
            raise exception.IllegalArgumentError("a_model_index is None.")
        
        # </editor-fold>

        tmp_type = a_model_index.data(enums.ModelEnum.TYPE_ROLE)
        # on protein pair node
        if tmp_type == "protein_pair":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).child(0, 0)
        # on protein node
        elif tmp_type == "protein":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().child(0, 0)
        # on header node (Scenes)
        elif tmp_type == "header" and a_model_index.data(Qt.DisplayRole) == "Scenes":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index)
        # on scene node
        elif tmp_type == "scene":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent()
        # on header node (Chains)
        elif tmp_type == "header" and a_model_index.data(Qt.DisplayRole) == "Chains":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().parent().child(0, 0)
        # on chain node
        elif tmp_type == "chain":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().parent().parent().child(0, 0)
        else:
            raise ValueError("Wrong type!")
        for i in range(tmp_scenes_header_item.rowCount()):
            if tmp_scenes_header_item.child(i, 0).data(Qt.DisplayRole) == "_scratch_":
                return True
            i += 1
        return False

    def add_scene(
            self,
            a_model_index: QtCore.QModelIndex,
            the_scene_item_to_add: QtGui.QStandardItem,
    ) -> None:
        """Add a scene item to the tree model.

        Args:
            a_model_index (QtCore.QModelIndex): The index of the tree model where the scene item should be added.
            the_scene_item_to_add (QtGui.QStandardItem): The scene item to be added.

        Raises:
            exception.IllegalArgumentError: If any of the arguments are None.
            ValueError: If the type of the model index is invalid.
        """
        # <editor-fold desc="Checks">
        if a_model_index is None:
            logger.error("a_model_index is None.")
            raise exception.IllegalArgumentError("a_model_index is None.")
        if the_scene_item_to_add is None:
            logger.error("the_scene_item_to_add is None.")
            raise exception.IllegalArgumentError("the_scene_item_to_add is None.")
        
        # </editor-fold>
        tmp_type = a_model_index.data(enums.ModelEnum.TYPE_ROLE)
        # on protein pair node
        if tmp_type == "protein_pair":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).child(0, 0)
        # on protein node
        elif tmp_type == "protein":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().child(0, 0)
        # on header node (Scenes)
        elif tmp_type == "header" and a_model_index.data(Qt.DisplayRole) == "Scenes":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index)
        # on scene node
        elif tmp_type == "scene":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent()
        # on header node (Chains)
        elif tmp_type == "header" and a_model_index.data(Qt.DisplayRole) == "Chains":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().parent().child(0, 0)
        # on chain node
        elif tmp_type == "chain":
            tmp_scenes_header_item = self.itemFromIndex(a_model_index).parent().parent().parent().child(0, 0)
        else:
            raise ValueError("Wrong type!")
        tmp_scenes_header_item.appendRow(the_scene_item_to_add)

    def remove_scene(self, the_model_index_of_the_scene: QtCore.QModelIndex) -> None:
        """Removes a scene from the model.

        Args:
            the_model_index_of_the_scene: The index of the scene to be removed.

        Returns:
            None

        Raises:
            exception.IllegalArgumentError: If the_model_index_of_the_scene is None.
            ValueError: If the type of the item at the given index is not `scene`.
        """
        # <editor-fold desc="Checks">
        if the_model_index_of_the_scene is None:
            logger.error("the_model_index_of_the_scene is None.")
            raise exception.IllegalArgumentError("the_model_index_of_the_scene is None.")
        if self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE) != "scene":
            raise ValueError("Wrong type!")

        # </editor-fold>

        tmp_scene_item = self.itemFromIndex(the_model_index_of_the_scene)
        tmp_scenes_item = tmp_scene_item.parent()
        tmp_scenes_item.removeRow(tmp_scene_item.row())

    def add_protein_pair(self,
                         a_protein_pair: "protein_pair.ProteinPair",
                         the_main_socket: zmq.Socket,
                         a_socket: zmq.Socket) -> None:
        """Add a protein pair to the data model.

        Args:
            a_protein_pair (protein_pair.ProteinPair): The protein pair to add.
            the_main_socket (zmq.Socket): The main socket used for communication.
            a_socket (zmq.Socket): A socket used for communication.
            
        Raises:
            exception.IllegalArgumentError: If any of the arguments are None.
        """
        # <editor-fold desc="Checks">
        if a_protein_pair is None:
            logger.error("a_protein_pair is None.")
            raise exception.IllegalArgumentError("a_protein_pair is None.")
        if the_main_socket is None:
            logger.error("the_main_socket is None.")
            raise exception.IllegalArgumentError("the_main_socket is None.")
        if a_socket is None:
            logger.error("a_socket is None.")
            raise exception.IllegalArgumentError("a_socket is None.")
        
        # </editor-fold>
    
        tmp_job_description = job.GeneralPurposeJobDescription(
            enums.JobShortDescription.GET_ALL_SCENES_OF_SESSION,
        )
        tmp_job_description.setup_dict(
            {enums.JobDescriptionKeys.PYMOL_SESSION.value: str(a_protein_pair.pymol_session)})
        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket, a_socket, tmp_job_description,
        )
        tmp_all_scenes = tmp_reply["data"]

        tmp_protein_pair_item = QtGui.QStandardItem(a_protein_pair.name)
        tmp_protein_pair_item.setData(a_protein_pair, enums.ModelEnum.OBJECT_ROLE)
        tmp_protein_pair_item.setData("protein_pair", enums.ModelEnum.TYPE_ROLE)
        self.appendRow(tmp_protein_pair_item)
        tmp_scenes_item = QtGui.QStandardItem("Scenes")
        tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_pair_item.appendRow(tmp_scenes_item)
        for tmp_scene in tmp_all_scenes:
            tmp_scene_item = QtGui.QStandardItem(tmp_scene)
            tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
            tmp_scenes_item.appendRow(tmp_scene_item)
        # Create protein 1 item
        tmp_protein_item_1 = QtGui.QStandardItem(a_protein_pair.protein_1.get_molecule_object())
        tmp_protein_item_1.setData(a_protein_pair.protein_1, enums.ModelEnum.OBJECT_ROLE)
        tmp_protein_item_1.setData("protein", enums.ModelEnum.TYPE_ROLE)
        tmp_chains_item_1 = QtGui.QStandardItem("Chains")
        tmp_chains_item_1.setData("header", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_item_1.appendRow(tmp_chains_item_1)
        for tmp_chain in a_protein_pair.protein_1.chains:
            tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
            tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
            tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
            tmp_chains_item_1.appendRow(tmp_chain_item)
        # Create protein 2 item
        tmp_protein_item_2 = QtGui.QStandardItem(a_protein_pair.protein_2.get_molecule_object())
        tmp_protein_item_2.setData(a_protein_pair.protein_2, enums.ModelEnum.OBJECT_ROLE)
        tmp_protein_item_2.setData("protein", enums.ModelEnum.TYPE_ROLE)
        tmp_chains_item_2 = QtGui.QStandardItem("Chains")
        tmp_chains_item_2.setData("header", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_item_2.appendRow(tmp_chains_item_2)
        for tmp_chain in a_protein_pair.protein_2.chains:
            tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
            tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
            tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
            tmp_chains_item_2.appendRow(tmp_chain_item)
        tmp_protein_pair_item.appendRow(tmp_protein_item_1)
        tmp_protein_pair_item.appendRow(tmp_protein_item_2)

    def remove_protein_pair(self, the_model_index_of_the_scene: QtCore.QModelIndex) -> None:
        """Removes a protein pair from the model.

        Args:
            the_model_index_of_the_scene (QtCore.QModelIndex): The index of the protein pair in the scene.

        Raises:
            exception.IllegalArgumentError: If the_model_index_of_the_scene is None.
            ValueError: If the type of the item at the given index is not "protein".
        """
        # <editor-fold desc="Checks">
        if the_model_index_of_the_scene is None:
            logger.error("the_model_index_of_the_scene is None.")
            raise exception.IllegalArgumentError("the_model_index_of_the_scene is None.")
        if self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE) != "protein_pair":
            raise ValueError("Wrong type!")

        # </editor-fold>

        tmp_protein_pair_item = self.itemFromIndex(the_model_index_of_the_scene)
        self.removeRow(tmp_protein_pair_item.row())
