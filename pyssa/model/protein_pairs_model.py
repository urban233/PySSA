from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from pyssa.internal.data_structures import protein, protein_pair
from pyssa.internal.data_structures.protein import Protein
from pyssa.internal.portal import pymol_io
from pyssa.util import enums


class ProteinPairsModel(QtGui.QStandardItemModel):

    def __init__(self):
        super().__init__()

    def build_model_from_scratch(self, the_protein_pair_objects: list["protein_pair.ProteinPair"]):
        """Builds the model from scratch."""
        tmp_root_item = self.invisibleRootItem()
        for tmp_protein_pair in the_protein_pair_objects:
            tmp_protein_pair_item = QtGui.QStandardItem(tmp_protein_pair.name)
            tmp_protein_pair_item.setData(tmp_protein_pair, enums.ModelEnum.OBJECT_ROLE)
            tmp_protein_pair_item.setData("protein_pair", enums.ModelEnum.TYPE_ROLE)
            tmp_scenes_item = QtGui.QStandardItem("Scenes")
            tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            tmp_protein_pair_item.appendRow(tmp_scenes_item)
            tmp_protein_pair.load_pymol_session()
            for tmp_scene in pymol_io.get_all_scenes_from_pymol_session():
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

    def add_scene(
            self,
            a_model_index: QtCore.QModelIndex,
            the_scene_item_to_add: QtGui.QStandardItem
    ) -> None:
        """Adds a scene to the model."""
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
        """Removes a scene from the model."""
        # <editor-fold desc="Checks">
        if self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE) != "scene":
            raise ValueError("Wrong type!")

        # </editor-fold>

        tmp_scene_item = self.itemFromIndex(the_model_index_of_the_scene)
        tmp_scenes_item = tmp_scene_item.parent()
        tmp_scenes_item.removeRow(tmp_scene_item.row())

    def add_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair"):
        """Adds a protein pair to the model."""
        tmp_protein_pair_item = QtGui.QStandardItem(a_protein_pair.name)
        tmp_protein_pair_item.setData(a_protein_pair, enums.ModelEnum.OBJECT_ROLE)
        tmp_protein_pair_item.setData("protein_pair", enums.ModelEnum.TYPE_ROLE)
        self.appendRow(tmp_protein_pair_item)
        tmp_scenes_item = QtGui.QStandardItem("Scenes")
        tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
        tmp_protein_pair_item.appendRow(tmp_scenes_item)
        a_protein_pair.load_pymol_session()
        for tmp_scene in pymol_io.get_all_scenes_from_pymol_session():
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

    def remove_protein_pair(self, the_model_index_of_the_scene: QtCore.QModelIndex):
        """Removes a protein from the model."""
        # <editor-fold desc="Checks">
        if self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE) != "protein_pair":
            tmp_value = self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE)
            raise ValueError("Wrong type!")

        # </editor-fold>

        tmp_protein_pair_item = self.itemFromIndex(the_model_index_of_the_scene)
        self.removeRow(tmp_protein_pair_item.row())
