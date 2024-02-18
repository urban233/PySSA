from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from pyssa.internal.data_structures import protein
from pyssa.internal.portal import pymol_io
from pyssa.util import enums


class ProteinsModel(QtGui.QStandardItemModel):

    def __init__(self):
        super().__init__()

    def build_model_from_scratch(self, the_protein_objects: list["protein.Protein"]):
        tmp_root_item = self.invisibleRootItem()
        for tmp_protein in the_protein_objects:
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
            tmp_protein.load_protein_pymol_session()
            for tmp_scene in pymol_io.get_all_scenes_from_pymol_session():
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
                tmp_chains_item.appendRow(tmp_chain_item)

    def add_scene(
            self,
            a_model_index: QtCore.QModelIndex,
            the_scene_item_to_add: QtGui.QStandardItem
    ) -> None:
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
        """Removes a scene from the model."""
        # <editor-fold desc="Checks">
        if self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE) != "scene":
            raise ValueError("Wrong type!")

        # </editor-fold>

        tmp_scene_item = self.itemFromIndex(the_model_index_of_the_scene)
        tmp_scenes_item = tmp_scene_item.parent()
        tmp_scenes_item.removeRow(tmp_scene_item.row())
