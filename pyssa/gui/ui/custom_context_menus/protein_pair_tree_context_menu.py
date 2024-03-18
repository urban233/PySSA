from PyQt5 import QtWidgets


class ProteinPairTreeContextMenu:
    """A custom context menu wrapper for the protein pair tree view of the main dialog."""

    def __init__(self):
        self._context_menu = QtWidgets.QMenu()
        self._open_results_summary_action = QtWidgets.QAction("Open Results Summary")
        self._color_based_on_rmsd_action = QtWidgets.QAction("Color Based On RMSD")
        self._help_action = QtWidgets.QAction("Get Help")
        
        self._context_menu.addAction(self._open_results_summary_action)
        self._context_menu.addAction(self._color_based_on_rmsd_action)
        self._context_menu.addAction(self._help_action)

        self._open_results_summary_action.triggered.connect(self._generic_action)
        self._color_based_on_rmsd_action.triggered.connect(self._generic_action)
        self._help_action.triggered.connect(self._generic_action)

    # <editor-fold desc="Connect actions from outside">
    def connect_open_results_summary_action(self, the_function_to_connect):
        self._open_results_summary_action.triggered.disconnect()
        self._open_results_summary_action.triggered.connect(the_function_to_connect)
    
    def connect_color_based_on_rmsd_action(self, the_function_to_connect):
        self._color_based_on_rmsd_action.triggered.disconnect()
        self._color_based_on_rmsd_action.triggered.connect(the_function_to_connect)

    def connect_help_action(self, the_function_to_connect):
        self._help_action.triggered.disconnect()
        self._help_action.triggered.connect(the_function_to_connect)
    # </editor-fold>
    
    def _generic_action(self):
        raise NotImplementedError("The called action is not connected to a function!")
    
    def get_context_menu(self,
                         the_selected_indexes,
                         is_protein_pair_in_current_session_flag: bool):
        # <editor-fold desc="Checks">
        if len(the_selected_indexes) > 0:
            level = 0
            index = the_selected_indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1
        else:
            # No protein pair is selected
            self._help_action.setVisible(True)
            self._open_results_summary_action.setVisible(False)
            self._color_based_on_rmsd_action.setVisible(False)
            return self._context_menu
        # </editor-fold>

        self._help_action.setVisible(False)
        if level == 0:
            # protein pair level
            self._open_results_summary_action.setVisible(True)
            self._color_based_on_rmsd_action.setVisible(True)

            if is_protein_pair_in_current_session_flag:
                self._color_based_on_rmsd_action.setEnabled(True)
            else:
                self._color_based_on_rmsd_action.setEnabled(False)
        else:
            self._open_results_summary_action.setVisible(False)
            self._color_based_on_rmsd_action.setVisible(False)
        return self._context_menu
