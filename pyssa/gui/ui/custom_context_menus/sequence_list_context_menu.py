from PyQt5 import QtWidgets


class SequenceListContextMenu:
    """A custom context menu wrapper for the sequence list view of the main dialog."""

    def __init__(self):
        self._context_menu = QtWidgets.QMenu()
        self._rename_sequence_action = QtWidgets.QAction("Rename Selected Protein")
        self._help_action = QtWidgets.QAction("Get Help")
        
        self._context_menu.addAction(self._rename_sequence_action)
        self._context_menu.addAction(self._help_action)

        self._rename_sequence_action.triggered.connect(self._generic_action)
        self._help_action.triggered.connect(self._generic_action)

    # <editor-fold desc="Connect actions from outside">
    def connect_rename_sequence_action(self, the_function_to_connect):
        self._rename_sequence_action.triggered.disconnect()
        self._rename_sequence_action.triggered.connect(the_function_to_connect)

    def connect_help_action(self, the_function_to_connect):
        self._help_action.triggered.disconnect()
        self._help_action.triggered.connect(the_function_to_connect)
    # </editor-fold>
    
    def _generic_action(self):
        raise NotImplementedError("The called action is not connected to a function!")
    
    def get_context_menu(self, the_selected_indexes):
        if len(the_selected_indexes) == 1:
            level = 0
            index = the_selected_indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1
        elif len(the_selected_indexes) > 1:
            # More than one sequence is selected
            self._help_action.setVisible(False)
            self._rename_sequence_action.setVisible(False)
            return self._context_menu
        else:
            # No sequences are selected
            self._help_action.setVisible(False)
            self._rename_sequence_action.setVisible(False)
            return self._context_menu
        # One sequence is selected
        self._help_action.setVisible(False)
        self._rename_sequence_action.setVisible(True)

        return self._context_menu
