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
import subprocess
import os
import logging
from pymol import Qt
from pyssa.gui.utilities import gui_utils
from pyssa.gui.ui.forms.auto_generated.auto_dialog_settings_global import Ui_Dialog
#from main_alt import global_var_settings_obj
from pyssa.gui.utilities.global_variables import global_var_settings_obj

# setup logger
logging.basicConfig(level=logging.DEBUG)


class DialogSettingsGlobal(Qt.QtWidgets.QDialog):
    """This class opens a settings customization dialog.

    """
    """This variable is for controlling whether the dialog opens or not"""
    ERROR = False

    def __init__(self, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.ui.txt_workspace_dir.setEnabled(False)
        self.ui.txt_zip_storage_dir.setEnabled(False)

        # try:
        #     self.xmlFile = self.xmlObj.load_xml_in_memory()
        # except IsADirectoryError:
        #     logging.error("There is only a directory and not a file.")
        #     gui_utils.error_dialog("There is only a directory.", "More information in the log.")
        #     # is used to stop opening this dialog
        #     self.ERROR = True
        # except FileNotFoundError:
        #     logging.error("The settings.xml file was not found and could therefore not be opened.")
        #     gui_utils.error_dialog("The settings.xml file was not found and could therefore not be opened.", "")
        #     # is used to stop opening this dialog
        #     self.ERROR = True

        logging.info("Settings dialog was opened.")
        # loading information from the settings.xml
        self.ui.txt_workspace_dir.setText(global_var_settings_obj.get_workspace_path())
        self.ui.txt_zip_storage_dir.setText(global_var_settings_obj.get_prediction_path())
        self.ui.spb_cycles.setValue(int(global_var_settings_obj.get_cycles()))
        self.ui.dspb_cutoff.setValue(float(global_var_settings_obj.get_cutoff()))
        logging.info("Loading values from settings.xml was successful.")
        # customize spin boxes
        self.ui.spb_cycles.setMinimum(0)
        # self.ui.spb_cycles.setMaximum(20) # is a maximum needed?
        self.ui.spb_cycles.setSingleStep(1)
        self.ui.dspb_cutoff.setMinimum(0.00)
        self.ui.dspb_cutoff.setMaximum(20.00)
        self.ui.dspb_cutoff.setSingleStep(0.1)

        # connect elements with function
        self.ui.btn_workspace_dir.clicked.connect(self.chooseWorkspaceDir)
        self.ui.btn_zip_storage_dir.clicked.connect(self.chooseZipStorageDir)
        self.ui.btn_cancel.clicked.connect(self.cancelDialog)
        self.ui.btn_ok.clicked.connect(self.okDialog)
        self.ui.btn_install_local_prediction.clicked.connect(self.install_local_colabfold)

        self.setWindowTitle("Global Settings")

    # @SLOT()
    def chooseWorkspaceDir(self):
        gui_utils.choose_directory(self, self.ui.txt_workspace_dir)

    def chooseZipStorageDir(self):
        gui_utils.choose_directory(self, self.ui.txt_zip_storage_dir)

    def cancelDialog(self):
        self.close()

    def okDialog(self):
        global_var_settings_obj.set_workspace_path(self.ui.txt_workspace_dir.text())
        global_var_settings_obj.set_prediction_path(self.ui.txt_zip_storage_dir.text())
        global_var_settings_obj.set_cycles(str(self.ui.spb_cycles.value()))
        global_var_settings_obj.set_cutoff(str(self.ui.dspb_cutoff.value()))
        global_var_settings_obj.save_settings_to_xml()
        self.close()

    def install_local_colabfold(self):
        print(subprocess.run(f"{os.environ['HOME']}/github_repos/tmpPySSA/pyssa/scripts/installation_colabfold.sh"))
