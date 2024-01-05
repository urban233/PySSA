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
"""Module for the main view of the pyssa plugin."""
import pathlib

from PyQt5 import QtWidgets
from PyQt5 import QtGui

from pyssa.gui.ui.forms.auto_generated import auto_main_window
from pyssa.util import constants


class MainView(QtWidgets.QMainWindow):
    """Class representing the main view of PySSA."""

    """
    The status bar of the main view.
    """

    def __init__(self) -> None:
        """Constructor."""
        super().__init__()
        # build ui object
        self.ui = auto_main_window.Ui_MainWindow()
        self.ui.setupUi(self)
        self.status_bar = QtWidgets.QStatusBar()
        self.initialize_ui()

    def initialize_ui(self) -> None:
        """Initialize the UI elements."""
        self.ui.lbl_page_title.setText("Home")
        self.setMinimumWidth(580)
        self.setMinimumHeight(200)
        # sets additional parameters
        pixmapi = QtWidgets.QStyle.SP_MessageBoxQuestion
        icon = self.style().standardIcon(pixmapi)
        self.ui.btn_info.setIcon(icon)
        self.ui.btn_info.setText("")
        self.ui.btn_info.setFixedWidth(50)

        self.status_bar = QtWidgets.QStatusBar()

        # create tooltips
        self._create_all_tooltips()

        self.ui.lbl_logo.setPixmap(
            QtGui.QPixmap(str(pathlib.Path(f"{constants.PLUGIN_ROOT_PATH}/assets/images/pyssa_logo.png"))),
        )
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("PySSA")
        constants.PYSSA_LOGGER.info(f"PySSA started with version {constants.VERSION_NUMBER}.")
        constants.PYSSA_LOGGER.info("Successful initialization of basic UI.")

    def _create_all_tooltips(self) -> None:
        """Creates all tooltips for the gui elements."""
        self.status_bar.setToolTip("Status information: Current process")
        # new project page
        self.ui.btn_new_choose_reference.setToolTip("Click to add a .pdb file")
        # sidebar
        self.ui.lbl_current_project_name.setToolTip("Name of the current project")
        # edit page
        self.ui.btn_edit_project_save.setToolTip("Save as a .pdb file")
        # view page
        # fixme: is this important? self.ui.list_view_project_proteins.setToolTip("Proteins of the current project")
        self.ui.txtedit_view_sequence.setToolTip("Protein sequence of the selected protein")
        # use page
        self.ui.txt_use_search.setToolTip("Enter a protein name to search in your current workspace")
        # prediction Monomer
        self.ui.table_pred_mono_prot_to_predict.setToolTip("Protein monomers which get predicted")
        self.ui.btn_pred_mono_seq_to_predict.setToolTip("Set up a protein which can be used for a prediction")
        self.ui.table_pred_analysis_mono_prot_to_predict.setToolTip("Protein monomers which get predicted")
        self.ui.list_pred_analysis_mono_overview.setToolTip("Protein pairs which get analyzed")
        self.ui.btn_pred_analysis_mono_seq_to_predict.setToolTip("Set up a protein which can be used for a prediction")
        # prediction Multimer
        self.ui.table_pred_multi_prot_to_predict.setToolTip("Protein multimers which get predicted")
        self.ui.btn_pred_multi_prot_to_predict_add.setToolTip("Set up a protein which can be used for a prediction")
        self.ui.table_pred_analysis_multi_prot_to_predict.setToolTip("Protein multimers which get predicted")
        self.ui.list_pred_analysis_multi_overview.setToolTip("Protein pairs which get analyzed")
        self.ui.btn_pred_analysis_multi_prot_to_predict_add.setToolTip(
            "Set up a protein which can be used for a prediction",
        )
        # image page
        self.ui.btn_save_scene.setToolTip("Create new PyMOL scene")
        self.ui.btn_update_scene.setToolTip("Overwrite current scene")
        self.ui.btn_preview_image.setToolTip("Preview current viewpoint")
        self.ui.btn_save_image.setToolTip("Save current viewpoint as png file")
        self.ui.cb_ray_tracing.setToolTip("Enable ray-tracing")
        self.ui.cb_transparent_bg.setToolTip("Enable transparent background")
        self.ui.box_representation.setToolTip("Choose a representation")
        self.ui.box_bg_color.setToolTip("Choose a background color")
        self.ui.box_renderer.setToolTip("Choose a ray-tracing renderer")
        self.ui.box_ray_trace_mode.setToolTip("Choose a ray-trace mode")
