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
"""Module for the histogram properties view."""
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.gui.ui.custom_widgets import custom_line_edit
from pyssa.gui.ui.forms.auto_generated import auto_histogram_properties_view
from pyssa.gui.ui.styles import styles
from pyssa.util import constants, enums, gui_utils


class HistogramPropertiesView(QtWidgets.QDialog):
    """Class representing a Create dialog."""
    new_properties = QtCore.pyqtSignal(tuple)

    def __init__(self, the_current_properties: dict) -> None:
        """Constructor.
        """
        QtWidgets.QDialog.__init__(self)
        # build ui object
        self.ui = auto_histogram_properties_view.Ui_Dialog()
        self.ui.setupUi(self)
        self._initialize_ui()
        self._fill_combo_box()
        self._set_current_properties_in_ui_elements(the_current_properties)
        self._connect_all_signals()
        self.ui.btn_cancel.clicked.connect(self.close)
        self.setModal(True)

    def _initialize_ui(self) -> None:
        """Initialize the UI elements."""
        self.ui.le_units_x_axis.setStyleSheet("""QLineEdit {min-width: 50px; max-width: 50px;}""")
        #self.ui.le_units_x_axis = custom_line_edit.CustomLineEditForEnteringNumbers()
        self.ui.cb_distance_interval.setStyleSheet("""QComboBox {min-width: 55px; max-width: 55px;}""")

        self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
        self.ui.btn_help.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help.setText("")
        styles.color_bottom_frame_button(self.ui.btn_save)
        styles.set_stylesheet(self)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Histogram Properties")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

    def _set_current_properties_in_ui_elements(self, the_current_properties: dict):
        self.ui.le_units_x_axis.setText(str(the_current_properties[enums.HistogramPropertiesEnum.X_AXIS_UNITS]))
        tmp_index = self.ui.cb_distance_interval.findText(str(the_current_properties[enums.HistogramPropertiesEnum.DISTANCE_INTERVAL]))
        self.ui.cb_distance_interval.setCurrentIndex(tmp_index)

    def _fill_combo_box(self):
        tmp_items = [
            "0.25",
            "0.5",
            "0.75",
            "1.0",
            "2.0",
            "5.0",
        ]
        gui_utils.fill_combo_box(self.ui.cb_distance_interval, tmp_items)

    def _connect_all_signals(self):
        self.ui.le_units_x_axis.textChanged.connect(self.validate_x_axis_units_input)
        self.ui.btn_save.clicked.connect(self.apply_properties)

    def validate_x_axis_units_input(self, text):
        print(text)
        allowed_chars = set("0123456789")
        new_text = ''.join(char for char in text if char in allowed_chars)
        self.ui.le_units_x_axis.setText(new_text)
        if new_text == '':
            self.ui.btn_save.setEnabled(False)
            self.ui.le_units_x_axis.setToolTip("A number is needed!")
            self.ui.le_units_x_axis.setStyleSheet(
                """QLineEdit {min-width: 50px; max-width: 50px; color: #ba1a1a; border-color: #ba1a1a;}"""
            )
        elif new_text == '0' or len(new_text) > 3:
            self.ui.btn_save.setEnabled(False)
            self.ui.le_units_x_axis.setToolTip("It is not allowed to enter a number greater 999!")
            self.ui.le_units_x_axis.setStyleSheet(
                """QLineEdit {min-width: 50px; max-width: 50px; color: #ba1a1a; border-color: #ba1a1a;}"""
            )
        elif new_text[0] == '0':
            self.ui.btn_save.setEnabled(False)
            self.ui.le_units_x_axis.setToolTip("It is not allowed to start with a 0!")
            self.ui.le_units_x_axis.setStyleSheet(
                """QLineEdit {min-width: 50px; max-width: 50px; color: #ba1a1a; border-color: #ba1a1a;}"""
            )
        elif len(new_text) == 3:
            self.ui.btn_save.setEnabled(True)
            self.ui.le_units_x_axis.setToolTip("You have entered a unit greater 100. Are you sure that this is correct?")
            self.ui.le_units_x_axis.setStyleSheet(
                """QLineEdit {min-width: 50px; max-width: 50px; color: #ff9500; border-color: #ff9500;}"""
            )
        else:
            self.ui.btn_save.setEnabled(True)
            self.ui.le_units_x_axis.setToolTip("")
            self.ui.le_units_x_axis.setStyleSheet(
                """QLineEdit {min-width: 50px; max-width: 50px; color: #000000; border-color: #DCDBE3;}"""
            )

    def apply_properties(self):
        self.close()
        self.new_properties.emit((self.ui.le_units_x_axis.text(), self.ui.cb_distance_interval.currentText()))
