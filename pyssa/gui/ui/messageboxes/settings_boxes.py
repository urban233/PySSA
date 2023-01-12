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
"""Module for all message boxes which occur in the settings dialog."""

from PyQt5.QtWidgets import QMessageBox


def installation_local_colabfold_accept() -> bool:
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Question)
    msg.setText("Are you sure that you want to install Local Colabfold?")
    msg.setWindowTitle("Local Colabfold installation")
    btn_installation_local_colabfold_accept_yes = msg.addButton("Yes", QMessageBox.ActionRole)
    btn_installation_local_colabfold_accept_no = msg.addButton("No", QMessageBox.ActionRole)
    msg.exec_()

    if msg.clickedButton() == btn_installation_local_colabfold_accept_yes:
        return True
        # you can better return a True, which gets evaluated in the dialog_settings_global.py through an
    else:
        msg.close()
        return False


def installation_local_colabfold_failed(message):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText(message)
    msg.setWindowTitle("Local Colabfold installation")
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec_()


def installation_local_colabfold_end() -> bool:
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Information)
    msg.setText("Installation is finished!")
    msg.setWindowTitle("Local Colabfold installation")
    msg.setStandardButtons(QMessageBox.Ok)

    if msg.clickedButton() == msg.setStandardButtons(QMessageBox.Ok):
        return True
    else:
        msg.close()
        return False


def installation_local_colabfold_remove() -> bool:
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Question)
    msg.setText("Are you sure that you want to remove Local Colabfold from your system?")
    msg.setWindowTitle("Remove Local Colabfold")
    btn_installation_local_colabfold_accept_yes = msg.addButton("Yes", QMessageBox.ActionRole)
    btn_installation_local_colabfold_accept_no = msg.addButton("No", QMessageBox.ActionRole)
    msg.exec_()

    if msg.clickedButton() == btn_installation_local_colabfold_accept_yes:

        return True
    else:
        msg.close()
        return False
