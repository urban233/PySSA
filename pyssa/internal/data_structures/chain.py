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
"""Module contains the chain class."""
import logging
from PyQt5 import QtCore

from pyssa.controller import pymol_session_manager
from pyssa.internal.portal import graphic_operations
from pyssa.util import enums
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures import sequence

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Chain:
    """This class contains information about a single chain."""

    # <editor-fold desc="Class attributes">
    _id: int
    """
    a letter of the chain
    """
    chain_letter: str
    """
    a sequence of the chain
    """
    chain_sequence: sequence.Sequence
    """
    a type of the chain, whether it is a protein, or nuclein acid chain or something different
    """
    chain_type: str
    """
    a dict of parameters that can be changed in pymol
    """
    pymol_parameters: dict
    """
    the protein id from the project database
    """
    db_protein_id: int

    # </editor-fold>

    def __init__(self, chain_letter: str, chain_sequence: sequence.Sequence, chain_type: str) -> None:
        """Constructor.

        Args:
            chain_letter:
                letter of the chain
            chain_sequence:
                sequence of the chain
            chain_type:
                type of the chain, whether it is a protein, or nuclein acid chain or something different
        """
        self.chain_letter = chain_letter
        self.chain_sequence = chain_sequence
        self.chain_type = chain_type
        self.pymol_parameters = {
            enums.PymolParameterEnum.COLOR.value: "green",
            enums.PymolParameterEnum.REPRESENTATION.value: "cartoon",
        }

    def get_id(self):
        return self._id

    def set_id(self, value):
        self._id = value

    def get_color(self, a_selection_string, the_pymol_session_manager: "pymol_session_manager.PymolSessionManager") -> tuple[str, bool]:
        self.pymol_parameters[enums.PymolParameterEnum.COLOR.value], _, tmp_is_colored_by_elements = the_pymol_session_manager.get_chain_color(a_selection_string, self.chain_letter)
        return self.pymol_parameters[enums.PymolParameterEnum.COLOR.value], tmp_is_colored_by_elements

    def get_representation_state(self, a_selection_string) -> dict:
        return graphic_operations.get_chain_repr_state(a_selection_string, self.chain_letter)

    def serialize(self, an_xml_writer: QtCore.QXmlStreamWriter):
        an_xml_writer.writeStartElement('chain')
        an_xml_writer.writeAttribute('chain_letter', str(self.chain_letter))
        an_xml_writer.writeAttribute('chain_sequence', str(self.chain_sequence.sequence))
        an_xml_writer.writeAttribute('chain_type', str(self.chain_type))
        an_xml_writer.writeStartElement('pymol_parameters')
        an_xml_writer.writeAttribute('chain_color', str(self.pymol_parameters['chain_color']))
        an_xml_writer.writeAttribute('chain_representation', str(self.pymol_parameters['chain_representation']))
        an_xml_writer.writeEndElement()
        an_xml_writer.writeEndElement()
