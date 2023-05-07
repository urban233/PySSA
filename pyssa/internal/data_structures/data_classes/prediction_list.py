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
from dataclasses import dataclass
from pyssa.util import constants


@dataclass
class PredictionList:
    protein_name: str
    protein_sequence: list[str]

    def get_tuple_notation(self):
        tuple_notation: tuple[str, list[str]] = (self.protein_name, self.protein_sequence)
        return tuple_notation

    def write_fasta_file(self):
        fasta_file = open(f"{constants.PREDICTION_FASTA_DIR}/{self.protein_name}.fasta", "w")
        fasta_file.write(f">{self.protein_name}\n")
        i = 0
        for sequence in self.protein_sequence:
            if i == len(self.protein_sequence) - 1:
                # should be the last entry
                fasta_file.write(sequence)
            else:
                fasta_file.write(f"{sequence}:")
            i += 1

        fasta_file.close()
