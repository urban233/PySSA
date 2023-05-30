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


@dataclass
class HotspotInfo:
    prot_1: str
    prot_1_seq_start: int
    prot_1_seq_end: int
    prot_2: str
    prot_2_seq_start: int
    prot_2_seq_end: int

    def __init__(self, prot_1: str, prot_1_seq_start: int, prot_1_seq_end: int, prot_2: str, prot_2_seq_start: int,
                 prot_2_seq_end: int):
        self.prot_1 = prot_1
        self.prot_1_seq_start = prot_1_seq_start
        self.prot_1_seq_end = prot_1_seq_end
        self.prot_2 = prot_2
        self.prot_2_seq_start = prot_2_seq_start
        self.prot_2_seq_end = prot_2_seq_end

    def get_prot_1(self):
        return self.prot_1

    def get_prot_1_seq_start(self):
        return self.prot_1_seq_start

    def get_prot_1_seq_end(self):
        return self.prot_1_seq_end

    def get_prot_2(self):
        return self.prot_2

    def get_prot_2_seq_start(self):
        return self.prot_2_seq_start

    def get_prot_2_seq_end(self):
        return self.prot_2_seq_end
