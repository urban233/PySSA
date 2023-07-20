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
import os

from pymol import cmd


if __name__ == '__main__':
    os.startfile("C:\\Users\\martin\\scratch\\test.html")

    # cmd.fetch("3bmp")
    # cmd.fetch("6omn")
    # cmd.align("/3bmp//A//CA", "/6omn//E//CA", cycles=0, object="aln")
    # raw_aln = cmd.get_raw_alignment(name="aln")
    # # print residue pairs (atom index)
    # for idx1, idx2 in raw_aln:
    #     print('%s`%d -> %s`%d' % tuple(idx1 + idx2))
    #
    # #idx2resi = {}
    # #cmd.iterate('aln', 'idx2resi[model, index] = resi', space={'idx2resi': idx2resi})
    # idx2resi = []
    # cmd.iterate('aln', 'idx2resi.append((model, chain, resi, resn))', space={'idx2resi': idx2resi})
    # print(idx2resi)
    # prot_1_indices = []
    # prot_1_name = "3bmp"
    # prot_2_indices = []
    # prot_2_name = "6omn"
    # for tmp_prot_atom in idx2resi:
    #     if tmp_prot_atom[0] == prot_1_name:
    #         prot_1_indices.append(tmp_prot_atom[1])
    #     if tmp_prot_atom[0] == prot_2_name:
    #         prot_2_indices.append(tmp_prot_atom[1])
    # print(prot_1_indices)
    # print(prot_2_indices)
    # # calculate the distance between the alpha-C atoms
    # for resi_no in range(len(prot_1_indices)):
    #     atom1 = f"/{prot_1_name}//A/{prot_1_indices[resi_no]}/CA"
    #     atom2 = f"/{prot_2_name}//E/{prot_2_indices[resi_no]}/CA"
    #     distance = round(cmd.get_distance(atom1, atom2), 2)
    #     print(distance)

