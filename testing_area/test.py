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
import json
import pathlib


def search_filesystem(protein_name_seq_tuples):
    filenames = []
    filename = []
    for tmp_prediction in protein_name_seq_tuples:
        for tmp_filename in os.listdir(pathlib.Path("C:/Users/martin/.pyssa/scratch/local_predictions/pdb")):
            test = tmp_filename.find(f"{tmp_prediction[0]}_relaxed_rank_001")
            if tmp_filename.find(f"{tmp_prediction[0]}_relaxed_rank_001") != -1:
                filename.append(tmp_filename)
        if len(filename) == 1:
            filenames.append((tmp_prediction, filename[0]))
    return filenames


if __name__ == '__main__':
    prot = [('re2', 'GFGGFTAGT')]
    print(search_filesystem(prot))
