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



if __name__ == '__main__':
    dict_te = {
        "list": [1,2,3,455,5]
    }
    tmp_object_file = open(("test.json"), "w", encoding="utf-8")
    json.dump(dict_te, tmp_object_file, indent=4)
    tmp_object_file.close()

    tmp_object_file = open("test.json", "r", encoding="utf-8")
    object_dict = json.load(tmp_object_file)
    tmp_object_file.close()

    print(object_dict.get("list")[4])
