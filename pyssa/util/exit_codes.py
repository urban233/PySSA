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
"""Module for exit code tuples."""
EXIT_CODE_ZERO = (0, "Process finished without errors.")
EXIT_CODE_ONE_UNKNOWN_ERROR = (1, "Unknown error occurred.")
ERROR_WRITING_FASTA_FILES = (99, "Error while writing FASTA file(s).")
ERROR_FASTA_FILES_NOT_FOUND = (98, "Error not found FASTA file(s)")
ERROR_PREDICTION_FAILED = (97, "Error while running prediction.")
ERROR_COLABFOLD_MODEL_NOT_FOUND = (96, "Error not found colabfold rank 1 model.")

ERROR_DISTANCE_ANALYSIS_FAILED = (87, "Error while running the distance analysis.")
