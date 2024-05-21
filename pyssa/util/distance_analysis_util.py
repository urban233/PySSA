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
"""Collection of helper functions for the distance analysis."""
import logging
import os
from typing import TYPE_CHECKING
from pyssa.logging_pyssa import log_handlers
from pyssa.util import protein_pair_util, constant_messages, constants
from pyssa.util import exception
from pyssa.internal.data_structures import results

if TYPE_CHECKING:
  from pyssa.internal.data_structures import protein_pair

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

#
# def create_scene_of_protein_pair(
#     a_protein_pair: "protein_pair.ProteinPair",
#     selection: str = "",
#     ray_shadows: bool = False,
#     opaque_background: int = 0,
#     take_images: bool = False,
#     filename: str = "",
# ) -> None:
#     """This function takes an image of the whole protein pair.
#
#     Args:
#         a_protein_pair: defines the type of molecular representation like cartoon or ribbon.
#         selection: the atoms which MUST NOT displayed in the image.
#         ray_shadows (optional): false if no shadows, true if shadows should be displayed.
#         opaque_background (optional): 0 for a transparent background and 1 for a white background.
#         take_images (optional): Is a boolean, True relates to creating images, False no images will be made.
#         filename (optional): name of the png image file without the extension.
#
#     Raises:
#         ValueError: If opaque_background is not 0 or 1.
#         UnableToColorProteinPairError: If the protein pair is not able to colored.
#     """
#     tmp_representation = "cartoon"
#     cmd.show(tmp_representation)
#     if selection != "":
#         cmd.hide(tmp_representation, selection)
#     cmd.hide("cgo", "all")
#     cmd.orient()
#     cmd.center()
#     # set image parameters
#     graphic_operations.setup_default_image_graphic_settings(ray_shadows, opaque_background)
#     graphic_operations.setup_default_session_graphic_settings()
#
#     cmd.scene(
#         key=f"{a_protein_pair.protein_1.get_molecule_object()}-{a_protein_pair.protein_2.get_molecule_object()}",
#         action="store",
#     )
#
#     if take_images is True:
#         if not os.path.exists(f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}"):
#             os.mkdir(f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}")
#         # save image as 300 dpi png image
#         cmd.ray(2400, 2400, renderer=0)
#         cmd.png(f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/{filename}.png", dpi=300)
#
#
# def create_scenes_of_interesting_regions(
#     the_distance_data: dict[str, np.ndarray],
#     protein_1_molecule_object: str,
#     protein_2_molecule_object: str,
#     cutoff: float,
#     ray_shadows: bool = False,
#     take_images: bool = False,
#     filename: str = "",
# ) -> None:
#     """Creates scenes of interesting regions of the alignment and optionally ray-traced images.
#
#     Args:
#         the_distance_data: Dictionary of distance data.
#         protein_1_molecule_object: Name of the protein 1.
#         protein_2_molecule_object: Name of the protein 2.
#         cutoff (float): Defines a border of which the specific regions begin, if the distance is greater than
#                         the cutoff, the amino acid is categorized as "interesting".
#         filename (str): Is the name of the png image file.
#         ray_shadows (bool): Is false if no shadows, true if shadows should be displayed.
#         opaque_background (int): 0 for a transparent background and 1 for a white background.
#         take_images (bool): Flag if images should be made or not.
#
#     Raises:
#         IllegalArgumentError: If var take_images is true but the filename is not set.
#         UnableToColorProteinPairError: If the protein pair is not able to colored.
#     """
#     # <editor-fold desc="Checks">
#     if take_images and filename == "":
#         tmp_msg: str = "Image creation is activated but no filename provided!"
#         logger.error(tmp_msg)
#         raise exception.IllegalArgumentError(tmp_msg)
#
#     # </editor-fold>
#
#     # set default parameters
#     graphic_operations.setup_default_session_graphic_settings()
#     graphic_operations.setup_default_graphic_settings_for_interesting_regions()
#     j: int = 0
#     i: int = 0
#     for distance_value in the_distance_data.get("distance"):
#         if float(distance_value) > float(cutoff):
#             # here algorithm for image
#             ref_pos_array: np.ndarray = the_distance_data.get("ref_pos")
#             ref_pos: int = ref_pos_array[i]
#             ref_chain_array: np.ndarray = the_distance_data.get("ref_chain")
#             ref_chain: str = ref_chain_array[i]
#             model_pos_array: np.ndarray = the_distance_data.get("model_pos")
#             model_pos: int = model_pos_array[i]
#             model_chain_array: np.ndarray = the_distance_data.get("model_chain")
#             model_chain: str = model_chain_array[i]
#
#             measurement_obj: str = f"measure{j}"
#             # create two atoms for the get_distance command
#             atom1: str = f"/{protein_1_molecule_object}//{ref_chain}/{ref_pos}/CA"
#             atom2: str = f"/{protein_2_molecule_object}//{model_chain}/{model_pos}/CA"
#             # zoom to reference amino acid
#             cmd.zoom(f"/{protein_1_molecule_object}//{ref_chain}/{ref_pos}", 10)
#             # create distance object with get_distance command
#             cmd.distance(measurement_obj, atom1, atom2)
#             cmd.label(atom1, "'%s-%s' % (resn, resi)")
#             cmd.label(atom2, "'%s-%s' % (resn, resi)")
#             cmd.set("label_position", (0, 0, 10))
#             # set image parameters
#             graphic_operations.setup_default_image_graphic_settings(ray_shadows)
#             cmd.scene(key=f"{ref_pos}-{model_pos}", action="store")
#             if take_images is True:
#                 if not os.path.exists(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
#                     os.mkdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR)
#                 # save image as 300 dpi png image
#                 cmd.ray(2400, 2400, renderer=0)
#                 cmd.png(
#                     f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{filename}_{ref_pos}.png",
#                     dpi=300,
#                 )
#             # hide created labels
#             cmd.hide("labels", atom1)
#             cmd.hide("labels", atom2)
#             cmd.hide("labels", measurement_obj)
#             cmd.hide("dashes", measurement_obj)
#             i += 1
#             j += 1
#         else:
#             i += 1
