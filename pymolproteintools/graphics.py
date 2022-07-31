# ----------------------------------------------------------------------- #
# Python Package: pymolproteintools
# Version 1.0 for Python 3.9
# -----------------------------------------------------------------------#
# Authors: Martin Urban & Hannah Kullik
# Westfaelische Hochschule
# Recklinghausen, Germany
#
# Copyright 2022 Martin Urban & Hannah Kullik
#
# Citation:
# Martin Urban & Hannah Kullik, pymolproteintools (PPT),
# Version 1.0, Recklinghausen, Germany, 2022.
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License (LGPL) as
# published by the Free Software Foundation, version 3 of the License.
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License (LGPL) for more details.
# You should have received a copy of the GNU Lesser General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

# -----------------------------------------------------------------------#
import os
from pymol import cmd
from typing import Dict, Tuple
from pymolproteintools import utility
from pymolproteintools import core
import matplotlib.pyplot as plt
import numpy as np


class graphics:

    def __init__(self, protein_pair: core.proteinpair,
                 results_hashtable: Dict[str, np.ndarray],
                 figure_size: Tuple[float, float]):
        """Constructor

        Args:
             protein_pair (core.proteinpair):
                protein pair object
             results_hashtable (dict[str, np.ndarray]):
                hash table which contains the results of the distance
                calculations
             figure_size (tuple[float, float]):
                figure size of the distance plot and distance histogram
        """
        self.protein_pair: core.proteinpair = protein_pair
        self.results_hashtable: Dict[str, np.ndarray] = results_hashtable
        self.figure_size: Tuple[float, float] = figure_size

    def take_image_of_protein_pair(self,
                                   alignment_filename: str,
                                   representation: str,
                                   filename: str,
                                   selection: str = "",
                                   ray_shadows: bool = False,
                                   opaque_background: int = 0) -> None:
        """This function takes an image of the whole protein/protein pair.

        Note:
            The png file will be saved under the relative path
            (if export_data_dir = "data/results"):
            ``data/results/images``

        Args:
            alignment_filename (str):
                name of the alignment file from the structure alignment
            representation (str):
                defines the type of molecular representation
                like cartoon or ribbon
            filename (str):
                name of the png image file
            selection (str, optional):
                the atoms which MUST NOT displayed in the image
            ray_shadows (bool, optional):
                false if no shadows, true if shadows should be displayed
            opaque_background (int, optional):
                0 for a transparent background and 1 for a white background

        Raises:
            ValueError: If opaque_background is not 0 or 1.
        """
        # argument test
        # if opaque_background != 0 or opaque_background != 1:
        #     raise Exception(
        #         "ValueError: The value for opaque_background MUST be 0 or 1!")

        # determine the option for ray_shadows
        if ray_shadows == False:
            opt_ray_shadows: str = "off"
        else:
            opt_ray_shadows: str = "on"

        REPRESENTATION: str = "cartoon"
        cmd.show(REPRESENTATION)

        if selection != "":
            cmd.hide(representation, selection)

        aln_obj_representation: str = "cgo"
        cmd.hide(aln_obj_representation, alignment_filename)
        cmd.orient()
        cmd.center()
        # set image parameters
        cmd.bg_color("white")
        cmd.set("ray_trace_mode", 3)
        cmd.set("antialias", 2)
        cmd.set("ray_shadows", opt_ray_shadows)
        cmd.set('ray_opaque_background', opaque_background)
        cmd.ray(2400, 2400, renderer=0)

        cmd.scene(key=f"{self.protein_pair.ref_obj.molecule_object}-{self.protein_pair.model_obj.molecule_object}", action="store")

        # check if path exists where the data will be exported,
        # if not the directory will be created
        if not os.path.exists(f"{self.protein_pair.results_dir}/images"):
            os.mkdir(f"{self.protein_pair.results_dir}/images")

        # save image as 300 dpi png image
        cmd.png(f'{self.protein_pair.results_dir}/images/{filename}.png',
                dpi=300)

    def take_image_of_interesting_regions(self,
                                          cutoff: float,
                                          filename: str,
                                          ray_shadows: bool = False,
                                          opaque_background: int = 0):
        """This function takes images of interesting regions of the alignment

        Args:
            cutoff (float):
                defines a border of which the specific regions begin,
                if the distance is greater than the cutoff, the amino acid is
                categorized as "interesting"
            filename (str):
                name of the png image file
            ray_shadows (bool, optional):
                false if no shadows, true if shadows should be displayed
            opaque_background (int, optional):
                0 for a transparent background and 1 for a white background

        """

        # set default parameters
        cmd.set("label_size", 14)
        cmd.set("label_font_id", 13)
        cmd.set("label_color", "hotpink")
        cmd.set("depth_cue", 0)
        # cmd.set("fog_start", 0.6)
        cmd.bg_color("white")

        j: int = 0
        measurement_obj: str = f"measure{j}"
        REPRESENTATION: str = "ribbon"

        cmd.hide("cartoon", "all")
        cmd.show(REPRESENTATION, "all")

        i: int = 0
        for distance_value in self.results_hashtable.get("distance"):
            if distance_value > cutoff:
                # here algorithm for image
                ref_pos_array: np.ndarray = self.results_hashtable.get("ref_pos")
                ref_pos: int = ref_pos_array[i]
                ref_chain_array: np.ndarray = self.results_hashtable.get("ref_chain")
                ref_chain: str = ref_chain_array[i]
                model_pos_array: np.ndarray = self.results_hashtable.get("model_pos")
                model_pos: int = model_pos_array[i]
                model_chain_array: np.ndarray = self.results_hashtable.get("model_chain")
                model_chain: str = model_chain_array[i]

                measurement_obj: str = f"measure{j}"
                # create two atoms for the get_distance command
                atom1: str = f"/{self.protein_pair.ref_obj.molecule_object}//" \
                             f"{ref_chain}/{ref_pos}/CA"
                atom2: str = f"/{self.protein_pair.model_obj.molecule_object}//" \
                             f"{model_chain}/{model_pos}/CA"
                # zoom to reference amino acid
                cmd.zoom(f"/{self.protein_pair.ref_obj.molecule_object}//"
                         f"{ref_chain}/{ref_pos}", 10)
                # create distance object with get_distance command
                cmd.distance(measurement_obj, atom1, atom2)
                cmd.label(atom1, "'%s-%s' % (resn, resi)")
                cmd.label(atom2, "'%s-%s' % (resn, resi)")
                cmd.set("label_position", (0,0,10))
                # determine the option for ray_shadows
                if ray_shadows == False:
                    opt_ray_shadows: str = "off"
                else:
                    opt_ray_shadows: str = "on"
                # set image parameters
                cmd.bg_color("white")
                cmd.set("ray_trace_mode", 3)
                cmd.set("antialias", 2)
                cmd.set("ray_shadows", opt_ray_shadows)
                cmd.set('ray_opaque_background', opaque_background)
                cmd.ray(2400, 2400, renderer=0)

                cmd.scene(key=f"{ref_pos}-{model_pos}", action="store")

                # check if path exists where the data will be exported,
                # if not the directory will be created
                if not os.path.exists(
                        f"{self.protein_pair.results_dir}/images"):
                    os.mkdir(f"{self.protein_pair.results_dir}/images")

                # save image as 300 dpi png image
                cmd.png(
                    f'{self.protein_pair.results_dir}/images/{filename}_{ref_pos}.png',
                    dpi=300)
                # hide created labels
                cmd.hide("labels", atom1)
                cmd.hide("labels", atom2)
                cmd.hide("labels",measurement_obj)
                cmd.hide("dashes",measurement_obj)
                i += 1
                j += 1
            else:
                i += 1

    def create_distance_plot(self, distance_label: str, cutoff: float,
                             limit_y_axis: Tuple[float, float] = (0, 0)):
        """This function creates a distance plot.

        Args:
            distance_label (str):
                label which defines the atoms of which the distance was
                calculated
            cutoff (float):
                value which was used for the structure alignment
            limit_y_axis (tuple[float, float], optional):
                defines the min and max value for the y-axis
        """
        x: np.ndarray = self.results_hashtable.get("index")
        y: np.ndarray = self.results_hashtable.get("distance")

        # max distance value
        max_distance = np.amax(y)

        # create an array for cutoff
        y_cutoff: np.ndarray = np.ndarray((x.size, 1), float)
        i: int = 0
        for j in y_cutoff:
            y_cutoff.flat[i] = cutoff
            i += 1

        # create an empty figure with no Axes
        fig = plt.figure()
        # create a figure with a single Axes
        fig, ax = plt.subplots(figsize=self.figure_size)
        # creates a basic scatter plot
        ax.scatter(x, y,
                   label=f"Distance {distance_label} pair")
        # creates a basic line plot
        ax.plot(x, y, linewidth=2.0)
        # creates a basic line plot for the cutoff
        ax.plot(x, y_cutoff, color='red',
                label="Cutoff = " + str(cutoff) + " $\mathring{A}$")

        # sets the label for the x-axis
        ax.set_xlabel('Residue no.')
        # sets the label for the y-axis
        ax.set_ylabel("Distance [$\mathring{A}$ngstrom]")
        # sets ticks for the x-axis
        ax.set_xticks(np.arange(0, x.size, 20))

        if limit_y_axis != (0,0):
            ax.set_ylim(limit_y_axis)
        else:
            # sets the y-axis limit
            ax.set_ylim(0, max_distance + 1)
        # creates legend for all "label" parameters
        ax.legend()
        # sets grid
        ax.grid(True)
        return fig

    def create_distance_histogram(self) -> None:
        """This function creates a histogram of the distances calculated with
        the function calculate_distance_between_ca_atoms()

        """
        y: np.ndarray = self.results_hashtable.get("distance")

        # max distance value
        max_distance = np.amax(y)

        # calculate figure size for y direction
        y_size: float = len(np.arange(0, max_distance, 0.25)) / 1.2
        FIGURE_SIZE: Tuple[float, float] = (11.0, y_size)
        # create an empty figure with no Axes
        fig = plt.figure()
        # create a figure with a single Axes
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        # creates a basic histogram
        counts, bins, patches = ax.hist(y, bins=np.arange(0, max_distance, 0.25), orientation="horizontal")
        # sets the label for the x-axis
        ax.set_xlabel("Frequency of $\\alpha$-C atoms distance")
        # sets the label for the y-axis
        ax.set_ylabel("Distance [$\mathring{A}$ngstrom]")
        # set coordinates for y-axis label
        ax.yaxis.set_label_coords(-0.12, 0.5)
        # hide y-ticks through empty list
        ax.set_yticks([])
        # create label bin position
        bins_centers = 0.15 * np.diff(bins) + bins[:-1]

        i: int = 0
        for count, y in zip(counts, bins_centers):
            # define bin label
            bin_name: str = f"{round(bins[i], 2)} - {round(bins[i + 1], 2)}"
            # set bin label through annotation
            ax.annotate(bin_name, xy=(0, y), xytext=(-70, y), textcoords="offset points")
            i += 1
        # sets grid
        ax.grid(True, axis="both")

    def color_by_distance(self, alignment_filename: str):

        cutoff_1 = 0.5
        cutoff_2 = 1.0
        cutoff_3 = 2
        cutoff_4 = 4
        cutoff_5 = 6

        color_1 = "br0"
        color_2 = "br2"
        color_3 = "br4"
        color_4 = "br6"
        color_5 = "br8"
        color_6 = "red"

        cmd.color("hydrogen", self.protein_pair.model_obj.molecule_object)

        i: int = 0
        for distance_value in self.results_hashtable.get("distance"):
            if distance_value <= cutoff_1:
                atom_info = utility.get_chain_and_position(self.results_hashtable, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{self.protein_pair.ref_obj.molecule_object}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{self.protein_pair.model_obj.molecule_object}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_1, atom1)
                cmd.color(color_1, atom2)
                # cmd.do(f"color {color_1}, {atom1}")
                # cmd.do(f"color {color_1}, {atom2}")
                i += 1

            elif distance_value <= cutoff_2:
                atom_info = utility.get_chain_and_position(
                    self.results_hashtable, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{self.protein_pair.ref_obj.molecule_object}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{self.protein_pair.model_obj.molecule_object}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_2, atom1)
                cmd.color(color_2, atom2)
                i += 1

            elif distance_value <= cutoff_3:
                atom_info = utility.get_chain_and_position(
                    self.results_hashtable, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{self.protein_pair.ref_obj.molecule_object}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}/CA"
                atom2: str = f"/{self.protein_pair.model_obj.molecule_object}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}/CA"
                # coloring
                cmd.color(color_3, atom1)
                cmd.color(color_3, atom2)
                i += 1

            elif distance_value <= cutoff_4:
                atom_info = utility.get_chain_and_position(
                    self.results_hashtable, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{self.protein_pair.ref_obj.molecule_object}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{self.protein_pair.model_obj.molecule_object}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_4, atom1)
                cmd.color(color_4, atom2)
                i += 1

            elif distance_value <= cutoff_5:
                atom_info = utility.get_chain_and_position(
                    self.results_hashtable, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{self.protein_pair.ref_obj.molecule_object}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{self.protein_pair.model_obj.molecule_object}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_5, atom1)
                cmd.color(color_5, atom2)
                i += 1

            elif distance_value > cutoff_5:
                atom_info = utility.get_chain_and_position(
                    self.results_hashtable, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{self.protein_pair.ref_obj.molecule_object}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{self.protein_pair.model_obj.molecule_object}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_6, f"({atom1})")
                cmd.color(color_6, f"({atom2})")
                # cmd.do(f"color {color_6}, {atom1}")
                # cmd.do(f"color {color_6}, {atom2}")
                i += 1

        # hide unnecessary representations
        cmd.hide("cartoon", self.protein_pair.ref_obj.molecule_object)
        cmd.hide("cartoon", f"{self.protein_pair.ref_obj.molecule_object}_CA")
        cmd.hide("cartoon", f"{self.protein_pair.model_obj.molecule_object}_CA")
        cmd.hide("cgo", alignment_filename)

    def create_gif(self):
        """This function creates a gif.

        """
        # create first scene
        cmd.scene("001", "store")
        # turn camera 180 degrees
        cmd.turn("y", 180)
        # create a second scene
        cmd.scene("002", "store")
        # define length of gif
        cmd.mset("1x240")
        # store first scene at frame 1
        cmd.mview("store", 1, scene="001")
        # store second scene at frame 120
        cmd.mview("store", 120, scene="002")
        # store first scene at frame 240
        cmd.mview("store", 240, scene="001")

        # check if path exists where the data will be exported,
        # if not the directory will be created
        if not os.path.exists(
                f"{self.protein_pair.results_dir}/gifs"):
            os.mkdir(f"{self.protein_pair.results_dir}/gifs")
        # save movie as gif
        cmd.movie.produce("rotation.gif", quality=1.0)
