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

import os
import sys
import time
from pymolproteintools_legacy import core
from pymolproteintools_legacy import graphics
from pymol import cmd
import matplotlib.pyplot as plt

# if upload_pdb == True:
#     # define reference Protein if a .pdb was uploaded
#     REFERENCE_OBJ_NAME: str = ref_name
# else:
#     # define reference Protein if a pdb id was given
#     REFERENCE_OBJ_NAME: str = pdb_id
#     cmd.set("fetch_type_default", "pdb")
#     cmd.fetch(REFERENCE_OBJ_NAME)
#     # move the fetched pdb file to the data directory
#     src_path: str = f"/content/{pdb_name}"
#     dest_path: str = f"/content/data/{pdb_name}"
#     shutil.move(src_path, dest_path)

if __name__ == "__main__":
    # import of reference Protein name from command line option
    REFERENCE_OBJ_NAME = sys.argv[1]

    t0 = time.time()
    # REFERENCE_OBJ_NAME = "pymol_6omn"
    MODEL_OBJ_NAME = "selected_prediction"
    # define a directory where all data is stored which needs to get imported
    IMPORT_SUBDIR = "data"
    # define a parent directory where all results will be saved
    EXPORT_SUBDIR = "data/results"

    # define parameters for the align command
    CYCLES = 1
    CUTOFF = 1.0
    # define a filename for the rmsd and aligned amino acids CSV file
    FILE_NAME = "rmsd_and_aligned_AA_form_alignment"

    # constant variables for plotting
    FIGURE_SIZE = (11.0, 6.0)
    # defines the atoms of which the distance gets computed
    DISTANCE_LABEL = "$\\alpha$-C"

    # define variable filenames for the different results
    ALIGNMENT_FILE_NAME = "alignment_file_%s" % MODEL_OBJ_NAME
    DISTANCE_FILE_NAME = "distance_file_%s" % MODEL_OBJ_NAME
    SESSION_FILE_NAME = "session_file_%s" % MODEL_OBJ_NAME

    # comment the below code in, if you wish to have a fixed y-axis for ALL plots
    # LIMIT_Y: tuple[float, float] = (0, 7)

    # create the Protein object for the reference
    reference_protein = core.Protein(REFERENCE_OBJ_NAME, IMPORT_SUBDIR)
    reference_selection = "/%s//G//CA, /%s//H//CA" % (
        reference_protein.molecule_object,
        reference_protein.molecule_object,
    )
    reference_protein.set_selection(reference_selection)

    # create model Protein object
    model_protein = core.Protein(MODEL_OBJ_NAME, IMPORT_SUBDIR)
    # set selection for the model, for the align command
    model_selection = "/%s//A//CA, /%s//B//CA" % (model_protein.molecule_object, model_protein.molecule_object)
    model_protein.set_selection(model_selection)
    # create Protein pair object
    bmp2_pair = core.ProteinPair(reference_protein, model_protein, EXPORT_SUBDIR)
    # reinitialize pymol session
    cmd.reinitialize()
    # load both proteins in the pymol session
    bmp2_pair.load_protein_pair()
    print("Finished loading proteins.")
    # color Protein with default colors; ref: green, model: blue
    bmp2_pair.color_protein_pair()
    print("Finished coloring proteins.")
    # do the structure alignment
    align_results = bmp2_pair.align_protein_pair(CYCLES, CUTOFF, ALIGNMENT_FILE_NAME)
    print("Finished aligning proteins.")
    # do the distance calculation
    distance_results = bmp2_pair.calculate_distance_between_ca_atoms(ALIGNMENT_FILE_NAME)
    bmp2_pair.export_distance_between_ca_atoms(distance_results)
    print("Finished distance calculations.")
    # create an instance of the Graphics class
    graphics_instance = graphics.Graphics(bmp2_pair, distance_results, FIGURE_SIZE)
    # create distance plot
    fig = graphics_instance.create_distance_plot(DISTANCE_LABEL, CUTOFF)
    print("Finished creation of distance plot.")
    # save distance plot
    if not os.path.exists("%s/plots" % bmp2_pair.results_dir):
        os.mkdir("%s/plots" % bmp2_pair.results_dir)
    if not os.path.exists("%s/plots/distance_plot" % bmp2_pair.results_dir):
        os.mkdir("%s/plots/distance_plot" % bmp2_pair.results_dir)
    plt.savefig("%s/plots/distance_plot/distance_plot_%s.svg" % (bmp2_pair.results_dir, MODEL_OBJ_NAME))
    plt.close(fig)

    # create distance histogram
    graphics_instance.create_distance_histogram()
    print("Finished creation of distance histogram.")
    # save distance histogram
    if not os.path.exists("%s/plots" % bmp2_pair.results_dir):
        os.mkdir("%s/plots" % bmp2_pair.results_dir)
    if not os.path.exists("%s/plots/distance_histogram" % bmp2_pair.results_dir):
        os.mkdir("%s/plots/distance_histogram" % bmp2_pair.results_dir)
    plt.savefig("%s/plots/distance_plot/distance_histogram_%s.svg" % (bmp2_pair.results_dir, MODEL_OBJ_NAME))

    # take image of whole structure alignment
    # graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME, "cartoon",
    #                                              "test")
    print("Finished creation of image which shows the whole structure " "alignment.")
    # take image of interesting regions
    # graphics_instance.take_image_of_interesting_regions(5.0,
    #                                                     "interesting_region",
    #                                                     opaque_background=1)
    # print(f"Finished creation of images which show interesting regions with"
    #       f"a distance greater than {CUTOFF} angstrom.")

    # color residues by distance
    graphics_instance.color_by_distance(ALIGNMENT_FILE_NAME)
    print("Finished coloring of prediction with color_by_distance.")
    # graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME, "cartoon",
    #                                              "coloredByRMSD")
    graphics_instance.create_gif()
    print("Finished creation of gif which shows the whole predicted " "structure colored by distance.")
    # save pymol session
    bmp2_pair.save_session_of_protein_pair(SESSION_FILE_NAME)
    print("Finished saving of pymol session file.")
    t1 = time.time()
    duration = round(t1 - t0, 2)
    print("The code took %f seconds." % duration)
