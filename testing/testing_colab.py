import os
import sys
import time
from pymolproteintools import core
from pymolproteintools import graphics
from pymol import cmd
import matplotlib.pyplot as plt

# if upload_pdb == True:
#     # define reference protein if a .pdb was uploaded
#     REFERENCE_OBJ_NAME: str = ref_name
# else:
#     # define reference protein if a pdb id was given
#     REFERENCE_OBJ_NAME: str = pdb_id
#     cmd.set("fetch_type_default", "pdb")
#     cmd.fetch(REFERENCE_OBJ_NAME)
#     # move the fetched pdb file to the data directory
#     src_path: str = f"/content/{pdb_name}"
#     dest_path: str = f"/content/data/{pdb_name}"
#     shutil.move(src_path, dest_path)

if __name__ == "__main__":
    # import of reference protein name from command line option
    REFERENCE_OBJ_NAME = sys.argv[1]

    t0 = time.time()
    REFERENCE_OBJ_NAME = "pymol_6omn"
    MODEL_OBJ_NAME: str = "selected_prediction"
    # define a directory where all data is stored which needs to get imported
    IMPORT_SUBDIR: str = "data"
    # define a parent directory where all results will be saved
    EXPORT_SUBDIR: str = "data/results"

    # define parameters for the align command
    CYCLES: int = 1
    CUTOFF: float = 1.0
    # define a filename for the rmsd and aligned amino acids CSV file
    FILE_NAME: str = "rmsd_and_aligned_AA_form_alignment"

    # constant variables for plotting
    FIGURE_SIZE = (11.0, 6.0)
    # defines the atoms of which the distance gets computed
    DISTANCE_LABEL: str = "$\\alpha$-C"

    # define variable filenames for the different results
    ALIGNMENT_FILE_NAME:str = f"alignment_file_{MODEL_OBJ_NAME}"
    DISTANCE_FILE_NAME:str = f"distance_file_{MODEL_OBJ_NAME}"
    SESSION_FILE_NAME: str = f"session_file_{MODEL_OBJ_NAME}"

    # comment the below code in, if you wish to have a fixed y-axis for ALL plots
    # LIMIT_Y: tuple[float, float] = (0, 7)

    # create the protein object for the reference
    reference_protein: core.protein = core.protein(REFERENCE_OBJ_NAME,
                                                   IMPORT_SUBDIR)
    reference_protein.set_selection(f"/{reference_protein.molecule_object}//G//CA, "
                                    f"/{reference_protein.molecule_object}//H//CA")

    # create model protein object
    model_protein: core.protein = core.protein(MODEL_OBJ_NAME, IMPORT_SUBDIR)
    # set selection for the model, for the align command
    model_protein.set_selection(f"/{model_protein.molecule_object}//A//CA, "
                                f"/{model_protein.molecule_object}//B//CA")
    # create protein pair object
    bmp2_pair: core.ProteinPair = core.ProteinPair(reference_protein,
                                                   model_protein,
                                                   EXPORT_SUBDIR)
    # reinitialize pymol session
    cmd.reinitialize()
    # load both proteins in the pymol session
    bmp2_pair.load_protein_pair()
    print(f"Finished loading proteins.")
    # color protein with default colors; ref: green, model: blue
    bmp2_pair.color_protein_pair()
    print(f"Finished coloring proteins.")
    # do the structure alignment
    align_results = bmp2_pair.align_protein_pair(CYCLES, CUTOFF,
                                                 ALIGNMENT_FILE_NAME)
    print(f"Finished aligning proteins.")
    # do the distance calculation
    distance_results = bmp2_pair.calculate_distance_between_ca_atoms(ALIGNMENT_FILE_NAME)
    bmp2_pair.export_distance_between_ca_atoms(distance_results)
    print(f"Finished distance calculations.")
    # create an instance of the graphics class
    graphics_instance: graphics.graphics = graphics.graphics(bmp2_pair,
                                                             distance_results,
                                                             FIGURE_SIZE)
    # create distance plot
    fig = graphics_instance.create_distance_plot(DISTANCE_LABEL, CUTOFF)
    print(f"Finished creation of distance plot.")
    # save distance plot
    if not os.path.exists(f"{bmp2_pair.results_dir}/plots"):
        os.mkdir(f"{bmp2_pair.results_dir}/plots")
    if not os.path.exists(f"{bmp2_pair.results_dir}/plots/distance_plot"):
        os.mkdir(f"{bmp2_pair.results_dir}/plots/distance_plot")
    plt.savefig(f"{bmp2_pair.results_dir}/plots/distance_plot/"
                f"distance_plot_{MODEL_OBJ_NAME}.svg")
    plt.close(fig)

    # create distance histogram
    graphics_instance.create_distance_histogram()
    print(f"Finished creation of distance histogram.")
    # save distance histogram
    if not os.path.exists(f"{bmp2_pair.results_dir}/plots"):
        os.mkdir(f"{bmp2_pair.results_dir}/plots")
    if not os.path.exists(f"{bmp2_pair.results_dir}/plots/distance_histogram"):
        os.mkdir(f"{bmp2_pair.results_dir}/plots/distance_histogram")
    plt.savefig(f"{bmp2_pair.results_dir}/plots/distance_histogram"
                f"/distance_histogram_{MODEL_OBJ_NAME}.svg")

    # take image of whole structure alignment
    # graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME, "cartoon",
    #                                              "test")
    print(f"Finished creation of image which shows the whole structure "
          f"alignment.")
    # take image of interesting regions
    # graphics_instance.take_image_of_interesting_regions(5.0,
    #                                                     "interesting_region",
    #                                                     opaque_background=1)
    # print(f"Finished creation of images which show interesting regions with"
    #       f"a distance greater than {CUTOFF} angstrom.")

    # color residues by distance
    graphics_instance.color_by_distance(ALIGNMENT_FILE_NAME)
    print(f"Finished coloring of prediction with color_by_distance.")
    # graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME, "cartoon",
    #                                              "coloredByRMSD")
    graphics_instance.create_gif()
    print(f"Finished creation of gif which shows the whole predicted "
          f"structure colored by distance.")
    # save pymol session
    bmp2_pair.save_session_of_protein_pair(SESSION_FILE_NAME)
    print(f"Finished saving of pymol session file.")
    t1 = time.time()
    duration = t1 - t0
    print(f"The code took {round(duration,2)} seconds.")
