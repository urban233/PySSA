#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Contains workers that are handling the requests to the auxiliary PyMOL."""
import queue
import zmq
import auxiliary_pymol_base
import local_enums
__docformat__ = "google"


def handle_request(a_job_type: str, a_queue: queue.Queue, a_socket: zmq.Socket) -> bool:
    """Handles a request to the auxiliary PyMOL server.
    
    Args:
        a_job_type (str): A string representing the type of job to handle.
        a_queue (queue.Queue): A queue.Queue object representing the queue of tasks to handle.
        a_socket (zmq.Socket): A zmq.Socket object representing the socket used for communication.
    
    Returns:
        A boolean value indicating success or failure.
    """
    # <editor-fold desc="Checks">
    if a_job_type is None:
        return False
    if a_queue is None:
        print("a_queue is None.")
        return False
    if a_socket is None:
        return False
    
    # </editor-fold>
    
    tmp_queue_is_empty = False
    try:
        while tmp_queue_is_empty is False:
            print(f"Checking queue for job type {a_job_type}...")
            # Get a task from the queue
            data = a_queue.get()
            try:
                if a_job_type == local_enums.JobType.PREDICTION.value:
                    base64_string = auxiliary_pymol_base.AuxiliaryPyMOL.create_pymol_session_for_protein(
                        data[local_enums.JobDescriptionKeys.PDB_FILEPATH.value],
                    )
                    print("Structure prediction finished.")
                    response = {"result": "", "data": (base64_string, )}
                elif a_job_type == local_enums.JobType.DISTANCE_ANALYSIS.value:
                    print("Starting distance analysis ...")
                    distance_analysis_results_object_values, base64_string = auxiliary_pymol_base.AuxiliaryPyMOL.do_distance_analysis(
                        data["the_protein_pair_name"],
                        data["a_protein_1_pdb_cache_filepath"],
                        data["a_protein_2_pdb_cache_filepath"],
                        data["a_protein_1_pymol_selection_string"],
                        data["a_protein_2_pymol_selection_string"],
                        data["a_cutoff"],
                        data["the_cycles"],
                    )
                    print("Distance analysis finished.")
                    response = {"result": "", "data": (distance_analysis_results_object_values, base64_string)}
                elif a_job_type == local_enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS.value:
                    pass
                elif a_job_type == local_enums.JobType.RAY_TRACING.value:
                    auxiliary_pymol_base.AuxiliaryPyMOL.create_ray_traced_image(
                        data[local_enums.JobDescriptionKeys.IMAGE_DESTINATION_FILEPATH.value],
                        data[local_enums.JobDescriptionKeys.CACHED_SESSION_FILEPATH.value],
                        data[local_enums.JobDescriptionKeys.RAY_TRACE_MODE.value],
                        data[local_enums.JobDescriptionKeys.RAY_TEXTURE.value],
                        data[local_enums.JobDescriptionKeys.RAY_TRACING_RENDERER.value],
                    )
                    response = {"result": "", "data": (data[
                                                           local_enums.JobDescriptionKeys.IMAGE_DESTINATION_FILEPATH.value],)}
                elif a_job_type == local_enums.JobType.GENERAL_PURPOSE.value:
                    if data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.CREATE_NEW_PROTEIN_PYMOL_SESSION.value:
                        # Create a new protein pymol session
                        base64_string = auxiliary_pymol_base.AuxiliaryPyMOL.create_pymol_session_for_protein(
                            data[local_enums.JobDescriptionKeys.PDB_FILEPATH.value],
                        )
                        print("Create PyMOL protein session finished.")
                        response = {"result": "", "data": (base64_string,)}
    
                    elif data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.GET_ALL_CHAINS_OF_GIVEN_PROTEIN.value:
                        tmp_chain_object_values = auxiliary_pymol_base.AuxiliaryPyMOL.get_protein_chains(
                            data[local_enums.JobDescriptionKeys.PDB_FILEPATH.value],
                        )
                        print("Get all chains of protein finished.")
                        response = {"result": "", "data": tmp_chain_object_values}
    
                    elif data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.CONSOLIDATE_MOLECULE_OBJECT_TO_FIRST_STATE.value:
                        tmp_new_pdb_filepath = auxiliary_pymol_base.AuxiliaryPyMOL.consolidate_molecule_object_to_first_state(
                            data[local_enums.JobDescriptionKeys.PDB_FILEPATH.value],
                        )
                        print("Consolidation of molecule object states finished.")
                        response = {"result": "", "data": tmp_new_pdb_filepath}
    
                    elif data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.GET_ALL_SCENES_OF_SESSION.value:
                        tmp_all_scenes = auxiliary_pymol_base.AuxiliaryPyMOL.get_all_scenes_of_session(
                            data[local_enums.JobDescriptionKeys.PYMOL_SESSION.value],
                        )
                        print("Get all scenes of pymol session finished.")
                        response = {"result": "", "data": tmp_all_scenes}
    
                    elif data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.CLEAN_PROTEIN_UPDATE_STRUCTURE.value:
                        tmp_pymol_session, tmp_new_pdb_filepath = auxiliary_pymol_base.AuxiliaryPyMOL.clean_protein_update(
                            data[local_enums.JobDescriptionKeys.PYMOL_SESSION.value], data[
                                local_enums.JobDescriptionKeys.PROTEIN_NAME.value],
                        )
                        print("Cleaning given protein finished.")
                        response = {"result": "", "data": (tmp_pymol_session, tmp_new_pdb_filepath)}
            except AttributeError as e:
                if str(e).find("Error: Selection 1: More than one atom found") != -1:
                    response = {"result": "error", "data": "Malformed pdb file."}
                else:
                    response = {"result": "error", "data": str(e)}
            except Exception as e:
                response = {"result": "error", "data": str(e)}
            # Indicate that the task is done
            a_queue.task_done()
            # Send the response back to the client
            #print("This is the response: {}".format(response))
            a_socket.recv_json()
            a_socket.send_json(response)  # Send JSON-encoded response
    except Exception as e:
        print(e)
        return False
    else:
        return True
