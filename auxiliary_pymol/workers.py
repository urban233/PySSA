import queue
import auxiliary_pymol_base
import local_enums


def handle_request(job_type, a_queue: queue.Queue, socket):
    while True:
        print(f"Checking queue for job type {job_type}...")
        # Get a task from the queue
        data = a_queue.get()
        try:
            if job_type == local_enums.JobType.PREDICTION.value:
                base64_string = auxiliary_pymol_base.AuxiliaryPyMOL.create_pymol_session_for_protein(
                    data[local_enums.JobDescriptionKeys.PDB_FILEPATH.value]
                )
                print("Structure prediction finished.")
                response = {"result": "", "data": (base64_string, )}
            elif job_type == local_enums.JobType.DISTANCE_ANALYSIS.value:
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
            elif job_type == local_enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS.value:
                pass
            elif job_type == local_enums.JobType.RAY_TRACING.value:
                auxiliary_pymol_base.AuxiliaryPyMOL.create_ray_traced_image(
                    data[local_enums.JobDescriptionKeys.IMAGE_DESTINATION_FILEPATH.value],
                    data[local_enums.JobDescriptionKeys.CACHED_SESSION_FILEPATH.value],
                    data[local_enums.JobDescriptionKeys.RAY_TRACE_MODE.value],
                    data[local_enums.JobDescriptionKeys.RAY_TEXTURE.value],
                    data[local_enums.JobDescriptionKeys.RAY_TRACING_RENDERER.value],
                )
                response = {"result": "", "data": (data[local_enums.JobDescriptionKeys.IMAGE_DESTINATION_FILEPATH.value],)}
            elif job_type == local_enums.JobType.GENERAL_PURPOSE.value:
                if data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.CREATE_NEW_PROTEIN_PYMOL_SESSION.value:
                    # Create a new protein pymol session
                    base64_string = auxiliary_pymol_base.AuxiliaryPyMOL.create_pymol_session_for_protein(
                        data[local_enums.JobDescriptionKeys.PDB_FILEPATH.value]
                    )
                    print("Create PyMOL protein session finished.")
                    response = {"result": "", "data": (base64_string,)}

                elif data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.GET_ALL_CHAINS_OF_GIVEN_PROTEIN.value:
                    tmp_chain_object_values = auxiliary_pymol_base.AuxiliaryPyMOL.get_protein_chains(
                        data[local_enums.JobDescriptionKeys.PDB_FILEPATH.value]
                    )
                    print("Get all chains of protein finished.")
                    response = {"result": "", "data": tmp_chain_object_values}

                elif data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.CONSOLIDATE_MOLECULE_OBJECT_TO_FIRST_STATE.value:
                    tmp_new_pdb_filepath = auxiliary_pymol_base.AuxiliaryPyMOL.consolidate_molecule_object_to_first_state(
                        data[local_enums.JobDescriptionKeys.PDB_FILEPATH.value]
                    )
                    print("Consolidation of molecule object states finished.")
                    response = {"result": "", "data": tmp_new_pdb_filepath}

                elif data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.GET_ALL_SCENES_OF_SESSION.value:
                    tmp_all_scenes = auxiliary_pymol_base.AuxiliaryPyMOL.get_all_scenes_of_session(
                        data[local_enums.JobDescriptionKeys.PYMOL_SESSION.value]
                    )
                    print("Get all scenes of pymol session finished.")
                    response = {"result": "", "data": tmp_all_scenes}

                elif data[local_enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value] == local_enums.JobShortDescription.CLEAN_PROTEIN_UPDATE_STRUCTURE.value:
                    tmp_pymol_session, tmp_new_pdb_filepath = auxiliary_pymol_base.AuxiliaryPyMOL.clean_protein_update(
                        data[local_enums.JobDescriptionKeys.PYMOL_SESSION.value]
                    )
                    print("Cleaning given protein finished.")
                    response = {"result": "", "data": (tmp_pymol_session, tmp_new_pdb_filepath)}
        except Exception as e:
            response = {"result": str(e), "data": ""}
        # Indicate that the task is done
        a_queue.task_done()
        # Send the response back to the client
        print("This is the response: {}".format(response))
        socket.recv_json()
        socket.send_json(response)  # Send JSON-encoded response
