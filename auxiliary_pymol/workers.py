import queue
import auxiliary_pymol_base


def handle_request(job_type, a_queue: queue.Queue, socket):
    while True:
        print(f"Checking queue for job type {job_type}...")
        # Get a task from the queue
        data = a_queue.get()
        print("Got data, start with rendering job.")
        try:
            if job_type == "Structure Prediction":
                base64_string = auxiliary_pymol_base.AuxiliaryPyMOL.create_pymol_session_for_protein(
                    data["a_pdb_filepath"]
                )
                print("Structure prediction finished.")
                response = {"result": "", "data": (base64_string, )}
            elif job_type == "Distance Analysis":
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
            elif job_type == "Structure Prediction and Distance Analysis":
                pass
            elif job_type == "Ray-tracing":
                auxiliary_pymol_base.AuxiliaryPyMOL.create_ray_traced_image(
                    data["dest"],
                    data["cached"],
                    data["mode"],
                    data["texture"],
                    data["renderer"]
                )
                response = {"result": "", "data": (data["dest"],)}
        except Exception as e:
            response = {"result": str(e), "data": ""}
        # Indicate that the task is done
        a_queue.task_done()
        # Send the response back to the client
        print("This is the response: {}".format(response))
        socket.recv_json()
        socket.send_json(response)  # Send JSON-encoded response
