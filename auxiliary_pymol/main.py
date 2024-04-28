import threading
import queue
import workers
import zmq
import local_enums


if __name__ == "__main__":
    try:
        context = zmq.Context()
        main_socket = context.socket(zmq.REP)
        main_socket.bind("tcp://127.0.0.1:8070")  # Binding to localhost on port 7070

        # <editor-fold desc="Queue and socket definitions">
        prediction_queue = queue.Queue()
        prediction_socket = context.socket(zmq.REP)
        prediction_socket.bind("tcp://127.0.0.1:8071")
        prediction_thread = threading.Thread(
            target=workers.handle_request,
            args=(
                local_enums.JobType.PREDICTION.value, prediction_queue, prediction_socket
            )
        )
        prediction_thread.daemon = True
        prediction_thread.start()

        distance_analysis_queue = queue.Queue()
        distance_analysis_socket = context.socket(zmq.REP)
        distance_analysis_socket.bind("tcp://127.0.0.1:8072")
        distance_analysis_thread = threading.Thread(
            target=workers.handle_request,
            args=(
                local_enums.JobType.DISTANCE_ANALYSIS.value, distance_analysis_queue, distance_analysis_socket
            )
        )
        distance_analysis_thread.daemon = True
        distance_analysis_thread.start()

        prediction_and_distance_analysis_queue = queue.Queue()
        prediction_and_distance_analysis_socket = context.socket(zmq.REP)
        prediction_and_distance_analysis_socket.bind("tcp://127.0.0.1:8073")
        prediction_and_distance_analysis_thread = threading.Thread(
            target=workers.handle_request,
            args=(
                local_enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS.value,
                prediction_and_distance_analysis_queue,
                prediction_and_distance_analysis_socket
            )
        )
        prediction_and_distance_analysis_thread.daemon = True
        prediction_and_distance_analysis_thread.start()

        ray_tracing_queue = queue.Queue()
        ray_tracing_socket = context.socket(zmq.REP)
        ray_tracing_socket.bind("tcp://127.0.0.1:8074")
        ray_tracing_thread = threading.Thread(
            target=workers.handle_request,
            args=(
                local_enums.JobType.RAY_TRACING.value, ray_tracing_queue, ray_tracing_socket
            )
        )
        ray_tracing_thread.daemon = True
        ray_tracing_thread.start()

        general_purpose_queue = queue.Queue()
        general_purpose_socket = context.socket(zmq.REP)
        general_purpose_socket.bind("tcp://127.0.0.1:8075")
        general_purpose_thread = threading.Thread(
            target=workers.handle_request,
            args=(
                local_enums.JobType.GENERAL_PURPOSE.value, general_purpose_queue, general_purpose_socket
            )
        )
        general_purpose_thread.daemon = True
        general_purpose_thread.start()
        # </editor-fold>
        # Group queues and sockets
        queues = {
            local_enums.JobType.PREDICTION.value: prediction_queue,
            local_enums.JobType.DISTANCE_ANALYSIS.value: distance_analysis_queue,
            local_enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS.value: prediction_and_distance_analysis_queue,
            local_enums.JobType.RAY_TRACING.value: ray_tracing_queue,
            local_enums.JobType.GENERAL_PURPOSE.value: general_purpose_queue
        }
        sockets = [prediction_socket, distance_analysis_socket, prediction_and_distance_analysis_socket, ray_tracing_socket, general_purpose_socket]
        threads = [prediction_thread, distance_analysis_thread, prediction_and_distance_analysis_thread, ray_tracing_thread, general_purpose_thread]
        # Bind sockets
        # i = 0
        # for tmp_socket in sockets:
        #     port = 7071 + i
        #     tmp_socket.bind(f"tcp://127.0.0.1:{port}")
        #     i += 1
        # for tmp_thread in threads:
        #     tmp_thread.daemon = True
        #     tmp_thread.start()
    except Exception as e:
        print("Exception raised during startup!")
        print(e)
        exit(1)

    print("Entering while loop ...")
    while True:
        # Wait for a request from the client
        try:
            print("Waiting to receive message ...")
            message = main_socket.recv_string()
            print("Received message: ", message)
            main_socket.send_string("I am ready to receive data.")
            print("Waiting for data from the client ...")
            data = main_socket.recv_json()
            print(data)
            task_type = data["job_type"]
            print(task_type)
            if task_type == "Abort":
                main_socket.send_string("Auxiliary PyMOL will now exit.")
                exit(0)
            else:
                queues[task_type].put(data)
                main_socket.send_string("Received data and added job to queue.")
        except Exception as e:
            response = {"result": str(e)}
            main_socket.send_json(response)  # Send JSON-encoded response


    # while True:
    #     # Wait for a request from the client
    #     message = socket.recv_string()
    #     socket.send_string("I am ready to receive data.")
    #     try:
    #         data = socket.recv_json()  # Receive JSON-encoded data
    #         """
    #         message = {
    #             "message": "Ray-tracing",
    #             "dest": str(self.dest_image_filepath),
    #             "cached": str(self.cached_session_filepath),
    #             "mode": self.image_ray_trace_mode,
    #             "texture": self.image_ray_texture,
    #             "renderer": self.image_renderer
    #         }
    #         """
    #         if data["message"] == "Ray-tracing":
    #             auxiliary_pymol_base.AuxiliaryPyMOL.create_ray_traced_image(
    #                 data["dest"],
    #                 data["cached"],
    #                 data["mode"],
    #                 data["texture"],
    #                 data["renderer"]
    #             )
    #         # Process the arguments (here, just echoing them back)
    #         response = {"result": ""}
    #     except Exception as e:
    #         response = {"result": e}
    #     # Send the response back to the client
    #     socket.send_json(response)  # Send JSON-encoded response
