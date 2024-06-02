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
"""Module to start Auxiliary PyMOL."""
import threading
import queue
import workers
import zmq
import local_enums
__docformat__ = "google"


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
                local_enums.JobType.PREDICTION.value, prediction_queue, prediction_socket,
            ),
        )
        prediction_thread.daemon = True
        prediction_thread.start()

        distance_analysis_queue = queue.Queue()
        distance_analysis_socket = context.socket(zmq.REP)
        distance_analysis_socket.bind("tcp://127.0.0.1:8072")
        distance_analysis_thread = threading.Thread(
            target=workers.handle_request,
            args=(
                local_enums.JobType.DISTANCE_ANALYSIS.value, distance_analysis_queue, distance_analysis_socket,
            ),
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
                prediction_and_distance_analysis_socket,
            ),
        )
        prediction_and_distance_analysis_thread.daemon = True
        prediction_and_distance_analysis_thread.start()

        ray_tracing_queue = queue.Queue()
        ray_tracing_socket = context.socket(zmq.REP)
        ray_tracing_socket.bind("tcp://127.0.0.1:8074")
        ray_tracing_thread = threading.Thread(
            target=workers.handle_request,
            args=(
                local_enums.JobType.RAY_TRACING.value, ray_tracing_queue, ray_tracing_socket,
            ),
        )
        ray_tracing_thread.daemon = True
        ray_tracing_thread.start()

        general_purpose_queue = queue.Queue()
        general_purpose_socket = context.socket(zmq.REP)
        general_purpose_socket.bind("tcp://127.0.0.1:8075")
        general_purpose_thread = threading.Thread(
            target=workers.handle_request,
            args=(
                local_enums.JobType.GENERAL_PURPOSE.value, general_purpose_queue, general_purpose_socket,
            ),
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
            local_enums.JobType.GENERAL_PURPOSE.value: general_purpose_queue,
        }
        sockets = [prediction_socket, distance_analysis_socket, prediction_and_distance_analysis_socket, ray_tracing_socket, general_purpose_socket]
        threads = [prediction_thread, distance_analysis_thread, prediction_and_distance_analysis_thread, ray_tracing_thread, general_purpose_thread]
    except Exception as e:
        print("Exception raised during startup!")
        print(e)
        exit(1)

    tmp_auxiliary_pymol_should_be_closed = False
    try:
        while tmp_auxiliary_pymol_should_be_closed is False:
            # Wait for a request from the client
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
                tmp_auxiliary_pymol_should_be_closed = True
            else:
                queues[task_type].put(data)
                main_socket.send_string("Received data and added job to queue.")
    except Exception as e:
        print("Caught general exception.")
        response = {"result": "error", "data": str(e)}
        main_socket.send_json(response)  # Send JSON-encoded response
    exit(0)
