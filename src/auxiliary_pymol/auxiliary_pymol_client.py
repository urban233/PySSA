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
"""Module for the interaction with the auxiliary pymol."""
from typing import Union
from typing import TYPE_CHECKING

import zmq

from src.pyssa.util import enums

__docformat__ = "google"

if TYPE_CHECKING:
    from src.pyssa.internal.data_structures import job


def send_request_to_auxiliary_pymol(
        the_main_socket: zmq.Socket,
        a_socket: zmq.Socket,
        a_job_description: Union["job.PredictionJobDescription", "job.DistanceAnalysisJobDescription", "job.RayTracingJobDescription", "job.GeneralPurposeJobDescription"],
) -> dict:
    """Sends a request to the auxiliary PyMOL socket and receives the results.

    Args:
        the_main_socket (zmq.Socket): The main socket for communication.
        a_socket (zmq.Socket): The auxiliary socket for communication.
        a_job_description (Union[job.PredictionJobDescription, job.DistanceAnalysisJobDescription, job.RayTracingJobDescription, job.GeneralPurposeJobDescription]): The job description.

    Returns:
        A dict of results from the auxiliary socket or an empty dict if an argument is None.
    """
    # <editor-fold desc="Checks">
    if the_main_socket is None:
        print("the_main_socket is None.")  # TODO: substitute with logger
        return {}
    if a_socket is None:
        print("a_socket is None.")  # TODO: substitute with logger
        return {}
    if a_job_description is None:
        print("a_job_description is None.")  # TODO: substitute with logger
        return {}
    
    # </editor-fold>
    
    # First check if main socket can receive messages (fixme: There is a better way to test this)
    the_main_socket.send_string("Check availability ...")
    the_main_socket.recv_string()
    the_main_socket.send_json(a_job_description.get_dict())
    the_main_socket.recv_string()
    # Send a json message to activate the socket
    if a_job_description.type == enums.JobType.PREDICTION:
        a_socket.send_json(
            {
                "job_type": enums.JobType.PREDICTION.value,
                "job_short_description": enums.JobShortDescription.RUN_STRUCTURE_PREDICTION.value,
            },
        )
    elif a_job_description.type == enums.JobType.DISTANCE_ANALYSIS:
        a_socket.send_json(
            {
                "job_type": enums.JobType.DISTANCE_ANALYSIS.value,
                "job_short_description": enums.JobShortDescription.RUN_DISTANCE_ANALYSIS.value,
            },
        )
    elif a_job_description.type == enums.JobType.RAY_TRACING:
        a_socket.send_json(
            {
                "job_type": enums.JobType.RAY_TRACING.value,
                "job_short_description": enums.JobShortDescription.CREATE_RAY_TRACED_IMAGE.value,
            },
        )
    elif a_job_description.type == enums.JobType.GENERAL_PURPOSE:
        a_socket.send_json(
            {
                "job_type": enums.JobType.GENERAL_PURPOSE.value,
                "job_short_description": "generic",
            },
        )

    # Receive results
    return a_socket.recv_json()
