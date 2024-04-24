from typing import Union
from typing import TYPE_CHECKING

from pyssa.util import enums

if TYPE_CHECKING:
    from pyssa.internal.data_structures import job


def send_request_to_auxiliary_pymol(
        the_main_socket,
        a_socket,
        a_job_description: Union["job.PredictionJobDescription", "job.DistanceAnalysisJobDescription", "job.RayTracingJobDescription", "job.GeneralPurposeJobDescription"]
) -> dict:
    """Sends a job request to the auxiliary pymol process."""
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
            }
        )
    elif a_job_description.type == enums.JobType.DISTANCE_ANALYSIS:
        a_socket.send_json(
            {
                "job_type": enums.JobType.DISTANCE_ANALYSIS.value,
                "job_short_description": enums.JobShortDescription.RUN_DISTANCE_ANALYSIS.value,
            }
        )
    elif a_job_description.type == enums.JobType.RAY_TRACING:
        a_socket.send_json(
            {
                "job_type": enums.JobType.RAY_TRACING.value,
                "job_short_description": enums.JobShortDescription.CREATE_RAY_TRACED_IMAGE.value,
            }
        )
    elif a_job_description.type == enums.JobType.GENERAL_PURPOSE:
        a_socket.send_json(
            {
                "job_type": enums.JobType.GENERAL_PURPOSE.value,
                "job_short_description": "generic",
            }
        )

    # Receive results
    return a_socket.recv_json()
