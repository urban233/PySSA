import os.path

import zmq
import colabfold_run


def read_log_file():
    log_filepath = "/home/rhel_user/scratch/local_predictions/pdb/log.txt"
    if os.path.exists(log_filepath):
        log_file = open(log_filepath)
        tmp_log_content: str = log_file.read()
    else:
        tmp_log_content: str = "No log file found."
    return tmp_log_content


def start_computation(fasta_dir: str, pdb_dir: str, use_amber: bool, use_templates: bool) -> dict:
    """Start structure prediction process based on custom arguments.

    Args:
        use_amber:  true if protein should be relaxed through amber
        use_templates: true if pdb70 templates should be used for prediction
    """
    try:
        colabfold_run.run_prediction(fasta_dir, pdb_dir, use_amber, use_templates)
    except Exception as e:
        return {"error": str(e), "log": read_log_file()}
    success = "Prediction finished without errors."
    return {"success": success, "log": read_log_file()}


if __name__ == "__main__":
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    # Bind to the WSL2 IP address and port
    socket.bind("tcp://127.0.0.1:7016")

    print("Server is ready to receive messages.")

    while True:
        # Wait for a message from the client
        message = socket.recv_json()

        # Perform the computation and send back the result
        result = start_computation(message["fasta_dir"], message["pdb_dir"], message["use_amber"], message["use_templates"])
        socket.send_json(result)
