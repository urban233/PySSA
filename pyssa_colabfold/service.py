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
"""Module that runs the colabfold micro service."""
import os.path
import shutil
import zmq
import colabfold_run
import subprocess


def change_ownership_recursive(directory_path):
    try:
        # Run the chown command recursively on the specified directory
        subprocess.run(["sudo", "chown", "-R", "rhel_user", directory_path], check=True)
        print(f"Ownership changed successfully for {directory_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")


def delete_scratch_directory_in_wsl2() -> bool:
    """Deletes the scratch directory in WSL2.

    Raises:
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    tmp_scratch_path: str = "/home/rhel_user/scratch"
    try:
        if os.path.exists(tmp_scratch_path):
            shutil.rmtree(tmp_scratch_path)
            return True
    except Exception as e:
        raise OSError(f"Scratch directory could not be deleted! {e}")
    return False


def create_fasta_directory_in_wsl2() -> None:
    """Creates the fasta directory inside the WSL2.

    Raises:
        IllegalArgumentError: If the argument is None.
        SubprocessExecutionError: If return code of subprocess is non-zero.
    """
    try:
        os.makedirs("/home/rhel_user/scratch/local_predictions/fasta")
    except Exception as e:
        raise OSError(f"Fasta directory could not be created! {e}")


def create_pdb_directory_in_wsl2() -> None:
    """Creates the pdb directory inside the WSL2.

    Raises:
        IllegalArgumentError: If the argument is None.
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    try:
        os.makedirs("/home/rhel_user/scratch/local_predictions/pdb")
    except Exception as e:
        raise OSError(f"Pdb directory could not be created! {e}")


def copy_fasta_files_from_windows_to_wsl2(the_fasta_path: str) -> None:
    """Copies fasta files from Windows host to WSL2.

    Raises:
        IllegalArgumentError: If an argument is None.
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    try:
        shutil.copytree(f"{the_fasta_path}", "/home/rhel_user/scratch/local_predictions/fasta/")
    except Exception as e:
        raise OSError(f"Fasta files from Windows host could not be copied to WSL2! {e}")


def disable_cuda_device_usage():
    import os
    os.environ['CUDA_VISIBLE_DEVICES'] = ''


def prepare_computation_environment(the_scratch_fasta_dir_of_windows_host: str) -> None:
    """Prepares the environment needed for the computation inside the WSL2."""
    str_conversion_1 = the_scratch_fasta_dir_of_windows_host.replace("\\", "/")
    str_conversion_2 = str_conversion_1.replace(":", "")
    str_conversion_3 = str_conversion_2.replace("C", "c")
    the_scratch_fasta_dir_of_windows_host = f"/mnt/{str_conversion_3}"
    try:
        disable_cuda_device_usage()
        change_ownership_recursive("/home/rhel_user/")
        delete_scratch_directory_in_wsl2()
        create_pdb_directory_in_wsl2()
        change_ownership_recursive("/home/rhel_user/scratch")
        copy_fasta_files_from_windows_to_wsl2(the_scratch_fasta_dir_of_windows_host)
    except OSError as e:
        raise EnvironmentError(e)


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
        prepare_computation_environment(fasta_dir)
        tmp_fasta_dir: str = "/home/rhel_user/scratch/local_predictions/fasta"
        tmp_pdb_dir: str = "/home/rhel_user/scratch/local_predictions/pdb"
        colabfold_run.run_prediction(tmp_fasta_dir, tmp_pdb_dir, use_amber, use_templates)
    except RuntimeError as e:
        return {"error": str(e), "log": read_log_file(), "exit_code": 2}
    except Exception as e:
        return {"error": str(e), "log": read_log_file(), "exit_code": 1}
    success = "Prediction finished without errors."
    return {"success": success, "log": read_log_file(), "exit_code": 0}


if __name__ == "__main__":
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    # Bind to the WSL2 IP address and port
    tmp_port = 7016
    socket.bind(f"tcp://127.0.0.1:{str(tmp_port)}")
    print("Server is ready to receive messages.")

    while True:
        # Wait for a message from the client
        message = socket.recv_json()

        # Perform the computation and send back the result
        result = start_computation(
            message["fasta_dir"], message["pdb_dir"], message["use_amber"], message["use_templates"],
        )
        socket.send_json(result)


# def delete_original_batch_py_file() -> None:
#     """Deletes the original batch.py file of colabfold.
#
#     Raises:
#         SubprocessExecutionError: If return code of subprocess is non-zero
#     """
#     # Remove original batch.py of Colabfold
#     tmp_batch_py_filepath: str = (
#         "/home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages/colabfold/batch.py"
#     )
#     try:
#         if os.path.exists(tmp_batch_py_filepath):
#             os.remove(tmp_batch_py_filepath)
#     except Exception as e:
#         raise OSError(f"Original batch.py file could not be deleted! {e}")
#
#
# def copy_modified_batch_py_file() -> None:
#     """Copies the modified batch.py file to the WSL2.
#
#     Raises:
#         SubprocessExecutionError: If return code of subprocess is non-zero
#     """
#     plugin_path = "/mnt/c/ProgramData/pyssa/mambaforge_pyssa/pyssa-mamba-env/Lib/site-packages/pymol/pymol_path/data/startup/PySSA"
#     tmp_batch_py_filepath_wsl: str = (
#         "/home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages/colabfold/batch.py"
#     )
#     tmp_batch_py_filepath_windows: str = f"{plugin_path}/pyssa_colabfold/colabfold_sub/batch.py"
#     try:
#         subprocess.run(["cp", tmp_batch_py_filepath_windows, tmp_batch_py_filepath_wsl])
#     except Exception as e:
#         raise OSError(f"Modified batch.py file could not be copied! {e}")