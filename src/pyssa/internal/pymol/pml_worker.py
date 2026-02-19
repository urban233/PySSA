from multiprocessing.connection import PipeConnection

import pymol2
from src.pyssa.gui.user_pymol import UserPyMOL
from src.pyssa.internal.pymol import worker_command


def start_pml_worker(pipe_connection: PipeConnection):
  tmp_pymol = pymol2.PyMOL()
  tmp_pymol.start()
  while not pipe_connection.closed:
    tmp_worker_command: worker_command.WorkerCommand = pipe_connection.recv()
    if tmp_worker_command.command == "shutdown":
      break
    print(f"Processing: {tmp_worker_command.command}")
    tmp_pymol.cmd.load(tmp_worker_command.session_path)
    if tmp_worker_command.command == "get_model":
      tmp_output = tmp_pymol.cmd.get_model("all")
    else:
      tmp_pymol.cmd.do(tmp_worker_command.get_command_with_args())
      tmp_output = f"Processed: {tmp_worker_command.command}"
    # pipe_connection.send(f"Processing: {tmp_worker_command.command}")
    if tmp_worker_command.sync:
      pipe_connection.send(tmp_output)
  tmp_pymol.stop()


def cache_session(user_pymol: UserPyMOL, session_filepath):
  user_pymol.get_cmd_module().save(session_filepath)
