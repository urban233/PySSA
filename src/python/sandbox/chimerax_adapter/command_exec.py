import zmq


class ChimeraXCommandExec:
  def __init__(self):
    pass

  @staticmethod
  def run_single_command(a_socket, command: str, args: str):
    a_socket.send_json({
      "command": command,
      "arguments": args,
    })
