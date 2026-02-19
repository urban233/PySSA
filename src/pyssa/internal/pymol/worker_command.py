from dataclasses import dataclass


@dataclass
class WorkerCommand:
    session_path: str
    command: str
    args: tuple
    sync: bool = False

    def get_command_with_args(self):
        tmp_args = ""
        if len(self.args) == 1:
            tmp_args = str(self.args[0])
        else:
            for tmp_arg in self.args:
                tmp_args += f"{tmp_arg}, "
        return f"{self.command} {tmp_args}"
