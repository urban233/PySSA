import build_inno_setup


AUTOMATION_TREE = {
  "build": {
    "help": "Build targets",
    "subcommands": {
      "inno_setup": {
        "help": "Builds the Inno Setup EXE file.",
        "func": build_inno_setup.build_full_setup_exe()
      },
      "update_inno_setup": {
        "help": "Builds the update Inno Setup EXE file.",
        "func": build_inno_setup.build_update_setup_exe()
      }
    }
  }
}


if __name__ == '__main__':
  from task_automator import automator
  automator.Automator(AUTOMATION_TREE).run()
