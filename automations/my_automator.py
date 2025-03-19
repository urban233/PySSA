import build_inno_setup
import make_end_user_docs


AUTOMATION_TREE = {
  "setup": {
    "help": "Temporary placeholder"
  },
  "build": {
    "help": "Build targets",
    "subcommands": {
      "inno_setup": {
        "help": "Builds the Inno Setup EXE file.",
        "func": build_inno_setup.build_full_setup_exe
      },
      "update_inno_setup": {
        "help": "Builds the update Inno Setup EXE file.",
        "func": build_inno_setup.build_update_setup_exe
      }
    }
  },
  "make": {
    "help": "Creates different resources",
    "subcommands": {
      "docs": {
        "help": "Creates the HTML files for the End-User documentation",
        "func": make_end_user_docs.make_docs
      }
    }
  }
}


if __name__ == '__main__':
  from task_automator import automator
  automator.Automator(AUTOMATION_TREE).run()
