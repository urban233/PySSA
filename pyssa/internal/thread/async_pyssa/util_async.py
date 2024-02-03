import os
import subprocess

import zmq
import pygetwindow

from pyssa.util import constants


def open_documentation_on_certain_page(
    a_page_name: str,
    placeholder: int
) -> tuple:
    """Opens the pyssa documentation on the given page name.

    Args:
        a_page_name: a name of a page from the documentation.

    Returns:
        a tuple with ("result", response)
    """
    window_title_to_check = "PySSA - Documentation Center"
    if len(pygetwindow.getWindowsWithTitle(window_title_to_check)):
        print("Update currently open documentation.")
    else:
        try:
            os.chdir(constants.PLUGIN_DOCS_PATH)
            subprocess.Popen(['C:\ProgramData\pyssa\mambaforge_pyssa\pyssa-mamba-env\Scripts\mkdocs.exe', 'serve'])
            subprocess.Popen([f"{constants.PLUGIN_EXTRA_TOOLS_PATH}\\browser.exe"])
        except subprocess.CalledProcessError as e:
            print(f"Error starting mkdocs serve: {e}")

    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.connect("tcp://127.0.0.1:7017")
    socket.send_string(f"update_url http://127.0.0.1:8000/{a_page_name}")
    response = socket.recv_string()
    print(f"Response from server: {response}")
    socket.close()

    return "result", response
