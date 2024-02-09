import logging
import os
import subprocess

import zmq
import pygetwindow

from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


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
        logger.info("Update currently open documentation.")
    else:
        logger.info("Trying to run the mkdocs serve command ...")
        try:
            os.chdir(constants.PLUGIN_DOCS_PATH)
            subprocess.Popen(['C:\ProgramData\pyssa\mambaforge_pyssa\pyssa-mamba-env\Scripts\mkdocs.exe', 'serve', '-a', '127.0.0.1:7018'])
            subprocess.Popen([r"C:\ProgramData\pyssa\extra_tools\browser.exe"])
        except subprocess.CalledProcessError as e:
            logger.error(f"Error starting mkdocs serve: {e}")
        else:
            logger.info("Running the mkdocs serve command finished without errors.")

    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.connect("tcp://127.0.0.1:7017")
    socket.send_string(f"update_url http://127.0.0.1:7018/{a_page_name}")
    response = socket.recv_string()
    logger.debug(f"Response from server: {response}")
    socket.close()

    return "result", response
