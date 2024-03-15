import logging
import os
import subprocess

import zmq
import pygetwindow

from pyssa.controller import interface_manager, pymol_session_manager
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, enums
from pyssa.util import globals

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def open_documentation_on_certain_page(
    a_page_name: str,
    the_docs_window
) -> tuple:
    """Opens the pyssa documentation on the given page name.

    Args:
        a_page_name: a name of a page from the documentation.

    Returns:
        a tuple with ("result", response)
    """
    if the_docs_window is not None:
        tmp_is_docs_open_flag = True
    else:
        tmp_is_docs_open_flag = False

    tmp_docs_window = None
    if tmp_is_docs_open_flag:
        # Docs window is already open
        logger.info("Update currently open documentation.")
        tmp_docs_window = the_docs_window
    elif globals.g_server_status == enums.DocsServerStatus.PENDING:
        # Docs are getting build
        logger.warning("Docs need to be build. The build process was started from another process!")
        flag = False
        while flag is not True:
            if len(pygetwindow.getWindowsWithTitle(
                    constants.WINDOW_TITLE_OF_HELP_CENTER)) == 1 or globals.g_server_status == enums.DocsServerStatus.ACTIVE:
                flag = True
                globals.g_server_status = enums.DocsServerStatus.ACTIVE
                tmp_docs_window = pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[0]
                tmp_docs_window.minimize()
    else:
        # Docs are need to be built in this process
        logger.info("Trying to run the mkdocs serve command ...")
        try:
            os.chdir(constants.PLUGIN_DOCS_PATH)
            subprocess.Popen(['C:\ProgramData\pyssa\mambaforge_pyssa\pyssa-mamba-env\Scripts\mkdocs.exe', 'serve', '-a', '127.0.0.1:7018'])
            subprocess.Popen([r"C:\ProgramData\pyssa\extra_tools\browser.exe"])
        except subprocess.CalledProcessError as e:
            logger.error(f"Error starting mkdocs serve: {e}")
            tmp_docs_window = None
        else:
            logger.info("Running the mkdocs serve command finished without errors.")
            globals.g_server_status = enums.DocsServerStatus.PENDING
            flag = False
            while flag is not True:
                if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) == 1:
                    flag = True
                    globals.g_server_status = enums.DocsServerStatus.ACTIVE
                    tmp_docs_window = pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[0]
                    tmp_docs_window.minimize()

    # After the docs window is ready, it can be accessed by ZeroMQ
    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.connect("tcp://127.0.0.1:7017")
    socket.send_string(f"update_url http://127.0.0.1:7018/{a_page_name}")
    response = socket.recv_string()
    logger.debug(f"Response from server: {response}")
    socket.close()
    return "result", response, tmp_docs_window


def start_documentation_server(
    placeholder_1: int,
    placeholder_2: int
) -> tuple:
    """Opens the pyssa documentation on the given page name.

    Args:
        placeholder_1: just a placeholder for the first argument
        placeholder_2: just a placeholder for the second argument

    Returns:
        a tuple with ("result", response)
    """
    logger.info("Trying to run the mkdocs serve command ...")
    try:
        os.chdir(constants.PLUGIN_DOCS_PATH)
        subprocess.Popen(['C:\ProgramData\pyssa\mambaforge_pyssa\pyssa-mamba-env\Scripts\mkdocs.exe', 'serve', '-a', '127.0.0.1:7018'])
        subprocess.Popen([r"C:\ProgramData\pyssa\extra_tools\browser.exe"])
    except subprocess.CalledProcessError as e:
        logger.error(f"Error starting mkdocs serve: {e}")
    else:
        logger.info("Running the mkdocs serve command finished without errors.")
        globals.g_server_status = enums.DocsServerStatus.PENDING

    flag = False
    while flag is not True:
        if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) == 1:
            flag = True
            globals.g_server_status = enums.DocsServerStatus.ACTIVE
    tmp_docs_window = pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[0]
    tmp_docs_window.minimize()

    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.connect("tcp://127.0.0.1:7017")
    socket.send_string("update_url http://127.0.0.1:7018/help/")
    response = socket.recv_string()
    logger.debug(f"Response from server: {response}")
    socket.close()

    return "result", tmp_docs_window


def unfreeze_pymol_session(the_pymol_session_manager: "pymol_session_manager.PymolSessionManager", placeholder_1: int):
    """Refreshes the protein model of the interface manager."""
    if the_pymol_session_manager.frozen_protein_object is not None:
        logger.info("Unfreezing protein pymol session ...")
        the_pymol_session_manager.unfreeze_current_protein_pymol_session()
        logger.info("Unfreezing protein pymol session finished.")
    elif the_pymol_session_manager.frozen_protein_pair_object is not None:
        logger.info("Unfreezing protein pair pymol session ...")
        the_pymol_session_manager.unfreeze_current_protein_pair_pymol_session()
        logger.info("Unfreezing protein pair pymol session finished.")
    else:
        logger.info("There is no pymol session to unfreeze.")
    return "result", the_pymol_session_manager
