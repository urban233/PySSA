# vim: set expandtab shiftwidth=4 softtabstop=4:
import pathlib
import sys

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

import os

# Add the script's directory to the Python path
script_dir = os.path.dirname(os.path.abspath(__file__))
script_path = pathlib.Path(script_dir)
sys.path.insert(0, str(script_path))
print(sys.path)

from chimerax.core.tools import ToolInstance
from gui import adapter_frame


class AdapterTool(ToolInstance):

    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    help = "help:user/tools/tutorial.html"  # Let ChimeraX know about our help page

    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "PyDD-Adapter"

        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        #
        # Note that by default, tool windows are only hidden rather than
        # destroyed when the user clicks the window's close button.  To change
        # this behavior, specify 'close_destroys=True' in the MainToolWindow
        # constructor.
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

        # Our user interface is simple enough that we could probably inline
        # the code right here, but for any kind of even moderately complex
        # interface, it is probably better to put the code in a method so
        # that this __init__ method remains readable.
        self._build_ui()

    def _build_ui(self):
        # Put our widgets in the tool window
        self.adapter_frame = adapter_frame.AdapterFrame(self, self.tool_window, self.session)
        self.tool_window.ui_area.setLayout(self.adapter_frame.ui.master_layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')

        # from PyQt6.QtWidgets import QMenuBar, QMenu
        # from PyQt6.QtGui import QAction
        # from pyssa.units.core.gui.forms.auto import auto_x_main_view
        # from pyssa.units.core.gui.forms.auto import auto_test_prototype
        # ui = auto_x_main_view.Ui_Dialog()
        # ui = auto_test_prototype.Ui_Dialog()
        # ui.setupUi(self.tool_window.ui_area)

    def delete(self):
        # Perform your cleanup tasks here
        print("Cleaning up resources...")
        # Example: if you need to close files, stop threads, etc.
        # Call the superclass's delete method to ensure base class cleanup
        #super().delete()

    def open_custom_help(self):
        # This method opens the docs of PySSA in their QWebEngineView
        from PyQt6.QtWebEngineWidgets import QWebEngineView
        from PyQt6.QtCore import QUrl
        self.browser = QWebEngineView()
        self.browser.setUrl(QUrl("https://www.youtube.com/watch?v=NIwRJjk-LzI"))
        self.browser.show()

    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.)
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.line_edit.clear())
        menu.addAction(clear_action)

    def take_snapshot(self, session, flags):
        return {
            'version': 1,
            'current text': self.line_edit.text()
        }

    @classmethod
    def restore_snapshot(class_obj, session, data):
        # Instead of using a fixed string when calling the constructor below, we could
        # have saved the tool name during take_snapshot() (from self.tool_name, inherited
        # from ToolInstance) and used that saved tool name.  There are pros and cons to
        # both approaches.
        inst = class_obj(session, "Tutorial (Qt)")
        inst.line_edit.setText(data['current text'])
        return inst
