import sys
from PyQt6 import QtWidgets, QtCore

from pyssa.gui.main import main_frame
#from pyssa.gui.main import main_frame_controller

from pyssa.model.util.gui_style import styles_utils
#from pyssa.extension.mv_wrappers.chimerax_wrapper import chimerax_pydd_client


if __name__ == "__main__":
  app = QtWidgets.QApplication(sys.argv)
  tmp_main_frame = main_frame.MainFrame()
  #tmp_main_frame_controller = main_frame_controller.MainFrameController(tmp_main_frame)
  # tmp_main_frame_controller.mv_client = chimerax_pydd_client.ChimeraxPyddClient()
  # tmp_main_frame_controller.mv_client.start_viewer()
  # tmp_main_frame_controller.mv_client.connect_to_viewer()
  styles_utils.set_stylesheet(app)
  #tmp_main_frame_controller.main_frame.show()
  tmp_main_frame.show()
  sys.exit(app.exec())
