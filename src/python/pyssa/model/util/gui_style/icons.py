import logging
import pathlib
from PyQt6 import QtGui, QtCore, QtWidgets
from pyssa.model.preference import model_definitions
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


# Define the root path for icons
_icons_root_path = model_definitions.ModelDefinitions.ICONS_PATH

# Create a dictionary to store icon paths
ICON_PATHS = {
  "HELP": pathlib.Path(f"{_icons_root_path}/help_w200.svg"),
  "OPEN_IN_NEW": pathlib.Path(f"{_icons_root_path}/open_in_new_w200.svg"),
  "OPEN_IN_NEW_DISABLED": pathlib.Path(f"{_icons_root_path}/open_in_new_disabled_w200.svg"),
  "ADD": pathlib.Path(f"{_icons_root_path}/add_w200.svg"),
  "ADD_DISABLED": pathlib.Path(f"{_icons_root_path}/add_disabled_w200.svg"),
  "ADD_CIRCLE": pathlib.Path(f"{_icons_root_path}/add_circle_w200.svg"),
  "ADD_CIRCLE_DISABLED": pathlib.Path(f"{_icons_root_path}/add_circle_disabled_w200.svg"),
  "NOTE_ADD": pathlib.Path(f"{_icons_root_path}/note_add_w200.svg"),
  "NOTE_ADD_DISABLED": pathlib.Path(f"{_icons_root_path}/note_add_disabled_w200.svg"),
  "CHANGE_CIRCLE": pathlib.Path(f"{_icons_root_path}/change_circle_w200.svg"),
  "CHANGE_CIRCLE_DISABLED": pathlib.Path(f"{_icons_root_path}/change_circle_disabled_w200.svg"),
  "UPLOAD_FILE": pathlib.Path(f"{_icons_root_path}/upload_file_w200.svg"),
  "UPLOAD_FILE_DISABLED": pathlib.Path(f"{_icons_root_path}/upload_file_disabled_w200.svg"),
  "FILE_SAVE": pathlib.Path(f"{_icons_root_path}/file_save_w200.svg"),
  "FILE_SAVE_DISABLED": pathlib.Path(f"{_icons_root_path}/file_save_disabled_w200.svg"),
  "SCAN_DELETE": pathlib.Path(f"{_icons_root_path}/scan_delete_w200.svg"),
  "SCAN_DELETE_DISABLED": pathlib.Path(f"{_icons_root_path}/scan_delete_disabled_w200.svg"),
  "CANCEL": pathlib.Path(f"{_icons_root_path}/cancel_w200.svg"),
  "CANCEL_DISABLED": pathlib.Path(f"{_icons_root_path}/cancel_disabled_w200.svg"),
  "DANGEROUS": pathlib.Path(f"{_icons_root_path}/dangerous_w200.svg"),
  "ERROR": pathlib.Path(f"{_icons_root_path}/error_w200.svg"),
  "INFO": pathlib.Path(f"{_icons_root_path}/info_w200.svg"),
  "WARNING": pathlib.Path(f"{_icons_root_path}/warning_w200.svg"),
  "KEYBOARD_ARROW_RIGHT": pathlib.Path(f"{_icons_root_path}/keyboard_arrow_right_w400.svg"),
  "KEYBOARD_ARROW_DOWN": pathlib.Path(f"{_icons_root_path}/keyboard_arrow_down_w400.svg"),
  "ARROW_DROP_UP": pathlib.Path(f"{_icons_root_path}/arrow_drop_up_w400.svg"),
  "ARROW_RIGHT": pathlib.Path(f"{_icons_root_path}/arrow_right_w400.svg"),
  "DONE_ROUND_EDGES": pathlib.Path(f"{_icons_root_path}/done_round_edges_w200_g200.svg"),
  "DO_NOT_DISTURB_ON": pathlib.Path(f"{_icons_root_path}/do_not_disturb_on_w200.svg"),
  "NOTIFICATIONS": pathlib.Path(f"{_icons_root_path}/notifications_w200.svg"),
  "NOTIFICATIONS_UNREAD": pathlib.Path(f"{_icons_root_path}/notifications_unread_w200.svg"),
  "PLAY_CIRCLE": pathlib.Path(f"{_icons_root_path}/play_circle_w200.svg"),
  "PLAY_CIRCLE_RUN": pathlib.Path(f"{_icons_root_path}/play_circle_run_w200_blue.svg"),
  "EXPAND_ALL": pathlib.Path(f"{_icons_root_path}/expand_all_w200.svg"),
  "COLLAPSE_ALL": pathlib.Path(f"{_icons_root_path}/collapse_all_w200.svg"),
  "DRAFT": pathlib.Path(f"{_icons_root_path}/draft_w200.svg"),
  "BROWSER_UPDATED": pathlib.Path(f"{_icons_root_path}/browser_updated_w200.svg"),
  "FOLDER_OPEN": pathlib.Path(f"{_icons_root_path}/folder_open_w200.svg"),
  "SHARE_WINDOWS": pathlib.Path(f"{_icons_root_path}/share_windows_w200.svg"),
  "CLOSE": pathlib.Path(f"{_icons_root_path}/close_w200.svg"),
  "DELETE": pathlib.Path(f"{_icons_root_path}/delete_w200.svg"),
  "ARROW_BACK": pathlib.Path(f"{_icons_root_path}/arrow_back_w200.svg"),
  "ANGSTROM_DISTANCE": pathlib.Path(f"{_icons_root_path}/angstrom_distance_w200.svg"),
  "CARTOON_REPR": pathlib.Path(f"{_icons_root_path}/cartoon_repr.png"),
  "STICKS_REPR": pathlib.Path(f"{_icons_root_path}/sticks_repr.png"),
  "RIBBON_REPR": pathlib.Path(f"{_icons_root_path}/ribbon_repr.png"),
  "SPHERES_REPR": pathlib.Path(f"{_icons_root_path}/spheres_repr.png"),
  "MESH_REPR": pathlib.Path(f"{_icons_root_path}/mesh_repr.png"),
  "DOTS_REPR": pathlib.Path(f"{_icons_root_path}/dots_repr.png"),
  "SURFACE_REPR": pathlib.Path(f"{_icons_root_path}/surface_repr.png"),
  "MORE_VERT": pathlib.Path(f"{_icons_root_path}/more_vert_w200.svg"),
  "GRID_VIEW": pathlib.Path(f"{_icons_root_path}/grid_view_w200.svg"),
  "PHOTO_FRAME": pathlib.Path(f"{_icons_root_path}/photo_frame_w200.svg"),
  "IMAGE": pathlib.Path(f"{_icons_root_path}/image_w200.svg"),
  "EDIT": pathlib.Path(f"{_icons_root_path}/edit_w200.svg"),
  "MONOMER_PROTEIN": pathlib.Path(f"{_icons_root_path}/monomer.png"),
  "MULTIMER_PROTEIN": pathlib.Path(f"{_icons_root_path}/multimer.png"),
  "DOCKING": pathlib.Path(f"{_icons_root_path}/docking.png"),
  "OPEN_IN_NEW_DOWN": pathlib.Path(f"{_icons_root_path}/open_in_new_down_w200.svg"),
  "HOME": pathlib.Path(f"{_icons_root_path}/home_w200.svg"),
}


def get_icon(
        icon_name: model_definitions.IconsEnum,
        disabled_icon_name: model_definitions.IconsEnum = None
):
  """
  Function to get an icon given its name and set its size.

  Args:
    icon_name (IconsEnum): The name of the icon.
    disabled_icon_name (IconsEnum): The name of the disabled icon, if any. (Default: None)

  Notes:
    If you don't want to specify the size of the icon, use None as the size of the icon.

  Returns:
    QtGui.QIcon: The corresponding QIcon object.

  Raises:
    exception.NoneValueError: If `icon_name` is None.
  """
  # <editor-fold desc="Checks">
  if icon_name is None:
    default_logging.append_to_log_file(logger, "icon_name is None.", logging.ERROR)
    raise exception.NoneValueError("icon_name is None.")
  # </editor-fold>
  icon = QtGui.QIcon(QtGui.QPixmap(str(ICON_PATHS[icon_name.value])))
  if disabled_icon_name:
    icon.addPixmap(
      QtGui.QPixmap(str(ICON_PATHS[disabled_icon_name.value])),
      mode=QtGui.QIcon.Mode.Disabled,
    )
  return icon


def set_icon(
      a_button: QtWidgets.QPushButton,
      icon_name: model_definitions.IconsEnum,
      disabled_icon_name: model_definitions.IconsEnum = None,
      size: QtCore.QSize = QtCore.QSize(32, 32)
):
  """
  Function to set an icon given its widget, its name and set its size.

  Args:
    a_button (QWidget): The widget to set the icon for
    icon_name (IconsEnum): The name of the icon.
    disabled_icon_name (IconsEnum): The name of the disabled icon, if any.
    size (QtCore.QSize): The desired size for the icon. Defaults to QtCore.QSize(30, 30)

  Notes:
    If you don't want to specify the size of the icon, use None as the size of the icon.

  Returns:
    QtGui.QIcon: The corresponding QIcon object.

  Raises:
    exception.NoneValueError: If `a_button` ,`icon_name` or `size` is None.
  """
  # <editor-fold desc="Checks">
  if a_button is None:
    default_logging.append_to_log_file(logger, "a_button is None.", logging.ERROR)
    raise exception.NoneValueError("a_button is None.")
  if icon_name is None:
    default_logging.append_to_log_file(logger, "icon_name is None.", logging.ERROR)
    raise exception.NoneValueError("icon_name is None.")
  if size is None:
    default_logging.append_to_log_file(logger, "size is None.", logging.ERROR)
    raise exception.NoneValueError("size is None.")
  # </editor-fold>
  a_button.setIcon(get_icon(icon_name, disabled_icon_name))
  a_button.setIconSize(
    a_button.icon().actualSize(size)
  )
  a_button.setText("")
