import datetime
import enum
import os
import pathlib


class ModelDefinitions:
  """Global basic definitions for the core component.

  Should be used as singleton class!
  """
  PROGRAM_ROOT_PATH = pathlib.Path(__file__).parent.parent.parent.parent.parent.parent
  print(PROGRAM_ROOT_PATH)
  PROGRAM_BIN_ROOT_PATH = "C:\\ProgramData\\IBCI\\PyDD\\bin\\PyDD"  # Until now, it's only valid for windows
  """Path to the root of the program."""
  #PROGRAM_SRC_PATH = pathlib.Path(__file__).parent.parent.parent.parent.parent.parent.parent.parent # TODO: Must be switched to if ready for first tests f"{PROGRAM_BIN_ROOT_PATH}\\src"
  PROGRAM_SRC_PATH = pathlib.Path(
    __file__).parent.parent.parent.parent.parent  # TODO: Must be switched to if ready for first tests f"{PROGRAM_BIN_ROOT_PATH}\\src"
  print(PROGRAM_SRC_PATH)
  """Path to the root of the program sources."""
  PLUGIN_LOGO_FILEPATH = str(
    pathlib.Path(f'{PROGRAM_ROOT_PATH}/assets/images/pyssa_logo.png')
  )
  """Filepath to the logo."""
  PLUGIN_LOGO_WITH_CAPTION_FILEPATH = str(
    pathlib.Path(f'{PROGRAM_ROOT_PATH}/assets/images/app_home_logo.png')
  )
  """Filepath to the logo displayed on the homepage."""
  DEFAULT_WORKSPACE_PATH: pathlib.Path = pathlib.Path(
    f"{os.path.expanduser('~')}/.pydd/default_workspace"
  )
  """Default workspace path"""
  SETTINGS_FILENAME = 'settings.json'
  """Default settings filename"""
  DEFAULT_SETTINGS_PATH: pathlib.Path = pathlib.Path(
    f"{os.path.expanduser('~')}/.pydd/"
  )
  DEFAULT_SETTINGS_FILEPATH: pathlib.Path = pathlib.Path(
    f"{os.path.expanduser('~')}/.pydd/{SETTINGS_FILENAME}"
  )
  """Default settings path"""
  DEFAULT_LOG_PATH: pathlib.Path = pathlib.Path(f'{DEFAULT_SETTINGS_PATH}/logs')
  """Default logging path"""
  LOG_FILENAME = f'{datetime.datetime.now().year}-{datetime.datetime.now().month:02d}-{datetime.datetime.now().day:02d}_{datetime.datetime.now().hour:02d}-{datetime.datetime.now().minute:02d}.log'  # noqa: E501
  """Log filename"""
  LOG_FILEPATH = pathlib.Path(f'{DEFAULT_SETTINGS_PATH}/logs/{LOG_FILENAME}')
  """Complete log filepath"""
  ASSETS_PATH = pathlib.Path(f"{PROGRAM_ROOT_PATH}/assets")
  """Path to the assets folder"""
  ICONS_PATH = pathlib.Path(f"{ASSETS_PATH}/icons")
  """Path to the icons folder"""
  AMINO_ACID_CODE = {
    'CYS': 'C',
    'ASP': 'D',
    'SER': 'S',
    'GLN': 'Q',
    'LYS': 'K',
    'ILE': 'I',
    'PRO': 'P',
    'THR': 'T',
    'PHE': 'F',
    'ASN': 'N',
    'GLY': 'G',
    'HIS': 'H',
    'LEU': 'L',
    'ARG': 'R',
    'TRP': 'W',
    'ALA': 'A',
    'VAL': 'V',
    'GLU': 'E',
    'TYR': 'Y',
    'MET': 'M',
  }
  """The three and one letter code of all human relevant amino acids."""
  STATUS_MESSAGE_TIMEOUT = 5000  # value in msec
  """The standard status bar message timeout"""

# <editor-fold desc="Model roles">
class RolesEnum(enum.IntEnum):
  """Enumeration for storing possible model roles."""
  OBJECT_ROLE = 1003
  """Model role to store an object."""
  TYPE_ROLE = 1004
  """Model role to store a type."""
  FILEPATH_ROLE = 1005
  """Model role to store a filepath."""
  CHAIN_COLOR_ROLE = 1006
  """Model role to store a chain color."""
# </editor-fold>


# <editor-fold desc="Model types">
class TypesEnum(enum.StrEnum):
  """Enumeration for storing possible model types."""
  PROJECT_TYPE = "project"
  """Model type for storing a project."""
  SEQUENCE_TYPE = "sequence"
  """Model type for storing a sequence."""
  PROTEIN_TYPE = "protein"
  """Model type for storing a protein."""
  CHAIN_TYPE = "chain"
  """Model type for storing a chain."""
  RESIDUE_TYPE = "residue"
  """Model type for storing a residue."""
  ATOM_TYPE = "atom"
  """Model type for storing an atom."""
  LIGAND_TYPE = "ligand"
  """Model type for storing an ligand."""
# </editor-fold>


class ComponentsEnum(enum.IntEnum):
  """Enumeration for storing all components."""
  CREATE_PROJECT = 1001
  OPEN_PROJECT = 1002

  IMPORT_SEQUENCE = 1101
  IMPORT_PROTEIN = 1201
  SAVE_PROTEIN = 1202
  IMPORT_LIGAND = 1301


class IconsEnum(enum.StrEnum):
  """Class for storing the icon filenames with their usage name."""
  HELP = "HELP"
  ADD_SEQUENCE = "NOTE_ADD"
  ADD_SEQUENCE_DISABLED = "NOTE_ADD_DISABLED"
  IMPORT_SEQUENCE = "UPLOAD_FILE"
  IMPORT_SEQUENCE_DISABLED = "UPLOAD_FILE_DISABLED"
  SAVE_SEQUENCE = "FILE_SAVE"
  SAVE_SEQUENCE_DISABLED = "FILE_SAVE_DISABLED"
  DELETE_SEQUENCE = "SCAN_DELETE"
  DELETE_SEQUENCE_DISABLED = "SCAN_DELETE_DISABLED"
  EXPAND_ALL = "EXPAND_ALL"
  COLLAPSE_ALL = "COLLAPSE_ALL"
  IMPORT_PROTEIN = "UPLOAD_FILE"
  IMPORT_PROTEIN_DISABLED = "UPLOAD_FILE_DISABLED"
  SAVE_PROTEIN = "FILE_SAVE"
  SAVE_PROTEIN_DISABLED = "FILE_SAVE_DISABLED"
  DELETE_PROTEIN = "SCAN_DELETE"
  DELETE_PROTEIN_DISABLED = "SCAN_DELETE_DISABLED"
  OPEN_SESSION = "OPEN_IN_NEW"
  OPEN_SESSION_DISABLED = "OPEN_IN_NEW_DISABLED"
  CREATE_SESSION_SCENE = "ADD_CIRCLE"
  CREATE_SESSION_SCENE_DISABLED = "ADD_CIRCLE_DISABLED"
  UPDATE_SESSION_SCENE = "CHANGE_CIRCLE"
  UPDATE_SESSION_SCENE_DISABLED = "CHANGE_CIRCLE_DISABLED"
  DELETE_SESSION_SCENE = "CANCEL"
  DELETE_SESSION_SCENE_DISABLED = "CANCEL_DISABLED"
  DELETE_PROTEIN_PAIR = "SCAN_DELETE"
  DELETE_PROTEIN_PAIR_DISABLED = "SCAN_DELETE_DISABLED"
  JOBS = "PLAY_CIRCLE"
  JOBS_RUNNING = "PLAY_CIRCLE_RUN"
  NOTIFY = "NOTIFICATIONS"
  NOTIFY_UNREAD = "NOTIFICATIONS_UNREAD"
  NEW = "DRAFT"
  OPEN = "FOLDER_OPEN"
  EXPORT = "SHARE_WINDOWS"
  IMPORT = "BROWSER_UPDATED"
  CLOSE = "CLOSE"
  DELETE = "DELETE"
  BACK = "ARROW_BACK"
  DISTANCE_ANALYSIS = "ANGSTROM_DISTANCE"
  CARTOON_REPR = "CARTOON_REPR"
  STICKS_REPR = "STICKS_REPR"
  RIBBON_REPR = "RIBBON_REPR"
  SPHERES_REPR = "SPHERES_REPR"
  LINES_REPR = "LINES_REPR"
  DOTS_REPR = "DOTS_REPR"
  MESH_REPR = "MESH_REPR"
  SURFACE_REPR = "SURFACE_REPR"
  MORE = "MORE_VERT"
  COLOR_GRID = "GRID_VIEW"
  IMAGE = "IMAGE"
  RAY_TRACED_IMAGE = "PHOTO_FRAME"
  EDIT = "EDIT"
  MONOMER_PROTEIN = "MONOMER_PROTEIN"
  MULTIMER_PROTEIN = "MULTIMER_PROTEIN"
  DOCKING = "DOCKING"
  RIBBON_PANEL_SETTINGS = "OPEN_IN_NEW_DOWN"
  HOME = "HOME"
