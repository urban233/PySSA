
from utils import settings_utils
from utils import project_constants

global_var_settings_obj = settings_utils.Settings(project_constants.SETTINGS_DIR, project_constants.SETTINGS_FILE)
global_var_settings_obj.load_settings_from_xml()
