import logging
from pyssa.internal.data_structures import settings
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class SettingsManager:

    settings: "settings.Settings"

    def __init__(self) -> None:
        self.settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        try:
            self.settings = self.settings.deserialize_settings()
        except ValueError:
            logging.error("Settings could not be loaded because of an unexpected error!")
        except TypeError:
            logging.warning("Settings are either damaged or outdated!")
