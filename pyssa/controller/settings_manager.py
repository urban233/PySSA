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
