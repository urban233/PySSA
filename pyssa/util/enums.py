import enum


class ApplicationModelEnum(enum.Enum):
    """An enum for all application model keys."""
    PROJECT = "current_project"


class ModelEnum(enum.IntEnum):
    """An enum for all model keys."""
    OBJECT_ROLE = 1003
    TYPE_ROLE = 1004
    FILEPATH_ROLE = 1005
