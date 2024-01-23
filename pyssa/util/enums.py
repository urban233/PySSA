import enum


class ApplicationModelEnum(enum.Enum):
    """An enum for all application model keys."""
    PROJECT = "current_project"


class ModelEnum(enum.IntEnum):
    """An enum for all model keys."""
    OBJECT_ROLE = 1003
    TYPE_ROLE = 1004
    FILEPATH_ROLE = 1005


class DatabaseEnum(enum.Enum):
    """An enum for all database fields."""
    PROJECT_NAME = "name"
    PROJECT_OS = "os"
    PROTEIN_NAME = "pymol_molecule_object"
    PROTEIN_PDB_FILEPATH = "pdb_filepath"
    PROTEIN_FASTA_FILEPATH = "fasta_filepath"
    PROTEIN_EXPORT_DIRNAME = "export_dirname"
    PROTEIN_PYMOL_SESSION_FILEPATH = "pymol_session_filepath"
    PROTEIN_PYMOL_SESSION = "pymol_session"
    # add more for the rest of the tables


class PymolParameterEnum(enum.Enum):
    """An enum for all pymol parameters."""
    COLOR = "chain_color"
    REPRESENTATION = "chain_representation"


class SQLQueryType(enum.Enum):
    """An enum for all possible sql queries for the database thread."""
    INSERT_NEW_PROTEIN = 'insert_new_protein'
    DELETE_EXISTING_PROTEIN = 'delete_existing_protein'

