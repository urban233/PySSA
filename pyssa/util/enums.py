import enum


class ApplicationModelEnum(enum.Enum):
    """An enum for all application model keys."""
    PROJECT = "current_project"


class ModelEnum(enum.IntEnum):
    """An enum for all model keys."""
    OBJECT_ROLE = 1003
    TYPE_ROLE = 1004
    FILEPATH_ROLE = 1005


class ModelTypeEnum(enum.Enum):
    """An enum for all model types."""
    MONOMER_SEQ = "monomer_seq"
    MULTIMER_SEQ = "multimer_seq"


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


class PdbAtomEnum(enum.Enum):
    """An enum for the PdbAtom table of the database"""
    RECORD_TYPE = "record_type"
    ATOM_NUMBER = "atom_number"
    ATOM_NAME = "atom_name"
    ALTERNATE_LOCATION_INDICATOR = "alternate_location_indicator"
    RESIDUE_NAME = "residue_name"
    CHAIN_IDENTIFIER = "chain_identifier"
    RESIDUE_SEQUENCE_NUMBER = "residue_sequence_number"
    CODE_FOR_INSERTIONS_OF_RESIDUES = "code_for_insertions_of_residues"
    X_COORD = "x_coord"
    Y_COORD = "y_coord"
    Z_COORD = "z_coord"
    OCCUPANCY = "occupancy"
    TEMPERATURE_FACTOR = "temperature_factor"
    SEGMENT_IDENTIFIER = "segment_identifier"
    ELEMENT_SYMBOL = "element_symbol"
    CHARGE = "charge"


class PymolParameterEnum(enum.Enum):
    """An enum for all pymol parameters."""
    COLOR = "chain_color"
    REPRESENTATION = "chain_representation"


class SQLQueryType(enum.Enum):
    """An enum for all possible sql queries for the database thread."""
    CLOSE_PROJECT = "close_project"
    INSERT_NEW_PROTEIN = 'insert_new_protein'
    DELETE_EXISTING_PROTEIN = 'delete_existing_protein'
    INSERT_NEW_PROTEIN_PAIR = 'insert_new_protein_pair'
    DELETE_EXISTING_PROTEIN_PAIR = 'delete_existing_protein_pair'
    UPDATE_PYMOL_SESSION_PROTEIN = 'update_pymol_session_protein'
    UPDATE_PYMOL_SESSION_PROTEIN_PAIR = 'update_pymol_session_protein_pair'
    INSERT_NEW_SEQUENCE = 'insert_new_sequence'
    DELETE_EXISTING_SEQUENCE = 'delete_existing_sequence'
    UPDATE_SEQUENCE_NAME = 'update_sequence_name'


class HistogramPropertiesEnum(enum.Enum):
    """An enum for all possible histogram properties."""
    X_AXIS_UNITS = 'x_axis_units'
    DISTANCE_INTERVAL = 'distance_interval'


class DocsServerStatus(enum.IntEnum):
    """An enum for all possible docs server status."""
    INACTIVE = 0
    ACTIVE = 1
    PENDING = 2


class PyMOLRepresentation(enum.Enum):
    """An enum for possible representations."""
    CARTOON = 'cartoon'
    STICKS = 'sticks'
    RIBBON = 'ribbon'
    LINES = 'lines'
    SPHERES = 'spheres'
    DOTS = 'dots'
    MESH = 'mesh'
    SURFACE = 'surface'


class StatusMessages(enum.Enum):
    """An enum for all status messages, displayed in the statusbar."""
    PREDICTION_IS_RUNNING = "A structure prediction is running ..."
    PREDICTION_IS_FINALIZING = "The structure prediction is finalizing ..."
    DISTANCE_ANALYSIS_IS_RUNNING = "A distance analysis is running ..."

    OPENING_PROJECT = "Opening existing project ..."
    OPENING_PROJECT_FINISHED = "Opening existing project finished."
    OPENING_PROJECT_FAILED = "Opening existing project failed!"


class JobType(enum.Enum):
    """An enum for all job types."""
    PREDICTION = "structure prediction"
    DISTANCE_ANALYSIS = "distance analysis"
    PREDICTION_AND_DISTANCE_ANALYSIS = "prediction and distance analysis"
    RAY_TRACING = "ray-tracing"


class JobProgress(enum.Enum):
    """An enum for all job progress states."""
    WAITING = "waiting"
    RUNNING = "running"
    FINISHED = "finished"
    FAILED = "failed"
