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


class SQLQueryStatement(enum.Enum):
    """An enum for all possible sql query statements."""
    # <editor-fold desc="Util statements">
    GET_LAST_INSERT_ROW_ID = """SELECT last_insert_rowid();"""
    GET_LATEST_ID_OF_PROTEIN_TABLE = """SELECT MAX(id) FROM Protein"""
    GET_LATEST_ID_OF_CHAIN_TABLE = """SELECT MAX(id) FROM Chain"""
    # </editor-fold>

    # <editor-fold desc="Project statements">
    INSERT_NEW_PROJECT = """   
        INSERT INTO Project(name, os)
        VALUES(?, ?)
    """
    GET_PROJECT_ID = """
        SELECT id 
        FROM Project 
        WHERE name = ?
    """
    UPDATE_PROJECT_NAME = """   
        UPDATE Project 
        SET name = ?
        WHERE id = ?
    """
    # </editor-fold>

    # <editor-fold desc="Sequence statements">
    INSERT_SEQUENCE = """
        INSERT INTO SeqRecord(seq_id, seq, name, project_id)
        VALUES (?, ?, ?, ?)
    """
    DELETE_SEQUENCE = """
        DELETE FROM SeqRecord
        WHERE name = ?
    """
    GET_SEQUENCES = """SELECT seq_id, seq, name FROM SeqRecord WHERE project_id = ?"""
    UPDATE_SEQUENCE_NAME = """   
        UPDATE SeqRecord 
        SET name = ?
        WHERE name = ? and seq = ?
    """
    # </editor-fold>

    # <editor-fold desc="Protein statements">
    INSERT_PROTEIN = """   
        INSERT INTO Protein(pymol_molecule_object, pymol_session, project_id)
        VALUES (?, ?, ?)
    """
    INSERT_CHAIN = """   
        INSERT INTO Chain(protein_id, chain_identifier, chain_type, chain_sequence)
        VALUES (?, ?, ?, ?)
    """
    INSERT_PYMOL_PARAMETER = """   
        INSERT INTO PyMOLParameter(color, representation, chain_id)
        VALUES (?, ?, ?)
    """
    INSERT_PYMOL_SELECTION = """   
        INSERT INTO PyMOLSelection(selection_string, protein_id)
        VALUES (?, ?)
    """
    INSERT_PDB_ATOM = """   INSERT INTO PdbAtom(record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
        chain_identifier, residue_sequence_number, code_for_insertions_of_residues,
        x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol,
        charge, protein_id)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    DELETE_PROTEIN = """   
        DELETE FROM Protein
        WHERE id = ?
    """
    DELETE_CHAINS = """   
        DELETE FROM Chain
        WHERE protein_id = ?
    """
    DELETE_PYMOL_PARAMETER = """   
        DELETE FROM PyMOLParameter
        WHERE chain_id = ?
    """
    DELETE_PYMOL_SELECTION = """   
        DELETE FROM PyMOLSelection
        WHERE protein_id = ?
    """
    DELETE_PDB_ATOM = """   
        DELETE FROM PdbAtom
        WHERE protein_id = ?
    """
    GET_PROTEINS = """SELECT id, pymol_molecule_object, pymol_session FROM Protein WHERE project_id = ?"""
    GET_PYMOL_SELECTIONS = """SELECT selection_string FROM PyMOLSelection WHERE protein_id = ?"""
    GET_CHAINS = """SELECT id, chain_identifier, chain_type, chain_sequence FROM Chain WHERE protein_id = ?"""
    GET_PYMOL_PARAMETERS = """SELECT color, representation FROM PyMOLParameter WHERE chain_id = ?"""
    GET_PROTEIN_NAME_BY_ID = """SELECT pymol_molecule_object FROM Protein WHERE id = ?"""
    GET_PDB_ATOMS_OF_PROTEIN = """   
        SELECT record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
        chain_identifier, residue_sequence_number, code_for_insertions_of_residues,
        x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol,
        charge, protein_id 
        FROM PdbAtom 
        WHERE protein_id = ?
    """
    GET_PROTEIN_BY_NAME = """SELECT id, pymol_molecule_object, pymol_session FROM Protein WHERE pymol_molecule_object = ?"""
    UPDATE_PROTEIN_CHAIN_COLOR = """
        UPDATE PyMOLParameter 
        SET color = ?
        WHERE chain_id = ?
    """
    UPDATE_PROTEIN_CHAIN_REPRESENTATION = """
        UPDATE PyMOLParameter 
        SET representation = ?
        WHERE chain_id = ?
    """
    UPDATE_PROTEIN_NAME = """
        UPDATE Protein 
        SET pymol_molecule_object = ?
        WHERE pymol_molecule_object = ? and id = ?
    """
    UPDATE_PYMOL_SESSION_OF_PROTEIN = """   
        UPDATE Protein 
        SET pymol_session = ?
        WHERE id = ?
    """
    # </editor-fold>

    # <editor-fold desc="Protein pair statements">
    INSERT_PROTEIN_PAIR = """   
        INSERT INTO ProteinPair(protein_1_id, protein_2_id, pymol_session, project_id, name)
        VALUES (?, ?, ?, ?, ?)
    """
    INSERT_PYMOL_PARAMETER_PROTEIN_PAIR = """   
        INSERT INTO PyMOLParameterProteinPair(protein_id, chain_letter, parameter_name, parameter_value, protein_pair_id)
        VALUES (?, ?, ?, ?, ?)
    """
    INSERT_DISTANCE_ANALYSIS = """   
        INSERT INTO DistanceAnalysis(name, cutoff, cycles, protein_pair_id, figure_size_x, figure_size_y)
        VALUES (?, ?, ?, ?, ?, ?)
    """
    INSERT_DISTANCE_ANALYSIS_RESULTS = """   
        INSERT INTO DistanceAnalysisResults(pymol_session, rmsd, aligned_aa, distance_analysis_id)
        VALUES (?, ?, ?, ?)
    """
    INSERT_DISTANCE_DATA_RECORDS = """   
        INSERT INTO DistanceAnalysisResultData(my_index, protein_1_chain, protein_1_position, protein_1_residue,
                                               protein_2_chain, protein_2_position, protein_2_residue, 
                                               distances, distance_analysis_results_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    DELETE_PROTEIN_PAIR = """   
        DELETE FROM ProteinPair
        WHERE id = ?
    """
    DELETE_PYMOL_PARAMETERS_PROTEIN_PAIR = """   
        DELETE FROM PyMOLParameterProteinPair
        WHERE protein_pair_id = ?
    """
    DELETE_DISTANCE_ANALYSIS = """   
        DELETE FROM DistanceAnalysis
        WHERE protein_pair_id = ?
    """
    DELETE_DISTANCE_ANALYSIS_RESULTS = """   
        DELETE FROM DistanceAnalysisResults
        WHERE distance_analysis_id = ?
    """
    DELETE_DISTANCE_ANALYSIS_RESULT_DATA = """   
        DELETE FROM DistanceAnalysisResultData
        WHERE distance_analysis_results_id = ?
    """
    GET_PROTEIN_PAIRS = """SELECT id, protein_1_id, protein_2_id, pymol_session, name FROM ProteinPair WHERE project_id = ?"""
    GET_PROTEIN_PAIR_BY_NAME = """SELECT id, protein_1_id, protein_2_id, pymol_session, name FROM ProteinPair WHERE name = ?"""
    GET_DISTANCE_ANALYSIS = """SELECT id, name, cutoff, cycles, figure_size_x, figure_size_y FROM DistanceAnalysis WHERE protein_pair_id = ?"""
    GET_DISTANCE_ANALYSIS_RESULTS = """SELECT id, pymol_session, rmsd, aligned_aa FROM DistanceAnalysisResults WHERE distance_analysis_id = ?"""
    GET_DISTANCE_ANALYSIS_RESULT_DATA = """   
        SELECT my_index, protein_1_chain, protein_1_position, protein_1_residue,
               protein_2_chain, protein_2_position, protein_2_residue, distances
        FROM DistanceAnalysisResultData WHERE distance_analysis_results_id = ?
    """
    GET_PYMOL_PARAMETER_FOR_PROTEIN_CHAIN_IN_PROTEIN_PAIR = """   
        SELECT parameter_value
        FROM PyMOLParameterProteinPair WHERE protein_id = ? and chain_letter = ? and protein_pair_id = ? and parameter_name = ?
    """
    UPDATE_PYMOL_PARAMETER_FOR_PROTEIN_CHAIN_IN_PROTEIN_PAIR = """   
        UPDATE PyMOLParameterProteinPair 
        SET parameter_value = ?
        WHERE protein_id = ? and chain_letter = ? and protein_pair_id = ? and parameter_name = ?
    """
    UPDATE_PYMOL_SESSION_OF_PROTEIN_PAIR = """   
        UPDATE ProteinPair 
        SET pymol_session = ?
        WHERE id = ?
    """

    # </editor-fold>


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

    OPENING_PROJECT = "Opening project ..."
    OPENING_PROJECT_FINISHED = "Opening project finished."
    OPENING_PROJECT_FAILED = "Opening project failed!"


class JobType(enum.Enum):
    """An enum for all job types."""
    PREDICTION = "structure prediction"
    DISTANCE_ANALYSIS = "distance analysis"
    PREDICTION_AND_DISTANCE_ANALYSIS = "prediction and distance analysis"
    RAY_TRACING = "ray-tracing"
    GENERAL_PURPOSE = "general-purpose"


class JobProgress(enum.Enum):
    """An enum for all job progress states."""
    WAITING = "waiting"
    RUNNING = "running"
    FINISHED = "finished"
    FAILED = "failed"


class JobDescriptionKeys(enum.Enum):
    """An enum for all job description keys.

    Notes:
        IMPORTANT: always use the .value of the enum! Otherwise, the connection to the auxiliary pymol would fail!
    """
    JOB_TYPE = "job_type"
    # Prediction keys
    PDB_FILEPATH = "pdb_filepath"
    # Distance analysis keys
    PROTEIN_PAIR_NAME = "the_protein_pair_name"
    PROTEIN_1_PDB_CACHE_FILEPATH = "a_protein_1_pdb_cache_filepath"
    PROTEIN_2_PDB_CACHE_FILEPATH = "a_protein_2_pdb_cache_filepath"
    PROTEIN_1_PYMOL_SELECTION_STRING = "a_protein_1_pymol_selection_string"
    PROTEIN_2_PYMOL_SELECTION_STRING = "a_protein_2_pymol_selection_string"
    CUTOFF = "a_cutoff"
    CYCLES = "the_cycles"
    # Ray-tracing keys
    IMAGE_DESTINATION_FILEPATH = "dest"
    CACHED_SESSION_FILEPATH = "cached"
    RAY_TRACE_MODE = "mode"
    RAY_TEXTURE = "texture"
    RAY_TRACING_RENDERER = "renderer"
    # General purpose
    JOB_SHORT_DESCRIPTION = "job_short_description"
    PYMOL_SESSION = "pymol_session"
    PROTEIN_NAME = "protein_name"


class JobShortDescription(enum.Enum):
    """An enum for all job short descriptions.

    Notes:
        IMPORTANT: always use the .value of the enum! Otherwise, the connection to the auxiliary pymol would fail!
    """
    RUN_STRUCTURE_PREDICTION = "Run ColabFold structure prediction."
    RUN_DISTANCE_ANALYSIS = "Run distance analysis."
    CREATE_RAY_TRACED_IMAGE = "Create new ray traced image."
    CREATE_NEW_PROTEIN_PYMOL_SESSION = "Create new protein pymol session."
    GET_ALL_CHAINS_OF_GIVEN_PROTEIN = "Get all chains of a given protein."
    GET_ALL_SCENES_OF_SESSION = "Get all scenes of a given pymol session."
    CONSOLIDATE_MOLECULE_OBJECT_TO_FIRST_STATE = "Consolidate molecule object to first state."
    CLEAN_PROTEIN_UPDATE_STRUCTURE = "Clean the existing protein structure."
