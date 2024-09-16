

class CommonQueryStatements:
  """Class for storing common sql queries for the database manager."""
  # <editor-fold desc="CREATE TABLE statements">
  CREATE_TABLE_PROJECT: str = """-- Project definition
      CREATE TABLE Project (
          id INTEGER NOT NULL,
          name TEXT NOT NULL,
          CONSTRAINT Project_PK PRIMARY KEY (id)
      );
  """
  CREATE_TABLE_SEQUENCE: str = """-- Sequence definition
      CREATE TABLE Sequence (
          id INTEGER NOT NULL,
          name TEXT,
          chain_sequence_map BLOB,
          type TEXT,
          project_id INTEGER,
          CONSTRAINT Sequence_PK PRIMARY KEY (id),
          CONSTRAINT Sequence_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
      );
  """
  CREATE_TABLE_PROTEIN: str = """-- Protein definition
      CREATE TABLE Protein (
          id INTEGER NOT NULL,
          name TEXT,
          structure BLOB,
          chain_sequence_map BLOB,
          session BLOB,
          project_id INTEGER,
          CONSTRAINT Protein_PK PRIMARY KEY (id),
          CONSTRAINT Protein_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
      );
  """
  CREATE_TABLE_LIGAND: str = """-- Ligand definition
      CREATE TABLE Ligand (
          id INTEGER NOT NULL,
          name TEXT,
          structure BLOB,
          session BLOB,
          project_id INTEGER,
          CONSTRAINT Ligand_PK PRIMARY KEY (id),
          CONSTRAINT Ligand_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
      );
  """
  CREATE_TABLE_PROTEIN_LIGAND_COMPLEX: str = """-- Protein-ligand complex definition
      CREATE TABLE ProteinLigandComplex (
          id INTEGER NOT NULL,
          name TEXT,
          protein_id INTEGER,
          ligand_id INTEGER,
          results BLOB,
          session BLOB,
          project_id INTEGER,
          CONSTRAINT ProteinLigandComplex_PK PRIMARY KEY (id),
          CONSTRAINT ProteinLigandComplex_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id)
          CONSTRAINT ProteinLigandComplex_Ligand_FK FOREIGN KEY (ligand_id) REFERENCES Ligand(id)
          CONSTRAINT ProteinLigandComplex_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
      );
  """
  CREATE_TABLE_PROTEIN_PROTEIN_COMPLEX: str = """-- Protein-protein complex definition
      CREATE TABLE ProteinProteinComplex (
          id INTEGER NOT NULL,
          name TEXT,
          protein_1_id INTEGER,
          protein_2_id INTEGER,
          results BLOB,
          session BLOB,
          project_id INTEGER,
          CONSTRAINT ProteinProteinComplex_PK PRIMARY KEY (id),
          CONSTRAINT ProteinProteinComplex_Protein_1_FK FOREIGN KEY (protein_1_id) REFERENCES Protein(id)
          CONSTRAINT ProteinProteinComplex_Protein_2_FK FOREIGN KEY (protein_2_id) REFERENCES Protein(id)
          CONSTRAINT ProteinProteinComplex_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
      );
  """
  CREATE_TABLE_DUMMY_PROTEIN: str = """-- DummyProtein definition
      CREATE TABLE DummyProtein (
          id INTEGER NOT NULL,
          name TEXT,
          chain_sequence_map BLOB,
          project_id INTEGER,
          CONSTRAINT Protein_PK PRIMARY KEY (id),
          CONSTRAINT Protein_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
      );
  """
  CREATE_TABLE_DUMMY_PROTEIN_LIGAND_COMPLEX: str = """-- Dummy protein-ligand complex definition
      CREATE TABLE DummyProteinLigandComplex (
          id INTEGER NOT NULL,
          name TEXT,
          project_id INTEGER,
          CONSTRAINT ProteinLigandComplex_PK PRIMARY KEY (id),
          CONSTRAINT ProteinLigandComplex_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
      );
  """
  CREATE_TABLE_DUMMY_PROTEIN_PROTEIN_COMPLEX: str = """-- Dummy protein-protein complex definition
      CREATE TABLE DummyProteinProteinComplex (
          id INTEGER NOT NULL,
          name TEXT,
          project_id INTEGER,
          CONSTRAINT ProteinProteinComplex_PK PRIMARY KEY (id),
          CONSTRAINT ProteinProteinComplex_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
      );
  """
  CREATE_TABLE_JOB: str = """-- Job definition
        CREATE TABLE Job (
            id INTEGER NOT NULL,
            name TEXT,
            type TEXT,
            status INTEGER,
            dummy_id INTEGER,
            dummy_type TEXT,
            project_id INTEGER,
            CONSTRAINT Jobs_PK PRIMARY KEY (id),
            CONSTRAINT Jobs_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
        );
    """
  # </editor-fold>

  # <editor-fold desc="Project statements">
  INSERT_PROJECT = """   
        INSERT INTO Project(name)
        VALUES (?)
  """
  # </editor-fold>

  # <editor-fold desc="Sequence statements">
  INSERT_SEQUENCE = """   
        INSERT INTO Sequence(name, seq_record, type, project_id)
        VALUES (?, ?, ?, ?)
  """
  UPDATE_CHAIN_SEQUENCE_MAP = """
        UPDATE Sequence 
        SET chain_sequence_map = ?
        WHERE name = ?
  """
  DELETE_SEQUENCE = """   
        DELETE FROM Sequence
        WHERE name = ?
  """
  GET_SEQUENCE = """SELECT id, name, chain_sequence_map, type FROM Sequence WHERE name = ?"""
  GET_ALL_SEQUENCES = """SELECT id, name, chain_sequence_map, type FROM Sequence WHERE project_id = ?"""
  # </editor-fold>

  # <editor-fold desc="Protein statements">
  INSERT_PROTEIN = """   
      INSERT INTO Protein(name, structure, chain_sequence_map, session, project_id)
      VALUES (?, ?, ?, ?, ?)
  """
  UPDATE_PROTEIN_STRUCTURE = """
      UPDATE Protein 
      SET structure = ?
      WHERE name = ?
  """
  UPDATE_PROTEIN_CHAIN_SEQ_MAP = """
      UPDATE Protein 
      SET chain_sequence_map = ?
      WHERE name = ?
  """
  UPDATE_PROTEIN_SESSION = """
      UPDATE Protein 
      SET session = ?
      WHERE name = ?
  """
  DELETE_PROTEIN = """   
      DELETE FROM Protein
      WHERE name = ?
  """
  GET_PROTEIN = """SELECT id, structure, chain_sequence_map, session FROM Protein WHERE name = ?"""
  GET_PROTEIN_BY_ID = """SELECT name, structure, chain_sequence_map, session FROM Protein WHERE id = ?"""
  GET_ALL_PROTEINS = """SELECT id, name, structure, chain_sequence_map, session FROM Protein WHERE project_id = ?"""
  # </editor-fold>

  # <editor-fold desc="Ligand statements">
  INSERT_LIGAND = """   
      INSERT INTO Ligand(name, structure, session, project_id)
      VALUES (?, ?, ?, ?)
  """
  UPDATE_LIGAND_STRUCTURE = """
      UPDATE Ligand 
      SET structure = ?
      WHERE name = ?
  """
  UPDATE_LIGAND_SESSION = """
      UPDATE Ligand 
      SET session = ?
      WHERE name = ?
  """
  DELETE_LIGAND = """   
      DELETE FROM Ligand
      WHERE name = ?
  """
  GET_LIGAND = """SELECT id, structure, session, project_id FROM Ligand WHERE name = ?"""
  GET_LIGANDS = """SELECT id, name, structure, session FROM Ligand WHERE project_id = ?"""
  # </editor-fold>

  # <editor-fold desc="Protein-ligand complex statements">
  INSERT_PROTEIN_LIGAND_COMPLEX = """   
      INSERT INTO ProteinLigandComplex(name, protein_id, ligand_id, results, session, project_id)
      VALUES (?, ?, ?, ?, ?, ?)
  """
  UPDATE_PROTEIN_LIGAND_COMPLEX_SESSION = """
      UPDATE ProteinLigandComplex 
      SET session = ?
      WHERE name = ?
  """
  DELETE_PROTEIN_LIGAND_COMPLEX = """   
      DELETE FROM ProteinLigandComplex
      WHERE name = ?
  """
  GET_PROTEIN_LIGAND_COMPLEX = """SELECT id, protein_id, ligand_id, results, session, project_id  FROM ProteinLigandComplex WHERE name = ?"""
  GET_PROTEIN_LIGAND_COMPLEXES = """SELECT id, name, protein_id, ligand_id, results, session FROM ProteinLigandComplex WHERE project_id = ?"""
  # </editor-fold>

  # <editor-fold desc="Protein-protein complex statements">
  INSERT_PROTEIN_PROTEIN_COMPLEX = """   
      INSERT INTO ProteinProteinComplex(name, protein_id, protein_id, results, session, project_id)
      VALUES (?, ?, ?, ?, ?, ?)
  """
  UPDATE_PROTEIN_PROTEIN_COMPLEX_SESSION = """
      UPDATE ProteinProteinComplex 
      SET session = ?
      WHERE name = ?
  """
  DELETE_PROTEIN_PROTEIN_COMPLEX = """   
      DELETE FROM ProteinProteinComplex
      WHERE name = ?
  """
  GET_PROTEIN_PROTEIN_COMPLEX = """SELECT id, protein_1_id, protein_2_id, results, session, project_id  FROM ProteinProteinComplex WHERE name = ?"""
  GET_PROTEIN_PROTEIN_COMPLEXES = """SELECT id, name, protein_1_id, protein_2_id, results, session FROM ProteinProteinComplex WHERE project_id = ?"""
  # </editor-fold>

  # <editor-fold desc="DummyProtein statements">
  INSERT_DUMMY_PROTEIN = """   
      INSERT INTO DummyProtein(name, chain_sequence_map, project_id)
      VALUES (?, ?, ?)
  """
  DELETE_DUMMY_PROTEIN = """   
      DELETE FROM DummyProtein
      WHERE name = ?
  """
  GET_DUMMY_PROTEIN = """SELECT id, chain_sequence_map, project_id FROM DummyProtein WHERE name = ?"""
  GET_DUMMY_PROTEINS = """SELECT id, name, chain_sequence_map FROM DummyProtein WHERE project_id = ?"""
  # </editor-fold>

  # <editor-fold desc="DummyProteinLigandComplex statements">
  INSERT_DUMMY_PROTEIN_LIGAND_COMPLEX = """   
      INSERT INTO DummyProteinLigandComplex(name, project_id)
      VALUES (?, ?)
  """
  DELETE_DUMMY_PROTEIN_LIGAND_COMPLEX = """   
      DELETE FROM DummyProteinLigandComplex
      WHERE name = ?
  """
  GET_DUMMY_PROTEIN_LIGAND_COMPLEX = """SELECT id, project_id FROM DummyProteinLigandComplex WHERE name = ?"""
  GET_DUMMY_PROTEIN_LIGAND_COMPLEXES = """SELECT id, name FROM DummyProteinLigandComplex WHERE project_id = ?"""
  # </editor-fold>

  # <editor-fold desc="DummyProteinLigandComplex statements">
  INSERT_DUMMY_PROTEIN_PROTEIN_COMPLEX = """   
      INSERT INTO DummyProteinProteinComplex(name, project_id)
      VALUES (?, ?)
  """
  DELETE_DUMMY_PROTEIN_PROTEIN_COMPLEX = """   
      DELETE FROM DummyProteinProteinComplex
      WHERE name = ?
  """
  GET_DUMMY_PROTEIN_PROTEIN_COMPLEX = """SELECT id, project_id FROM DummyProteinProteinComplex WHERE name = ?"""
  GET_DUMMY_PROTEIN_PROTEIN_COMPLEXES = """SELECT id, name FROM DummyProteinProteinComplex WHERE project_id = ?"""
  # </editor-fold>

  # <editor-fold desc="Job statements">
  INSERT_JOB = """   
      INSERT INTO Job(name, type, status, dummy_id, dummy_type, project_id)
      VALUES (?, ?, ?, ?, ?, ?)
  """
  UPDATE_JOB_STATUS = """
      UPDATE Job 
      SET status = ?
      WHERE name = ?
  """
  DELETE_JOB = """   
      DELETE FROM Job
      WHERE name = ?
  """
  GET_JOB = """SELECT id, type, status, dummy_id, dummy_type, project_id FROM Job WHERE name = ?"""
  GET_JOBS = """SELECT id, name, type, status, dummy_id, dummy_type FROM Job WHERE project_id = ?"""
  # </editor-fold>
