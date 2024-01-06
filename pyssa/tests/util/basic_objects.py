import os
import pathlib

from pyssa.io_pyssa import path_util
from pyssa.tests.util import io_test
from pyssa.internal.data_structures import protein


def create_protein_obj_from_pdb_id(pdb_id: str) -> "protein.Protein":
    tmp_filepath: pathlib.Path = pathlib.Path(io_test.download_pdb_file_by_id(pdb_id, pdb_id))
    tmp_protein: "protein.Protein" = protein.Protein(pdb_id, path_util.FilePath(tmp_filepath))
    if os.path.exists(tmp_filepath):
        os.remove(tmp_filepath)
    return tmp_protein


def create_protein_obj_from_pdb_file(url: str, filename: str) -> "protein.Protein":
    tmp_filepath: pathlib.Path = pathlib.Path(io_test.download_file_from_url(url, filename))
    tmp_protein: "protein.Protein" = protein.Protein(filename.replace(".pdb", ""), path_util.FilePath(tmp_filepath))
    if os.path.exists(tmp_filepath):
        os.remove(tmp_filepath)
    return tmp_protein
