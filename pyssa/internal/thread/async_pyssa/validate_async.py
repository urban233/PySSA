import os
import pathlib

from pyssa.io_pyssa import bio_data
from pyssa.util import constants


def validate_add_protein_view_input(the_entered_text: str, placeholder: int):
    # checks if a pdb id was entered
    if len(the_entered_text) == 4:
        # pdb id is used
        pdb_id = the_entered_text.upper()
        tmp_filepath: str = str(pathlib.Path(f"{constants.SCRATCH_DIR}/{pdb_id}.pdb"))
        bio_data.download_pdb_file(pdb_id, tmp_filepath)
        if os.path.exists(tmp_filepath):
            os.remove(tmp_filepath)
            return 1, True, pdb_id
        else:
            return 1, False, pdb_id
    else:
        if os.path.exists(the_entered_text):
            return 2, True, pathlib.Path(the_entered_text).name.replace(".pdb", "")
        else:
            return 2, False, pathlib.Path(the_entered_text).name.replace(".pdb", "")
