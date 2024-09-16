import pathlib
from dataclasses import dataclass
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Draw, AllChem


@dataclass
class LigandPreview:
  """A dataclass storing a preview of a ligand object."""
  structure: Chem.Mol
  filepath: pathlib.Path
  name: str = ""
  image_filepath: Optional[pathlib.Path] = None

  def create_molecule_image(self) -> None:
    """Creates an image of the ligand."""
    AllChem.Compute2DCoords(self.structure)
    img = Draw.MolToFile(self.structure, "test.png", size=(400, 400))
    self.image_filepath = pathlib.Path("test.png")

  def get_ligand_name(self) -> str:
    """Returns the name of the ligand based on either the molecule or the filepath."""
    if self.structure.GetProp("_Name") == "":
      return self.filepath.name
    return self.structure.GetProp("_Name")
