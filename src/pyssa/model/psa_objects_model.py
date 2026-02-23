"""Module contains PySSAProjectModel — the single model for the project QTreeView.

This model combines sequences, standalone proteins, and protein pairs into one
``QStandardItemModel`` so that a single ``QTreeView`` can display all project
data.  Three permanent top-level section nodes act as visual dividers and
define the fixed display order.

Tree structure
--------------
(invisible root)
├── Sequences                          [TYPE_SECTION]
│   └── <sequence name>                [TYPE_SEQUENCE]
├── Proteins                           [TYPE_SECTION]
│   └── <protein name>                 [TYPE_PROTEIN]
│       ├── Scenes                     [TYPE_HEADER]
│       │   └── <scene name> …         [TYPE_SCENE]
│       └── Chains                     [TYPE_HEADER]
│           └── <chain id>             [TYPE_CHAIN]
│               └── <resi> - <resn>    [TYPE_RESIDUE]
│                   └── <atom name>    [TYPE_ATOM]
└── Protein Pairs                      [TYPE_SECTION]
    └── <pair name>                    [TYPE_PROTEIN_PAIR]
        ├── Scenes                     [TYPE_HEADER]
        │   └── <scene name> …         [TYPE_SCENE]
        ├── <protein 1 name>           [TYPE_PROTEIN]
        │   └── Chains                 [TYPE_HEADER]
        │       └── …
        └── <protein 2 name>           [TYPE_PROTEIN]
            └── Chains                 [TYPE_HEADER]
                └── …

Design rationale
----------------
Section nodes are created in the constructor in the required display order
(Sequences → Proteins → Protein Pairs) so the ordering is always correct
regardless of when or in what order data is added at runtime.

All per-domain logic is delegated to the three focused sub-models, keeping
this class small and easy to follow.  The ``_with_root`` helper temporarily
redirects each sub-model's root node to the matching section node so its
``add_node`` calls build inside the correct section of this model's tree
without any changes being required in the sub-model classes.
"""
import logging

import zmq
from Bio.SeqRecord import SeqRecord

from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.qt import QtCore

from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.model import base_tree_model
from src.pyssa.model.selection_snapshot import SelectionSnapshot
from src.pyssa.model.protein_subtree_mixin import (
  TYPE_PROTEIN,
  TYPE_PROTEIN_PAIR,
  TYPE_SEQUENCE,
  TYPE_SECTION,
  TYPE_HEADER,
  TYPE_SCENE,
  TYPE_CHAIN,
  TYPE_RESIDUE,
  TYPE_ATOM,
  LABEL_SEQUENCES,
  LABEL_PROTEINS,
  LABEL_PROTEIN_PAIRS,
  LABEL_SCENES,
  LABEL_CHAINS,
)
from src.pyssa.model import psa_sequence_model
from src.pyssa.model import psa_protein_model
from src.pyssa.model import psa_protein_pair_model
from src.pyssa.internal.data_structures import protein, protein_pair, project
from src.pyssa.util import enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

__docformat__ = "google"


class PSAObjectsModel(base_tree_model.BaseTreeModel):
  """Single tree model for displaying sequences, proteins, and protein pairs in one QTreeView.

  The model owns three permanent section nodes created in the order
  Sequences → Proteins → Protein Pairs, and delegates all subtree-building
  and management operations to the three specialised sub-models.

  Typical usage
  -------------
  ::

      model = PySSAProjectModel()
      model.build_model(sequences, proteins, protein_pairs, main_socket, socket)
      tree_view.setModel(model)
  """

  def __init__(self) -> None:
    """Constructor.

    Creates the three permanent section nodes in display order and
    initialises the sub-models.
    """
    super().__init__()
    self.create_root_node()

    # Section nodes are appended in the required display order.
    self._sequences_section: QtGui.QStandardItem = self._append_section_node(LABEL_SEQUENCES)
    self._proteins_section: QtGui.QStandardItem = self._append_section_node(LABEL_PROTEINS)
    self._protein_pairs_section: QtGui.QStandardItem = self._append_section_node(LABEL_PROTEIN_PAIRS)

    # Focused sub-models that own the building and management logic.
    # Their root_node attributes are temporarily redirected to the matching
    # section node when a build/add call is made (see _with_root).
    self._sequences_model = psa_sequence_model.PSASequenceModel()
    self._proteins_model = psa_protein_model.PSAProteinModel()
    self._protein_pairs_model = psa_protein_pair_model.PSAProteinPairModel()

  # ------------------------------------------------------------------
  # Public API — building the model
  # ------------------------------------------------------------------

  def build_model(
          self,
          the_project: "project.Project",
  ) -> None:
    """Populate the model from all project data in one call.

    Args:
        the_project: The project that contains all PySSA objects.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if the_project is None:
      logger.error("the_project is None.")
      raise exception.IllegalArgumentError("the_project is None.")

    for sequence in the_project.sequences:
      self.add_sequence(sequence.name)

    for tmp_protein in the_project.proteins:
      self.add_protein(tmp_protein)

    for tmp_pair in the_project.protein_pairs:
      self.add_protein_pair(tmp_pair)

  # ------------------------------------------------------------------
  # Public API — adding items
  # ------------------------------------------------------------------

  def add_sequence(self, sequence: str | SeqRecord) -> None:
    """Add an amino acid sequence under the "Sequences" section.

    Args:
        sequence: Either a string sequence (used as display name and value)
                  or a SeqRecord object (name and sequence extracted automatically).

    Raises:
        exception.IllegalArgumentError: If ``sequence`` is ``None`` or empty.
    """
    if sequence is None:
      logger.error("sequence is None.")
      raise exception.IllegalArgumentError("sequence is None.")

    if isinstance(sequence, str):
      if not sequence:
        logger.error("sequence is an empty string.")
        raise exception.IllegalArgumentError("sequence is an empty string.")
      self._with_root(
        self._sequences_model,
        self._sequences_section,
        lambda: self._sequences_model.add_sequence(sequence),
      )
    elif isinstance(sequence, SeqRecord):
      name = sequence.id or str(sequence.seq)
      seq_str = str(sequence.seq)
      if not seq_str:
        logger.error("SeqRecord has empty sequence.")
        raise exception.IllegalArgumentError("SeqRecord has empty sequence.")
      self._with_root(
        self._sequences_model,
        self._sequences_section,
        lambda: self._sequences_model.add_named_sequence(name, seq_str),
      )
    else:
      logger.error("sequence must be a string or SeqRecord.")
      raise exception.IllegalArgumentError("sequence must be a string or SeqRecord.")

  def add_protein(
          self,
          a_protein: "protein.Protein"
  ) -> None:
    """Add a standalone protein under the "Proteins" section.

    Args:
        a_protein: The protein to add.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    self._with_root(
      self._proteins_model,
      self._proteins_section,
      lambda: self._proteins_model.add_protein_from_protein_object(a_protein),
    )

  def update_protein(
          self,
          a_protein: "protein.Protein"
  ) -> None:
    """Update a standalone protein under the "Proteins" section in-place.

    Args:
        a_protein: The protein to update.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    target_id = a_protein.get_id()
    for row in range(self._proteins_section.rowCount()):
      item = self._proteins_section.child(row)
      stored = item.data(enums.ModelEnum.OBJECT_ROLE)
      if stored is not None and stored.get_id() == target_id:
        self._with_root(
          self._proteins_model,
          self._proteins_section,
          lambda: self._proteins_model.update_protein_node(item, a_protein),
        )
        return
    logger.warning("Protein to update not found in the model.")

  def add_protein_pair(
          self,
          a_protein_pair: "protein_pair.ProteinPair",
  ) -> None:
    """Add a protein pair under the "Protein Pairs" section.

    Args:
        a_protein_pair: The protein pair to add.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    self._with_root(
      self._protein_pairs_model,
      self._protein_pairs_section,
      lambda: self._protein_pairs_model.add_protein_pair(a_protein_pair),
    )

  def add_protein_pair_from_proteins(
          self,
          protein_1_name: str,
          protein_2_name: str,
  ) -> "protein_pair.ProteinPair | None":
    """Create and add a protein pair from two existing standalone proteins.

    This method looks up both proteins by name in the Proteins section,
    creates a deep copy of each for the pair, creates the ProteinPair
    object, and adds it to the model.

    Args:
        protein_1_name: The name of the first protein.
        protein_2_name: The name of the second protein.

    Returns:
        The newly created ProteinPair object, or None if either protein
        was not found.

    Raises:
        exception.IllegalArgumentError: If either protein name is None or empty.
    """
    if not protein_1_name:
      logger.error("protein_1_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
        "protein_1_name is either None or an empty string."
      )
    if not protein_2_name:
      logger.error("protein_2_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
        "protein_2_name is either None or an empty string."
      )

    # Look up both proteins
    protein_1 = self.get_protein_by_name(protein_1_name)
    protein_2 = self.get_protein_by_name(protein_2_name)

    if protein_1 is None:
      logger.error(f"Protein '{protein_1_name}' not found in Proteins section.")
      return None
    if protein_2 is None:
      logger.error(f"Protein '{protein_2_name}' not found in Proteins section.")
      return None

    # Create the protein pair (ProteinPair constructor handles deep copying)
    new_pair = protein_pair.ProteinPair(protein_1, protein_2)

    # Add to model
    self.add_protein_pair(new_pair)

    return new_pair

  def add_temporary_protein(
          self,
          a_protein: "protein.Protein"
  ) -> None:
    """Add a standalone protein minimally (without PyMOL query).
    
    Delegates to the Proteins sub-model to append a simple layout
    of just the protein and its chains for temporary UIs.

    Args:
        a_protein: The protein to add.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    self._with_root(
      self._proteins_model,
      self._proteins_section,
      lambda: self._proteins_model.add_temporary_protein(a_protein),
    )

  # ------------------------------------------------------------------
  # Public API — removing items
  # ------------------------------------------------------------------

  def remove_sequence(self, sequence: str) -> None:
    """Remove a sequence from the "Sequences" section.

    Args:
        sequence: The sequence string to remove.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None`` or empty.
        ValueError: If the sequence is not found in the model.
    """
    if not sequence:
      logger.error("sequence is either None or an empty string.")
      raise exception.IllegalArgumentError("sequence is either None or an empty string.")

    # Find the sequence in the Sequences section
    for row in range(self._sequences_section.rowCount()):
      item = self._sequences_section.child(row)
      if item.data(enums.ModelEnum.OBJECT_ROLE) == sequence:
        self._sequences_section.removeRow(row)
        return

    logger.warning(f"Sequence not found in model during removal.")

  def remove_protein(self, a_protein: "protein.Protein") -> None:
    """Remove a protein from the "Proteins" section.

    Args:
        a_protein: The protein to remove.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the protein is not found in the model.
    """
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    target_name = a_protein.get_molecule_object()
    target_id = a_protein.get_id()
    for row in range(self._proteins_section.rowCount()):
      item = self._proteins_section.child(row)
      stored = item.data(enums.ModelEnum.OBJECT_ROLE)
      if stored is not None and stored.get_id() == target_id:
        self._proteins_section.removeRow(row)
        return

    logger.warning("Protein '%s' not found in model during removal.", target_name)

  def remove_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair") -> None:
    """Remove a protein pair from the "Protein Pairs" section.

    Args:
        a_protein_pair: The protein pair to remove.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the protein pair is not found in the model.
    """
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    target_name = a_protein_pair.name
    target_id = a_protein_pair.get_id()
    for row in range(self._protein_pairs_section.rowCount()):
      item = self._protein_pairs_section.child(row)
      stored = item.data(enums.ModelEnum.OBJECT_ROLE)
      if stored is not None and stored.get_id() == target_id:
        self._protein_pairs_section.removeRow(row)
        return

    logger.warning("Protein pair '%s' not found in model during removal.", target_name)

  # ------------------------------------------------------------------
  # Public API — scene management (delegates to protein/pair sub-models)
  # ------------------------------------------------------------------

  def add_scene(
          self,
          scene_name: str,
          target: "protein.Protein | protein_pair.ProteinPair",
  ) -> None:
    """Add a scene to a protein or protein pair.

    Args:
        scene_name: The name of the scene to add.
        target: The protein or protein pair to add the scene to.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None`` or scene_name is empty.
        ValueError: If the target object is not found in the model.
    """
    if not scene_name:
      logger.error("scene_name is either None or an empty string.")
      raise exception.IllegalArgumentError("scene_name is either None or an empty string.")
    if target is None:
      logger.error("target is None.")
      raise exception.IllegalArgumentError("target is None.")

    # Create the scene item
    scene_item = QtGui.QStandardItem(scene_name)
    scene_item.setData(TYPE_SCENE, enums.ModelEnum.TYPE_ROLE)

    # Find the target in the model
    target_index = self._find_object_index(target)
    if target_index is None or not target_index.isValid():
      logger.error(f"Target object not found in model.")
      raise ValueError("Target object not found in model.")

    # Find the Scenes header and add the scene item directly
    target_item = self.itemFromIndex(target_index)
    scenes_header = self._find_scenes_header_item(target_item)
    if scenes_header is None:
      logger.error("Could not find Scenes header for target.")
      raise ValueError("Could not find Scenes header for target.")

    scenes_header.appendRow(scene_item)

  def remove_scene(
          self,
          scene_name: str,
          target: "protein.Protein | protein_pair.ProteinPair",
  ) -> None:
    """Remove a scene from a protein or protein pair.

    Args:
        scene_name: The name of the scene to remove.
        target: The protein or protein pair to remove the scene from.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None`` or scene_name is empty.
        ValueError: If the target object or scene is not found in the model.
    """
    if not scene_name:
      logger.error("scene_name is either None or an empty string.")
      raise exception.IllegalArgumentError("scene_name is either None or an empty string.")
    if target is None:
      logger.error("target is None.")
      raise exception.IllegalArgumentError("target is None.")

    # Find the target in the model
    target_index = self._find_object_index(target)
    if target_index is None or not target_index.isValid():
      logger.error(f"Target object not found in model.")
      raise ValueError("Target object not found in model.")

    # Find the Scenes header
    target_item = self.itemFromIndex(target_index)
    scenes_header = self._find_scenes_header_item(target_item)
    if scenes_header is None:
      logger.error("Could not find Scenes header for target.")
      raise ValueError("Could not find Scenes header for target.")

    # Find and remove the scene by name
    for row in range(scenes_header.rowCount()):
      scene_item = scenes_header.child(row)
      if scene_item.text() == scene_name:
        scenes_header.removeRow(row)
        return

    logger.error(f"Scene '{scene_name}' not found for target.")
    raise ValueError(f"Scene '{scene_name}' not found for target.")

  def check_if_scratch_scene_exists(
          self, target: "protein.Protein | protein_pair.ProteinPair"
  ) -> bool:
    """Return ``True`` if a ``_scratch_`` scene exists for the given protein or protein pair.

    Args:
        target: The protein or protein pair to check.

    Raises:
        exception.IllegalArgumentError: If ``target`` is ``None``.
        ValueError: If the target object is not found in the model.
    """
    if target is None:
      logger.error("target is None.")
      raise exception.IllegalArgumentError("target is None.")

    # Find the target in the model
    target_index = self._find_object_index(target)
    if target_index is None or not target_index.isValid():
      logger.error(f"Target object not found in model.")
      raise ValueError("Target object not found in model.")

    sub_model = self._scene_capable_sub_model_for_index(target_index)
    return sub_model.check_if_scratch_scene_exists(target_index)

  # ------------------------------------------------------------------
  # Public API — protein lookup
  # ------------------------------------------------------------------

  def get_protein_by_name(self, protein_name: str) -> "protein.Protein | None":
    """Search for a protein by name in the Proteins section.

    Args:
        protein_name: The name of the protein to find.

    Returns:
        The Protein object if found, None otherwise.

    Raises:
        exception.IllegalArgumentError: If protein_name is None or empty.
    """
    if not protein_name:
      logger.error("protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
        "protein_name is either None or an empty string."
      )

    # Iterate through the Proteins section
    for row in range(self._proteins_section.rowCount()):
      protein_item = self._proteins_section.child(row)
      if protein_item.data(QtCore.Qt.ItemDataRole.DisplayRole) == protein_name:
        return protein_item.data(enums.ModelEnum.OBJECT_ROLE)
    return None

  # ------------------------------------------------------------------
  # Public API — sequence data access
  # ------------------------------------------------------------------

  def get_sequence_string(self, a_model_index: QtCore.QModelIndex) -> str:
    """Return the amino acid sequence string stored at *a_model_index*.

    Args:
        a_model_index: Index of a sequence node.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
        ValueError: If the node is not a sequence node.
    """
    return self._sequences_model.get_sequence_string(a_model_index)

  # ------------------------------------------------------------------
  # Public API — section index accessors
  # ------------------------------------------------------------------

  def get_sequences_section_index(self) -> QtCore.QModelIndex:
    """Return the model index of the "Sequences" section node."""
    return self.indexFromItem(self._sequences_section)

  def get_proteins_section_index(self) -> QtCore.QModelIndex:
    """Return the model index of the "Proteins" section node."""
    return self.indexFromItem(self._proteins_section)

  def get_protein_pairs_section_index(self) -> QtCore.QModelIndex:
    """Return the model index of the "Protein Pairs" section node."""
    return self.indexFromItem(self._protein_pairs_section)

  # ------------------------------------------------------------------
  # Public API — PyMOL selection string construction
  # ------------------------------------------------------------------

  def construct_selection_string(self, a_model_index: QtCore.QModelIndex) -> str:
    """Build a PyMOL selection expression from a tree node.

    Walks the node's ``TYPE_ROLE`` and its ancestors to produce a valid
    PyMOL selection string.  The returned string is suitable for use with
    ``cmd.select("sele", <string>)``.

    Supported node types and their resulting formats:

    * ``TYPE_PROTEIN`` / ``TYPE_PROTEIN_PAIR``
        → ``"<molecule_object_name>"``
    * ``TYPE_CHAIN``
        → ``"(<protein> and chain <letter>)"``
    * ``TYPE_RESIDUE``
        → ``"(<protein> and chain <letter> and resi <number>)"``
    * ``TYPE_ATOM``
        → ``"(<protein> and chain <letter> and resi <number> and name <atom>)"``

    All other node types (sections, headers, scenes, sequences) return an
    empty string because they have no direct PyMOL 3-D representation.

    Args:
        a_model_index: Index of a node in this model.

    Returns:
        A PyMOL selection expression, or ``""`` for non-selectable node types.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
    """
    if a_model_index is None:
      logger.error("a_model_index is None.")
      raise exception.IllegalArgumentError("a_model_index is None.")

    node_type = a_model_index.data(enums.ModelEnum.TYPE_ROLE)

    if node_type == TYPE_PROTEIN:
      return f"/{a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole)}"

    if node_type == TYPE_PROTEIN_PAIR:
      pair_object = a_model_index.data(enums.ModelEnum.OBJECT_ROLE)
      if pair_object is not None:
        return (
          f"/{pair_object.protein_1.get_molecule_object()} or "
          f"/{pair_object.protein_2.get_molecule_object()}"
        )
      return f"/{a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole)}"

    if node_type == TYPE_CHAIN:
      # Chain → parent is "Chains" header → parent is protein node.
      protein_index = a_model_index.parent().parent()
      protein_name = protein_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      chain_letter = a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      return f"/{protein_name}//{chain_letter}"

    if node_type == TYPE_RESIDUE:
      chain_index = a_model_index.parent()
      protein_index = chain_index.parent().parent()
      protein_name = protein_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      chain_letter = chain_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      resi = a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole).split(" - ")[0]
      return f"/{protein_name}//{chain_letter}/{resi}"

    if node_type == TYPE_ATOM:
      residue_index = a_model_index.parent()
      chain_index = residue_index.parent()
      protein_index = chain_index.parent().parent()
      protein_name = protein_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      chain_letter = chain_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      resi = residue_index.data(QtCore.Qt.ItemDataRole.DisplayRole).split(" - ")[0]
      atom_name = a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      return f"/{protein_name}//{chain_letter}/{resi}/{atom_name}"

    return ""

  # ------------------------------------------------------------------
  # Public API — selection resolution
  # ------------------------------------------------------------------

  def resolve_selection(
          self, selected_indexes: list[QtCore.QModelIndex]
  ) -> SelectionSnapshot:
    """Parse a list of selected QModelIndex elements into an immutable SelectionSnapshot.
    
    This method iterates over the provided indices, extracting the underlying objects based
    on the TYPE_ROLE. It also resolves parent references recursively to ensure that if a 
    child (e.g., atom) is selected, its parent Protein is also added to the snapshot so that
    downstream logic works seamlessly.
    
    Args:
        selected_indexes: A list of selected model indices from the view.
        
    Returns:
        An immutable SelectionSnapshot populated with the exact semantic selection state.
    """
    raw_sequences = set()
    raw_standalone_proteins = set()
    raw_protein_pair_children = set()
    raw_protein_pairs = set()
    raw_chains = set()
    raw_residues = set()
    raw_atoms = set()
    raw_scenes = set()

    # We need a helper to safely extract the upper-level Protein and its parent TYPE from any node in the tree.
    def _resolve_protein_and_parent_type(index: QtCore.QModelIndex) -> tuple["protein.Protein | None", int | None]:
        current_idx = index
        while current_idx.isValid():
            node_type = current_idx.data(enums.ModelEnum.TYPE_ROLE)
            if node_type == TYPE_PROTEIN:
                parent_idx = current_idx.parent()
                parent_type = parent_idx.data(enums.ModelEnum.TYPE_ROLE) if parent_idx.isValid() else None
                return current_idx.data(enums.ModelEnum.OBJECT_ROLE), parent_type
            current_idx = current_idx.parent()
        return None, None

    for index in selected_indexes:
        if not index.isValid():
            continue

        node_type = index.data(enums.ModelEnum.TYPE_ROLE)
        obj_data = index.data(enums.ModelEnum.OBJECT_ROLE)

        if node_type == TYPE_SEQUENCE:
            # For sequences, OBJECT_ROLE stores the string.
            raw_sequences.add(obj_data)
        elif node_type == TYPE_PROTEIN:
            parent_idx = index.parent()
            parent_type = parent_idx.data(enums.ModelEnum.TYPE_ROLE) if parent_idx.isValid() else None
            if parent_type == TYPE_PROTEIN_PAIR:
                raw_protein_pair_children.add(obj_data)
            else:
                raw_standalone_proteins.add(obj_data)
        elif node_type == TYPE_PROTEIN_PAIR:
            raw_protein_pairs.add(obj_data)
        elif node_type == TYPE_CHAIN:
            raw_chains.add(obj_data)
        elif node_type == TYPE_RESIDUE:
            raw_residues.add(obj_data)
        elif node_type == TYPE_ATOM:
            raw_atoms.add(obj_data)
        elif node_type == TYPE_SCENE:
            # Sequences don't have scenes; the tree text defines the scene.
            raw_scenes.add(index.data(QtCore.Qt.ItemDataRole.DisplayRole))

        # IMPORTANT: Regardless of what child object is selected (chain, residue, atom), 
        # we resolve the host Protein and place it into the appropriate set (standalone vs pair-child). 
        # This gives the snapshot its power: a user can select a single atom, 
        # and downstream code knows immediately which protein it belongs to and its context.
        if node_type in (TYPE_CHAIN, TYPE_RESIDUE, TYPE_ATOM, TYPE_SCENE):
            parent_prot, parent_prot_parent_type = _resolve_protein_and_parent_type(index)
            if parent_prot:
                if parent_prot_parent_type == TYPE_PROTEIN_PAIR:
                    raw_protein_pair_children.add(parent_prot)
                else:
                    raw_standalone_proteins.add(parent_prot)

    # Build the combined PyMOL selection string from all selectable indexes.
    selection_fragments: list[str] = []
    for index in selected_indexes:
        if not index.isValid():
            continue
        fragment = self.construct_selection_string(index)
        if fragment:
            selection_fragments.append(fragment)
    pymol_selection = " or ".join(selection_fragments)

    return SelectionSnapshot(
        raw_sequences=raw_sequences,
        raw_standalone_proteins=raw_standalone_proteins,
        raw_protein_pair_children=raw_protein_pair_children,
        raw_protein_pairs=raw_protein_pairs,
        raw_chains=raw_chains,
        raw_residues=raw_residues,
        raw_atoms=raw_atoms,
        raw_scenes=raw_scenes,
        pymol_selection_string=pymol_selection,
    )

  # ------------------------------------------------------------------
  # Public API — reverse selection (PyMOL → tree indexes)
  # ------------------------------------------------------------------

  def find_indexes_for_pymol_atoms(
          self,
          a_chempy_model: "Indexed",
  ) -> list[QtCore.QModelIndex]:
    """Find tree-view indexes that match atoms from a PyMOL selection.

    Given a chempy ``Indexed`` model (typically from ``cmd.get_model("sele")``),
    this method walks the Proteins and Protein Pairs sections of the tree to
    locate the corresponding chain → residue → atom nodes.

    **Upward collapsing:** if every atom of a residue is selected, the method
    returns the residue's ``QModelIndex`` instead of the individual atom
    indexes.  This mirrors how the legacy demo rolled up full-residue
    selections.

    Args:
        a_chempy_model: The chempy ``Indexed`` model containing the selected
            atoms.  May be ``None`` or have an empty ``atom`` list, in which
            case an empty list is returned.

    Returns:
        A list of ``QModelIndex`` objects suitable for passing to a
        ``QItemSelectionModel.select()`` call.
    """
    if a_chempy_model is None:
      return []
    atoms_list = getattr(a_chempy_model, "atom", [])
    if not atoms_list:
      return []

    # Step 1 — Build a hierarchy map from the incoming atoms:
    #   { object_name → { chain_id → { (resi, resn) → { atom_names } } } }
    hierarchy: dict[str, dict[str, dict[tuple[str, str], set[str]]]] = {}
    for atom in atoms_list:
      obj_name = getattr(atom, "model", None) or ""
      chain_id = getattr(atom, "chain", None) or ""
      resi = str(getattr(atom, "resi", ""))
      resn = str(getattr(atom, "resn", ""))
      atom_name = str(getattr(atom, "name", ""))
      (
        hierarchy
        .setdefault(obj_name, {})
        .setdefault(chain_id, {})
        .setdefault((resi, resn), set())
        .add(atom_name)
      )

    # Step 2 — Walk sections that can contain protein nodes.
    result_indexes: list[QtCore.QModelIndex] = []
    sections = [self._proteins_section, self._protein_pairs_section]

    for section in sections:
      for sec_row in range(section.rowCount()):
        top_item = section.child(sec_row)
        top_type = top_item.data(enums.ModelEnum.TYPE_ROLE)

        # Collect protein items to inspect.
        protein_items: list[QtGui.QStandardItem] = []
        if top_type == TYPE_PROTEIN:
          protein_items.append(top_item)
        elif top_type == TYPE_PROTEIN_PAIR:
          # A pair's children include headers and two protein nodes.
          for pair_child_row in range(top_item.rowCount()):
            pair_child = top_item.child(pair_child_row)
            if pair_child.data(enums.ModelEnum.TYPE_ROLE) == TYPE_PROTEIN:
              protein_items.append(pair_child)

        for protein_item in protein_items:
          protein_display = protein_item.text()
          if protein_display not in hierarchy:
            continue
          chains_for_protein = hierarchy[protein_display]

          # Find the "Chains" header under this protein node.
          chains_header = self._find_chains_header_item(protein_item)
          if chains_header is None:
            continue

          self._match_chain_hierarchy(
            chains_header, chains_for_protein, result_indexes
          )

    return result_indexes

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  def _append_section_node(self, a_label: str) -> QtGui.QStandardItem:
    """Append a permanent, non-selectable top-level section divider and return it."""
    section_item = QtGui.QStandardItem(a_label)
    section_item.setData(TYPE_SECTION, enums.ModelEnum.TYPE_ROLE)
    section_item.setFlags(
      section_item.flags() & ~QtCore.Qt.ItemFlag.ItemIsSelectable
    )
    self.root_node.appendRow(section_item)
    return section_item

  def _with_root(
          self,
          a_sub_model: base_tree_model.BaseTreeModel,
          a_section_node: QtGui.QStandardItem,
          a_callable,
  ) -> None:
    """Run *a_callable* with *a_sub_model*'s root temporarily set to *a_section_node*.

    This lets the sub-model's ``add_node`` calls build directly inside the
    correct section of this model's item tree, with no changes needed in
    the sub-model classes.  The original root is always restored, even if
    *a_callable* raises.

    Args:
        a_sub_model: The sub-model whose root to redirect.
        a_section_node: The section node to use as the temporary root.
        a_callable: Zero-argument callable to invoke with the root redirected.
    """
    original_root = a_sub_model.root_node
    a_sub_model.root_node = a_section_node
    try:
      a_callable()
    finally:
      a_sub_model.root_node = original_root

  def _scene_capable_sub_model_for_index(
          self, a_model_index: QtCore.QModelIndex
  ) -> "psa_protein_model.PSAProteinModel | psa_protein_pair_model.PSAProteinPairModel":
    """Return the protein or protein-pair sub-model that owns *a_model_index*.

    Sequences have no scenes, so this raises ``ValueError`` if the index
    belongs to the Sequences section.

    Raises:
        ValueError: If the section is "Sequences" or cannot be found in
            the ancestor chain.
    """
    item = self.itemFromIndex(a_model_index)
    while item is not None:
      if item.data(enums.ModelEnum.TYPE_ROLE) == TYPE_SECTION:
        label = item.text()
        if label == LABEL_SEQUENCES:
          raise ValueError(
            "Scene operations are not supported for the Sequences section."
          )
        if label == LABEL_PROTEINS:
          return self._proteins_model
        if label == LABEL_PROTEIN_PAIRS:
          return self._protein_pairs_model
      item = item.parent()
    raise ValueError("Could not determine the section for the given index.")

  def _find_object_index(
          self, target_object: "protein.Protein | protein_pair.ProteinPair"
  ) -> QtCore.QModelIndex | None:
    """Find the QModelIndex for a given protein or protein pair object.

    Args:
        target_object: The protein or protein pair to find.

    Returns:
        The QModelIndex if found, None otherwise.
    """
    # Determine which section to search
    if isinstance(target_object, protein.Protein):
      # Search in Proteins section
      for row in range(self._proteins_section.rowCount()):
        item = self._proteins_section.child(row)
        if item.data(enums.ModelEnum.OBJECT_ROLE) is target_object:
          return self.indexFromItem(item)
    elif isinstance(target_object, protein_pair.ProteinPair):
      # Search in Protein Pairs section
      for row in range(self._protein_pairs_section.rowCount()):
        item = self._protein_pairs_section.child(row)
        if item.data(enums.ModelEnum.OBJECT_ROLE) is target_object:
          return self.indexFromItem(item)
    return None

  def _find_scenes_header_item(
          self, protein_or_pair_item: QtGui.QStandardItem
  ) -> QtGui.QStandardItem | None:
    """Find the Scenes header item under a protein or protein pair item.

    Args:
        protein_or_pair_item: The protein or protein pair item.

    Returns:
        The Scenes header item if found, None otherwise.
    """
    # The Scenes header is always the first child of a protein or protein pair node
    for row in range(protein_or_pair_item.rowCount()):
      child = protein_or_pair_item.child(row)
      if (child.data(enums.ModelEnum.TYPE_ROLE) == TYPE_HEADER and
          child.text() == LABEL_SCENES):
        return child
    return None

  def _find_chains_header_item(
          self, protein_item: QtGui.QStandardItem
  ) -> QtGui.QStandardItem | None:
    """Find the Chains header item under a protein item.

    Args:
        protein_item: The protein item to search within.

    Returns:
        The Chains header ``QStandardItem`` if found, ``None`` otherwise.
    """
    for row in range(protein_item.rowCount()):
      child = protein_item.child(row)
      if (child.data(enums.ModelEnum.TYPE_ROLE) == TYPE_HEADER and
          child.text() == LABEL_CHAINS):
        return child
    return None

  def _match_chain_hierarchy(
          self,
          chains_header: QtGui.QStandardItem,
          chains_map: dict[str, dict[tuple[str, str], set[str]]],
          result_indexes: list[QtCore.QModelIndex],
  ) -> None:
    """Match a PyMOL atom hierarchy against chain→residue→atom tree nodes.

    For each chain in *chains_map*, this method locates the corresponding
    chain node under *chains_header*, then walks its residue children.
    If **all** atoms of a residue are present in the map, the residue's
    ``QModelIndex`` is appended (upward collapse); otherwise only matching
    atom indexes are appended.

    Args:
        chains_header: The "Chains" header ``QStandardItem`` to walk.
        chains_map: Per-chain data ``{chain_id → {(resi, resn) → {atom_names}}}``.
        result_indexes: Accumulator list; matched indexes are appended in place.
    """
    for chain_row in range(chains_header.rowCount()):
      chain_item = chains_header.child(chain_row)
      chain_letter = chain_item.text()
      if chain_letter not in chains_map:
        continue
      residues_map = chains_map[chain_letter]

      for resi_row in range(chain_item.rowCount()):
        residue_item = chain_item.child(resi_row)
        residue_label = residue_item.text()

        # Parse residue label "resi - resn" back into the (resi, resn) key.
        parts = residue_label.split(" - ", 1)
        if len(parts) != 2:
          continue
        resi_key = (parts[0], parts[1])
        if resi_key not in residues_map:
          continue

        selected_atoms = residues_map[resi_key]
        total_atoms_in_tree = residue_item.rowCount()

        if len(selected_atoms) >= total_atoms_in_tree and total_atoms_in_tree > 0:
          # All atoms of this residue are selected → collapse to residue node.
          result_indexes.append(self.indexFromItem(residue_item))
        else:
          # Only some atoms are selected → append individual atom indexes.
          for atom_row in range(total_atoms_in_tree):
            atom_item = residue_item.child(atom_row)
            if atom_item.text() in selected_atoms:
              result_indexes.append(self.indexFromItem(atom_item))
