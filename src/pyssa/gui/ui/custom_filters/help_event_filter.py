from src.pyssa.gui.qt import QtCore, QtGui
from src.pyssa.gui.ui.views import help_panel
import logging

logger = logging.getLogger(__name__)

# Import type constants for tree item identification.
# These are resolved lazily to avoid circular imports.
_TREE_TYPES_LOADED = False
_TYPE_SECTION = None
_TYPE_SEQUENCE = None
_TYPE_PROTEIN = None
_TYPE_PROTEIN_PAIR = None
_TYPE_HEADER = None
_TYPE_SCENE = None
_TYPE_CHAIN = None
_TYPE_RESIDUE = None
_TYPE_ATOM = None
_LABEL_SEQUENCES = None
_LABEL_PROTEINS = None
_LABEL_PROTEIN_PAIRS = None
_LABEL_SCENES = None
_LABEL_CHAINS = None


def _load_tree_types():
  global _TREE_TYPES_LOADED
  global _TYPE_SECTION, _TYPE_SEQUENCE, _TYPE_PROTEIN, _TYPE_PROTEIN_PAIR
  global _TYPE_HEADER, _TYPE_SCENE, _TYPE_CHAIN, _TYPE_RESIDUE, _TYPE_ATOM
  global _LABEL_SEQUENCES, _LABEL_PROTEINS, _LABEL_PROTEIN_PAIRS
  global _LABEL_SCENES, _LABEL_CHAINS
  if _TREE_TYPES_LOADED:
    return
  try:
    from src.pyssa.model.protein_subtree_mixin import (
      TYPE_SECTION, TYPE_SEQUENCE, TYPE_PROTEIN, TYPE_PROTEIN_PAIR,
      TYPE_HEADER, TYPE_SCENE, TYPE_CHAIN, TYPE_RESIDUE, TYPE_ATOM,
      LABEL_SEQUENCES, LABEL_PROTEINS, LABEL_PROTEIN_PAIRS,
      LABEL_SCENES, LABEL_CHAINS,
    )
    from src.pyssa.util import enums as _enums
    _TYPE_SECTION = TYPE_SECTION
    _TYPE_SEQUENCE = TYPE_SEQUENCE
    _TYPE_PROTEIN = TYPE_PROTEIN
    _TYPE_PROTEIN_PAIR = TYPE_PROTEIN_PAIR
    _TYPE_HEADER = TYPE_HEADER
    _TYPE_SCENE = TYPE_SCENE
    _TYPE_CHAIN = TYPE_CHAIN
    _TYPE_RESIDUE = TYPE_RESIDUE
    _TYPE_ATOM = TYPE_ATOM
    _LABEL_SEQUENCES = LABEL_SEQUENCES
    _LABEL_PROTEINS = LABEL_PROTEINS
    _LABEL_PROTEIN_PAIRS = LABEL_PROTEIN_PAIRS
    _LABEL_SCENES = LABEL_SCENES
    _LABEL_CHAINS = LABEL_CHAINS
    _TREE_TYPES_LOADED = True
  except Exception as exc:
    logger.warning(f"Could not load tree type constants: {exc}")


class HelpEventFilter(QtCore.QObject):
  def __init__(self, help_browser, help_map, help_panel=None):
    super().__init__()
    self.help_browser = help_browser
    self.help_map = help_map
    self.help_panel: "help_panel.HelpPanel" = help_panel

  def eventFilter(self, obj, event):
    if event.type() == QtCore.QEvent.Type.Enter:
      object_name = obj.objectName()
      logger.debug(f"Enter event on object: {object_name}")
      help_text = self.help_map.get(object_name, "No help available")
      logger.debug(f"Setting help text: {help_text[:50]}...")
      self.help_browser.setHtml(help_text)

    elif event.type() == QtCore.QEvent.Type.Leave:
      logger.debug(f"Leave event on object: {obj.objectName()}")
      self.help_browser.clear()

    return super().eventFilter(obj, event)

  def handle_menu_action_hovered(self, action: QtGui.QAction):
    """Handle when a menu action is hovered."""
    object_name = action.objectName()
    logger.debug(f"Menu action hovered: {object_name}")
    help_text = self.help_map.get(object_name, "")
    if help_text:
      self.help_browser.setHtml(help_text)
      self.help_panel.set_help_panel_header()

  def handle_menu_about_to_hide(self):
    """Handle when  menu is about to close."""
    logger.debug("Menu about to hide, clearing help")
    self.help_browser.clear()

  def handle_tree_item_entered(self, index: QtCore.QModelIndex):
    """Show help text for the tree item the cursor just entered.

    Called by connecting the QTreeView's ``entered(QModelIndex)`` signal
    to this slot.  The help key is resolved from the item's TYPE_ROLE
    (and, for section/header nodes, from the display text as well).
    """
    _load_tree_types()
    if not _TREE_TYPES_LOADED or not index.isValid():
      return

    try:
      from src.pyssa.util.enums import ModelEnum
      node_type = index.data(ModelEnum.TYPE_ROLE)
      display_text = index.data(QtCore.Qt.ItemDataRole.DisplayRole)

      help_key = None

      if node_type == _TYPE_SECTION:
        if display_text == _LABEL_SEQUENCES:
          help_key = "section_sequences"
        elif display_text == _LABEL_PROTEINS:
          help_key = "section_proteins"
        elif display_text == _LABEL_PROTEIN_PAIRS:
          help_key = "section_protein_pairs"

      elif node_type == _TYPE_SEQUENCE:
        help_key = "sequence_item"

      elif node_type == _TYPE_PROTEIN:
        # Distinguish standalone protein from protein inside a protein pair.
        parent_idx = index.parent()
        parent_type = parent_idx.data(ModelEnum.TYPE_ROLE) if parent_idx.isValid() else None
        if parent_type == _TYPE_PROTEIN_PAIR:
          help_key = "protein_item_in_pair"
        else:
          help_key = "protein_item"

      elif node_type == _TYPE_PROTEIN_PAIR:
        help_key = "protein_pair_item"

      elif node_type == _TYPE_HEADER:
        if display_text == _LABEL_SCENES:
          help_key = "header_scenes"
        elif display_text == _LABEL_CHAINS:
          help_key = "header_chains"

      elif node_type == _TYPE_SCENE:
        help_key = "scene_item"

      elif node_type == _TYPE_CHAIN:
        help_key = "chain_item"

      elif node_type == _TYPE_RESIDUE:
        help_key = "residue_item"

      elif node_type == _TYPE_ATOM:
        help_key = "atom_item"

      if help_key:
        help_text = self.help_map.get(help_key, "")
        if help_text:
          logger.debug(f"Tree item entered, key={help_key}")
          self.help_browser.setHtml(help_text)
        else:
          self.help_browser.clear()
      else:
        self.help_browser.clear()

    except Exception as exc:
      logger.warning(f"handle_tree_item_entered: {exc}")

  def handle_tree_left(self):
    """Clear the help panel when the cursor leaves the tree viewport."""
    self.help_browser.clear()
