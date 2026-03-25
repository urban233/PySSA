"""
Utility module to load hover help texts from external HTML files.

This module reads help content from HTML files stored in the
docs/pyssa-documentation/docs/help/help directory structure.
"""
import logging
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# Cache for loaded help texts (loaded once per application lifetime)
_help_map_cache: Optional[Dict[str, str]] = None

# Color definitions for generating help text programmatically
_COLOR_DEFINITIONS = {
  # Reds
  "color_red": ("red", "#ff0000", "pure red"),
  "color_tv_red": ("tv_red", "#ff3333", "TV red"),
  "color_salmon": ("salmon", "#ff9999", "salmon"),
  "color_raspberry": ("raspberry", "#b24c66", "raspberry"),
  # Greens
  "color_green": ("green", "#00ff00", "pure green"),
  "color_tv_green": ("tv_green", "#33ff33", "TV green"),
  "color_palegreen": ("palegreen", "#a5e5a5", "pale green"),
  "color_forest": ("forest", "#339933", "forest green"),
  # Blues
  "color_blue": ("blue", "#0000ff", "pure blue"),
  "color_tv_blue": ("tv_blue", "#4c4cff", "TV blue"),
  "color_lightblue": ("lightblue", "#bfbfff", "light blue"),
  "color_skyblue": ("skyblue", "#337fcc", "sky blue"),
  # Yellows
  "color_yellow": ("yellow", "#ffff00", "pure yellow"),
  "color_tv_yellow": ("tv_yellow", "#ffff33", "TV yellow"),
  "color_paleyellow": ("paleyellow", "#ffff7f", "pale yellow"),
  "color_sand": ("sand", "#b78c4c", "sand"),
  # Magentas
  "color_magenta": ("magenta", "#ff00ff", "magenta"),
  "color_purple": ("purple", "#bf00bf", "purple"),
  "color_pink": ("pink", "#ffa5d8", "pink"),
  "color_hotpink": ("hotpink", "#ff007f", "hot pink"),
  # Cyans
  "color_cyan": ("cyan", "#00ffff", "cyan"),
  "color_aquamarine": ("aquamarine", "#7fffff", "aquamarine"),
  "color_palecyan": ("palecyan", "#ccffff", "pale cyan"),
  "color_teal": ("teal", "#00bfbf", "teal"),
  # Oranges
  "color_orange": ("orange", "#ff7f00", "orange"),
  "color_tv_orange": ("tv_orange", "#ff8c26", "TV orange"),
  "color_lightorange": ("lightorange", "#ffcc7f", "light orange"),
  "color_olive": ("olive", "#c4b200", "olive"),
  # Greys / Black / White
  "color_white": ("white", "#ffffff", "white"),
  "color_grey70": ("grey70", "#b2b2b2", "grey70"),
  "color_grey30": ("grey30", "#4c4c4c", "grey30"),
  "color_black": ("black", "#000000", "black"),
}


def _generate_color_help_text(color_name: str, hex_code: str, description: str) -> str:
  """Generate help text for a color button."""
  return f"<p>Color: <b>{color_name}</b> ({hex_code}) &mdash; Click to apply {description} to the selected protein or structure.</p>"


def _load_help_texts_from_files(help_base_dir: Path) -> Dict[str, str]:
  """Load help texts from external HTML files."""
  help_map = {}

  # Define the subdirectories to search
  subdirs = [
    "panels",
    "viewer_toolbar",
    "popup_menus",
    "objects_panel_toolbar",
    "tree_nodes",
    "dialogs",
    "menu_actions"
  ]

  # Load help texts from each subdirectory
  for subdir in subdirs:
    subdir_path = help_base_dir / subdir
    if not subdir_path.exists():
      logger.debug(f"Help subdirectory not found: {subdir_path}")
      continue

    # Read all HTML files in this subdirectory
    for html_file in subdir_path.glob("*.html"):
      try:
        with open(html_file, 'r', encoding='utf-8') as f:
          content = f.read().strip()

        # Use the filename (without extension) as the key
        key = html_file.stem
        help_map[key] = content

        logger.debug(f"Loaded help text for: {key}")

      except Exception as exc:
        logger.error(f"Failed to load help text from {html_file}: {exc}")

  return help_map


def load_help_texts() -> Dict[str, str]:
  """
    Load all hover help texts from external HTML files with caching.

    This function loads help texts once and caches them for subsequent calls.
    Color button help texts are generated programmatically to avoid duplication.

    Returns:
        Dict[str, str]: A dictionary mapping help keys to their HTML content.
    """
  global _help_map_cache

  # Return a cached version if already loaded
  if _help_map_cache is not None:
    return _help_map_cache

  help_map = {}

  # Determine the base path for help files
  current_file = Path(__file__)
  pyssa_root = current_file.parent.parent
  help_base_dir = pyssa_root / "data" / "help"

  if help_base_dir.exists():
    # Load help texts from files
    help_map = _load_help_texts_from_files(help_base_dir)
    logger.info(f"Loaded {len(help_map)} help texts from {help_base_dir}")
  else:
    logger.warning(f"Help directory not found: {help_base_dir}")

  # Generate color button help texts programmatically (no need for 36 separate files)
  for key, (color_name, hex_code, description) in _COLOR_DEFINITIONS.items():
    # Only generate if not already loaded from file (file takes precedence)
    if key not in help_map:
      help_map[key] = _generate_color_help_text(color_name, hex_code, description)

  # Cache the loaded help texts
  _help_map_cache = help_map

  logger.info(f"Total help texts available: {len(help_map)}")

  return help_map
