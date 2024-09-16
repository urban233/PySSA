import enum


# <editor-fold desc="Model types">
class SeqTypesEnum(enum.StrEnum):
  """Enumeration for storing possible model types."""
  PROTEIN = "protein"
  """Seq type for storing a protein sequence."""
  NON_PROTEIN = "non_protein"
  """Seq type for storing a non protein sequence."""
# </editor-fold>
