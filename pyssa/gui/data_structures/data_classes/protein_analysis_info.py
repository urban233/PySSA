from dataclasses import dataclass


@dataclass
class ProteinAnalysisInfo:
    protein_name: str
    protein_chains: str
    analysis_name: str

    def get_tuple_notation(self):
        tuple_notation: tuple[str, str] = (self.protein_name, self.protein_chains)
        return tuple_notation
