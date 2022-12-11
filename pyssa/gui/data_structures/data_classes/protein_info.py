from dataclasses import dataclass


@dataclass
class ProteinInfo:
    name: str
    path: str

    def get_tuple_notation(self):
        tuple_notation: tuple[str, str] = (self.name, self.path)
        return tuple_notation
