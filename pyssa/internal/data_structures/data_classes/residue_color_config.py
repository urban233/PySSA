from dataclasses import dataclass


@dataclass
class ResidueColorConfig:
    carbon_color: str
    nitrogen_color: str
    oxygen_color: str

    def atoms_are_colored_by_elements(self) -> bool:
        if self.carbon_color == "grey70" and self.nitrogen_color == "N-blue" and self.oxygen_color == "O-red":
            return True
        return False
