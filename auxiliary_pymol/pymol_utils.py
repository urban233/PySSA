

def remove_solvent_molecules_in_protein() -> None:
    """Removes solvent molecules in a protein."""
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.remove("solvent")
    except pymol.CmdException:
        print("No solvent molecules needs to be removed.")


def remove_organic_molecules_in_protein() -> None:
    """Removes organic molecules in a protein."""
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.remove("organic")
    except pymol.CmdException:
        print("No organic molecules needs to be removed.")


