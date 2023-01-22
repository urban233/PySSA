import pathlib

from internal.data_structures import protein
from internal.data_structures.data_classes import protein_analysis_info


class DataTransformer:
    
    def __init__(self, ui):
        self.ui = ui
    
    def transform_to_analysis(self, project):
        prot_1_name = self.ui.lbl_analysis_prot_struct_1.text().replace(".pdb", "")
        prot_2_name = self.ui.lbl_analysis_prot_struct_2.text().replace(".pdb", "")
        prot_1: protein.Protein = project.search_protein(prot_1_name)
        if prot_1_name == prot_2_name:
            prot_2: protein.Protein = prot_1.duplicate_protein()
            prot_1.molecule_object = f"{prot_1.molecule_object}_1"
            prot_2.molecule_object = f"{prot_2.molecule_object}_2"
        else:
            prot_2: protein.Protein = project.search_protein(prot_2_name)

        prot_1_chains_selected = self.ui.list_analysis_ref_chains.selectedItems()
        prot_1_chains = []
        for tmp_chain in prot_1_chains_selected:
            prot_1_chains.append(tmp_chain.text())
        prot_1.set_chains(prot_1_chains)

        prot_2_chains_selected = self.ui.list_analysis_model_chains.selectedItems()
        prot_2_chains = []
        for tmp_chain in prot_2_chains_selected:
            prot_2_chains.append(tmp_chain.text())
        prot_2.set_chains(prot_2_chains)

        if len(prot_1.chains) != 0:
            analysis_name = f"{prot_1.molecule_object};{prot_1_chains}_vs_{prot_2.molecule_object};{prot_2_chains}"
            analysis_name = analysis_name.replace(";", "_")
            analysis_name = analysis_name.replace(",", "_")
            analysis_name = analysis_name.replace("[", "")
            analysis_name = analysis_name.replace("]", "")
            analysis_name = analysis_name.replace("'", "")
        else:
            analysis_name = f"{prot_1.molecule_object}_vs_{prot_2.molecule_object}"
        export_dir = pathlib.Path(
            f"{project.get_results_path()}/{analysis_name}")
        return prot_1, prot_2, export_dir, analysis_name

    def transform_data_for_analysis(self, project, protein_list: list[tuple[
        protein_analysis_info.ProteinAnalysisInfo, protein_analysis_info.ProteinAnalysisInfo]]):
        tmp_data = []
        for protein_pair_tuple in protein_list:
            prot_1_name = protein_pair_tuple[0].protein_name.replace(".pdb", "")
            prot_2_name = protein_pair_tuple[1].protein_name.replace(".pdb", "")

            prot_1 = project.search_protein(prot_1_name)
            if prot_1_name == prot_2_name:
                prot_2 = prot_1.duplicate_protein()
                prot_1 = prot_2.duplicate_protein()
                prot_1.molecule_object = f"{prot_1.molecule_object}_1"
                prot_2.molecule_object = f"{prot_2.molecule_object}_2"
            else:
                prot_2 = project.search_protein(prot_2_name)

            prot_1_chains_selected = protein_pair_tuple[0].protein_chains
            if prot_1_chains_selected is not None:
                prot_1_chains = []
                for tmp_chain in prot_1_chains_selected:
                    prot_1_chains.append(tmp_chain)
                prot_1.set_chains(prot_1_chains)

            prot_2_chains_selected = protein_pair_tuple[1].protein_chains
            if prot_2_chains_selected is not None:
                prot_2_chains = []
                for tmp_chain in prot_2_chains_selected:
                    prot_2_chains.append(tmp_chain)
                prot_2.set_chains(prot_2_chains)

            if len(prot_1.chains) != 0:
                analysis_name = f"{prot_1.molecule_object};{prot_1_chains}_vs_{prot_2.molecule_object};{prot_2_chains}"
                analysis_name = analysis_name.replace(";", "_")
                analysis_name = analysis_name.replace(",", "_")
                analysis_name = analysis_name.replace("[", "")
                analysis_name = analysis_name.replace("]", "")
                analysis_name = analysis_name.replace("'", "")
            else:
                analysis_name = f"{prot_1.molecule_object}_vs_{prot_2.molecule_object}"
            export_dir = pathlib.Path(
                f"{project.get_results_path()}/{analysis_name}")
            transformed_data: tuple = prot_1, prot_2, export_dir, analysis_name
            tmp_data.append(transformed_data)
        return tmp_data
