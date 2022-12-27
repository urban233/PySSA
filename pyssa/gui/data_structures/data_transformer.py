import pathlib


class DataTransformer:
    
    def __init__(self, ui):
        self.ui = ui
    
    def transform_to_analysis(self, project):
        prot_1_name = self.ui.lbl_analysis_prot_struct_1.text().replace(".pdb", "")
        prot_2_name = self.ui.lbl_analysis_prot_struct_2.text().replace(".pdb", "")
        prot_1 = project.search_protein(prot_1_name)
        if prot_1_name == prot_2_name:
            prot_2 = prot_1.duplicate_protein()
            prot_1.molecule_object = f"{prot_1.molecule_object}_1"
            prot_2.molecule_object = f"{prot_2.molecule_object}_2"
        else:
            prot_2 = project.search_protein(prot_2_name)

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

        export_dir = pathlib.Path(
            f"{project.get_results_path()}/{prot_1.molecule_object}_vs_{prot_2.molecule_object}")
        
        return prot_1, prot_2, export_dir
