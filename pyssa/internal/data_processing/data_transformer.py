#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import pathlib
from pyssa.util import types


class DataTransformer:
    
    def __init__(self,
                 project,
                 proteins_for_analysis: list[tuple[types.PROTEIN_ANALYSIS_INFO, types.PROTEIN_ANALYSIS_INFO]]):
        self.project = project
        self.proteins_for_analysis = proteins_for_analysis
        print(self.project.proteins)
    # def transform_to_analysis(self, project):
    #     prot_1_name = self.ui.lbl_analysis_prot_struct_1.text().replace(".pdb", "")
    #     prot_2_name = self.ui.lbl_analysis_prot_struct_2.text().replace(".pdb", "")
    #     prot_1: protein.Protein = project.search_protein(prot_1_name)
    #     if prot_1_name == prot_2_name:
    #         prot_2: protein.Protein = prot_1.duplicate_protein()
    #         prot_1.molecule_object = f"{prot_1.molecule_object}_1"
    #         prot_2.molecule_object = f"{prot_2.molecule_object}_2"
    #     else:
    #         prot_2: protein.Protein = project.search_protein(prot_2_name)
    #
    #     prot_1_chains_selected = self.ui.list_analysis_ref_chains.selectedItems()
    #     prot_1_chains = []
    #     for tmp_chain in prot_1_chains_selected:
    #         prot_1_chains.append(tmp_chain.text())
    #     prot_1.set_chains(prot_1_chains)
    #
    #     prot_2_chains_selected = self.ui.list_analysis_model_chains.selectedItems()
    #     prot_2_chains = []
    #     for tmp_chain in prot_2_chains_selected:
    #         prot_2_chains.append(tmp_chain.text())
    #     prot_2.set_chains(prot_2_chains)
    #
    #     if len(prot_1.chains) != 0:
    #         analysis_name = f"{prot_1.molecule_object};{prot_1_chains}_vs_{prot_2.molecule_object};{prot_2_chains}"
    #         analysis_name = analysis_name.replace(";", "_")
    #         analysis_name = analysis_name.replace(",", "_")
    #         analysis_name = analysis_name.replace("[", "")
    #         analysis_name = analysis_name.replace("]", "")
    #         analysis_name = analysis_name.replace("'", "")
    #     else:
    #         analysis_name = f"{prot_1.molecule_object}_vs_{prot_2.molecule_object}"
    #     export_dir = pathlib.Path(
    #         f"{project.get_results_path()}/{analysis_name}")
    #     return prot_1, prot_2, export_dir, analysis_name

    def transform_data_for_analysis(self) -> list[tuple[types.PROTEIN, types.PROTEIN, pathlib.Path, str]]:
        """This function transforms the data from a protein list to a usable format for the analysis algorithm

        Returns:
            a list which consists of the two proteins, the export path for the results and the name of the analysis run
        """
        tmp_data = []
        for protein_pair_tuple in self.proteins_for_analysis:
            prot_1_name = protein_pair_tuple[0].protein_name.replace(".pdb", "")
            prot_2_name = protein_pair_tuple[1].protein_name.replace(".pdb", "")

            prot_1: types.PROTEIN = self.project.search_protein(prot_1_name)
            if prot_1_name == prot_2_name:
                prot_2: types.PROTEIN = prot_1.duplicate_protein()
                prot_1: types.PROTEIN = prot_2.duplicate_protein()
                prot_1.molecule_object = f"{prot_1.molecule_object}_1"
                prot_2.molecule_object = f"{prot_2.molecule_object}_2"
            else:
                prot_2: types.PROTEIN = self.project.search_protein(prot_2_name)

            prot_1_chains_selected = protein_pair_tuple[0].protein_chains
            if prot_1_chains_selected is not None:
                prot_1_chains = []
                for tmp_chain in prot_1_chains_selected:
                    prot_1_chains.append(tmp_chain)
                prot_1.set_chains(prot_1_chains)
            else:
                prot_1_chains = []

            prot_2_chains_selected = protein_pair_tuple[1].protein_chains
            if prot_2_chains_selected is not None:
                prot_2_chains = []
                for tmp_chain in prot_2_chains_selected:
                    prot_2_chains.append(tmp_chain)
                prot_2.set_chains(prot_2_chains)
            else:
                prot_2_chains = []

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
                f"{self.project.get_results_path()}/{analysis_name}")
            transformed_data: tuple = prot_1, prot_2, export_dir, analysis_name
            tmp_data.append(transformed_data)
        return tmp_data
