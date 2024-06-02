import sys

from colabfold.batch import get_queries, run
from colabfold.download import default_data_dir
from colabfold.utils import setup_logging
from pathlib import Path


def run_prediction(fasta_dir: str, pdb_dir: str, use_amber: bool, use_templates: bool):
    input_dir = "/home/rhel_user/pyssa_colabfold/fasta"  # @param {type:"string"}
    result_dir = "/home/rhel_user/pyssa_colabfold/pdb"  # @param {type:"string"}

    input_dir = fasta_dir
    result_dir = pdb_dir

    # number of models to use
    # @markdown ---
    # @markdown ### Advanced settings
    msa_mode = "MMseqs2 (UniRef+Environmental)"  # @param ["MMseqs2 (UniRef+Environmental)", "MMseqs2 (UniRef only)","single_sequence","custom"]
    num_models = 5  # @param [1,2,3,4,5] {type:"raw"}
    num_recycles = 3  # @param [1,3,6,12,24,48] {type:"raw"}
    stop_at_score = 100  # @param {type:"string"}
    # @markdown - early stop computing models once score > threshold (avg. plddt for "structures" and ptmscore for "complexes")
    use_custom_msa = False
    num_relax = 1  # @param [0, 1, 5] {type:"raw"}
    # use_amber = num_relax > 0
    relax_max_iterations = 200  # @param [0,200,2000] {type:"raw"}
    # use_templates = False  # @param {type:"boolean"}
    do_not_overwrite_results = True  # @param {type:"boolean"}
    zip_results = False  # @param {type:"boolean"}

    # For some reason we need that to get pdbfixer to import
    if use_amber and f"/home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages" not in sys.path:
        sys.path.insert(0, f"/home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages")

    setup_logging(Path(result_dir).joinpath("log.txt"))

    queries, is_complex = get_queries(input_dir)
    run(
        queries=queries,
        result_dir=result_dir,
        use_templates=use_templates,
        num_relax=num_relax,
        relax_max_iterations=relax_max_iterations,
        msa_mode=msa_mode,
        model_type="auto",
        num_models=num_models,
        num_recycles=num_recycles,
        model_order=[1, 2, 3, 4, 5],
        is_complex=is_complex,
        data_dir=default_data_dir,
        keep_existing_results=do_not_overwrite_results,
        rank_by="auto",
        pair_mode="unpaired+paired",
        stop_at_score=stop_at_score,
        zip_results=zip_results,
        user_agent="colabfold/pyssa",
    )
