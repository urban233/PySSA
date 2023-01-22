# '/path/to/your/colabfold_batch' should be substituted to your path, e.g. '/home/moriwaki/Desktop/colabfold_batch'
# install GPU-supported jaxlib
COLABFOLDDIR="/home/$USER/.pyssa/colabfold_batch/"
${COLABFOLDDIR}/colabfold-conda/bin/python3.7 -m pip uninstall "colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold" -y
${COLABFOLDDIR}/colabfold-conda/bin/python3.7 -m pip uninstall alphafold-colabfold -y
${COLABFOLDDIR}/colabfold-conda/bin/python3.7 -m pip install "colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold"
${COLABFOLDDIR}/colabfold-conda/bin/python3.7 -m pip install https://storage.googleapis.com/jax-releases/cuda11/jaxlib-0.3.25+cuda11.cudnn82-cp37-cp37m-manylinux2014_x86_64.whl
${COLABFOLDDIR}/colabfold-conda/bin/python3.7 -m pip install jax==0.3.25 biopython==1.79
# fix jax.tree_(un)flatten warnings (ad hoc)
sed -i -e "s/jax.tree_flatten/jax.tree_util.tree_flatten/g" ${COLABFOLDDIR}/colabfold-conda/lib/python3.7/site-packages/alphafold/model/mapping.py
sed -i -e "s/jax.tree_unflatten/jax.tree_util.tree_unflatten/g" ${COLABFOLDDIR}/colabfold-conda/lib/python3.7/site-packages/alphafold/model/mapping.py
