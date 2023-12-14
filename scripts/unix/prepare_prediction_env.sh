#!/usr/bin/env bash

echo 'Preparation started.'
echo '1 - Change to home directory'
cd /home/rhel_user

echo '2 - Remove old directories'
rm -r /home/rhel_user/pyssa_colabfold
rm -r /home/rhel_user/scratch
echo '3 - Copy prediction service to WSL2 filesystem'
cp -r /mnt/c/ProgramData/pyssa/mambaforge_pyssa/pyssa-mamba-env/Lib/site-packages/pymol/pymol_path/data/startup/PySSA/pyssa_colabfold /home/rhel_user/
echo '4 - Remove original batch.py of colabfold'
sudo rm /home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages/colabfold/batch.py
echo '5 - Copy modified batch.py of colabfold'
sudo cp /mnt/c/ProgramData/pyssa/mambaforge_pyssa/pyssa-mamba-env/Lib/site-packages/pymol/pymol_path/data/startup/PySSA/pyssa_colabfold/colabfold_sub/batch.py /home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages/colabfold/batch.py
echo '6 - Create scratch directories'
mkdir /home/rhel_user/scratch
mkdir /home/rhel_user/scratch/local_predictions
mkdir /home/rhel_user/scratch/local_predictions/fasta
mkdir /home/rhel_user/scratch/local_predictions/pdb
echo 'Preparation finished.'
