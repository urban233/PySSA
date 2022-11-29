#!/usr/bin/env bash

# LocalColabfold
# check if curl, git and wget is installed on PC

# TODO: check operating system
type curl || { sudo dnf install curl -y;}
type git  || { sudo dnf install git -y;}
type wget || { sudo dnf install wget -y;}

# with GPU # TODO: write in installer documentation with source

# check if GNU compiler is version 4.9 or later
gcc_v=$(gcc -dumpversion)
if [ $gcc_v < "4.9" ]; then
    sudo dnf install gcc-10 -y
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10
fi

# copy install_colabbatch_linux.sh
cd scripts
pwd
cp install_colabbatch_linux.sh $HOME/.pyssa
# remove install_colabbatch_linux.sh in origin dir
# rm install_colabbatch_linux.sh
# run script
cd $HOME/.pyssa
chmod +x install_colabbatch_linux.sh
bash install_colabbatch_linux.sh

echo 'Do not move this directory after the installation! Keep the network unblocked!'

# Add environment variable PATH
# For bash or zsh
export PATH="$HOME/.pyssa/colabfold_batch/bin:$PATH"
