##!/usr/bin/env bash

# LocalColabfold
# check if curl, git and wget is installed on PC

type curl || { sudo apt install curl -y;}
type git  || { sudo apt install git -y;}
type wget || { sudo apt install wget -y;}

# with GPU # TODO: write in installer documentation with source

# check if GNU compiler is version 4.9 or later
#gcc_v=$(gcc -dumpversion)
#if [ $gcc_v < "4.9" ]; then
#    sudo dnf install gcc-10 -y
#    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10
#fi

# copy install_colabbatch_linux.sh
#cd scripts
echo $HOME
cd ../scripts/unix
pwd
cp install_colabbatch_linux.sh $HOME/.pyssa
cp post_colabfold_installation.sh $HOME/.pyssa
cp update.sh $HOME/.pyssa
# remove install_colabbatch_linux.sh in origin dir
# rm install_colabbatch_linux.sh
# run script
cd $HOME/.pyssa
pwd
chmod +x install_colabbatch_linux.sh
chmod +x post_colabfold_installation.sh
chmod +x update.sh
