#!/bin/bash

export PATH="$HOME/.pyssa/localcolabfold/colabfold-conda/bin:$PATH"
cd /home/ubuntu_colabfold/.pyssa/localcolabfold/colabfold-conda/bin
pwd
./colabfold_batch --help
pwd
colabfold_batch --amber --num-recycle 3 $1 $2
