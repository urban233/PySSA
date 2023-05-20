#!/bin/bash

export PATH="$HOME/.pyssa/colabfold_batch/bin:$PATH"
colabfold_batch --help
pwd
colabfold_batch --num-recycle 3 $1 $2
