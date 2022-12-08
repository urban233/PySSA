# Add environment variable PATH
# For bash or zsh
export PATH="$HOME/.pyssa/colabfold_batch/bin:$PATH"

export TF_FORCE_UNIFIED_MEMORY="1"
export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
export TF_FORCE_GPU_ALLOW_GROWTH="true"
