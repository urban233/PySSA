FROM ubuntu:22.04

# Update the package repository and install sudo
RUN apt-get update && apt-get install -y sudo

# Create a non-root user (ubuntu_colabfold)
RUN useradd -m ubuntu_colabfold
RUN echo "ubuntu_colabfold ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/ubuntu_colabfold

# Switch to the newly created user
USER ubuntu_colabfold

# Install required packages (including curl, git, and wget)
RUN sudo apt-get update && \
    sudo apt-get install -y curl git wget

# Create a directory in the user's home directory
RUN mkdir -p /home/ubuntu_colabfold/.pyssa

# Download the script to the user's home directory
RUN wget https://w-hs.sciebo.de/s/01E5g76sjsuwowG/download -O /home/ubuntu_colabfold/install_colabbatch_linux.sh

# Make the script executable
RUN sudo chmod +x /home/ubuntu_colabfold/install_colabbatch_linux.sh

# Execute the script
RUN cd /home/ubuntu_colabfold && sudo ./install_colabbatch_linux.sh

CMD ["/bin/bash"]
