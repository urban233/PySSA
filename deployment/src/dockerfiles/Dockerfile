# This is the latest DOCKERFILE that was used to build the rootfs that is
# located under deployment/offline_resources -> alma-colabfold-9-rootfs.tar

FROM almalinux/9-base:9

# Set environment variables
ENV PATH="/home/alma_user/localcolabfold/colabfold-conda/bin:$PATH"

# Install required dependencies
RUN dnf -y update && \
    dnf -y install git wget gcc gcc-c++ make sudo && \
    dnf clean all

# Ensure GCC version is at least 9.0
RUN gcc --version || dnf -y install gcc-toolset-9 && \
    ln -sf /opt/rh/gcc-toolset-9/root/bin/gcc /usr/bin/gcc && \
    ln -sf /opt/rh/gcc-toolset-9/root/bin/g++ /usr/bin/g++

# Create a non-root user with sudo privileges
RUN useradd -m -s /bin/bash alma_user && \
    echo 'alma_user ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Set working directory
WORKDIR /home/alma_user

# Switch to non-root user
USER alma_user

# Download and install localcolabfold (CPU version)
RUN wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh && \
    bash install_colabbatch_linux.sh && \
    rm install_colabbatch_linux.sh

# Ensure the path is set for colabfold
RUN echo 'export PATH="/home/alma_user/localcolabfold/colabfold-conda/bin:$PATH"' >> ~/.bashrc

RUN /home/alma_user/localcolabfold/colabfold-conda/bin/pip install pyzmq

# Define entrypoint
ENTRYPOINT ["/bin/bash"]
