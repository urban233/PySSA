Podman Build Guide
==================

Build process
-------------
podman machine start
podman build -f .\Dockerfile -t alma_colabfold_9:1.0.0.0
podman run --name almaColabfold alma_colabfold_9:1.0.0.0
podman export -o alma-colabfold-9-rootfs.tar almaColabfold
podman rm almaColabfold

podman rmi alma_colabfold_9:1.0.0.0

Useful podman commands
----------------------
List all podman containers
podman ps -a
