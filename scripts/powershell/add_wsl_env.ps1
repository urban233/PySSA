#$plugin_name = "tmpPySSA"
#$default_username = "ubuntu_colabfold"

Write-Output "Start importing ..."

wsl --import UbuntuColabfold C:\Users\$env:username\.pyssa\wsl C:\Users\$env:username\.pyssa\UbuntuColabfold.tar

Write-Output "Finished importing. Shutting down ..."

wsl -s UbuntuColabfold

wsl cp /mnt/c/Users/$env:username/AppData/Roaming/pymol/startup/$plugin_name/config/wsl/wsl.conf /etc

#wsl --shutdown

# wait for 10s so that all distros are shutdown

#wsl mkdir /home/$default_username/.pyssa

#wsl /mnt/c/Users/$env:username/AppData/Roaming/pymol/startup/$plugin_name/scripts/unix/installation_colabfold.sh

#wsl cd /home/$USER/.pyssa && wget -q -P . https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
