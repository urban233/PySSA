wsl --shutdown
Remove-Item -Path C:\Users\$env:username\.pyssa\wsl -recurse
wsl --unregister UbuntuColabfold
