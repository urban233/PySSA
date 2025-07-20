REM Setup WSL2 distro
if not exist "C:\ProgramData\IBCI\wsl2\PySSA" mkdir "C:\ProgramData\IBCI\wsl2\PySSA"
if exist C:\ProgramData\IBCI\wsl2\PySSA\ext4.vhdx exit 0
wsl --import almaColabfold9 C:\ProgramData\IBCI\wsl2\PySSA C:\ProgramData\IBCI\PySSA\tmp\alma-colabfold-9-rootfs.tar
if not exist C:\ProgramData\IBCI\wsl2\PySSA\ext4.vhdx exit -1
del /F C:\ProgramData\IBCI\PySSA\tmp\alma-colabfold-9-rootfs.tar
