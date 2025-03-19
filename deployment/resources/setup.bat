REM Setup WSL2 distro
wsl --import almaColabfold9 C:\ProgramData\IBCI\wsl2\PySSA C:\ProgramData\IBCI\PySSA\tmp\alma-colabfold-9-rootfs.tar
if not exist C:\ProgramData\IBCI\wsl2\PySSA\ext4.vhdx exit -1
del /F C:\ProgramData\IBCI\PySSA\tmp\alma-colabfold-9-rootfs.tar
