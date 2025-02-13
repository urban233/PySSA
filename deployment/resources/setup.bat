REM Setup WSL2 distro
wsl --import almaColabfold9 C:\ProgramData\IBCI\wsl2\PySSA C:\ProgramData\IBCI\PySSA\tmp\alma-colabfold-9-rootfs.tar
if not exist C:\ProgramData\IBCI\wsl2\PySSA\ext4.vhdx exit -1
REM Setup python environment
if not exist C:\ProgramData\IBCI\PySSA\bin\.venv C:\ProgramData\IBCI\PySSA\bin\python.exe -m venv C:\ProgramData\IBCI\PySSA\bin\.venv
C:\ProgramData\IBCI\PySSA\bin\.venv\Scripts\pip.exe install -r "C:\ProgramData\IBCI\PySSA\tmp\requirements.txt" --no-index --find-links "C:\ProgramData\IBCI\PySSA\tmp\wheelfiles"
C:\ProgramData\IBCI\PySSA\bin\.venv\Scripts\pip.exe install C:\ProgramData\IBCI\PySSA\tmp\pymol-3.1.0a0-py3-none-any.whl
