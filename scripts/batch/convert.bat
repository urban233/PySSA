@echo off

set "directory_path=%USERPROFILE%\.pyssa\scripts\unix"
set "dos2unixExe=C:\ProgramData\pyssa\plugin\Miniconda3\envs\pyssa_colab\Lib\site-packages\pymol\pymol_path\data\startup\tmpPySSA\externals\dos2unix.exe"


for /r "%directory_path%" %%F in (*) do (
	echo File: "%%F"
    rem Do something with the file
	%dos2unixExe% "%%F"
)
pause
