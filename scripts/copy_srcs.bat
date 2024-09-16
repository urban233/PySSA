@echo off
rem Copying the new src files to the plugin directory for testing purposes
set "sourceDirectory=C:\Users\hannah\github_repos\PyDD_Prototype\src"
set "destinationDirectory=C:\Users\hannah\AppData\Local\UCSF\ChimeraX\1.7\Python311\site-packages\chimerax"
set "copy_ignore_filepath=C:\Users\hannah\github_repos\PyDD_Prototype\scripts\.copy_ignore.txt"

rmdir /s /q "%destinationDirectory%"
xcopy "%sourceDirectory%" "%destinationDirectory%" /s /e /i /h /EXCLUDE:%copy_ignore_filepath%
