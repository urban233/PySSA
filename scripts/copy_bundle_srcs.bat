@echo off
rem Copying the new src files to the plugin directory for testing purposes
C:\ProgramData\IBCI\PyDD\bin\ChimeraX\bin\ChimeraX-console.exe --nogui --script %USERPROFILE%\github_repos\PyDD_Prototype\scripts\update_bundle_src.cxc

rem set "sourceDirectory=%USERPROFILE%\github_repos\PyDD_Prototype\adapter_bundle"
rem set "destinationDirectory=%USERPROFILE%\AppData\Local\UCSF\ChimeraX\1.8\Python311\site-packages\PyDDAdapter-package"
rem set "copy_ignore_filepath=%USERPROFILE%\github_repos\PyDD_Prototype\scripts\.copy_ignore_bundle.txt"

rem rmdir /s /q "%destinationDirectory%"
rem xcopy "%sourceDirectory%" "%destinationDirectory%" /s /e /i /h /EXCLUDE:%copy_ignore_filepath%
