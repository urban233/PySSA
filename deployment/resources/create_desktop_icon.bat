@echo off
setlocal

:: Set the path to the EXE file
set "exePath=C:\ProgramData\IBCI\PySSA\pyssa.exe"

:: Set the name of the shortcut
set "shortcutName=PySSA 2"

:: Set the path to the ICO file
set "icoPath=C:\ProgramData\IBCI\PySSA\assets\logo.ico"

:: Get the path to the desktop
set "desktopPath=%USERPROFILE%\Desktop"

:: Create a VBScript to generate the shortcut
set "vbsPath=%temp%\create_shortcut.vbs"

(
echo Set oWS = WScript.CreateObject("WScript.Shell"^)
echo sLinkFile = "%desktopPath%\%shortcutName%.lnk"
echo Set oLink = oWS.CreateShortcut(sLinkFile^)
echo oLink.TargetPath = "%exePath%"
echo oLink.WorkingDirectory = "%~dp0"
echo oLink.IconLocation = "%icoPath%"
echo oLink.Save
) > "%vbsPath%"

:: Execute the VBScript
cscript //nologo "%vbsPath%"

:: Clean up the VBScript
del "%vbsPath%"

echo Shortcut created on the desktop.
