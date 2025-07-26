@echo off
setlocal enabledelayedexpansion

echo Searching for IBCI\PySSA in all user AppData folders...
echo.

set found=0

REM Search in all user profiles
for /d %%u in ("C:\Users\*") do (
    set "username=%%~nu"

    REM Skip system accounts and common non-user folders
    if /i not "!username!"=="Public" (
        if /i not "!username!"=="Default" (
            if /i not "!username!"=="All Users" (

                REM Check Local AppData
                set "localpath=%%u\AppData\Local\IBCI\PySSA"
                if exist "!localpath!" (
                    echo [FOUND] Local AppData: !localpath!
                    set /a found+=1

                    REM Add your custom logic here for when found in Local AppData
                    REM Example: call :ProcessLocalAppData "!localpath!" "!username!"
                )

                REM Check Roaming AppData
                set "roamingpath=%%u\AppData\Roaming\IBCI\PySSA"
                if exist "!roamingpath!" (
                    echo [FOUND] Roaming AppData: !roamingpath!
                    set /a found+=1

                    REM Add your custom logic here for when found in Roaming AppData
                    REM Example: call :ProcessRoamingAppData "!roamingpath!" "!username!"
                )
            )
        )
    )
)

echo.
REM Conditional logic based on search results
if !found! gtr 1 (
    echo [RESULT] IBCI\PySSA directories were found in user AppData folders ^(!found! occurrence^(s^)^).

    REM Add your custom logic here for when any directories are found
    REM Example: call :HandleFoundDirectories

) else (
    echo [RESULT] No IBCI\PySSA directories found in any user AppData folders.
    wsl --unregister almaColabfold9
    REM Add your custom logic here for when no directories are found
    REM Example: call :HandleNotFound
)

echo.
echo Search completed.
goto :end

REM Example subroutines - uncomment and modify as needed
REM :ProcessLocalAppData
REM echo Processing Local AppData for user %~2: %~1
REM REM Add your specific logic here
REM goto :eof

REM :ProcessRoamingAppData
REM echo Processing Roaming AppData for user %~2: %~1
REM REM Add your specific logic here
REM goto :eof

REM :HandleFoundDirectories
REM echo Executing logic for when directories are found...
REM REM Add your specific logic here
REM goto :eof

REM :HandleNotFound
REM echo Executing logic for when no directories are found...
REM REM Add your specific logic here
REM goto :eof

:end
endlocal
pause
