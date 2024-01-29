# Define paths and file names
$batchScriptPath1 = "C:\Users\martin\Desktop\copy_pyssa_src.bat"
$batchScriptPath2 = "C:\ProgramData\pyssa\win_start\vb_script\PySSA.bat"

# Run the batch script using Start-Process
Start-Process -FilePath "cmd.exe" -ArgumentList "/c", $batchScriptPath1 -Wait
# Run the batch script using Start-Process
Start-Process -FilePath "cmd.exe" -ArgumentList "/c", $batchScriptPath2 -Wait

# Pause execution at the end
Read-Host "Press Enter to exit"
