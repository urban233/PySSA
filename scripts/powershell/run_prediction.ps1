$processInfo = New-Object System.Diagnostics.ProcessStartInfo
$processInfo.FileName = "wsl.exe"
$processInfo.Arguments = "--list --verbose"
$processInfo.RedirectStandardOutput = $true
$processInfo.UseShellExecute = $false

$process = New-Object System.Diagnostics.Process
$process.StartInfo = $processInfo
$process.Start() | Out-Null
$output = $process.StandardOutput.ReadToEnd()
$process.WaitForExit()

$output
