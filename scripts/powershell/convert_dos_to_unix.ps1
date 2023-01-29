
$powershell_scripts_path = $PSScriptRoot
$scripts_path = Split-Path -Path $powershell_scripts_path -Parent
$unix_scripts_path = Join-Path -Path $scripts_path -ChildPath "unix"
$tmp_PySSA_path = Split-Path -Path $scripts_path -Parent
$dos_to_unix_path = Join-Path -Path $tmp_PySSA_path -ChildPath "externals\dos2unix.exe"
$wsl_config_path = Join-Path -Path $tmp_PySSA_path -ChildPath "config\wsl\"

# convert unix scripts
$directory_content = Get-ChildItem -Path $unix_scripts_path -Recurse
for ($i = 0; $i -lt $directory_content.Count; $i++) {
    $script_file = $directory_content[$i]
    $script_path = Join-Path -Path $unix_scripts_path -ChildPath $script_file
    & $dos_to_unix_path $script_path
}

# convert config files
$directory_content_extras = Get-ChildItem -Path $wsl_config_path -Recurse
for ($i = 0; $i -lt $directory_content_extras.Count; $i++) {
    $script_file = $directory_content_extras[$i]
    $script_path = Join-Path -Path $wsl_config_path -ChildPath $script_file
    & $dos_to_unix_path $script_path
}
