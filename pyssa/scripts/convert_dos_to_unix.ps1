$directory_content = Get-ChildItem -Path .\scripts -Recurse

for ($i = 0; $i -lt $directory_content.Count; $i++) {
    $script_file = $directory_content[$i]
    .\externals\dos2unix.exe .\scripts\$script_file
}

$directory_content_extras = Get-ChildItem -Path .\wsl_extras -Recurse
for ($i = 0; $i -lt $directory_content_extras.Count; $i++) {
    $script_file = $directory_content_extras[$i]
    .\externals\dos2unix.exe .\wsl_extras\$script_file
}
