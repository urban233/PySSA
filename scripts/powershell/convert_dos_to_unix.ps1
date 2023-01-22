
# convert unix scripts
$directory_content = Get-ChildItem -Path ..\unix -Recurse
for ($i = 0; $i -lt $directory_content.Count; $i++) {
    $script_file = $directory_content[$i]
    ..\..\externals\dos2unix.exe ..\unix\$script_file
}

# convert config files
$directory_content_extras = Get-ChildItem -Path ..\..\config\wsl -Recurse
for ($i = 0; $i -lt $directory_content_extras.Count; $i++) {
    $script_file = $directory_content_extras[$i]
    ..\..\externals\dos2unix.exe ..\..\config\wsl\$script_file
}
