$directory_content = Get-ChildItem -Path .\pyssa\scripts -Recurse

for ($i = 0; $i -lt $directory_content.Count; $i++) {
    $script_file = $directory_content[$i]
    .\pyssa\externals\dos2unix.exe .\pyssa\scripts\$script_file
}
