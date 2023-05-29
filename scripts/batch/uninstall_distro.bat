
set DISTRO_NAME=UbuntuColabfold
set STORAGE_LOCATION=C:\Users\%USERNAME%\.pyssa\wsl\%DISTRO_NAME%

wsl --shutdown
wsl --unregister %DISTRO_NAME%
rmdir /s /q "%STORAGE_LOCATION%"
