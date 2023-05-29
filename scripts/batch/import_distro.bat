
set DISTRO_NAME=UbuntuColabfold
set TAR_FILE=C:\Users\%USERNAME%\.pyssa\UbuntuColabfold.tar
set STORAGE_LOCATION=C:\Users\%USERNAME%\.pyssa\wsl\%DISTRO_NAME%

wsl --import %DISTRO_NAME% %STORAGE_LOCATION% %TAR_FILE%
