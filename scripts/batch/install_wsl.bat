@echo off

echo Enabling Windows Optional Feature: Virtual Machine Platform
DISM /Online /Enable-Feature /FeatureName:VirtualMachinePlatform /All

echo Setting WSL default version to 2
wsl --set-default-version 2

echo Installing Ubuntu WSL distribution
wsl --install -d Ubuntu
