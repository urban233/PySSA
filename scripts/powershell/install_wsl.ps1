Enable-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform -All
wsl --set-default-version 2
wsl --install -d Ubuntu
