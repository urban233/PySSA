; -- CodePrepareToInstall.iss --
;
; This script shows how the PrepareToInstall event function can be used to
; install prerequisites and handle any reboots in between, while remembering
; user selections across reboots.

[Setup]
WizardStyle=modern
AppName=PySSA
AppVersion=1.1.0
AppCopyright=Martin Urban, Hannah Kullik, IBCI
AppId={{192F52C3-D86D-4735-9929-C7DF599CB538}
DefaultDirName={commonappdata}\IBCI\PySSA
AppPublisher=IBCI
VersionInfoProductName=PySSA
MinVersion=10.0.19045
PrivilegesRequired=lowest
OutputDir=..\..\dist
OutputBaseFilename=pyssa_src_update_1.1.0
DisableReadyPage=True
DisableWelcomePage=False
DisableDirPage=True
DisableProgramGroupPage=True
UninstallDisplayName=PySSA
UninstallDisplayIcon={app}\assets\logo.ico
ArchitecturesInstallIn64BitMode=x64
LicenseFile=LICENSE.txt
; This is necessary because the setup will exceed 2 GB (due to almalinux rootfs)
DiskSpanning=no
;DiskSliceSize=2100000000


[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Dirs]
Name: "{commonappdata}\IBCI\wsl2\PySSA"
Name: "{app}"
Name: "{app}\assets"
Name: "{app}\bin"
Name: "{app}\third_party"

[Files]
Source: "..\..\..\inno-build-release\inno-sources\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs;
Source: "..\..\..\inno-build-release\inno-assets\logo.ico"; DestDir: "{app}\assets"; Flags: ignoreversion recursesubdirs createallsubdirs;

[Icons]
Name: "{userdesktop}\PySSA"; Filename: "{app}\pyssa.exe"; IconFilename: "{app}\assets\logo.ico"
Name: "{userstartmenu}\PySSA"; Filename: "{app}\pyssa.exe"; IconFilename: "{app}\assets\logo.ico"

[Run]
Filename: "{app}\third_party\VC_redist.x64.exe"; Parameters: "/quiet /norestart"; Flags: runhidden waituntilterminated

[UninstallRun]
Filename: "{cmd}"; Parameters: "/C wsl --unregister almaColabfold9"

[UninstallDelete]
Type: filesandordirs; Name: "{app}"

[Code]
const
  (*** Customize the following to your own name. ***)
  RunOnceName = 'My Program Setup restart';

  QuitMessageReboot = 'To complete the installation of WSL2 which is a prerequisite you will need to restart your computer. After restarting your computer, the setup will continue.';
  QuitMessageError = 'An error occurred during the WSL2 installation. Please try again.';

var
  Restarted: Boolean;
    
function Quote(const S: String): String;
begin
  Result := '"' + S + '"';
end;

function AddParam(const S, P, V: String): String;
begin
  if V <> '""' then
    Result := S + ' /' + P + '=' + V;
end;

function AddSimpleParam(const S, P: String): String;
begin
 Result := S + ' /' + P;
end;

procedure CreateRunOnceEntry;
var
  RunOnceData: String;
begin
  RunOnceData := Quote(ExpandConstant('{srcexe}')) + ' /restart=1';
  RunOnceData := AddParam(RunOnceData, 'LANG', ExpandConstant('{language}'));
  RunOnceData := AddParam(RunOnceData, 'DIR', Quote(WizardDirValue));
  RunOnceData := AddParam(RunOnceData, 'GROUP', Quote(WizardGroupValue));
  if WizardNoIcons then
    RunOnceData := AddSimpleParam(RunOnceData, 'NOICONS');
  RunOnceData := AddParam(RunOnceData, 'TYPE', Quote(WizardSetupType(False)));
  RunOnceData := AddParam(RunOnceData, 'COMPONENTS', Quote(WizardSelectedComponents(False)));
  RunOnceData := AddParam(RunOnceData, 'TASKS', Quote(WizardSelectedTasks(False)));
  
  RegWriteStringValue(HKCU, 'Software\Microsoft\Windows\CurrentVersion\RunOnce', RunOnceName, RunOnceData);
end;

function IsVirtualMachinePlatformEnabled: Boolean;
var
  Command: string;
  Output: Integer;
  FilePath: string;
begin
  Command := '"$vm = Get-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform | Select-Object -ExpandProperty State; $fileName = ''vm_platform_'' + $vm + ''.txt''; Out-File -FilePath (Join-Path ' + ExpandConstant('{tmp}') + ' $fileName);"'
  ShellExec('runas', 'powershell.exe', Command, '', SW_HIDE, ewWaitUntilTerminated, Output);
  
  FilePath := ExpandConstant('{tmp}') + '\vm_platform_Enabled.txt';
  Result := FileExists(FilePath);
end;

function DetectAndInstallPrerequisites: Boolean;
var
  ResultCode: Integer;
begin
  Result := True;
  
  if not IsVirtualMachinePlatformEnabled() then
  begin
    if MsgBox('WSL2 is not installed. Do you want to install WSL2 now?' + #13#10 + 'IMPORTANT: The WSL2 will integrate into the Windows OS and be a system component that cannot be uninstalled!', mbConfirmation, MB_YESNO) = IDNO then
    begin
      MsgBox('WSL2 is required for this installation. Setup will now exit.', mbError, MB_OK);
      Result := False;
      Exit;
    end;

    if ShellExec('runas', 'powershell.exe', 'wsl --install; Out-File -FilePath (Join-Path ' + ExpandConstant('{tmp}') + ' ''wsl2_just_installed.txt'');', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then
    begin
      if ResultCode = 0 then
      begin
        if not Restarted then
        begin
          // Create RunOnce entry and request restart
          CreateRunOnceEntry;
          Result := True;
        end
        else
        begin
          Result := IsVirtualMachinePlatformEnabled();
          if not Result then
            MsgBox('WSL2 installation did not complete successfully after reboot.', mbError, MB_OK);
        end;
      end
      else
      begin
        MsgBox('WSL2 installation failed. Error code: ' + IntToStr(ResultCode), mbError, MB_OK);
        Result := False;
      end;
    end
    else
    begin
      MsgBox('Failed to start WSL2 installation.', mbError, MB_OK);
      Result := False;
    end;
  end;
end;

function InitializeSetup(): Boolean;
begin
  Restarted := ExpandConstant('{param:restart|0}') = '1';

  if not Restarted then begin
    Result := not (RegValueExists(HKCU, 'Software\Microsoft\Windows\CurrentVersion\RunOnce', RunOnceName));
    if not Result then
      MsgBox(QuitMessageReboot, mbError, mb_Ok);
  end else
    Result := True;
end;

function PrepareToInstall(var NeedsRestart: Boolean): String;
var
  FilePath: string;

begin
  Result := '';
  
  if not DetectAndInstallPrerequisites then 
  begin
    Result := QuitMessageError;
    Exit;
  end;
  
  // If we're not restarted and WSL2 was just installed, we need to restart
  FilePath := ExpandConstant('{tmp}') + '\wsl2_just_installed.txt';
  if not Restarted and FileExists(FilePath) then
  begin
    NeedsRestart := True;
    Result := QuitMessageReboot;
  end;
end;

function ShouldSkipPage(PageID: Integer): Boolean;
begin
  Result := Restarted;
end;
