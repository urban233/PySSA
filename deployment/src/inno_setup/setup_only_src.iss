; -- CodePrepareToInstall.iss --
;
; This script shows how the PrepareToInstall event function can be used to
; install prerequisites and handle any reboots in between, while remembering
; user selections across reboots.

[Setup]
WizardStyle=modern
AppName=PySSA
AppVersion=1.0.3
AppCopyright=Martin Urban, Hannah Kullik, IBCI
AppId={{192F52C3-D86D-4735-9929-C7DF599CB538}
DefaultDirName={commonappdata}\IBCI\PySSA
AppPublisher=IBCI
VersionInfoProductName=PySSA
MinVersion=10.0.19045
OutputDir=..\..\dist
OutputBaseFilename=pyssa_src_update_1.0.3
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
; Place any prerequisite files here, for example:
Source: "..\..\..\inno-build-release\inno-sources\prerequisite\WindowsTasks.exe"; Flags: dontcopy;
; Place any regular files here, so *after* all your prerequisites.
Source: "..\..\..\inno-build-release\inno-sources\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs;
Source: "..\..\..\inno-build-release\inno-assets\logo.ico"; DestDir: "{app}\assets"; Flags: ignoreversion recursesubdirs createallsubdirs;

[Icons]
Name: "{commondesktop}\PySSA"; Filename: "{app}\win_start\vb_script\window_arrangement.exe"; IconFilename: "{app}\assets\logo.ico"
Name: "{commonstartmenu}\PySSA"; Filename: "{app}\win_start\vb_script\window_arrangement.exe"; IconFilename: "{app}\assets\logo.ico"

[Run]
Filename: "{app}\third_party\VC_redist.x64.exe"; Parameters: "/quiet /norestart"; Flags: runhidden waituntilterminated
; Filename: "{app}\tmp\setup.bat"; Flags: runhidden waituntilterminated

[UninstallRun]
Filename: "{cmd}"; Parameters: "/C wsl --unregister almaColabfold9"

[UninstallDelete]
Type: filesandordirs; Name: "{commonappdata}\IBCI\PySSA"

[Code]
const
  (*** Customize the following to your own name. ***)
  RunOnceName = 'My Program Setup restart';

  QuitMessageReboot = 'To complete the installation of WSL2 which is a prerequisite you will need to restart your computer. After restarting your computer, Setup will continue next time an administrator logs in.';
  QuitMessageError = 'An error occurred during the WSL2 installation. Please try again.';

var
  Restarted: Boolean;

function InitializeSetup(): Boolean;
begin
  Restarted := ExpandConstant('{param:restart|0}') = '1';

  if not Restarted then begin
    Result := not RegValueExists(HKA, 'Software\Microsoft\Windows\CurrentVersion\RunOnce', RunOnceName);
    if not Result then
      MsgBox(QuitMessageReboot, mbError, mb_Ok);
  end else
    Result := True;
end;

function IsWSL2Installed(): Boolean;
var
  ResultCode: Integer;
begin
  Result := False;

  ExtractTemporaryFile('WindowsTasks.exe');

  if ShellExec(
    'runas',
    ExpandConstant('{tmp}\WindowsTasks.exe'),
    '',
    '',
    SW_HIDE,  // Consider hiding the window if you don't need user interaction
    ewWaitUntilTerminated,
    ResultCode) then
  begin
    // Return True only if the exit code is 0
    Result := (ResultCode = 0);
  end
  else
  begin
    // Show error if execution fails
    MsgBox('Failed to launch WindowsTasks.exe with elevated privileges: ' +
      SysErrorMessage(ResultCode), mbError, MB_OK);
    Result := False;
  end;
end;

function DetectAndInstallPrerequisites: Boolean;
var
  ResultCode: Integer;
begin
  Result := True;
  if not IsWSL2Installed() then
  begin
    // Ask user if they want to install WSL2
    if MsgBox('WSL2 is not installed. Do you want to install WSL2 now?' + #13#10 + 'IMPORTANT: The WSL2 will integrate into the Windows OS and be a system component that cannot be uninstalled!', mbConfirmation, MB_YESNO) = IDNO then
    begin
      MsgBox('WSL2 is required for this installation. Setup will now exit.', mbError, MB_OK);
      Result := False;
      Exit;
    end;

    // Proceed with installation if user selects Yes
    if Exec('wsl.exe', '--install --no-distribution', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then
    begin
      if ResultCode = 0 then
      begin
        // Schedule reboot if needed
        if not Restarted then
        begin
          RestartReplace(ParamStr(0), '');
          Result := True; // Requires reboot
        end
        else
        begin
          // Verify after reboot
          Result := IsWSL2Installed();
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

  (*** Place any custom user selection you want to remember below. ***)

  //<your code here>

  RegWriteStringValue(HKA, 'Software\Microsoft\Windows\CurrentVersion\RunOnce', RunOnceName, RunOnceData);
end;

function PrepareToInstall(var NeedsRestart: Boolean): String;
var
  ChecksumBefore, ChecksumAfter: String;
begin
  ChecksumBefore := MakePendingFileRenameOperationsChecksum;
  if DetectAndInstallPrerequisites then begin
    ChecksumAfter := MakePendingFileRenameOperationsChecksum;
    if ChecksumBefore <> ChecksumAfter then begin
      CreateRunOnceEntry;
      NeedsRestart := True;
      Result := QuitMessageReboot;
    end;
  end else
    Result := QuitMessageError;
end;

function ShouldSkipPage(PageID: Integer): Boolean;
begin
  Result := Restarted;
end;
